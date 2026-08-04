# The durable v1 codec is deliberately semantic rather than a generic float
# codec: Z, EAF and SE are the sufficient GWAS values, while beta and P are
# reconstructed on read. Codes are ordinary integer columns in Parquet, so
# the store remains inspectable from Arrow/DuckDB/Polars.

q_profile_metadata <- function(se) {
  good <- is.finite(se) & se > 0
  if (!any(good)) return(list(se_log_min = -10, se_log_max = 2))
  limits <- stats::quantile(log2(se[good]), c(0.005, 0.995), names = FALSE, na.rm = TRUE)
  if (limits[1] == limits[2]) limits[2] <- limits[1] + 1
  list(se_log_min = unname(limits[1]), se_log_max = unname(limits[2]))
}

q_encode <- function(beta, se, eaf, z = NULL, z_bits = 9L, se_bits = 6L,
                     eaf_bits = 8L, block_rows = 65536L) {
  n <- length(beta)
  if (is.null(z)) z <- beta / se
  z <- as.numeric(z)
  se <- as.numeric(se)
  eaf <- as.numeric(eaf)
  z_count <- 2^as.integer(z_bits) - 2L
  se_count <- 2^as.integer(se_bits) - 2L
  eaf_max <- 2^as.integer(eaf_bits) - 1L
  z_missing <- z_count + 1L
  se_missing <- se_count + 1L
  eaf_missing <- eaf_max
  z_min <- -3.5
  z_max <- 3.5
  z_step <- (z_max - z_min) / z_count

  z_code <- rep.int(z_missing, n)
  z_ok <- is.finite(z) & z >= z_min & z <= z_max
  z_code[z_ok] <- pmin(z_count - 1L, pmax(0L, as.integer(floor((z[z_ok] - z_min) / z_step))))

  eaf_code <- rep.int(eaf_missing, n)
  eaf_ok <- is.finite(eaf) & eaf >= 0 & eaf <= 1
  eaf_code[eaf_ok] <- pmin(eaf_max - 1L,
                           pmax(0L, as.integer(round(asin(sqrt(eaf[eaf_ok])) /
                             (pi / 2) * (eaf_max - 1L)))))
  decoded_eaf <- rep(NA_real_, n)
  decoded_eaf[eaf_ok] <- sin(as.numeric(eaf_code[eaf_ok]) / (eaf_max - 1L) * (pi / 2))^2
  safe_eaf <- pmin(1 - 1e-12, pmax(1e-12, ifelse(is.finite(decoded_eaf), decoded_eaf, 0.5)))
  safe_var <- 2 * safe_eaf * (1 - safe_eaf)
  target <- ifelse(is.finite(se) & se > 0, log2(se) + 0.5 * log2(safe_var), NA_real_)

  blocks <- if (n) ceiling(n / as.integer(block_rows)) else 0L
  centers <- numeric(blocks)
  se_code <- rep.int(se_missing, n)
  se_ok <- is.finite(target) & is.finite(se) & se > 0
  for (block in seq_len(blocks)) {
    start <- (block - 1L) * as.integer(block_rows) + 1L
    stop <- min(n, block * as.integer(block_rows))
    inside <- start:stop
    good <- inside[se_ok[inside]]
    centers[block] <- if (length(good)) stats::median(target[good], na.rm = TRUE) else 0
    local <- se_ok[inside]
    residual <- (target[inside] - centers[block] + 1) / 2 * se_count - 0.5
    values <- pmin(se_count - 1L, pmax(0L, as.integer(floor(residual))))
    se_code[inside[local]] <- values[local]
  }

  # Exact exceptions are sparse and field-selective. They are used for values
  # outside the compact central ranges and for missing/invalid EAF; ordinary
  # in-range values remain quantised by design.
  z_exception <- !z_ok & is.finite(z)
  se_exception <- !se_ok & is.finite(se) & se > 0
  eaf_exception <- !eaf_ok
  exception_mask <- z_exception | se_exception | eaf_exception
  exceptions <- data.frame(
    row = which(exception_mask) - 1L,
    z_value = ifelse(z_exception[exception_mask], z[exception_mask], NA_real_),
    se_value = ifelse(se_exception[exception_mask], se[exception_mask], NA_real_),
    eaf_value = ifelse(eaf_exception[exception_mask], eaf[exception_mask], NA_real_),
    flags = as.integer(z_exception[exception_mask]) + 2L * as.integer(se_exception[exception_mask]) +
      4L * as.integer(eaf_exception[exception_mask])
  )
  exceptions <- exceptions[exceptions$flags > 0L, , drop = FALSE]
  main <- data.frame(row = seq_len(n) - 1L,
                     z_code = as.integer(z_code),
                     se_code = as.integer(se_code),
                     eaf_code = as.integer(eaf_code))
  list(main = main, exceptions = exceptions, metadata = list(
    name = "semantic_z9_eaf8_se6",
    version = "1.0.0",
    z_bits = as.integer(z_bits), eaf_bits = as.integer(eaf_bits), se_bits = as.integer(se_bits),
    z_min = z_min, z_max = z_max, z_count = z_count,
    eaf_count = eaf_max - 1L, se_count = se_count,
    se_residual_min = -1, se_residual_max = 1,
    block_rows = as.integer(block_rows), block_centers_log2_residual = unname(centers),
    z_exceptions = sum(z_exception), se_exceptions = sum(se_exception),
    eaf_exceptions = sum(eaf_exception), exception_rows = nrow(exceptions),
    p_storage = "omitted; derived as 2*pnorm(-abs(z))"
  ))
}

q_decode <- function(main, exceptions, metadata, include_beta = TRUE, include_p = FALSE) {
  n <- nrow(main)
  z_count <- as.integer(metadata$z_count %||% (2^as.integer(metadata$z_bits) - 2L))
  se_count <- as.integer(metadata$se_count %||% (2^as.integer(metadata$se_bits) - 2L))
  eaf_count <- as.integer(metadata$eaf_count %||% (2^as.integer(metadata$eaf_bits) - 1L))

  # Full scans dominate genome-wide readers and can be decoded without
  # materialising intermediate R vectors.  Regional/subset reads retain the
  # validated R path below until their row-mapping behaviour has its own
  # native benchmark.  This is a safe acceleration only: the returned names
  # and values are identical to the public decoder contract.
  full_scan <- n == nrow(main) && n > 0L &&
    identical(as.integer(main$row), seq_len(n) - 1L)
  native_available <- isTRUE(getOption("CompreSSoR.native_decode", TRUE)) &&
    is.loaded("compressor_decode_native", PACKAGE = "CompreSSoR")
  if (full_scan && isTRUE(native_available)) {
    native_exceptions <- if (nrow(exceptions)) exceptions else data.frame(
      row = integer(), z_value = numeric(), se_value = numeric(),
      eaf_value = numeric(), flags = integer()
    )
    return(.Call("compressor_decode_native",
                 as.integer(main$z_code), as.integer(main$se_code),
                 as.integer(main$eaf_code),
                 as.numeric(metadata$z_min), as.numeric(metadata$z_max),
                 z_count, se_count, eaf_count,
                 as.integer(metadata$z_bits), as.integer(metadata$se_bits),
                 as.integer(metadata$eaf_bits),
                 as.integer(metadata$block_rows %||% n),
                 as.numeric(unlist(metadata$block_centers_log2_residual %||% 0)),
                 as.integer(native_exceptions$row),
                 as.numeric(native_exceptions$z_value),
                 as.numeric(native_exceptions$se_value),
                 as.numeric(native_exceptions$eaf_value),
                 as.integer(native_exceptions$flags),
                 isTRUE(include_beta), isTRUE(include_p),
                 as.numeric(metadata$se_residual_min %||% -1),
                 as.numeric(metadata$se_residual_max %||% 1),
                 PACKAGE = "CompreSSoR"))
  }

  # The code domains are tiny.  Build the nonlinear transforms once per
  # decode and use indexed loads for every row.  This keeps the public R
  # format unchanged while avoiding row-wise sin/log2/2^ calculations.
  z_table <- rep(NA_real_, 2^as.integer(metadata$z_bits))
  z_step <- (as.numeric(metadata$z_max) - as.numeric(metadata$z_min)) / z_count
  z_table[seq_len(z_count)] <- as.numeric(metadata$z_min) +
    (seq_len(z_count) - 0.5) * z_step
  z_index <- pmin(pmax(as.integer(main$z_code) + 1L, 1L), length(z_table))
  z <- z_table[z_index]
  p_table <- 2 * stats::pnorm(-abs(z_table))

  eaf_table <- rep(NA_real_, 2^as.integer(metadata$eaf_bits))
  eaf_codes <- seq_len(eaf_count + 1L) - 1L
  eaf_table[eaf_codes + 1L] <- sin(eaf_codes / eaf_count * (pi / 2))^2
  eaf_index <- pmin(pmax(as.integer(main$eaf_code) + 1L, 1L), length(eaf_table))
  eaf <- eaf_table[eaf_index]

  block_rows <- as.integer(metadata$block_rows %||% n)
  centers <- as.numeric(unlist(metadata$block_centers_log2_residual %||% 0))
  block_id <- floor(as.numeric(main$row) / block_rows) + 1L

  # Include the fallback EAF used by the original decoder for missing EAF,
  # because SE can still be valid when EAF is not.  The table stores the
  # residual SE factor; the block centre contributes one multiplication.
  se_eaf <- eaf_table
  se_eaf[is.na(se_eaf)] <- 0.5
  safe_eaf <- pmin(1 - 1e-12, pmax(1e-12, se_eaf))
  se_min <- as.numeric(metadata$se_residual_min %||% -1)
  se_max <- as.numeric(metadata$se_residual_max %||% 1)
  se_residual <- se_min + (seq_len(se_count) - 0.5) * ((se_max - se_min) / se_count)
  eaf_correction <- -0.5 * log2(2 * safe_eaf * (1 - safe_eaf))
  se_table <- matrix(NA_real_, nrow = 2^as.integer(metadata$se_bits),
                     ncol = length(se_eaf))
  se_table[seq_len(se_count), ] <- matrix(
    2^(matrix(rep(se_residual, length(se_eaf)), nrow = se_count,
              ncol = length(se_eaf)) +
      matrix(rep(eaf_correction, each = se_count), nrow = se_count,
             ncol = length(se_eaf))),
    nrow = se_count, ncol = length(se_eaf)
  )
  se_index <- pmin(pmax(as.integer(main$se_code) + 1L, 1L), nrow(se_table))
  se <- rep(NA_real_, n)
  se_ok <- main$se_code < se_count
  if (any(se_ok) && length(centers)) {
    centre_factor <- 2^centers[pmin(block_id[se_ok], length(centers))]
    se[se_ok] <- se_table[cbind(se_index[se_ok], eaf_index[se_ok])] * centre_factor
  }

  if (nrow(exceptions)) {
    full_scan <- n == nrow(main) && n > 0L &&
      identical(as.integer(main$row), seq_len(n) - 1L)
    pos <- if (full_scan) as.integer(exceptions$row) + 1L else match(exceptions$row, main$row)
    good <- !is.na(pos)
    good <- good & pos >= 1L & pos <= n
    pos <- pos[good]
    ex <- exceptions[good, , drop = FALSE]
    z_keep <- bitwAnd(as.integer(ex$flags), 1L) != 0L
    se_keep <- bitwAnd(as.integer(ex$flags), 2L) != 0L
    eaf_keep <- bitwAnd(as.integer(ex$flags), 4L) != 0L
    z[pos[z_keep]] <- ex$z_value[z_keep]
    se[pos[se_keep]] <- ex$se_value[se_keep]
    eaf[pos[eaf_keep]] <- ex$eaf_value[eaf_keep]
  }
  out <- list(z = z, standard_error = se,
              effect_allele_frequency = eaf)
  if (isTRUE(include_beta)) out$beta <- z * se
  if (isTRUE(include_p)) {
    p <- p_table[z_index]
    if (nrow(exceptions)) {
      full_scan <- n == nrow(main) && n > 0L &&
        identical(as.integer(main$row), seq_len(n) - 1L)
      pos <- if (full_scan) as.integer(exceptions$row) + 1L else match(exceptions$row, main$row)
      good <- !is.na(pos) & pos >= 1L & pos <= n
      pos <- pos[good]
      ex <- exceptions[good, , drop = FALSE]
      z_keep <- bitwAnd(as.integer(ex$flags), 1L) != 0L
      if (any(z_keep)) p[pos[z_keep]] <- 2 * stats::pnorm(-abs(ex$z_value[z_keep]))
    }
    out$p_value <- p
  }
  out[c("z", if (isTRUE(include_beta)) "beta", "standard_error",
        "effect_allele_frequency", if (isTRUE(include_p)) "p_value")]
}
