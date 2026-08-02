q_encode <- function(beta, se, eaf, z_bits = 9L, se_bits = 10L, se_lo, se_hi) {
  n <- length(beta)
  z <- beta / se
  zmax <- 2^(z_bits - 1L) - 1L
  semax <- 2^(se_bits - 1L) - 1L
  sentinel <- -32768L
  z_code <- rep.int(sentinel, n)
  se_code <- rep.int(sentinel, n)
  z_ok <- is.finite(z) & z >= -3.5 & z <= 3.5
  se_log <- log(se)
  se_ok <- is.finite(se_log) & se_log >= se_lo & se_log <= se_hi
  z_code[z_ok] <- as.integer(round((z[z_ok] + 3.5) * (2 * zmax) / 7 - zmax))
  se_code[se_ok] <- as.integer(round((se_log[se_ok] - se_lo) * (2 * semax) / (se_hi - se_lo) - semax))
  eaf_code <- rep.int(65535L, n)
  eaf_ok <- is.finite(eaf) & eaf >= 0 & eaf <= 1
  eaf_code[eaf_ok] <- as.integer(round(asin(sqrt(eaf[eaf_ok])) / (pi / 2) * 4094))
  main <- data.frame(row = seq_len(n) - 1L, z_code = z_code, se_code = se_code, eaf_code = eaf_code)
  exceptions <- data.frame(
    row = seq_len(n) - 1L,
    z_value = as.numeric(ifelse(is.finite(z) & !z_ok, z, NA_real_)),
    se_value = as.numeric(ifelse(is.finite(se) & se > 0 & !se_ok, se, NA_real_))
  )
  exceptions <- exceptions[is.finite(exceptions$z_value) | is.finite(exceptions$se_value), , drop = FALSE]
  list(main = main, exceptions = exceptions, metadata = list(
    z_bits = z_bits, se_bits = se_bits, z_min = -3.5, z_max = 3.5,
    se_log_min = se_lo, se_log_max = se_hi, eaf_bits = 12L,
    z_exceptions = sum(is.finite(z) & !z_ok), se_exceptions = sum(is.finite(se) & se > 0 & !se_ok)
  ))
}

q_decode <- function(main, exceptions, metadata) {
  zmax <- 2^(metadata$z_bits - 1L) - 1L
  semax <- 2^(metadata$se_bits - 1L) - 1L
  z <- ((as.numeric(main$z_code) + zmax) * 7 / (2 * zmax)) - 3.5
  z[main$z_code == -32768L] <- NA_real_
  se <- exp(((as.numeric(main$se_code) + semax) * (metadata$se_log_max - metadata$se_log_min) / (2 * semax)) + metadata$se_log_min)
  se[main$se_code == -32768L] <- NA_real_
  if (nrow(exceptions)) {
    pos <- match(exceptions$row, main$row)
    z_keep <- !is.na(exceptions$z_value) & !is.na(pos)
    se_keep <- !is.na(exceptions$se_value) & !is.na(pos)
    z[pos[z_keep]] <- exceptions$z_value[z_keep]
    se[pos[se_keep]] <- exceptions$se_value[se_keep]
  }
  eaf <- rep(NA_real_, nrow(main))
  ok <- main$eaf_code != 65535L
  eaf[ok] <- sin(as.numeric(main$eaf_code[ok]) / 4094 * (pi / 2))^2
  list(beta = z * se, standard_error = se, effect_allele_frequency = eaf)
}

q_profile_metadata <- function(se) {
  good <- is.finite(se) & se > 0
  if (!any(good)) return(list(se_log_min = -10, se_log_max = 2))
  limits <- stats::quantile(log(se[good]), c(0.005, 0.995), names = FALSE, na.rm = TRUE)
  if (limits[1] == limits[2]) limits[2] <- limits[1] + 1
  list(se_log_min = unname(limits[1]), se_log_max = unname(limits[2]))
}
