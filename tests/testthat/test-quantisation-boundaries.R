next_down <- function(x) {
  # Largest double below x (x > 0 or x < 0, not 0).
  y <- x - abs(x) * .Machine$double.eps
  while ((x + y) / 2 < x && (x + y) / 2 > y) y <- (x + y) / 2
  y
}
next_up <- function(x) -next_down(-x)

test_that("adjacent-double helpers step by exactly one ulp", {
  expect_true(next_down(3.5) < 3.5)
  expect_identical((next_down(3.5) + 3.5) / 2, 3.5)
  expect_identical(next_down(3.5), 0.0875 / 0.025)
  expect_true(next_up(3.5) > 3.5)
  expect_true(next_up(-3.5) > -3.5)
  expect_true(next_down(-3.5) < -3.5)
})

test_that("central bin codes never collide with the missing sentinel", {
  bin <- CompreSSoR:::pcodec_native_bin_code
  z_step <- 7 / 510
  # Every Z bin edge, the values either side of it, and the range ends.
  edges <- -3.5 + (0:510) * z_step
  probe <- c(edges, vapply(edges, next_down, numeric(1)),
             vapply(edges, next_up, numeric(1)))
  probe <- probe[probe >= -3.5 & probe < 3.5]
  codes <- bin(probe, -3.5, z_step, 510L)
  expect_true(all(codes >= 0L & codes <= 509L))
  # The raw formula does hit the sentinel just below 3.5; the clamp keeps it
  # in the last bin, which contains it.
  expect_identical(as.integer(floor((next_down(3.5) + 3.5) / z_step)), 510L)
  expect_identical(bin(next_down(3.5), -3.5, z_step, 510L), 509L)
  expect_identical(bin(-3.5, -3.5, z_step, 510L), 0L)
  expect_identical(bin(next_up(-3.5), -3.5, z_step, 510L), 0L)

  se_step <- 2 / 62
  edges <- -1 + (0:62) * se_step
  probe <- c(edges, vapply(edges[edges != 0], next_down, numeric(1)),
             vapply(edges[edges != 0], next_up, numeric(1)))
  probe <- probe[probe >= -1 & probe < 1]
  codes <- bin(probe, -1, se_step, 62L)
  expect_true(all(codes >= 0L & codes <= 61L))
  expect_identical(bin(next_down(1), -1, se_step, 62L), 61L)
  expect_identical(bin(-1, -1, se_step, 62L), 0L)
})

test_that("the Pcodec quantiser keeps every finite Z at the boundary", {
  z <- c(3.5, next_down(3.5), next_up(3.5), -3.5, next_up(-3.5),
         next_down(-3.5), 0.0875 / 0.025, 0.0574 / 0.0164)
  data <- data.frame(z = z, standard_error = rep(0.02, length(z)),
                     effect_allele_frequency = rep(0.3, length(z)))
  q <- CompreSSoR:::pcodec_native_quantise(data)
  expect_false(any(q$z == 510L))
  expect_identical(q$z, c(511L, 509L, 511L, 0L, 0L, 511L, 509L, 509L))
  expect_identical(q$exceptions$row, c(0L, 2L, 5L))
  expect_identical(q$exceptions$flags, c(1L, 1L, 1L))
})

test_that("the Pcodec quantiser keeps SE residuals at the bin boundary", {
  # One SE-centre block whose median residual is the middle row. With these
  # values the third row's residual delta is the double just below +1, the
  # upper end of the SE range, and the unclamped bin formula gives 62: the SE
  # missing sentinel.
  se <- c(0.7071 / 4, 0.7071, 1.4142)
  data <- data.frame(z = c(0, 0, 0), standard_error = se,
                     effect_allele_frequency = rep(0.5, 3L))
  q <- CompreSSoR:::pcodec_native_quantise(data)
  expect_false(any(q$se == CompreSSoR:::PCODEC_NATIVE_SE_MISSING_CODE))
  expect_identical(q$se[3L], CompreSSoR:::PCODEC_NATIVE_SE_COUNT - 1L)

  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(3L)
  input$standard_error <- se
  input$beta <- c(0.01, -0.02, 0.03)
  input$effect_allele_frequency <- 0.5
  input$p_value <- 2 * pnorm(-abs(input$beta / se))
  store <- compress_sumstats(input, tempfile("se-boundary-"), overwrite = TRUE)
  full <- read_sumstats(store)
  expect_false(anyNA(full$standard_error))
  expect_false(anyNA(full$beta))
  expect_equal(full$standard_error, se, tolerance = 0.04)
})

test_that("Z, beta and p survive a store round trip at |z| = 3.5", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  z_true <- c(3.5, next_down(3.5), next_up(3.5), -3.5, next_up(-3.5),
              next_down(-3.5), 0.0875 / 0.025, -0.0875 / 0.025, 1, -1)
  n <- 64L
  input <- make_fixture(n)
  input$standard_error <- 1
  input$beta <- rep(c(0.1, -0.2), length.out = n)
  input$beta[seq_along(z_true)] <- z_true
  input$p_value <- 2 * pnorm(-abs(input$beta))
  store <- compress_sumstats(input, tempfile("z-boundary-"), overwrite = TRUE)

  for (threads in c(1L, 2L)) {
    full <- read_sumstats(store, threads = threads)
    expect_false(anyNA(full$z))
    expect_false(anyNA(full$beta))
    expect_false(anyNA(full$p_value))
    expect_false(anyNA(full$standard_error))
    got <- full$z[seq_along(z_true)]
    # The central range is [-3.5, 3.5): values outside it are exact float32
    # exceptions, the rest are quantised to within half a bin.
    outside <- z_true >= 3.5 | z_true < -3.5
    expect_equal(got[outside], z_true[outside], tolerance = 1e-6)
    expect_true(all(abs(got[!outside] - z_true[!outside]) <= 7 / 510 / 2 + 1e-12))
    expect_true(all(sign(got) == sign(z_true)))
  }

  # Row-selective reads go through a different decoder.
  by_row <- read_sumstats(store, variants = seq_along(z_true) - 1L,
                          columns = c("z", "beta", "p_value"))
  expect_false(anyNA(by_row$z))
  expect_false(anyNA(by_row$p_value))
  keys <- compressor_variant_key(input$chromosome, input$base_pair_location,
                                 input$other_allele, input$effect_allele)
  by_key <- read_sumstats(store, variants = keys[seq_along(z_true)],
                          columns = c("z", "beta", "p_value"))
  expect_false(anyNA(by_key$z))
  expect_false(anyNA(by_key$beta))
  regional <- read_sumstats(store, region = "chr1:100001-100010",
                            columns = c("z", "beta", "p_value"))
  expect_equal(nrow(regional), length(z_true))
  expect_false(anyNA(regional$z))
})

test_that("EAF round trips at 0, 1 and the adjacent doubles", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  eaf_true <- c(0, 1, next_down(1), 4.940656458412465e-324, .Machine$double.xmin, 0.5,
                next_down(0.5), next_up(0.5))
  n <- 32L
  input <- make_fixture(n)
  input$effect_allele_frequency[seq_along(eaf_true)] <- eaf_true
  store <- compress_sumstats(input, tempfile("eaf-boundary-"), overwrite = TRUE)
  full <- read_sumstats(store)
  expect_false(anyNA(full$effect_allele_frequency))
  expect_false(anyNA(full$standard_error))
  expect_equal(full$effect_allele_frequency[seq_along(eaf_true)], eaf_true,
               tolerance = 0.006)
  expect_identical(full$effect_allele_frequency[1:2], c(0, 1))
  q <- CompreSSoR:::pcodec_native_quantise(data.frame(
    z = 0, standard_error = 0.02, effect_allele_frequency = eaf_true
  ))
  expect_true(all(q$eaf >= 0L & q$eaf <= CompreSSoR:::PCODEC_NATIVE_EAF_COUNT))
})

test_that("the Parquet semantic codec also clamps the top bin", {
  encoded <- CompreSSoR:::q_encode(
    beta = c(next_down(3.5), 3.5, -3.5, next_up(-3.5)), se = rep(1, 4),
    eaf = c(0, 1, next_down(1), 0.5)
  )
  expect_identical(encoded$main$z_code, c(509L, 509L, 0L, 0L))
  expect_true(all(encoded$main$eaf_code <= 254L))
  expect_true(all(encoded$main$se_code <= 61L))
})
