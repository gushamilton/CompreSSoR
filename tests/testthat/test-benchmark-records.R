test_that("current packaged benchmark records remain readable", {
  headline <- benchmark_table()
  expect_s3_class(headline, "data.frame")
  expect_true(nrow(headline) > 0L)

  threads <- benchmark_table("pcodec_10m_threads")
  expect_identical(threads$threads, c(1L, 4L))
  expect_equal(threads$median_seconds, c(0.552, 0.451), tolerance = 1e-12)
  expect_identical(threads$rows, c(10000000L, 10000000L))
})
