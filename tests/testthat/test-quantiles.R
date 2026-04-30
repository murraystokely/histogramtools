test_that("ApproxMean computes the bin-midpoint-weighted mean", {
  x <- hist(c(1, 2, 3, 4), plot = FALSE)
  # x$breaks = 1 2 3 4 ; x$counts = 2 1 1
  expect_equal(ApproxMean(x), 2.25)
})

test_that("ApproxQuantile is monotonic and brackets the median bin", {
  x <- hist(c(1, 2, 3, 4), plot = FALSE)
  # The bucketing means a reasonable approximation for the median is in [1.5, 2.5].
  expect_true(unname(ApproxQuantile(x, .5)) >= 1.5)
  expect_true(unname(ApproxQuantile(x, .5)) <= 2.5)
  expect_true(ApproxQuantile(x, .5) < ApproxQuantile(x, .51))
  # The midpoint of the last bucket is the approximation of the 100th percentile.
  expect_equal(unname(ApproxQuantile(x, 1)), max(x$mids))
})
