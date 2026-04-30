test_that("HistToEcdf inverse round-trips at the inverse function's knots", {
  set.seed(0)
  x <- rexp(100)
  h <- hist(x, plot = FALSE)
  f <- HistToEcdf(h)
  f.inv <- HistToEcdf(h, inverse = TRUE)
  expect_identical(
    knots(f.inv),
    sapply(knots(f.inv), function(x) f(f.inv(x)))
  )
})

test_that("HistToEcdf brackets the data range and matches cumulative counts", {
  set.seed(0)
  x <- rexp(100)
  h <- hist(x, plot = FALSE)
  f <- HistToEcdf(h)
  expect_true(min(knots(f)) < min(range(x)))
  expect_true(max(knots(f)) > max(range(x)))

  expect_equal(f(h$breaks[3]), (cumsum(h$counts) / sum(h$counts))[2])
})
