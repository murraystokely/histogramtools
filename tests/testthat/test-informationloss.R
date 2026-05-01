expect_metric_invariants <- function(metric) {
  expect_true(metric >= 0)
  expect_true(metric <= 1)
}

test_that("KSDCC respects the [0,1] bound and decreases as the histogram refines", {
  set.seed(0)
  x <- rexp(100)
  h1 <- hist(x, plot = FALSE)
  h2 <- hist(x, breaks = seq(0, round(max(x) + 1), by = 0.1), plot = FALSE)

  ksdcc.1 <- KSDCC(h1)
  ksdcc.2 <- KSDCC(h2)

  expect_metric_invariants(ksdcc.1)
  expect_metric_invariants(ksdcc.2)
  expect_true(ksdcc.1 >= ksdcc.2)

  # KSDCC equals the KS statistic between the under- and over-estimate
  # representations of the binned data.
  x1.min <- rep(head(h1$breaks, -1), h1$counts)
  x1.max <- rep(tail(h1$breaks, -1), h1$counts)
  expect_equal(unname(ks.test(x1.min, x1.max, exact = FALSE)$statistic), KSDCC(h1))

  x2.min <- rep(head(h2$breaks, -1), h2$counts)
  x2.max <- rep(tail(h2$breaks, -1), h2$counts)
  expect_equal(unname(ks.test(x2.min, x2.max, exact = FALSE)$statistic), KSDCC(h2))
})

test_that("EMDCC respects [0,1] and matches emdist::emd when emdist is available", {
  set.seed(0)
  x <- rexp(100)
  h1 <- hist(x, plot = FALSE)
  h2 <- hist(x, breaks = seq(0, round(max(x) + 1), by = 0.1), plot = FALSE)

  emdcc.1 <- EMDCC(h1)
  emdcc.2 <- EMDCC(h2)

  expect_metric_invariants(emdcc.1)
  expect_metric_invariants(emdcc.2)
  expect_true(emdcc.1 >= emdcc.2)

  skip_if_not_installed("emdist")

  MinEcdf <- HistToEcdf(h1, f = 0)
  MaxEcdf <- HistToEcdf(h1, f = 1)

  A1 <- matrix(c(rep(1, length(h1$counts)),
                 h1$mids, MaxEcdf(tail(knots(MinEcdf), -1))), ncol = 3)
  A2 <- matrix(c(rep(1, length(h1$counts)),
                 h1$mids, MinEcdf(head(knots(MinEcdf), -1))), ncol = 3)
  # emdist uses single-precision floats internally, hence the 2^-23 tolerance.
  expect_equal(emdist::emd(A1, A2), emdcc.1, tolerance = 2^-23)
})
