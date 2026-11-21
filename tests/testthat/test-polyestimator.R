# Tests for the CDF and quantile estimators of polynomial histograms:
# the HistToEcdf, ApproxQuantile and ApproxMean methods in R/polynomial.R.

# L1 distance between a CDF function and the ECDF of x over the breaks'
# range, normalized by that range like EMDCC: exact integration of the
# absolute difference against the step ECDF, with the reconstruction
# interpolated linearly between a fine grid plus the data values.
l1_cdf <- function(F, x, breaks, n.grid = 20000) {
  grid <- sort(unique(c(seq(min(breaks), max(breaks), length.out = n.grid), x, breaks)))
  d <- abs(F(grid) - ecdf(x)(grid))
  sum(diff(grid) * (d[-1] + d[-length(d)]) / 2) / diff(range(breaks))
}

log_breaks <- function(x, K) {
  br <- exp(seq(log(min(x)), log(max(x)), length.out = K + 1))
  br[1] <- min(x)
  br[K + 1] <- max(x)
  br
}

# ----- Count-only methods are unchanged by the generics ------------------------

test_that("the count-only HistToEcdf, ApproxQuantile and ApproxMean are unchanged", {
  set.seed(0)
  x <- rexp(200)
  h <- hist(x, plot = FALSE)
  e <- HistToEcdf(h)
  expect_s3_class(e, c("ecdf", "stepfun", "function"), exact = TRUE)
  expect_equal(e(h$breaks), c(0, cumsum(h$counts) / sum(h$counts)))
  expect_equal(HistToEcdf(h, inverse = TRUE)(1), max(h$breaks))
  expect_equal(HistToEcdf(h, method = "linear")(mean(h$breaks[1:2])),
               h$counts[1] / 2 / sum(h$counts))
  expect_equal(ApproxMean(h), weighted.mean(h$mids, h$counts))
  expect_equal(suppressWarnings(ApproxQuantile(h, 0.5)),
               suppressWarnings(Hmisc::wtd.quantile(h$mids, weights = h$counts, probs = 0.5)))
  expect_warning(ApproxQuantile(h, 0.5), "fewer than 100 buckets")
  expect_error(HistToEcdf(1:3))
  expect_error(ApproxQuantile(1:3, 0.5))
  expect_error(ApproxMean(1:3))
})

# ----- CDF -------------------------------------------------------------------

test_that("the reconstructed CDF is a CDF that passes through the cumulative counts", {
  set.seed(1)
  x <- rlnorm(3000)
  br <- log_breaks(x, 10)
  for (m in 1:2) {
    h <- PolynomialHistogram(x, br, moments = m)
    F <- HistToEcdf(h)
    expect_type(F, "closure")
    expect_equal(F(c(-Inf, br[1] - 1)), c(0, 0))
    expect_equal(F(c(br[11], br[11] + 1, Inf)), c(1, 1, 1))
    g <- seq(br[1], br[11], length.out = 5001)
    expect_true(all(diff(F(g)) >= -1e-12))
    # No bucket here has an atom on an edge, so F at every break is the
    # cumulative count there.
    expect_equal(F(br), c(0, cumsum(h$counts) / sum(h$counts)))
    expect_equal(F(c(NA_real_, 1, NA_real_))[c(1, 3)], c(NA_real_, NA_real_))
  }
})

test_that("the reconstructed CDF has the stored moments in every bucket", {
  set.seed(2)
  x <- c(rgamma(4000, 0.7), rep(5, 200))
  br <- c(0, 0.1, 0.5, 1, 2, 4, 6, 10, max(x))
  h <- PolynomialHistogram(x, br, moments = 2)
  F <- HistToEcdf(h)
  N <- sum(h$counts)
  cs <- c(0, cumsum(h$counts))
  for (j in which(h$counts > 0)) {
    a <- br[j]
    b <- br[j + 1]
    u <- seq(a, b, length.out = 20001)
    # Mass in the bucket as a CDF on [a, b]: subtract the mass below, use
    # the left limit at a so an atom on the lower edge is included.
    Fb <- (F(u) - cs[j] / N) / (h$counts[j] / N)
    Fb[1] <- 0
    m <- a + sum(diff(u) * (2 - Fb[-1] - Fb[-length(Fb)]) / 2)
    expect_equal(m, h$means[j], tolerance = 1e-4)
    g <- 2 * (u - m) * (1 - Fb)
    v <- (a - m)^2 + sum(diff(u) * (g[-1] + g[-length(g)]) / 2)
    expect_equal(v, h$variances[j], tolerance = 1e-3)
  }
})

test_that("the measured L1 distance of the reconstruction never exceeds the EMDCC bound", {
  # The bound is a theorem about every CDF consistent with the stored
  # statistics; the reconstruction has them, so it can never measure
  # above the bound.  Smooth data, heavy tails, and streams with atoms
  # at and inside bucket edges.
  set.seed(3)
  streams <- list(
    lognormal = rlnorm(5000),
    pareto = 1 / runif(5000)^(1 / 1.5),
    atoms = c(rlnorm(2000), rep(c(1, 2, 4, 7), c(500, 300, 200, 100))),
    discrete = rpois(3000, 20) + 1,
    uniform = runif(4000, 1, 100))
  for (nm in names(streams)) {
    x <- streams[[nm]]
    for (K in c(6, 24)) {
      br <- log_breaks(x, K)
      for (m in 0:2) {
        h <- PolynomialHistogram(x, br, moments = m, right = FALSE)
        F <- if (m == 0) HistToEcdf(h, method = "linear") else HistToEcdf(h)
        measured <- l1_cdf(F, x, br)
        expect_lte(measured, EMDCC(h) * (1 + 1e-9) + 1e-12,
                   label = sprintf("%s K=%d moments=%d measured", nm, K, m))
      }
    }
  }
})

test_that("a bucket holding a single value or edge atoms is reconstructed exactly", {
  # Values on the lower edges of [a, b) buckets: with the mean on an
  # edge even one moment is exact, and the CDF is the data's ECDF.
  x <- rep(c(1, 2, 4, 8), c(10, 20, 30, 40))
  br <- c(1, 2, 4, 8, 16)
  for (m in 1:2) {
    h <- PolynomialHistogram(x, br, moments = m, right = FALSE)
    F <- HistToEcdf(h)
    expect_equal(F(c(0.5, 1, 1.5, 2, 3, 4, 7.9, 8, 100)), ecdf(x)(c(0.5, 1, 1.5, 2, 3, 4, 7.9, 8, 100)))
    expect_equal(EMDCC(h), 0)
    expect_equal(unname(ApproxQuantile(h, c(0, 0.1, 0.1001, 0.3, 0.5, 0.95, 1))),
                 c(1, 1, 2, 2, 4, 8, 8))
  }
  # Values in the interior: exact with two moments only.
  y <- rep(c(1.5, 3, 5), c(5, 3, 2))
  br <- c(1, 2, 4, 8)
  h2 <- PolynomialHistogram(y, br, moments = 2)
  expect_equal(HistToEcdf(h2)(c(1, 1.49, 1.5, 2.9, 3, 4.99, 5)), c(0, 0, 0.5, 0.5, 0.8, 0.8, 1))
  expect_equal(unname(ApproxQuantile(h2, c(0.2, 0.5, 0.7, 0.9, 1))), c(1.5, 1.5, 3, 5, 5))
  h1 <- PolynomialHistogram(y, br, moments = 1)
  expect_gt(abs(HistToEcdf(h1)(1.49) - 0), 0.01)
})

# ----- Quantiles ---------------------------------------------------------------

test_that("ApproxQuantile is monotone, inverts the CDF, and stays inside the bound", {
  set.seed(4)
  x <- c(rlnorm(5000, 0, 1.5), rep(3, 400))
  br <- log_breaks(x, 16)
  p <- c(0, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1)
  for (m in 1:2) {
    h <- PolynomialHistogram(x, br, moments = m)
    q <- ApproxQuantile(h, p)
    expect_named(q, paste0(formatC(100 * p, format = "fg", width = 1,
                                   digits = max(2L, getOption("digits"))), "%"))
    expect_true(all(diff(q) >= 0))
    expect_equal(unname(q[c(1, length(p))]), c(br[1], br[17]))
    F <- HistToEcdf(h)
    Q <- HistToEcdf(h, inverse = TRUE)
    expect_equal(unname(q), Q(p))
    # Where the reconstruction is continuous, F(Q(p)) = p; at an atom it
    # can only overshoot, by at most the atom's mass.
    expect_true(all(F(q) >= p - 1e-9))
    inner <- p > 0 & p < 1
    expect_true(all(F(q[inner]) - p[inner] <= 400 / length(x) + 1e-9))
    # Rank error against the data never exceeds the envelope width,
    # which is at most the mass of the bucket holding the quantile.
    N <- length(x)
    for (i in which(inner)) {
      j <- findInterval(q[i], br, rightmost.closed = TRUE)
      expect_lte(abs(ecdf(x)(q[i]) - p[i]), h$counts[j] / N + 1e-9)
    }
  }
})

test_that("stored moments sharpen the tail quantiles on skewed data", {
  # Fixed seed and wide buckets: the count-only estimate interpolates
  # linearly across a bucket that holds a decaying tail, while the
  # moment-annotated fits follow it.  Checked at the 90th and 99th
  # percentiles as a median over 20 replicates.
  set.seed(5)
  br <- c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 64, 256)
  err <- replicate(20, {
    x <- rlnorm(5000)
    p <- c(0.9, 0.99)
    tr <- quantile(x, p)
    c(count = mean(abs(unname(suppressWarnings(
            ApproxQuantile(PolynomialHistogram(x, br, moments = 0), p))) - tr) / tr),
      one = mean(abs(unname(ApproxQuantile(PolynomialHistogram(x, br, moments = 1), p)) - tr) / tr),
      two = mean(abs(unname(ApproxQuantile(PolynomialHistogram(x, br, moments = 2), p)) - tr) / tr))
  })
  med <- apply(err, 1, median)
  expect_lt(med["one"], med["count"])
  expect_lt(med["two"], med["one"])
})

test_that("ApproxQuantile handles empty buckets, NA probabilities and bad input", {
  h <- PreBinnedPolynomialHistogram(c(0, 1, 2, 3, 4), c(0, 5, 0, 5), c(NA, 1.5, NA, 3.5),
                                    c(NA, 0.05, NA, 0.05))
  q <- ApproxQuantile(h, c(0, 0.5, 0.5000001, 1, NA))
  expect_equal(unname(q[1]), 1)                      # the first non-empty bucket
  expect_lte(unname(q[2]), 2)
  expect_gte(unname(q[3]), 3)                        # skips the empty bucket
  expect_equal(unname(q[4]), 4)
  expect_equal(unname(q[5]), NA_real_)
  F <- HistToEcdf(h)
  expect_equal(F(c(0.5, 2, 2.5, 3)), c(0, 0.5, 0.5, 0.5))
  expect_error(ApproxQuantile(h, 1.5))
  expect_error(ApproxQuantile(h, -0.1))
  expect_error(HistToEcdf(PreBinnedPolynomialHistogram(c(0, 1), 0, NA)))
})

# ----- Mean ------------------------------------------------------------------

test_that("ApproxMean of a polynomial histogram is the exact data mean", {
  set.seed(6)
  x <- rgamma(1000, 3) * 1e6
  br <- c(0, 1e6, 2e6, 4e6, 8e6, max(x))
  h <- PolynomialHistogram(x, br, moments = 1)
  expect_equal(ApproxMean(h), mean(x))
  expect_gt(abs(ApproxMean(hist(x, breaks = br, plot = FALSE)) - mean(x)), 1e4)
  h.empty <- PreBinnedPolynomialHistogram(c(0, 1, 2), c(0, 4), c(NA, 1.5))
  expect_equal(ApproxMean(h.empty), 1.5)
})
