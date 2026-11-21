# Tests for the maximum-entropy in-bucket reconstruction (R/maxent.R).
# These exercise the internal solver directly; the user-facing methods
# are covered in test-polynomial.R.

# ----- Gauss-Legendre table ---------------------------------------------------

test_that("the tabulated 16-point Gauss-Legendre rule on [0, 1] is exact to degree 31", {
  expect_length(.GL16.X, 16)
  expect_length(.GL16.W, 16)
  expect_equal(sum(.GL16.W), 1, tolerance = 1e-15)
  # Symmetric about 1/2.
  expect_equal(.GL16.X + rev(.GL16.X), rep(1, 16), tolerance = 1e-15)
  expect_equal(.GL16.W, rev(.GL16.W), tolerance = 1e-15)
  # Integrates u^k over [0, 1] exactly for k <= 2n - 1 = 31.
  for (k in c(1, 2, 7, 16, 25, 31)) {
    expect_equal(sum(.GL16.W * .GL16.X^k), 1 / (k + 1), tolerance = 1e-14)
  }
})

test_that("the tabulated rule agrees with statmod::gauss.quad", {
  skip_if_not_installed("statmod")
  g <- statmod::gauss.quad(16L, kind = "legendre")
  # statmod's eigen-solver output is symmetric only to about 1e-16; the
  # table was symmetrized, so agreement is to that level, not exact.
  expect_equal(.GL16.X, 0.5 * (g$nodes + 1), tolerance = 1e-14)
  expect_equal(.GL16.W, 0.5 * g$weights, tolerance = 1e-14)
})

# ----- One moment: closed forms ----------------------------------------------

test_that("lambda(alpha) round-trips through the mean of exp(l u) / Z", {
  alphas <- c(1e-9, 1e-5, 0.01, 0.2, 0.5, 0.8, 0.99, 1 - 1e-5, 1 - 1e-9)
  for (a in alphas) {
    l <- .LambdaOfAlpha(a)
    expect_true(is.finite(l))
    # Absolute: the mean is 1 / |l| + O(e^-|l|) at the edges, where double
    # precision resolves it to about 1e-16, not to a relative 1e-12.
    expect_lt(abs(.AlphaOfLambda(l) - a), 1e-15)
  }
  expect_equal(.LambdaOfAlpha(0.5), 0)
  # Antisymmetric about 1/2.
  expect_equal(.LambdaOfAlpha(0.3), -.LambdaOfAlpha(0.7))
})

test_that("the one-moment CDF is a valid CDF with the stored mean, and its inverse inverts it", {
  u <- seq(0, 1, length.out = 1001)
  for (l in c(-2000, -50, -3, -1e-12, 0, 0.5, 20, 700)) {
    F <- .MaxEnt1Cdf(u, l)
    expect_equal(F[1], 0)
    expect_equal(F[length(F)], 1)
    expect_true(all(diff(F) >= 0))
    # Mean of the density as the integral of 1 - F.
    m <- integrate(function(t) 1 - .MaxEnt1Cdf(t, l), 0, 1, rel.tol = 1e-10,
                   subdivisions = 2000L)$value
    expect_equal(m, .AlphaOfLambda(l), tolerance = 1e-8)
    p <- c(1e-6, 0.1, 0.5, 0.9, 1 - 1e-6)
    expect_equal(.MaxEnt1Cdf(.MaxEnt1Quantile(p, l), l), p, tolerance = 1e-10)
  }
})

# ----- Two moments: solver ----------------------------------------------------

test_that("the two-moment fit matches the stored mean and variance across the feasible region", {
  # alpha from a near point mass at the edge to the center; tau as a
  # fraction of its maximum alpha (1 - alpha) from a near point mass to
  # a near pair of edge atoms.  Mirrors the sweep in the paper's
  # verification script.
  alphas <- c(1e-10, 1e-6, 1e-3, 0.1, 0.3, 0.5)
  fracs <- c(1e-12, 1e-8, 1e-4, 1e-2, 0.3, 0.9, 0.999, 1 - 1e-8)
  u <- seq(0, 1, length.out = 2001)
  for (a in alphas) {
    for (f in fracs) {
      tau <- f * a * (1 - a)
      fit <- expect_no_warning(.MaxEnt2Solve(a, tau))
      r <- .MaxEnt2Residual(fit, a, tau)
      floor <- .MaxEnt2Floor(fit)
      expect_lt(abs(r[1]), max(.ME.TOL.MEAN, floor))
      expect_lt(abs(r[2]), max(.ME.TOL.VAR, floor))
      F <- .MaxEnt2Cdf(fit, u)
      expect_equal(F[1], 0)
      expect_equal(F[length(F)], 1)
      expect_true(all(diff(F) >= -1e-12))
    }
  }
})

test_that("the two-moment CDF equals the truncated-normal closed form when l2 < 0", {
  # f(u) = exp(l2 (u - v)^2) / Z is a normal with mean v and
  # sd 1 / sqrt(-2 l2) truncated to [0, 1].
  u <- seq(0, 1, length.out = 501)
  for (case in list(c(0.3, 0.01), c(0.5, 0.05), c(0.1, 0.002), c(0.45, 0.0005))) {
    a <- case[1]
    tau <- case[2]
    fit <- .MaxEnt2Solve(a, tau)
    expect_equal(fit$form, "vertex")
    v <- fit$p1
    s <- 1 / sqrt(-2 * fit$l2)
    Ftn <- (pnorm((u - v) / s) - pnorm(-v / s)) / (pnorm((1 - v) / s) - pnorm(-v / s))
    expect_equal(.MaxEnt2Cdf(fit, u), Ftn, tolerance = 1e-12)
  }
})

test_that("the two-moment CDF agrees with adaptive quadrature when l2 > 0", {
  # Large variance relative to the one-moment fit gives exponential
  # layers at both edges (l2 > 0), where no closed form exists.
  u <- c(0.001, 0.05, 0.3, 0.5, 0.8, 0.99)
  for (case in list(c(0.3, 0.17), c(0.5, 0.2), c(0.2, 0.12))) {
    a <- case[1]
    tau <- case[2]
    fit <- .MaxEnt2Solve(a, tau)
    expect_equal(fit$form, "layers")
    dens <- function(t) exp(.MaxEnt2Exponent(fit, t))
    Z <- integrate(dens, 0, 1, rel.tol = 1e-12)$value
    Fq <- vapply(u, function(x) integrate(dens, 0, x, rel.tol = 1e-12)$value / Z,
                 numeric(1))
    expect_equal(.MaxEnt2Cdf(fit, u), Fq, tolerance = 1e-9)
  }
})

test_that("a two-moment fit with the one-moment variance reduces to the one-moment density", {
  a <- 0.3
  f0 <- .MaxEnt2Fit(.LambdaOfAlpha(a), 0)
  fit <- .MaxEnt2Solve(a, f0$c2)
  u <- seq(0, 1, length.out = 201)
  expect_equal(.MaxEnt2Cdf(fit, u), .MaxEnt1Cdf(u, .LambdaOfAlpha(a)), tolerance = 1e-7)
})

# ----- Per-bucket dispatch ----------------------------------------------------

make.h <- function(means, variances = NULL, counts = NULL) {
  # Two buckets [0, 1) and [1, 3): the second is the one under test.
  if (is.null(counts)) counts <- c(10, 10)
  list(breaks = c(0, 1, 3), counts = counts, means = c(0.5, means),
       variances = if (is.null(variances)) NULL else c(1 / 12, variances))
}

test_that("degenerate stored moments are reconstructed as atoms, others as densities", {
  # Mean at an edge: a point mass there, with or without a variance.
  d <- .PolyBucketDescription(make.h(1), 2)
  expect_equal(d$kind, "atom")
  expect_equal(d$loc, 1)
  d <- .PolyBucketDescription(make.h(3, 0), 2)
  expect_equal(d$kind, "atom")
  expect_equal(d$loc, 3)
  # Zero variance: a point mass at the mean.
  d <- .PolyBucketDescription(make.h(2.2, 0), 2)
  expect_equal(d, list(kind = "atom", loc = 2.2, wt = 1))
  # Maximum variance alpha (1 - alpha) w^2: atoms at both edges with
  # masses 1 - alpha and alpha.
  alpha <- 0.25
  d <- .PolyBucketDescription(make.h(1 + 2 * alpha, alpha * (1 - alpha) * 4), 2)
  expect_equal(d$kind, "atom")
  expect_equal(d$loc, c(1, 3))
  expect_equal(d$wt, c(1 - alpha, alpha))
  # Interior mean, no variance: a one-moment density with the right sign.
  d <- .PolyBucketDescription(make.h(1.5), 2)
  expect_equal(d$kind, "cont")
  expect_lt(d$l, 0)
  d <- .PolyBucketDescription(make.h(2.5), 2)
  expect_gt(d$l, 0)
  # Interior mean and variance: a two-moment fit, mirrored when the mean
  # is in the upper half.
  d <- .PolyBucketDescription(make.h(2.5, 0.1), 2)
  expect_equal(d$kind, "cont")
  expect_true(isTRUE(d$fit$mirror))
  d <- .PolyBucketDescription(make.h(1.5, 0.1), 2)
  expect_false(isTRUE(d$fit$mirror))
  # An empty bucket is filled uniformly.
  d <- .PolyBucketDescription(make.h(NA, NA, counts = c(10, 0)), 2)
  expect_equal(d$kind, "cont")
  expect_equal(.PolyBucketCdf(d, c(0, 0.25, 1), 1, 3), c(0, 0.25, 1))
})

test_that("the within-bucket CDF has the stored moments and the quantile inverts it", {
  a <- 1
  b <- 3
  u <- seq(0, 1, length.out = 4001)
  for (h in list(make.h(1.4), make.h(2.7), make.h(1.4, 0.05), make.h(2.1, 0.3))) {
    d <- .PolyBucketDescription(h, 2)
    F <- .PolyBucketCdf(d, u, a, b)
    expect_equal(F[1], 0)
    expect_equal(F[length(F)], 1)
    # Mean from the CDF by the trapezoid rule, in x units.
    m <- a + (b - a) * sum(diff(u) * (2 - F[-1] - F[-length(F)]) / 2)
    expect_equal(m, h$means[2], tolerance = 1e-5)
    if (!is.null(h$variances)) {
      # Variance by parts: E[(x - m)^2] = (a - m)^2 + int 2 (x - m) (1 - F) dx.
      x <- a + (b - a) * u
      g <- 2 * (x - m) * (1 - F)
      v <- (a - m)^2 + sum(diff(x) * (g[-1] + g[-length(g)]) / 2)
      expect_equal(v, h$variances[2], tolerance = 1e-4)
    }
    for (f in c(1e-6, 0.05, 0.5, 0.95, 1 - 1e-6)) {
      q <- .PolyBucketQuantile(d, f, a, b)
      expect_true(q >= a && q <= b)
      expect_lt(abs(.PolyBucketCdf(d, (q - a) / (b - a), a, b) - f), 1e-9)
    }
    expect_equal(.PolyBucketQuantile(d, 0, a, b), a)
    expect_equal(.PolyBucketQuantile(d, 1, a, b), b)
  }
})

test_that("atom buckets give step CDFs and the atom locations as quantiles", {
  alpha <- 0.25
  h <- make.h(1 + 2 * alpha, alpha * (1 - alpha) * 4)
  d <- .PolyBucketDescription(h, 2)
  expect_equal(.PolyBucketCdf(d, c(0, 0.5, 1), 1, 3), c(1 - alpha, 1 - alpha, 1))
  expect_equal(.PolyBucketQuantile(d, 0.5, 1, 3), 1)
  expect_equal(.PolyBucketQuantile(d, 0.9, 1, 3), 3)
  d <- .PolyBucketDescription(make.h(2.2, 0), 2)
  expect_equal(.PolyBucketCdf(d, c(0, 0.59, 0.6, 1), 1, 3), c(0, 0, 1, 1))
  expect_equal(.PolyBucketQuantile(d, 0.3, 1, 3), 2.2)
})

test_that("a steep one-moment bucket returns finite quantiles at both ends", {
  # Mean within 1e-9 of the upper edge but outside the degenerate
  # tolerance: |lambda| is in the hundreds of millions.
  h <- make.h(3 - 2e-9)
  d <- .PolyBucketDescription(h, 2)
  expect_equal(d$kind, "cont")
  expect_equal(.PolyBucketQuantile(d, 1, 1, 3), 3)
  expect_equal(.PolyBucketQuantile(d, 0, 1, 3), 1)
  expect_true(is.finite(.PolyBucketQuantile(d, 0.5, 1, 3)))
})
