# Tests for polynomial histograms (R/polynomial.R): the class, its two
# constructors, the per-bucket closed-form EMDCC, and the EMDCC generic.

with_null_pdf <- function(expr) {
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file)
  on.exit({
    dev.off()
    unlink(pdf_file)
  })
  force(expr)
}

# Population variance of each bucket's observations, by hand.
bucket_stat <- function(x, k, K, f) {
  out <- rep(NA_real_, K)
  s <- tapply(x, factor(k, levels = seq_len(K)), f)
  out[!is.na(s)] <- s[!is.na(s)]
  out
}

# ----- Constructors from raw data --------------------------------------------

test_that("PolynomialHistogram with moments = 0 is exactly hist()", {
  set.seed(1)
  x <- rexp(500)
  br <- c(0, 0.5, 1, 2, 4, 8)
  expect_identical(PolynomialHistogram(x, br, moments = 0),
                   hist(x, breaks = br, plot = FALSE))
  expect_identical(PolynomialHistogram(x, moments = 0),
                   hist(x, plot = FALSE))
  expect_identical(PolynomialHistogram(x, br, moments = 0, right = FALSE),
                   hist(x, breaks = br, right = FALSE, plot = FALSE))
})

test_that("PolynomialHistogram stores the per-bucket mean and variance of hist()'s buckets", {
  set.seed(2)
  x <- c(rlnorm(2000), 7, 7, 7)
  br <- c(0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)
  K <- length(br) - 1
  for (right in c(TRUE, FALSE)) {
    h <- PolynomialHistogram(x, br, moments = 2, right = right)
    h0 <- hist(x, breaks = br, right = right, plot = FALSE)
    expect_s3_class(h, c("polynomial.histogram", "histogram"), exact = TRUE)
    expect_equal(h$breaks, h0$breaks)
    expect_equal(h$counts, h0$counts)
    expect_equal(h$density, h0$density)
    expect_equal(h$mids, h0$mids)
    expect_equal(h$xname, "x")
    # Bucket membership by the convention hist() used.
    k <- if (right) findInterval(x, br, left.open = TRUE, rightmost.closed = TRUE)
         else findInterval(x, br, rightmost.closed = TRUE)
    expect_equal(h$means, bucket_stat(x, k, K, mean))
    expect_equal(h$variances,
                 bucket_stat(x, k, K, function(v) mean((v - mean(v))^2)))
    expect_true(all(h$means >= head(br, -1) & h$means <= tail(br, -1), na.rm = TRUE))
    # The pooled bucket means give the exact data mean.
    expect_equal(sum(h$counts * h$means, na.rm = TRUE) / sum(h$counts), mean(x))
  }
  h1 <- PolynomialHistogram(x, br, moments = 1)
  expect_null(h1$variances)
  expect_equal(h1$means, PolynomialHistogram(x, br, moments = 2)$means)
})

test_that("observations on a break are assigned as hist() assigns them", {
  # hist() bins on fuzzed breaks; the moments must follow the counts.
  x <- c(1, 1, 2, 2, 2, 3, 4, 4)
  br <- 1:4
  for (right in c(TRUE, FALSE)) {
    for (il in c(TRUE, FALSE)) {
      if (!il && ((right && any(x == 1)) || (!right && any(x == 4)))) next
      h0 <- hist(x, breaks = br, right = right, include.lowest = il, plot = FALSE)
      h <- PolynomialHistogram(x, br, moments = 2, right = right, include.lowest = il)
      expect_equal(h$counts, h0$counts)
      # Every bucket's mean is the mean of the values hist() counted there:
      # with (a, b] buckets the 2s are in bucket 1 and the 3 in bucket 2.
      if (right) {
        expect_equal(h$means, c(mean(c(1, 1, 2, 2, 2)), 3, 4))
        expect_equal(h$variances, c(mean((c(1, 1, 2, 2, 2) - 1.6)^2), 0, 0))
      } else {
        expect_equal(h$means, c(1, 2, mean(c(3, 4, 4))))
      }
    }
  }
})

test_that("a value that hist()'s fuzz counts in the bucket below its break gets an edge mean", {
  # 1 - 0.9 is 2e-17 below 0.1, so 0.1 lies just past that break; hist()
  # counts it in the lower bucket (its fuzz is 1e-7 of a bucket width).
  br <- c(0, 1 - 0.9, 0.2, 0.3)
  x <- c(0.1, 0.1, 0.25)
  h0 <- hist(x, breaks = br, plot = FALSE)
  expect_equal(h0$counts, c(2, 0, 1))
  h <- PolynomialHistogram(x, br, moments = 2)
  expect_equal(h$counts, h0$counts)
  expect_equal(h$means, c(br[2], NA, 0.25))
  expect_equal(h$variances, c(0, NA, 0))
  expect_equal(EMDCC(h), 0)
})

test_that("an empty bucket has NA moments and non-finite observations are dropped", {
  x <- c(0.5, 0.5, NA, Inf, 3.5)
  h <- PolynomialHistogram(x, breaks = 0:4, moments = 2)
  expect_equal(h$counts, c(2, 0, 0, 1))
  expect_equal(h$means, c(0.5, NA, NA, 3.5))
  expect_equal(h$variances, c(0, NA, NA, 0))
  expect_equal(Count(h), 3)
})

test_that("PolynomialHistogram validates its arguments", {
  expect_error(PolynomialHistogram("a"))
  expect_error(PolynomialHistogram(1:10, moments = 3))
  expect_error(PolynomialHistogram(1:10, breaks = c(2, 5), moments = 1))
})

# ----- Constructor from aggregates -------------------------------------------

test_that("PreBinnedPolynomialHistogram builds the same object as the raw-data constructor", {
  set.seed(3)
  x <- rgamma(1000, 2)
  br <- c(0, 1, 2, 4, 8, 16)
  h <- PolynomialHistogram(x, br, moments = 2, xname = "x")
  p <- PreBinnedPolynomialHistogram(h$breaks, h$counts, h$means, h$variances, xname = "x")
  expect_identical(p, h)
  p1 <- PreBinnedPolynomialHistogram(h$breaks, h$counts, h$means, xname = "x")
  expect_identical(p1, PolynomialHistogram(x, br, moments = 1, xname = "x"))
})

test_that("PreBinnedPolynomialHistogram validates means and variances", {
  br <- c(0, 1, 2, 4)
  n <- c(5, 0, 3)
  expect_error(PreBinnedPolynomialHistogram(br, n, means = c(1.5, NA, 3)),
               "within its bucket")
  expect_error(PreBinnedPolynomialHistogram(br, n, means = c(NA, NA, 3)),
               "must have a mean")
  expect_error(PreBinnedPolynomialHistogram(br, n, means = c(0.5, 1.5)))
  expect_error(PreBinnedPolynomialHistogram(br, n, means = c(0.5, NA, 3),
                                            variances = c(-1, NA, 0)),
               "non-negative")
  expect_error(PreBinnedPolynomialHistogram(br, n, means = c(0.5, NA, 3),
                                            variances = c(0.1, NA, NA)),
               "must have a variance")
  # The variance cannot exceed (mean - a) (b - mean): here 1 * 1 = 1.
  expect_error(PreBinnedPolynomialHistogram(br, n, means = c(0.5, NA, 3),
                                            variances = c(0.1, NA, 1.01)),
               "at most")
  # Within rounding of that maximum it is clamped to it.
  p <- PreBinnedPolynomialHistogram(br, n, means = c(0.5, NA, 3),
                                    variances = c(0.1, NA, 1 + 1e-12))
  expect_equal(p$variances[3], 1)
  # The mean and variance given for an empty bucket are stored as NA.
  p <- PreBinnedPolynomialHistogram(br, n, means = c(0.5, 1.5, 3),
                                    variances = c(0.1, 0.2, 0.5))
  expect_equal(p$means[2], NA_real_)
  expect_equal(p$variances[2], NA_real_)
  # A mean on a bucket edge is allowed; its variance must then be 0.
  p <- PreBinnedPolynomialHistogram(br, n, means = c(1, NA, 2), variances = c(0, NA, 0))
  expect_equal(p$variances, c(0, NA, 0))
  expect_error(PreBinnedPolynomialHistogram(br, n, means = c(1, NA, 2),
                                            variances = c(0.01, NA, 0)),
               "at most")
})

test_that("a polynomial histogram works with the count-only functions", {
  set.seed(4)
  x <- rexp(300)
  h <- PolynomialHistogram(x, breaks = seq(0, 10, by = 0.5), moments = 2)
  h0 <- hist(x, breaks = seq(0, 10, by = 0.5), plot = FALSE)
  expect_equal(Count(h), length(x))
  expect_equal(ApproxMean(h), mean(x))              # exact from the bucket means
  expect_equal(ApproxMean(h0), weighted.mean(h0$mids, h0$counts))
  expect_equal(KSDCC(h), KSDCC(h0))
  expect_equal(HistToEcdf(h)(h$breaks), HistToEcdf(h0)(h0$breaks))
  expect_no_error(with_null_pdf(plot(h)))
  expect_no_error(with_null_pdf(PlotKSDCC(h)))
})

# ----- Per-bucket closed forms -----------------------------------------------

test_that("EMDCCPerBucket with no moments is 1 and with one moment is the entropy curve", {
  expect_equal(EMDCCPerBucket(), 1)
  a <- seq(0, 1, by = 0.05)
  e1 <- EMDCCPerBucket(a)
  expect_length(e1, length(a))
  expect_equal(e1, rev(e1))                         # symmetric about 1/2
  expect_equal(e1[c(1, length(a))], c(0, 0))        # a mean at an edge is a point mass
  expect_equal(max(e1), log(2))                     # largest at the center
  expect_equal(EMDCCPerBucket(0.5), log(2))
  expect_true(all(e1 >= 0 & e1 <= 1))
  expect_equal(EMDCCPerBucket(0.3), -0.3 * log(0.3) - 0.7 * log(0.7))
  expect_equal(EMDCCPerBucket(c(0.3, NA)), c(EMDCCPerBucket(0.3), NA))
  expect_error(EMDCCPerBucket(1.5), "0, 1")
})

test_that("the one-moment bound crosses the matched-memory count-only bound near alpha = 0.2", {
  # Observation 1 of the paper: a stored mean beats doubling the buckets
  # (per-bucket bound 1/2) only when the mean lies in the outer fifth.
  crossing <- uniroot(function(a) EMDCCPerBucket(a) - 1 / 2, c(0.1, 0.3), tol = 1e-10)$root
  expect_equal(round(crossing, 3), 0.200)
})

test_that("EMDCCPerBucket with two moments matches a numerical integral of the envelope", {
  gap <- function(x, a, t) {
    # Fmax - Fmin at x, Theorems 3.2 and 3.3.
    c1 <- a - t / (1 - a)
    c2 <- a + t / a
    ifelse(x < c1 | x > c2, t / (t + (x - a)^2), (a * (1 - a) - t) / (x * (1 - x)))
  }
  numerical <- function(a, t) {
    k <- c(0, a - t / (1 - a), a + t / a, 1)
    sum(vapply(1:3, function(i) {
      integrate(gap, k[i], k[i + 1], a = a, t = t, rel.tol = 1e-10)$value
    }, numeric(1)))
  }
  for (a in c(0.05, 0.3, 0.5, 0.8)) {
    for (f in c(0.01, 0.3, 0.6, 0.9, 0.99)) {
      t <- f * a * (1 - a)
      expect_equal(EMDCCPerBucket(a, t), numerical(a, t), tolerance = 1e-8)
    }
  }
})

test_that("EMDCCPerBucket with two moments reproduces the Figure 5 labels and is vectorized", {
  # Figure 5 of the paper: at mu = 0.5 and sigma^2 = 0.025, 0.075, 0.15,
  # 0.225 the per-bucket EMDCC is about 0.39, 0.53, 0.47 and 0.20; only
  # the last lies below the H(3K, 0) baseline of 1/3.
  mu <- rep(0.5, 4)
  s2 <- c(0.025, 0.075, 0.15, 0.225)
  e <- EMDCCPerBucket(mu, s2)
  expect_equal(round(e, 3), c(0.393, 0.528, 0.473, 0.197), tolerance = 1e-3)
  expect_true(all(e[1:3] > 1 / 3) && e[4] < 1 / 3)
  expect_equal(e, vapply(1:4, function(i) EMDCCPerBucket(mu[i], s2[i]), numeric(1)))
})

test_that("EMDCCPerBucket with two moments is 0 at the determined configurations and never above one moment", {
  a <- c(0.1, 0.3, 0.5, 0.7, 0.95)
  expect_equal(EMDCCPerBucket(a, rep(0, 5)), rep(0, 5))         # point mass at the mean
  expect_equal(EMDCCPerBucket(a, a * (1 - a)), rep(0, 5))       # atoms at both edges
  expect_equal(EMDCCPerBucket(c(0, 1), c(0, 0)), c(0, 0))       # mean at an edge
  for (f in c(0.001, 0.1, 0.5, 0.9, 0.999)) {
    e2 <- EMDCCPerBucket(a, f * a * (1 - a))
    expect_true(all(e2 > 0 & e2 < EMDCCPerBucket(a)))
  }
  expect_equal(EMDCCPerBucket(c(0.3, NA, 0.5), c(0.1, 0.1, NA)),
               c(EMDCCPerBucket(0.3, 0.1), NA, NA))
  expect_error(EMDCCPerBucket(0.3, 0.3), "alpha")
  expect_error(EMDCCPerBucket(0.3, -0.1), "alpha")
  expect_error(EMDCCPerBucket(c(0.3, 0.4), 0.1))
})

# ----- EMDCC generic ----------------------------------------------------------

test_that("EMDCC on a count-only histogram is unchanged by the generic", {
  set.seed(0)
  h <- hist(rexp(100), plot = FALSE)
  expect_equal(EMDCC(h), sum(diff(h$breaks) * h$counts) / sum(h$counts) / diff(range(h$breaks)))
  expect_equal(EMDCC(h), EMDCC.histogram(h))
  expect_error(EMDCC(1:3))
})

test_that("EMDCC of a polynomial histogram is the count-weighted sum of the per-bucket bounds", {
  set.seed(5)
  x <- rlnorm(5000)
  br <- exp(seq(log(min(x)), log(max(x)), length.out = 13))
  br[1] <- min(x)
  br[13] <- max(x)
  h0 <- hist(x, breaks = br, plot = FALSE)
  h1 <- PolynomialHistogram(x, br, moments = 1)
  h2 <- PolynomialHistogram(x, br, moments = 2)
  w <- diff(br)
  alpha <- (h1$means - head(br, -1)) / w
  tau <- h2$variances / w^2
  expect_equal(EMDCC(h1),
               sum(w * h1$counts * EMDCCPerBucket(alpha)) / sum(h1$counts) / diff(range(br)))
  expect_equal(EMDCC(h2),
               sum(w * h2$counts * EMDCCPerBucket(alpha, tau)) / sum(h2$counts) / diff(range(br)))
  # More stored moments on the same grid never widen the envelope.
  expect_lt(EMDCC(h1), EMDCC(h0))
  expect_lt(EMDCC(h2), EMDCC(h1))
  expect_true(EMDCC(h2) > 0)
})

test_that("EMDCC of a polynomial histogram is exact where the data are quantization-aligned", {
  # Every bucket holds a single value in its interior: the variance is
  # zero, so two moments determine the CDF exactly.  One moment does
  # not (a mean at the bucket center is the worst case for it) unless
  # the value sits on a bucket edge.
  x <- rep(c(1, 2, 3), c(5, 3, 2))
  br <- c(0.5, 1.5, 2.5, 3.5)
  expect_gt(EMDCC(hist(x, breaks = br, plot = FALSE)), 0)
  expect_equal(EMDCC(PolynomialHistogram(x, br, moments = 1)), log(2) / 3)
  expect_equal(EMDCC(PolynomialHistogram(x, br, moments = 2)), 0)
  expect_equal(EMDCC(PolynomialHistogram(x, breaks = 0:3, moments = 1)), 0)
  # Every bucket holds atoms at both of its edges (as an aggregate, since
  # a single bucket convention cannot put both edges in one bucket): the
  # variance is at its maximum, so two moments determine the CDF and
  # one does not.
  br <- c(0, 1, 2, 3)
  means <- c(0.25, 1.5, 2.9)
  vmax <- (means - head(br, -1)) * (tail(br, -1) - means)
  h1 <- PreBinnedPolynomialHistogram(br, c(4, 6, 10), means)
  h2 <- PreBinnedPolynomialHistogram(br, c(4, 6, 10), means, vmax)
  expect_gt(EMDCC(h1), 0)
  expect_equal(EMDCC(h2), 0)
})

test_that("EMDCC ignores empty buckets", {
  h <- PreBinnedPolynomialHistogram(c(0, 1, 2, 3), c(4, 0, 6), c(0.5, NA, 2.5),
                                    c(0.1, NA, 0.2))
  h.dense <- PreBinnedPolynomialHistogram(c(0, 1, 2, 3), c(4, 0, 6), c(0.5, 1.5, 2.5),
                                          c(0.1, 0.1, 0.2))
  expect_equal(EMDCC(h), EMDCC(h.dense))
  expect_equal(EMDCC(h), (4 * EMDCCPerBucket(0.5, 0.1) + 6 * EMDCCPerBucket(0.5, 0.2)) / 10 / 3)
})
