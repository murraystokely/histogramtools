# Tests that the histogram manipulation functions keep a polynomial
# histogram's per-bucket moments consistent: AddHistograms pools them
# exactly, MergeBuckets pools adjacent buckets, SubsetHistogram and
# TrimHistogram carry them along, ScaleHistogram leaves them alone, and
# IntersectHistograms returns a count-only histogram.

expect_same_moments <- function(h, ref) {
  expect_s3_class(h, c("polynomial.histogram", "histogram"), exact = TRUE)
  expect_equal(h$breaks, ref$breaks)
  expect_equal(h$counts, ref$counts)
  expect_equal(h$means, ref$means)
  expect_equal(h$variances, ref$variances)
  expect_equal(h$density, ref$density)
  expect_equal(h$mids, ref$mids)
}

set.seed(10)
x <- c(rlnorm(3000), rep(2, 100), rep(7, 50))
br <- c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 64)

# ----- AddHistograms -----------------------------------------------------------

test_that("AddHistograms pools polynomial histograms to the histogram of the pooled data", {
  i <- seq_along(x) %% 3
  for (m in 1:2) {
    parts <- lapply(0:2, function(k) PolynomialHistogram(x[i == k], br, moments = m))
    whole <- PolynomialHistogram(x, br, moments = m)
    merged <- AddHistograms(x = parts)
    expect_same_moments(merged, whole)
    expect_equal(EMDCC(merged), EMDCC(whole))
    expect_equal(ApproxMean(merged), mean(x))
    expect_equal(merged$xname, "Merge of 3 histograms")
    two <- AddHistograms(parts[[1]], parts[[2]])
    expect_equal(two$counts, parts[[1]]$counts + parts[[2]]$counts)
    expect_equal(two$xname, "Merge of x[i == k] and x[i == k]")
  }
})

test_that("AddHistograms pools empty buckets and single-input groups correctly", {
  a <- PreBinnedPolynomialHistogram(0:3, c(5, 0, 2), c(0.4, NA, 2.5), c(0.02, NA, 0.1))
  b <- PreBinnedPolynomialHistogram(0:3, c(0, 3, 2), c(NA, 1.5, 2.9), c(NA, 0.05, 0.01))
  s <- AddHistograms(a, b)
  expect_equal(s$counts, c(5, 3, 4))
  expect_equal(s$means, c(0.4, 1.5, 2.7))
  # Pooled variance: mean within-variance plus spread of the means.
  expect_equal(s$variances[3], (2 * (0.1 + 0.2^2) + 2 * (0.01 + 0.2^2)) / 4)
  expect_equal(s$variances[1:2], c(0.02, 0.05))
  # One-moment inputs give a one-moment result.
  a1 <- PreBinnedPolynomialHistogram(0:3, c(5, 0, 2), c(0.4, NA, 2.5))
  b1 <- PreBinnedPolynomialHistogram(0:3, c(0, 3, 2), c(NA, 1.5, 2.9))
  s1 <- AddHistograms(a1, b1)
  expect_null(s1$variances)
  expect_equal(s1$means, s$means)
})

test_that("AddHistograms refuses to mix moment orders or count-only histograms", {
  h0 <- hist(x, breaks = br, plot = FALSE)
  h1 <- PolynomialHistogram(x, br, moments = 1)
  h2 <- PolynomialHistogram(x, br, moments = 2)
  expect_error(AddHistograms(h0, h1), "same number of per-bucket moments")
  expect_error(AddHistograms(h1, h2), "same number of per-bucket moments")
  expect_error(AddHistograms(h2, PolynomialHistogram(x, c(0, 1, 64), moments = 2)),
               "identical breakpoints")
  # Count-only addition is unchanged.
  s <- AddHistograms(h0, h0)
  expect_s3_class(s, "histogram", exact = TRUE)
  expect_equal(s$counts, 2 * h0$counts)
})

# ----- MergeBuckets ------------------------------------------------------------

test_that("MergeBuckets pools adjacent buckets of a polynomial histogram", {
  for (m in 1:2) {
    h <- PolynomialHistogram(x, br, moments = m)
    # Pairs of adjacent buckets: the same as building on every other break.
    coarse <- PolynomialHistogram(x, br[c(TRUE, FALSE)], moments = m)
    expect_same_moments(MergeBuckets(h, adj.buckets = 2), coarse)
    expect_same_moments(MergeBuckets(h, breaks = br[c(TRUE, FALSE)]), coarse)
    # An uneven grouping by a break list.
    nb <- c(0, 1, 8, 64)
    expect_same_moments(MergeBuckets(h, breaks = nb), PolynomialHistogram(x, nb, moments = m))
    # A number of buckets that does not divide evenly: the last group is short.
    merged <- MergeBuckets(h, breaks = 3)
    expect_equal(length(merged$counts), 3)
    expect_equal(sum(merged$counts), sum(h$counts))
    expect_equal(ApproxMean(merged), mean(x))
    expect_lte(EMDCC(h), EMDCC(merged))
  }
  h2 <- PolynomialHistogram(x, br, moments = 2)
  expect_error(MergeBuckets(h2, adj.buckets = 2, FUN = max), "FUN=sum")
  # Count-only merging is unchanged, including FUN.
  h0 <- hist(x, breaks = br, plot = FALSE)
  expect_equal(MergeBuckets(h0, adj.buckets = 2)$counts,
               as.vector(tapply(h0$counts, rep(1:4, each = 2), sum)))
  expect_equal(MergeBuckets(h0, adj.buckets = 2, FUN = max)$counts,
               as.vector(tapply(h0$counts, rep(1:4, each = 2), max)))
})

# ----- SubsetHistogram and TrimHistogram -----------------------------------------

test_that("SubsetHistogram keeps the moments of the retained buckets", {
  for (m in 1:2) {
    h <- PolynomialHistogram(x, br, moments = m)
    s <- SubsetHistogram(h, minbreak = 0.5, maxbreak = 8)
    keep <- 3:6
    expect_s3_class(s, c("polynomial.histogram", "histogram"), exact = TRUE)
    expect_equal(s$breaks, br[3:7])
    expect_equal(s$counts, h$counts[keep])
    expect_equal(s$means, h$means[keep])
    expect_equal(s$variances, h$variances[keep])
    expect_equal(s$density, s$counts / (sum(s$counts) * diff(s$breaks)))
    expect_equal(SubsetHistogram(h, minbreak = 0.5)$means, h$means[-(1:2)])
    expect_equal(SubsetHistogram(h, maxbreak = 8)$variances, h$variances[1:6])
  }
})

test_that("TrimHistogram keeps the moments of the retained buckets", {
  h <- PreBinnedPolynomialHistogram(0:7, c(0, 0, 4, 6, 0, 2, 0),
                                    c(NA, NA, 2.5, 3.2, NA, 5.5, NA),
                                    c(NA, NA, 0.1, 0.08, NA, 0.2, NA))
  t <- TrimHistogram(h)
  expect_s3_class(t, c("polynomial.histogram", "histogram"), exact = TRUE)
  expect_equal(t$breaks, 2:6)
  expect_equal(t$counts, c(4, 6, 0, 2))
  expect_equal(t$means, c(2.5, 3.2, NA, 5.5))
  expect_equal(t$variances, c(0.1, 0.08, NA, 0.2))
  expect_equal(t$mids, 2.5:5.5)
  expect_equal(TrimHistogram(h, right = FALSE)$means, c(2.5, 3.2, NA, 5.5, NA))
  expect_equal(TrimHistogram(h, left = FALSE)$variances, c(NA, NA, 0.1, 0.08, NA, 0.2))
  expect_equal(EMDCC(t) * diff(range(t$breaks)), EMDCC(h) * diff(range(h$breaks)))
  expect_equal(HistToEcdf(t)(3), HistToEcdf(h)(3))
})

# ----- ScaleHistogram and IntersectHistograms -------------------------------------

test_that("ScaleHistogram leaves the per-bucket moments alone", {
  h <- PolynomialHistogram(x, br, moments = 2)
  s <- ScaleHistogram(h)
  expect_equal(sum(s$counts), 1)
  expect_equal(s$means, h$means)
  expect_equal(s$variances, h$variances)
  expect_equal(EMDCC(s), EMDCC(h))
  expect_equal(HistToEcdf(s)(c(1, 3)), HistToEcdf(h)(c(1, 3)))
})

test_that("IntersectHistograms of polynomial histograms is count-only", {
  a <- PolynomialHistogram(x[1:1500], br, moments = 2)
  b <- PolynomialHistogram(x[1501:3150], br, moments = 2)
  i <- IntersectHistograms(a, b)
  expect_s3_class(i, "histogram", exact = TRUE)
  expect_null(i$means)
  expect_equal(i$counts, pmin(a$counts, b$counts))
})

test_that(".UpdateHistogram catches moments that no longer match the counts", {
  h <- PolynomialHistogram(x, br, moments = 1)
  h$counts <- h$counts[-1]
  h$breaks <- h$breaks[-1]
  expect_error(.UpdateHistogram(h))
})

test_that("converting a polynomial histogram to a protocol buffer warns that moments are dropped", {
  h <- PolynomialHistogram(x, br, moments = 1)
  if (requireNamespace("RProtoBuf", quietly = TRUE)) {
    expect_warning(as.Message(h), "moments are not stored")
  } else {
    expect_warning(try(as.Message(h), silent = TRUE), "moments are not stored")
  }
})
