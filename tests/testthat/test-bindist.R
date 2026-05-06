# Tests for the two-histogram comparison functions in R/bindist.R:
# IntersectHistograms, minkowski.dist, intersect.dist, kl.divergence,
# jeffrey.divergence.

# A handful of helpers used across the tests below. Building histograms
# with hand-picked counts keeps the expected values easy to verify by
# hand from the formulas.
hist_with_counts <- function(counts, breaks = seq(0, length(counts))) {
  PreBinnedHistogram(breaks = breaks, counts = counts)
}

# ----- IntersectHistograms ---------------------------------------------------

test_that("IntersectHistograms takes elementwise pmin of counts", {
  h1 <- hist_with_counts(c(3, 2, 1))
  h2 <- hist_with_counts(c(1, 2, 3))
  out <- IntersectHistograms(h1, h2)
  expect_s3_class(out, "histogram")
  expect_equal(out$breaks, h1$breaks)
  expect_equal(out$counts, c(1, 2, 1))
})

test_that("IntersectHistograms errors on mismatched breakpoints", {
  h1 <- hist_with_counts(c(1, 1, 1), breaks = 0:3)
  h2 <- hist_with_counts(c(1, 1, 1), breaks = 1:4)
  expect_error(IntersectHistograms(h1, h2))
})

test_that("IntersectHistograms errors when an input is not a histogram", {
  h <- hist_with_counts(c(1, 2, 3))
  expect_error(IntersectHistograms(h, list()))
  expect_error(IntersectHistograms(c(1, 2, 3), h))
})

# ----- minkowski.dist --------------------------------------------------------

test_that("minkowski.dist returns 0 when histograms are identical", {
  h <- hist_with_counts(c(1, 2, 3))
  expect_equal(minkowski.dist(h, h, p = 1), 0)
  expect_equal(minkowski.dist(h, h, p = 2), 0)
})

test_that("minkowski.dist computes L1 distance with p=1", {
  # diff in counts: (4, -1, -1); sum of |diff| = 6
  h1 <- hist_with_counts(c(5, 0, 0))
  h2 <- hist_with_counts(c(1, 1, 1))
  expect_equal(minkowski.dist(h1, h2, p = 1), 6)
})

test_that("minkowski.dist computes L2 distance with p=2", {
  # diff in counts: (3, -4); sum of squares = 25; sqrt = 5
  h1 <- hist_with_counts(c(3, 0))
  h2 <- hist_with_counts(c(0, 4))
  expect_equal(minkowski.dist(h1, h2, p = 2), 5)
})

test_that("minkowski.dist errors on non-positive or non-scalar p", {
  h <- hist_with_counts(c(1, 2, 3))
  expect_error(minkowski.dist(h, h, p = 0))
  expect_error(minkowski.dist(h, h, p = -1))
  expect_error(minkowski.dist(h, h, p = c(1, 2)))
})

test_that("minkowski.dist errors on mismatched breakpoints", {
  h1 <- hist_with_counts(c(1, 1, 1), breaks = 0:3)
  h2 <- hist_with_counts(c(1, 1, 1), breaks = 1:4)
  expect_error(minkowski.dist(h1, h2, p = 1))
})

# ----- intersect.dist --------------------------------------------------------

test_that("intersect.dist is 0 when histograms are identical", {
  h <- hist_with_counts(c(1, 2, 3))
  expect_equal(intersect.dist(h, h), 0)
})

test_that("intersect.dist is 1 when histograms are completely disjoint", {
  # pmin on disjoint supports is all zeros => intersect / sum(h2) = 0
  h1 <- hist_with_counts(c(5, 0))
  h2 <- hist_with_counts(c(0, 5))
  expect_equal(intersect.dist(h1, h2), 1)
})

test_that("intersect.dist matches a hand-computed partial overlap", {
  # pmin((5,5), (3,7)) = (3,5); sum = 8; sum(h2) = 10; 1 - 8/10 = 0.2
  h1 <- hist_with_counts(c(5, 5))
  h2 <- hist_with_counts(c(3, 7))
  expect_equal(intersect.dist(h1, h2), 0.2)
})

# ----- kl.divergence --------------------------------------------------------

test_that("kl.divergence is 0 when histograms are identical (no zero counts)", {
  h <- hist_with_counts(c(1, 2, 3))
  expect_equal(kl.divergence(h, h), 0)
})

test_that("kl.divergence matches a hand-computed two-bin example", {
  # h1 normalized = (0.75, 0.25); h2 normalized = (0.25, 0.75)
  # KL = 0.75 * log(0.75/0.25) + 0.25 * log(0.25/0.75) = 0.5 * log(3)
  h1 <- hist_with_counts(c(3, 1))
  h2 <- hist_with_counts(c(1, 3))
  expect_equal(kl.divergence(h1, h2), 0.5 * log(3))
})

test_that("kl.divergence is asymmetric in general (KL(h1,h2) != KL(h2,h1))", {
  # Note: swap-symmetric distributions like (0.8, 0.2) and (0.2, 0.8) happen
  # to give the same KL in either direction, so we deliberately use a pair
  # that isn't swap-symmetric.
  h1 <- hist_with_counts(c(3, 1))   # normalized (0.75, 0.25)
  h2 <- hist_with_counts(c(2, 2))   # normalized (0.5, 0.5)
  expect_false(isTRUE(all.equal(
    kl.divergence(h1, h2),
    kl.divergence(h2, h1)
  )))
})

# ----- jeffrey.divergence ----------------------------------------------------

test_that("jeffrey.divergence is 0 when histograms are identical", {
  h <- hist_with_counts(c(1, 2, 3))
  expect_equal(jeffrey.divergence(h, h), 0)
})

test_that("jeffrey.divergence is symmetric: J(h1,h2) == J(h2,h1)", {
  h1 <- hist_with_counts(c(3, 1))
  h2 <- hist_with_counts(c(1, 3))
  expect_equal(jeffrey.divergence(h1, h2), jeffrey.divergence(h2, h1))
})

test_that("jeffrey.divergence matches a hand-computed two-bin example", {
  # h1 normalized = (0.75, 0.25); h2 normalized = (0.25, 0.75)
  # m = (0.5, 0.5)
  # J = 0.75 * log(0.75/0.5) + 0.25 * log(0.25/0.5)         # h1 contribution
  #   + 0.25 * log(0.25/0.5) + 0.75 * log(0.75/0.5)         # h2 contribution
  #   = 1.5 * log(1.5) + 0.5 * log(0.5)
  h1 <- hist_with_counts(c(3, 1))
  h2 <- hist_with_counts(c(1, 3))
  expect_equal(jeffrey.divergence(h1, h2), 1.5 * log(1.5) + 0.5 * log(0.5))
})
