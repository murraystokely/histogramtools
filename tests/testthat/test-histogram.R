test_that("AddHistograms sums two histograms with identical breaks", {
  hist.1 <- hist(c(1, 2, 3, 4), plot = FALSE)
  hist.2 <- hist(c(1, 2, 2, 4), plot = FALSE)
  hist.sum <- AddHistograms(hist.1, hist.2)
  expect_equal(hist.sum$breaks, c(1, 2, 3, 4))
  expect_equal(hist.sum$counts, c(5, 1, 2))
})

test_that("AddHistograms sums many histograms in either calling form", {
  hist.1 <- hist(c(1, 2, 3), breaks = 0:9, plot = FALSE)
  hist.2 <- hist(c(1, 2, 3), breaks = 0:9, plot = FALSE)
  hist.3 <- hist(c(4, 5, 6), breaks = 0:9, plot = FALSE)

  hist.sum <- AddHistograms(hist.1, hist.2, hist.3)
  expect_equal(hist.sum$breaks, 0:9)
  expect_equal(hist.sum$counts, c(2, 2, 2, 1, 1, 1, 0, 0, 0))

  hist.sum <- AddHistograms(x = list(hist.1, hist.2, hist.3))
  expect_equal(hist.sum$breaks, 0:9)
  expect_equal(hist.sum$counts, c(2, 2, 2, 1, 1, 1, 0, 0, 0))
})

test_that("MergeBuckets coarsens by adjacency, by count, and by explicit breaks", {
  hist.1 <- hist(c(1, 2, 3), breaks = 0:9, plot = FALSE)

  hist.2 <- MergeBuckets(hist.1, adj = 2)
  expect_equal(hist.2$breaks, c(0, 2, 4, 6, 8, 9))
  expect_equal(hist.2$counts, c(2, 1, 0, 0, 0))

  hist.3 <- MergeBuckets(hist.1, breaks = 3)
  expect_equal(hist.3$breaks, c(0, 3, 6, 9))
  expect_equal(hist.3$counts, c(3, 0, 0))

  # Breakpoint list that is not a subset of the existing list errors out.
  expect_error(MergeBuckets(hist.1, breaks = c(17, 34)))

  hist.4 <- MergeBuckets(hist.1, breaks = c(0, 3, 6, 9))
  expect_equal(hist.4$breaks, c(0, 3, 6, 9))
  expect_equal(hist.4$counts, c(3, 0, 0))

  # A breakpoint range that would drop existing buckets errors out.
  expect_error(MergeBuckets(hist.1, breaks = 2:9))

  hist.5 <- hist(c(1, 2, 3), breaks = 0:10, plot = FALSE)
  hist.6 <- MergeBuckets(hist.5, adj = 2)
  expect_equal(hist.6$breaks, c(0, 2, 4, 6, 8, 10))
  expect_equal(hist.6$counts, c(2, 1, 0, 0, 0))
})
