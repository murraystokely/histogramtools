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

# ----- Count -----------------------------------------------------------------

test_that("Count returns the sum of counts", {
  h <- PreBinnedHistogram(breaks = 0:4, counts = c(7, 3, 0, 2))
  expect_equal(Count(h), 12)
})

test_that("Count errors on non-histogram input", {
  expect_error(Count(c(1, 2, 3)))
})

# ----- ScaleHistogram --------------------------------------------------------

test_that("ScaleHistogram with explicit factor scales counts linearly", {
  h <- PreBinnedHistogram(breaks = 0:3, counts = c(1, 2, 3))
  out <- ScaleHistogram(h, factor = 10)
  expect_equal(out$counts, c(10, 20, 30))
  expect_equal(out$breaks, h$breaks)
})

test_that("ScaleHistogram default factor (1/Count) normalizes counts to sum to 1", {
  h <- PreBinnedHistogram(breaks = 0:3, counts = c(1, 2, 3))
  out <- ScaleHistogram(h)
  expect_equal(sum(out$counts), 1)
})

test_that("ScaleHistogram errors on non-histogram input or non-scalar factor", {
  h <- PreBinnedHistogram(breaks = 0:3, counts = c(1, 2, 3))
  expect_error(ScaleHistogram(c(1, 2, 3), factor = 2))
  expect_error(ScaleHistogram(h, factor = c(1, 2)))
  expect_error(ScaleHistogram(h, factor = "two"))
})

# ----- PreBinnedHistogram ----------------------------------------------------

test_that("PreBinnedHistogram constructs a valid histogram object", {
  h <- PreBinnedHistogram(breaks = 0:5, counts = c(1, 2, 3, 4, 5), xname = "demo")
  expect_s3_class(h, "histogram")
  expect_equal(h$breaks, 0:5)
  expect_equal(h$counts, c(1, 2, 3, 4, 5))
  expect_equal(h$mids, c(0.5, 1.5, 2.5, 3.5, 4.5))
  expect_true(h$equidist)
  expect_equal(h$xname, "demo")
})

test_that("PreBinnedHistogram detects non-equidist breaks", {
  h <- PreBinnedHistogram(breaks = c(0, 1, 3, 7), counts = c(2, 5, 1))
  expect_false(h$equidist)
})

# ----- AddHistograms edge cases ----------------------------------------------

test_that("AddHistograms with a single histogram returns it (length-1 path)", {
  h <- PreBinnedHistogram(breaks = 0:3, counts = c(1, 2, 3))
  out <- AddHistograms(h)
  # The length(x) == 1 branch returns the list itself.
  expect_equal(out, list(h))
})

test_that("AddHistograms errors when histograms have mismatched breakpoints", {
  h1 <- PreBinnedHistogram(breaks = 0:3, counts = c(1, 1, 1))
  h2 <- PreBinnedHistogram(breaks = 1:4, counts = c(1, 1, 1))
  expect_error(AddHistograms(h1, h2), "identical breakpoints")
})

# ----- Protocol Buffer round-trip --------------------------------------------

test_that("histogram <-> protobuf round-trips identically when RProtoBuf is available", {
  skip_if_not_installed("RProtoBuf")

  set.seed(0)
  h <- hist(rexp(100), plot = FALSE)
  msg <- as.Message(h)
  expect_s4_class(msg, "Message")

  back <- as.histogram(msg)
  expect_s3_class(back, "histogram")
  expect_equal(back$breaks, h$breaks)
  expect_equal(back$counts, h$counts)
  expect_equal(back$mids, h$mids)
})

test_that("as.histogram on a non-HistogramState message errors with a clear message", {
  skip_if_not_installed("RProtoBuf")

  # Construct a message of an unrelated descriptor type. RProtoBuf ships
  # tutorial.proto with descriptor `tutorial.Person`, which we use here
  # as a foreign type.
  if (!exists("tutorial.Person",
              where = "RProtoBuf:DescriptorPool",
              inherits = FALSE)) {
    skip("tutorial.Person descriptor not in RProtoBuf descriptor pool")
  }
  foreign <- RProtoBuf::new(RProtoBuf::P("tutorial.Person"), id = 1, name = "x")
  expect_error(as.histogram(foreign), "Unknown protocol message type")
})
