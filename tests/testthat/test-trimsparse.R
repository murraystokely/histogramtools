test_that("TrimHistogram drops trailing/leading empty buckets but preserves total count", {
  hist.1 <- hist(c(1, 2, 3), breaks = 0:9, plot = FALSE)
  hist.trimmed <- TrimHistogram(hist.1)
  expect_equal(sum(hist.1$counts), sum(hist.trimmed$counts))
  expect_true(length(hist.1$counts) > length(hist.trimmed$counts))

  hist.2 <- hist(c(4, 5, 6), breaks = 0:9, plot = FALSE)
  hist.trimmed <- TrimHistogram(hist.2)
  expect_equal(sum(hist.2$counts), sum(hist.trimmed$counts))
  expect_true(length(hist.2$counts) > length(hist.trimmed$counts))
})

test_that("TrimHistogram is a no-op when every bucket is empty", {
  zero.hist <- suppressWarnings(
    hist(numeric(), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), plot = FALSE)
  )
  zero.trimmed <- TrimHistogram(zero.hist)
  expect_equal(length(zero.hist$counts), length(zero.trimmed$counts))
})
