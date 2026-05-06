test_that("SubsetHistogram trims to the requested range", {
  hist.1 <- hist(c(1, 1, 2, 2, 7), breaks = 0:9, plot = FALSE)
  hist.subset <- SubsetHistogram(hist.1, maxbreak = 3)
  expect_equal(hist.subset$breaks, c(0, 1, 2, 3))
  expect_equal(hist.subset$counts, c(2, 2, 0))
})

test_that("SubsetHistogram errors on out-of-range or non-existent breakpoints", {
  hist.1 <- hist(c(1, 1, 2, 2, 7), breaks = 0:9, plot = FALSE)

  # Out of range
  expect_error(SubsetHistogram(hist.1, minbreak = -100))
  expect_error(SubsetHistogram(hist.1, maxbreak = 100))

  # In range but not an existing breakpoint
  expect_error(SubsetHistogram(hist.1, minbreak = 0.5))
  expect_error(SubsetHistogram(hist.1, maxbreak = 6.5))
})

test_that("SubsetHistogram is a no-op when bounds match the existing endpoints", {
  hist.1 <- hist(c(1, 1, 2, 2, 7), breaks = 0:9, plot = FALSE)

  hist.nosubset <- SubsetHistogram(hist.1, minbreak = 0)
  expect_equal(hist.nosubset$breaks, hist.1$breaks)
  expect_equal(hist.nosubset$counts, hist.1$counts)

  hist.nosubset <- SubsetHistogram(hist.1, maxbreak = 9)
  expect_equal(hist.nosubset$breaks, hist.1$breaks)
  expect_equal(hist.nosubset$counts, hist.1$counts)
})

test_that("SubsetHistogram with both bounds NULL warns and returns the input", {
  hist.1 <- hist(c(1, 1, 2, 2, 7), breaks = 0:9, plot = FALSE)
  expect_warning(
    out <- SubsetHistogram(hist.1),
    "No new breakpoints specified"
  )
  expect_equal(out, hist.1)
})
