# Tests for HistToASH (R/ash.R), a thin wrapper around ash::ash1.
# We don't re-test the ash package itself, only that our wrapper hands
# off correctly and validates its inputs.

test_that("HistToASH returns an ash-shaped result on an equidist histogram", {
  set.seed(0)
  h <- hist(rexp(1000), plot = FALSE)
  out <- HistToASH(h)
  # ash1 returns a list with at least 'x' and 'y' components.
  expect_true(is.list(out))
  expect_true(all(c("x", "y") %in% names(out)))
  expect_equal(length(out$x), length(out$y))
})

test_that("HistToASH errors on non-histogram input", {
  expect_error(HistToASH(c(1, 2, 3)))
  expect_error(HistToASH(list()))
})

test_that("HistToASH errors on a non-equidist histogram", {
  # PreBinnedHistogram with deliberately uneven breaks.
  h <- PreBinnedHistogram(breaks = c(0, 1, 3, 7), counts = c(2, 5, 1))
  expect_false(h$equidist)
  expect_error(HistToASH(h), "equidist")
})
