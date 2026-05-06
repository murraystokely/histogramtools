# Smoke tests for the plot wrappers in R/plot.R. These verify that the
# functions run to completion without erroring on representative input;
# visual correctness is out of scope (use vdiffr if that ever matters).
#
# Each test opens a null PDF device so we don't try to render to the
# active device during R CMD check.

with_null_pdf <- function(expr) {
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file)
  on.exit({
    dev.off()
    unlink(pdf_file)
  })
  force(expr)
}

# ----- PlotRelativeFrequency -------------------------------------------------

test_that("PlotRelativeFrequency runs without error on a simple histogram", {
  set.seed(0)
  h <- hist(runif(100), plot = FALSE)
  expect_no_error(with_null_pdf(PlotRelativeFrequency(h)))
})

# ----- PlotLog2ByteEcdf ------------------------------------------------------

# A histogram whose breaks land on powers of 2 from 2^0 .. 2^10. That gives
# 11 power-of-2 boundaries -- well above the >= 3 the function requires.
.pow2_histogram <- function() {
  breaks <- 2^(0:10)
  counts <- rep(10, length(breaks) - 1)
  PreBinnedHistogram(breaks = breaks, counts = counts)
}

test_that("PlotLog2ByteEcdf runs on a histogram with power-of-2 breaks", {
  expect_no_error(with_null_pdf(PlotLog2ByteEcdf(.pow2_histogram())))
})

test_that("PlotLog2ByteEcdf accepts a pre-computed ecdf", {
  e <- HistToEcdf(.pow2_histogram())
  expect_no_error(with_null_pdf(PlotLog2ByteEcdf(e)))
})

test_that("PlotLog2ByteEcdf errors when the ecdf has fewer than 3 powers-of-2 in its knots", {
  # Breaks deliberately not aligned to powers of two -> 0 power-of-2 knots.
  h <- PreBinnedHistogram(breaks = c(3, 5, 7, 9), counts = c(1, 2, 3))
  expect_error(with_null_pdf(PlotLog2ByteEcdf(h)), "Insufficient powers of 2")
})

# ----- PlotLogTimeDurationEcdf ----------------------------------------------

test_that("PlotLogTimeDurationEcdf runs on a histogram with time-duration-shaped breaks", {
  # Breaks roughly at 1s, 1m, 1h, 1d, 1w, 1mo, 1y in seconds.
  breaks <- c(1, 60, 3600, 86400, 86400 * 7, 86400 * 30, 86400 * 365)
  counts <- rep(5, length(breaks) - 1)
  h <- PreBinnedHistogram(breaks = breaks, counts = counts)
  expect_no_error(with_null_pdf(PlotLogTimeDurationEcdf(h)))
})

test_that("PlotLogTimeDurationEcdf accepts a pre-computed ecdf", {
  breaks <- c(1, 60, 3600, 86400, 86400 * 7, 86400 * 30, 86400 * 365)
  counts <- rep(5, length(breaks) - 1)
  e <- HistToEcdf(PreBinnedHistogram(breaks = breaks, counts = counts))
  expect_no_error(with_null_pdf(PlotLogTimeDurationEcdf(e)))
})
