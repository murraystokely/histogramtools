test_that("ReadHistogramsFromDtraceOutputFile parses 54 histograms from the FreeBSD fixture", {
  filename <- system.file(
    "extdata/buildkernel-readsize-dtrace.txt",
    package = "HistogramTools"
  )
  dtrace.hists <- ReadHistogramsFromDtraceOutputFile(filename)
  expect_equal(length(dtrace.hists), 54)
  expect_true(all(sapply(dtrace.hists, inherits, "histogram")))
})
