## HistogramTools

<!-- badges: start -->

[![R-CMD-check](https://github.com/murraystokely/histogramtools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/murraystokely/histogramtools/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/HistogramTools)](https://cran.r-project.org/package=HistogramTools)
[![Codecov test coverage](https://codecov.io/gh/murraystokely/histogramtools/graph/badge.svg)](https://app.codecov.io/gh/murraystokely/histogramtools)
<!-- badges: end -->

This package provides a number of utility functions useful for manipulating large histograms. This includes methods to trim, subset, merge buckets, merge histograms, convert to CDF, and calculate information loss due to binning. It also provides a protocol buffer representation of R's native histogram class to allow histograms over large data sets to be computed and combined in distributed analytical pipelines.

## Installation

You can either install from source via this repo, or install
[the CRAN package](https://cran.r-project.org/package=HistogramTools)
the usual way from [R](https://www.r-project.org).

## More Info

[RProtoBuf & HistogramTools: Statistical Analysis Tools for Large Data Sets](https://opensource.googleblog.com/2013/10/rprotobuf-histogramtools-statistical_10.html) Google Open Source Blog, October 10, 2013

## History

This package was originally developed at Google between 2011 and 2015.
It is now independently maintained by the original author.

## Authors

Murray Stokely

## License

Apache 2.0


