## HistogramTools 0.5.0

This release adds polynomial histograms (histograms that store the mean
and variance of each bucket alongside its count), the closed-form
worst-case information-loss bound for them, and CDF and quantile
estimation from the stored moments.  See `NEWS.md`.

* `EMDCC()`, `HistToEcdf()`, `ApproxQuantile()` and `ApproxMean()` are
  now S3 generics.  Their behaviour on ordinary `"histogram"` objects is
  unchanged and covered by the existing tests.
* New dependency in `Suggests`: `statmod`, used by one test to
  regenerate a tabulated quadrature rule.
* Tim Hesterberg is added as author and copyright holder.

## Test environments

* Local: macOS (Apple Silicon), R 4.6.0: `R CMD check --as-cran`
  0 errors, 0 warnings.
* GitHub Actions matrix (R-CMD-check): macOS-latest, Windows-latest,
  Ubuntu-latest x R release / devel / oldrel-1.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies to check.
