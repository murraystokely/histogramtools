# HistogramTools 0.5.0 (2026-11-21)

Adds polynomial histograms: histograms whose buckets store the mean and,
optionally, the variance of their observations alongside the count.

* New S3 class `polynomial.histogram`, inheriting from `histogram`, with
  components `means` and (for two moments) `variances`.  Every function
  that reads the breaks and counts of a histogram works on it unchanged.
* `PolynomialHistogram(x, breaks, moments, right, include.lowest)` builds
  one from raw data on the buckets `hist()` would use, computing the
  moments on exactly the bucket assignment that produced the counts;
  `moments=0` returns `hist()`'s object.  `PreBinnedPolynomialHistogram()`
  builds one from pre-aggregated per-bucket statistics and validates them.
* `EMDCC()` is now a generic.  For a polynomial histogram it returns the
  worst-case information loss from the closed forms for one and two
  stored moments; `EMDCCPerBucket(alpha, tau)` exposes the per-bucket
  forms on the unit bucket.  Count-only histograms are unchanged.
* `HistToEcdf()`, `ApproxQuantile()` and `ApproxMean()` are now generics.
  For a polynomial histogram the CDF and quantiles come from filling each
  bucket with the maximum-entropy density that has its stored moments, so
  the estimate never leaves the EMDCC envelope; buckets whose moments
  determine their contents (an edge mean, zero variance or the maximum
  variance) are reconstructed exactly.  `ApproxMean()` of a polynomial
  histogram is the exact data mean.
* `AddHistograms()` merges polynomial histograms exactly, pooling counts,
  means and variances.  `MergeBuckets()`, `SubsetHistogram()` and
  `TrimHistogram()` keep the per-bucket moments aligned with the
  retained or merged buckets.  `IntersectHistograms()` returns a
  count-only histogram.  `as.Message()` warns that the protocol buffer
  representation does not carry the moments.
* `AddHistograms()` with a single histogram now returns that histogram;
  through 0.4.1 it returned the list wrapping it.
* `.UpdateHistogram()` stops if a histogram's per-bucket moments are not
  the same length as its counts.
* New vignette section "Polynomial Histograms".
* Tim Hesterberg joins as author.

# HistogramTools 0.4.1 (2026-05-06)

Addresses CRAN reviewer feedback on the 0.4.0 resubmission.

* `DESCRIPTION`: cite the methods implemented by the package
  (Rubner, Tomasi & Guibas 2000; Swain & Ballard 1991; Puzicha, Hofmann &
  Buhmann 1997; Scott 2015) using the `<doi:...>` / `ISBN:...` form CRAN
  expects.
* `man/*.Rd`: add a `\value` section to every exported function, giving
  the class and meaning of the returned object (or noting that the
  function is called for its side effect).
* `man/dtrace.Rd`: replace the `\dontrun{}` example, which required a
  live DTrace install, with a runnable one that parses the bundled
  `inst/extdata/buildkernel-readsize-dtrace.txt` sample.
* `man/informationloss.Rd`: unwrap the `PlotKSDCC` / `PlotEMDCC` calls
  from `\dontrun{}`; they run in well under five seconds.
* `vignettes/HistogramTools.Rnw`, `vignettes/HistogramTools-quickref.Rnw`:
  capture `par()` and `options()` at the top of each vignette and
  restore them at the end, so building the vignettes leaves the user's
  graphics parameters and options unchanged.

# HistogramTools 0.4.0 (2026-04-29)

Returns the package to CRAN after its 2024-04-20 archival.

* Fixes the package alias to match the package name (`HistogramTools-package`),
  the NOTE that triggered the archival.
* Switches `Suggests`/`Enhances` examples from `if (require(RProtoBuf))` to
  `if (requireNamespace("RProtoBuf", quietly = TRUE))` per CRAN policy 3.1.4,
  with `RProtoBuf::P()` fully qualified so example bodies don't depend on
  attachment.
* Declares `Encoding: UTF-8`.
* Replaces the dead R-Forge URL with the GitHub repository
  (<https://github.com/murraystokely/histogramtools>) and adds `BugReports`.
* Drops `stringr` from `Imports`; the single `str_trim` call is replaced by
  base R `trimws()`.
* `inst/NEWS.Rd` and `ChangeLog` are replaced by this `NEWS.md`. The detailed
  per-commit history remains in `git log`.

# HistogramTools 0.3.2 (2015-07-29)

* `NAMESPACE`: add explicit `importFrom` for symbols in `graphics`, `stats`,
  and `utils` that R 3.3+ requires to be declared.
* `R/histogram.R`: replace `stopifnot(require(RProtoBuf))` with
  `requireNamespace()` and an explicit error; use `RProtoBuf::new()` for
  message construction since `new` is generic in RProtoBuf.
* `DESCRIPTION`: drop the leading "This package" from the `Description`
  field per CRAN policy.

# HistogramTools 0.3.1 (2014-08-26)

* `DESCRIPTION`: add `Copyright` and `Authors@R` fields.
* `R/plot.R`: fully qualify `gdata::humanReadable` so `R CMD check` on R-devel
  is clean when the conditionally-suggested `gdata` package is available.

# HistogramTools 0.3 (2013-12-11)

* Move `Hmisc` from `Depends` to `Imports`.
* Significantly improved introduction vignette.
* Added `ScaleHistograms()`.
* Added `PlotRelativeFrequency()` for plotting relative-frequency histograms.
* Added two-histogram distance and divergence measures: `minkowski.dist`,
  `intersect.dist`, `kl.divergence`, `jeffrey.divergence`.
* Added `PreBinnedHistogram()` for constructing histogram objects directly
  from a vector of breaks and counts.

# HistogramTools 0.2 (2013-10-02)

* Moved `RProtoBuf` from `Suggests` to `Enhances` and made all uses
  conditional, so the package can be built and checked without it.
* Improved plot titles for `PlotKSDCC` and `PlotEMDCC`.
* Added a quick-reference vignette.
* Implemented histogram smoothing via Average Shifted Histograms from the
  `ash` package.
* Added `PlotLog2ByteEcdf()` and `PlotLogTimeDurationEcdf()` for ECDFs
  of data binned on power-of-two byte or log-scale time-duration boundaries.
* Removed the unused `methods` dependency and tightened up `NAMESPACE`
  imports/exports.

# HistogramTools 0.1 (2013-09-13)

* Initial CRAN release.
