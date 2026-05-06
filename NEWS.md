# HistogramTools 0.4.1

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

# HistogramTools 0.4.0

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
