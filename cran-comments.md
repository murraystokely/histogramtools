## Resubmission

This is a resubmission of HistogramTools 0.4.1. It addresses the four
points raised by the CRAN reviewer on the 0.4.0 resubmission:

* `DESCRIPTION` now cites the methods implemented by the package
  (Rubner, Tomasi & Guibas 2000; Swain & Ballard 1991; Puzicha,
  Hofmann & Buhmann 1997; Scott 2015) using the `<doi:...>` and
  `ISBN:...` form requested by the cookbook.

* Every exported function listed in the reviewer's note now has a
  `\value` section describing the class and meaning of the returned
  object. Plotting helpers explicitly state "No return value, called
  for side effects".

* `man/dtrace.Rd` no longer needs `\dontrun{}` — the example now
  parses the bundled `inst/extdata/buildkernel-readsize-dtrace.txt`
  sample so it runs without a live DTrace install.
  `man/informationloss.Rd` unwraps the `PlotKSDCC` / `PlotEMDCC`
  calls; they complete in well under five seconds. No `\dontrun{}`
  remains in the package.

* Both Sweave vignettes (`HistogramTools.Rnw`,
  `HistogramTools-quickref.Rnw`) capture `par()` and `options()` at
  the top and restore them at the end, so building the vignettes no
  longer mutates the user's graphics parameters or options.

## Test environments

* Local: Ubuntu 24.04, R 4.3.3 — `R CMD check --as-cran`: 0 errors,
  0 warnings, 3 notes (all environment-only: time-sync,
  RProtoBuf headers not installed locally — `RProtoBuf` is an
  `Enhances:` dependency).
* GitHub Actions matrix (R-CMD-check): macOS-latest, Windows-latest,
  Ubuntu-latest × R release / devel / oldrel-1 — all green.
* win-builder: `devtools::check_win_devel()` submitted 2026-05-06.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Murray Stokely <murray@stokely.org>'
  New submission
  Package was archived on CRAN

  This is the resubmission of an archived package. The 0.4.0
  resubmission addressed the original archival cause (the
  package-alias mismatch). The current 0.4.1 resubmission addresses
  the additional reviewer feedback summarised above.

## Reverse dependencies

There are no reverse dependencies to check.
