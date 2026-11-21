% Copyright 2013 Google Inc. All Rights Reserved.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Author: mstokely@google.com (Murray Stokely)

\name{HistToEcdf}
\alias{HistToEcdf}
\alias{HistToEcdf.histogram}
\alias{HistToEcdf.polynomial.histogram}
\title{Empirical Cumulative Distribution Function From a Histogram.}
\description{
Computes an approximate empirical cumulative distribution function of a
data set given a binned histogram representation of that dataset, using
the per-bucket moments when the histogram stores them.
}
\usage{
HistToEcdf(h, \dots)
\method{HistToEcdf}{histogram}(h, method="constant", f=0, inverse=FALSE, \dots)
\method{HistToEcdf}{polynomial.histogram}(h, inverse=FALSE, \dots)
}
\arguments{
  \item{h}{A histogram object (created by \code{\link[graphics]{hist}})
    representing a pre-binned dataset on which we'd like to calculate
    an ECDF, or a polynomial histogram (created by
    \code{\link{PolynomialHistogram}}) that also stores per-bucket
    moments.}
  \item{method}{specifies the interpolation method to be used in call to
    \code{\link[stats]{approxfun}()}.  Choices are \code{"linear"} or \code{"constant"}.}
  \item{f}{if \code{method="constant"} a number between 0 and 1 inclusive,
    indicating a compromise between left- and right-continuous
    step functions.  See \code{?approxfun}}
  \item{inverse}{if \code{TRUE} return the inverse Ecdf.}
  \item{\dots}{Unused.}
}
\details{
  This function approximates the e.c.d.f. (empirical cumulative
  distribution function) of a data set given a binned
  histogram representation of that data set.

  For an ordinary histogram the e.c.d.f. passes through the cumulative
  counts at the breaks and is constant (or, with
  \code{method="linear"}, linear) in between.

  For a polynomial histogram each bucket is filled with the
  maximum-entropy density that has the bucket's stored mean and, when
  present, variance, so the estimate assumes nothing beyond the stored
  moments.  A bucket whose stored moments determine its contents
  exactly (a mean on a bucket edge or a variance of zero, which is a
  point mass, or the largest feasible variance, which is two atoms on
  the edges) is reconstructed as those atoms.  Because the
  reconstruction has the stored moments, it lies within the envelope of
  CDFs consistent with the histogram, so its \eqn{L_1} distance from
  the true CDF never exceeds \code{\link{EMDCC}(h)}.  With one stored
  moment the per-bucket fit and its inverse are in closed form; with
  two, the two exponents of \eqn{\exp(\lambda_1 u + \lambda_2 u^2)} are
  solved by damped Newton iteration and the inverse by Brent's method.
}
\value{
  For an ordinary histogram, a function of class
  \code{c("ecdf", "stepfun", "function")}, as returned by
  \code{\link[stats]{ecdf}}.  Calling it with a numeric vector \code{x}
  returns the approximate cumulative probability at each value of
  \code{x}.  If \code{inverse=TRUE} the returned function instead maps
  a probability in \code{[0,1]} to the corresponding approximate
  quantile.

  For a polynomial histogram, a plain function with the same calling
  conventions.  It is not a step function, so it is plotted with
  \code{\link[graphics]{curve}} or \code{plot(F, from, to)} rather than
  \code{plot.ecdf}.
}
\references{
  Stokely, M., Hesterberg, T., Coehlo, N., and Merchant, A.  Worst-Case
  Information Loss Bounds for Moment-Augmented Histograms in Cloud
  Telemetry.  Manuscript, 2026.
}
\author{Murray Stokely \email{murray@stokely.org}}
\seealso{
\code{\link{HistogramTools-package}},
\code{\link{PolynomialHistogram}},
\code{\link{ApproxQuantile}},
\code{\link{ecdf}},
\code{\link[stats]{approxfun}}, and
\code{\link[graphics]{hist}}.
}
\examples{
h <- hist(runif(100), plot=FALSE)
plot(HistToEcdf(h))

# With per-bucket moments the reconstructed CDF follows the data inside
# each bucket instead of stepping at the breaks.
set.seed(1)
x <- rlnorm(5000)
breaks <- c(0, 0.5, 1, 2, 4, 8, 16, 64)
h2 <- PolynomialHistogram(x, breaks, moments=2)
F <- HistToEcdf(h2)
F(c(0.5, 1, 1.5, 3))
plot(ecdf(x), xlim=c(0, 8), do.points=FALSE, main="Reconstructed CDF")
curve(F, 0, 8, add=TRUE, col="red", n=1001)
plot(HistToEcdf(hist(x, breaks, plot=FALSE)), add=TRUE, col="blue",
     do.points=FALSE)
Q <- HistToEcdf(h2, inverse=TRUE)
Q(c(0.5, 0.9, 0.99))
}
\keyword{manip}
\keyword{utilities}
