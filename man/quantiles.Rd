% Copyright 2012 Google Inc. All Rights Reserved.
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

\name{quantile}
\alias{ApproxQuantile}
\alias{ApproxQuantile.histogram}
\alias{ApproxQuantile.polynomial.histogram}
\alias{ApproxMean}
\alias{ApproxMean.histogram}
\alias{ApproxMean.polynomial.histogram}
\alias{Count}
\title{Histogram Approximate Quantiles.}
\description{
Approximate the quantiles of the underlying distribution for which only
a histogram is available.
}
\usage{
ApproxQuantile(x, probs, \dots)
\method{ApproxQuantile}{histogram}(x, probs, \dots)
\method{ApproxQuantile}{polynomial.histogram}(x, probs, \dots)
ApproxMean(x, \dots)
\method{ApproxMean}{histogram}(x, \dots)
\method{ApproxMean}{polynomial.histogram}(x, \dots)
Count(x)
}
\arguments{
  \item{x}{A histogram object (created by \code{\link[graphics]{hist}})
    or a polynomial histogram (created by
    \code{\link{PolynomialHistogram}}) that also stores per-bucket
    moments.}
  \item{probs}{Numeric vector of probabilities with values in [0,1].}
  \item{\dots}{For ordinary histograms, any other arguments to pass to
    \code{\link[Hmisc]{wtd.quantile}}; unused otherwise.}
}
\details{
Many data analysis pipelines write out histogram protocol buffers with
thousands of buckets so as to be applicable in a wide range of
contexts.  This function provides a way to transform the histogram into
approximations of the quantile, median, mean, etc of the underlying
distribution.  \code{Count(x)} returns the number of
observations in the histogram.

For an ordinary histogram, \code{ApproxMean} is the count-weighted
mean of the bucket midpoints and \code{ApproxQuantile} the weighted
quantile of the midpoints; both are accurate only with many narrow
buckets.

For a polynomial histogram the stored bucket means pool to the exact
data mean, so \code{ApproxMean} is exact up to rounding.
\code{ApproxQuantile} reads each quantile in three steps: locate the
bucket holding that rank from the cumulative counts, fill that bucket
with the maximum-entropy density that has its stored mean and
variance (see \code{\link{HistToEcdf}}), and invert.  Only the bucket
holding the quantile is fitted.  The rank error of the result never
exceeds the width of the EMDCC envelope at that point, and a bucket
whose moments determine its contents exactly (a single repeated value,
or values on both bucket edges) returns the exact quantile.
}
\value{
  \code{Count} returns a single numeric value giving the total number
  of observations summarized by the histogram (the sum of
  \code{x$counts}).  \code{ApproxMean} returns a single numeric value
  giving an approximation of the mean of the underlying distribution
  computed as a weighted mean of the bucket midpoints.
  \code{ApproxQuantile} returns a numeric vector the same length as
  \code{probs} containing approximate quantiles of the underlying
  distribution, named as \code{\link[stats]{quantile}} names them.
}
\keyword{methods}
\author{Murray Stokely \email{murray@stokely.org}}
\seealso{
\code{\link{HistogramTools-package}},
\code{\link{PolynomialHistogram}},
\code{\link{HistToEcdf}}, and
\code{\link[graphics]{hist}}.
}
\examples{
x <- hist(c(1,2,3), breaks=c(0,1,2,3,4,5,6,7,8,9), plot=FALSE)
Count(x)
ApproxMean(x)
ApproxQuantile(x, .5)
ApproxQuantile(x, c(.05, .95))

# Eight wide buckets with stored moments against the true quantiles.
set.seed(1)
y <- rlnorm(10000)
breaks <- c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 64)
p <- c(0.5, 0.9, 0.99)
rbind(true=quantile(y, p),
      count.only=ApproxQuantile(hist(y, breaks, plot=FALSE), p),
      one.moment=ApproxQuantile(PolynomialHistogram(y, breaks, moments=1), p),
      two.moments=ApproxQuantile(PolynomialHistogram(y, breaks, moments=2), p))
c(mean(y), ApproxMean(PolynomialHistogram(y, breaks)))
}
\keyword{manip}
\keyword{utilities}
