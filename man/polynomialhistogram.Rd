% Copyright 2026 Murray Stokely and Tim Hesterberg. All Rights Reserved.
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
% Authors: murray@stokely.org (Murray Stokely)
%          timhesterberg@gmail.com (Tim Hesterberg)

\name{PolynomialHistogram}
\alias{PolynomialHistogram}
\alias{PreBinnedPolynomialHistogram}
\alias{polynomial.histogram}
\title{Polynomial Histograms: Histograms with Per-Bucket Moments}
\description{
Build a polynomial histogram, a histogram whose buckets store the first
one or two moments of their observations alongside the count, either
from raw data or from pre-aggregated per-bucket statistics.
}
\usage{
PolynomialHistogram(x, breaks="Sturges", moments=1, right=TRUE,
                    include.lowest=TRUE, xname=deparse1(substitute(x)))
PreBinnedPolynomialHistogram(breaks, counts, means, variances=NULL,
                             xname="")
}
\arguments{
  \item{x}{A numeric vector of observations.  Non-finite values are
    dropped, as \code{\link[graphics]{hist}} does.}
  \item{breaks}{For \code{PolynomialHistogram}, anything
    \code{\link[graphics]{hist}} accepts: a vector of breakpoints, a
    single number of buckets, or the name of a binning algorithm.  For
    \code{PreBinnedPolynomialHistogram}, a numeric vector of \eqn{K + 1}
    breakpoints.}
  \item{moments}{The number of per-bucket moments to store: 0, 1 or 2.
    With 0 the result is an ordinary \code{"histogram"} object.}
  \item{right}{As in \code{\link[graphics]{hist}}: \code{TRUE} for
    buckets of the form \eqn{(a, b]}, \code{FALSE} for \eqn{[a, b)}.}
  \item{include.lowest}{As in \code{\link[graphics]{hist}}: whether an
    observation on the outermost breakpoint belongs to the end bucket.}
  \item{xname}{A character string naming the data, used in plots.}
  \item{counts}{A numeric vector of \eqn{K} bucket counts.}
  \item{means}{A numeric vector of \eqn{K} per-bucket means.  Each must
    lie within its bucket.  The value given for an empty bucket is
    ignored and stored as \code{NA}.}
  \item{variances}{An optional numeric vector of \eqn{K} per-bucket
    population variances (the mean squared deviation from the bucket
    mean, not the unbiased sample variance).  Each must lie between 0
    and the largest variance any distribution on the bucket with that
    mean can have, \eqn{(\mu - a)(b - \mu)} for a bucket \eqn{[a, b]}
    with mean \eqn{\mu}.  A value above that limit by no more than
    rounding error is clamped to it.}
}
\details{
A polynomial histogram is a fixed-grid histogram that keeps, for each
bucket, the sample mean of the observations that fell in it and
optionally their variance.  Telemetry pipelines produce these
statistics at the same cost as the count, and they tighten the
worst-case bound on what the histogram can say about the underlying
distribution (see \code{\link{EMDCC}}) and sharpen the quantiles read
from it.  In the notation of the reference below, a histogram with
\eqn{K} buckets and \eqn{p} stored moments is \eqn{H(K, p)}; an
ordinary histogram is \eqn{H(K, 0)}.

The result is an object of class \code{c("polynomial.histogram",
"histogram")}: an ordinary R histogram with the additional components
\code{means} and, when two moments are stored, \code{variances}.  Every
function in this package and in base R that reads the breaks and
counts of a histogram works on it unchanged.  The number of stored
moments is 2 when \code{variances} is present and 1 otherwise.

\code{PolynomialHistogram} obtains its breaks and counts from
\code{\link[graphics]{hist}} with the same bucket-convention arguments,
and computes the moments on exactly the bucket assignment that produced
the counts, so that \code{moments=0} returns precisely what
\code{hist(x, breaks, right, include.lowest, plot=FALSE)} returns.
\code{right=FALSE} gives the \eqn{[a, b)} buckets, with the last one
closed, used in the reference.

\code{PreBinnedPolynomialHistogram} is the constructor for statistics
aggregated elsewhere, as \code{\link{PreBinnedHistogram}} is for
counts.  It validates that each mean lies within its bucket and each
variance within the feasible range for its mean.
}
\value{
  For \code{PolynomialHistogram} with \code{moments=0}, an object of
  class \code{"histogram"}.  Otherwise an object of class
  \code{c("polynomial.histogram", "histogram")}: a list with the
  components of a histogram (\code{breaks}, \code{counts},
  \code{density}, \code{mids}, \code{xname}, \code{equidist}) and
  \code{means}, a numeric vector of per-bucket means with \code{NA} for
  empty buckets, and for two moments \code{variances}, a numeric vector
  of per-bucket population variances with \code{NA} for empty buckets.
}
\references{
  Stokely, M., Hesterberg, T., Coehlo, N., and Merchant, A.  Worst-Case
  Information Loss Bounds for Moment-Augmented Histograms in Cloud
  Telemetry.  Manuscript, 2026.
}
\author{Murray Stokely \email{murray@stokely.org} and Tim Hesterberg}
\seealso{
\code{\link{EMDCC}} for the worst-case information loss of a polynomial
histogram,
\code{\link{PreBinnedHistogram}},
\code{\link{HistogramTools-package}}, and
\code{\link[graphics]{hist}}.
}
\examples{
set.seed(1)
x <- rlnorm(10000)
breaks <- c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)

# The same buckets with zero, one and two stored moments.
h0 <- PolynomialHistogram(x, breaks, moments=0)   # identical to hist()
h1 <- PolynomialHistogram(x, breaks, moments=1)
h2 <- PolynomialHistogram(x, breaks, moments=2)
class(h2)
h2$counts
h2$means
h2$variances

# More stored moments on the same grid never widen the worst-case
# envelope of CDFs consistent with the histogram.
c(EMDCC(h0), EMDCC(h1), EMDCC(h2))

# From statistics aggregated upstream, e.g. by a telemetry pipeline.
p <- PreBinnedPolynomialHistogram(breaks=c(0, 1, 2, 4),
                                  counts=c(120, 80, 30),
                                  means=c(0.6, 1.4, 2.5),
                                  variances=c(0.05, 0.08, 0.3),
                                  xname="latency (ms)")
EMDCC(p)
plot(p)
}
\keyword{manip}
\keyword{utilities}
