% Copyright 2011 Google Inc. All Rights Reserved.
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

\name{AddHistograms}
\alias{AddHistograms}
\alias{.NewHistogramName}
\title{Aggregate histograms that have identical breaks.}
\description{
Aggregate histogram objects that have identical breaks.
}
\usage{
  AddHistograms(\dots, x=list(\dots), main=.NewHistogramName(x))
  .NewHistogramName(x)
}
\arguments{
  \item{\dots}{Histogram objects (created by \code{\link[graphics]{hist}}),
    or polynomial histograms (created by
    \code{\link{PolynomialHistogram}}) that all store the same number
    of per-bucket moments.}
  \item{x}{A list of histogram objects. If \code{x} is supplied then the
  \ldots are ignored.}
  \item{main}{A title to add to the aggregated/merged histogram.  By default, if
    two histograms are provided a title will be created that includes
    the names of the original histograms.  If more histograms are
    provided the title will simply include the number of aggregated histograms.
  }
}
\details{
This function adds the buckets of the provided histograms to return a single
aggregated histogram.

Polynomial histograms merge exactly.  Their counts add, their bucket
means combine count-weighted, and their bucket variances combine by the
parallel-variance formula (the count-weighted mean of the within-bucket
variances plus the count-weighted spread of the bucket means about the
pooled mean), so the result carries the statistics of the union of the
observations and its \code{\link{EMDCC}} is the worst-case bound for
the merged data.  This is the property that lets per-process
histograms be summed into a fleet-wide one with no loss.  Count-only
and polynomial histograms, or polynomial histograms with different
numbers of stored moments, cannot be added together.

\code{.NewHistogramName} is a utility that takes the list of histogram
objects to be aggregated and returns a name for the new merged
histogram.  It is normally hidden, but can be viewed using
\code{HistogramTools:::.NewHistogramName}.
}
\value{
  \code{AddHistograms} returns an object of class \code{"histogram"}
  whose \code{breaks} match those of the inputs and whose \code{counts}
  are the elementwise sum of the counts of the supplied histograms.
  For polynomial histograms the result is of class
  \code{c("polynomial.histogram", "histogram")} with the pooled
  \code{means} and \code{variances}.
  \code{.NewHistogramName} returns a single character string suitable
  for use as the \code{xname} of the merged histogram.
}
\author{Murray Stokely \email{murray@stokely.org}}
\seealso{
\code{\link{HistogramTools-package}} and
\code{\link[graphics]{hist}}.
}
\examples{
hist.1 <- hist(c(1,2,3,4), plot=FALSE)
hist.2 <- hist(c(1,2,2,4), plot=FALSE)
hist.sum <- AddHistograms(hist.1, hist.2)
hist.3 <- hist(c(1,2,2,4), plot=FALSE)
hist.sum <- AddHistograms(hist.1, hist.2, hist.3)

# Polynomial histograms built on two halves of a data set merge to the
# histogram of the whole.
set.seed(1)
x <- rlnorm(2000)
breaks <- c(0, 0.5, 1, 2, 4, 8, 64)
p1 <- PolynomialHistogram(x[1:1000], breaks, moments=2)
p2 <- PolynomialHistogram(x[1001:2000], breaks, moments=2)
p12 <- AddHistograms(p1, p2)
all.equal(p12$means, PolynomialHistogram(x, breaks, moments=2)$means)
c(EMDCC(p1), EMDCC(p2), EMDCC(p12))
}
\keyword{manip}
\keyword{utilities}
