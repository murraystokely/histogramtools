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

\name{InformationLoss}
\alias{KSDCC}
\alias{EMDCC}
\alias{EMDCC.histogram}
\alias{EMDCC.polynomial.histogram}
\alias{EMDCCPerBucket}
\alias{PlotKSDCC}
\alias{PlotEMDCC}
\title{Information Loss Metrics for Histograms}
\description{
Computes a metric between 0 and 1 of the amount of information lost
about the underlying distribution of data for a given histogram, with
or without per-bucket moments.
}
\usage{
KSDCC(h)
EMDCC(h, \dots)
\method{EMDCC}{histogram}(h, \dots)
\method{EMDCC}{polynomial.histogram}(h, \dots)
EMDCCPerBucket(alpha=NULL, tau=NULL)
PlotKSDCC(h, arrow.size.scale=1, main=paste("KSDCC =", KSDCC(h)), \dots)
PlotEMDCC(h, main=paste("EMDCC =", EMDCC(h)), \dots)
}
\arguments{
  \item{h}{A \code{"histogram"} object (created by
    \code{\link[graphics]{hist}}) representing a pre-binned dataset on
    which we'd like to calculate the information loss due to binning,
    or for \code{EMDCC} a \code{"polynomial.histogram"} object (created
    by \code{\link{PolynomialHistogram}}) that also stores per-bucket
    moments.}
  \item{alpha}{Numeric vector of bucket means, each rescaled to the
    unit interval as \eqn{(\mu - a) / (b - a)} for a bucket
    \eqn{[a, b]}, or \code{NULL} for buckets with no stored moments.}
  \item{tau}{Numeric vector of bucket variances, each rescaled as
    \eqn{\sigma^2 / (b - a)^2} and at most \eqn{\alpha (1 - \alpha)}, or
    \code{NULL} for buckets with only the mean stored.}
  \item{arrow.size.scale}{specifies a size scaling factor for the arrow
    illustrating the point of Kolmogorov-Smirnov distance between the
    two e.c.d.fs}
  \item{main}{if 'method="constant"' a number between 0 and 1 inclusive,
    indicating a compromise between left- and right-continuous
    step functions.  See ?approxfun}
  \item{\ldots}{For the plotting functions, any other arguments to pass
    to \code{\link[graphics]{plot}}; unused by \code{EMDCC}.}
}
\details{
The \code{KSDCC} (Kolmogorov-Smirnov Distance of the Cumulative Curves)
function provides the Kolmogorov-Smirnov distance between the empirical
distribution functions of the smallest and largest datasets that could
be represented by the binned data in the provided histogram.  This
quantity is also called the Maximum Displacement of the Cumulative
Curves (MDCC) in the computer science performance evaluation community (see
references).

The \code{EMDCC} (Earth Mover's Distance of the Cumulative Curves)
function is like the Kolmogorov-Smirnov statistic, but uses an integral
to capture the difference across all points of the curve rather than
just the maximum difference.  This is also known as Mallows distance, or
Wasserstein distance with \eqn{p=1}.

\code{EMDCC} is generic.  For an ordinary histogram every bucket's
mass can sit anywhere within the bucket, so the two extreme CDFs differ
by the bucket's full mass across its full width and the metric reduces
to the count-weighted mean bucket width over the value range.  For a
\code{"polynomial.histogram"}, which stores the mean and optionally the
variance of each bucket, only distributions with those moments are
consistent with the data, and the envelope between the highest and
lowest consistent CDFs is narrower.  Its width integrates to a closed
form per bucket (the reference, Theorems 3.1 to 3.3), implemented by
\code{EMDCCPerBucket} on the bucket rescaled to \eqn{[0, 1]}:
\itemize{
  \item no moments: 1;
  \item one moment, mean \eqn{\alpha}:
    \eqn{-\alpha \log \alpha - (1 - \alpha) \log (1 - \alpha)}, which is 0
    for a mean on a bucket edge (a point mass) and \eqn{\log 2} at the
    center;
  \item two moments, mean \eqn{\alpha} and variance \eqn{\tau}: the
    integral of \eqn{\tau / (\tau + (x - \alpha)^2)} outside
    \eqn{[c_1, c_2]} and of \eqn{(\alpha (1 - \alpha) - \tau) / (x (1 - x))}
    inside it, where \eqn{c_1 = \alpha - \tau / (1 - \alpha)} and
    \eqn{c_2 = \alpha + \tau / \alpha}, which is 0 both for zero variance
    (a point mass at the mean) and for the largest feasible variance
    \eqn{\alpha (1 - \alpha)} (two atoms on the edges).
}
The EMDCC of a polynomial histogram is the count-weighted sum of these
per-bucket values times the bucket widths, over the total count and the
value range, exactly as for ordinary histograms.  It is a bound on the
\eqn{L_1} distance between the true CDF and any CDF reconstructed from
the stored statistics, independent of the data's shape.

The \code{PlotKSDCC} and \code{PlotEMDCC} functions take a histogram and
generate a plot showing a geometric representation of the information
loss metrics for the provided histogram.
}
\value{
  \code{KSDCC} and \code{EMDCC} each return a single numeric value
  between 0 and 1 quantifying the information loss for the supplied
  histogram.  \code{EMDCCPerBucket} returns a numeric vector the same
  length as \code{alpha} (or the single value 1 when \code{alpha} is
  \code{NULL}) of per-bucket worst-case information loss on the unit
  bucket, with \code{NA} where \code{alpha} or \code{tau} is \code{NA}.  \code{PlotKSDCC} and \code{PlotEMDCC} have no return
  value and are called for the side effect of drawing an annotated
  ECDF plot on the current graphics device.
}
\references{Douceur, John R., and William J. Bolosky. "A large-scale
  study of file-system contents." \emph{ACM SIGMETRICS Performance
  Evaluation Review} \bold{27.1} (1999): 59-70.

  Stokely, M., Hesterberg, T., Coehlo, N., and Merchant, A.  Worst-Case
  Information Loss Bounds for Moment-Augmented Histograms in Cloud
  Telemetry.  Manuscript, 2026.}
\author{Murray Stokely \email{murray@stokely.org}; the polynomial
  histogram bounds with Tim Hesterberg}
\seealso{
\code{\link{HistogramTools-package}},
\code{\link{PolynomialHistogram}},
\code{\link{ecdf}}, and
\code{\link[graphics]{hist}}.
}
\examples{
x <- rexp(1000)
h <- hist(x, breaks=c(0,1,2,3,4,8,16,32), plot=FALSE)
KSDCC(h)

# For small enough data sets we can construct the two extreme data sets
# that can be constructed from a histogram.  One assuming every data point
# is on the left boundary of its bucket, and another assuming every data
# point is on the right boundary of its bucket.  Our KSDCC metric for
# histograms is equivalent to the ks.test statistics for these two
# extreme data sets.

x.min <- rep(head(h$breaks, -1), h$counts)
x.max <- rep(tail(h$breaks, -1), h$counts)
ks.test(x.min, x.max, exact=FALSE)

PlotKSDCC(h)

EMDCC(h)
PlotEMDCC(h)

# Storing each bucket's mean, or its mean and variance, narrows the
# envelope of CDFs consistent with the histogram.
h1 <- PolynomialHistogram(x, breaks=c(0,1,2,3,4,8,16,32), moments=1)
h2 <- PolynomialHistogram(x, breaks=c(0,1,2,3,4,8,16,32), moments=2)
c(EMDCC(h), EMDCC(h1), EMDCC(h2))

# The per-bucket closed forms on the unit bucket: a stored mean beats
# doubling the buckets (per-bucket bound 1/2) only when the mean lies in
# the outer fifth of the bucket.
alpha <- seq(0, 1, by=0.1)
EMDCCPerBucket(alpha)
EMDCCPerBucket(alpha=rep(0.5, 5), tau=c(0, 0.025, 0.1, 0.225, 0.25))
}
