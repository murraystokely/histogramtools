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

\name{MergeBuckets}
\alias{MergeBuckets}
\title{Merge adjacent buckets in a histogram to create a new histogram.}
\description{
Merge adjacent buckets in a histogram and return a new histogram.
}
\usage{
  MergeBuckets(x, adj.buckets=NULL, breaks=NULL, FUN=sum)
}
\arguments{
  \item{x}{A histogram object (created by \code{\link[graphics]{hist}}).}
  \item{adj.buckets}{The number of adjacent buckets to merge together.}
  \item{breaks}{If \code{adj.buckets} is equal to \code{NULL}, then this
    argument either specifies a vector giving the breakpoints between cells, or a
    single number giving the total number of cells in the new
    histogram.  The vector of new buckets must have the same range as
    the original break list.  If a single number is provided it must be
    less than \code{length(x$breaks)}.}
  \item{FUN}{The user defined function that should be run to merge the
    counts of adjacent buckets in the histogram.}
}
\details{
Many data analysis pipelines write out histogram protocol buffers with
thousands of buckets so as to be applicable in a wide range of
contexts.  This function provides a way to transform the histogram into
one with fewer buckets.
}
\keyword{methods}
\author{Murray Stokely \email{mstokely@google.com}}
\seealso{
\code{\link{histogramtools-package}} and
\code{\link[graphics]{hist}}.
}
\examples{
hist.1 <- hist(c(1,2,3), breaks=c(0,1,2,3,4,5,6,7,8,9), plot=FALSE)
hist.2 <- MergeBuckets(hist.1, adj.buckets=2)

hist.1
hist.2
}
\keyword{manip}
\keyword{utilities}
