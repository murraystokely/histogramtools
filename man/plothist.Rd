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

\name{PlotRelativeFrequency}
\alias{PlotRelativeFrequency}
\title{Plot Relative Frequency Histogram}
\description{
Produces a relative frequency histogram.
}
\usage{
PlotRelativeFrequency(x, ylab="Relative Frequency", \dots)
}
\arguments{
  \item{x}{A \code{"histogram"} object (created by
    \code{\link[graphics]{hist}}) representing a pre-binned dataset.}
  \item{ylab}{y-axis label for the plot.}
  \item{\dots}{Additional parameters are passed to
    \code{\link[graphics]{plot}()}.}
}
\details{
The default \code{\link[graphics]{plot.histogram}} function supports Frequency or
Density plots, but does not provide a way to produce a relative
frequency histogram.  This function plots this type of histogram.
}
\author{Murray Stokely \email{mstokely@google.com}}
\seealso{
\code{\link{histogramtools-package}},
\code{\link[graphics]{hist}},
\code{\link[graphics]{plot.histogram}}.
}
\examples{
x <- runif(100)
h <- hist(x, plot=FALSE)
PlotRelativeFrequency(h)
}
