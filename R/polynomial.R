# Copyright 2026 Murray Stokely and Tim Hesterberg. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Authors: murray@stokely.org (Murray Stokely)
#          timhesterberg@gmail.com (Tim Hesterberg)

# Polynomial histograms.
#
# A polynomial histogram is a fixed-grid histogram whose buckets store,
# alongside the count, the first one or two moments of the observations
# that fell in them.  It is represented here as an object of class
# c("polynomial.histogram", "histogram"): an ordinary R histogram (so
# every function that reads breaks and counts keeps working) with the
# additional fields
#
#   means      per-bucket sample mean, NA where the bucket is empty;
#   variances  per-bucket population variance, NA where the bucket is
#              empty; absent when only the mean is stored.
#
# The number of stored moments is 2 when variances is present and 1
# otherwise.  In the notation of Stokely, Hesterberg, Coehlo and
# Merchant (2027) these are H(K, 1) and H(K, 2); a plain histogram is
# H(K, 0).
#
# The worst-case information loss of a polynomial histogram, EMDCC, is
# the L1 distance between the highest and lowest CDFs consistent with
# the stored statistics.  It is additive over buckets and has a closed
# form per bucket (EMDCCPerBucket below).  The CDF and quantile
# estimators that use the stored moments are the HistToEcdf,
# ApproxQuantile and ApproxMean methods at the end of this file, built
# on the maximum-entropy in-bucket reconstruction in R/maxent.R.

# ----- Constructors --------------------------------------------------------

.PolynomialBucketIndex <- function(x, breaks, right, include.lowest) {
  # Assigns observations to buckets exactly as hist.default() does.
  #
  # hist.default() does not bin on the breaks themselves but on breaks
  # nudged by a small fuzz, so that an observation lying on a break by
  # floating-point accident is binned robustly.  The per-bucket moments
  # have to be computed on the same assignment as the counts, so the
  # fuzz is reproduced here (fuzz = 1e-7 is hist.default's default).
  #
  # Args:
  #   x:              Numeric vector of observations.
  #   breaks:         The histogram's breaks.
  #   right:          As in hist(): TRUE for (a, b] buckets, FALSE for [a, b).
  #   include.lowest: As in hist(): whether the outer edge is included.
  #
  # Returns:
  #   An integer vector of bucket indices in 1:(length(breaks) - 1).
  nB <- length(breaks)
  h <- diff(breaks)
  diddle <- 1e-7 * if (nB > 5) median(h) else if (nB <= 3) diff(range(x)) else min(h[h > 0])
  fuzz <- if (right) {
    c(if (include.lowest) -diddle else diddle, rep.int(diddle, nB - 1L))
  } else {
    c(rep.int(-diddle, nB - 1L), if (include.lowest) diddle else -diddle)
  }
  k <- findInterval(x, breaks + fuzz, left.open = right)
  if (any(k < 1L | k >= nB)) {
    stop("some 'x' not counted; maybe 'breaks' do not span range of 'x'")
  }
  k
}

PolynomialHistogram <- function(x, breaks="Sturges", moments=1, right=TRUE,
                                include.lowest=TRUE,
                                xname=deparse1(substitute(x))) {
  # Builds a polynomial histogram from raw data.
  #
  # The breaks and counts come from hist(), with the same bucket
  # convention arguments, so that moments=0 returns exactly what hist()
  # returns.  right=FALSE gives the [a, b) buckets with the last one
  # closed used in the paper.
  #
  # Args:
  #   x:              Numeric vector of observations.  Non-finite values
  #                   are dropped, as hist() does.
  #   breaks:         Anything hist() accepts: a vector of breaks, a
  #                   number of buckets, or the name of an algorithm.
  #   moments:        0, 1 or 2: the number of per-bucket moments to store.
  #   right:          TRUE for (a, b] buckets, FALSE for [a, b).
  #   include.lowest: Whether the outer edge belongs to the end bucket.
  #   xname:          A name for the data, used in plots.
  #
  # Returns:
  #   For moments=0 an object of class "histogram"; otherwise an object
  #   of class c("polynomial.histogram", "histogram") with per-bucket
  #   means and, for moments=2, variances.
  xname <- xname                      # capture the name before x is filtered
  stopifnot(is.numeric(x))
  stopifnot(length(moments) == 1, moments %in% 0:2)
  x <- x[is.finite(x)]
  h <- hist(x, breaks = breaks, right = right, include.lowest = include.lowest,
            plot = FALSE)
  h$xname <- xname
  if (moments == 0) return(h)
  K <- length(h$counts)
  k <- .PolynomialBucketIndex(x, h$breaks, right, include.lowest)
  n <- tabulate(k, K)
  if (!identical(n, h$counts)) {
    stop("internal error: bucket assignment disagrees with hist()")
  }
  sums <- numeric(K)
  agg <- rowsum(x, k)
  sums[as.integer(rownames(agg))] <- agg[, 1]
  means <- ifelse(n > 0, sums / n, NA_real_)
  variances <- NULL
  if (moments == 2) {
    # Two passes, about the bucket mean, so that a bucket far from the
    # origin keeps full precision in its variance.
    agg <- rowsum((x - means[k])^2, k)
    ss <- numeric(K)
    ss[as.integer(rownames(agg))] <- agg[, 1]
    variances <- ifelse(n > 0, ss / n, NA_real_)
  }
  # The fuzz lets hist() count an observation lying within 1e-7 of a
  # bucket width past a break in the bucket it was meant for, so a
  # bucket's mean can sit outside the bucket by that much.  Such a mean
  # is on the edge to the precision of the stored breaks; store it there.
  means <- pmin(pmax(means, head(h$breaks, -1)), tail(h$breaks, -1))
  PreBinnedPolynomialHistogram(h$breaks, h$counts, means, variances, xname)
}

# A stored variance may exceed the largest variance any distribution in
# the bucket can have, alpha (1 - alpha) w^2, by rounding: pipelines
# compute it as E[x^2] - mean^2 in double precision.  Excess within this
# relative tolerance is clamped; more is an error.  It is the solver's
# own variance tolerance (.ME.TOL.VAR in maxent.R).
.VARIANCE.TOL <- 1e-9

PreBinnedPolynomialHistogram <- function(breaks, counts, means, variances=NULL,
                                         xname="") {
  # Builds a polynomial histogram from pre-aggregated per-bucket
  # statistics, as telemetry pipelines deliver them.
  #
  # Args:
  #   breaks:    Numeric vector of K + 1 breaks.
  #   counts:    Numeric vector of K counts.
  #   means:     Numeric vector of K per-bucket means.  Each must lie
  #              within its bucket; the value for an empty bucket is
  #              ignored and stored as NA.
  #   variances: Optional numeric vector of K per-bucket population
  #              variances, each in [0, alpha (1 - alpha) w^2] where
  #              alpha is the mean's position in the bucket of width w.
  #   xname:     A name for the data, used in plots.
  #
  # Returns:
  #   An object of class c("polynomial.histogram", "histogram").
  h <- PreBinnedHistogram(breaks, counts, xname)
  stopifnot(is.numeric(means), length(means) == length(counts))
  stopifnot(all(diff(breaks) > 0))
  empty <- counts == 0
  lo <- head(breaks, -1)
  hi <- tail(breaks, -1)
  if (any(is.na(means[!empty]))) {
    stop("Each non-empty bucket must have a mean.")
  }
  if (any(means[!empty] < lo[!empty] | means[!empty] > hi[!empty])) {
    stop("Each mean must lie within its bucket.")
  }
  means[empty] <- NA_real_
  if (!is.null(variances)) {
    stopifnot(is.numeric(variances), length(variances) == length(counts))
    if (any(is.na(variances[!empty]))) {
      stop("Each non-empty bucket must have a variance.")
    }
    if (any(variances[!empty] < 0)) {
      stop("Each variance must be non-negative.")
    }
    vmax <- (means - lo) * (hi - means)
    over <- !empty & variances > vmax
    if (any(variances[over] > vmax[over] + .VARIANCE.TOL * (hi - lo)[over]^2)) {
      stop("Each variance must be at most (mean - a) * (b - mean) for its bucket [a, b].")
    }
    variances[over] <- vmax[over]
    variances[empty] <- NA_real_
  }
  h$means <- means
  if (!is.null(variances)) h$variances <- variances
  class(h) <- c("polynomial.histogram", class(h))
  h
}

# ----- Merging ---------------------------------------------------------------

.PoolBucketMoments <- function(counts, means, variances=NULL, group) {
  # Pools per-bucket statistics over groups of buckets: the counts add,
  # the means combine count-weighted, and the variances combine as the
  # count-weighted mean of the within-group variances plus the
  # count-weighted spread of the group means about the pooled mean
  # (the parallel-variance formula of Chan, Golub and LeVeque, 1979).
  # The pooled statistics are exactly those of the union of the
  # observations, so merging polynomial histograms loses nothing and
  # the EMDCC of the result is the bound for the merged data.
  #
  # Args:
  #   counts:    Numeric vector of bucket counts.
  #   means:     Numeric vector of bucket means (NA where empty).
  #   variances: Optional numeric vector of bucket variances.
  #   group:     Integer vector, the group of each bucket, in 1:G.
  #
  # Returns:
  #   A list with counts, means and (if given) variances of length G.
  G <- max(group)
  n <- counts
  m <- ifelse(n > 0, means, 0)
  N <- as.vector(rowsum(n, group, reorder = TRUE))
  present <- sort(unique(group))
  N.full <- numeric(G)
  N.full[present] <- N
  M <- numeric(G)
  M[present] <- as.vector(rowsum(n * m, group, reorder = TRUE))
  M <- ifelse(N.full > 0, M / N.full, NA_real_)
  out <- list(counts = N.full, means = M)
  if (!is.null(variances)) {
    v <- ifelse(n > 0, variances, 0)
    V <- numeric(G)
    V[present] <- as.vector(rowsum(n * (v + (m - ifelse(is.na(M), 0, M)[group])^2),
                                   group, reorder = TRUE))
    out$variances <- ifelse(N.full > 0, V / N.full, NA_real_)
  }
  out
}

.MomentOrder <- function(h) {
  # 0 for a count-only histogram, 1 or 2 for a polynomial histogram.
  if (!inherits(h, "polynomial.histogram")) return(0L)
  if (is.null(h$variances)) 1L else 2L
}

# ----- Worst-case information loss ------------------------------------------

EMDCCPerBucket <- function(alpha=NULL, tau=NULL) {
  # Worst-case EMDCC of one bucket rescaled to [0, 1], given the bucket
  # mean alpha and variance tau in those coordinates.  Vectorized over
  # buckets; the number of moments is inferred from the arguments
  # supplied.  EMDCC is additive across buckets, so each is treated on
  # its own with endpoints 0 and 1, and the value is the integral over
  # [0, 1] of Fmax - Fmin, the envelope of all distributions matching
  # the stored moments.
  #
  # No moments: with no information beyond the count, Fmax puts the
  #   whole mass at the left edge and Fmin at the right edge, so the
  #   envelope is 1 across the unit bucket.
  # One moment (Theorem 3.1 of the paper): the envelope is
  #   (1 - alpha) / (1 - x) below the mean and alpha / x above it;
  #   integrating gives
  #       -(1 - alpha) log(1 - alpha) - alpha log(alpha).
  #   A mean at an edge is a point mass, so 0.
  # Two moments (Theorems 3.2 and 3.3): the envelope is
  #   p3 = tau / (tau + (x - alpha)^2) for x outside [c1, c2] and
  #   p4 = (alpha (1 - alpha) - tau) / (x (1 - x)) inside it, where
  #   c1 = alpha - tau / (1 - alpha) and c2 = alpha + tau / alpha.  With
  #   antiderivatives P3 = sqrt(tau) atan((x - alpha) / sqrt(tau)) and
  #   P4 = (alpha (1 - alpha) - tau) log(x / (1 - x)),
  #       EMDCC = [P3]_0^c1 + [P4]_c1^c2 + [P3]_c2^1.
  #   Zero variance is a point mass and variance alpha (1 - alpha) is
  #   the two-point distribution on the edges; both are determined
  #   exactly, so the bound is 0 there.
  #
  # Args:
  #   alpha: Bucket means rescaled to [0, 1], or NULL for no moments.
  #   tau:   Bucket variances rescaled to [0, 1] (at most
  #          alpha (1 - alpha)), or NULL for one moment.
  #
  # Returns:
  #   The per-bucket EMDCC, a number in [0, 1], one per element of
  #   alpha (NA where alpha or tau is NA).
  if (is.null(alpha)) return(1)
  stopifnot(is.numeric(alpha))
  if (any(alpha < 0 | alpha > 1, na.rm = TRUE)) {
    stop("'alpha' must lie in [0, 1].")
  }
  e1 <- -alpha * log(alpha) - (1 - alpha) * log1p(-alpha)
  e1[alpha <= 0 | alpha >= 1] <- 0
  if (is.null(tau)) return(e1)
  stopifnot(is.numeric(tau), length(tau) == length(alpha))
  s2max <- alpha * (1 - alpha)
  if (any(tau < 0 | tau > s2max, na.rm = TRUE)) {
    stop("'tau' must lie in [0, alpha * (1 - alpha)].")
  }
  e2 <- rep(NA_real_, length(alpha))
  ok <- !is.na(alpha) & !is.na(tau)
  e2[ok] <- 0
  i <- which(ok & tau > 0 & tau < s2max)
  if (length(i)) {
    a <- alpha[i]
    t <- tau[i]
    s <- sqrt(t)
    c1 <- a - t / (1 - a)
    c2 <- a + t / a
    P3 <- function(x) s * atan((x - a) / s)
    P4 <- function(x) (a * (1 - a) - t) * log(x / (1 - x))
    e2[i] <- (P3(c1) - P3(0)) + (P4(c2) - P4(c1)) + (P3(1) - P3(c2))
  }
  e2
}

EMDCC.polynomial.histogram <- function(h, ...) {
  # Earth Mover's Distance of the Cumulative Curves for a polynomial
  # histogram: the L1 distance between the highest and lowest CDFs
  # consistent with the stored counts and moments, normalized by the
  # value range as for count-only histograms.
  #
  # Args:
  #   h: A polynomial.histogram object.
  # Return Value:
  #   A number between 0 and 1 giving the EMDCC.
  stopifnot(inherits(h, "polynomial.histogram"))
  w <- diff(h$breaks)
  n <- h$counts
  alpha <- (h$means - head(h$breaks, -1)) / w
  tau <- if (is.null(h$variances)) NULL else h$variances / w^2
  per <- EMDCCPerBucket(alpha, tau)
  per[n == 0] <- 0                     # empty buckets carry no mass
  return(sum(w * n * per) / sum(n) / diff(range(h$breaks)))
}

# ----- CDF and quantile estimation -----------------------------------------
#
# Each bucket is filled with the maximum-entropy density that has the
# bucket's stored mean and variance (R/maxent.R).  Because that density
# has the stored moments, the reconstructed CDF lies inside the EMDCC
# envelope, so its L1 distance from the true CDF never exceeds EMDCC(h)
# and its rank error at any quantile never exceeds the envelope width
# there.  Some stored configurations determine the CDF exactly and need
# no fit: a mean on a bucket edge or a zero variance is a point mass at
# the mean, and the largest feasible variance is two atoms on the edges.

.PolyBucketIndexOfRank <- function(h, p) {
  # The bucket holding the quantile at probability p, and the fraction
  # of that bucket's mass below the quantile.
  cs <- cumsum(h$counts)
  N <- cs[length(cs)]
  target <- p * N
  idx <- which(cs >= target & h$counts > 0)[1]
  if (is.na(idx)) idx <- max(which(h$counts > 0))     # p = 1 on a trailing empty bucket
  cum.before <- if (idx == 1L) 0 else cs[idx - 1L]
  f <- (target - cum.before) / h$counts[idx]
  list(idx = idx, f = min(max(f, 0), 1))
}

HistToEcdf.polynomial.histogram <- function(h, inverse=FALSE, ...) {
  # The CDF reconstructed from a polynomial histogram, or its inverse.
  #
  # Every bucket is fitted once; the returned closure evaluates the CDF
  # (right-continuous, with the atoms of exactly determined buckets as
  # jumps) at any vector of values.
  #
  # Args:
  #   h:       A polynomial.histogram object.
  #   inverse: If TRUE, return the quantile function instead.
  #   ...:     Unused.
  #
  # Returns:
  #   A function.  Called with a numeric vector x it returns the
  #   estimated cumulative probability at each x; with inverse=TRUE,
  #   called with probabilities in [0, 1] it returns the estimated
  #   quantiles.
  stopifnot(inherits(h, "polynomial.histogram"))
  stopifnot(sum(h$counts) > 0)
  K <- length(h$counts)
  desc <- lapply(seq_len(K), function(j) .PolyBucketDescription(h, j))
  breaks <- h$breaks
  width <- diff(breaks)
  cs <- c(0, cumsum(h$counts))
  N <- cs[K + 1L]
  counts <- h$counts
  if (inverse) {
    return(function(p) {
      stopifnot(is.numeric(p), all(p >= 0 & p <= 1, na.rm = TRUE))
      vapply(p, function(pp) {
        if (is.na(pp)) return(NA_real_)
        r <- .PolyBucketIndexOfRank(h, pp)
        .PolyBucketQuantile(desc[[r$idx]], r$f, breaks[r$idx], breaks[r$idx + 1L])
      }, numeric(1))
    })
  }
  function(x) {
    stopifnot(is.numeric(x))
    out <- rep(NA_real_, length(x))
    ok <- which(!is.na(x))
    if (!length(ok)) return(out)
    xx <- x[ok]
    k <- findInterval(xx, breaks, rightmost.closed = TRUE)
    k <- pmin(pmax(k, 1L), K)
    u <- pmin(pmax((xx - breaks[k]) / width[k], 0), 1)
    within <- numeric(length(xx))
    for (j in unique(k)) {
      i <- which(k == j)
      within[i] <- .PolyBucketCdf(desc[[j]], u[i], breaks[j], breaks[j + 1L])
    }
    val <- (cs[k] + within * counts[k]) / N
    val[xx < breaks[1L]] <- 0
    val[xx >= breaks[K + 1L]] <- 1
    out[ok] <- val
    out
  }
}

ApproxQuantile.polynomial.histogram <- function(x, probs, ...) {
  # Quantiles read off a polynomial histogram.
  #
  # For each probability: locate the bucket holding that rank by the
  # cumulative counts, fit that bucket's maximum-entropy density to its
  # stored moments, and invert.  Only the bucket holding the quantile is
  # fitted, so the cost does not grow with the number of buckets.
  #
  # Args:
  #   x:     A polynomial.histogram object.
  #   probs: Numeric vector of probabilities in [0, 1].
  #   ...:   Unused.
  #
  # Returns:
  #   A named numeric vector of quantiles, one per probability, named
  #   as quantile() names them.
  stopifnot(inherits(x, "polynomial.histogram"))
  stopifnot(is.numeric(probs), all(probs >= 0 & probs <= 1, na.rm = TRUE))
  stopifnot(sum(x$counts) > 0)
  q <- vapply(probs, function(p) {
    if (is.na(p)) return(NA_real_)
    r <- .PolyBucketIndexOfRank(x, p)
    .PolyBucketQuantile(.PolyBucketDescription(x, r$idx), r$f,
                        x$breaks[r$idx], x$breaks[r$idx + 1L])
  }, numeric(1))
  names(q) <- paste0(formatC(100 * probs, format = "fg", width = 1,
                             digits = max(2L, getOption("digits"))), "%")
  q
}

ApproxMean.polynomial.histogram <- function(x, ...) {
  # The mean of the data summarized by a polynomial histogram.  The
  # stored bucket means pool to the exact data mean, so nothing is
  # approximated beyond floating-point rounding.
  stopifnot(inherits(x, "polynomial.histogram"))
  return(sum(x$counts * x$means, na.rm = TRUE) / sum(x$counts))
}
