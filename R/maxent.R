# Copyright 2026 Murray Stokely. All Rights Reserved.
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
# Author: murray@stokely.org (Murray Stokely)

# Maximum-entropy in-bucket reconstruction for polynomial histograms.
#
# A polynomial histogram stores, for each bucket, the count and the
# first one or two moments of the observations that fell in it.  To
# estimate a CDF or a quantile the bucket has to be filled with some
# distribution consistent with those statistics.  The functions in this
# file fill it with the maximum-entropy density, which assumes nothing
# beyond the stored moments.  Within a bucket rescaled to u in [0, 1]
# the density is
#
#     f(u) = exp(l1 u + l2 u^2) / Z
#
# with (l1, l2) solved so that f has the stored mean alpha and variance
# tau (l2 = 0 when only the mean is stored).  Because f matches the
# stored moments, the reconstructed CDF lies inside the EMDCC envelope
# and its L1 distance from the true CDF never exceeds EMDCC(h).
#
# The exponents are unbounded: the family has to reach a point mass
# (tau -> 0, or a mean at an edge) and a pair of atoms at the edges
# (tau -> alpha (1 - alpha)).  A bucket whose stored moments sit within
# .ME.DEGENERATE.TOL of one of those configurations is reconstructed as
# that configuration exactly, with no fit.
#
# Numerics.  The exponent is shifted by its maximum over [0, 1] before
# exponentiation.  Integrals use a 16-point Gauss-Legendre rule on
# panels whose boundaries are placed dyadically at the density's own
# scales (sd = 1 / sqrt(-2 l2) around the vertex when l2 < 0, 1 / |slope|
# from an edge the density decays away from), so a peak 1e-12 wide is
# resolved with the same ~60 panels as a flat one.  Moments are taken
# about the fitted mean so the variance keeps full relative precision
# however small it is.  The same quadrature defines the fitted moments
# and the reconstructed CDF.
#
# None of these functions are exported; the user-facing methods are
# HistToEcdf.polynomial.histogram and ApproxQuantile.polynomial.histogram
# in polynomial.R.

# 16-point Gauss-Legendre nodes and weights on [0, 1].  These are
# mathematical constants (the roots of the Legendre polynomial P_16 and
# the matching weights).  Generated 2026-11-20 with
# statmod::gauss.quad(16, "legendre") (Golub-Welsch, Netlib gausq2),
# mapped from [-1, 1] to [0, 1] and symmetrized about 1/2; the weights
# sum to 1 within 1e-15 and the rule integrates u^k exactly for
# k <= 31 to the same precision.  Tabulated here so the package does not
# import a compiled orthogonal-polynomial library for one fixed list of
# 16 numbers; test-maxent.R regenerates the table with statmod when it
# is installed.
.GL16.X <- c(0.0052995325041749197, 0.027712488463383755, 0.067184398806084095,
             0.12229779582249839, 0.19106187779867806, 0.27099161117138632,
             0.35919822461037054, 0.45249374508118134, 0.54750625491881877,
             0.64080177538962957, 0.72900838882861363, 0.80893812220132189,
             0.87770220417750155, 0.93281560119391593, 0.9722875115366163,
             0.99470046749582508)
.GL16.W <- c(0.013576229705877092, 0.031126761969323902, 0.047579255841246372,
             0.062314485627767348, 0.074797994408288396, 0.084578259697501379,
             0.091301707522461653, 0.094725305227534251, 0.094725305227534251,
             0.091301707522461653, 0.084578259697501379, 0.074797994408288396,
             0.062314485627767348, 0.047579255841246372, 0.031126761969323902,
             0.013576229705877092)

# Relative precision of alpha = (mean - a) / width and tau = var / width^2
# as computed from stored sums in double precision: the subtraction
# E[x^2] - mean^2 loses about 1e-16 (mean / width)^2, which is at most
# about 1e-12 for the log-spaced grids telemetry systems use.  A bucket
# whose moments sit within this of a degenerate configuration is
# reconstructed as that configuration.
.ME.DEGENERATE.TOL <- 1e-12
# Solver tolerances on |fitted mean - alpha| / alpha and
# |log(fitted var / tau)|, each floored at the rounding limit of the
# exponent (see .MaxEnt2Floor).
.ME.TOL.MEAN <- 1e-10
.ME.TOL.VAR <- 1e-9

# ---- One moment: f(u) = exp(l u) / Z, everything in closed form ----------

.AlphaOfLambda <- function(l) {
  # Mean of exp(l u) / Z on [0, 1]; stable for any l.
  if (abs(l) < 1e-9) return(0.5)
  if (l > 0) 1 / (-expm1(-l)) - 1 / l else 1 - .AlphaOfLambda(-l)
}

.LambdaOfAlpha <- function(alpha) {
  # Inverse of .AlphaOfLambda on (0, 1).  For l > 0 the mean satisfies
  # alpha(l) >= 1 - 1 / l, which gives the bracket for uniroot.
  if (abs(alpha - 0.5) < 1e-12) return(0)
  if (alpha < 0.5) return(-.LambdaOfAlpha(1 - alpha))
  hi <- 1 / (1 - alpha)
  l <- uniroot(function(l) .AlphaOfLambda(l) - alpha, lower = 0, upper = hi,
               tol = 1e-13 * max(1, hi))$root
  # Newton polish where the derivative (the variance) is well conditioned.
  if (l > 1e-3) {
    for (i in 1:3) {
      v <- 1 / l^2 - 1 / (4 * sinh(l / 2)^2)
      if (!is.finite(v) || v <= 0) break
      l <- l - (.AlphaOfLambda(l) - alpha) / v
    }
  }
  l
}

.MaxEnt1Cdf <- function(u, l) {
  # F(u) = (exp(l u) - 1) / (exp(l) - 1), written without overflow.
  if (abs(l) < 1e-9) return(u)
  if (l > 0) exp(l * (u - 1)) * expm1(-l * u) / expm1(-l)
  else expm1(l * u) / expm1(l)
}

.MaxEnt1Quantile <- function(p, l) {
  # Inverse of .MaxEnt1Cdf.
  if (abs(l) < 1e-9) return(p)
  if (l > 0) 1 + log(p + (1 - p) * exp(-l)) / l
  else log1p(p * expm1(l)) / l
}

.MaxEnt1Panels <- function(l) {
  # Dyadic points at the exponential's scale 1 / |l| from the edge it
  # decays away from, as fractions of the bucket.  Used by callers that
  # integrate the reconstructed CDF numerically.
  if (abs(l) < 1e-9) return(numeric(0))
  s <- 1 / abs(l)
  pts <- s * 2^(-3:60)
  pts <- pts[pts < 1]
  if (l < 0) pts else 1 - pts
}

# ---- Two moments: f(u) = exp(l1 u + l2 u^2) / Z, adaptive quadrature -----
#
# A fit is stored in the form that keeps its parameters exact where the
# density is extreme, and the exponent is evaluated from that form:
#   "lin"     l2 = 0, parameter l1:         E(t) = l1 t
#   "vertex"  l2 < 0, parameters (v, l2):   E(t) = l2 (t - v)^2 + const
#             (a Gaussian peak at v with sd 1 / sqrt(-2 l2))
#   "layers"  l2 > 0, parameters (q, l2):   E(t) = l2 t (t - 1) - q t
#             (exponential layers at both edges; the mass ratio of the
#             layer at 1 to the layer at 0 is about exp(-q) times the
#             ratio of their rates, so q is the exact split variable)
# with l1 = -2 l2 v, respectively l1 = -l2 - q.  Everything downstream
# uses E(t) minus its maximum over [0, 1], computed from the form so that
# no large terms cancel.

.MaxEnt2Params <- function(form, p1, p2) {
  switch(form,
         lin = c(l1 = p1, l2 = 0),
         vertex = c(l1 = -2 * p2 * p1, l2 = p2),
         layers = c(l1 = -p2 - p1, l2 = p2))
}

.MaxEnt2Bounds <- function(form, p1, p2, l1, l2) {
  # Panel boundaries: dyadic points at every scale of the density, and
  # no panel wider than 1/8.
  cand <- numeric(0)
  j <- 2^(-3:60)
  if (form == "vertex") {
    s <- 1 / sqrt(-2 * l2)
    cand <- c(cand, p1 + s * j, p1 - s * j)
  }
  if (form == "layers") cand <- c(cand, 0.5 + p1 / (2 * l2))  # the minimum
  if (l1 < 0) cand <- c(cand, (1 / -l1) * j)                    # decays from 0
  if (l1 + 2 * l2 > 0) cand <- c(cand, 1 - (1 / (l1 + 2 * l2)) * j)  # from 1
  bd <- sort(unique(c(0, 1, cand[cand > 0 & cand < 1])))
  wide <- which(diff(bd) > 1 / 8)
  if (length(wide)) {
    extra <- unlist(lapply(wide, function(i) {
      m <- ceiling(diff(bd)[i] * 8)
      bd[i] + diff(bd)[i] * seq_len(m - 1) / m
    }))
    bd <- sort(unique(c(bd, extra)))
  }
  bd[c(TRUE, diff(bd) > 1e-15)]
}

.MaxEnt2Exponent <- function(fit, t) {
  # E(t) - max over [0, 1] of E, evaluated from the stored form.
  switch(fit$form,
         lin = fit$l1 * (t - as.numeric(fit$l1 > 0)),
         vertex = {
           v <- fit$p1
           c <- if (v <= 0) 0 else if (v >= 1) 1 else v
           fit$l2 * (t - c) * (t + c - 2 * v)
         },
         layers = {
           q <- fit$p1
           if (q >= 0) fit$l2 * t * (t - 1) - q * t
           else (t - 1) * (fit$l2 * t - q)
         })
}

.MaxEnt2Make <- function(form, p1, p2 = NULL) {
  # Builds a fit object: the form and its parameters, the panel
  # boundaries, the normalizer, the cumulative mass at every boundary,
  # and the first four central moments.
  L <- .MaxEnt2Params(form, p1, p2)
  l1 <- L[["l1"]]
  l2 <- L[["l2"]]
  fit <- list(form = form, p1 = p1, p2 = p2, l1 = l1, l2 = l2)
  bd <- .MaxEnt2Bounds(form, p1, p2, l1, l2)
  np <- length(bd) - 1L
  h <- diff(bd)
  t <- outer(h, .GL16.X) + bd[-(np + 1L)]            # np x 16 nodes
  w <- outer(h, .GL16.W)
  wg <- w * exp(.MaxEnt2Exponent(fit, t))
  Z <- sum(wg)
  m1 <- sum(wg * t) / Z
  d <- t - m1
  fit$bd <- bd
  fit$Z <- Z
  fit$cum <- c(0, cumsum(rowSums(wg)) / Z)
  fit$m1 <- m1
  fit$c2 <- sum(wg * d^2) / Z
  fit$c3 <- sum(wg * d^3) / Z
  fit$c4 <- sum(wg * d^4) / Z
  fit
}

.MaxEnt2Fit <- function(l1, l2) {
  # A fit from (l1, l2), choosing the form by the sign of l2.
  if (l2 == 0) .MaxEnt2Make("lin", l1)
  else if (l2 < 0) .MaxEnt2Make("vertex", -l1 / (2 * l2), l2)
  else .MaxEnt2Make("layers", -(l1 + l2), l2)
}

.MaxEnt2CdfRaw <- function(fit, u) {
  out <- numeric(length(u))
  out[u >= 1] <- 1
  mid <- which(u > 0 & u < 1)
  if (length(mid)) {
    p <- findInterval(u[mid], fit$bd, rightmost.closed = TRUE)
    # Chunked so the outer() below stays a modest matrix for long u.
    for (st in seq(1L, length(mid), by = 100000L)) {
      ii <- st:min(st + 99999L, length(mid))
      uu <- u[mid[ii]]
      pp <- p[ii]
      lo <- fit$bd[pp]
      hh <- uu - lo
      tt <- outer(hh, .GL16.X) + lo
      gg <- exp(.MaxEnt2Exponent(fit, tt))
      out[mid[ii]] <- fit$cum[pp] + hh * as.numeric(gg %*% .GL16.W) / fit$Z
    }
  }
  pmin(pmax(out, 0), 1)
}

# A bucket whose mean is nearer the upper edge is solved in the mirrored
# coordinate u' = 1 - u, so that the mean is always resolved from the
# nearer edge in double precision, and carries fit$mirror = TRUE.
.MaxEnt2Cdf <- function(fit, u) {
  if (isTRUE(fit$mirror)) 1 - .MaxEnt2CdfRaw(fit, 1 - u)
  else .MaxEnt2CdfRaw(fit, u)
}

.MaxEnt2PanelPoints <- function(fit) {
  if (isTRUE(fit$mirror)) 1 - rev(fit$bd) else fit$bd
}

# Residual: mean error relative to alpha (the solver always has
# alpha <= 1/2, the distance to the nearer edge) and the log of the
# variance ratio.
.MaxEnt2Residual <- function(fit, alpha, tau) {
  c((fit$m1 - alpha) / alpha, log(fit$c2 / tau))
}

# Both tolerances are floored at the rounding limit of the exponent,
# about eps (|l1| + |l2|): a layer thinner than 1 / (|l1| + |l2|) of the
# bucket cannot be placed more precisely than that in double precision.
.MaxEnt2Floor <- function(fit) {
  4 * .Machine$double.eps * (abs(fit$l1) + abs(fit$l2))
}

.MaxEnt2TolVar <- function(fit) max(.ME.TOL.VAR, .MaxEnt2Floor(fit))

.MaxEnt2Converged <- function(fit, alpha, tau) {
  r <- .MaxEnt2Residual(fit, alpha, tau)
  all(is.finite(r)) &&
    abs(r[1]) < max(.ME.TOL.MEAN, .MaxEnt2Floor(fit)) &&
    abs(r[2]) < .MaxEnt2TolVar(fit)
}

.MaxEnt2Jacobian <- function(fit, alpha) {
  # Jacobian of (relative mean, log var) with respect to (l1, l2), from
  # the central moments: d mean = (c2, c3 + 2 m1 c2),
  # d var = (c3, c4 - c2^2 + 2 m1 c3).
  matrix(c(fit$c2 / alpha, fit$c3 / fit$c2,
           (fit$c3 + 2 * fit$m1 * fit$c2) / alpha,
           (fit$c4 - fit$c2^2 + 2 * fit$m1 * fit$c3) / fit$c2), 2, 2)
}

.Solve2x2 <- function(J, r) {
  det <- J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]
  if (!is.finite(det) || det == 0) return(NULL)
  s <- c(J[2, 2] * r[1] - J[1, 2] * r[2], -J[2, 1] * r[1] + J[1, 1] * r[2]) / det
  if (any(!is.finite(s))) NULL else s
}

.MaxEnt2Newton <- function(fit, alpha, tau, max.iter = 40L) {
  # Damped Newton in (l1, l2): at most a 100% change of either per step.
  for (iter in seq_len(max.iter)) {
    if (.MaxEnt2Converged(fit, alpha, tau)) return(list(fit = fit, ok = TRUE))
    r <- .MaxEnt2Residual(fit, alpha, tau)
    if (!all(is.finite(r))) break
    step <- .Solve2x2(.MaxEnt2Jacobian(fit, alpha), r)
    if (is.null(step)) break
    L <- c(fit$l1, fit$l2)
    fac <- min(1, 1 / max(abs(step) / pmax(1, abs(L))))
    fit <- .MaxEnt2Fit(L[1] - fac * step[1], L[2] - fac * step[2])
  }
  list(fit = fit, ok = FALSE)
}

.MaxEnt2NewtonNatural <- function(fit, alpha, tau, max.iter = 60L) {
  # Damped Newton in the form's own parameters (p1, theta = log |l2|),
  # where (l1, l2) is ill conditioned: vertex (v, theta) with
  # d(l1, l2)/d(v, theta) = [[-2 l2, l1], [0, l2]], layers (q, theta)
  # with [[-1, -l2], [0, l2]].  Steps capped at 0.25 in v, 2 in q, 1 in
  # theta.
  if (fit$form == "lin") {
    return(list(fit = fit, ok = .MaxEnt2Converged(fit, alpha, tau)))
  }
  cap1 <- if (fit$form == "vertex") 0.25 else 2
  for (iter in seq_len(max.iter)) {
    if (.MaxEnt2Converged(fit, alpha, tau)) return(list(fit = fit, ok = TRUE))
    r <- .MaxEnt2Residual(fit, alpha, tau)
    if (!all(is.finite(r))) break
    D <- if (fit$form == "vertex") {
      matrix(c(-2 * fit$l2, 0, fit$l1, fit$l2), 2, 2)
    } else {
      matrix(c(-1, 0, -fit$l2, fit$l2), 2, 2)
    }
    step <- .Solve2x2(.MaxEnt2Jacobian(fit, alpha) %*% D, r)
    if (is.null(step)) break
    fac <- min(1, cap1 / abs(step[1]), 1 / abs(step[2]))
    p1 <- fit$p1 - fac * step[1]
    l2 <- sign(fit$l2) * exp(log(abs(fit$l2)) - fac * step[2])
    fit <- .MaxEnt2Make(fit$form, p1, l2)
  }
  list(fit = fit, ok = FALSE)
}

.MaxEnt2LayerStart <- function(alpha, tau) {
  # Start for the layers regime: two exponential layers, rate r0 = l2 + q
  # at u = 0 and r1 = l2 - q at u = 1, the layer at 1 holding mass p with
  # p / (1 - p) = exp(-q) r0 / r1.  Mixture mean and variance are closed
  # forms in (q, l2); the best point of a grid (scored on the mean, the
  # variance and its deficit from alpha (1 - alpha)) seeds Newton.
  g <- expand.grid(q = seq(-60, 60, by = 0.25), L = seq(-2, 60, by = 0.1))
  l2 <- exp(g$L)
  q <- g$q
  r0 <- l2 + q
  r1 <- l2 - q
  p <- 1 / (1 + exp(q) * r1 / r0)
  m <- (1 - p) / r0 + p * (1 - 1 / r1)
  s2 <- (1 - p) / r0^2 + p / r1^2 + p * (1 - p) * (1 - 1 / r0 - 1 / r1)^2
  s2max <- alpha * (1 - alpha)
  ok <- is.finite(m) & is.finite(s2) & r0 > 2 & r1 > 2 & s2 > 0 & s2 < s2max
  sc <- suppressWarnings(abs(m / alpha - 1) + abs(log(s2 / tau)) +
                         abs(log((s2max - s2) / (s2max - tau))))
  sc[!ok] <- Inf
  i <- which.min(sc)
  if (!is.finite(sc[i])) return(NULL)
  .MaxEnt2Make("layers", q[i], l2[i])
}

.MaxEnt2Solve <- function(alpha, tau) {
  # Solves (l1, l2) for mean alpha <= 1/2 and variance
  # 0 < tau < alpha (1 - alpha).
  #
  # The regime is fixed by tau against the one-moment variance v0: below
  # it the solution is a vertex form, above it a layers form.  Newton in
  # (l1, l2) from the Gaussian start or the one-moment solution; then
  # Newton in the form's own parameters from the Gaussian start or the
  # two-layer closed-form start; last, a coarse grid on the exact
  # moments.
  #
  # Args:
  #   alpha: Bucket mean rescaled to [0, 1], at most 1/2.
  #   tau:   Bucket variance rescaled to [0, 1], strictly between 0 and
  #          alpha (1 - alpha).
  #
  # Returns:
  #   A fit object (see .MaxEnt2Make).  Warns if the fit did not reach
  #   the solver tolerances.
  f0 <- .MaxEnt2Fit(.LambdaOfAlpha(alpha), 0)
  v0 <- f0$c2
  vertex <- tau < v0
  gauss <- if (vertex) .MaxEnt2Make("vertex", alpha, -1 / (2 * tau)) else NULL
  nw <- .MaxEnt2Newton(if (vertex) gauss else f0, alpha, tau)
  if (nw$ok && (nw$fit$form != "lin")) return(nw$fit)
  st <- if (vertex) gauss else .MaxEnt2LayerStart(alpha, tau)
  if (!is.null(st)) {
    nw <- .MaxEnt2NewtonNatural(st, alpha, tau)
    if (nw$ok) return(nw$fit)
  }
  best <- NULL
  best.r <- Inf
  p1.grid <- if (vertex) seq(0, 0.5, by = 0.025) else seq(-60, 60, by = 2)
  for (p1 in p1.grid) {
    for (th in seq(-3, 60, by = 1)) {
      f <- .MaxEnt2Make(if (vertex) "vertex" else "layers", p1,
                        if (vertex) -exp(th) else exp(th))
      r <- .MaxEnt2Residual(f, alpha, tau)
      if (!all(is.finite(r))) next
      sc <- abs(r[1]) + abs(r[2])
      if (sc < best.r) {
        best.r <- sc
        best <- f
      }
    }
  }
  if (!is.null(best)) nw <- .MaxEnt2NewtonNatural(best, alpha, tau)
  if (!nw$ok) {
    r <- .MaxEnt2Residual(nw$fit, alpha, tau)
    warning(sprintf("max-entropy fit did not converge: alpha=%g tau=%g resid=%g,%g",
                    alpha, tau, r[1], r[2]))
  }
  nw$fit
}

# ---- Per-bucket description shared by the one- and two-moment cases -----
#
# A description is a list with kind "atom" (fields loc and wt, the atom
# locations in the bucket's own x units and their weights summing to 1)
# or kind "cont" with either l (one moment) or fit (two moments).

.PolyBucketDescription <- function(h, j) {
  # Describes bucket j of a polynomial histogram for reconstruction.
  #
  # Args:
  #   h: A polynomial.histogram object.
  #   j: A bucket index.
  #
  # Returns:
  #   A bucket description as above.
  a <- h$breaks[j]
  b <- h$breaks[j + 1L]
  two <- !is.null(h$variances)
  if (h$counts[j] == 0 || is.na(h$means[j])) {
    # An empty bucket is filled uniformly; it carries no mass anyway.
    return(if (two) list(kind = "cont", fit = .MaxEnt2Fit(0, 0))
           else list(kind = "cont", l = 0))
  }
  mu <- h$means[j]
  alpha <- (mu - a) / (b - a)
  beta <- (b - mu) / (b - a)
  # A mean at an edge is a point mass there (an atom sits at the stored
  # mean itself; a log-spaced break can be an ulp off the data value).
  if (alpha <= .ME.DEGENERATE.TOL || beta <= .ME.DEGENERATE.TOL) {
    return(list(kind = "atom", loc = mu, wt = 1))
  }
  if (!two) {
    return(list(kind = "cont",
                l = if (alpha <= beta) .LambdaOfAlpha(alpha)
                    else -.LambdaOfAlpha(beta)))
  }
  tau <- max(h$variances[j], 0) / (b - a)^2
  s2max <- alpha * beta
  if (tau <= .ME.DEGENERATE.TOL * s2max) {
    return(list(kind = "atom", loc = mu, wt = 1))
  }
  if (tau >= s2max * (1 - .ME.DEGENERATE.TOL)) {
    return(list(kind = "atom", loc = c(a, b), wt = c(beta, alpha)))
  }
  # Solve from the nearer edge (see .MaxEnt2Cdf).
  if (alpha <= beta) {
    fit <- .MaxEnt2Solve(alpha, tau)
  } else {
    fit <- .MaxEnt2Solve(beta, tau)
    fit$mirror <- TRUE
  }
  list(kind = "cont", fit = fit)
}

.PolyBucketCdf <- function(d, u, a, b) {
  # Within-bucket CDF at the rescaled positions u in [0, 1].
  #
  # Args:
  #   d:    A bucket description.
  #   u:    Numeric vector of positions rescaled to the bucket.
  #   a, b: The bucket's edges, needed to place atoms.
  #
  # Returns:
  #   The fraction of the bucket's mass at or below each u.
  if (d$kind == "atom") {
    f <- numeric(length(u))
    loc.u <- (d$loc - a) / (b - a)
    for (i in seq_along(d$loc)) {
      f <- f + d$wt[i] * (u >= loc.u[i] - .ME.DEGENERATE.TOL)
    }
    return(pmin(f, 1))
  }
  if (!is.null(d$fit)) return(.MaxEnt2Cdf(d$fit, u))
  .MaxEnt1Cdf(u, d$l)
}

.PolyBucketQuantile <- function(d, f, a, b) {
  # Within-bucket quantile at mass fraction f, in the bucket's x units.
  if (d$kind == "atom") {
    cw <- cumsum(d$wt)
    return(d$loc[which(cw >= f - 1e-12)[1]])
  }
  # The edges are returned directly: at f = 1 the closed-form inverse of
  # a steep exponential (|l| > 700) overflows, and uniroot needs a sign
  # change.
  if (f <= 0) return(a)
  if (f >= 1) return(b)
  if (!is.null(d$fit)) {
    u <- uniroot(function(u) .MaxEnt2Cdf(d$fit, u) - f, c(0, 1),
                 tol = 1e-12)$root
    return(a + u * (b - a))
  }
  a + .MaxEnt1Quantile(f, d$l) * (b - a)
}
