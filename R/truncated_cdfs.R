#' Zero-truncated CDF wrappers for signed-support delay families
#'
#' Internal helpers used by \code{kerlikelihood()} to interpret the gaussian
#' and skew-normal families as non-negative delay distributions. Each wrapper
#' returns a proper CDF on the non-negative half-line by clipping at zero and
#' renormalising the mass to the right of zero. They are passed into
#' \code{primarycensored::dprimarycensored()} as \code{pdist}, so the
#' primary-censoring integral is evaluated on the truncated underlying delay
#' distribution rather than on the full real line.
#'
#' The wrappers are CDFs in the formal sense: monotone non-decreasing, right-
#' continuous, 0 at \code{q = 0} and 1 at \code{q = Inf}.
#'
#' @name truncated-cdfs
#' @keywords internal
NULL

# Clamp to [0, 1] and enforce monotone non-decreasing output on sorted inputs.
# pskewnorm uses Owen's T, which can produce values a hair outside [0, 1] in
# extreme parameter regions, and rounding noise around the tails can leave
# the cummax a hair non-monotone. This is reused by both the raw skew-normal
# CDF wrapper inside kerlikelihood() and the truncated wrapper below so the
# two paths agree on the same tail handling.
clamp_monotone_cdf <- function(q, v) {
  v <- pmin(1, pmax(0, v))
  ord <- order(q)
  v[ord] <- cummax(v[ord])
  v
}

#' @rdname truncated-cdfs
#' @keywords internal
ptruncnorm_nonneg <- function(q, mean, sd) {
  f0 <- stats::pnorm(0, mean = mean, sd = sd)
  # Deep negative means push f0 towards 1; the renormaliser then blows up.
  # Return 0 in that regime so downstream dprimarycensored sees a valid
  # (though degenerate) CDF rather than NaN.
  denom <- 1 - f0
  out <- rep(0, length(q))
  pos <- q > 0
  if (any(pos) && denom > 0) {
    out[pos] <- (stats::pnorm(q[pos], mean = mean, sd = sd) - f0) / denom
  }
  out[is.infinite(q) & q > 0] <- 1
  clamp_monotone_cdf(q, out)
}

#' @rdname truncated-cdfs
#' @keywords internal
ptruncskewnorm_nonneg <- function(q, location, scale, slant) {
  f0 <- pskewnorm(x = 0, par1 = location, par2 = scale, par3 = slant)
  denom <- 1 - f0
  out <- rep(0, length(q))
  pos <- q > 0
  if (any(pos) && denom > 0) {
    raw <- pskewnorm(x = q[pos], par1 = location, par2 = scale, par3 = slant)
    out[pos] <- (raw - f0) / denom
  }
  out[is.infinite(q) & q > 0] <- 1
  clamp_monotone_cdf(q, out)
}
