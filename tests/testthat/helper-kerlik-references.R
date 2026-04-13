# Test-local reference oracles for kerlikelihood() equivalence tests.
#
# Two flavours are provided:
#
# - kerlik_integrate_reference_truncated() takes whatever pdist / pars are
#   passed in and reconstructs the doubly-interval-censored inner integral
#   via stats::integrate(). For the signed-support families (gaussian,
#   skewnorm) the caller is expected to pass the zero-truncated wrapper so
#   this oracle matches the production path in kerlikelihood(). For the
#   positive-support families (gamma, lognormal, weibull) the raw CDF is
#   identical to the truncated CDF because F(q) = 0 for q < 0, and both
#   flavours agree to numerical precision.
#
# - kerlik_integrate_reference_upstream() is a verbatim copy of the
#   pre-primarycensored-swap integrate() inner loop. It applies NO
#   truncation and NO renormalisation to the supplied pdist. This is what
#   upstream oswaldogressani/EpiDelays kerlikelihood() computed before the
#   dprimarycensored swap, and is kept as a permanent anchor for
#   reproducing upstream's numbers when debugging.
#
# Both oracles are CDF-agnostic: the same family_cases entry drives them,
# keyed on whatever pdist is passed in.

# nolint start: object_length_linter, line_length_linter.
kerlik_integrate_reference_truncated <- function(x, pdist, pars) {
  n <- nrow(x)
  z <- 0
  for (i in seq_len(n)) {
    h <- function(t1) {
      do.call(pdist, c(list(x$x2r[i] - t1), pars)) -
        do.call(pdist, c(list(x$x2l[i] - t1), pars))
    }
    logint <- log(
      stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
    )
    z <- z + logint - log(x$x1r[i] - x$x1l[i])
  }
  z
}

# Permanent anchor for upstream's pre-swap shape. Do not add truncation or
# renormalisation here — this is the shape of what upstream
# oswaldogressani/EpiDelays kerlikelihood() computes before the
# dprimarycensored swap.
kerlik_integrate_reference_upstream <- function(x, pdist, pars) {
  n <- nrow(x)
  z <- 0
  for (i in seq_len(n)) {
    h <- function(t1) {
      do.call(pdist, c(list(x$x2r[i] - t1), pars)) -
        do.call(pdist, c(list(x$x2l[i] - t1), pars))
    }
    logint <- log(
      stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
    )
    z <- z + logint - log(x$x1r[i] - x$x1l[i])
  }
  z
}

# Idiomatic primarycensored oracle. CDF-agnostic: whatever pdist and pars
# are supplied feed straight into dprimarycensored().
kerlik_loglik_via_dprimarycensored <- function(x, pdist, pars) {
  sum(vapply(seq_len(nrow(x)), function(i) {
    do.call(
      primarycensored::dprimarycensored,
      c(
        list(
          x = x$x2l[i] - x$x1l[i],
          pdist = pdist,
          pwindow = x$x1r[i] - x$x1l[i],
          swindow = x$x2r[i] - x$x2l[i],
          log = TRUE
        ),
        pars
      )
    )
  }, numeric(1)))
}
# nolint end
