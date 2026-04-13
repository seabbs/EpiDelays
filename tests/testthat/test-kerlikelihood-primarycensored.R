# Equivalence tests: kerlikelihood() vs primarycensored.
#
# Both the numerical-integration (`ni`) branch of kerlikelihood() and
# primarycensored::pprimarycensored() compute the same doubly-interval-censored
# likelihood. These tests pin that equivalence down numerically so we can
# later swap the internals without changing semantics.
#
# Derivation:
#
#   kerlikelihood inner integral (per row i):
#     (1 / (x1r - x1l)) * int_{x1l}^{x1r} [F(x2r - t1) - F(x2l - t1)] dt1
#
# Substitute s = t1 - x1l, so s ~ Uniform(0, pwindow) with pwindow = x1r - x1l:
#     = E_s[ F((x2r - x1l) - s) - F((x2l - x1l) - s) ]
#
# pprimarycensored(q, pwindow) computes E_s[F(q - s)], therefore the inner
# integral equals
#     pprimarycensored(x2r - x1l, pwindow) - pprimarycensored(x2l - x1l, pwindow).
#
# And by definition dprimarycensored(x2l - x1l, pwindow, swindow = x2r - x2l)
# returns the same difference.

make_double_data <- function(seed = 1L, n = 5L) {
  set.seed(seed)
  x1l <- sort(stats::runif(n, 0, 3))
  x1r <- x1l + stats::runif(n, 0.5, 1.5)
  x2l <- x1r + stats::runif(n, 0.2, 4)
  x2r <- x2l + stats::runif(n, 0.3, 1.5)
  data.frame(x1l = x1l, x1r = x1r, x2l = x2l, x2r = x2r)
}

kerlik_inner_via_primarycensored <- function(x, pdist, ...) {
  # Per-row contribution to the kerlikelihood() log-likelihood, computed via
  # primarycensored. Each element equals exp(per-row log-lik contribution):
  #   F_cens(x2r - x1l) - F_cens(x2l - x1l)
  # where F_cens is pprimarycensored with pwindow = x1r - x1l. The
  # normalisation by pwindow is built into pprimarycensored, which averages
  # F(q - s) with s ~ Uniform(0, pwindow).
  dots <- list(...)
  vapply(seq_len(nrow(x)), function(i) {
    pwindow <- x$x1r[i] - x$x1l[i]
    qr <- x$x2r[i] - x$x1l[i]
    ql <- x$x2l[i] - x$x1l[i]
    args <- c(list(q = c(ql, qr), pdist = pdist, pwindow = pwindow), dots)
    cdfs <- do.call(primarycensored::pprimarycensored, args)
    cdfs[2] - cdfs[1]
  }, numeric(1))
}

test_that("gaussian double-interval ni matches pprimarycensored", {
  skip_if_no_primarycensored()
  x <- make_double_data()
  v <- c(1.5, log(0.8))  # unbounded-space params: mean=1.5, sd=0.8

  m <- kerlikelihood(x = x, family = "gaussian", likapprox = "ni")
  ker_value <- m$loglik(v, x)

  inner <- kerlik_inner_via_primarycensored(
    x, pdist = stats::pnorm, mean = v[1], sd = exp(v[2])
  )
  expected <- sum(log(inner))

  expect_equal(ker_value, expected, tolerance = 1e-8)
})
