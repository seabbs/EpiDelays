# Runnable version of the smoke-test reprex posted to
# https://github.com/oswaldogressani/EpiDelays/issues/1
#
# Confirms equivalence between kerlikelihood() and a vectorised primarycensored
# call on integer-day doubly-interval-censored data, for both an analytical
# family (gamma) and a numerical-only family (gaussian). The dedup pattern for
# the gaussian case is the same one the issue comment uses to show primary-
# censored still wins when there is no closed form, because realistic
# integer-day data has few unique delay quantiles.

make_integer_day_data <- function(seed = 1L, n = 100L) {
  set.seed(seed)
  x1l <- sample(0:14, n, replace = TRUE)
  data.frame(
    x1l = x1l,
    x1r = x1l + 1L,
    x2l = x1l + sample(2:8, n, replace = TRUE),
    x2r = NA_integer_
  ) |>
    transform(x2r = x2l + 1L)
}

pc_loglik <- function(v, x, pdist, par_map) {
  n <- nrow(x)
  qs_all <- c(x$x2l - x$x1l, x$x2r - x$x1l)
  uq <- unique(qs_all)
  args <- c(list(q = uq, pdist = pdist, pwindow = 1), par_map(v))
  cdfs_u <- do.call(primarycensored::pprimarycensored, args)
  cdfs <- cdfs_u[match(qs_all, uq)]
  sum(log(cdfs[seq.int(n + 1L, 2L * n)] - cdfs[seq_len(n)]))
}

test_that("reprex (issue #1): gamma kerlikelihood == vectorised pprimarycensored", {
  skip_if_no_primarycensored()

  x <- make_integer_day_data()
  v <- c(log(3), log(1))  # shape = 3, rate = 1

  m <- kerlikelihood(x, family = "gamma", likapprox = "ni")
  expect_equal(
    pc_loglik(v, x, stats::pgamma,
              function(v) list(shape = exp(v[1]), rate = exp(v[2]))),
    m$loglik(v, x),
    tolerance = 1e-8
  )
})

test_that("reprex (issue #1): gaussian kerlikelihood == vectorised pprimarycensored", {
  skip_if_no_primarycensored()

  x <- make_integer_day_data()
  v <- c(2.5, log(1))  # mean = 2.5, sd = 1

  m <- kerlikelihood(x, family = "gaussian", likapprox = "ni")
  expect_equal(
    pc_loglik(v, x, stats::pnorm,
              function(v) list(mean = v[1], sd = exp(v[2]))),
    m$loglik(v, x),
    tolerance = 1e-8
  )
})
