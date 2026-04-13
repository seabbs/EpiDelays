# Runnable version of the smoke-test reprex posted to
# https://github.com/oswaldogressani/EpiDelays/issues/1
#
# Confirms that kerlikelihood()'s doubly-interval-censored log-likelihood for
# integer-day data agrees with a single vectorised call to
# primarycensored::dprimarycensored(). dprimarycensored already dedupes its
# inputs via an internal lookup table over `unique(c(x, x + swindow))`, so a
# correct client call is a one-liner — no hand-rolled match() wrapper needed.

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
  sum(do.call(
    primarycensored::dprimarycensored,
    c(
      list(
        x = x$x2l - x$x1l,
        pdist = pdist,
        pwindow = 1,
        swindow = 1,
        log = TRUE
      ),
      par_map(v)
    )
  ))
}

test_that("reprex (issue #1): gamma kerlikelihood == dprimarycensored", {
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

test_that("reprex (issue #1): gaussian kerlikelihood == dprimarycensored", {
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
