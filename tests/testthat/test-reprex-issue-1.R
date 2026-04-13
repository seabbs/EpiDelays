# Runnable version of the smoke-test reprex posted to
# https://github.com/oswaldogressani/EpiDelays/issues/1
#
# Confirms per-row that kerlikelihood()'s doubly-interval-censored contribution
# for the gaussian family matches primarycensored::dprimarycensored() with a
# uniform primary event window. Kept deliberately small and dependency-light so
# anyone can paste it into an R session and check for themselves.

test_that("reprex (issue #1): kerlikelihood == dprimarycensored per row", {
  skip_if_no_primarycensored()

  x <- data.frame(
    x1l = c(0.5, 1.0, 2.0),
    x1r = c(1.5, 2.5, 3.2),
    x2l = c(2.0, 3.5, 5.0),
    x2r = c(3.0, 4.5, 6.0)
  )
  mean_par <- 1.5
  sd_par <- 0.8

  per_row_ker <- vapply(seq_len(nrow(x)), function(i) {
    h <- function(t1) {
      stats::pnorm(x$x2r[i] - t1, mean_par, sd_par) -
        stats::pnorm(x$x2l[i] - t1, mean_par, sd_par)
    }
    stats::integrate(h, x$x1l[i], x$x1r[i])$value / (x$x1r[i] - x$x1l[i])
  }, numeric(1))

  per_row_pc <- vapply(seq_len(nrow(x)), function(i) {
    primarycensored::dprimarycensored(
      x = x$x2l[i] - x$x1l[i],
      pdist = stats::pnorm,
      pwindow = x$x1r[i] - x$x1l[i],
      swindow = x$x2r[i] - x$x2l[i],
      mean = mean_par,
      sd = sd_par
    )
  }, numeric(1))

  expect_equal(per_row_pc, per_row_ker, tolerance = 1e-8)
})
