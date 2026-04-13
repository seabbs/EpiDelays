# Equivalence tests: kerlikelihood() vs primarycensored.
#
# kerlikelihood()'s doubly-interval-censored ni branch is implemented on top
# of primarycensored::dprimarycensored(). These tests pin that equivalence
# down numerically and lock the current implementation against a
# reconstruction of the previous stats::integrate() based algorithm so that
# future refactors cannot drift silently.
#
# Derivation:
#
#   kerlikelihood inner integral (per row i):
#     (1 / (x1r - x1l)) * int_{x1l}^{x1r} [F(x2r - t1) - F(x2l - t1)] dt1
#
# Substitute s = t1 - x1l, so s ~ Uniform(0, pwindow) with pwindow = x1r - x1l:
#     = E_s[ F((x2r - x1l) - s) - F((x2l - x1l) - s) ]
#
# By definition dprimarycensored(x2l - x1l, pwindow, swindow = x2r - x2l)
# returns the log of that difference, so the total log-likelihood is a single
# vectorised call to dprimarycensored with log = TRUE.

make_double_data <- function(seed = 1L, n = 5L) {
  set.seed(seed)
  x1l <- sort(stats::runif(n, 0, 3))
  x1r <- x1l + stats::runif(n, 0.5, 1.5)
  x2l <- x1r + stats::runif(n, 0.2, 4)
  x2r <- x2l + stats::runif(n, 0.3, 1.5)
  data.frame(x1l = x1l, x1r = x1r, x2l = x2l, x2r = x2r)
}

# Idiomatic primarycensored oracle: sum the per-row log-density from
# dprimarycensored() directly. primarycensored::dprimarycensored requires
# scalar pwindow/swindow, so rows are iterated one at a time — the call
# shape still mirrors the one end users should reach for when writing a
# doubly-interval-censored likelihood by hand.
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

# Reconstruction of the pre-primarycensored stats::integrate() inner loop for
# gaussian doubly interval-censored data. Kept test-local so the old
# algorithm lives on in exactly one place as a regression lock; if the
# production dprimarycensored path ever drifts from the numerical integral,
# this test will fail.
kerlik_integrate_reference_gaussian <- function(v, x) {
  mean_par <- v[1]
  sd_par <- exp(v[2])
  n <- nrow(x)
  z <- 0
  for (i in seq_len(n)) {
    h <- function(t1) {
      stats::pnorm(x$x2r[i] - t1, mean = mean_par, sd = sd_par) -
        stats::pnorm(x$x2l[i] - t1, mean = mean_par, sd = sd_par)
    }
    logint <- log(
      stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
    )
    z <- z + logint - log(x$x1r[i] - x$x1l[i])
  }
  z
}

# CDF wrapper for pskewnorm with a signature compatible with
# primarycensored::dprimarycensored (which calls pdist(q, ...)).
pskewnorm_cdf <- function(q, location, scale, slant) {
  EpiDelays::pskewnorm(x = q, par1 = location, par2 = scale, par3 = slant)
}

# Parametric family configurations. Each entry describes the family-specific
# bits needed to (a) evaluate kerlikelihood and (b) call dprimarycensored with
# the matching CDF and parameters.
family_cases <- list(
  gaussian = list(
    v = c(1.5, log(0.8)),
    pdist = stats::pnorm,
    pars = list(mean = 1.5, sd = 0.8)
  ),
  gamma = list(
    v = c(log(2), log(0.5)),
    pdist = stats::pgamma,
    pars = list(shape = 2, rate = 0.5)
  ),
  lognormal = list(
    v = c(0.5, log(0.4)),
    pdist = stats::plnorm,
    pars = list(meanlog = 0.5, sdlog = 0.4)
  ),
  weibull = list(
    v = c(log(2), log(2)),
    pdist = stats::pweibull,
    pars = list(shape = 2, scale = 2)
  ),
  skewnorm = list(
    v = c(1, log(1), 2),
    pdist = pskewnorm_cdf,
    pars = list(location = 1, scale = 1, slant = 2)
  )
)

for (fam in names(family_cases)) {
  local({
    family <- fam
    case <- family_cases[[family]]

    test_that(sprintf(
      "%s double-interval ni matches dprimarycensored oracle", family
    ), {
      skip_if_no_primarycensored()
      x <- make_double_data()

      m <- kerlikelihood(x = x, family = family, likapprox = "ni")
      ker_value <- m$loglik(case$v, x)

      expected <- kerlik_loglik_via_dprimarycensored(
        x = x, pdist = case$pdist, pars = case$pars
      )

      expect_equal(ker_value, expected, tolerance = 1e-8)
    })

    test_that(sprintf(
      "%s double-interval mc matches dprimarycensored oracle", family
    ), {
      skip_if_no_primarycensored()
      x <- make_double_data()

      expected <- kerlik_loglik_via_dprimarycensored(
        x = x, pdist = case$pdist, pars = case$pars
      )

      # kerlikelihood's mc closure uses stats::runif internally; seeding the
      # global RNG immediately before the call pins the Monte Carlo draws and
      # keeps the test reproducible. With M = 1000 samples per row and 5 rows
      # the summed log-likelihood sits well within 5e-2 of the ni target for
      # all families tested here.
      m <- kerlikelihood(x = x, family = family, likapprox = "mc")
      set.seed(20260413L)
      ker_value <- m$loglik(case$v, x)

      expect_equal(ker_value, expected, tolerance = 5e-2)
    })

    test_that(sprintf(
      "%s single-interval matches F(xr) - F(xl)", family
    ), {
      x_double <- make_double_data()
      # Reuse the existing interval for the secondary event as a single
      # interval-censored observation (xl, xr) — positive throughout so the
      # gamma/lognormal/weibull supports are respected.
      x <- data.frame(xl = x_double$x2l, xr = x_double$x2r)

      m <- kerlikelihood(x = x, family = family, likapprox = "ni")
      ker_value <- m$loglik(case$v, x)

      Fl <- do.call(case$pdist, c(list(x$xl), case$pars))
      Fr <- do.call(case$pdist, c(list(x$xr), case$pars))
      expected <- sum(log(Fr - Fl))

      expect_equal(ker_value, expected, tolerance = 1e-10)
    })
  })
}

test_that("gaussian ni matches the legacy integrate-based reference", {
  # Regression lock: the production dprimarycensored path must agree with a
  # direct reconstruction of the removed stats::integrate() inner loop. This
  # keeps the previous algorithm reproducible from one test-local helper and
  # flags any drift in a future refactor. Gaussian is used because it has no
  # closed-form pcens_cdf method inside primarycensored, so the agreement is
  # strictly numerical.
  skip_if_no_primarycensored()
  x <- make_double_data()
  v <- c(1.5, log(0.8))

  m <- kerlikelihood(x = x, family = "gaussian", likapprox = "ni")
  expect_equal(
    m$loglik(v, x),
    kerlik_integrate_reference_gaussian(v, x),
    tolerance = 1e-8
  )
})

test_that("parfitml fits gamma end-to-end on doubly censored data", {
  set.seed(42L)
  n <- 30L
  # Simulate doubly interval-censored gamma data. Use narrow primary windows
  # and small secondary windows so the fit converges quickly with few
  # bootstrap iterations.
  x1l <- sort(stats::runif(n, 0, 2))
  x1r <- x1l + stats::runif(n, 0.2, 0.8)
  true_draw <- stats::rgamma(n, shape = 3, rate = 1)
  x2l <- x1l + true_draw
  x2r <- x2l + stats::runif(n, 0.2, 0.8)
  x <- data.frame(x1l = x1l, x1r = x1r, x2l = x2l, x2r = x2r)

  fit <- parfitml(x = x, family = "gamma", Bboot = 10L, pgbar = FALSE)
  expect_true(fit$mleconv)
  expect_equal(fit$censtype, "double")
})

test_that("kerlikelihood handles mixed pwindow rows", {
  skip_if_no_primarycensored()
  # Mix two distinct primary windows (1 and 2) so the group-by logic in
  # build_pc_loglik() has to dispatch more than one dprimarycensored call.
  x <- data.frame(
    x1l = c(0, 1, 2, 3, 4, 5),
    x1r = c(1, 2, 3, 5, 6, 7),
    x2l = c(3, 4, 5, 6, 7, 8),
    x2r = c(4, 5, 6, 7, 8, 9)
  )
  stopifnot(length(unique(x$x1r - x$x1l)) == 2L)

  v <- c(1.5, log(0.8))
  m <- kerlikelihood(x = x, family = "gaussian", likapprox = "ni")

  expected <- kerlik_loglik_via_dprimarycensored(
    x = x, pdist = stats::pnorm,
    pars = list(mean = v[1], sd = exp(v[2]))
  )
  expect_equal(m$loglik(v, x), expected, tolerance = 1e-8)
})
