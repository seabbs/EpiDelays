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

# CDF wrapper for pskewnorm with a signature compatible with
# primarycensored::pprimarycensored (which calls pdist(q, ...)).
pskewnorm_cdf <- function(q, location, scale, slant) {
  EpiDelays::pskewnorm(x = q, par1 = location, par2 = scale, par3 = slant)
}

# Parametric family configurations. Each entry describes the family-specific
# bits needed to (a) evaluate kerlikelihood and (b) call pprimarycensored with
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

    test_that(sprintf("%s double-interval ni matches pprimarycensored", family), {
      skip_if_no_primarycensored()
      x <- make_double_data()

      m <- kerlikelihood(x = x, family = family, likapprox = "ni")
      ker_value <- m$loglik(case$v, x)

      inner <- do.call(
        kerlik_inner_via_primarycensored,
        c(list(x = x, pdist = case$pdist), case$pars)
      )
      expected <- sum(log(inner))

      expect_equal(ker_value, expected, tolerance = 1e-8)
    })

    test_that(sprintf("%s double-interval mc matches pprimarycensored", family), {
      skip_if_no_primarycensored()
      x <- make_double_data()

      inner <- do.call(
        kerlik_inner_via_primarycensored,
        c(list(x = x, pdist = case$pdist), case$pars)
      )
      expected <- sum(log(inner))

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

    test_that(sprintf("%s single-interval matches F(xr) - F(xl)", family), {
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

for (fam in names(family_cases)) {
  local({
    family <- fam
    case <- family_cases[[family]]

    test_that(sprintf(
      "%s kerlikelihood engine primarycensored matches integrate", family
    ), {
      skip_if_no_primarycensored()
      x <- make_double_data()

      m_int <- kerlikelihood(
        x = x, family = family, likapprox = "ni", engine = "integrate"
      )
      m_pc <- kerlikelihood(
        x = x, family = family, likapprox = "ni", engine = "primarycensored"
      )

      expect_equal(
        m_pc$loglik(case$v, x),
        m_int$loglik(case$v, x),
        tolerance = 1e-8
      )
    })
  })
}

test_that("kerlikelihood engine defaults to integrate", {
  x <- make_double_data()
  m_default <- kerlikelihood(x = x, family = "gaussian", likapprox = "ni")
  m_int <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "ni", engine = "integrate"
  )
  v <- c(1.5, log(0.8))
  expect_equal(m_default$loglik(v, x), m_int$loglik(v, x))
})

test_that("kerlikelihood engine silently falls back for single interval", {
  # nc == 2: primarycensored engine must be ignored (no primary window).
  x <- data.frame(xl = c(1, 2, 3), xr = c(2, 3, 4))
  v <- c(1.5, log(0.8))
  m_int <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "ni", engine = "integrate"
  )
  m_pc <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "ni", engine = "primarycensored"
  )
  expect_equal(m_pc$loglik(v, x), m_int$loglik(v, x))
})

test_that("kerlikelihood engine silently falls back for mc approximation", {
  # mc branch is not wired to primarycensored; the engine argument should be
  # accepted without error and produce the same result as the integrate path.
  x <- make_double_data()
  v <- c(1.5, log(0.8))
  m_int <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "mc", engine = "integrate"
  )
  m_pc <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "mc", engine = "primarycensored"
  )
  set.seed(20260413L)
  a <- m_int$loglik(v, x)
  set.seed(20260413L)
  b <- m_pc$loglik(v, x)
  expect_equal(a, b)
})

test_that("kerlikelihood engine argument is validated via match.arg", {
  x <- make_double_data()
  expect_error(
    kerlikelihood(
      x = x, family = "gaussian", likapprox = "ni", engine = "nonsense"
    )
  )
})

test_that("parfitml default engine still fits gamma end-to-end", {
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

test_that("primarycensored engine handles mixed pwindow rows", {
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
  m_int <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "ni", engine = "integrate"
  )
  m_pc <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "ni", engine = "primarycensored"
  )
  expect_equal(m_pc$loglik(v, x), m_int$loglik(v, x), tolerance = 1e-8)
})

test_that("gaussian dprimarycensored matches per-row kerlikelihood contribution", {
  skip_if_no_primarycensored()
  x <- make_double_data()
  v <- c(1.5, log(0.8))
  mean_par <- v[1]
  sd_par <- exp(v[2])

  # Per-row contribution from kerlikelihood: recompute the inner integral by
  # hand so we have a per-row vector (m$loglik only returns the sum).
  n <- nrow(x)
  per_row_kerlik <- vapply(seq_len(n), function(i) {
    h <- function(t1) {
      stats::pnorm(x$x2r[i] - t1, mean = mean_par, sd = sd_par) -
        stats::pnorm(x$x2l[i] - t1, mean = mean_par, sd = sd_par)
    }
    integral <- stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
    integral / (x$x1r[i] - x$x1l[i])
  }, numeric(1))

  per_row_dprim <- vapply(seq_len(n), function(i) {
    primarycensored::dprimarycensored(
      x = x$x2l[i] - x$x1l[i],
      pdist = stats::pnorm,
      pwindow = x$x1r[i] - x$x1l[i],
      swindow = x$x2r[i] - x$x2l[i],
      mean = mean_par,
      sd = sd_par
    )
  }, numeric(1))

  expect_equal(per_row_dprim, per_row_kerlik, tolerance = 1e-8)
})
