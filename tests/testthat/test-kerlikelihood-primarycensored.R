# Equivalence tests: kerlikelihood() vs primarycensored.
#
# kerlikelihood()'s doubly-interval-censored ni branch is implemented on top
# of primarycensored::dprimarycensored(). These tests pin that equivalence
# down numerically and lock the current implementation against a direct
# integrate()-based reconstruction of the same inner integral, so that
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
  # Continuous primary onset calendar baseline plus row-varying primary and
  # secondary observation windows. Delays are drawn via
  # primarycensored::rprimarycensored() so the simulator demonstrates the
  # primarycensored idiom on the sampling side. rprimarycensored quantises
  # the total delay to a multiple of swindow, so the secondary observation
  # interval is [x1l + delay, x1l + delay + swidth]. The evaluation family in
  # the tests below is not required to match the sampling distribution — the
  # simulator just needs to produce plausible doubly-interval-censored rows.
  x1l <- sort(stats::runif(n, 0, 3))
  pwidth <- stats::runif(n, 0.5, 1.5)
  swidth <- stats::runif(n, 0.3, 1.5)
  delays <- vapply(seq_len(n), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rnorm,
      pwindow = pwidth[i],
      swindow = swidth[i],
      mean = 2.5,
      sd = 1
    )
  }, numeric(1))
  data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )
}

# Reference oracles (kerlik_loglik_via_dprimarycensored,
# kerlik_integrate_reference_truncated, kerlik_integrate_reference_upstream)
# are defined in tests/testthat/helper-kerlik-references.R. They are
# CDF-agnostic, so each family_cases entry supplies its own `pdist`. The
# signed-support families (gaussian, skewnorm) use the zero-truncated
# wrappers that kerlikelihood() itself passes into dprimarycensored, so the
# production path and the oracles share the same underlying delay CDF.
family_cases <- list(
  gaussian = list(
    v = c(1.5, log(0.8)),
    pdist = ptruncnorm_nonneg,
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
    pdist = ptruncskewnorm_nonneg,
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
      "%s ni matches the truncated integrate-based reference", family
    ), {
      # Regression lock: the production dprimarycensored path must agree
      # with the direct integrate()-based reconstruction in
      # kerlik_integrate_reference_truncated(), flagging any drift in a
      # future refactor. For the signed-support families the `pdist` in
      # family_cases is already a zero-truncated wrapper, so the oracle
      # and the production path agree on the semantics of the underlying
      # delay distribution.
      skip_if_no_primarycensored()
      x <- make_double_data()

      m <- kerlikelihood(x = x, family = family, likapprox = "ni")
      expect_equal(
        m$loglik(case$v, x),
        kerlik_integrate_reference_truncated(
          x = x, pdist = case$pdist, pars = case$pars
        ),
        tolerance = 1e-8
      )
    })

    test_that(sprintf(
      "%s single-interval matches F(xr) - F(xl)", family
    ), {
      skip_if_no_primarycensored()
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

test_that("parfitml fits gamma end-to-end on doubly censored data", {
  skip_if_no_primarycensored()
  set.seed(42L)
  n <- 30L
  # Simulate doubly interval-censored gamma data via
  # primarycensored::rprimarycensored(). Narrow primary and secondary windows
  # keep the fit quick with few bootstrap iterations. rprimarycensored needs
  # scalar pwindow/swindow, so draw row-by-row.
  x1l <- sort(stats::runif(n, 0, 2))
  pwidth <- stats::runif(n, 0.2, 0.8)
  swidth <- stats::runif(n, 0.2, 0.8)
  delays <- vapply(seq_len(n), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rgamma,
      pwindow = pwidth[i],
      swindow = swidth[i],
      shape = 3,
      rate = 1
    )
  }, numeric(1))
  x <- data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )

  fit <- parfitml(x = x, family = "gamma", Bboot = 10L, pgbar = FALSE)
  expect_true(fit$mleconv)
  expect_equal(fit$censtype, "double")
})

test_that("gaussian left-tail: delay = 0 row matches truncated oracle", {
  # Regression pin for the delay = 0 left-tail bug. At mean = 0, sd = 1, a
  # row with x2l - x1l = 0 must not fall back to primarycensored's d <= 0
  # short-circuit in pcens_cdf.default() (which hardcodes F_cens(0) = 0 on
  # the raw pnorm CDF and leaves an incoherent subprobability on the
  # interval [0, swindow]). The correct interpretation is that the delay is
  # a non-negative truncated normal; with the wrapper in place this row
  # agrees with an inline truncated-normal oracle.
  skip_if_no_primarycensored()
  x <- data.frame(
    x1l = 0,
    x1r = 1,
    x2l = 0,
    x2r = 1
  )
  v <- c(0, log(1))  # mean = 0, sd = 1

  # Inline zero-truncated normal CDF used purely as a regression oracle.
  ptrunc <- function(q, mean, sd) {
    f0 <- stats::pnorm(0, mean, sd)
    out <- (stats::pnorm(q, mean, sd) - f0) / (1 - f0)
    out[q <= 0] <- 0
    out
  }
  expected <- primarycensored::dprimarycensored(
    x = x$x2l - x$x1l, pdist = ptrunc,
    pwindow = x$x1r - x$x1l, swindow = x$x2r - x$x2l,
    mean = v[1], sd = exp(v[2]), log = TRUE
  )

  m <- kerlikelihood(x = x, family = "gaussian", likapprox = "ni")
  expect_equal(m$loglik(v, x), sum(expected), tolerance = 1e-8)
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
    x = x, pdist = ptruncnorm_nonneg,
    pars = list(mean = v[1], sd = exp(v[2]))
  )
  expect_equal(m$loglik(v, x), expected, tolerance = 1e-8)
})

# Delta-documenting tests.
#
# These pin the intentional semantic delta between the upstream
# integrate() reference (which uses the raw CDF on the full real line) and
# the truncated reference (which treats the delay as non-negative). For
# positive-support families the two agree trivially; for signed-support
# families the two differ by design.

test_that("upstream vs truncated oracles agree on positive-support families", {
  skip_if_no_primarycensored()
  x <- make_double_data()
  for (family in c("gamma", "lognormal", "weibull")) {
    case <- family_cases[[family]]
    upstream <- kerlik_integrate_reference_upstream(
      x = x, pdist = case$pdist, pars = case$pars
    )
    truncated <- kerlik_integrate_reference_truncated(
      x = x, pdist = case$pdist, pars = case$pars
    )
    # F(q) = 0 for q < 0 so the renormaliser is 1 and both oracles collapse
    # to the same value.
    expect_equal(upstream, truncated, tolerance = 1e-10, info = family)
  }
})

test_that("upstream vs truncated oracles differ for gaussian at mean=0,sd=1", {
  # Signed-support families at saturating params sit deep in the regime
  # where the untruncated and truncated semantics disagree. At mean = 0,
  # sd = 1 on the test frame the upstream integrate() oracle is smaller
  # (more negative loglik) than the truncated oracle by roughly 3-5 nats
  # because half of each row's mass under the untruncated pnorm sits on
  # the negative half-line and is lost to the d <= 0 short-circuit.
  skip_if_no_primarycensored()
  x <- make_double_data()
  pars <- list(mean = 0, sd = 1)
  upstream <- kerlik_integrate_reference_upstream(
    x = x, pdist = stats::pnorm, pars = pars
  )
  truncated <- kerlik_integrate_reference_truncated(
    x = x, pdist = ptruncnorm_nonneg, pars = pars
  )
  expect_false(isTRUE(all.equal(upstream, truncated)))
  # Upstream has more mass on the unconstrained real line and therefore
  # assigns less density to each observed interval on the positive side.
  expect_lt(upstream, truncated)
})

test_that("upstream vs truncated differ for skewnorm at saturating params", {
  skip_if_no_primarycensored()
  x <- make_double_data()
  # slant = 5 puts substantial mass to the right of zero already; at
  # location = 0, scale = 1 the zero truncation still redistributes a
  # non-negligible chunk of probability and the two oracles diverge.
  pars <- list(location = 0, scale = 1, slant = 5)
  raw_pskewnorm <- function(q, location, scale, slant) {
    EpiDelays::pskewnorm(x = q, par1 = location, par2 = scale, par3 = slant)
  }
  upstream <- kerlik_integrate_reference_upstream(
    x = x, pdist = raw_pskewnorm, pars = pars
  )
  truncated <- kerlik_integrate_reference_truncated(
    x = x, pdist = ptruncskewnorm_nonneg, pars = pars
  )
  expect_false(isTRUE(all.equal(upstream, truncated)))
  expect_lt(upstream, truncated)
})

test_that("upstream vs truncated agree for gaussian in the safe regime", {
  # At mean = 2.5, sd = 1 on make_double_data(), pnorm(0, mean, sd) is
  # tiny (~0.006) and the renormaliser is effectively 1. Upstream and
  # truncated agree within roughly 3e-2 on the summed loglik across 5
  # rows — useful context for anyone debugging the benchmark or wondering
  # why the default test params don't expose the delta. The remaining gap
  # comes from the zero-truncation redistributing the ~0.6% of mass below
  # zero onto the positive half-line, which slightly inflates each row's
  # truncated density relative to the upstream full-real oracle.
  skip_if_no_primarycensored()
  x <- make_double_data()
  pars <- list(mean = 2.5, sd = 1)
  upstream <- kerlik_integrate_reference_upstream(
    x = x, pdist = stats::pnorm, pars = pars
  )
  truncated <- kerlik_integrate_reference_truncated(
    x = x, pdist = ptruncnorm_nonneg, pars = pars
  )
  expect_equal(upstream, truncated, tolerance = 5e-2)
})
