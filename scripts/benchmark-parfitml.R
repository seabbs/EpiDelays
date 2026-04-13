# End-to-end parfitml() benchmark across the five parametric families.
#
# Times a single representative parfitml() fit for each
# (family, n) cell on doubly-interval-censored data simulated via
# primarycensored::rprimarycensored(). Ground-truth parameters per family
# mirror those in tests/testthat/test-kerlikelihood-primarycensored.R so
# the benchmark stays consistent with the existing test oracles.
#
# Not a unit test. Excluded from the package build via .Rbuildignore.
#
# Reproduce from the repo root with:
#   Rscript scripts/benchmark-parfitml.R

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

set.seed(20260413L)

# Ground-truth parameters per family. Matched to the `family_cases` list
# in tests/testthat/test-kerlikelihood-primarycensored.R so that the
# benchmark sims stay consistent with existing oracles.
family_truth <- list(
  gaussian = list(
    rdist = stats::rnorm,
    pars = list(mean = 1.5, sd = 0.8)
  ),
  gamma = list(
    rdist = stats::rgamma,
    pars = list(shape = 2, rate = 0.5)
  ),
  lognormal = list(
    rdist = stats::rlnorm,
    pars = list(meanlog = 0.5, sdlog = 0.4)
  ),
  weibull = list(
    rdist = stats::rweibull,
    pars = list(shape = 2, scale = 2)
  ),
  skewnorm = list(
    # sn::rsn would be the natural sampler but sn is not a dependency.
    # Inverse-CDF sampling via EpiDelays::qskewnorm keeps the benchmark
    # dependency-free and matches the (location, scale, slant) truth
    # used in tests/testthat/test-kerlikelihood-primarycensored.R.
    rdist = function(n) EpiDelays::qskewnorm(
      stats::runif(n), par1 = 1, par2 = 1, par3 = 2
    ),
    pars = list()
  )
)

simulate_double <- function(n, rdist, pars, pwindow = 1, swindow = 1,
                            seed = 1L, nonneg = FALSE) {
  set.seed(seed)
  # Scalar primary and secondary windows reflect the common integer-day
  # case and avoid the per-row vapply that the tests use for mixed windows.
  # For non-negative-support families parfitml rejects rows where
  # x2l - x1r < 0 (i.e. delay < pwindow); over-sample and filter so the
  # benchmark cell always hits the requested n.
  oversample <- if (nonneg) ceiling(n * 4) + 200L else n
  delays <- do.call(
    primarycensored::rprimarycensored,
    c(
      list(
        n = oversample,
        rdist = rdist,
        pwindow = pwindow,
        swindow = swindow
      ),
      pars
    )
  )
  if (nonneg) {
    delays <- delays[delays >= pwindow]
    if (length(delays) < n) {
      stop("insufficient non-negative delays simulated; bump oversample")
    }
    delays <- delays[seq_len(n)]
  }
  x1l <- sort(stats::runif(n, 0, 5))
  data.frame(
    x1l = x1l,
    x1r = x1l + pwindow,
    x2l = x1l + delays,
    x2r = x1l + delays + swindow
  )
}

nonneg_families <- c("gamma", "lognormal", "weibull")

ns <- c(50L, 200L, 1000L)
families <- names(family_truth)
Bboot <- 100L

rows <- list()
total_tic <- proc.time()

for (family in families) {
  truth <- family_truth[[family]]
  for (n in ns) {
    x <- simulate_double(
      n = n, rdist = truth$rdist, pars = truth$pars, seed = n,
      nonneg = family %in% nonneg_families
    )

    fit <- tryCatch(
      {
        tic <- proc.time()
        f <- parfitml(
          x = x, family = family, Bboot = Bboot, pgbar = FALSE
        )
        toc <- proc.time() - tic
        list(
          ok = TRUE,
          elapsed = as.numeric(toc[3]),
          converged = isTRUE(f$mleconv),
          pars = vapply(
            f$parfit,
            function(p) as.numeric(p$point),
            numeric(1)
          )
        )
      },
      error = function(e) {
        list(
          ok = FALSE,
          elapsed = NA_real_,
          converged = FALSE,
          pars = NA_real_,
          msg = conditionMessage(e)
        )
      }
    )

    pars_str <- if (fit$ok) {
      paste(sprintf("%.3f", fit$pars), collapse = ", ")
    } else {
      "DNF"
    }
    cat(sprintf(
      "%-9s n=%4d Bboot=%d  elapsed=%7.2fs  conv=%s  pars=[%s]\n",
      family, n, Bboot,
      if (is.null(fit$elapsed) || is.na(fit$elapsed)) NA_real_
      else fit$elapsed,
      if (!fit$ok) "DNF" else as.character(fit$converged),
      pars_str
    ))

    rows[[length(rows) + 1L]] <- data.frame(
      family = family,
      n = n,
      Bboot = Bboot,
      elapsed_s = if (fit$ok) round(fit$elapsed, 2) else NA_real_,
      converged = if (fit$ok) fit$converged else NA,
      pars = pars_str,
      stringsAsFactors = FALSE
    )
  }
}

total_toc <- proc.time() - total_tic

bench <- do.call(rbind, rows)
cat("\n--- benchmark results ---\n")
print(bench, row.names = FALSE)
cat(sprintf(
  "\nTotal wall-clock: %.1fs\n", as.numeric(total_toc[3])
))
cat(sprintf(
  "R %s on %s\n",
  paste(R.version$major, R.version$minor, sep = "."),
  R.version$platform
))

invisible(bench)
