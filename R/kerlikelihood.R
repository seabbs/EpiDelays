#' Likelihood function for single and doubly interval-censored data
#'
#' @description
#' The likelihood function plays a key role in maximum likelihood estimation and
#' Bayesian inference. This routine admits as data input a data frame \code{x}
#' with either two columns (named \code{xl} and \code{xr}) representing the left
#' and right bound, respectively, of the delay variable, or four columns
#' (named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}) representing the left
#' and right bound, respectively, of the primary and secondary events of the
#' delay variable. The naming convention of the columns is strict and different
#' namings are not allowed. When data frame \code{x} has two columns, a single
#' interval-censored likelihood function is used. Note that the left bound
#' should be strictly smaller than the right bound, i.e. \code{xl < xr} must
#' be satisfied for all rows in \code{x}. When data frame \code{x} has four
#' columns, a doubly interval-censored likelihood function is used. In that case
#' \code{x1l < x1r} and \code{x2l < x2r} must hold for all rows in \code{x}.
#' Moreover, \code{NA} values are not allowed in data frame \code{x}.
#'
#' @details
#' Returns the log-likelihood function of parameters that have been transformed
#' to live in an unbounded parameter space. This makes the
#' maximization/evaluation of the likelihood function numerically more stable.
#' The doubly interval-censored likelihood function marginalises out the
#' primary event time. By default (\code{likapprox = "ni"}), this
#' marginalisation is delegated to
#' \code{primarycensored::dprimarycensored()}, which chooses the most
#' efficient method for the given CDF and windows. Another option is
#' \code{likapprox = "mc"}, which uses crude Monte Carlo sampling to evaluate
#' the inner integral in the doubly interval-censored likelihood function.
#' Note that Monte Carlo sampling results in much longer computation times to
#' evaluate the likelihood function. The number of Monte Carlo samples is
#' fixed at 1000.
#'
#' @param x A data frame with either two columns named \code{xl} and \code{xr},
#' or four columns named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}. See
#' description for constraints imposed on the columns.
#' @param family A character string specifying the name of the parametric
#' family.
#' @param likapprox Approximation of the likelihood function in case of doubly
#' interval-censored data. Default is \code{"ni"}, which delegates evaluation
#' to \code{primarycensored::dprimarycensored()}. Another possibility is
#' \code{"mc"} for crude Monte Carlo sampling with 1000 samples.
#' @param L Lower truncation point applied to the underlying delay
#' distribution. Defaults to \code{0}. When \code{L > 0} the per-row
#' contribution is rescaled by the truncated CDF.
#' @param D Upper truncation point applied to the underlying delay
#' distribution. Defaults to \code{Inf}. Together with \code{L}, rows whose
#' observed interval is incompatible with \code{[L, D]} are rejected at
#' input time.
#' @param dprimary Primary event density function. Used only when \code{x} has
#' four columns (doubly interval-censored data) and \code{likapprox = "ni"}.
#' Must take a vector \code{x} and the arguments \code{min} and \code{max},
#' return a density normalised to integrate to 1 on \code{[min, max]}, and be
#' compatible with \code{primarycensored::dprimarycensored()}. Defaults to
#' \code{stats::dunif} (uniform primary onset within the observation window),
#' which reproduces the behaviour of earlier EpiDelays versions. Non-uniform
#' choices include \code{primarycensored::dexpgrowth} for exponential growth
#' during an outbreak. Ignored when \code{x} has two columns (single
#' interval-censored data) but must equal the default in that case.
#' @param dprimary_args A named list of additional arguments passed to
#' \code{dprimary}, mirroring the \code{dprimary_args} argument of
#' \code{primarycensored::dprimarycensored()}. Example: \code{list(r = 0.1)}
#' for \code{primarycensored::dexpgrowth}. Defaults to an empty list. Must be
#' empty unless \code{x} has four columns and \code{likapprox = "ni"}.
#'
#' @return A list containing information on the chosen parametric family,
#' the log-likelihood function, a function that transforms back the parameters
#' in their original scale, the censoring type, the approximation of the
#' likelihood function used in case of doubly interval-censored data, and the
#' primary event density settings.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerlikelihood <- function(x, family, likapprox = "ni", L = 0, D = Inf,
                          dprimary = stats::dunif,
                          dprimary_args = list()) {
  if(!(likapprox %in% c("ni", "mc"))) {
    stop("likapprox should either be ni or mc")
  }
  # Truncation bound validation mirrors
  # primarycensored::.check_truncation_bounds.
  if (!is.numeric(L) || length(L) != 1L || is.na(L) || L < 0) {
    stop("L must be a non-negative scalar.")
  }
  if (!is.numeric(D) || length(D) != 1L || is.na(D) || L >= D) {
    stop("L must be less than D.")
  }
  # Validate dprimary / dprimary_args shape. Identity-check against the
  # default so users writing `dprimary = stats::dunif` explicitly still hit
  # the uniform-primary fast path.
  if(!is.function(dprimary)) {
    stop("dprimary must be a function")
  }
  if(!is.list(dprimary_args) ||
     (length(dprimary_args) > 0 && is.null(names(dprimary_args)))) {
    stop("dprimary_args must be a named list")
  }
  dprimary_default <- identical(dprimary, stats::dunif) &&
    length(dprimary_args) == 0
  # Input checks
  dfck <- kerdata_check(x = x) # data frame check
    if (dfck$result == "fail") {
      stop(dfck$message)
    }
  famck <- kerfamily_check(x = family) # family check
    if (famck$result == "fail") {
      stop(famck$message)
    }
  domck <- kerdomain_check(x = x, family = family) # Domain check
  if (domck$result == "fail") {
    stop(domck$message)
  }
  fset <- kerfamilies()
  fnames <- sapply(fset, "[[", "fname")
  famdesc <- fset[[match(family, fnames)]]
  n <- nrow(x)
  nc <- ncol(x)
  if(nc == 2) {
    censtype <- "single"
    likapprox <- "None"
    # Non-uniform primary has no meaning without a primary event window.
    if(!dprimary_default) {
      stop(
        "dprimary only applies to doubly interval-censored data (four ",
        "columns); x has two columns so no primary event window is ",
        "modelled"
      )
    }
    # Reject rows whose observed interval cannot intersect [L, D].
    if ((L > 0 || is.finite(D)) &&
        (any(x$xr <= L) || any(x$xl >= D))) {
      stop(
        "Some rows of x are incompatible with the requested truncation ",
        "bounds [L, D]: need xr > L and xl < D for every row."
      )
    }
  } else if (nc == 4) {
    censtype <- "double"
    M <- 1000 #  Monte carlo samples for doubly interval-censored likelihood
    # The mc branch samples primary onsets from a uniform; a non-uniform
    # dprimary would silently bias the Monte Carlo estimate of F_cens.
    if(likapprox == "mc" && !dprimary_default) {
      stop(
        "likapprox = \"mc\" does not support non-uniform dprimary; ",
        "use likapprox = \"ni\" (the default)"
      )
    }
    # Reject rows whose observed interval cannot intersect [L, D]. Rows that
    # straddle the boundary are retained and clamped inside the loglik
    # closures below.
    if (L > 0 || is.finite(D)) {
      lowers <- x$x2l - x$x1l
      uppers <- x$x2r - x$x1l
      if (any(uppers <= L) || any(lowers >= D)) {
        stop(
          "Some rows of x are incompatible with the requested truncation ",
          "bounds [L, D]: need x2r - x1l > L and x2l - x1l < D for every ",
          "row."
        )
      }
    }
  }
  # Helper: sum log-contributions for the single-interval (nc == 2) branch,
  # applying the truncation correction when L > 0 or D < Inf. Reduces to
  # sum(log(Fr - Fl)) in the default case, so the pre-existing code path is
  # preserved bit-for-bit.
  single_interval_sum <- function(Fl, Fr, FL, FD) {
    if (L <= 0 && is.infinite(D)) {
      return(sum(log(Fr - Fl)))
    }
    num <- pmin(Fr, FD) - pmax(Fl, FL)
    den <- FD - FL
    sum(log(num)) - length(Fl) * log(den)
  }
  # Helper: per-row mc truncation correction. Given the primary draws x1s
  # for a single row, subtract log(mean(F(D - x1s) - F(L - x1s))). Returns 0
  # in the default case, keeping the untruncated mc path unchanged.
  mc_row_correction <- function(x1s, pdist, pars) {
    if (L <= 0 && is.infinite(D)) {
      return(0)
    }
    if (is.finite(D)) {
      Fd <- do.call(pdist, c(list(D - x1s), pars))
    } else {
      Fd <- rep(1, length(x1s))
    }
    if (L > 0) {
      Fl <- do.call(pdist, c(list(L - x1s), pars))
    } else {
      Fl <- rep(0, length(x1s))
    }
    log(mean(Fd - Fl))
  }
  # Shared loglik builder for the doubly interval-censored ni branch. Given a
  # CDF and a function that extracts a named parameter list from `v`, returns
  # a closure that sums the per-row log-density from
  # primarycensored::dprimarycensored(). Rows are grouped by (pwindow,
  # swindow) because both are scalar arguments in dprimarycensored. L and D
  # are forwarded unchanged: the default bounds (L = 0, D = Inf) reduce to
  # the untruncated call and non-default bounds pick up the
  # F_cens(D) - F_cens(L) renormaliser inside primarycensored. The
  # user-supplied dprimary / dprimary_args are forwarded too.
  build_pc_loglik <- function(pdist, pars_fn) {
    force(pdist)
    force(pars_fn)
    # Close over the primary event density so parfitml() and its bootstrap
    # loop automatically re-use whatever dprimary the user passed in.
    dprimary_local <- dprimary
    dprimary_args_local <- dprimary_args
    function(v, x) {
      pars <- pars_fn(v)
      pwindows <- x$x1r - x$x1l
      swindows <- x$x2r - x$x2l
      lowers <- x$x2l - x$x1l
      groups <- split(
        seq_along(lowers), list(pwindows, swindows), drop = TRUE
      )
      z <- 0
      for (idx in groups) {
        pw <- pwindows[idx[1]]
        sw <- swindows[idx[1]]
        logd <- do.call(
          primarycensored::dprimarycensored,
          c(
            list(
              x = lowers[idx], pdist = pdist,
              pwindow = pw, swindow = sw,
              L = L, D = D,
              dprimary = dprimary_local,
              dprimary_args = dprimary_args_local,
              log = TRUE
            ),
            pars
          )
        )
        z <- z + sum(logd)
      }
      z
    }
  }
  if (family == "gaussian") {
    # Gaussian and skewnorm are interpreted as zero-truncated on the
    # non-negative half-line. The truncated wrapper becomes the underlying
    # delay CDF, so every code path (single-interval, mc, ni) routes through
    # the same ptruncnorm_nonneg without a bespoke F_cens(0) correction.
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- ptruncnorm_nonneg(q = x$xl, mean = par1, sd = par2)
        Fr  <- ptruncnorm_nonneg(q = x$xr, mean = par1, sd = par2)
        FL  <- ptruncnorm_nonneg(q = L, mean = par1, sd = par2)
        FD  <- ptruncnorm_nonneg(q = D, mean = par1, sd = par2)
        single_interval_sum(Fl, Fr, FL, FD)
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          pars <- list(mean = par1, sd = par2)
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- ptruncnorm_nonneg(q = x$x2l[i] - x1s, mean = par1, sd = par2)
            Fr  <- ptruncnorm_nonneg(q = x$x2r[i] - x1s, mean = par1, sd = par2)
            z <- z + log(mean(Fr - Fl)) -
              mc_row_correction(x1s, ptruncnorm_nonneg, pars)
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- build_pc_loglik(
          pdist = ptruncnorm_nonneg,
          pars_fn = function(v) list(mean = v[1], sd = exp(v[2]))
        )
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      z
    }
  }else if (family == "skewnorm") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        par3 <- v[3]
        Fl <- ptruncskewnorm_nonneg(
          q = x$xl, location = par1, scale = par2, slant = par3
        )
        Fr <- ptruncskewnorm_nonneg(
          q = x$xr, location = par1, scale = par2, slant = par3
        )
        FL <- ptruncskewnorm_nonneg(
          q = L, location = par1, scale = par2, slant = par3
        )
        FD <- ptruncskewnorm_nonneg(
          q = D, location = par1, scale = par2, slant = par3
        )
        single_interval_sum(Fl, Fr, FL, FD)
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          par3 <- v[3]
          pars <- list(location = par1, scale = par2, slant = par3)
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl <- ptruncskewnorm_nonneg(
              q = x$x2l[i] - x1s,
              location = par1, scale = par2, slant = par3
            )
            Fr <- ptruncskewnorm_nonneg(
              q = x$x2r[i] - x1s,
              location = par1, scale = par2, slant = par3
            )
            z <- z + log(mean(Fr - Fl)) -
              mc_row_correction(x1s, ptruncskewnorm_nonneg, pars)
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- build_pc_loglik(
          pdist = ptruncskewnorm_nonneg,
          pars_fn = function(v) list(
            location = v[1], scale = exp(v[2]), slant = v[3]
          )
        )
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]), v[3])
      colnames(z) <- c(famdesc$par1, famdesc$par2, famdesc$par3)
      z
    }
  } else if (family == "gamma") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- exp(v[1])
        par2 <- exp(v[2])
        Fl  <- stats::pgamma(q = x$xl, shape = par1, rate = par2)
        Fr  <- stats::pgamma(q = x$xr, shape = par1, rate = par2)
        FL  <- stats::pgamma(q = L, shape = par1, rate = par2)
        FD  <- stats::pgamma(q = D, shape = par1, rate = par2)
        single_interval_sum(Fl, Fr, FL, FD)
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          pars <- list(shape = par1, rate = par2)
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::pgamma(q = x$x2l[i] - x1s, shape = par1, rate = par2)
            Fr  <- stats::pgamma(q = x$x2r[i] - x1s, shape = par1, rate = par2)
            z <- z + log(mean(Fr - Fl)) -
              mc_row_correction(x1s, stats::pgamma, pars)
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- build_pc_loglik(
          pdist = stats::pgamma,
          pars_fn = function(v) list(shape = exp(v[1]), rate = exp(v[2]))
        )
      }
    }
    originscale <- function(v) {
      z <- data.frame(exp(v[1]), exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      z
    }
  } else if (family == "lognormal") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- stats::plnorm(q = x$xl, meanlog = par1, sdlog = par2)
        Fr  <- stats::plnorm(q = x$xr, meanlog = par1, sdlog = par2)
        FL  <- stats::plnorm(q = L, meanlog = par1, sdlog = par2)
        FD  <- stats::plnorm(q = D, meanlog = par1, sdlog = par2)
        single_interval_sum(Fl, Fr, FL, FD)
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          pars <- list(meanlog = par1, sdlog = par2)
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::plnorm(q = x$x2l[i] - x1s, meanlog = par1, sdlog = par2)
            Fr  <- stats::plnorm(q = x$x2r[i] - x1s, meanlog = par1, sdlog = par2)
            z <- z + log(mean(Fr - Fl)) -
              mc_row_correction(x1s, stats::plnorm, pars)
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- build_pc_loglik(
          pdist = stats::plnorm,
          pars_fn = function(v) list(meanlog = v[1], sdlog = exp(v[2]))
        )
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      z
    }
  } else if (family == "weibull") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- exp(v[1])
        par2 <- exp(v[2])
        Fl  <- stats::pweibull(q = x$xl, shape = par1, scale = par2)
        Fr  <- stats::pweibull(q = x$xr, shape = par1, scale = par2)
        FL  <- stats::pweibull(q = L, shape = par1, scale = par2)
        FD  <- stats::pweibull(q = D, shape = par1, scale = par2)
        single_interval_sum(Fl, Fr, FL, FD)
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          pars <- list(shape = par1, scale = par2)
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::pweibull(q = x$x2l[i] - x1s, shape = par1, scale = par2)
            Fr  <- stats::pweibull(q = x$x2r[i] - x1s, shape = par1, scale = par2)
            z <- z + log(mean(Fr - Fl)) -
              mc_row_correction(x1s, stats::pweibull, pars)
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- build_pc_loglik(
          pdist = stats::pweibull,
          pars_fn = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
        )
      }
    }
    originscale <- function(v) {
      z <- data.frame(exp(v[1]), exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      z
    }
  }
  o <- c(famdesc, list(loglik = loglik, originscale = originscale,
           censtype = censtype, likapprox = likapprox,
           dprimary = dprimary, dprimary_args = dprimary_args))
  return(o)
}
