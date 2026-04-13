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
#' The doubly interval-censored likelihood function requires to evaluate an
#' inner integral. By default, this is implemented via numerical integration
#' (ni) with \code{likapprox = "ni"}, where the integral is evaluated with
#' \code{primarycensored::dprimarycensored()}. Another option is
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
#' interval-censored data. Default is \code{"ni"} for numerical integration.
#' Another possibility is \code{"mc"} for crude Monte Carlo sampling with 1000
#' samples.
#'
#' @return A list containing information on the chosen parametric family,
#' the log-likelihood function, a function that transforms back the parameters
#' in their original scale, the censoring type and the approximation of the
#' likelihood function used in case of doubly interval-censored data.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerlikelihood <- function(x, family, likapprox = "ni") {
  if(!(likapprox %in% c("ni", "mc"))) {
    stop("likapprox should either be ni or mc")
  }
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
  } else if (nc == 4) {
    censtype <- "double"
    M <- 1000 #  Monte carlo samples for doubly interval-censored likelihood
  }
  # CDF wrapper exposing pskewnorm with a primarycensored-friendly signature.
  pskewnorm_cdf <- function(q, location, scale, slant) {
    pskewnorm(x = q, par1 = location, par2 = scale, par3 = slant)
  }
  # Shared loglik builder for the doubly interval-censored ni branch: given a
  # CDF and a function that extracts a named parameter list from `v`, returns
  # a closure that sums the per-row log-density from
  # primarycensored::dprimarycensored(). Rows are grouped by (pwindow, swindow)
  # so dprimarycensored's internal lookup table is reused within each group.
  build_pc_loglik <- function(pdist, pars_fn) {
    force(pdist)
    force(pars_fn)
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
              pwindow = pw, swindow = sw, log = TRUE
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
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- stats::pnorm(q = x$xl, mean = par1, sd = par2)
        Fr  <- stats::pnorm(q = x$xr, mean = par1, sd = par2)
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::pnorm(q = x$x2l[i] - x1s , mean = par1, sd = par2)
            Fr  <- stats::pnorm(q = x$x2r[i] - x1s , mean = par1, sd = par2)
            z <- z + log(mean(Fr - Fl))
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- build_pc_loglik(
          pdist = stats::pnorm,
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
        Fl <- pskewnorm(x = x$xl, par1 = par1, par2 = par2, par3 = par3)
        Fr <- pskewnorm(x = x$xr, par1 = par1, par2 = par2, par3 = par3)
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          par3 <- v[3]
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl <- pskewnorm(x = x$x2l[i] - x1s, par1 = par1, par2 = par2, par3 = par3)
            Fr <- pskewnorm(x = x$x2r[i] - x1s, par1 = par1, par2 = par2, par3 = par3)
            z <- z + log(mean(Fr - Fl))
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- build_pc_loglik(
          pdist = pskewnorm_cdf,
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
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::pgamma(q = x$x2l[i] - x1s, shape = par1, rate = par2)
            Fr  <- stats::pgamma(q = x$x2r[i] - x1s, shape = par1, rate = par2)
            z <- z + log(mean(Fr - Fl))
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
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::plnorm(q = x$x2l[i] - x1s, meanlog = par1, sdlog = par2)
            Fr  <- stats::plnorm(q = x$x2r[i] - x1s, meanlog = par1, sdlog = par2)
            z <- z + log(mean(Fr - Fl))
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
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::pweibull(q = x$x2l[i] - x1s, shape = par1, scale = par2)
            Fr  <- stats::pweibull(q = x$x2r[i] - x1s, shape = par1, scale = par2)
            z <- z + log(mean(Fr - Fl))
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
           censtype = censtype, likapprox = likapprox))
  return(o)
}
