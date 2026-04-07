#' Description of available parametric families
#'
#' @description
#' Provides a summary of available parametric families.
#'
#' @returns
#' A list summarizing available parametric families, the number of
#' parameters involved, the parameter names and the support.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerfamilies <- function() {
  f1 <- list(
    fname = "gaussian",
    npars = 2,
    par1 = "mean",
    par2 = "sd",
    support = "real"
  )
  f2 <- list(
    fname = "gamma",
    npars = 2,
    par1 = "shape",
    par2 = "rate",
    support = "non-neg"
  )
  f3 <- list(
    fname = "lognormal",
    npars = 2,
    par1 = "location",
    par2 = "scale",
    support = "non-neg"
  )
  f4 <- list(
    fname = "weibull",
    npars = 2,
    par1 = "shape",
    par2 = "scale",
    support = "non-neg"
  )
  f5 <- list(
    fname = "skewnorm",
    npars = 3,
    par1 = "location",
    par2 = "scale",
    par3 = "slant",
    support = "real"
  )
  return(list(f1, f2, f3, f4, f5))
}
