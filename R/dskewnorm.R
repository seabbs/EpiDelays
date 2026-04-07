#' Density function of the Skew-Normal Distribution
#'
#' @param x A vector of quantiles.
#' @param par1 Location parameter of the skew-normal distribution.
#' @param par2 Scale parameter of the skew-normal distribution.
#' @param par3 Slant parameter of the skew-normal distribution.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references
#' Owen, D. B. (1956). Tables for Computing Bivariate Normal Probabilities.
#' \emph{The Annals of Mathematical Statistics}, \strong{27}(4), 1075–1090.
#'
#' Azzalini, A., and Capitanio, A. (2014). The Skew-Normal and Related Families.
#' \emph{Cambridge University Press}.
#'
#' @examples
#' dskewnorm(x = c(-0.5,0.9))
#'
#' @export

dskewnorm <- function(x, par1 = 0, par2 = 1, par3 = 0) {
  z <- (x - par1) / par2
  o <- (2 / par2) * stats::dnorm(z) * stats::pnorm(z * par3)
  return(o)
}
