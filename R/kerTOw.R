#' Evaluation of Owen's \eqn{T(h,a)} function
#'
#' @description
#' Evaluates the function \eqn{T(h,a)} described in Owen (1956).
#'
#' @param h A real number.
#' @param a A real number.
#' @param kmax An integer to truncate the infinite series expression of \eqn{T(h,a)}.
#'
#' @returns
#' A real number.
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
#' @keywords internal

kerTOw <- function(h, a, kmax = 50) {
  TOwseries <- function(hs, as) {
    c <- 0.5 * hs^2
    s <- seq(0, kmax, by = 1)
    gs <- 2 * s + 1
    ck <- (((-1)^s * as^gs) / gs) * (1 - exp(-c) * cumsum(c^s / gamma(s + 1)))
    v <- (1 / (2 * pi)) * (atan(as) - sum(ck, na.rm = TRUE))
    return(v)
  }
  if (a == 0 || is.infinite(h))
    o <- 0
  else if (abs(a) <= 1) {
    o <- TOwseries(hs = h, as = a)
  } else{
    o <- 0.5 * (stats::pnorm(h) + stats::pnorm(a * h)) -
      stats::pnorm(h) * stats::pnorm(a * h) -
      TOwseries(hs = (a * h), as = (1 / a)) - 0.5 * as.numeric(a < 0)
  }
  return(o)
}
