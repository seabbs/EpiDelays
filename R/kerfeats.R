#' Computes features for a given parametric family
#'
#' @description
#' Computes the mean, variance, standard deviation and different quantiles
#' given a parametric family and parameter values.
#'
#' @param family A characther string specifying the parametric family.
#' @param par A vector of parameter values.
#'
#' @returns
#' A data frame of features.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerfeats <- function(family, par) {
  pfeats <- c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)
  if (family == "gaussian") {
    mean <- par[1]
    var <- par[2]^2
    sd <- par[2]
    quants <- stats::qnorm(p = pfeats, mean = par[1], sd = par[2])
  } else if (family == "skewnorm") {
    d3 <- par[3] / sqrt(1 + par[3]^2)
    mean <- par[1] + par[2] * d3 * sqrt(2 / pi)
    var <- par[2]^2 * (1 - (2 / pi) * d3^2)
    sd <- sqrt(var)
    quants <- qskewnorm(p = pfeats, par1 = par[1], par2 = par[2], par3 = par[3])
  } else if (family == "gamma") {
    mean <- par[1] / par[2]
    var <- par[1] / (par[2]^2)
    sd <- sqrt(var)
    quants <- stats::qgamma(p = pfeats, shape = par[1], rate = par[2])
  } else if (family == "lognormal") {
    mean <- exp(par[1] + 0.5 * par[2]^2)
    var <- exp(2 * par[1] + par[2]^2) * (exp(par[2]^2) - 1)
    sd <- sqrt(var)
    quants <- stats::qlnorm(p = pfeats, meanlog = par[1], sdlog = par[2])
  } else if (family == "weibull") {
    mean <- par[2] * gamma(1 + 1 / par[1])
    var <- par[2]^2 * (gamma(1 + 2 / par[1]) - (gamma(1 + 1 / par[1])^2))
    sd <- sqrt(var)
    quants <- stats::qweibull(p = pfeats, shape = par[1], scale = par[2])
  }
  o <- data.frame(mean = mean, var = var, sd = sd,
                  q1 = quants[1], q5 = quants[2], q25 = quants[3],
                  q50 = quants[4], q75 = quants[5], q95 = quants[6],
                  q99 = quants[7])
  return(o)
}






