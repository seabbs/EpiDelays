#' Bootstrap resampling routine
#'
#' @description
#' Resample with replacement rows of a data frame.
#'
#' @param x A data frame.
#'
#' @returns
#' A data frame representing a bootstrap sample of `x`. The number of
#' bootstrap samples equals the number of rows of `x`.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerboot <- function(x) {
  n <- nrow(x)
  rowid <- seq_len(n)
  return(x[sample(rowid, size = n, replace = TRUE), ])
}
