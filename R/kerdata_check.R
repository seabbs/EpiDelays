#' Internal check of data frame provided by user
#'
#' @description
#' Checks if the data frame provided by user is in correct format.
#'
#' @param x A data frame.
#'
#' @returns
#' A list with "pass" (if `x` is in correct format) or "fail" if
#' (`x` is not in correct format) and a message.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerdata_check <- function(x) {
  res <- "pass"
  msg <- "Data frame in correct format"
  if (!is.data.frame(x)) {
    res <- "fail"
    msg <- "Input data should be a data frame"
  } else if (anyNA(x)) {
    res <- "fail"
    msg <- "Data frame cannot contain missing values"
  } else {
    nc <- ncol(x)
    if (!(nc %in% c(2, 4))) {
      res <- "fail"
      msg <- "Data frame must have either 2 or 4 columns"
    } else {
      cnames <- colnames(x)
      if (nc == 2) {
        if (!identical(cnames, c("xl", "xr"))) {
          res <- "fail"
          msg <- "Column names of data frame should be xl and xr"
        } else if (any(x$xl >= x$xr)) {
          res <- "fail"
          msg <- "xl < xr not satisfied in all rows of data frame"
        }
      } else if (nc == 4) {
        if (!identical(cnames, c("x1l", "x1r", "x2l", "x2r"))) {
          res <- "fail"
          msg <- "Column names of data frame should be x1l, x1r, x2l, x2r"
        } else if (any(x$x1l >= x$x1r | x$x2l >= x$x2r)) {
          res <- "fail"
          msg <- "x1l < x1r and x2l < x2r not satisfied in all rows of data frame"
        }
      }
    }
  }
  return(list(result = res, message = msg))
}

