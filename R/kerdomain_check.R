#' Check if the sample space of parametric family is compatible with data
#'
#' @param x A data frame with either two columns named \code{xl} and \code{xr},
#' or four columns named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}.
#' @param family A character string specifying the name of the parametric
#' family.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#'
#' @keywords internal

kerdomain_check <- function(x, family){
  fset <- kerfamilies()
  fnames <- sapply(fset, "[[", "fname")
  famdesc <- fset[[match(family, fnames)]]
  nc <- ncol(x)
  res <- "pass"
  msg <- "Parametric family sample space compatible with data"
  if(nc == 2) {
    if((famdesc$support == "non-neg") & any(x$xl < 0)) {
      res <- "fail"
      msg <- "Chosen parametric family has a non-negative sample space but data frame x contains at least one negative observed left bound xl. The data suggests that negative values for the delay variable cannot be excluded. Choose a parametric family with a sample space covering negative values (e.g. family = 'gaussian' or family = 'skewnorm')."
    }
  } else if(nc == 4) {
    xl <- x$x2l - x$x1r
    if((famdesc$support == "non-neg") & any(xl < 0)) {
      res <- "fail"
      msg <- "Chosen parametric family has a non-negative sample space but data frame x contains at least one negative observed left bound xl = x$x2l - x$x1r. The data suggests that negative values for the delay variable cannot be excluded. Choose a parametric family with a sample space covering negative values (e.g. family = 'gaussian' or family = 'skewnorm')."
    }
  }
  o <- list(result = res, message = msg)
  return(o)
}
