#' Internal check of parametric family name provided by user
#'
#' @description
#' Checks if the parametric family name specified by user is available.
#'
#' @param x A string representing a parametric family name.
#'
#' @returns
#' A list with "pass" (if `x` is an available family) or "fail" if
#' (`x` is not an available family) and a message.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerfamily_check <- function(x) {
  fset <- kerfamilies()
  fnames <- sapply(fset, "[[", "fname")
  res <- "pass"
  msg <- "Family check passed"
  if (!is.character(x)) {
    res <- "fail"
    msg <- "Family must be of type character"
  }
  else if (!(x %in% fnames)) {
    res <- "fail"
    msg <- paste0("Specify a family among: ", paste(fnames, collapse = ", "))
  }
  return(list(result = res, message = msg))
}
