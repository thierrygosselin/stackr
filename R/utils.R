#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' @importFrom utils packageDescription
#' @importFrom stringi stri_join

.onAttach <- function(libname, pkgname) {
  stackr.version <- utils::packageDescription("stackr", fields = "Version")

  startup.message <- stringi::stri_join("
******************************* IMPORTANT NOTICE *******************************\n",
"stackr v.", stackr.version, "\n",
"stackr was modified heavily and now focuses exclusively on running stacks
software in R and reading particular files it produces.\n
For filters and other functions that was in stackr,
please see my new package called radiator.
https://github.com/thierrygosselin/radiator
********************************************************************************",
sep = "")
packageStartupMessage(startup.message)
}


#' @title split_vec_row
#' @description Split input into chunk for parallel processing
#' @rdname split_vec_row
#' @keywords internal
#' @export
split_vec_row <- function(x, cpu.rounds, parallel.core = parallel::detectCores() - 1) {
  n.row <- nrow(x)
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))
  return(split.vec)
}#End split_vec_row


.onUnload <- function(libpath) {
  library.dynam.unload("stackr", libpath)
}
