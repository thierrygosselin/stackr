#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

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
