# magrittr ---------------------------------------------------------------------
# Forward-pipe operator --------------------------------------------------------
#' @title Forward-pipe operator
#' @description magrittr forward-pipe operator
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

# Exposition pipe-operator ------------------------------------------------------
#' @title Exposition pipe-operator
#' @description magrittr Exposition pipe-operator
#' @name %$%
#' @rdname Exposition_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %$%
#' @usage lhs \%$\% rhs
NULL

# Compound assignment pipe operator --------------------------------------------
#' @title compound assignment pipe operator
#' @description magrittr compound assignment pipe operator
#' @name %<>%
#' @rdname compound_assignment_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL

# dplyr n ----------------------------------------------------------------------
# The number of observations in the current group.
#' @title The number of observations in the current group.
#' @description Check dplyr
#' @name n
#' @rdname n
#' @keywords internal
#' @export
#' @importFrom dplyr n
#' @usage n()
NULL

# importFrom -------------------------------------------------------------------
#' @importFrom pbmcapply pbmclapply
#' @importFrom stringdist stringdist
#' @importFrom dplyr group_by tally ungroup summarise summarise_if mutate filter distinct select bind_rows n_distinct arrange
#' @importFrom ShortRead readFastq


# .onAttach <- function(libname, pkgname) {
#   stackr.version <- utils::packageDescription("stackr", fields = "Version")
#
#   startup.message <- stringi::stri_join("
# ******************************* IMPORTANT NOTICE *******************************\n",
# "stackr v.", stackr.version, "\n",
# "stackr was modified heavily and now focuses exclusively on running stacks
# software in R and reading particular files it produces.\n
# For filters and other functions that were in stackr,
# please see my new package called radiator.
# https://github.com/thierrygosselin/radiator
#
#
# Reproducibility note: all zenodo DOI and stackr versions were left intact.
# ********************************************************************************",
# sep = "")
# packageStartupMessage(startup.message)
# }


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


# .onUnload <- function(libpath) {
#   library.dynam.unload("stackr", libpath)
# }



#' @title split_tibble_rows
#' @description Split rows of tibble for parallel processing
#' @rdname split_tibble_rows
#' @keywords internal
#' @export
split_tibble_rows <- function(
  x,
  lines.cpu = 1000, #lines per CPU rounds
  parallel.core = parallel::detectCores() - 1,
  group.split = TRUE # does dplyr: group_by and group_split
) {
  n.row <- nrow(x)
  n.cores <- parallel::detectCores()
  if (parallel.core > n.cores) parallel.core <- n.cores
  if (n.row < parallel.core) parallel.core <- n.row
  if (lines.cpu > n.row) lines.cpu <- n.row
  lines.rounds <- parallel.core * lines.cpu
  x$SPLIT_VEC <- sort(rep_len(x = 1:floor(n.row / lines.rounds), length.out = n.row))
  if (group.split) {
    x %<>%
      dplyr::group_by(SPLIT_VEC) %>%
      dplyr::group_split(.tbl = ., keep = FALSE)
  }
  return(x)
}#End split_tibble_rows
