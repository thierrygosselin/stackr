#' @name build_stackr_workflow_dir
#' @title Build stacks workflow directories
#' @description This prepares the different folder required to run the
#' \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS} workflow
#' proposed in stackr.

#' @param main.folder.name Use this argument to customize the name of the main
#' folder for workflow. e.g. \code{main.folder.name = "stacks_lobster"}
#' Default: \code{main.folder.name = "stacks_run"}

#' @param date (logical) Should the current date be appended on the main folder
#' name? e.g. \code{stacks_lobster_20160630}
#' Default: \code{date = TRUE}


#' @rdname build_stackr_workflow_dir
#' @export
#' @import stringi
#' @import dplyr
#' @import readr

#' @return A main folder containing different folders used during the stacks
#' workflow proposed here.

#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' build_stackr_workflow_dir()
#' # that's it ! Now if you have your own main folder name:
#' build_stackr_workflow_dir(main.folder.name = "stacks_whiteshark", date = TRUE)
#' }

# # required to pass the R CMD check and have 'no visible binding for global variable'
# if (getRversion() >= "2.15.1") {
#   utils::globalVariables(
#     c("INDIVIDUALS_REP")
#   )
# }

build_stackr_workflow_dir <- function (
  main.folder.name = "stacks_run",
  date = TRUE
){

  # Get date and time to have unique filenaming --------------------------------
  file.date <- stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")

  # Main folder ----------------------------------------------------------------
  if (date) {
    # main.folder.name <- "stacks_test"
    main.folder.name <- stri_paste(main.folder.name, file.date, sep = "_")
  }
  dir.create(path = main.folder.name)

  # Workflows folder -----------------------------------------------------------
  dir.create(path = stri_paste(main.folder.name, "/01_scripts"))
  dir.create(path = stri_paste(main.folder.name, "/02_project_info"))
  dir.create(path = stri_paste(main.folder.name, "/03_sequencing_lanes"))
  dir.create(path = stri_paste(main.folder.name, "/04_process_radtags"))
  dir.create(path = stri_paste(main.folder.name, "/05_paralogy_clustering_mismatch"))
  dir.create(path = stri_paste(main.folder.name, "/06_ustacks_cstacks_sstacks"))
  dir.create(path = stri_paste(main.folder.name, "/07_rxstacks_cstacks_sstacks_populations"))
  dir.create(path = stri_paste(main.folder.name, "/08_stacks_results"))
  dir.create(path = stri_paste(main.folder.name, "/09_log_files"))
  dir.create(path = stri_paste(main.folder.name, "/10_filters"))
} # end build_stackr_workflow_dir
