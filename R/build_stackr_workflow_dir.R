#' @name build_stackr_workflow_dir
#' @title Automatically build the stacks/stackr workflow directories
#' @description This prepares the different folders required to run the
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
#' @return A main folder containing different folders used during the stacks
#' workflow proposed here.
#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' stackr::build_stackr_workflow_dir()
#' # that's it ! Now if you have your own main folder name:
#' build_stackr_workflow_dir(main.folder.name = "stacks_whiteshark", date = TRUE)
#' }


build_stackr_workflow_dir <- function(
  main.folder.name = "stacks_run",
  date = TRUE
){

  # Get date and time to have unique filenaming --------------------------------
  file.date <- stringi::stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")

  # Main folder ----------------------------------------------------------------
  if (date) {
    main.folder.name <- stringi::stri_join(main.folder.name, file.date, sep = "_")
  }
  dir.create(path = main.folder.name)

  # Workflows folder -----------------------------------------------------------
  dir.create(path = file.path(main.folder.name, "01_scripts"))
  dir.create(path = file.path(main.folder.name, "02_project_info"))
  dir.create(path = file.path(main.folder.name, "03_sequencing_lanes"))
  dir.create(path = file.path(main.folder.name, "04_process_radtags"))
  dir.create(path = file.path(main.folder.name, "05_clustering_mismatches"))
  dir.create(path = file.path(main.folder.name, "06_ustacks_2_gstacks"))
  dir.create(path = file.path(main.folder.name, "07_populations"))
  dir.create(path = file.path(main.folder.name, "08_stacks_results"))
  dir.create(path = file.path(main.folder.name, "09_log_files"))
  dir.create(path = file.path(main.folder.name, "10_filters"))
} # end build_stackr_workflow_dir
