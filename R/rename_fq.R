#' @name rename_fq
#' @title Rename fastq files
#' @description Change of heart concerning the names of your fastq files?
#'
#' Renaming them shouldn't be a complicated and dangerous venture, even for people with no
#' UNIX bash experience.
#'
#' This function allows you to rename your fastq files without
#' messing things with your new data. Plus, it will fit nicely in your reproducible workflow!
#'
#' The function is parallelize and can also move the renamed files in a different folder.
#' @rdname rename_fq
#' @export

#' @param change.fq (object, path to a file) Data frame in the global environment or a tab separated file.
#' The dataframe as 2 columns:
#' \code{OLD_FQ} and \code{NEW_FQ} that show the id change you want.
#' If \code{NEW_FQ} column contains a different path, this will also move the files.
#' See example below.

#' @param parallel.core Enable parallel execution using multiple core.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @return The function returns nothing in the global environment. The function
#' just renames the fq files.

#' @examples
#' \dontrun{
#' # change.fq = "fq.new.naming.tsv". This files contains a data frame
#' # with 2 columns separated by a TAB: OLD_FQ\tNEW_FQ
#'
#' # The first line after column headers could be:
#' # 04_process_radtags/124134249.fastq.gz\t04_process_radtags/EST-NEL-ADU-2013-0001.fastq.gz
#' # and so on.
#' # To run:
#' rename_fq(change.fq = "fq.new.naming.tsv", parallel.core = 12)
#'
#'
#' # The function can also rename and move files at the same time.
#' # The first line after column headers could be:
#' # Downloads/124134249.fastq.gz\t04_process_radtags/EST-NEL-ADU-2013-0001.fastq.gz
#' # The second column points to a different and existing folder.
#' # To run, same as above:
#' rename_fq(change.fq = "fq.new.naming.tsv", parallel.core = 12)
#' # the copied and renamed fastq files will be in: 04_process_radtags
#'
#' # Now, how to generate the file fq.new.naming.tsv ? It's actually quite easy within R.
#'
#' # Below, give the path to the folder containing the fastq files and the fq file extension
#' fq.path <- "04_process_radtags"
#' fq.ext <- ".fastq.gz"
#'
#' # Generate a data frame with the column OLD_FQ and NEW_FQ
#' # Here, the NEW_FQ column is identical to OLD_FQ (see next step)
#' fq.new.naming <- tibble::as_tibble(list(
#' OLD_FQ = list.files(path = fq.path, pattern = fq.ext, full.names = TRUE))) %>%
#' dplyr::mutate(NEW_FQ = OLD_FQ)
#'
#'
#' # Save the file
#' readr::write_tsv(x = fq.new.naming, path = "fq.new.naming.tsv")
#'
#' # Then, you could edit the NEW_FQ column, by hand, in MS EXCEL,
#' # or inside R using dplyr package and mutate function.
#' # All this can be done with very different tools. The stringi package can
#' become very handy to parse the OLD_FQ column, etc.
#'
#' fq.new.naming <- fq.new.naming %>%
#' dplyr::mutate(NEW_FQ = bla bla bla)
#'
#' fq.new.naming <- fq.new.naming %>%
#' dplyr::mutate(NEW_FQ = stringi::stri_replace_all_fixed(
#' str = OLD_FQ,
#' pattern = "blablabla",
#' replacement = "newblablabla",
#' vectorize_all = FALSE))
#'
#' }


rename_fq <- function(change.fq, parallel.core = parallel::detectCores() - 1) {
  opt.change <- getOption("width")
  options(width = 70)
  cat("#######################################################################\n")
  cat("########################## stackr::rename_fq ##########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()


  # Checking for missing and/or default arguments-------------------------------
  if (missing(change.fq)) stop("Input file argument is missing")

  if (is.vector(change.fq)) {# for file in the working directory
    change.fq <- readr::read_tsv(file = change.fq, col_types = readr::cols(.default = readr::col_character()))
  }

  if (!tibble::has_name(change.fq, "OLD_FQ") || !tibble::has_name(change.fq, "NEW_FQ")) {
    stop("Verify data frame column naming")
  }

  message("Renaming ", nrow(change.fq), " fastq files...")

  rename_one_fq <- function(x, change.fq) {
    change.fq <- dplyr::filter(change.fq, OLD_FQ == x)
    file.rename(from = change.fq$OLD_FQ, to = change.fq$NEW_FQ)
  }

  suppressWarnings(
    .stackr_parallel_mc(
      X = change.fq$OLD_FQ,
      FUN = rename_one_fq,
      mc.cores = parallel.core,
      change.fq = change.fq))

  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  options(width = opt.change)
}#End rename_fq
