#' @name run_radproc
#' @title Run RADproc
#' @description Runs \href{https://github.com/beiko-lab/RADProc}{RADproc}.
#' The approach \emph{replaces ustacks and cstacks} steps. Read the paper
#' for more information.
#' Read \href{http://thierrygosselin.github.io/stackr/articles/stackr.html}{stackr} vignette.

#' @param file.type (character) Input file Type.
#' Supported types: fasta, fastq, gzfasta, or gzfastq.
#' Default: \code{file.type = "gzfastq"}.

#' @param f (path) Input file path. Usually,
#' the stacks process_radtags output folder.
#' Default: \code{f = "04_process_radtags"}.

#' @param o (path) Output path to write results.
#' Default: \code{o = "06_ustacks_2_gstacks"}.

#' @param a (logical) Enable parameter sweep mode.
#' Default: \code{a = FALSE}.
#' Please use \code{run_ustacks} with mismatch testing.

#' @param M Maximum distance (in nucleotides) allowed between stacks to form
#' network.
#' Default: \code{M = 2}.

#' @param m Minimum depth of coverage.
#' Default: \code{m = 3}.

#' @param n Maximum distance (in nucleotides) allowed between catalog loci to merge
#' Default: \code{n = 2}.

#' @param parallel.core (integer) Enable parallel execution with num_threads
#' threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param x (integer) Maximum number stacks per locus.
#' Default: \code{x = 3}.

#' @param S (percentage) Minimum sample percentage.
#' Default: \code{S = 30}.

#' @param D (percentage) Minimum average coverage depth.
#' Default: \code{D = 7}.


#' @param cmd.path (character, path) Provide the FULL path to RADProc
#' program. See details on how to install RADProc in
#' \href{http://thierrygosselin.github.io/stackr/articles/stackr.html#ustackscstacks-alternatives}{stackr vignette}.
#' Default: \code{cmd.path = "/usr/local/bin/RADProc"}.

#' @rdname run_radproc
#' @export


#' @return Returns 3 files per samples: \code{.snps.tsv}, \code{.tags.tsv},
#' \code{.alleles.tsv}. Also returns 3 catalog files: \code{catalog.snps.tsv,
#' catalog.tags.tsv, catalog.alleles.tsv}.

#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' u <- stackr::run_radproc() # that's it !
#' }

#' @seealso
#' \href{https://github.com/beiko-lab/RADProc}{RADproc}.

#' @references Ravindran, P., Bentzen, P., Bradbury, I., Beiko, R. (2019).
#' RADProc: A computationally efficient de novo locus assembler for population
#' studies using RADseq data.
#' Molecular Ecology Resources 19(1), 272-282.
#' https://dx.doi.org/10.1111/1755-0998.12954

run_radproc <- function(
  file.type = "gzfastq",
  f = "04_process_radtags",
  o = "06_ustacks_2_gstacks",
  a = FALSE,
  M = 2,
  m = 3,
  n = 2,
  parallel.core = parallel::detectCores() - 1,
  x = 3,
  S = 2,
  D = 7,
  cmd.path = "/usr/local/bin/RADProc"
  ) {

  # testing

  cat("#######################################################################\n")
  cat("######################## stackr::run_radproc ##########################\n")
  cat("#######################################################################\n")

  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists(f)) rlang::abort("Input file path (f) not found... ")
  if (!dir.exists("09_log_files")) dir.create("09_log_files")
  if (!dir.exists("08_stacks_results")) dir.create("08_stacks_results")
  if (!dir.exists(o)) dir.create(o)



  # check that RADProc is installed -----------------------0---------------------
  if (!file.exists(cmd.path)) stop("Path to RADProc is not valid")


  # file data and time ---------------------------------------------------------
  file.date.time <- format(Sys.time(), "%Y%m%d@%H%M")

  # logs file ------------------------------------------------------------------
  radproc.log.file <- stringi::stri_join("09_log_files/radproc_", file.date.time,".log")
  message("For progress, look in the log file:\n", radproc.log.file)

  # RADProc arguments ----------------------------------------------------------

  # file type
  file.type <- stringi::stri_join("-t ", shQuote(file.type))
  f <- stringi::stri_join("-f ", shQuote(f))

  # output dir
  o.bk <- o
  o <- stringi::stri_join("-o ", shQuote(o))

  if (a) {
    a <- "-a "
  } else {
    a <- ""
  }
  M <- stringi::stri_join("-M ", M)
  m <- stringi::stri_join("-m ", m)
  n <- stringi::stri_join("-n ", n)
  parallel.core.bk <- parallel.core
  parallel.core <- stringi::stri_join("-p ", parallel.core)
  x <- stringi::stri_join("-x ", x)
  S <- stringi::stri_join("-S ", S)
  D <- stringi::stri_join("-D ", D)


  # command args ---------------------------------------------------------------
  command.arguments <- paste(file.type, f, o, a, M, m, n, parallel.core, x, S, D)

  # run command ----------------------------------------------------------------
  system2(command = "RADProc", args = command.arguments, stderr = radproc.log.file)

  stats.filename <- file.path("08_stacks_results", stringi::stri_join("De_novo_stats", file.date.time, sep = "_"))
  file.rename(from = file.path(o.bk, "De_novo_stats"), to = stats.filename)


  # We have to rename the files and move them out of the RADProc folder
  # detect RADProc folder created

  rp.folder <- readr::read_tsv(file = stats.filename, comment = "#", col_types = "ccddd", n_max = 3) %>%
    dplyr::distinct(Parameters) %$%
    Parameters

  move.tibble <- tibble::tibble(
    OLD_FQ = list.files(path = file.path(o.bk, rp.folder), full.names = TRUE),
    SHORT_NAMES = list.files(path = file.path(o.bk, rp.folder), full.names = FALSE)
  ) %>%
    dplyr::mutate(
      NEW_FQ = file.path(o.bk, SHORT_NAMES),
      SHORT_NAMES = NULL
    )

  message("Moved ", nrow(move.tibble), " files from folder ", rp.folder, " to output folder ", basename(o.bk))

  stackr::rename_fq(change.fq = move.tibble, parallel.core = parallel.core.bk, verbose = FALSE)
  file.remove(file.path(o.bk, rp.folder))

  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("########################## RADProc completed ##########################\n")
  return(radproc.log.file)
}# end run_radproc

