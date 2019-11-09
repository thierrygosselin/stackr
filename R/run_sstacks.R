#' @name run_sstacks
#' @title Run STACKS sstacks module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}
#' module inside R!
#' Inside the \code{P} folder (where the \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks} and the \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks} files are),
#' you should have:
#' \itemize{
#'   \item \strong{3 Catalog files:} the files created in
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks}
#' and usually looking like this:
#' \code{batch_1.catalog.alleles.tsv.gz,
#' batch_1.catalog.snps.tsv.gz,
#' batch_1.catalog.tags.tsv.gz}
#'   \item \strong{3 files for each samples:} The sample name is the prefix of
#'   the files ending with:
#' \code{.alleles.tsv.gz, .snps.tsv.gz, .tags.tsv.gz}.
#' Those files are created in the
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' module.
#' }

#' @param P (character) Path to the directory containing STACKS files.
#' Contrary to \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks},
#' this argument can be use with \code{sample.list} (the \code{s} in \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks})
#' Default: \code{P = "06_ustacks_cstacks_sstacks"}.

#' @param M (character) Path to a population map file from which to take sample
#' names (this argument won't work if \code{sample.list}).
#' Advice: don't modify the default (please see details).
#' Default: \code{M = NULL}.

#' @param sample.list This is for the \code{s} option in
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}.
#' \code{s: Filename prefix from which to load sample loci}.
#' Here, you have 2 choices: 1. you leave empty and let the function use the
#' default:\code{sample.list = NULL} which will scan for the files in the
#' \code{P} folder given above. 2. you supply a character string of the
#' samples. This could come from the \code{INDIVIDUALS_REP} column of the
#' project info file, e.g. \code{sample.list = project.info$INDIVIDUALS_REP}.
#' (Please see details).

#' @param c (character, path) Path to the catalog.
#' Default: \code{c = "06_ustacks_cstacks_sstacks"}

#' @param p (Integer) Enable parallel execution with num_threads threads.
#' Default: \code{p = parallel::detectCores() - 1}

#' @param o output path to write results.
#' Default: \code{o = "06_ustacks_cstacks_sstacks"}

#' @param x Don't verify haplotype of matching locus.
#' Default: \code{x = FALSE}

#' @param v Print program version.
#' Default: \code{v = FALSE}

#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @param disable.gapped Disable gapped alignments between stacks.
#' Default: \code{disable.gapped = FALSE} (use gapped alignments).

#' @rdname run_sstacks
#' @export
#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}
#' returns a \code{.matches.tsv.gz file for each sample}.
#' If \code{lnl_dist = TRUE}, the function will also return a
#' summary of catalog loci log-likelihood.


#' @details \strong{Computer or server problem during the sstacks ?}
#' Just launch again the same command, the function will start again, but only
#' with the unmatched samples!
#'
#' \strong{Some argument can't be use together:}
#'
#' Don't use these arguments together: \code{M} and \code{sample.list/c}.
#' You can't add samples to an existing catalog using a population map.
#' If you want to match new samples to an existing catalog, very easy,
#' just add the samples in the folder and let \code{run_sstacks} function find
#' the samples that were not matched to the catalog.

#' @examples
#' \dontrun{
#' # The simplest form of the function when using the stackr workflow:
#' run_sstacks()
#' # that's it !
#'
#'
#' Now if you have your own workflow folders, etc. Enter them like this:
#' sstacks <- run_sstacks (P = "/my/input/path", p = 32, b = 2,
#' sample.list = c("ind1", "ind2", "..."), o = "/my/output/path",
#' x = FALSE)
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

# sstacks ----------------------------------------------------------------------
run_sstacks <- function(
  P = "06_ustacks_cstacks_sstacks",
  M = NULL,
  sample.list = NULL,
  c = "06_ustacks_cstacks_sstacks",
  p = parallel::detectCores() - 1,
  o = "06_ustacks_cstacks_sstacks",
  x = FALSE,
  disable.gapped = FALSE,
  # lnl_dist = TRUE,
  v = FALSE,
  h = FALSE
) {

  cat("#######################################################################\n")
  cat("######################## stackr::run_sstacks ##########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  res <- list() # return results in this list

  # Check directory ------------------------------------------------------------
  if (!dir.exists(P)) dir.create(P)
  if (!dir.exists("09_log_files")) dir.create("09_log_files")

  # sstacks options ------------------------------------------------------------

  # b: MySQL ID of this batch.
  # if (b == "guess") {
  #   b <- ""
  # } else {
  #   b <- stringi::stri_join("-b ", b)
  # }

  # p: enable parallel execution with num_threads threads.
  parallel.core <- p # backup to use later
  p <- stringi::stri_join("-p ", p)

  # o: output path to write results.
  ustacks.folder <- o
  o <- stringi::stri_join("-o ", shQuote(o))

  # g: base matching on genomic location, not sequence identity.
  # if (g) {
  #   g <- stringi::stri_join("-g ")
  # } else {
  #   g <- ""
  # }

  # x: don't verify haplotype of matching locus.
  if (x) {
    x <- stringi::stri_join("-x ")
  } else {
    x <- ""
  }

  # v: print program version.
  if (v) {
    v <- stringi::stri_join("-v ")
  } else {
    v <- ""
  }

  # h: display this help messsage.
  if (h) {
    h <- stringi::stri_join("-h ")
  } else {
    h <- ""
  }

  # Gapped assembly options ---------------------------------------------------
  if (disable.gapped) {
    disable.gapped <- stringi::stri_join("--disable-gapped ")
  } else {
    disable.gapped <- ""
  }



  # Pop map --------------------------------------------------------------------
  if (is.null(M)) {
    c <- stringi::stri_join("-c ", c)

    # s: filename prefix from which to load sample loci---------------------------

    # Get the samples in the folder
    samples.in.folder <- tibble::tibble(INDIVIDUALS_REP = list.files(P)) %>%
      dplyr::filter(!grepl("catalog", INDIVIDUALS_REP))

    # Search for those already matched
    samples.matched <- samples.in.folder %>%
      dplyr::filter(grepl("matches", INDIVIDUALS_REP)) %>%
      dplyr::mutate(INDIVIDUALS_REP = stringi::stri_replace_all_fixed(
        str = INDIVIDUALS_REP,
        pattern = ".matches.tsv.gz", replacement = "",
        vectorized_all = FALSE)
      ) %>%
      dplyr::distinct(INDIVIDUALS_REP)

    message(stringi::stri_join("Samples in folder already matched to the catalog: ",
                               length(samples.matched$INDIVIDUALS_REP)))

    # Get the name of samples that need to be match to the catalog
    samples.to.match <- samples.in.folder %>%
      dplyr::filter(grepl("alleles", INDIVIDUALS_REP)) %>%
      dplyr::mutate(INDIVIDUALS_REP = stringi::stri_replace_all_fixed(
        str = INDIVIDUALS_REP,
        pattern = ".alleles.tsv.gz", replacement = "",
        vectorized_all = FALSE)
      ) %>%
      dplyr::distinct(INDIVIDUALS_REP) %>%
      dplyr::filter(!INDIVIDUALS_REP %in% samples.matched$INDIVIDUALS_REP)

    message(stringi::stri_join("Samples in folder not matched to the catalog: ", length(samples.to.match$INDIVIDUALS_REP)))


    if (is.null(sample.list)) {
      s <- stringi::stri_join("-s ", shQuote(stringi::stri_join(P, "/", samples.to.match$INDIVIDUALS_REP)))
      message(stringi::stri_join("Matching ", length(samples.to.match$INDIVIDUALS_REP), " sample(s) to the catalog..."))
    } else {
      sample.list.before <- sample.list
      sample.list <- purrr::keep(.x = sample.list, .p = sample.list %in% samples.to.match$INDIVIDUALS_REP)
      sample.list <- stringi::stri_join(P, "/", sample.list)
      s <- stringi::stri_join("-s ", shQuote(sample.list))

      samples.to.remove <- length(sample.list) - length(sample.list.before)
      if (samples.to.remove != 0) {
        message(stringi::stri_join("Removed ", samples.to.remove, " problematic sample(s) from the samples list to match"))
      }
      message(stringi::stri_join("Matching ", length(sample.list), " sample(s) to the catalog..."))
    }
    # remove unwanted // (better solution ? Please tell me!)
    s <- stringi::stri_replace_all_fixed(
      str = s,
      pattern = "//", replacement = "/", vectorize_all = FALSE)

    M <- ""
    P <- ""
  } else {
    M <- stringi::stri_join("-M ", M)
    s <- ""
    c <- ""
    o <- ""
    message("Automatically turning off sstacks arguments: -s, -c, -o")
  }

  # logs files -----------------------------------------------------------------
  file.date.time <- format(Sys.time(), "%Y%m%d@%H%M")
  log.file <- stringi::stri_join("09_log_files/sstacks_", file.date.time,".log")
  message(stringi::stri_join("For progress, look in the log file: ", log.file))

  # command --------------------------------------------------------------------
  command.arguments <- c(P, M, s, c, p, o, x, disable.gapped, v, h)

  # command
  system2(
    command = "sstacks",
    args = command.arguments,
    # stdout = log.file,
    stderr = log.file
  )
  # Summary sstacks ------------------------------------------------------------
  sum <- stackr::summary_sstacks(
    sstacks.log = log.file,
    verbose = FALSE
  )
  sum <- NULL

  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)

}# end run_sstacks
