#' @name run_cstacks
#' @title Run STACKS cstacks module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks}
#' module inside R!

#' @param b database/batch ID for this catalog. Advice: don't modify.
#' Default: \code{b = 1}.

#' @param P path to the directory containing STACKS files.
#' Default: \code{P = "06_ustacks_cstacks_sstacks"}.
#' Inside the folder \code{06_ustacks_cstacks_sstacks}, you should have:
#' \itemize{
#'   \item \strong{4 files for each samples:} The sample name is the prefix of
#'   the files ending with:
#' \code{.alleles.tsv.gz, .models.tsv.gz, .snps.tsv.gz, .tags.tsv.gz}.
#' Those files are created in the
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' module.
#' }

#' @param M path to a population map file (Required when P is used).
#' Default: \code{M = "06_ustacks_cstacks_sstacks/population.map.catalog.tsv"}.

#' @param g base catalog construction on alignment position, not sequence identity.
#' Advice: don't modify.
#' Default: \code{g = FALSE}.

#' @param n number of mismatches allowed between sample loci when build the catalog.
#' Default: \code{n = 1}

#' @param p enable parallel execution with num_threads threads.
#' Default: \code{p = parallel::detectCores() - 1}

#' @param catalog.path This is for the "Catalog editing" part in cstacks where
#' you can provide the path to an existing catalog.
#' cstacks will add data to this existing catalog.
#' With default: \code{catalog.path = NULL} or with a supplied path, the function
#' The function scan automatically for the presence of a catalog inside the input folder.
#' If none is found, a new catalog is created.
#' If your catalog is not in the input folder, supply a path here.
#' e.g. \code{catalog.path = ~/catalog_folder}, the catalog files are inside the
#' P folder along the samples files and detected automatically.
#' If a catalog is detected in the input folder,
#' the samples in the \code{sample.list} argument
#' will be added in this catalog. The catalog is made of 3 files:
#' \code{batch_1.catalog.alleles.tsv.gz,
#' batch_1.catalog.snps.tsv.gz,
#' batch_1.catalog.tags.tsv.gz}

#' @param gapped Gapped assembly options: do you want to preform
#' gapped alignments between stacks.
#' Default: \code{gapped = TRUE}

#' @param max_gaps The number of gaps allowed between stacks before merging.
#' Default: \code{max_gaps = 2}

#' @param min_aln_len The minimum length of aligned sequence in a gapped
#' alignment.
#' Default: \code{min_aln_len = 0.8}

#' @param m Include tags in the catalog that match to more than one entry.
#' Advice: don't modify.
#' Default: \code{m = FALSE}

#' @param k_len Specify k-mer size for matching between between catalog loci
#' (automatically calculated by default).
#' Advice: don't modify.
#' Default: \code{k_len = NULL}

#' @param report_mmatches Report query loci that match more than one catalog locus.
#' Advice: don't modify.
#' Default: \code{report_mmatches = FALSE}

#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @rdname run_cstacks
#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed

#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}
#' returns a \code{.matches.tsv.gz file for each sample}

#' @details \strong{Computer or server problem during the cstacks ?} Look
#' in the log file to see which individuals remains to be included. Create a
#' new list of individuals to include and use the catalog.path argument to point
#' to the catalog created before the problem.

#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' run_cstacks()
#' # that's it ! Now if you have your own workflow folders, etc. See below.
#' Next example, let say you only want to include 10 individuals/pop and
#' include in the catalog samples with more than 2000000 reads. With the project
#' info file in the global environment:
#' library(tidyverse)
#' individuals.catalog <- project.info.file) %>%
#' filter(RETAINED > 2000000) %>%
#' group_by(POP_ID) %>%
#' sample_n(size = 10, replace = FALSE) %>%
#' ungroup %>%
#' arrange(desc(RETAINED)) %>%
#' distinct(INDIVIDUALS_REP, POP_ID)
#' # Write file to disk
#' readr::write_tsv(x = individuals.catalog,
#' path = "06_ustacks_cstacks_sstacks/population.map.catalog.tsv")
#' # The next line will give you the list of individuals to include
#' individuals.catalog <- individuals.catalog$INDIVIDUALS_REP
#'
#' # To keep your info file updated with this information:
#' project.info.file <- project.info.file %>%
#' mutate(CATALOG = if_else(INDIVIDUALS_REP %in% individuals.catalog,
#' true = "catalog", false = "not_catalog")
#' )
#' write_tsv(project.info.file, "project.info.catalog.tsv")
#'
#' # Then run the command this way:
#' run_cstacks (
#' P = "06_ustacks_cstacks_sstacks",
#' catalog.path = NULL,
#' b = 1,
#' g = FALSE,
#' m = FALSE,
#' n = 1,
#' p = 32,
#' h = FALSE,
#' gapped = TRUE, max_gaps = 2, min_aln_len = 0.8,
#' k_len = NULL, report_mmatches = FALSE
#' )
#' }

#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

run_cstacks <- function(
  b = 1,
  P = "06_ustacks_cstacks_sstacks",
  M = "06_ustacks_cstacks_sstacks/population.map.catalog.tsv",
  g = FALSE,
  n = 1,
  p = parallel::detectCores() - 1,
  catalog.path = NULL,
  gapped = TRUE, max_gaps = 2, min_aln_len = 0.8,
  m = FALSE, k_len = NULL, report_mmatches = FALSE,
  h = FALSE
  # , transfer.s3 = FALSE,
  # from.folder = NULL, destination.folder = NULL,
) {

  cat("#######################################################################\n")
  cat("######################## stackr::run_cstacks ##########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists(P)) stop("Missing P directory")
  if (!dir.exists("09_log_files")) dir.create("09_log_files")

  # Catalog editing ------------------------------------------------------------

  if (is.null(catalog.path)) { # no catalog path, searching in the input path...
    old.catalog <- list.files(path = P, pattern = "batch_")
    detect.rxstacks.lnls <- list.files(path = P, pattern = "rxstacks_lnls")
    if (length(detect.rxstacks.lnls) == 1) {
      old.catalog <- purrr::discard(.x = old.catalog, .p = old.catalog %in% detect.rxstacks.lnls)
    }
    if (length(old.catalog) > 0 & length(old.catalog) == 3) {
      message("Found a catalog in the input folder, using files: ")
      message(stringi::stri_join(old.catalog, "\n"))

      catalog.path <- stringi::stri_replace_all_fixed(
        str = old.catalog[1],
        pattern = ".catalog.alleles.tsv.gz",
        replacement = "",
        vectorize_all = FALSE
      )
      catalog.path <- stringi::stri_join(P, "/", catalog.path)
      catalog.path <- stringi::stri_join("--catalog ", shQuote(catalog.path))
    }
    if (length(old.catalog) > 0 & length(old.catalog) < 3) {
      stop("Incomplete catalog, 3 files are required, see argument documentation")
    }

    if (length(old.catalog) == 0) {
      message("Builing catalog for the first time")
      catalog.path <- ""
    }
  } else {
    old.catalog <- list.files(path = catalog.path, pattern = "batch_")
    if (length(old.catalog) > 0 & length(old.catalog) == 3) {
      message("Found the catalog in the catalog path using files: ")
      message(stringi::stri_join(old.catalog, "\n"))

      catalog.path <- stringi::stri_replace_all_fixed(
        str = old.catalog[1],
        pattern = ".catalog.alleles.tsv.gz",
        replacement = "",
        vectorize_all = FALSE
      )
      catalog.path <- stringi::stri_join(P, "/", catalog.path)
      catalog.path <- stringi::stri_join("--catalog ", shQuote(catalog.path))
    }
    if (length(old.catalog) > 0 & length(old.catalog) < 3) {
      stop("Incomplete catalog, 3 files are required, see argument documentation")
    }

    if (length(old.catalog) == 0) {
      message("Builing catalog for the first time")
      catalog.path <- ""
    }
  }


  # cstacks options ------------------------------------------------------------
  b <- stringi::stri_join("-b ", b)
  P <- stringi::stri_join("-P ", P)
  M <- stringi::stri_join("-M ", M)

  if (g) {
    g <- stringi::stri_join("-g ")
  } else {
    g <- ""
  }

  n <- stringi::stri_join("-n ", n)
  p <- stringi::stri_join("-p ", p)

  # gapped assembly options ---------------------------------------------------
  if (gapped) {
    gapped <- stringi::stri_join("--gapped ")
  } else {
    gapped <- ""
  }

  max_gaps <- stringi::stri_join("--max_gaps ", max_gaps)
  min_aln_len <- stringi::stri_join("--min_aln_len ", min_aln_len)


  # Advanced options -----------------------------------------------------------

  if (m) {
    m <- stringi::stri_join("-m ")
  } else {
    m <- ""
  }

  if (is.null(k_len)) {
    k_len <- ""
  } else {
    k_len <- stringi::stri_join("--k_len ", k_len)
  }

  if (report_mmatches) {
    report_mmatches <- stringi::stri_join("--report_mmatches ")
  } else {
    report_mmatches <- ""
  }

  # Help  ------------------------------------------------------------------------
  if (h) {
    h <- stringi::stri_join("-h ")
  } else {
    h <- ""
  }


  # logs files -----------------------------------------------------------------
  file.date.time <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
  file.date.time <- stringi::stri_replace_all_fixed(
    file.date.time,
    pattern = c("-", " ", ":"),
    replacement = c("", "@", ""),
    vectorize_all = FALSE
  )
  file.date.time <- stringi::stri_sub(file.date.time, from = 1, to = 13)

  cstacks.log.file <- stringi::stri_join("09_log_files/cstacks_", file.date.time,".log")
  message(stringi::stri_join("For progress, look in the log file: ", cstacks.log.file))


  # command args ---------------------------------------------------------------
  command.arguments <- paste(
    b, P, M, g, n, p, catalog.path,
    gapped, max_gaps, min_aln_len,
    m, k_len, report_mmatches,
    h
  )

  # command
  system2(command = "cstacks", args = command.arguments, stderr = cstacks.log.file)

  # # transfer back to s3
  # if (transfer.s3) {
  #   cstacks.files.to.s3 <- list.files(path = sample.list.path, pattern = individual, full.names = FALSE)
  #   purrr::walk(.x = cstacks.files.to.s3, .f = copy_s3, from.folder = from.folder, destination.folder = destination.folder)
  # }


timing <- proc.time() - timing
message("\nComputation time: ", round(timing[[3]]), " sec")
cat("############################## completed ##############################\n")
}# end run_cstacks
