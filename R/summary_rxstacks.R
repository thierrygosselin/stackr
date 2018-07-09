#' @name summary_rxstacks
#' @title Summarize STACKS rxstacks files
#' @description This function reads the output of STACKS
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/rxstacks.php}{rxstacks}
#' and summarise the information.

#' @param rxstacks.folder (character). The path to the rxstacks output files.
#' e.g. \code{rxstacks.folder = "07_rxstacks_cstacks_sstacks_populations"}

#' @param strata (optional) Path to a strata file to get
#' the rxstacks results by population. With default, the results is presented
#' overall samples.
#' The strata used here is a tab delimited file with a minimum of 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' Additionnally, if you include the corresponding \code{SQL_ID}
#' column header of your samples,
#' the summary of the rxstacks haplotypes and snps files
#' will include the user friendly sample id (easy to get if you've used the stackr pipeline).
#' The \code{STRATA} column can be any hierarchical grouping.
#' If you have already run
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data,
#' the strata file is similar to a stacks \emph{population map file},
#' make sure you
#' have the required column names (\code{INDIVIDUALS}, \code{STRATA} and optionally \code{SQL_ID}).
#' Default: \code{strata = NULL}.

#' @param parallel.core (integer, optional) The number of core used for parallel
#' execution.
#' Default: \code{parallel::detectCores() - 1}.

#' @param verbose (logical, optional) Make the function a little more chatty during
#' execution.
#' Default: \code{verbose = FALSE}.

#' @param ... Used internally in stackr.

#' @rdname summary_rxstacks
#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr group_by tally ungroup summarise summarise_if mutate filter distinct select bind_rows n_distinct arrange
#' @importFrom readr read_tsv
#' @importFrom stats cor
#' @importFrom tibble add_column
#' @importFrom rlang dots_list

#' @return The function returns a summary (data frame) containing:
#' \enumerate{
#' \item INDIVIDUALS: the sample id
#' \item LOCUS_NUMBER: the number of locus
#' \item BLACKLIST_rxstacks: the number of locus blacklisted by ustacks
#' \item FOR_CATALOG: the number of locus available to generate the catalog
#' \item BLACKLIST_ARTIFACT: the number of artifact genotypes (> 2 alleles, see details)
#' \item FILTERED: the number of locus after artifacts are removed
#' \item HOMOZYGOSITY: the number of homozygous genotypes
#' \item HETEROZYGOSITY: the number of heterozygous genotypes
#' \item MEAN_NUMBER_SNP_LOCUS: the mean number of SNP/locus (excluding the artifact locus)
#' \item MAX_NUMBER_SNP_LOCUS: the max number of SNP/locus observed for this individual (excluding the artifact locus)
#' \item NUMBER_LOCUS_4SNP: the number of locus with 4 or more SNP/locus (excluding the artifact locus)
#' }


#' @examples
#' \dontrun{
#' rxstacks.sum <- stackr::summary_rxstacks(
#' rxstacks.folder = "07_rxstacks_cstacks_sstacks_populations")
#'
#' # To get the number of snp/locus for a specific individual:
#' snp.info <- ustacks.summary %>%
#' dplyr::filter(INDIVIDUALS == "ID2") %>%
#' dplyr::select(INDIVIDUALS, SNP_LOCUS) %>%
#' tidyr::unnest(.)
#'
#' #Similarly, for the blacklisted locus (artifactual):
#' blacklist <- ustacks.summary %>%
#' dplyr::filter(INDIVIDUALS == "ID2") %>%
#' dplyr::select(INDIVIDUALS, BLACKLIST_ARTIFACT) %>%
#' tidyr::unnest(.)
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/rxstacks.php}{rxstacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

summary_rxstacks <- function(
  rxstacks.folder,
  strata = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...) {
  message("\n\n\nThis function will be deprecated in future version, please update your pipeline to stacks v.2.1\n\n\n")
  if (verbose) cat("#######################################################################\n")
  if (verbose) cat("##################### stackr::summary_rxstacks ########################\n")
  if (verbose) cat("#######################################################################\n")

  timing <- proc.time()
  res <- list() # return results in this list
  if (missing(rxstacks.folder)) stop("rxstacks.folder argument required")
  if (!dir.exists("09_log_files")) dir.create("09_log_files")

  # Date and time extension ----------------------------------------------------
  file.date.time <- rlang::dots_list(...)

  if (length(file.date.time) == 0) {
    file.date.time <- stringi::stri_replace_all_fixed(
      str = Sys.time(),
      pattern = " EDT", replacement = "") %>%
      stringi::stri_replace_all_fixed(
        str = .,
        pattern = c("-", " ", ":"),
        replacement = c("", "@", ""),
        vectorize_all = FALSE
      ) %>%
      stringi::stri_sub(str = ., from = 1, to = 13)
  }

  # rxstacks log file ---------------------------------------------------------
  log.file <- list.files(
    path = rxstacks.folder, pattern = "rxstacks.log", full.names = TRUE)
  new.log.file <- stringi::stri_join("09_log_files/rxstacks_correction_summary_",
                                     file.date.time,".log")
  file.rename(from = log.file, to = new.log.file)
  message("\nImporting and summarizing stacks log file:\n", new.log.file)

  res$rxstacks.correction.ind <- readr::read_tsv(
    new.log.file,
    comment = "#",
    col_types = "ciiiiiiiiiiiii",
    col_names = c("INDIVIDUALS", "TOTAL_NUCS", "TOTAL_NUCS_CONVERTED",
                  "UNK_TO_HOM", "UNK_TO_HET", "HOM_TO_UNK", "HET_TO_UNK",
                  "HOM_TO_HET", "HET_TO_HOM","BLACKLIST_CALL_HAPLOTYPE",
                  "BLACKLIST_CONFOUNDED_LOCI", "LNL_FILTERED_LOCI",
                  "PRUNED_HAPLOTYPES_RARE", "PRUNED_HAPLOTYPES_TREE")) %>%
    dplyr::mutate(
      PRUNED_HAPLOTYPES = (PRUNED_HAPLOTYPES_RARE + PRUNED_HAPLOTYPES_TREE),
      BLACKLISTED_TOTAL = BLACKLIST_CALL_HAPLOTYPE + BLACKLIST_CONFOUNDED_LOCI + LNL_FILTERED_LOCI
    )
  res$rxstacks.correction.overall <- dplyr::summarise_if(
    .tbl = res$rxstacks.correction.ind, .predicate = is.integer, .funs = mean) %>%
    tibble::add_column(.data = ., POP_ID = "OVERALL", .before = 1)

  if (!is.null(strata)) {
    # get the number of column in the strata file
    n.col <- ncol(suppressWarnings(suppressMessages(readr::read_tsv(file = strata, n_max = 1))))

    if (n.col < 3) {
      strata.info <- readr::read_tsv(file = strata, col_types = "cc") %>%
        dplyr::rename(POP_ID = STRATA)
    } else {
      strata.info <- readr::read_tsv(file = strata, col_types = "cci") %>%
        dplyr::rename(POP_ID = STRATA)
    }
    n.col <- NULL
    rxstacks.correction.ind.pop <- suppressWarnings(dplyr::left_join(
      res$rxstacks.correction.ind, dplyr::select(strata.info, INDIVIDUALS, POP_ID), by = "INDIVIDUALS") %>%
        dplyr::group_by(POP_ID))

    res$rxstacks.correction.overall <- dplyr::summarise_if(
      .tbl = rxstacks.correction.ind.pop, .predicate = is.integer, .funs = mean) %>%
      dplyr::bind_rows(res$rxstacks.correction.overall)
    rxstacks.correction.ind.pop <- NULL
  }

  res$rxstacks.correction.overall.proportion <- res$rxstacks.correction.overall %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::mutate(
      NUCS_CONVERTED = TOTAL_NUCS_CONVERTED / TOTAL_NUCS,
      BLACKLISTED = BLACKLISTED_TOTAL / TOTAL_NUCS,
      PRUNED_HAPLOTYPES = PRUNED_HAPLOTYPES / TOTAL_NUCS,
      NUC_NOT_CONVERTED = 1 - (NUCS_CONVERTED + BLACKLISTED+PRUNED_HAPLOTYPES),
      UNK_TO_HOM = UNK_TO_HOM / TOTAL_NUCS_CONVERTED,
      UNK_TO_HET = UNK_TO_HET / TOTAL_NUCS_CONVERTED,
      HOM_TO_UNK = HOM_TO_UNK / TOTAL_NUCS_CONVERTED,
      BLACKLIST_CALL_HAPLOTYPE = BLACKLIST_CALL_HAPLOTYPE / BLACKLISTED_TOTAL,
      BLACKLIST_CONFOUNDED_LOCI = BLACKLIST_CONFOUNDED_LOCI / BLACKLISTED_TOTAL,
      LNL_FILTERED_LOCI = LNL_FILTERED_LOCI / BLACKLISTED_TOTAL
    ) %>%
    dplyr::select(POP_ID, NUC_NOT_CONVERTED, NUCS_CONVERTED, BLACKLISTED, PRUNED_HAPLOTYPES, UNK_TO_HOM, UNK_TO_HET, HOM_TO_UNK, BLACKLIST_CALL_HAPLOTYPE, BLACKLIST_CONFOUNDED_LOCI, LNL_FILTERED_LOCI)

  # Percentages to make the figure

  fig.data <- tidyr::gather(res$rxstacks.correction.overall.proportion, RXSTACKS, PROPORTION, -POP_ID) %>%
    dplyr::mutate(
      GROUP = dplyr::if_else(
        RXSTACKS == "NUCS_CONVERTED" | RXSTACKS == "BLACKLISTED" | RXSTACKS == "PRUNED_HAPLOTYPES" | RXSTACKS == "NUC_NOT_CONVERTED", "correction",
        dplyr::if_else(RXSTACKS == "UNK_TO_HOM" | RXSTACKS == "UNK_TO_HET" | RXSTACKS == "HOM_TO_UNK", "converted", "blacklisted")),
      GROUP = factor(GROUP, levels = c("correction", "converted", "blacklisted"), ordered = TRUE)
    )

  res$rxstacks.correction.filled.bar.plot <- ggplot2::ggplot(fig.data, ggplot2::aes(x = POP_ID, y = PROPORTION, fill = RXSTACKS)) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(y = "Proportion") +
    ggplot2::labs(x = "Sampling sites") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0,size = 12, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")) +
    ggplot2::facet_wrap(~GROUP, nrow = 1, ncol = 3)


  res$rxstacks.correction.bar.plot <- ggplot2::ggplot(fig.data, ggplot2::aes(x = POP_ID, y = PROPORTION)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(y = "Proportion") +
    ggplot2::labs(x = "Sampling sites") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0,size = 12, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")) +
    ggplot2::facet_wrap(RXSTACKS~GROUP, scales = "free", nrow = 5, ncol = 2)
  fig.data <- NULL


  # rxstacks haplotype log file ------------------------------------------------
  haplo.log.file <- list.files(
    path = rxstacks.folder, pattern = "rxstacks.haplotypes.log", full.names = TRUE)
  if (length(haplo.log.file) > 0) {
    new.haplo.log.file <- stringi::stri_join("09_log_files/rxstacks.haplotypes_",
                                             file.date.time, ".log")
    file.rename(from = haplo.log.file, to = new.haplo.log.file)
    message("\nImporting and moving stacks rxstacks haplotypes log file:\n", new.haplo.log.file)
    message("Summarizing by individuals...")
    res$rxstacks.haplo <- readr::read_tsv(
      file = new.haplo.log.file,
      comment = "#",
      col_types = "iiiccccc",
      col_names = c("CATALOG_LOCUS", "SQL_ID", "SAMPLE_LOCUS", "SAMPLE_HAPLOTYPE", "CATALOG_HAPLOTYPE", "CORRECTED_SAMPLE_HAPLOTYPE", "CORRECTED_CATALOG_HAPLOTYPE", "ALGORITHM")) %>%
      dplyr::group_by(SQL_ID, ALGORITHM) %>%
      dplyr::tally(.) %>%
      dplyr::group_by(SQL_ID) %>%
      tidyr::spread(data = ., key = ALGORITHM, value = n) %>%
      dplyr::rename(RXSTACKS_MST = mst, RXSTACKS_RARE_STEP_1 = rare_step_1)

    if (tibble::has_name(strata.info, "SQL_ID")) {
      res$rxstacks.haplo <- dplyr::left_join(res$rxstacks.haplo, strata.info, by = "SQL_ID") %>%
        dplyr::select(INDIVIDUALS, SQL_ID, POP_ID, RXSTACKS_MST, RXSTACKS_RARE_STEP_1)
    }
  }

  # rxstacks snps log file ------------------------------------------------
  snps.log.file <- list.files(
    path = rxstacks.folder, pattern = "rxstacks.snps.log", full.names = TRUE)
  if (length(snps.log.file) > 0) {
    new.snps.log.file <- stringi::stri_join("09_log_files/rxstacks.snps_",
                                            file.date.time, ".log")
    file.rename(from = snps.log.file, to = new.snps.log.file)
    message("\nImporting and moving stacks rxstacks snps log file:\n", new.snps.log.file)
    message("Summarizing by individuals...")
    res$rxstacks.snps <- readr::read_tsv(
      file = new.snps.log.file,
      comment = "#",
      col_types = "iiicc",
      col_names = c("SQL_ID", "LOCUS_ID", "SNP_COL", "ORIG_VALUE", "CORR_VALUE")) %>%
      dplyr::group_by(SQL_ID, ORIG_VALUE, CORR_VALUE) %>%
      dplyr::tally(.) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(FROM_TO = stringi::stri_join(ORIG_VALUE, "_", CORR_VALUE)) %>%
      dplyr::select(-ORIG_VALUE, -CORR_VALUE) %>%
      dplyr::group_by(SQL_ID) %>%
      tidyr::spread(data = ., key = FROM_TO, value = n)

    if (tibble::has_name(strata.info, "SQL_ID")) {
      res$rxstacks.snps <- dplyr::left_join(res$rxstacks.snps, strata.info, by = "SQL_ID") %>%
        dplyr::select(INDIVIDUALS, SQL_ID, POP_ID, dplyr::everything(.))
    }
  }
  timing <- proc.time() - timing
  if (verbose) message("\nComputation time: ", round(timing[[3]]), " sec")
  if (verbose) cat("############################## completed ##############################\n")
  return(res)
}#End summary_rxstacks
