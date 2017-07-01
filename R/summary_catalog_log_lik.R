#' @name summary_catalog_log_lik
#' @title Summary of catalog log likelihood

#' @description This function reads inside the output of
#' STACKS ustasks-cstacks-sstacks folder
#' and and summarize the sstacks output files \code{.matches} to get
#' the mean and median log likelihood of catalog loci.

#' @param matches.folder (logical). Default: \code{mismatch.testing = FALSE}.

#' @param parallel.core (optional) The number of core used for parallel
#' execution.
#' Default: \code{parallel::detectCores() - 1}.

#' @param verbose (optional) Make the function a little more chatty during
#' execution.
#' Default: \code{verbose = FALSE}.

#' @rdname summary_catalog_log_lik
#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr group_by tally ungroup summarise summarise_if mutate filter distinct select bind_rows n_distinct arrange ntile
#' @importFrom readr read_tsv
#' @importFrom stats cor median
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot geom_histogram labs facet_wrap theme element_text


#' @return The function returns a summary (data frame) containing:
#' \enumerate{
#' \item INDIVIDUALS: the sample id
#' \item LOCUS_NUMBER: the number of locus
#' \item BLACKLIST_USTACKS: the number of locus blacklisted by ustacks
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
#' catalog.log.lik <- stackr::summary_catalog_log_lik(
#' matches.folder = "06_ustacks_cstacks_sstacks", verbose = TRUE)
#' # To get the summary after correction module STACKS rxstacks, use the built-in
#' feature of \code{run_rxstacks} or run this function with rxstacks outout folder selected.
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

summary_catalog_log_lik <- function(
  matches.folder,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE) {
  if (verbose) cat("#######################################################################\n")
  if (verbose) cat("################# stackr::summary_catalog_log_lik #####################\n")
  if (verbose) cat("#######################################################################\n")
  timing <- proc.time()

  res <- list() # return results in this list

  if (missing(matches.folder)) stop("matches.folder argument required")
  message("Printing the distribution of mean log likelihoods for catalog loci: selected")
  # message("After visualizing the distribution:")
  # message("    1. select the best threshold")
  # message("    2. use function run_rxstacks with the appropriate filter thresholds.")

  # Import matches
  matches.files <- list.files(
    path = matches.folder, pattern = "matches", full.names = TRUE)

  n.matches <- length(matches.files)

  if (n.matches == 0) stop("Missing matches files: generate the files for each samples with run_sstacks or sstacks in command line mode")
  message("Reading ", n.matches, " matches files\n    for the distribution of mean log likelihoods of catalog loci...")

  opt.change <- getOption("width")
  options(width = 70)

  # functions
  read_matches <- function(matches.files) {
    matches <- readr::read_tsv(
      file = matches.files,
      col_names = c("SQL_ID", "BATCH_ID", "CATALOG_ID", "SAMPLE_ID", "LOCUS_ID", "HAPLOTYPE", "STACK_DEPTH", "LOG_LIKELIHOOD"),
      col_types = "iiiiicid", na = "-",
      comment = "#")
    return(matches)
  }
  summarize_likelihood <- function(x) {
    lik <- dplyr::group_by(x, CATALOG_ID) %>%
      dplyr::summarise(
        MEAN = mean(LOG_LIKELIHOOD),
        MEDIAN = stats::median(LOG_LIKELIHOOD))
    return(lik)
  }

  res$log.likelihood.summary <- list()
  res$log.likelihood.summary <- .stackr_parallel(
    X = matches.files,
    FUN = read_matches,
    mc.cores = parallel.core
  ) %>%
    dplyr::bind_rows(.) %>%
    dplyr::distinct(CATALOG_ID, SAMPLE_ID, LOG_LIKELIHOOD) %>%
    dplyr::arrange(CATALOG_ID)

  message("Summarizing the log-likelihood....")
  res$log.likelihood.summary <- res$log.likelihood.summary %>%
    dplyr::left_join(
      dplyr::distinct(res$log.likelihood.summary, CATALOG_ID) %>%
        dplyr::mutate(
          SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3))
      , by = "CATALOG_ID") %>%
    split(x = ., f = .$SPLIT_VEC) %>%
    .stackr_parallel(
      X = .,
      FUN = summarize_likelihood,
      mc.cores = parallel.core
    ) %>%
    dplyr::bind_rows(.)

  options(width = opt.change) #restore the option width

  overall.mean <- round(mean(res$log.likelihood.summary$MEAN), 2)
  overall.median <- round(stats::median(res$log.likelihood.summary$MEDIAN), 2)

  res$log.likelihood.fig <- tidyr::gather(
    data = res$log.likelihood.summary, key = STATS, value = VALUE, -CATALOG_ID) %>%
    ggplot2::ggplot(data = ., ggplot2::aes(x = VALUE)) +
    ggplot2::geom_histogram() +
    ggplot2::labs(x = "Log likelihoods of catalog loci") +
    ggplot2::labs(y = "Number of catalog loci") +
    ggplot2::facet_wrap(~STATS, nrow = 1, ncol = 2) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold"))

  message("Distribution plot available inside the object in the Global environment")
  message("Overall mean log likelihood: ", overall.mean)
  message("Overall median log likelihood: ", overall.median)

  timing <- proc.time() - timing
  if (verbose) message("\nComputation time: ", round(timing[[3]]), " sec")
  if (verbose) cat("############################## completed ##############################\n")
  return(res)
}#End summary_catalog_log_lik
