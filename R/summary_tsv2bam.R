#' @name summary_tsv2bam
#' @title Summary of catalog log likelihood

#' @description This function reads the log file of
#' \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{tsv2bam}
#' module and produce a summary data frame and plots with loci distribution.
#' This function is run automatically with
#' \pkg{stackr} \code{\link{run_tsv2bam}} function.

#' @param tsv2bam.output (character, path). Path to tsv2bam output files.
#' Default: \code{P = "06_ustacks_cstacks_sstacks"}.

#' @param verbose (optional) Make the function a little more chatty during
#' execution.
#' Default: \code{verbose = TRUE}.

#' @rdname summary_tsv2bam
#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr group_by tally ungroup summarise summarise_if mutate filter distinct select bind_rows n_distinct arrange ntile
#' @importFrom readr read_tsv
#' @importFrom stats cor median
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot geom_histogram labs facet_wrap theme element_text


#' @return The function returns a list with the number of individuals, the batch ID number,
#' a summary data frame and a plot containing:
#' \enumerate{
#' \item INDIVIDUALS: the sample id
#' \item ALL_LOCUS: the total number of locus for the individual (shown in subplot A)
#' \item LOCUS: the number of locus with a one-to-one relationship (shown in subplot B)
#' with the catalog
#' \item MATCH_PERCENT: the percentage of locus with a one-to-one relationship
#' with the catalog (shown in subplot C)
#' }


#' @examples
#' \dontrun{
#' # with all the files in default folders, very simple:
#' sum <- stackr::summary_tsv2bam()
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{stacks Version 2.0Beta6}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

summary_tsv2bam <- function(tsv2bam.output = "06_ustacks_cstacks_sstacks",
  verbose = TRUE) {
  if (verbose) cat("#######################################################################\n")
  if (verbose) cat("##################### stackr::summary_tsv2bam #########################\n")
  if (verbose) cat("#######################################################################\n")
  timing <- proc.time()
  # opt.change <- getOption("width")
  # options(width = 70)

  res <- list() # return results in this list

  # Import log file
  if (verbose) message("Importing log file")

  log.file <- list.files(
    path = tsv2bam.output, pattern = "tsv2bam.log", full.names = TRUE)

  if (verbose) message("Log file located: ", log.file)
  tsv2bam.log <- readr::read_tsv(file = log.file, col_names = "INPUT",
                                 col_types = "c")

  n.ind.match <-  which(stringi::stri_detect_fixed(str = tsv2bam.log$INPUT,
                                                   pattern = "Num. samples: "))
  res$n.ind <- readr::read_lines(file = log.file, skip = n.ind.match - 1, n_max = 1) %>%
    stringi::stri_extract_all_regex(str = ., pattern = "[1-9]+", simplify = TRUE) %>%
    as.integer

  res$batch.id <- which(stringi::stri_detect_fixed(str = tsv2bam.log$INPUT,
                                                               pattern = "Batch ID: "))
  res$batch.id <- readr::read_lines(file = log.file, skip = res$batch.id - 1, n_max = 1) %>%
    stringi::stri_extract_all_regex(str = ., pattern = "[1-9]+", simplify = TRUE) %>%
    as.integer


  if (verbose) message("Number of individuals: ", res$n.ind)

  sample.indices <- which(stringi::stri_detect_fixed(str = tsv2bam.log$INPUT,
                                                     pattern = "Sample"))

  if (!identical(res$n.ind, length(sample.indices))) {
    stop("Contact author, problem with the file import workflow")
  }

  data <- readr::read_tsv(file = log.file, col_names = "INPUT", col_types = "c", skip = sample.indices[1]) %>%
    dplyr::filter(INPUT != "tsv2bam is done.")

  sample.indices <- which(stringi::stri_detect_fixed(str = data$INPUT,
                                                     pattern = "Sample"))

  if (verbose) message("Summarizing loci with a one-to-one relationship with the catalog")
  res$tsv2bam.summary <- dplyr::bind_cols(
    dplyr::slice(.data = data, sample.indices) %>%
      dplyr::rename(INDIVIDUALS = INPUT) %>%
      dplyr::mutate(INDIVIDUALS = stringi::stri_replace_all_fixed(
        str = INDIVIDUALS, pattern = c("Sample '", "':"), replacement = c("", ""), vectorize_all = FALSE)),
    dplyr::slice(.data = data, -c(1,sample.indices)) %>%
      dplyr::rename(MATCH = INPUT) %>%
      tidyr::separate(
        data = .,
        col = "MATCH",
        into = c("LOCUS", "delete_1", "delete_2", "MATCH_PERCENT", "delete_3", "ALL_LOCUS"),
        sep = " ", extra = "drop") %>%
      dplyr::select(-dplyr::contains("delete")) %>%
      dplyr::mutate(
        MATCH_PERCENT = stringi::stri_replace_all_fixed(str = MATCH_PERCENT, pattern = c("(", "%"), replacement = c("", ""), vectorize_all = FALSE),
        MATCH_PERCENT = as.numeric(MATCH_PERCENT),
        ALL_LOCUS = stringi::stri_replace_all_fixed(str = ALL_LOCUS, pattern = ")", replacement = "", vectorize_all = FALSE)
      ) %>%
      dplyr::mutate_at(.tbl = ., .vars = c("LOCUS", "ALL_LOCUS"), .funs = as.integer) %>%
      dplyr::select(ALL_LOCUS, LOCUS, MATCH_PERCENT)
  )

  sample.indices <- data <- tsv2bam.log <- NULL

  axis.theme.bold <- ggplot2::element_text(size = 12, family = "Helvetica", face = "bold")
  axis.theme.normal <- ggplot2::element_text(size = 12, family = "Helvetica")

  if (verbose) message("Generating distribution plots")
  all.sample.loci.plot <- ggplot2::ggplot(data = res$tsv2bam.summary, ggplot2::aes(x = ALL_LOCUS)) +
    ggplot2::geom_histogram() +
    ggplot2::labs(x = "Total number of loci") +
    ggplot2::labs(y = "Individuals (number)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = axis.theme.bold,
      axis.title.y = axis.theme.bold,
      legend.title = axis.theme.bold,
      legend.text = axis.theme.normal,
      axis.text.x = axis.theme.normal,
      axis.text.y = axis.theme.normal)

  sample.loci.plot <- ggplot2::ggplot(data = res$tsv2bam.summary, ggplot2::aes(x = LOCUS)) +
    ggplot2::geom_histogram() +
    ggplot2::labs(x = "Number of loci matching catalog") +
    ggplot2::labs(y = "Individuals (number)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = axis.theme.bold,
      axis.title.y = axis.theme.bold,
      legend.title = axis.theme.bold,
      legend.text = axis.theme.normal,
      axis.text.x = axis.theme.normal,
      axis.text.y = axis.theme.normal)

  percent.match.plot <- ggplot2::ggplot(data = res$tsv2bam.summary, ggplot2::aes(x = MATCH_PERCENT)) +
    ggplot2::geom_histogram() +
    ggplot2::labs(x = "Percentage of loci matching catalog") +
    ggplot2::labs(y = "Individuals (number)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = axis.theme.bold,
      axis.title.y = axis.theme.bold,
      legend.title = axis.theme.bold,
      legend.text = axis.theme.normal,
      axis.text.x = axis.theme.normal,
      axis.text.y = axis.theme.normal)

  res$distribution.tsv2bam.plots <- suppressMessages(
    cowplot::plot_grid(
    all.sample.loci.plot, sample.loci.plot, percent.match.plot,
    labels = c("A","B", "C"), align = "h", nrow = 1))
  # distribution.tsv2bam.plots

  ggplot2::ggsave(
    filename = stringi::stri_join(tsv2bam.output, "/distribution.tsv2bam.plots.pdf"),
    plot = res$distribution.tsv2bam.plots,
    width = 30, height = 10,
    dpi = 600, units = "cm", useDingbats = FALSE)

  message("Distribution plots written in the folder")
  message("Mean number of loci: ", round(mean(res$tsv2bam.summary$ALL_LOCUS)))
  message("Mean number of loci with a one-to-one relationship with the catalog: ",
          round(mean(res$tsv2bam.summary$LOCUS)), " (", round(mean(res$tsv2bam.summary$MATCH_PERCENT)), " %)")

  timing <- proc.time() - timing
  if (verbose) message("\nComputation time: ", round(timing[[3]]), " sec")
  if (verbose) cat("###################### summary_tsv2bam completed ######################\n")
  return(res)
}#End summary_tsv2bam
