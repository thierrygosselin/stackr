#' @name reads_length_distribution
#' @title Generate the read length distribution of an individual fastq file
#' @description This function reads the fastq file of an individual and
#' generate the read length distribution to help decide the threshold to cut the
#' reads to a specific length.
#'
#' @param fq.file (character, path). The path to the individual fastq file to check.
#' Default: \code{fq.file = "my-sample.fq.gz"}.

#' @param parallel.core (integer) Enable parallel execution with the number of threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @details
#'
#' coming soon, just try it in the meantime...
#'

#' @rdname reads_length_distribution
#' @export

#' @return The function returns a plot and a tibble with potential reads length
#' thresholds and associated number of reads.

#' @examples
#' \dontrun{
#' require(vroom)
#' reads.length.info <- stackr::reads_length_distribution(
#'   fq.file = "my-sample.fq.gz")
#' }

reads_length_distribution <- function(
  fq.file,
  parallel.core = parallel::detectCores() - 1
) {
  # Required package -----------------------------------------------------------
  if (!"vroom" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install vroom for this option:\n
                 install.packages("vroom")')
  }
  # testing
  # fq.file = "BET-ATL-60002-1682158.fq.gz"
  # parallel.core = parallel::detectCores() - 1

  # Extract sample name
  sample.clean <- stringi::stri_replace_all_fixed(
    str = fq.file,
    pattern = c(".fq.gz", ".fq", ".fasta", ".fastq", ".gzfasta", ".gzfastq", ".fastq.gz", ".FASTQ.gz", ".FASTQ.GZ"),
    replacement = c("", "", "", "", "", "", "", "", ""),
    vectorize_all = FALSE
  )
  message("Sample: ", sample.clean)

  fq <- vroom::vroom(
    file = fq.file,
    col_names = "READS",
    col_types = "c",
    num_threads = parallel.core,
    progress = TRUE
  ) %>%
    dplyr::mutate(
      INFO = seq.int(from = 1L, to = n()),
      SEQ = rep(1:4, n() / 4)
    ) %>%
    dplyr::filter(SEQ == 2L) %>%
    dplyr::select(READS)

  total.sequences <- nrow(fq)
  message("Number of reads: ", total.sequences)

  fq %<>%
    dplyr::mutate(LENGTH = stringi::stri_length(READS))


  cum_length <- function(threshold, x) x <- length(x$READS[x$LENGTH >= threshold])
  read.breaks <- read.seq <- seq(from = (max(min(fq$LENGTH), 50)), to = (max(min(fq$LENGTH), 200)), by = 10L)
  names(read.seq) <- read.seq
  reads.info <- tibble::tibble(
    READS_LENGTH = read.seq,
    N = purrr::map_int(.x = read.seq, .f = cum_length, x = fq)
    )

  reads.plot <- ggplot2::ggplot(
    data = reads.info,
    ggplot2::aes(x = READS_LENGTH, y = N)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
    ggplot2::scale_x_continuous(name = "Read length maximum size", breaks = read.breaks) +
    ggplot2::scale_y_continuous(name = "Number of reads")+
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8)
    )
  filename.plot <- stringi::stri_join(sample.clean, "_reads_length_dist.png")
  ggplot2::ggsave(
    plot = reads.plot,
    filename = filename.plot,
    width = 25,
    height = 15,
    dpi = 300,
    units = "cm"
  )
  return(list(reads.plot, reads.info))
}# End reads_length_distribution
