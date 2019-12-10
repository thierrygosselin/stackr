#' @name reads_length_distribution
#' @title Generate the read length distribution of a fastq file
#' @description This function reads the fastq file of an individual, lane or chip
#' and generate the read length distribution to help decide the threshold to cut the
#' reads to a specific length.
#'
#' @param fq.file (character, path). The path to the fastq file
#' (individal, lane or chip).
#' Default: \code{fq.file = "my-sample.fq.gz"}.

#' @param with.future (logical) When \code{TRUE} will use future package to run
#' the code in parallel. Set \code{parallel.core} to the number of physical, not
#' logical, cores. See example below.
#' Default: \code{with.future = FALSE}.

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
#' require(ShortRead)
#' reads.length.info <- stackr::reads_length_distribution(
#'   fq.file = "my-sample.fq.gz")
#'
#'  # with future package to get faster results:
#'  require(future)
#'  require(listenv)
#'  reads.length.info <- stackr::reads_length_distribution(
#'   fq.file = "my-sample.fq.gz",
#'   with.future = TRUE
#'   )
#' }


reads_length_distribution <- function(
  fq.file,
  parallel.core = parallel::detectCores() - 1,
  with.future = FALSE
) {
  timing <- proc.time()

  # Required package -----------------------------------------------------------
  # vroom turns out to be slower to do this kind of stuff...

  # if (!"vroom" %in% utils::installed.packages()[,"Package"]) {
  #   rlang::abort('Please install vroom for this option:\n
  #                install.packages("vroom")')
  # }
  if (with.future) {
    if (!"listenv" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install listenv for this option:\n
                 install.packages("listenv")')
    }
    if (!"future" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install future for this option:\n
                 install.packages("future")')
    }
  }

  # Extract sample name
  sample.clean <- stringi::stri_replace_all_fixed(
    str = fq.file,
    pattern = c(".fq.gz", ".gz",".fq", ".fasta", ".fastq", ".gzfasta", ".gzfastq", ".fastq.gz", ".FASTQ.gz", ".FASTQ.GZ"),
    replacement = c("", "", "", "", "", "", "", "", "", ""),
    vectorize_all = FALSE
  )
  message("Reading and summarizing read length for: ", sample.clean)
  # fq <- vroom::vroom(
  #   file = fq.file,
  #   col_names = "READS",
  #   col_types = "c",
  #   num_threads = parallel.core,
  #   progress = TRUE
  # ) %>%
  #   dplyr::mutate(
  #     INFO = seq.int(from = 1L, to = n()),
  #     SEQ = rep(1:4, n() / 4)
  #   ) %>%
  #   dplyr::filter(SEQ == 2L) %>%
  #   dplyr::select(READS)

  # no future
  if (with.future) {
    # with future
    future::plan(strategy = "multisession", workers = parallel.core)

    read_length <- function(fq.file, counter = 0L) {
      f <- ShortRead::FastqStreamer(fq.file, verbose = FALSE)
      res <- listenv::listenv()
      while (length(fq <- ShortRead::yield(f))) {
        new.counter <- counter + 1L
        res[[new.counter]] <- as.integer(fq@sread@ranges@width)
        counter <- new.counter
      }
      close(f)
      res <- tibble::tibble(READ_LENGTH = purrr::flatten_int(.x = as.list(res)))
      return(res)
    }#End read_length
  } else {
    read_length <- function(fq.file, counter = 0L) {
      f <- ShortRead::FastqStreamer(fq.file, verbose = FALSE)
      res <- list()
      while (length(fq <- ShortRead::yield(f))) {
        new.counter <- counter + 1L
        res[[new.counter]] <- as.integer(fq@sread@ranges@width)
        counter <- new.counter
      }
      close(f)
      res <- tibble::tibble(READ_LENGTH = purrr::flatten_int(.x = res))
      return(res)
    }#End read_length
  }


  fq <- read_length(fq.file = fq.file)
  message("Number of reads: ", nrow(fq))
  cum_length <- function(threshold, x) x <- length(x$READ_LENGTH[x$READ_LENGTH >= threshold])
  read.breaks <- read.seq <- seq(from = (max(min(fq$READ_LENGTH), 50)), to = (min(max(fq$READ_LENGTH), 200)), by = 10L)
  names(read.seq) <- read.seq
  reads.info <- tibble::tibble(
    READS_LENGTH = read.seq,
    N = purrr::map_int(.x = read.seq, .f = cum_length, x = fq)
  )
  fq <- NULL
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
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  return(list(reads.plot, reads.info))
}# End reads_length_distribution
