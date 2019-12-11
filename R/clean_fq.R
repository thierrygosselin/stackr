#' @name clean_fq
#' @title Removes the noise of an individual fastq file
#' @description This function reads the fastq file of an individual and clean it
#' by removing:
#' \itemize{
#' \item unique reads with high coverage (likely paralogs or TE)
#' \item distinct reads with low coverage
#' }
#'
#' @param fq.file (character, path). The path to the individual fastq file to check.
#' Default: \code{fq.file = "my-sample.fq.gz"}.

#' @param min.coverage.threshold (character, path). Minimum coverage threshold.
#' The function will remove distinct reads with coverage <= to the threshold.
#' Default: \code{min.coverage.threshold = 2L}.

#' @param remove.unique.reads (logical). Remove distinct unique reads with high
#' coverage. Likely paralogs or Transposable elements.
#' Default: \code{remove.unique.reads = TRUE}.

#' @param write.blacklist (logical). Write the blacklisted reads to a file.
#' Default: \code{write.blacklist = FALSE}.

#' @param write.blacklist.fasta (logical). Write the blacklisted reads to a
#' fasta file.
#' Default: \code{write.blacklist.fasta = FALSE}.

#' @param compress (logical) To compress the output files. If you have the disk
#' space, don't compress, it's way faster this way to write.
#' Default: \code{compress = FALSE}.

#' @param parallel.core (integer) Enable parallel execution with the number of threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @details
#'
#' coming soon, just try it in the meantime...
#'

#' @rdname clean_fq
#' @export

#' @return The function returns a cleaned fq file with the name of the sample and
#' \code{_cleaned} appended to the filename.

#' @examples
#' \dontrun{
#' require(vroom)
#' clean.id <- stackr::clean_fq(
#'   fq.file = "my-sample.fq.gz",
#'   min.coverage.threshold = 7,
#'   remove.unique.reads = TRUE,
#'   write.blacklist.fasta = TRUE
#'   )
#' }

clean_fq <- function (
  fq.file,
  min.coverage.threshold = 2L,
  remove.unique.reads = TRUE,
  write.blacklist = FALSE,
  write.blacklist.fasta = FALSE,
  compress = FALSE,
  parallel.core = parallel::detectCores() - 1
) {
  # Required package -----------------------------------------------------------
  if (!"vroom" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install vroom for this option:\n
                 install.packages("vroom")')
  }
  # testing
  # fq.file = "BET-ATL-60002-1682158.fq.gz"
  # min.coverage.threshold = 7L
  # remove.unique.reads = TRUE
  # write.blacklist = TRUE
  # write.blacklist.fasta = TRUE
  # parallel.core = parallel::detectCores() - 1

  # Extract sample name
  sample.clean <- stringi::stri_replace_all_fixed(
    str = fq.file,
    pattern = c(".fq.gz", ".fq", ".fasta", ".fastq", ".gzfasta", ".gzfastq", ".fastq.gz", ".FASTQ.gz", ".FASTQ.GZ"),
    replacement = c("", "", "", "", "", "", "", "", ""),
    vectorize_all = FALSE
  )
  message("Sample: ", sample.clean)

  # sample name and blacklist name
  if (compress) {
    if (write.blacklist) bl.filename <- stringi::stri_join(sample.clean, "_blacklisted_reads.fq.gz")
    if (write.blacklist.fasta) bl.fasta <- stringi::stri_join(sample.clean, "_blacklisted_reads.fasta.gz")
    clean.name <- stringi::stri_join(sample.clean, "_cleaned.fq.gz")
  } else {
    if (write.blacklist) bl.filename <- stringi::stri_join(sample.clean, "_blacklisted_reads.fq")
    if (write.blacklist.fasta) bl.fasta <- stringi::stri_join(sample.clean, "_blacklisted_reads.fasta")
    clean.name <- stringi::stri_join(sample.clean, "_cleaned.fq")
  }
  bl.fq <- NULL

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
    )

  fq.stats <- fq %>%
    dplyr::filter(SEQ == 2L) %>%
    dplyr::select(-SEQ)

  total.sequences <- nrow(fq.stats)
  message("Number of reads: ", total.sequences)

  fq.stats %<>%
    dplyr::group_by(READS) %>%
    dplyr::summarise(
      INFO = list(INFO),
      DEPTH = n()
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-READS) %>%
    dplyr::group_by(DEPTH) %>%
    dplyr::summarise(
      INFO = list(INFO),
      NUMBER_DISTINCT_READS = n()
    ) %>%
    dplyr::ungroup(.)

  # blacklist unique reads -----------------------------------------------------
  if (remove.unique.reads) {
    bl <- fq.stats %>%
      dplyr::filter(NUMBER_DISTINCT_READS == 1L) %>%
      dplyr::select(INFO) %$%
      INFO %>%
      unlist %>%
      sort

    message("Number of high coverage unique reads blacklisted: ", length(bl))

    if (write.blacklist.fasta) {
      fq %>%
        dplyr::filter(INFO %in% bl) %>%
        dplyr::select(READS) %>%
        vroom::vroom_write(
          x = .,
          path =  bl.fasta,
          col_names = FALSE,
          append = TRUE,
          num_threads = parallel.core,
          progress = TRUE
        )
    }

    bl <- sort(c(bl - 1L, bl, bl + 1L, bl + 2L))

    if (write.blacklist) {
      fq %>%
        dplyr::filter(INFO %in% bl) %>%
        dplyr::select(READS) %>%
        vroom::vroom_write(
          x = .,
          path =  bl.filename,
          col_names = FALSE,
          append = TRUE,
          num_threads = parallel.core,
          progress = TRUE
        )
    }

    fq %<>%
      dplyr::filter(!INFO %in% bl)
  }

  if (min.coverage.threshold > 0) {
    bl <- fq.stats %>%
      dplyr::filter(DEPTH <= min.coverage.threshold) %>%
      dplyr::select(INFO) %$%
      INFO %>%
      unlist %>%
      sort
    message("Number of distinct reads with low coverage blacklisted: ", length(bl))

    if (write.blacklist.fasta) {
      fq %>%
        dplyr::filter(INFO %in% bl) %>%
        dplyr::select(READS) %>%
        vroom::vroom_write(
          x = .,
          path =  bl.fasta,
          col_names = FALSE,
          append = TRUE,
          num_threads = parallel.core,
          progress = TRUE
        )
    }


    bl <- sort(c(bl - 1L, bl, bl + 1L, bl + 2L))

    if (write.blacklist) {
      fq %>%
        dplyr::filter(INFO %in% bl) %>%
        dplyr::select(READS) %>%
        vroom::vroom_write(
          x = .,
          path =  bl.filename,
          col_names = FALSE,
          append = TRUE,
          num_threads = parallel.core,
          progress = TRUE
        )
    }

    fq %<>%
      dplyr::filter(!INFO %in% bl)
  }

  after.cleaning <- nrow(dplyr::filter(fq, SEQ == 2L))
  message("Number of reads after cleaning: ", after.cleaning, " (kept ", round(after.cleaning / total.sequences, 2), ")")


  fq %<>%
    dplyr::select(READS) %>%
    vroom::vroom_write(
      x = .,
      path =  clean.name,
      col_names = FALSE,
      num_threads = parallel.core,
      progress = TRUE
    )
  return(clean.name)
}# End clean_fq
