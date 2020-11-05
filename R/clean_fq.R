#' @name clean_fq
#' @title Removes the noise of an individual fastq file
#' @description This function reads the fastq file of an individual and clean it
#' by removing:
#' \itemize{
#' \item unique reads with high coverage (likely paralogs or TE)
#' \item distinct reads with low coverage
#' }
#'
#' @param fq.files (character, path). The path to the individual fastq file to check.
#' Default: \code{fq.files = "my-sample.fq.gz"}.

#' @param min.coverage.threshold (integer). Minimum coverage threshold.
#' The function will remove distinct reads with coverage <= to the threshold.
#' To turn off, \code{min.coverage.threshold = NULL or 0L}.
#' Default: \code{min.coverage.threshold = 2L}.

#' @param max.coverage.threshold (integer, character). Maximum coverage threshold.
#' The function will remove distinct reads with coverage >= than this threshold.
#' To turn off, \code{max.coverage.threshold = NULL}.
#' The default, use the starting depth where high coverage unique reads are observed.
#' Default: \code{max.coverage.threshold = "high.coverage.unique.reads"}.

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

#' @param output.dir (path) Write the cleaned fq files in a specific directory.
#' Default: \code{output.dir = NULL}, uses the working directory.

#' @param parallel.core (integer) Enable parallel execution with the number of threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.


#' @details
#'
#' coming soon, just try it in the meantime...
#'

#' @rdname clean_fq
#' @export

#' @return The function returns a cleaned fq file with the name of the sample and
#' \code{-C} appended to the filename.

#' @examples
#' \dontrun{
#' require(vroom)
#'
#' # for one sample
#' clean.id <- stackr::clean_fq(
#'   fq.files = "my-sample.fq.gz",
#'   min.coverage.threshold = 7L,
#'   max.coverage.threshold = "high.coverage.unique.reads"
#'   )
#'
#' # for multiple samples in parallel
#' # require(progressr)
#'
#'  progressr::with_progress({
#'    clean <- stackr::clean_fq(
#'      fq.files = 04_process_radtags,
#'       min.coverage.threshold = 2L,
#'       max.coverage.threshold = "high.coverage.unique.reads",
#'       write.blacklist = TRUE,
#'       write.blacklist.fasta = TRUE,
#'       compress = FALSE,
#'       output.dir = "04_process_radtags/cleaned_fq"
#'  )
#'  })
#' }

clean_fq <- function(
  fq.files,
  min.coverage.threshold = 2L,
  max.coverage.threshold = "high.coverage.unique.reads",
  remove.unique.reads = TRUE,
  write.blacklist = TRUE,
  write.blacklist.fasta = TRUE,
  compress = FALSE,
  output.dir = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  opt.change <- getOption("width")
  options(width = 70)
  cat("#######################################################################\n")
  cat("######################### stackr::clean_fq ############################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Required package -----------------------------------------------------------
  if (!"vroom" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install vroom for this option:\n
                 install.packages("vroom")')
  }
  # testing
  # fq.files = "BET-ATL-60002-1682158.fq.gz"
  # min.coverage.threshold = 7L
  # max.coverage.threshold = "high.coverage.unique.reads",
  # remove.unique.reads = TRUE
  # write.blacklist = TRUE
  # write.blacklist.fasta = TRUE
  # parallel.core = parallel::detectCores() - 1
  # compress = FALSE,

  if (is.null(output.dir)) {
    output.dir <- getwd()
  } else {
    if (!dir.exists(output.dir)) dir.create(output.dir)
  }

  if (assertthat::is.string(fq.files) && assertthat::is.dir(fq.files)) {
    fq.files <- stackr::list_sample_file(f =  fq.files, full.path = TRUE, recursive = TRUE)
  }
  n.files <- length(fq.files)

  if (n.files > 1) {
    path.fq <- dirname(fq.files)
    if (length(path.fq) == 0) {
      rlang::abort("FQ files requires full path, look at the function doc examples...")
    }
    message("Cleaning ", n.files, " samples...")
  }
  if (n.files < parallel.core) parallel.core <- n.files

  if (n.files > 1) {
    future::plan(future::multisession, workers = parallel.core)

    p <- progressr::progressor(steps = n.files)
    res <- furrr::future_map(
      .x = fq.files,
      .f = clean,
      min.coverage.threshold = min.coverage.threshold,
      max.coverage.threshold = "high.coverage.unique.reads",
      remove.unique.reads = remove.unique.reads,
      write.blacklist = write.blacklist,
      write.blacklist.fasta = write.blacklist.fasta,
      compress = compress,
      output.dir = output.dir,
      parallel.core = parallel.core,
      p = p,
      .progress = FALSE
    )
  } else {
    res <- clean(
      fq.files = fq.files,
      min.coverage.threshold = min.coverage.threshold,
      max.coverage.threshold = "high.coverage.unique.reads",
      remove.unique.reads = remove.unique.reads,
      write.blacklist = write.blacklist,
      write.blacklist.fasta = write.blacklist.fasta,
      compress = compress,
      output.dir = output.dir,
      parallel.core = parallel.core,
      p = NULL
    )
  }

  timing <- proc.time() - timing
  options(width = opt.change)
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}#End clean_fq


# internal function ------------------------------------------------------------
clean <- function(
  fq.files,
  min.coverage.threshold = 2L,
  max.coverage.threshold = "high.coverage.unique.reads",
  remove.unique.reads = TRUE,
  write.blacklist = TRUE,
  write.blacklist.fasta = TRUE,
  compress = FALSE,
  output.dir = NULL,
  parallel.core = parallel::detectCores() - 1,
  p = NULL
) {

  if (!is.null(p)) p()
  # Extract sample name
  clean.names <- stackr::clean_fq_filename(basename(fq.files))
  message("Sample: ", clean.names)

  # sample name and blacklist name
  if (write.blacklist) {
    bl.min.filename <- stringi::stri_join(clean.names, "_blacklisted_min_reads.fq")
    bl.max.filename <- stringi::stri_join(clean.names, "_blacklisted_max_reads.fq")
    bl.high.unique.filename <- stringi::stri_join(clean.names, "_blacklisted_high_unique_reads.fq")
    if (compress) {
      bl.min.filename <- stringi::stri_join(clean.names, "_blacklisted_min_reads.fq.gz")
      bl.max.filename <- stringi::stri_join(clean.names, "_blacklisted_max_reads.fq.gz")
      bl.high.unique.filename <- stringi::stri_join(clean.names, "_blacklisted_high_unique_reads.fq.gz")
    }
  }
  if (write.blacklist.fasta) {
    bl.min.fasta <- stringi::stri_join(clean.names, "_blacklisted_min_reads.fasta")
    bl.max.fasta <- stringi::stri_join(clean.names, "_blacklisted_max_reads.fasta")
    bl.high.unique.fasta <- stringi::stri_join(clean.names, "_blacklisted_high_unique_reads.fasta")
    if (compress) {
      bl.min.fasta <- stringi::stri_join(clean.names, "_blacklisted_min_reads.fasta.gz")
      bl.max.fasta <- stringi::stri_join(clean.names, "_blacklisted_max_reads.fasta.gz")
      bl.high.unique.fasta <- stringi::stri_join(clean.names, "_blacklisted_high_unique_reads.fasta.gz")
    }
  }
  clean.names %<>% stringi::stri_join("-C.fq")
  if (compress) clean.names <- stringi::stri_join(clean.names, ".gz")

  bl.fq <- NULL

  fq <- vroom::vroom(
    file = fq.files,
    col_names = "READS",
    col_types = "c",
    num_threads = parallel.core,
    delim = "\t",
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
      DEPTH = n(),
      .groups = "drop"
    ) %>%
    dplyr::select(-READS) %>%
    dplyr::group_by(DEPTH) %>%
    dplyr::summarise(
      INFO = list(INFO),
      NUMBER_DISTINCT_READS = n(),
      .groups = "drop"
    )

  # a couple of required stats
  max.depth <- max(fq.stats$DEPTH, na.rm = TRUE)
  min.depth <- min(fq.stats$DEPTH, na.rm = TRUE)



  # max.coverage.threshold ---------------------------------------------------
  filter.high <- FALSE

  # max.coverage.threshold = "high.coverage.unique.reads"
  # max.coverage.threshold = 20L
  # max.coverage.threshold = 20000L
  # max.coverage.threshold = NULL
  if (!is.null(max.coverage.threshold)) {
    max.coverage.fig <- dplyr::filter(fq.stats, NUMBER_DISTINCT_READS == 1L)

    # determine high coverage unique reads thresholds
    if (nrow(max.coverage.fig) == 0) {
      high.coverage.unique.reads <- max.depth
    } else {
      high.coverage.unique.reads <- min(fq.stats$DEPTH[fq.stats$NUMBER_DISTINCT_READS == 1L]) - 1
    }

    if (max.coverage.threshold == "high.coverage.unique.reads") {
      max.coverage.threshold <- high.coverage.unique.reads
      if (max.coverage.threshold == max.depth) {
        message("The data as no high coverage unique reads")
        message("The max.coverage.threshold selected will not blacklist reads in the fq file")
        remove.unique.reads <- FALSE
        filter.high <- FALSE
      } else {
        message("The max.coverage.threshold selected will remove all high coverage unique reads")
        remove.unique.reads <- FALSE
        filter.high <- TRUE
      }
    } else {
      if (max.coverage.threshold < high.coverage.unique.reads) {
        message("The max.coverage.threshold selected will remove all high coverage unique reads")
        remove.unique.reads <- FALSE
        filter.high <- TRUE
      }
      if (max.coverage.threshold > max.depth) {
        message("The max.coverage.threshold selected is higher than the maximum depth of reads observed")
        message("The max.coverage.threshold selected will not blacklist reads in the fq file")
        remove.unique.reads <- FALSE
        filter.high <- FALSE
      }
    }

    if (filter.high) {
      bl <- fq.stats %>%
        dplyr::filter(DEPTH >= max.coverage.threshold) %$%
        INFO %>%
        unlist %>%
        sort
      message("Number of reads blacklisted with max.coverage.threshold: ", length(bl))

      if (length(bl) > 0) {
        if (write.blacklist.fasta) {
          fq %>%
            dplyr::filter(INFO %in% bl) %>%
            dplyr::distinct(READS) %>%
            vroom::vroom_write(
              x = .,
              path =  file.path(output.dir, bl.high.unique.fasta),
              col_names = FALSE,
              append = FALSE,
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
              path =  file.path(output.dir, bl.high.unique.filename),
              col_names = FALSE,
              append = FALSE,
              num_threads = parallel.core,
              progress = TRUE
            )
        }

        fq %<>% dplyr::filter(!INFO %in% bl)
      }
    }#End filter.high
  }


  # blacklist unique reads -----------------------------------------------------
  if (remove.unique.reads) {
    bl <- fq.stats %>%
      dplyr::filter(NUMBER_DISTINCT_READS == 1L) %$%
      INFO %>%
      unlist %>%
      sort
    message("Number of high coverage unique reads blacklisted: ", length(bl))

    if (length(bl) > 0) {
      if (write.blacklist.fasta) {
        fq %>%
          dplyr::filter(INFO %in% bl) %>%
          dplyr::distinct(READS) %>%
          vroom::vroom_write(
            x = .,
            path =  file.path(output.dir, bl.high.unique.fasta),
            col_names = FALSE,
            append = FALSE,
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
            path =  file.path(output.dir, bl.high.unique.filename),
            col_names = FALSE,
            append = FALSE,
            num_threads = parallel.core,
            progress = TRUE
          )
      }

      fq %<>% dplyr::filter(!INFO %in% bl)
    }
  } # remove.unique.reads

  if (is.null(min.coverage.threshold)) min.coverage.threshold <- 0L
  if (min.coverage.threshold <= min.depth) {
    message("The min.coverage.threshold selected is lower than the minimum depth of reads observed")
    message("The min.coverage.threshold selected will not blacklist reads in the fq file")
    min.coverage.threshold <- 0L
  }
  if (min.coverage.threshold > 0L) {
    bl <- fq.stats %>%
      dplyr::filter(DEPTH <= min.coverage.threshold) %$%
      INFO %>%
      unlist %>%
      sort
    message("Number of distinct reads with low coverage blacklisted: ", length(bl))

    if (length(bl) > 0) {
      if (write.blacklist.fasta) {
        fq %>%
          dplyr::filter(INFO %in% bl) %>%
          dplyr::distinct(READS) %>%
          vroom::vroom_write(
            x = .,
            path =  file.path(output.dir, bl.min.fasta),
            col_names = FALSE,
            append = FALSE,
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
            path =  file.path(output.dir, bl.min.filename),
            col_names = FALSE,
            append = FALSE,
            num_threads = parallel.core,
            progress = TRUE
          )
      }

      fq %<>% dplyr::filter(!INFO %in% bl)
    }
  }# min.coverage.threshold


  # Write the cleaned fq file---------------------------------------------------
  after.cleaning <- nrow(dplyr::filter(fq, SEQ == 2L))
  message("Number of reads after cleaning: ", after.cleaning, " (kept ", round(after.cleaning / total.sequences, 2), ")")


  fq %<>%
    dplyr::select(READS) %>%
    vroom::vroom_write(
      x = .,
      path =  file.path(output.dir, clean.names),
      col_names = FALSE,
      num_threads = parallel.core,
      progress = TRUE
    )
  return(clean.names)
}# End clean
