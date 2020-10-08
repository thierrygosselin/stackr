#' @name run_ustacks
#' @title Run STACKS ustacks module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' module inside R!
#' Inside the folder \code{04_process_radtags}, you should have all the
#' individual's fastq files. Those files are created in the
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php}{process_radtags}
#' module.

#' @param mismatch.testing (logical). Default: \code{mismatch.testing = FALSE}.

#' @param sample.list The default use all the samples in \code{f},
#' usually, the stacks process_radtags output folder, here defaulted to
#' \code{f = "04_process_radtags"}, see \code{f} argument below.
#' When using \code{mismatch.testing = TRUE}, only one sample is allowed inside the folder specified in \code{f}
#' (e.g. choose one with a mean number of read, or MB, file size).
#' Default: \code{sample.list = NULL}.
#' Power outage? no problem, see details below.

#' @param project.info When using the stackr pipeline,
#' a project info file is created. This file will be modified inside this function.
#' The file is in the working directory (given the path or in the global environment).
#' If no \code{project.info} file is provided, the function first look in the
#' working directory for file(s) with "project.info" in it's name. If several files are found,
#' the latest one created is used.
#' Default: \code{project.info = NULL}.
#' Power outage? no problem, see details below.

#' @param f Input file path. Usually,
#' the stacks process_radtags output folder.
#' Default: \code{f = "04_process_radtags"}.
#' @param o Output path to write results.
#' Default: \code{o = "06_ustacks_2_gstacks"}.
#' @param m Minimum depth of coverage required to create a stack.
#' Default: \code{m = 3}.
#' @param M Maximum distance (in nucleotides) allowed between stacks.
#' Default: \code{M = 2}.
#' @param N Maximum distance allowed to align secondary reads to primary stacks.
#' Default: \code{N = M + 2}.
#' @param t Input file type.
#' Supported types: fasta, fastq, gzfasta, gzfastq, fq.gz, fastq.gz.
#' Default: \code{t = "guess"}.
#' @param R Retain unused reads. Default: \code{R = FALSE}.
#' @param H Disable calling haplotypes from secondary reads.
#' Default: \code{H = TRUE}.
#' @param parallel.core Enable parallel execution with num_threads threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.
#' @param h Display this help messsage. Default: \code{h = FALSE}.

# Stack assembly options:
#' @param d Enable the Deleveraging algorithm, used for resolving over merged tags.
#' Default: \code{d = TRUE}.
#' @param keep.high.cov Disable the algorithm that removes highly-repetitive stacks and nearby errors.
#' Default: \code{keep.high.cov = FALSE}.

#' @param high.cov.thres (double) Highly-repetitive stacks threshold,
#' in standard deviation units.
#' Default: \code{high.cov.thres = 3.0}.


#' @param max.locus.stacks Maximum number of stacks at a single de novo locus.
#' Default: \code{max.locus.stacks = 3}.
#' @param k.len Specify k-mer size for matching between alleles and loci.
#' Default: \code{k.len = NULL}.


# Gapped assembly options:
#' @param max.gaps Number of gaps allowed between stacks before merging.
#' Default: \code{max.gaps = 2}.
#' @param min.aln.len Minimum length of aligned sequence in a gapped alignment.
#' Default: \code{min.aln.len = 0.8}.
#' @param disable.gapped (logical) do not preform gapped alignments between stacks
#' (default: gapped alignements enabled).
#' Default: \code{disable.gapped = FALSE}.

# Model options:
#' @param model.type Either 'snp' (default), 'bounded', or 'fixed'.
#' Default: \code{model.type = "snp"}.
#' @param alpha For the SNP or Bounded SNP model,
#' Chi square significance level required to call
#' a heterozygote or homozygote, either 0.1, 0.05.
#' Default: \code{alpha = 0.05}.
#' @param bound.low For the bounded SNP model, lower bound for epsilon,
#' the error rate, between 0 and 1.0.
#' Default: \code{bound.low = 0}.
#' @param bound.high For the bounded SNP model, upper bound for epsilon,
#' the error rate, between 0 and 1.0.
#' Default: \code{bound.high = 0.2}.
#' @param bc.err.freq For the fixed model, specify the barcode error frequency, between 0 and 1.0.
#' Default: \code{bc.err.freq = NULL}.

# @param transfer.s3 (todo) When working on Amazon CLOUD and S3.
# Default: \code{transfer.s3 = FALSE}.
# @param from.folder (todo) When working on Amazon CLOUD and S3.
# Default: \code{from.folder = NULL}.
# @param destination.folder When working on Amazon CLOUD and S3.
# (todo) Default: \code{destination.folder = NULL}.


#' @rdname run_ustacks
#' @export

#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' returns 4 files per samples: \code{.snps.tsv.gz}, \code{.tags.tsv.gz},
#' \code{.alleles.tsv.gz},  \code{.models.tsv.gz}. In the global environment,
#' the function returns a project info file (updated if one was provided) and
#' a summary of ustacks for each samples.

#' @details \code{-i} the unique integer ID to identify the sample (SQL ID),
#' is taken from the project info file. If no project info file is provided,
#' the id is created sequentially from the sample files. This id will be written
#' in the project info file (if no file is found or given, a new file is created in the working directory).
#'
#' \strong{Power outage? No problem:}
#'
#' Restart the function as it was. After the re-start the project info file created automatically during the previous run
#' will be used. This ensure: i) that the unique SQL ids are not duplicated and
#' ii) that ustacks can start at the sample is was assembling before the outage.



#' @examples
#' \dontrun{
#' # to test different mismatches (M values from 1 to 5) on a sample with an average read number:
#' mismatch <- stackr::run_ustacks(
#' mismatch.testing = TRUE,
#' f = “mismatch_testing_folder_path”,
#' M = 1:5,
#' p = 12)
#' # the other arguments: defaults
#' # Check out the list to figure out the best threshold to use on all samples
#' # A summary table is also written in the output folder.
#'
#' # To run ustacks on all samples, using defaults:
#' ustacks.lobster <- stackr::run_ustacks()
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}



# ustacks ----------------------------------------------------------------------

run_ustacks <- function(
  mismatch.testing = FALSE,
  sample.list = NULL,
  project.info = NULL,
  f = "04_process_radtags",
  o = "06_ustacks_2_gstacks",
  m = 3,
  M = 2,
  N = M + 2,
  t = "guess",
  R = FALSE,
  H = TRUE,
  parallel.core = parallel::detectCores() - 1,
  h = FALSE,
  d = TRUE,
  keep.high.cov = FALSE,
  high.cov.thres = 3.0,
  max.locus.stacks = 3,
  k.len = NULL,
  max.gaps = 2,
  min.aln.len = 0.8,
  disable.gapped = FALSE,
  model.type = "snp",
  alpha = 0.05,
  bound.low = 0,
  bound.high = 0.2,
  bc.err.freq = NULL
) {

  cat("#######################################################################\n")
  cat("######################## stackr::run_ustacks ##########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists(o) && !mismatch.testing) dir.create(o)
  if (!dir.exists("09_log_files")) dir.create("09_log_files")

  if (mismatch.testing && length(M) < 2) {
    stop("Mismatch testing requires a range of values for M argument. e.g. M = 1:5 to test 1, 2, 3, 4, 5 mismatches")
  }

  # Samples ---------------------------
  if (is.null(sample.list)) {
    sample.list <- list_sample_file(f = f)
    sample.list.path <- list_sample_file(f = f, full.path = TRUE)
  } else {
    sample.list.path <- stringi::stri_join(
      f, "/", sample.list) %>%
      stringi::stri_replace_all_fixed(
        str = ., pattern = "//", replacement = "/", vectorize_all = TRUE)
  }

  fq.file.type <- unique(fq_file_type(sample.list))


  # Project info file -------------------
  # sql id of the sample
  # Add SQL_ID column for ustacks if missing in project.info df
  if (is.null(project.info)) {
    potential.project.file <- list.files(
      path = getwd(),
      pattern = "project.info",
      ignore.case = TRUE,
      full.names = TRUE,
      recursive = TRUE,
      include.dirs = FALSE
    )
    if (length(potential.project.file) > 0) {
      mtime <- NULL # also in global_variables.R file...
      project.file.info <- file.info(potential.project.file) %>%
        tibble::rownames_to_column(.data = ., var = "FILE") %>%
        dplyr::filter(mtime == max(mtime))
      project.info <- suppressMessages(readr::read_tsv(file = project.file.info$FILE))
      message("Using project info file found in the working directory: ", project.file.info$FILE)
    } else {
      project.info <- tibble::tibble(
        INDIVIDUALS_REP = stringi::stri_replace_all_fixed(
          str = sample.list,
          pattern = fq.file.type,
          replacement = "", vectorize_all = FALSE),
        FQ_FILES = sample.list
      )
    }
    potential.project.file <- project.file.info <- NULL
  } else {
    project.info <- suppressMessages(readr::read_tsv(file = project.info))
  }

  if (!tibble::has_name(project.info, "SQL_ID")) {
    INDIVIDUALS_REP <- NULL
    project.info <- project.info %>%
      dplyr::arrange(INDIVIDUALS_REP) %>%
      dplyr::mutate(SQL_ID = seq(1, n()))

    if (!mismatch.testing) {
      readr::write_tsv(x = project.info, path = "project.info.sqlid.tsv")
      message("Unique id info was generated: project.info.sqlid.tsv")
    }
  }

  # Mismatch -------------------------------------------------------------------
  res <- list()
  if (mismatch.testing) {
    if (f == "05_clustering_mismatches") if (!dir.exists(f)) dir.create(f)
    o <- f

    # Map samples to ustacks ---------------------------------------------------
    res$mismatches <- purrr::pmap(
      .l = list(
        M = M,
        N = N,
        max.gaps = max.gaps,
        o = purrr::map(.x = M, mismatch_dir, o = o, sample.list = sample.list, fq.file.type = fq.file.type)
      ),
      .f = run_ustacks_one_sample,
      mismatch.testing = TRUE,
      sample.list = sample.list,
      project.info = project.info,
      f = f,
      m = m,
      t = t,
      R = R,
      H = H,
      parallel.core = parallel.core,
      h = h,
      d = d,
      keep.high.cov = keep.high.cov,
      max.locus.stacks = max.locus.stacks,
      k.len = k.len,
      min.aln.len = min.aln.len,
      disable.gapped = disable.gapped,
      model.type = model.type,
      alpha = alpha,
      bound.low = bound.low,
      bound.high = bound.high,
      bc.err.freq = bc.err.freq
    )

    # summaries the info
    mismatches.summary.list <- mismatch_fig(res$mismatches)
    res$mismatches.summary <- mismatches.summary.list$mismatches.summary
    res$mismatches.plot <- mismatches.summary.list$mismatches.plot
    mismatches.summary.list <- NULL
  } else {
    if (is.null(sample.list)) {
      sample.list <- list_sample_file(f = f, full.path = FALSE)
    }

    # subsample to get extension
    subsample.sample.list <- sample(x = sample.list, size = ceiling(length(sample.list) * 0.10))
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "fastq.gz"))) encoding <- ".gzfastq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "fastq"))) encoding <- ".fastq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "gzfastq"))) encoding <- ".gzfastq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "fq"))) encoding <- ".fq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "fq.gz"))) encoding <- ".fq.gz"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "FASTQ.gz"))) encoding <- ".gzfastq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "FASTQ.GZ"))) encoding <- ".gzfastq"

    # changing the encoding for the actual fq.file.type...
    sample.assembled <- stringi::stri_replace_all_fixed(
      str = list.files(path = o, pattern = c("alleles.tsv", "alleles.tsv.gz"), full.names = FALSE),
      pattern = ".alleles.tsv.gz", replacement = fq.file.type, vectorize_all = FALSE)

    n.sample.assembled <- length(sample.assembled)
    sample.assembled <- sample.assembled[-n.sample.assembled]

    sample.before <- length(sample.list)

    sample.after <- sample.before - n.sample.assembled
    if (n.sample.assembled > 0) {
      sample.list <- purrr::discard(.x = sample.list, .p = sample.list %in% sample.assembled)
      message("ustacks restarted, oups...")
      message("  Number of samples in the directory: ", sample.before)
      message("  Number of samples already assembled, minus last sample for safe recovery: ", n.sample.assembled)
      message("  Number of samples to perform de novo assembly: ", sample.after, "\n")
    } else {
      message("Number of samples in the directory: ", sample.before)
    }

    if (sample.after == 0) stop("No samples to assemble")

    purrr::walk(
      .x = sample.list,
      run_ustacks_one_sample,
      mismatch.testing = FALSE,
      project.info = project.info,
      f = f,
      o = o,
      m = m,
      M = M,
      N = N,
      t = t,
      R = R,
      H = H,
      parallel.core = parallel.core,
      h = h,
      d = d,
      keep.high.cov = keep.high.cov,
      high.cov.thres = high.cov.thres,
      max.locus.stacks = max.locus.stacks,
      k.len = k.len,
      max.gaps = max.gaps,
      min.aln.len = min.aln.len,
      disable.gapped = disable.gapped,
      model.type = model.type,
      alpha = alpha,
      bound.low = bound.low,
      bound.high = bound.high,
      bc.err.freq = bc.err.freq
    )
    res$project.info <- project.info

    res$summary.ustacks <- stackr::summary_ustacks(
      ustacks.folder = o,
      parallel.core = parallel.core,
      verbose = FALSE)
  }#no mismatch test


  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}# end run_ustacks

# Internal function ------------------------------------------------------------
#' @title Generate mismatch directory
#' @description Generate mismatch directory
#' @rdname mismatch_dir
#' @export
#' @keywords internal
mismatch_dir <- function(M, o, sample.list, fq.file.type) {
  o = stringi::stri_join(o, "/", stringi::stri_replace_all_fixed(str = sample.list, pattern = fq.file.type, replacement = "", vectorize_all = FALSE),"_mismatch_", M)
  return(o)
}



#' @title Read stacks ustacks log
#' @description Read stacks ustacks log
#' @rdname read_stacks_ustacks_log
#' @export
#' @keywords internal
read_stacks_ustacks_log <- function(
  log.file = NULL,
  ustacks.folder = NULL,
  parallel.core = parallel::detectCores() - 1
) {

  # log.file = "09_log_files/ustacks_STU-COD-ADU-001-R-cleaned_mismatch_1.log"
  # ustacks.folder = "mismatches_tests/STU-COD-ADU-001-R-cleaned_mismatch_1"
  # parallel.core = parallel::detectCores() - 1


  ustacks.log <- suppressMessages(readr::read_lines(file = log.file))

  # Parameters
  VALUE <- NULL
  mismatch <- suppressWarnings(suppressMessages(
    readr::read_delim(
      log.file,
      delim = ":",
      skip = 2,
      n_max = 11,
      col_names = c("PARAMETER", "VALUE")) %>%
      dplyr::mutate(VALUE = as.character(VALUE))))

  # when gapped is disabled...
  remove <- which(
    stringi::stri_detect_fixed(str = mismatch$PARAMETER,
                               pattern = "Load"))

  if (length(remove) > 0L) {
    mismatch %<>% dplyr::filter(PARAMETER < remove)
  }

  # n.radtags.start ------------------------------------------------------------
  n.radtags.start <- tibble::tibble(
    PARAMETER = "Number of RAD-Tags loaded",
    VALUE = stringi::stri_extract_all_charclass(
      str = readr::read_lines(
        log.file,
        skip = which(
          stringi::stri_detect_fixed(str = ustacks.log,
                                     pattern = "Loaded")) - 1,
        n_max = 1),
      pattern = "[0-9]")[1]) %>%
    dplyr::mutate(VALUE = as.character(VALUE))

  # coverage.info --------------------------------------------------------------
  coverage.info <- which(
    stringi::stri_detect_fixed(str = ustacks.log,
                               pattern = "Initial coverage mean"))
  if (length(coverage.info) == 0) {
    coverage <- tibble::tibble(
      PARAMETER = c(
        "Coverage_start_mean",
        "Coverage_start_sd",
        "Coverage_start_max",
        "Coverage_n_reads",
        "Coverage_primary_reads_percent"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(
          log.file,
          skip = which(
            stringi::stri_detect_fixed(str = ustacks.log,
                                       pattern = "Stack coverage")) - 1,
          n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist
    )
  } else {
    coverage <- tibble::tibble(
      PARAMETER = c("Coverage_start_mean", "Coverage_start_sd", "Coverage_start_max"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(log.file, skip = coverage.info - 1, n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist
    )
  }

  # Removing repetitive stacks--------------------------------------------------
  rep.stacks.info <- which(stringi::stri_detect_fixed(str = ustacks.log,
                                                      pattern = "Removed"))
  if (length(rep.stacks.info) == 0) {
    rep.stacks <- tibble::tibble(
      PARAMETER = "Repetitive stacks blacklisted",
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(
          log.file,
          skip = which(stringi::stri_detect_regex(
            str = ustacks.log,
            pattern = "^Removing repetitive stacks")),
          n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
  } else {
    rep.stacks <- tibble::tibble(
      PARAMETER = "Repetitive stacks blacklisted",
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(
          log.file,
          skip = rep.stacks.info - 1,
          n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
  }

  # post repeat removal coverage------------------------------------------------

  coverage.info <- which(stringi::stri_detect_fixed(str = ustacks.log,
                                                    pattern = "Post-Repeat Removal"))
  if (length(coverage.info) == 0) {
    coverage.post.repeat <- tibble::tibble(
      PARAMETER = c("Coverage_post_repeat_mean", "Coverage_post_repeat_sd", "Coverage_post_repeat_max", "Coverage_post_repeat_n_reads", "Coverage_post_repeat_primary_reads_percent"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(
          log.file,
          skip = which(
            stringi::stri_detect_fixed(
              str = ustacks.log,
              pattern = "Coverage after repeat removal")) - 1,
          n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist
    )
  } else {
    coverage.post.repeat <- tibble::tibble(
      PARAMETER = c("Coverage_post_repeat_mean", "Coverage_post_repeat_sd", "Coverage_post_repeat_max"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(log.file, skip = coverage.info - 1, n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist
    )
  }

  # merged stacks --------------------------------------------------------------
  merge.stacks.info <- which(stringi::stri_detect_fixed(str = ustacks.log,
                                                        pattern = "Merging stacks"))
  if (length(merge.stacks.info) == 0) {
    stacks.merge1 <- tibble::tibble(
      PARAMETER = c("Stacks_number_merged", "Loci", "Blacklisted_loci"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(
          log.file,
          skip = which(
            stringi::stri_detect_fixed(
              str = ustacks.log,
              pattern = "Assembled"))[[1]] - 1,
          n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist)
  } else {
    stacks.merge1 <- tibble::tibble(
      PARAMETER = c("Stacks_number_merged", "Loci", "Deleveraged_loci", "Blacklisted_loci"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(log.file, skip = merge.stacks.info, n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist)
  }


  # coverage after merging------------------------------------------------------
  coverage.info <- which(stringi::stri_detect_fixed(
    str = ustacks.log,
    pattern = "After merging, coverage depth"))
  if (length(coverage.info) == 0) {
    coverage.post.merging <- tibble::tibble(
      PARAMETER = c("Coverage_post_merging_mean", "Coverage_post_merging_sd", "Coverage_post_merging_max", "Coverage_post_merging_n_reads", "Coverage_post_merging_primary_reads_percent"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(
          log.file,
          skip = which(stringi::stri_detect_fixed(
            str = ustacks.log,
            pattern = "Coverage after assembling stacks")) - 1,
          n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
  } else {
    coverage.post.merging <- tibble::tibble(
      PARAMETER = c("Coverage_post_merging_mean", "Coverage_post_merging_sd", "Coverage_post_merging_max"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(log.file, skip = coverage.info - 1, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist)
  }
  # merging secondary reads-----------------------------------------------------
  merge.sec.stacks.info <- which(stringi::stri_detect_fixed(str = ustacks.log,
                                                            pattern = "Merging secondary stacks"))
  if (length(merge.stacks.info) == 0) {
    stacks.sec.merge <- tibble::tibble(
      PARAMETER = c("Secondary_reads_merged", "Secondary_reads_total", "Secondary_reads_percent", "Reads_merged_gap_align"),
      VALUE = c(NA, NA, NA, NA))
  } else {
    stacks.sec.merge <- tibble::tibble(
      PARAMETER = c("Secondary_reads_merged", "Secondary_reads_total", "Secondary_reads_percent", "Reads_merged_gap_align"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(log.file, skip = merge.sec.stacks.info, n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist)
  }

  coverage.info <- which(stringi::stri_detect_fixed(
    str = ustacks.log,
    pattern = "After remainders merged, coverage depth"))
  if (length(coverage.info) == 0) {
    coverage.post.merging.sec <- tibble::tibble(
      PARAMETER = c("Coverage_post_merging_sec_mean", "Coverage_post_merging_sec_sd", "Coverage_post_merging_sec_max", "Coverage_post_merging_sec_n_reads", "Coverage_post_merging_sec_primary_reads_percent"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(
          log.file,
          skip = which(stringi::stri_detect_fixed(
            str = ustacks.log,
            pattern = "Coverage after merging secondary stacks")) - 1,
          n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
  } else {
    coverage.post.merging.sec <- tibble::tibble(
      PARAMETER = c("Coverage_post_merging_sec_mean", "Coverage_post_merging_sec_sd", "Coverage_post_merging_sec_max"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(log.file, skip = coverage.info - 1, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
      )[1] %>% unlist)
  }

  # GAP info--------------------------------------------------------------------
  gap.info <- which(stringi::stri_detect_fixed(
    str = ustacks.log,
    pattern = "Searching for gaps between merged stacks"))

  # gapped ON
  if (length(remove) == 0) {
    if (length(coverage.info) == 0) {
      gap <- tibble::tibble(
        PARAMETER = c("Stacks_assembled_before_gap", "Stacks_assembled_after_gap"),
        VALUE = stringi::stri_match_all_regex(
          str = readr::read_lines(
            log.file,
            skip = which(
              stringi::stri_detect_fixed(
                str = ustacks.log,
                pattern = "Assembled"))[[2]] - 1,
            n_max = 1),
          pattern = "[0-9]*[.]*[0-9]+"
        )[1] %>% unlist)
    } else {
      gap  <- tibble::tibble(
        PARAMETER = c("Stacks_before_merged_gap", "Stacks_after_merged_gap", "Gapped_alignments"),
        VALUE = stringi::stri_match_all_regex(
          str = readr::read_lines(log.file,
                                  skip = gap.info + 1,
                                  n_max = 1),
          pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
    }
    # coverage.post.gap ----------------------------------------------------------
    coverage.info <- which(stringi::stri_detect_fixed(
      str = ustacks.log,
      pattern = "After gapped alignments, coverage depth"))

    if (length(coverage.info) == 0) {
      coverage.post.gap <- tibble::tibble(
        PARAMETER = c("Coverage_post_gapped_alignments_mean",
                      "Coverage_post_gapped_alignments_sd",
                      "Coverage_post_gapped_alignments_max",
                      "Coverage_post_gapped_alignments_n_reads",
                      "Coverage_post_gapped_alignments_primary_reads_percent"),
        VALUE = stringi::stri_match_all_regex(
          str = readr::read_lines(
            log.file,
            skip = which(stringi::stri_detect_fixed(
              str = ustacks.log,
              pattern = "Coverage after gapped assembly")) - 1,
            n_max = 1),
          pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
    } else {
      coverage.post.gap <- tibble::tibble(
        PARAMETER = c("Coverage_post_gapped_alignments_mean",
                      "Coverage_post_gapped_alignments_sd",
                      "Coverage_post_gapped_alignments_max"),
        VALUE = stringi::stri_match_all_regex(
          str = readr::read_lines(log.file, skip = coverage.info - 1, n_max = 1),
          pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
    }
  } else {# Gapped OFF
    gap <- coverage.post.gap <- NULL
  }




  # Final coverage -------------------------------------------------------------
  coverage.info <- which(stringi::stri_detect_fixed(
    str = ustacks.log,
    pattern = "Final coverage"))

  if (length(coverage.info) == 0) {
    coverage.final <- NULL
  } else {
    coverage.final <- tibble::tibble(
      PARAMETER = c("Coverage_final_mean",
                    "Coverage_final_sd",
                    "Coverage_final_max",
                    "Coverage_final_n_reads",
                    "Coverage_final_primary_reads_percent"),
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(log.file, skip = coverage.info - 1, n_max = 1),
        pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
      # VALUE = stringi::stri_match_all_regex(
      #   str = readr::read_lines(
      #     log.file,
      #     skip = which(stringi::stri_detect_fixed(
      #       str = ustacks.log,
      #       pattern = "Coverage after gapped assembly")) - 1,
      #     n_max = 1),
      #   pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
  }

  # reads.used.info ------------------------------------------------------------
  reads.used.info <- which(stringi::stri_detect_fixed(
    str = ustacks.log,
    pattern = "Number of utilized reads"))

  if (length(reads.used.info) == 0) {
    reads.used <- NULL
  } else {
    reads.used <- tibble::tibble(
      PARAMETER = "Reads_used",
      VALUE = stringi::stri_match_all_regex(
        str = readr::read_lines(
          log.file,
          skip = reads.used.info - 1,
          n_max = 1), pattern = "[0-9]*[.]*[0-9]+")[1] %>% unlist)
  }

  INDIVIDUALS <- PARAMETER <- NULL
  polymorphism <- stackr::summary_ustacks(
    ustacks.folder = ustacks.folder,
    parallel.core = parallel.core) %>%
    dplyr::select(-INDIVIDUALS) %>%
    tidyr::gather(data = ., key = PARAMETER, value = VALUE) %>%
    dplyr::mutate(VALUE = as.character(VALUE))

  res <- dplyr::bind_rows(mismatch,
                          n.radtags.start,
                          coverage,
                          rep.stacks,
                          coverage.post.repeat,
                          stacks.merge1,
                          coverage.post.merging,
                          stacks.sec.merge,
                          coverage.post.merging.sec,
                          gap,
                          coverage.post.gap)

  if (!is.null(coverage.final)) {
    res <- dplyr::bind_rows(res, coverage.final)
  }

  if (!is.null(reads.used)) {
    res <- dplyr::bind_rows(res, reads.used, polymorphism)
  } else {
    res <- dplyr::bind_rows(res, polymorphism)
  }

  return(res)

}#End read_stacks_ustacks_log



#' @title mismatch_fig
#' @description Summary of mismatches
#' @rdname mismatch_fig
#' @export
#' @keywords internal
mismatch_fig <- function(mismatch.run) {
  if (!dir.exists("05_clustering_mismatches")) {
    dir.create("clustering_mismatches")
    out.path <- "clustering_mismatches"
  } else {
    out.path <- "05_clustering_mismatches"
  }

  # mismatch.run <- res$mismatches#test
  res <- list() # store results
  remove_first_column <- function(x) {
    x <- x[, -1]
  }
  res$mismatches.summary <- mismatch.run[[1]][,1] %>%
    dplyr::bind_cols(purrr::map(.x = mismatch.run,.f = remove_first_column))

  # write in the working directory
  readr::write_tsv(x = res$mismatches.summary, path = file.path(out.path, "mismatches.summary.tsv"))
  message("Summary of all mismatches written in folder: mismatches.summary.tsv")


  # lets check the impact of mismatch thresholds on polymorphism discovery
  MISMATCH <- PARAMETER <- VALUE <- TOTAL <- PROP <- NULL
  mismatch.polymorphism <- res$mismatches.summary %>%
    tidyr::gather(data = ., key = MISMATCH, value = VALUE, -PARAMETER) %>%
    dplyr::filter(PARAMETER %in% c("HOMOZYGOSITY", "HETEROZYGOSITY", "BLACKLIST_ARTIFACT")) %>%
    dplyr::group_by(MISMATCH) %>%
    dplyr::mutate(
      VALUE = as.numeric(VALUE),
      TOTAL = sum(VALUE)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
      PROP = VALUE / TOTAL,
      PERCENT = PROP * 100,
      MISMATCH = as.integer(
        stringi::stri_replace_all_fixed(
          str = MISMATCH, pattern = "M_", replacement = "",
          vectorize_all = FALSE)),
      PARAMETER = factor(x = PARAMETER,
                         levels = c("HOMOZYGOSITY", "HETEROZYGOSITY", "BLACKLIST_ARTIFACT"), ordered = TRUE),
      SPECIES = rep("TA", n()))

  # Generate figure
  n.mismatch <- dplyr::n_distinct(mismatch.polymorphism$MISMATCH)

  # polymorphism.detail <- c("1 allele/cluster\n(homozygote)",
  #                          "2 alleles/cluster\n(heterozygote)",
  #                          ">= 3 alleles/cluster\n(artifacts)")

  polymorphism.detail <- c("1 allele/cluster (homozygote)",
                           "2 alleles/cluster (heterozygote)",
                           ">= 3 alleles/cluster (artifacts)")


  PERCENT <- NULL
  res$mismatches.plot <- ggplot2::ggplot(
    mismatch.polymorphism,
    ggplot2::aes(y = PERCENT, x = MISMATCH)) +
    ggplot2::geom_line(ggplot2::aes(colour = mismatch.polymorphism$PARAMETER)) +
    ggplot2::geom_point(ggplot2::aes(colour = mismatch.polymorphism$PARAMETER,
                                     size = mismatch.polymorphism$TOTAL)) +
    ggplot2::scale_colour_manual(
      name = "Haplotype Cluster:",
      labels = polymorphism.detail, values = c("blue","green","darkred")) +
    ggplot2::scale_size_continuous(name = "Number of locus") +
    ggplot2::scale_x_continuous(breaks = 1:n.mismatch) +
    ggplot2::labs(x = "Maximum divergence between reads within a cluster") +
    ggplot2::labs(y = "Total clusters (%)") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16, family = "Helvetica",face = "bold"),
      legend.title = ggplot2::element_text(size = 14,family = "Helvetica",face = "bold"),
      legend.text = ggplot2::element_text(size = 10,family = "Helvetica"), #,face = "bold"), #legend.position = "bottom",
      axis.text = ggplot2::element_text(size = 12, family = "Helvetica"))

  ggplot2::ggsave(
    filename = file.path(out.path, "mismatches.plot.pdf"),
    plot = res$mismatches.plot,
    height = 15, width = 30, dpi = 600, units = "cm", useDingbats = FALSE)
  return(res)
}#End mismatch_fig

#' @title run_ustacks_one_sample
#' @description run ustacks for 1 sample
#' @rdname run_ustacks_one_sample
#' @export
#' @keywords internal
run_ustacks_one_sample <- function(
  mismatch.testing = FALSE,
  sample.list = NULL,
  project.info = NULL,
  f = "04_process_radtags",
  o = "06_ustacks_2_gstacks",
  m = 3,
  M = 2,
  N = M + 2,
  t = "guess",
  R = FALSE,
  H = TRUE,
  parallel.core = parallel::detectCores() - 1,
  h = FALSE,
  d = TRUE,
  keep.high.cov = FALSE,
  high.cov.thres = 3.0,
  max.locus.stacks = 3,
  k.len = NULL,
  max.gaps = 2,
  min.aln.len = 0.8,
  disable.gapped = FALSE,
  model.type = "snp",
  alpha = 0.05,
  bound.low = 0,
  bound.high = 0.2,
  bc.err.freq = NULL
) {
  # ustacks common options ---------------------------------------------------
  if (mismatch.testing) {
    message("\nMismatch threshold testing: M = ", M, "\n")
  }

  if (is.null(sample.list)) {
    stop("No sample provided")
  } else {
    sample.list.path <- stringi::stri_join(f, "/", sample.list) %>% stringi::stri_replace_all_fixed(str = ., pattern = "//", replacement = "/", vectorize_all = TRUE)
  }

  if (mismatch.testing && length(sample.list.path) > 1) {
    stop("When testing mismatch threshold, 1 sample in a folder is required")
  }

  f <- stringi::stri_join("-f ", shQuote(sample.list.path))
  FQ_FILES <- NULL

  if (!tibble::has_name(project.info, "FQ_FILES")) {
    project.info %<>%
    dplyr::rename(FQ_FILES = FQ_FILES_F)
  }
  i <- dplyr::filter(.data = project.info, FQ_FILES %in% sample.list)

  individual <- i$INDIVIDUALS_REP

  message("\nGenerating de novo assembly of sample: ", individual)

  # don't erase file already present...
  if (stringi::stri_detect_fixed(str = o, pattern = getwd())) {
    individual2 <- stringi::stri_join(o, "/", individual,".snps.tsv.gz")
  } else {
    individual2 <- stringi::stri_join(getwd(), "/", o, "/", individual,".snps.tsv.gz")
  }

  if (file.exists(individual2)) {
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    individual <- stringi::stri_join(individual, "date", file.date, "param_MNm", M, N, m, sep = "_")
    o <- stringi::stri_join(o, "/date", "_", file.date, "_", "param_MNm", "_", M, "_", N, "_", m)
  }
  if (!dir.exists(o)) {
    dir.create(o)
    message("New folder to store results:\n  ", o)
  }

  individual2 <- NULL
  o.bk <- o
  o <- stringi::stri_join("-o ", shQuote(o))

  sql.id <- as.numeric(i$SQL_ID)
  message("SQL ID: ", sql.id, "\n")

  i <- stringi::stri_join("-i ", sql.id)

  # if(missing(y)) stop("output type, either 'fastq', 'gzfastq', 'fasta',
  # or 'gzfasta' (default is to match the input file type).")
  if (t == "guess") {
    t <- ""
  } else {
    t <- stringi::stri_join("-t ", shQuote(t))
  }

  if (R) {
    R <- stringi::stri_join("-R ")
  } else {
    R <- ""
  }

  p <- stringi::stri_join("-p ", parallel.core)

  if (h) {
    h <- stringi::stri_join("-h ")
  } else {
    h <- ""
  }

  # Stack assembly and genotyping options-------------------------------------

  # ustacks.sample.log.file
  if (mismatch.testing) {
    ustacks.sample.log.file <- stringi::stri_join("09_log_files/ustacks_", individual, "_mismatch_", M, ".log")
  } else {
    ustacks.sample.log.file <- stringi::stri_join("09_log_files/ustacks_", individual, ".log")
  }

  m <- stringi::stri_join("-m ", m)
  # Max distance allowed to align secondary reads
  N <- stringi::stri_join("-N ", N)
  M.bk <- M
  M <- stringi::stri_join("-M ", M)

  if (H) {
    H <- stringi::stri_join("-H ")
  } else {
    H <- ""
  }

  if (d) {
    d <- stringi::stri_join("-d ")
  } else {
    d <- ""
  }

  if (keep.high.cov) {
    keep.high.cov <- stringi::stri_join("--keep-high-cov ")
  } else {
    keep.high.cov <- ""
  }

  high.cov.thres <- stringi::stri_join("--high-cov-thres ", high.cov.thres)


  max.locus.stacks <- stringi::stri_join("--max-locus-stacks ", max.locus.stacks)

  if (is.null(k.len)) {
    k.len <- ""
  } else {
    k.len <- stringi::stri_join("--k-len ", k.len)
  }

  # gapped assembly options ---------------------------------------------------
  # preform gapped alignments between stacks
  # implements the Needleman–Wunsch algorithm to obtain gapped alignments
  # between putative alleles.

  # number of gaps allowed between stacks before merging (stacks default: 2)
  max.gaps <- stringi::stri_join("--max-gaps ", max.gaps)
  # minimum length of aligned sequence in a gapped alignment (default: 0.80)
  min.aln.len <- stringi::stri_join("--min-aln-len ", min.aln.len)

  if (disable.gapped) {
    disable.gapped <- "--disable-gapped "
  } else {
    disable.gapped <- ""
  }



  # Model options --------------------------------------------------------------
  alpha <- stringi::stri_join("--alpha ", alpha)

  if (model.type == "bounded") {
    bound.low <- stringi::stri_join("--bound-low ", bound.low)
    bound.high <- stringi::stri_join("--bound-high ", bound.high)
  } else {
    bound.low <- ""
    bound.high <- ""
  }

  if (model.type == "fixed") {
    bc.err.freq <- stringi::stri_join("--bc-err-freq ", bc.err.freq)
  } else {
    bc.err.freq <- ""
  }

  model.type <- stringi::stri_join("--model-type ", model.type)


  # command args
  command.arguments <- paste(
    f, i, o, M, m, N, p, t, R, H, h,
    d, keep.high.cov, high.cov.thres, max.locus.stacks, k.len,
    max.gaps, min.aln.len, disable.gapped,
    model.type, alpha, bound.low, bound.high, bc.err.freq
  )


  # command
  # stdout doesn't work need to use stderr instead to log the info
  system2(command = "ustacks", args = command.arguments, stderr = ustacks.sample.log.file)

  # transfer back to s3
  # if (transfer.s3) {
  #   message("in development")
  #   # ustacks.files.to.s3 <- list.files(path = sample.list.path, pattern = individual, full.names = FALSE)
  #   # purrr::walk(.x = ustacks.files.to.s3, .f = sunnier::copy_s3, from.folder = from.folder, destination.folder = destination.folder)
  # }

  # compile log when using mismatch testing
  if (mismatch.testing) {
    message("Reading and summarizing log file")

    if (!stringi::stri_detect_fixed(str = ustacks.sample.log.file, pattern = getwd())) {
      ustacks.sample.log.file <- stringi::stri_join(getwd(), "/", ustacks.sample.log.file)
    }

    if (!stringi::stri_detect_fixed(str = o.bk, pattern = getwd())) {
      o.bk <- stringi::stri_join(getwd(), "/", o.bk)
    }

    message("directory use for test: ", o.bk)
    res <- read_stacks_ustacks_log(
      log.file = ustacks.sample.log.file,
      ustacks.folder = o.bk,
      parallel.core = parallel.core)

    colnames(res) <- c("PARAMETER", stringi::stri_join("M_", M.bk))
  } else {
    res <- "ustacks files in output directory"
  }
  return(res)

} # end run_ustacks_one_sample
