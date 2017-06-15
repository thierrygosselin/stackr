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
#' Default: \code{o = "06_ustacks_cstacks_sstacks"}.
#' @param m Minimum depth of coverage required to create a stack.
#' Default: \code{m = 3}.
#' @param M Maximum distance (in nucleotides) allowed between stacks.
#' Default: \code{M = 2}.
#' @param N Maximum distance allowed to align secondary reads to primary stacks.
#' Default: \code{N = M + 2}.
#' @param t Input file type.
#' Supported types: fasta, fastq, gzfasta, gzfastq, fq.gz, fastq.gz.
#' Default: \code{t = "gzfastq"}.
#' @param R Retain unused reads. Default: \code{R = FALSE}.
#' @param H Disable calling haplotypes from secondary reads.
#' Default: \code{H = TRUE}.
#' @param p Enable parallel execution with num_threads threads.
#' Default: \code{p = parallel::detectCores() - 1}.
#' @param h Display this help messsage. Default: \code{h = FALSE}.


#' @param d Enable the Deleveraging algorithm, used for resolving over merged tags.
#' Default: \code{d = TRUE}.
#' @param keep_high_cov Disable the algorithm that removes highly-repetitive stacks and nearby errors.
#' Default: \code{keep_high_cov = FALSE}.
#' @param max_locus_stacks Maximum number of stacks at a single de novo locus.
#' Default: \code{max_locus_stacks = 3}.
#' @param k_len Specify k-mer size for matching between alleles and loci.
#' Default: \code{k_len = NULL}.
#' @param gapped Perform gapped alignments between stacks.
#' Default: \code{gapped = TRUE}.
#' @param max_gaps Number of gaps allowed between stacks before merging.
#' Default: \code{max_gaps = 2}.
#' @param min_aln_len Minimum length of aligned sequence in a gapped alignment.
#' Default: \code{min_aln_len = 0.8}.
#' @param model_type Either 'snp' (default), 'bounded', or 'fixed'.
#' Default: \code{model_type = "snp"}.
#' @param alpha For the SNP or Bounded SNP model,
#' Chi square significance level required to call
#' a heterozygote or homozygote, either 0.1, 0.05.
#' Default: \code{alpha = 0.05}.
#' @param bound_low For the bounded SNP model, lower bound for epsilon,
#' the error rate, between 0 and 1.0.
#' Default: \code{bound_low = 0}.
#' @param bound_high For the bounded SNP model, upper bound for epsilon,
#' the error rate, between 0 and 1.0.
#' Default: \code{bound_high = 0.2}.
#' @param bc_err_freq For the fixed model, specify the barcode error frequency, between 0 and 1.0.
#' Default: \code{bc_err_freq = NULL}.

#' @param transfer.s3 (todo) When working on Amazon CLOUD and S3.
#' Default: \code{transfer.s3 = FALSE}.
#' @param from.folder (todo) When working on Amazon CLOUD and S3.
#' Default: \code{from.folder = NULL}.
#' @param destination.folder When working on Amazon CLOUD and S3.
#' (todo) Default: \code{destination.folder = NULL}.


#' @rdname run_ustacks
#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr mutate filter distinct
#' @importFrom purrr keep walk pwalk pmap
#' @importFrom tibble data_frame

#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' returns 4 files per samples: \code{.snps.tsv.gz}, \code{.tags.tsv.gz},
#' \code{.alleles.tsv.gz},  \code{.models.tsv.gz}.

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
#' to do
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

# ustacks ----------------------------------------------------------------------

run_ustacks <- function(
  mismatch.testing = FALSE,
  sample.list = NULL,
  project.info = NULL,
  f = "04_process_radtags",
  o = "06_ustacks_cstacks_sstacks",
  m = 3,
  M = 2,
  N = M + 2,
  t = "gzfastq",
  R = FALSE,
  H = TRUE,
  p = parallel::detectCores() - 1,
  h = FALSE,
  d = TRUE,
  keep_high_cov = FALSE,
  max_locus_stacks = 3,
  k_len = NULL,
  gapped = TRUE,
  max_gaps = 2,
  min_aln_len = 0.8,
  model_type = "snp",
  alpha = 0.05,
  bound_low = 0,
  bound_high = 0.2,
  bc_err_freq = NULL,
  transfer.s3 = FALSE,
  from.folder = NULL,
  destination.folder = NULL
) {

  cat("#######################################################################\n")
  cat("######################## stackr::run_ustacks ##########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists("06_ustacks_cstacks_sstacks")) dir.create("06_ustacks_cstacks_sstacks")
  if (!dir.exists("09_log_files")) dir.create("09_log_files")

  if (mismatch.testing && length(M) < 2) {
    stop("Mismatch testing requires a range of values for M argument. e.g. M = 1:5 to test 1, 2, 3, 4, 5 mismatches")
  }

  parallel.core <- p

  run_ustacks_one_sample <- function(
    mismatch.testing = FALSE,
    sample.list = NULL,
    project.info = NULL,
    f = "04_process_radtags",
    o = "06_ustacks_cstacks_sstacks",
    m = 3,
    M = 2,
    N = M + 2,
    t = "gzfastq",
    R = FALSE,
    H = TRUE,
    p = parallel::detectCores() - 1,
    h = FALSE,
    d = TRUE,
    keep_high_cov = FALSE,
    max_locus_stacks = 3,
    k_len = NULL,
    gapped = TRUE,
    max_gaps = 2,
    min_aln_len = 0.8,
    model_type = "snp",
    alpha = 0.05,
    bound_low = 0,
    bound_high = 0.2,
    bc_err_freq = NULL,
    transfer.s3 = FALSE,
    from.folder = NULL,
    destination.folder = NULL
  ) {
    # ustacks common options ---------------------------------------------------
    if (mismatch.testing) {
      message("\nMismatch threshold testing: M = ", M, "\n")
    }

    if (is.null(sample.list)) {
      stop("No sample provided")
    } else {
      sample.list.path <- stringi::stri_join(f, "/", sample.list)
    }

    if (mismatch.testing && length(sample.list.path) > 1) {
      stop("When testing mismatch threshold, 1 sample is required")
    }

    f <- stringi::stri_join("-f ", shQuote(sample.list.path))

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
      file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
      file.date <- stringi::stri_replace_all_fixed(
        str = file.date,
        pattern = c("-", " ", ":"),
        replacement = c("", "_", ""),
        vectorize_all = FALSE
      )

      individual <- stringi::stri_join(individual, "date", file.date, "param_MNm", M, N, m, sep = "_")
      o <- stringi::stri_join(o, "/date", "_", file.date, "_", "param_MNm", "_", M, "_", N, "_", m)
    }
    if (!dir.exists(o)) {
      dir.create(o)
      message("New folder to store results:\n  ", o)
    }

    individual2 <- NULL
    if (mismatch.testing) {
      o.bk <- o
      # o.bk <- stringi::stri_join(
      #   o, "/", individual,"_mismatch_", M)
    }
    o <- stringi::stri_join("-o ", shQuote(o))

    sql.id <- as.numeric(i$SQL_ID)
    message("SQL ID: ", sql.id, "\n")

    i <- stringi::stri_join("-i ", sql.id)

    # if(missing(y)) stop("output type, either 'fastq', 'gzfastq', 'fasta',
    # or 'gzfasta' (default is to match the input file type).")
    t <- stringi::stri_join("-t ", shQuote(t))


    if (R) {
      R <- stringi::stri_join("-R ")
    } else {
      R <- ""
    }

    parallel.core <- p
    p <- stringi::stri_join("-p ", p)

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

    if (keep_high_cov) {
      keep_high_cov <- stringi::stri_join("--keep_high_cov ")
    } else {
      keep_high_cov <- ""
    }


    max_locus_stacks <- stringi::stri_join("--max_locus_stacks ", max_locus_stacks)

    if (is.null(k_len)) {
      k_len <- ""
    } else {
      k_len <- stringi::stri_join("--k_len ", k_len)
    }

    # gapped assembly options ---------------------------------------------------
    # preform gapped alignments between stacks
    # implements the Needleman–Wunsch algorithm to obtain gapped alignments
    # between putative alleles.
    if (gapped) {
      gapped <- stringi::stri_join("--gapped ")
    } else {
      gapped <- ""
    }

    # number of gaps allowed between stacks before merging (stacks default: 2)
    max_gaps <- stringi::stri_join("--max_gaps ", max_gaps)
    # minimum length of aligned sequence in a gapped alignment (default: 0.80)
    min_aln_len <- stringi::stri_join("--min_aln_len ", min_aln_len)

    # Model options --------------------------------------------------------------
    alpha <- stringi::stri_join("--alpha ", alpha)

    if (model_type == "bounded") {
      bound_low <- stringi::stri_join("--bound_low ", bound_low)
      bound_high <- stringi::stri_join("--bound_high ", bound_high)
    } else {
      bound_low <- ""
      bound_high <- ""
    }

    if (model_type == "fixed") {
      bc_err_freq <- stringi::stri_join("--bc_err_freq ", bc_err_freq)
    } else {
      bc_err_freq <- ""
    }

    model_type <- stringi::stri_join("--model_type ", model_type)


    # command args
    command.arguments <- paste(
      f, i, o, M, m, N, p, t, R, H, h,
      d, keep_high_cov, max_locus_stacks, k_len,
      gapped, max_gaps, min_aln_len,
      model_type, alpha, bound_low, bound_high, bc_err_freq
    )


    # command
    # stdout doesn't work need to use stderr instead to log the info
    system2(command = "ustacks", args = command.arguments, stderr = ustacks.sample.log.file)

    # transfer back to s3
    if (transfer.s3) {
      message("in development")
      # ustacks.files.to.s3 <- list.files(path = sample.list.path, pattern = individual, full.names = FALSE)
      # purrr::walk(.x = ustacks.files.to.s3, .f = sunnier::copy_s3, from.folder = from.folder, destination.folder = destination.folder)
    }

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

  # Samples ---------------------------
  if (is.null(sample.list)) {
    sample.list <- list.files(
      path = f,
      pattern = c("fq.gz", "fq", "fasta", "fastq", "gzfasta", "gzfastq", "fastq.gz"),
      full.names = FALSE)

    sample.list.path <- list.files(
      path = f,
      pattern = c("fq.gz", "fq", "fasta", "fastq", "gzfasta", "gzfastq", "fastq.gz"),
      full.names = TRUE)
  } else {
    sample.list.path <- stringi::stri_join(f, "/", sample.list)
  }


  # Project info file -------------------
  # sql id of the sample
  # Add SQL_ID column for ustacks if missing in project.info df
  if (is.null(project.info)) {
    potential.project.file <- list.files(path = getwd(), pattern = "project.info", ignore.case = TRUE)

    if (length(potential.project.file) > 0) {
      project.file.info <- file.info(potential.project.file) %>%
        tibble::rownames_to_column(df = ., var = "FILE") %>%
        dplyr::filter(mtime == max(mtime))
      project.info <- suppressMessages(readr::read_tsv(file = project.file.info$FILE))
      message("Using project info file found in the working directory: ", project.file.info$FILE)
    } else {
      project.info <- tibble::data_frame(
        INDIVIDUALS_REP = stringi::stri_replace_all_fixed(
          str = sample.list,
          pattern = c(".fq.gz", ".fq", ".fasta", ".fa", ".fa.gz", ".gzfastq", ".fastq", ".fastq.gz"),
          replacement = "", vectorize_all = FALSE),
        FQ_FILES = sample.list
      )
    }
    potential.project.file <- project.file.info <- NULL
  }

  if (!tibble::has_name(project.info, "SQL_ID")) {
    project.info <- project.info %>%
      dplyr::arrange(INDIVIDUALS_REP) %>%
      dplyr::mutate(SQL_ID = seq(1, n()))

    readr::write_tsv(x = project.info, path = "project.info.sqlid.tsv")
    message("Unique id info was generated: project.info.sqlid.tsv")
  }

  # Mismatch ---------------------------
  if (mismatch.testing) {
    mismatch.to.do <- length(dir(path = o, full.names = TRUE, recursive = FALSE))
    if (mismatch.to.do > 0) {
      M <- (mismatch.to.do - 1):length(M)
    }
    mismatch.to.do <- NULL

    mismatch_dir <- function(M, o) {
      o = stringi::stri_join(o, "/", stringi::stri_replace_all_fixed(str = sample.list, pattern = ".fq.gz", replacement = "", vectorize_all = FALSE),"_mismatch_", M)
      return(o)
    }

    # Map samples to ustacks --------------
    res <- purrr::pmap(
      .l = list(
        M = M,
        N = N,
        max_gaps = max_gaps,
        o = purrr::map(.x = M, mismatch_dir, o = o)
      ),
      .f = run_ustacks_one_sample,
      mismatch.testing = TRUE,
      sample.list = sample.list,
      project.info = project.info,
      transfer.s3 = transfer.s3,
      from.folder = from.folder,
      destination.folder = destination.folder,
      f = f,
      m = m,
      t = t,
      R = R,
      H = H,
      p = p,
      h = h,
      d = d,
      keep_high_cov = keep_high_cov,
      max_locus_stacks = max_locus_stacks,
      k_len = k_len,
      gapped = gapped,
      min_aln_len = min_aln_len,
      model_type = model_type,
      alpha = alpha,
      bound_low = bound_low,
      bound_high = bound_high,
      bc_err_freq = bc_err_freq
    )
  } else {
    if (is.null(sample.list)) {
      sample.list <- list.files(
        path = f,
        pattern = c("fq.gz", "fq", "fasta", "fastq", "gzfasta", "gzfastq", "fastq.gz"),
        full.names = FALSE)
    }

    # subsample to get extension
    subsample.sample.list <- sample(x = sample.list, size = ceiling(length(sample.list) * 0.10))
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "fastq.gz"))) encoding <- ".gzfastq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "fastq"))) encoding <- ".fastq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "gzfastq"))) encoding <- ".gzfastq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "fq"))) encoding <- ".fq"
    if (unique(stringi::stri_detect_fixed(str = subsample.sample.list, pattern = "fq.gz"))) encoding <- ".fq.gz"

    #Check if ustacks as already done some files
    sample.assembled <- stringi::stri_replace_all_fixed(
      str = list.files(path = o, pattern = c("alleles.tsv", "alleles.tsv.gz"), full.names = FALSE),
      pattern = ".alleles.tsv.gz", replacement = encoding, vectorize_all = FALSE)
    n.sample.assembled <- length(sample.assembled)

    sample.before <- length(sample.list)

    if (n.sample.assembled > 0) {
      sample.after <- sample.before - n.sample.assembled
      sample.list <- purrr::discard(.x = sample.list, .p = sample.list %in% sample.assembled)
      message("ustacks restarted, oups...")
      message("  Number of samples in the directory: ", sample.before)
      message("  Number of samples already assembled: ", n.sample.assembled)
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
      transfer.s3 = transfer.s3,
      from.folder = from.folder,
      destination.folder = destination.folder,
      f = f,
      o = o,
      m = m,
      M = M,
      N = N,
      t = t,
      R = R,
      H = H,
      p = p,
      h = h,
      d = d,
      keep_high_cov = keep_high_cov,
      max_locus_stacks = max_locus_stacks,
      k_len = k_len,
      gapped = gapped,
      max_gaps = max_gaps,
      min_aln_len = min_aln_len,
      model_type = model_type,
      alpha = alpha,
      bound_low = bound_low,
      bound_high = bound_high,
      bc_err_freq = bc_err_freq
    )
    res <- project.info
  }


  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}# end run_ustacks


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

  ustacks.log <- NULL
  ustacks.log <- suppressMessages(readr::read_lines(file = log.file))

  mismatch <- suppressWarnings(suppressMessages(
    readr::read_delim(
      log.file,
      delim = ":",
      skip = 2,
      n_max = 9,
      col_names = c("PARAMETER", "VALUE")) %>%
      dplyr::mutate(VALUE = as.character(VALUE))))

  n.radtags.start <- tibble::data_frame(
    PARAMETER = "Number of RAD-Tags loaded",
    VALUE = stringi::stri_extract_all_charclass(
      str = readr::read_lines(log.file, skip = 13, n_max = 1), pattern = "[0-9]"
    )[1]
  ) %>% dplyr::mutate(VALUE = as.character(VALUE))

  coverage <- tibble::data_frame(
    PARAMETER = c("Coverage_start_mean", "Coverage_start_sd", "Coverage_start_max"),
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 17, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )

  # Removing repetitive stacks
  rep.stacks <- tibble::data_frame(
    PARAMETER = "Repetitive stacks",
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 22, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )


  # post repeat removal coverage
  coverage.post.repeat <- tibble::data_frame(
    PARAMETER = c("Coverage_post_repeat_mean", "Coverage_post_repeat_sd", "Coverage_post_repeat_max"),
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 24, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )

  stacks.merge1 <- tibble::data_frame(
    PARAMETER = c("Stacks_number_merged", "Loci", "Deleveraged_loci", "Blacklisted_loci"),
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 28, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )


  # coverage after merging
  coverage.post.merging <- tibble::data_frame(
    PARAMETER = c("Coverage_post_merging_mean", "Coverage_post_merging_sd", "Coverage_post_merging_max"),
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 29, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )

  # merging secondary reads
  coverage.post.merging.sec <- tibble::data_frame(
    PARAMETER = c("Coverage_post_merging_sec_mean", "Coverage_post_merging_sec_sd", "Coverage_post_merging_sec_max"),
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 34, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )

  gap  <- tibble::data_frame(
    PARAMETER = c("Stacks_before_merged_gap", "Stacks_after_merged_gap", "Gapped_alignments"),
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 37, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )

  coverage.post.gap <- tibble::data_frame(
    PARAMETER = c("Coverage_post_gapped_alignments_mean", "Coverage_post_gapped_alignments_sd", "Coverage_post_gapped_alignments_max"),
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 38, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )

  reads.used <- tibble::data_frame(
    PARAMETER = "Reads_used",
    VALUE = stringi::stri_match_all_regex(
      str = readr::read_lines(log.file, skip = 40, n_max = 1), pattern = "[0-9]*[.]*[0-9]+"
    )[1] %>% unlist
  )

  # # Read Allele.tsv file
  # allele.summary <- suppressMessages(readr::read_tsv(
  #   file = list.files(path = ustacks.folder, pattern = ".alleles", full.names = TRUE),
  #   col_names = c("SQL_ID", "ID", "LOCUS", "HAPLOTYPE", "PERCENT", "COUNT"),
  #   comment = "#"))
  # n.locus <- dplyr::n_distinct(allele.summary$LOCUS)
  # message("Number of locus:  ", n.locus)

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
                          coverage.post.merging.sec,
                          gap,
                          coverage.post.gap,
                          reads.used,
                          polymorphism)
  return(res)

}#End read_stacks_ustacks_log
