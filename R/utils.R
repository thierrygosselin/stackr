# BADGES
#[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/stackr?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/stackr)
#[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/stackr.svg?branch=master)](https://travis-ci.org/thierrygosselin/stackr)

# magrittr ---------------------------------------------------------------------
# Forward-pipe operator --------------------------------------------------------
#' @title Forward-pipe operator
#' @description magrittr forward-pipe operator
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

# Exposition pipe-operator ------------------------------------------------------
#' @title Exposition pipe-operator
#' @description magrittr Exposition pipe-operator
#' @name %$%
#' @rdname Exposition_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %$%
#' @usage lhs \%$\% rhs
NULL

# Compound assignment pipe operator --------------------------------------------
#' @title compound assignment pipe operator
#' @description magrittr compound assignment pipe operator
#' @name %<>%
#' @rdname compound_assignment_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL

# dplyr n ----------------------------------------------------------------------
# The number of observations in the current group.
#' @title The number of observations in the current group.
#' @description Check dplyr
#' @name n
#' @rdname n
#' @keywords internal
#' @export
#' @importFrom dplyr n
#' @usage n()
NULL

# importFrom -------------------------------------------------------------------
#' @importFrom pbmcapply pbmclapply
#' @importFrom stringdist stringdist
#' @importFrom dplyr group_by tally ungroup summarise summarise_if mutate filter distinct select bind_rows n_distinct arrange
#' @importFrom ShortRead readFastq
#' @importFrom stats IQR


# .onAttach <- function(libname, pkgname) {
#   stackr.version <- utils::packageDescription("stackr", fields = "Version")
#
#   startup.message <- stringi::stri_join("
# ******************************* IMPORTANT NOTICE *******************************\n",
# "stackr v.", stackr.version, "\n",
# "stackr was modified heavily and now focuses exclusively on running stacks
# software in R and reading particular files it produces.\n
# For filters and other functions that were in stackr,
# please see my new package called radiator.
# https://github.com/thierrygosselin/radiator
#
#
# Reproducibility note: all zenodo DOI and stackr versions were left intact.
# ********************************************************************************",
# sep = "")
# packageStartupMessage(startup.message)
# }


#' @title split_vec_row
#' @description Split input into chunk for parallel processing
#' @rdname split_vec_row
#' @keywords internal
#' @export
split_vec_row <- function(x, cpu.rounds, parallel.core = parallel::detectCores() - 1) {
  n.row <- nrow(x)
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))
  return(split.vec)
}#End split_vec_row


# .onUnload <- function(libpath) {
#   library.dynam.unload("stackr", libpath)
# }



#' @title split_tibble_rows
#' @description Split rows of tibble for parallel processing
#' @rdname split_tibble_rows
#' @keywords internal
#' @export
split_tibble_rows <- function(
  x,
  lines.cpu = 1000, #lines per CPU rounds
  parallel.core = parallel::detectCores() - 1,
  group.split = TRUE # does dplyr: group_by and group_split
) {
  n.row <- nrow(x)
  n.cores <- parallel::detectCores()
  if (parallel.core > n.cores) parallel.core <- n.cores
  if (n.row < parallel.core) parallel.core <- n.row
  if (lines.cpu > n.row) lines.cpu <- n.row
  lines.rounds <- parallel.core * lines.cpu
  x$SPLIT_VEC <- sort(rep_len(x = 1:floor(n.row / lines.rounds), length.out = n.row))
  if (group.split) {
    x %<>%
      dplyr::group_by(SPLIT_VEC) %>%
      dplyr::group_split(.tbl = ., keep = FALSE)
  }
  return(x)
}#End split_tibble_rows

#' @title list_sample_file
#' @description List sample file in folder
#' @rdname list_sample_file
#' @export
#' @keywords internal
list_sample_file <- function(f, full.path = FALSE, recursive = FALSE, paired.end = FALSE) {
  sample_file <- function(x, f) {
    sample.file <- list.files(
      path = f,
      pattern = x,
      full.names = full.path,
      recursive = recursive
    )

    # fq files with .rem.
    not.wanted <- list.files(
      path = f,
      pattern = ".rem.",
      full.names = full.path,
      recursive = recursive
    )

    sample.file <- purrr::keep(.x = sample.file, .p = !sample.file %in% not.wanted)

    if (!paired.end) {
      # reverse paired-end files
      not.wanted <- list.files(
        path = f,
        pattern = "\\.2\\.",
        full.names = full.path,
        recursive = recursive
      )

      sample.file <- purrr::keep(.x = sample.file, .p = !sample.file %in% not.wanted)
    }

    if (length(sample.file) > 0) {
      return(sample.file)
    } else {
      return(NULL)
    }
  }
  sample.list <- purrr::map(
    .x = c("fq.gz", "fq", "fasta", "fastq", "gzfasta", "gzfastq", "fastq.gz", "FASTQ.gz", "FASTQ.GZ"),
    .f = sample_file, f = f) %>%
    purrr::flatten_chr(.) %>%
    unique
  return(sample.list)
}#End list_sample_file

#' @title fq_file_type
#' @description Detect fq file type
#' @rdname fq_file_type
#' @export
#' @keywords internal
fq_file_type <- function(x) {
  fq.file.type <-  suppressWarnings(
    stringi::stri_match_all_regex(
      str = x,
      omit_no_match = TRUE,
      pattern = c( ".fq", ".fq.gz", ".fasta", ".gzfasta", ".gzfastq", ".fastq", ".fastq.gz", ".FASTQ.gz", ".FASTQ.GZ")
    ) %>%
      purrr::flatten_chr(.)
  ) %>%
    unique

  # if (identical(x = c(".fastq", ".fastq.gz", ".FASTQ.gz", ".FASTQ.GZ"), y = fq.file.type)) fq.file.type <- ".fastq.gz"
  # if (identical(x = c(".FASTQ.GZ"), y = fq.file.type)) fq.file.type <- ".fastq.gz"
  # if (identical(x = c(".FASTQ.gz"), y = fq.file.type)) fq.file.type <- ".fastq.gz"
  if (identical(x = c(".fastq", ".fastq.gz"), y = fq.file.type)) fq.file.type <- ".fastq.gz"
  if (identical(x = c(".fq", ".fq.gz"), y = fq.file.type)) fq.file.type <- ".fq.gz"
  if (identical(x = c(".fasta", ".gzfasta"), y = fq.file.type)) fq.file.type <- ".gzfasta"
  return(fq.file.type)
}#End fq_file_type


#' @title clean_fq_filename
#' @description Removes the last part of the fq filename
#' @rdname clean_fq_filename
#' @export
#' @keywords internal
clean_fq_filename <- function(x) {
  x %<>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = unique(fq_file_type(.)),
      replacement = "",
      vectorize_all = FALSE
    )
}#End clean_fq_filename


#' @title stats_stackr
#' @description Generate useful stats
#' @rdname stats_stackr
#' @export
#' @keywords internal

stats_stackr <- function(data, x, group.by = NULL, digits = NULL) {

  if (!is.null(group.by)) {
    data <- dplyr::group_by(.data = data, .data[[group.by]])
  }

  s <- dplyr::summarise(
    .data = data,
    N = n(),
    SUM = sum(.data[[x]], na.rm = TRUE),
    MEAN = mean(.data[[x]], na.rm = TRUE),
    SE = sqrt(stats::var(.data[[x]]) / length(.data[[x]])),
    SD = stats::sd(.data[[x]], na.rm = TRUE),
    MEDIAN = stats::median(.data[[x]], na.rm = TRUE),
    Q25 = stats::quantile(.data[[x]], 0.25, na.rm = TRUE),
    Q75 = stats::quantile(.data[[x]], 0.75, na.rm = TRUE),
    IQR = stats::IQR(.data[[x]], na.rm = TRUE),
    # IQR = abs(diff(stats::quantile(.data[[x]], probs = c(0.25, 0.75), na.rm = TRUE))),
    MIN = min(.data[[x]], na.rm = TRUE),
    MAX = max(.data[[x]], na.rm = TRUE),
    OUTLIERS_LOW = Q25 - (1.5 * IQR),
    OUTLIERS_HIGH = Q75 + (1.5 * IQR),
    OUTLIERS_LOW = ifelse(OUTLIERS_LOW < MIN, MIN, OUTLIERS_LOW), # don'T use dplyr::if_else here... you don't want to preserve types
    OUTLIERS_LOW_N = length(.data[[x]][.data[[x]] < OUTLIERS_LOW]),
    OUTLIERS_HIGH = ifelse(OUTLIERS_HIGH > MAX, MAX, OUTLIERS_HIGH),
    OUTLIERS_HIGH_N = length(.data[[x]][.data[[x]] > OUTLIERS_HIGH]),
    .groups = "keep"
  )

  if (!is.null(digits)) {
    s %<>% dplyr::mutate(.data = ., dplyr::across(.cols = where(is.numeric), .fns = round, digits = digits))
  }

  return(s)
}#End stats_stackr


#' @title merge_fq
#' @description Merge FQ files (single-end only)
#' @rdname merge_fq
#' @export
#' @keywords internal
merge_fq <- function(
  x,
  output.dir = NULL,
  append.filename = "-R.fq",
  parallel.core = parallel::detectCores() - 1
) {

  opt.change <- getOption("width")
  options(width = 70)
  cat("#######################################################################\n")
  cat("########################## stackr::merge_fq ###########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()



  combine_replicates_fq <- function(
    x,
    output.dir = NULL,
    append.filename = "-R.fq"
  ) {

    # output directory
    if (is.null(output.dir)) output.dir <- getwd()

    # Create a new empty fq files
    new.fq <- file.path(output.dir, paste0(unique(x$INDIVIDUALS), append.filename))
    file.create(new.fq)

    # get the fq to combined in a vector
    append_fq <- function(fq.reps, new.fq) {
      file.append(file1 = new.fq, file2 = fq.reps)
    }

    purrr::walk(.x = x$FQ_FILES,
                .f = append_fq,
                new.fq = new.fq)
    return(new.fq)
  }# End combine_replicates_fq

  if (parallel.core > 1) {
    if (length(x) < parallel.core) parallel.core <- length(x)
    new.combined.fq <- .stackr_parallel(
      X = x,
      FUN = combine_replicates_fq,
      mc.cores = parallel.core,
      output.dir = output.dir,
      append.filename = append.filename
    )
  } else {
    new.combined.fq <- purrr::map(
      .x = x,
      .f = combine_replicates_fq,
      output.dir = output.dir,
      append.filename = append.filename
    )
  }
  res <- clean_fq_filename(basename(purrr::flatten_chr(.x = new.combined.fq)))

  timing <- proc.time() - timing
  options(width = opt.change)
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}#End merge_fq



# stackr_future --------------------------------------------------------------
#' @name stackr_future
#' @title radiator parallel function
#' @description Updating radiator to use future
# @inheritParams future::plan
# @inheritParams future::availableCores
#' @inheritParams future.apply::future_apply
#' @rdname stackr_future
#' @keywords internal
stackr_future <- function(
    .x,
    .f,
    flat.future = c("int", "chr", "dfr", "dfc", "walk", "drop"),
    split.vec = FALSE,
    split.with = NULL,
    split.chunks = 4L,
    parallel.core = parallel::detectCores() - 1,
    forking = FALSE,
    ...
) {
  os <- Sys.info()[['sysname']]
  if (os == "Windows") forking <- FALSE

  opt.change <- getOption("width")
  options(width = 70)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(if (parallel.core > 1L && !forking) future::plan(strategy = "sequential"), add = TRUE)

  # argument for flattening the results
  flat.future <- match.arg(
    arg = flat.future,
    choices = c("int", "chr", "dfr", "dfc", "walk", "drop"),
    several.ok = FALSE
  )

  # splitting into chunks-------------------------------------------------------
  if (split.vec && is.null(split.with)) {
    # d: data, data length, data size
    # sv: split vector
    d <- .x
    df <- FALSE
    if (any(class(d) %in% c("tbl_df","tbl","data.frame"))) {
      d <- nrow(d)
      df <- TRUE
    }
    if (length(d) > 1L) d <- length(d)
    stopifnot(is.integer(d))
    sv <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
    # sv <- as.integer(floor((parallel.core * cpu.rounds * (seq_len(d) - 1) / d) + 1))
    stopifnot(length(sv) == d)

    # split
    if (df) {
      .x$SPLIT_VEC <- sv
      .x %<>% dplyr::ungroup(.) %>% dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% split(x = ., f = sv)
    }
  }
  if (!is.null(split.with)) {
    # check
    if (length(split.with) != 1 || !is.character(split.with)) {
      rlang::abort(message = "Contact author: problem with parallel computation")
    }
    .data <- NULL
    stopifnot(rlang::has_name(.x, split.with))
    if (split.vec) {
      sv <- dplyr::distinct(.x, .data[[split.with]])
      d <- nrow(sv)
      sv$SPLIT_VEC <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
      .x %<>%
        dplyr::left_join(sv, by = split.with) %>%
        dplyr::ungroup(.) %>%
        dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% dplyr::ungroup(.) %>% dplyr::group_split(.tbl = ., .data[[split.with]], .keep = TRUE)
    }
  }


  if (parallel.core == 1L) {
    future::plan(strategy = "sequential")
  } else {
    # parallel.core <- parallel_core_opt(parallel.core = parallel.core)
    lx <- length(.x)
    if (lx < parallel.core) {
      future::plan(strategy = "multisession", workers = lx)
    } else {
      if (!forking) future::plan(strategy = "multisession", workers = parallel.core)
    }
  }

  # Run the function in parallel and account for dots-dots-dots argument

  if (forking) {
    if (length(list(...)) == 0) {
      rad_map <- switch(
        flat.future,
        int = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_int(.)
        },
        chr = {.x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_chr(.)
        },
        dfr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_dfr(.)
        },
        dfc = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_dfc(.)
        },
        walk = {furrr::future_walk},
        drop = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core)
        }
      )
    } else {
      rad_map <- switch(
        flat.future,
        int = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_int(.)
        },
        chr = {.x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_chr(.)
        },
        dfr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_dfr(.)
        },
        dfc = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_dfc(.)
        },
        walk = {furrr::future_walk},
        drop = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core)
        }
      )
    }
  } else {
    rad_map <- switch(flat.future,
                      int = {furrr::future_map_int},
                      chr = {furrr::future_map_chr},
                      dfr = {furrr::future_map_dfr},
                      dfc = {furrr::future_map_dfc},
                      walk = {furrr::future_walk},
                      drop = {furrr::future_map}
    )
    p <- NULL
    p <- progressr::progressor(along = .x)
    opts <- furrr::furrr_options(globals = FALSE, seed = TRUE)
    if (length(list(...)) == 0) {
      .x %<>% rad_map(.x = ., .f = .f, .options = opts)
    } else {
      .x %<>% rad_map(.x = ., .f = .f, ..., .options = opts)
    }
  }
  return(.x)
}#End stackr_future
