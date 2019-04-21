#' @name run_tsv2bam
#' @title Run STACKS tsv2bam and merges BAM files
#' @description Runs \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{tsv2bam}
#' module and additionnally, this function will also generate a summary of
#' stacks tsv2bam and will merge in parallel BAM sample files into a unique
#' BAM catalog file using SAMtools or Sambamba.
#' \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{tsv2bam}
#' converts the data (single-end or paired-end) from being organized by sample
#' into being organized by locus. This allows downstream improvements
#' (e.g. Bayesian SNP calling).

#' @param P (path, character) Path to the directory containing STACKS files.
#' Default: \code{P = "06_ustacks_cstacks_sstacks"}.
#' Inside the folder, you should have:
#' \itemize{
#'   \item \strong{the catalog files:} starting with \code{batch_} and ending with
#'   \code{.alleles.tsv.gz, .snps.tsv.gz, .tags.tsv.gz};
#'   \item \strong{3 files for each samples:} The sample name is the prefix for
#'   the files ending with:
#' \code{.alleles.tsv.gz, .snps.tsv.gz, .tags.tsv.gz}.
#' Those files are created in the
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks},
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks} and
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cxstacks.php}{cxstacks}
#' modules.
#' }

#' @param M (character, path) Path to a population map file.
#' Note that the \code{-s} option is not used inside \strong{stackr}.
#' Default: \code{M = "06_ustacks_cstacks_sstacks/population.map.tsv2bam.tsv"}.

#' @param R (path, character) Directory where to find the paired-end reads files
#' (in fastq/fasta/bam (gz) format).

#' @param t (integer) Enable parallel execution with the number of threads.
#' Default: \code{t = parallel::detectCores() - 1}

#' @param cmd.path (character, path) Provide the FULL path to SAMtools
#' program. See details on how to install SAMtools.
#' Default: \code{cmd.path = "/usr/local/bin/samtools"}.

#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @rdname run_tsv2bam
#' @export

#' @return \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{tsv2bam}
#' returns a set of \code{.matches.bam} files.
#'
#' The function \code{run_tsv2bam} returns a list with the number of individuals, the batch ID number,
#' a summary data frame and a plot containing:
#' \enumerate{
#' \item INDIVIDUALS: the sample id
#' \item ALL_LOCUS: the total number of locus for the individual (shown in subplot A)
#' \item LOCUS: the number of locus with a one-to-one relationship (shown in subplot B)
#' with the catalog
#' \item MATCH_PERCENT: the percentage of locus with a one-to-one relationship
#' with the catalog (shown in subplot C)
#'
#' Addtionally, the function returns a catalog.bam file, generated
#' by merging all the individual BAM files in parallel.
#' }

#' @details \strong{Install SAMtools}
#' \href{link to detailed instructions on how to install SAMtools}{http://gbs-cloud-tutorial.readthedocs.io/en/latest/07_start_from_new_image.html#install-gbs-radseq-software-time-required-30min}
#'
#'
#'

#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' bam.sum <- stackr::run_tsv2bam() # that's it !
#' }

#' @seealso
#'\href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#'
#' \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{stacks Version 2.0Beta6}
#'
#' \href{http://www.htslib.org}{SAMtools}
#'
#' \href{http://lomereiter.github.io/sambamba/}{Sambamba}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @references Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N.,
#' Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing
#' Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools.
#' Bioinformatics, 25, 2078-9.
#' @references Li H A statistical framework for SNP calling,
#' mutation discovery, association mapping and population genetical parameter
#' estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93.
#' @references A. Tarasov, A. J. Vilella, E. Cuppen, I. J. Nijman, and P. Prins.
#' Sambamba: fast processing of NGS alignment formats. Bioinformatics, 2015.

run_tsv2bam <- function(
  P = "06_ustacks_cstacks_sstacks",
  M = "06_ustacks_cstacks_sstacks/population.map.tsv2bam.tsv",
  R = NULL,
  t = parallel::detectCores() - 1,
  cmd.path = "/usr/local/bin/samtools",
  h = FALSE
) {

  cat("#######################################################################\n")
  cat("######################## stackr::run_tsv2bam ##########################\n")
  cat("#######################################################################\n")

  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists("06_ustacks_cstacks_sstacks")) dir.create("06_ustacks_cstacks_sstacks")
  if (!dir.exists("09_log_files")) dir.create("09_log_files")
  if (!dir.exists("08_stacks_results")) dir.create("08_stacks_results")


  # check SAMtools or Sambamba are installed -----------------------------------
  use.samtools <- stringi::stri_detect_fixed(str = cmd.path, pattern = "samtools")
  if (!file.exists(cmd.path)) {
    if (use.samtools) {
      stop("Path to SAMtools is not valid")
    } else {
      stop("Path to Sambamba is not valid")
    }
  }


  # file data and time ---------------------------------------------------------
  file.date.time <- format(Sys.time(), "%Y%m%d@%H%M")

  # logs file ------------------------------------------------------------------
  tsv2bam.log.file <- stringi::stri_join("09_log_files/tsv2bam_", file.date.time,".log")
  message("For progress, look in the log file:\n", tsv2bam.log.file)

  # tsv2bam arguments ----------------------------------------------------------

  # Input filder path
  output.folder <- P # keep a distinct copy for other use
  P <- stringi::stri_join("-P ", shQuote(P))

  # Population map path
  M <- stringi::stri_join("-M ", M)

  # Threads
  parallel.core <- t # keep a distinct copy for other use
  t <- stringi::stri_join("-t ", t)


  # Catalog batch ID
  # if (b == "guess") {
  #   b <- ""
  # } else {
  #   b <- stringi::stri_join("-b ", b)
  # }

  # paired-end path
  if (is.null(R)) {
    R <- ""
  } else {
    R <- stringi::stri_join("-R ", R)
  }

  # Help
  if (h) {
    h <- stringi::stri_join("-h ")
  } else {
    h <- ""
  }

  # command args ---------------------------------------------------------------
  command.arguments <- paste(P, M, R, t, h)

  # run command ----------------------------------------------------------------
  system2(command = "tsv2bam", args = command.arguments, stdout = tsv2bam.log.file)

  # summarize the log file -----------------------------------------------------
  message("tsv2bam completed")
  log.file <- list.files(
    path = output.folder, pattern = "tsv2bam.log", full.names = FALSE)
  new.log.file <- stringi::stri_join(
    "09_log_files/",
    stringi::stri_replace_all_fixed(
      str = log.file,
      pattern = ".log",
      replacement = stringi::stri_join("_", file.date.time, ".log"),
      vectorize_all = FALSE))
  log.file <- list.files(
    path = output.folder, pattern = "tsv2bam.log", full.names = TRUE)
  suppressMessages(transfer <- file.rename(from = log.file, to = new.log.file))

  message("\nMoving/Renaming stacks tsv2bam log file:\n", new.log.file)

  # plots.file <- list.files(
  #   path = output.folder, pattern = "distribution.tsv2bam.plots.pdf", full.names = FALSE)
  # new.plots.file <- stringi::stri_join(
  #   "08_stacks_results/",
  #   stringi::stri_replace_all_fixed(
  #     str = plots.file,
  #     pattern = ".pdf",
  #     replacement = stringi::stri_join("_", file.date.time, ".pdf"),
  #     vectorize_all = FALSE))
  # plots.file <- list.files(
  #   path = output.folder, pattern = "distribution.tsv2bam.plots.pdf", full.names = TRUE)
  # suppressMessages(transfer <- file.rename(from = plots.file, to = new.plots.file))
  # message("\nMoving/Renaming summary plots file:\n", new.plots.file)

  # Merging BAM files ----------------------------------------------------------
  if (use.samtools) {
    message("Merging BAM files with SAMtools to generate a catalog.bam file...")
  } else {
    message("Merging BAM files with Sambamba to generate a catalog.bam file...")
  }

  merge.res <- merge_parallel(
    cmd.path = cmd.path,
    output.folder = output.folder,
    parallel.core = parallel.core)

  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("########################## tsv2bam completed ##########################\n")
  return(merge.res)
}# end run_tsv2bam

# Internal nested Function -----------------------------------------------------
#' @title split_bam_list
#' @description Split a list of bam files for parallel merging
#' @rdname split_bam_list
#' @keywords internal
#' @export
split_bam_list <- function(x, parallel.core = parallel::detectCores() - 1) {
  n.bam <- length(x)
  if (n.bam < 10) {
    bam.core <- n.bam
  } else {
    bam.core <- 10
  }
  cpu.rounds <- n.bam/bam.core/parallel.core
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.bam - 1) / n.bam) + 1))
  # bam.split <- split(x, split.vec)
  # we want a df to have the list of bam and the iteration
  bam.split <- tibble::data_frame(BAM = x, SPLIT_VEC = split.vec) %>%
    split(., .$SPLIT_VEC)
  return(bam.split)
}#End split_bam_list

#' @title merge_bam
#' @description Function that merge an iteration of BAM files
#' @rdname merge_bam
#' @keywords internal
#' @export
merge_bam <- function(x, new.folder, cmd.path, output.folder, parallel.core) {
  use.samtools <- stringi::stri_detect_fixed(str = cmd.path, pattern = "samtools")

  i <- unique(x$SPLIT_VEC)
  bam.files <- stringi::stri_join(output.folder,"/", x$BAM)
  cmd.option <- "merge"
  catalog.bam.temp <- stringi::stri_join(new.folder,"/catalog_temp_", i, ".bam")

  if (use.samtools) {
    parallel.core <- stringi::stri_join("--threads ", shQuote(parallel.core))
  } else {
    parallel.core <- stringi::stri_join("-t ", shQuote(parallel.core))
  }
  # merge.log.file <- stringi::stri_join("09_log_files/stackr.mergebam_", file.date.time, ".log")
  command.arguments <- paste(c(cmd.option, parallel.core, catalog.bam.temp, bam.files), collapse = " ")
  system2(command = cmd.path, args = command.arguments)
  res <- stringi::stri_join("Finished merging iteration: ", i)
  return(res)
}# merge_bam


#' @title merge_parallel
#' @description Merge bam in parallel. This is the complete function to merge
#' in parallel BAM files
#' @rdname merge_parallel
#' @keywords internal
#' @export
merge_parallel <- function(cmd.path, output.folder, parallel.core) {
  timing <- proc.time()
  bam.list <- list.files(path = output.folder, pattern = ".matches.bam", full.names = FALSE)
  n.bam <- length(bam.list)
  message("Number of bam files to merge: ", n.bam)

  if (dir.exists("merge_bam_temp")) {
    new.folder <- "merge_bam_temp_2"
  } else {
    new.folder <- "merge_bam_temp"
  }
  dir.create(new.folder)
  merge.res <- split_bam_list(
    x = bam.list,
    parallel.core = parallel.core) %>%
    .stackr_parallel(
      X = ., FUN = merge_bam,
      mc.cores = parallel.core,
      new.folder = new.folder,
      cmd.path = cmd.path,
      output.folder = output.folder,
      parallel.core = parallel.core)

  # merge remaining bam
  use.samtools <- stringi::stri_detect_fixed(str = cmd.path, pattern = "samtools")

  # remove sambamba index bam
  if (!use.samtools) remove.bai <- file.remove(list.files(path = "merge_bam_temp", pattern = ".bai", full.names = TRUE))
  bam.files <- list.files(path = "merge_bam_temp", pattern = ".bam", full.names = TRUE)

  more.merge <- length(bam.files) > 30
  if (more.merge) {
    message("Another round of merging bam in parallel is required...")
    if (dir.exists("merge_bam_temp")) {
      new.folder <- "merge_bam_temp_2"
    } else {
      new.folder <- "merge_bam_temp"
    }
    dir.create(new.folder)
    merge.res <- split_bam_list(
      x = list.files(path = "merge_bam_temp", pattern = ".bam"),
      parallel.core = parallel.core) %>%
      .stackr_parallel(
        X = ., FUN = merge_bam,
        mc.cores = parallel.core,
        new.folder = new.folder,
        cmd.path = cmd.path,
        output.folder = "merge_bam_temp",
        parallel.core = parallel.core)

    # remove sambamba index bam
    if (!use.samtools) remove.bai <- file.remove(list.files(path = "merge_bam_temp_2", pattern = ".bai", full.names = TRUE))
    bam.files <- list.files(path = "merge_bam_temp_2", pattern = ".bam", full.names = TRUE)
    unlink(x = "merge_bam_temp", recursive = TRUE, force = TRUE)
  }

  cmd.option <- "merge"
  catalog.bam <- stringi::stri_join(output.folder, "/catalog.bam")

  if (use.samtools) {
    parallel.core <- stringi::stri_join("--threads ", shQuote(parallel.core))
  } else {
    parallel.core <- stringi::stri_join("-t ", shQuote(parallel.core))
  }
  command.arguments <- paste(c(cmd.option, parallel.core, catalog.bam, bam.files), collapse = " ")
  system2(command = cmd.path, args = command.arguments)

  # file.remove(bam.files)
  unlink(x = "merge_bam_temp", recursive = TRUE, force = TRUE)
  unlink(x = "merge_bam_temp_2", recursive = TRUE, force = TRUE)

  res <- stringi::stri_join("Finished merging, ", n.bam, " bam files")
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  return(res)
}#End merge_parallel
