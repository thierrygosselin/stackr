#' @name run_gstacks
#' @title Run STACKS new module called gstacks
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{gstacks}
#' module inside R!

#' @param denovo (logical) de novo vs reference-based mode.
#' Default: \code{denovo = TRUE}.
#' @param P (path, character) De novo mode.
#' Path to the directory containing STACKS files.
#' Default: \code{P = "06_ustacks_2_gstacks"}.
#' Inside the folder, you should have:
#' \itemize{
#'   \item \strong{the catalog files:} starting with \code{batch_} and ending with
#'   \code{.alleles.tsv.gz, .snps.tsv.gz, .tags.tsv.gz and .bam};
#'   \item \strong{5 files for each samples:} The sample name is the prefix for
#'   the files ending with:
#' \code{.alleles.tsv.gz, .models.tsv.gz, .snps.tsv.gz, .tags.tsv.gz and .bam}.
#' Those files are created in the
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks},
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks},
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cxstacks.php}{cxstacks},
#' \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{tsv2bam}
#' modules
#' }.

#' @param M (path, character) Path to a population map giving the list of samples.

# @param b (integer) De novo mode. Database/batch ID of the input catalog
# to consider. Advice: don't modify the default.
# Default: \code{b = "guess"}.

#' @param I (character, path) Reference-based mode.
#' Input directory containing BAM files.
#' Default: \code{I = NULL}.

#' @param B (character, path) Reference-based mode. Path to input BAM files.
#' The input BAM file(s) must be sorted by coordinate.
#' With \code{B}, records must be assigned to samples using BAM "reads groups"
#' (gstacks uses the ID/identifier and SM/sample name fields). Read groups
#' must be consistent if repeated different files. With \code{I},
#' read groups are unneeded and ignored.
#' Please refer to the gstacks manual page for information. stacks also provide
#' information about how to
#' generate such a BAM file with Samtools or Sambamba, and examples.
#' Default: \code{B = NULL}.
#
#' @param O (character, path) Path to output directory. \code{stackr} by default
#' will output in the input directory.
#' Default: \code{O = NULL} (same as input).


#' @param rm.pcr.duplicates (logical) remove all but one set of read pairs of
#' the same sample that have the same insert length
#' (implies \code{rm.unpaired.reads = TRUE}.).
#' Default: \code{rm.pcr.duplicates = FALSE}.

#' @param rm.unpaired.reads (logical) Discard unpaired reads.
#' Default: \code{rm.unpaired.reads = FALSE}.

#' @param ignore.pe.reads (logical) With default the function will
#' ignore paired-end reads even if present in the input.
#' Default: \code{ignore.pe.reads = TRUE}.

#' @param unpaired (logical) Ignore read pairing (only for paired-end GBS;
#' treat READ2's as if they were READ1's).
#' Default: \code{unpaired = TRUE}.

#' @param parallel.core (integer) Enable parallel execution with the number of threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param model (character) The model to use to call variants and genotypes;
#' one of \code{"marukilow"}, \code{"marukihigh"}, or \code{"snp"}.
#' See ref for more details on algorithms.
#' Default: \code{model = "marukilow"}.
#' @param var.alpha (double) Alpha threshold for discovering SNPs.
#' Default: \code{var.alpha = 0.01}.
#' @param gt.alpha (double) Alpha threshold for calling genotypes.
#' Default: \code{gt.alpha = 0.05}.


#' @param kmer.length (integer) De novo mode.
#' kmer length for the de Bruijn graph. For expert.
#' Default: \code{kmer.length = 31}.

#' @param max.debruijn.reads (integer) Maximum number of reads to use in the
#' de Bruijn graph. For expert.
#' Default: \code{max.debruijn.reads = 1000}.

#' @param min.kmer.cov (integer) De novo mode.
#' Minimum coverage to consider a kmer. For expert.
#' Default: \code{min.kmer.cov = 2}.

#' @param write.alignments (logical) Save read alignments (heavy BAM files).
#' Default: \code{write.alignments = FALSE}.

#' @param min.mapq (double) Reference-based mode.
#' Minimum PHRED-scaled mapping quality to consider a read.
#' Default: \code{min.mapq = 10}.
#' @param max.clipped (double) Reference-based mode.
#' Maximum soft-clipping level, in fraction of read length.
#' Default: \code{max.clipped = 0.20}.
#' @param max.insert.len (integer) Reference-based mode.
#' Maximum allowed sequencing insert length
#' Default: \code{max.insert.len = 1000}.

#' @param details (logical) With default the function will write a more detailed
#' output.
#' Default: \code{details = FALSE}.

#' @param phasing.cooccurrences.thr.range (integer) range of edge coverage thresholds to
#' iterate over when building the graph of allele cooccurrences for
#' SNP phasing.
#' Default: \code{phasing.cooccurrences.thr.range = c(1,2)}.

#' @param phasing.dont.prune.hets (logical) Don't try to ignore dubious heterozygote
#' genotypes during phasing. By default, during phasing, dubious het are ignored.
#' Default: \code{phasing.dont.prune.hets = FALSE}.

#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @rdname run_gstacks
#' @export

#' @return \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{tsv2bam}
#' returns a set of \code{.matches.bam} files.
#'
#' The function \code{run_gstacks} returns a list with the number of individuals, the batch ID number,
#' a summary data frame and a plot containing:
#' \enumerate{
#' \item INDIVIDUALS: the sample id
#' \item ALL_LOCUS: the total number of locus for the individual (shown in subplot A)
#' \item LOCUS: the number of locus with a one-to-one relationship (shown in subplot B)
#' with the catalog
#' \item MATCH_PERCENT: the percentage of locus with a one-to-one relationship
#' with the catalog (shown in subplot C)
#'
#' Addtionally, the function returns a batch_X.catalog.bam file that was generated
#' by merging all the individual BAM files in parallel.
#' }



#' @examples
#' \dontrun{
#' # The simplest form of the function with De novo data:
#' bam.sum <- stackr::run_gstacks() # that's it !
#' # This will use, by default, the same population map used in run_tsv2bam
#' }

#' @seealso
#'\href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}


#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @references Maruki T, Lynch M (2017)
#' Genotype Calling from Population-Genomic Sequencing Data. G3, 7, 1393-1404.

run_gstacks <- function(
  denovo = TRUE,
  P = "06_ustacks_2_gstacks",
  M = "02_project_info/population.map.tsv2bam.tsv",
  # b = "guess",
  I = NULL,
  B = NULL,
  O = NULL,
  parallel.core = parallel::detectCores() - 1,
  model = "marukilow",
  var.alpha = 0.01,
  gt.alpha = 0.05,
  rm.pcr.duplicates = FALSE,
  rm.unpaired.reads = FALSE,
  ignore.pe.reads = TRUE,
  unpaired = TRUE,
  kmer.length = 31,
  max.debruijn.reads = 1000,
  min.kmer.cov = 2,
  write.alignments = FALSE,
  min.mapq = 10,
  max.clipped = 0.20,
  max.insert.len = 1000,
  details = FALSE,
  phasing.cooccurrences.thr.range = c(1,2),
  phasing.dont.prune.hets = FALSE,
  h = FALSE
) {

  cat("#######################################################################\n")
  cat("######################## stackr::run_gstacks ##########################\n")
  cat("#######################################################################\n")

  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists("09_log_files")) dir.create("09_log_files")

  # file data and time ---------------------------------------------------------
  file.date.time <- format(Sys.time(), "%Y%m%d@%H%M")

  # logs file ------------------------------------------------------------------
  gstacks.log.file <- stringi::stri_join("09_log_files/gstacks_", file.date.time,".log")
  message("For progress, look in the log file:\n", gstacks.log.file)

  # gstacks arguments ----------------------------------------------------------
  # De novo approach -----------------------------------------------------------
  if (denovo) {
    if (!dir.exists(P)) dir.create(P)
    # Input folder path
    if (!is.null(O)) {
      O <- stringi::stri_join("-O ", O)
    } else {
      O <- stringi::stri_join("-O ", P)
    }
    P <- stringi::stri_join("-P ", shQuote(P))

    # Advanced options -------------------------------------------------------------
    kmer.length <- stringi::stri_join("--kmer-length ", kmer.length)
    max.debruijn.reads <- stringi::stri_join("--max-debruijn-reads ", max.debruijn.reads)
    min.kmer.cov <- stringi::stri_join("--min-kmer-cov ", min.kmer.cov)

    if (write.alignments) {
      write.alignments <- "--write-alignments "
    } else {
      write.alignments <- ""
    }

    # to counteract reference-based...
    B <- ""
    I <- ""
    min.mapq <- ""
    max.clipped <- ""
    max.insert.len <- ""
    phasing.cooccurrences.thr.range <- ""
    phasing.dont.prune.hets <- ""
    details <- ""


  } else {# reference-based approach -------------------------------------------
    B.bk <- B
    B <- stringi::stri_join("-B ", B)
    I <- stringi::stri_join("-I ", I)

    # Advanced options -------------------------------------------------------
    # Reference-based mode
    min.mapq <- stringi::stri_join("--min-mapq ", min.mapq)
    max.clipped <- stringi::stri_join("--max-clipped ", max.clipped)
    max.insert.len <- stringi::stri_join("--max-insert-len ", max.insert.len)
    phasing.cooccurrences.thr.range <- stringi::stri_join("--phasing-cooccurrences-thr-range ", stringi::stri_join(phasing.cooccurrences.thr.range, collapse = ","))

    if (!phasing.dont.prune.hets) {
      phasing.dont.prune.hets <- "--phasing-dont-prune-hets "
    } else {
      phasing.dont.prune.hets <- ""
    }
    # details
    if (details) {
      details <- "--details"
    }

    # counteract de novo
    O <- ""
    P <- ""
    kmer.length <- ""
    max.debruijn.reads <- ""
    min.kmer.cov <- ""
    write.alignments <- ""

  }# End reference-based mode


  # Shared options -------------------------------------------------------------
  # Population map path
  M <- stringi::stri_join("-M ", M)

  # Threads
  t <- stringi::stri_join("-t ", parallel.core)

  # SNP Model options ----------------------------------------------------------
  model <- stringi::stri_join("--model ", model)
  var.alpha <- stringi::stri_join("--var-alpha ", var.alpha)
  gt.alpha <- stringi::stri_join("--gt-alpha ", gt.alpha)


  # paired-end -----------------------------------------------------------------
  # ignore.pe.reads = TRUE
  if (ignore.pe.reads) {
    ignore.pe.reads <- "--ignore-pe-reads"
  } else {
    ignore.pe.reads <- ""
  }

  if (!unpaired) {
    unpaired <- ""
    if (rm.unpaired.reads) {
      rm.unpaired.reads <- "--rm-unpaired-reads"
      if (rm.pcr.duplicates) {
        rm.pcr.duplicates <- "--rm-pcr-duplicates"
      } else {
        rm.pcr.duplicates <- ""
      }
    } else {
      rm.unpaired.reads <- ""
      rm.pcr.duplicates <- ""
    }
  } else {
    unpaired <- "--unpaired"
    rm.unpaired.reads <- ""
    rm.pcr.duplicates <- ""
  }

  if (denovo) unpaired <- ""


  # Help
  if (h) {
    h <- stringi::stri_join("-h ")
  } else {
    h <- ""
  }

  # command args ---------------------------------------------------------------
  command.arguments <- paste(
    P,
    I, B,
    M, O, t,
    model, var.alpha, gt.alpha,
    rm.pcr.duplicates, rm.unpaired.reads, ignore.pe.reads, unpaired,
    kmer.length, max.debruijn.reads, min.kmer.cov, write.alignments,
    min.mapq, max.clipped, max.insert.len, details,
    phasing.cooccurrences.thr.range, phasing.dont.prune.hets,
    h
  )

  # run command ----------------------------------------------------------------
  system2(command = "gstacks", args = command.arguments, stderr = gstacks.log.file)
  # stdout = gstacks.log.file

  # summarize the log file -----------------------------------------------------
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("########################## gstacks completed ##########################\n")
  res <- "gstacks finished"
  return(res)
}# end run_gstacks
