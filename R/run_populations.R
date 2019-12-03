#' @name run_populations
#' @title Run STACKS Version >= 2.0 populations module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/populations.php}{populations}
#' module inside R!

#' @param P (path, character) Path to the directory containing all the STACKS files.
#' Default: \code{P = "06_ustacks_2_gstacks"}.

#' @param V (character) Path to an input VCF file. When this module is used to
#' filter an existing vcf file.
#' Default: \code{V = NULL}.

#' @param O (character) Path to a directory where to write the output files.
#' With default: \code{O = "07_populations"}, the function creates a folder inside
#' \code{07_populations}, with date and time appended to \code{populations_}.

#' @param M path to a population map file. The format is a tab-separated file,
#' with first column containing sample name and second column population id.
#' No heather (column name).
#' e.g. \code{M = "02_project_info/population.map.catalog.tsv"}

#' @param parallel.core (integer) enable parallel execution with num_threads threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}

#' @param batch_size (integer) The number of loci (de novo mode) or
#' chromosome (reference mode), to process in a batch.
#' Increase to speed analysis, uses more memory, decrease to save memory).
#' Default in de novo mode (loci/batch): \code{batch_size = 10000}.
#' Default in reference mode (chromosome/batch): \code{batch_size = 1}.

#' @param p (integer) Minimum number of populations a locus must be present in to process a locus.
#' Default: \code{p = 1}.
#' @param r (double) Minimum percentage of individuals in a population required to process a locus for that population.
#' Default: \code{r = 0.3}.
#' @param R (double) Minimum percentage of individuals across populations required to process a locus.
#' Default: \code{r = 0.3}.
#' @param H (logical) Apply the above filters haplotype wise
#' (unshared SNPs will be pruned to reduce haplotype-wise missing data).
#' Default: \code{H = TRUE}.

#' @param min.maf (double) Specify a minimum minor allele frequency required to process a nucleotide site at a locus (0 < min.maf < 0.5).
#' Default: \code{min.maf = NULL}. Recommendation: use my other R package RADIATOR to filter data based on MAF.
#' @param min.mac (integer) Specify a minimum minor allele count required to process a nucleotide site at a locus.
#' Default: \code{min.mac = NULL}.

#' @param max.obs.het (double) Specify a maximum observed heterozygosity required to process a nucleotide site at a locus.
#' Default: \code{max.obs.het = NULL}.  Recommendation: use my other R package RADIATOR to filter data based on heterozygosity.

#' @param write.single.snp (logical) Restrict data analysis to only the first SNP per locus (implies --no-haps).
#' Default: \code{write.single.snp = FALSE}.
#' @param write.random.snp (logical) Restrict data analysis to one random SNP per locus (implies --no-haps).
#' Default: \code{write.random.snp = FALSE}.
#' @param B (character) Path to a file containing Blacklisted markers to be excluded from the export.
#' Default: \code{B = NULL}.
#' @param W (character) Path to a file containing Whitelisted markers to include in the export.
#' Default: \code{W = NULL}.

#' @param e (character) Restriction enzyme name.
#' Default: \code{e = NULL}.
#' @param merge.sites (logical) Merge loci that were produced from the same restriction enzyme cutsite (requires reference-aligned data).
#' Default: \code{merge.sites = FALSE}.
#' @param merge.prune.lim (integer) When merging adjacent loci, if at least X% samples posses both loci prune the remaining samples out of the analysis.
#' Default: \code{merge.prune.lim = NULL}.

#' @param hwe (logical) Calculate divergence from Hardy-Weinberg equilibrium
#' using the exact test at the SNP level and Guo and Thompson MCMC algorithm at
#' the haplotype level.
#' Default: \code{hwe = FALSE}.

#' @param fstats (logical) Enable SNP and haplotype-based F statistics.
#' Default: \code{fstats = FALSE}.
#' @param fst.correction (character) Specify a correction to be applied to Fst values: 'p_value', 'bonferroni_win', or 'bonferroni_gen'.
#' Default: \code{fst.correction = NULL}.
#' @param p.value.cutoff (double) maximum p-value to keep an Fst measurement.
#' Also used as base for Bonferroni correction.
#' Default: \code{p.value.cutoff = 0.05}.

#' @param k (logical) Enable kernel-smoothed Pi, Fis, Fst, Fst', and Phi_st calculations.
#' Default: \code{k = FALSE}.
#' @param smooth.fstats (logical) Enable kernel-smoothed Fst, Fst', and Phi_st calculations.
#' Default: \code{smooth.fstats = FALSE}.
#' @param smooth.popstats (logical) Enable kernel-smoothed Pi and Fis calculations.
#' Default: \code{smooth.popstats = FALSE}.

#' @param sigma (integer) Standard deviation of the kernel smoothing weight distribution.
#' Default: \code{sigma = 150000} (150kb).
#' @param bootstrap (logical) Turn on boostrap resampling for all smoothed statistics.
#' Default: \code{bootstrap = FALSE}.
#' @param N (integer) Number of bootstrap resamplings to calculate.
#' Default: \code{N = 100}.
#' @param bootstrap.pifis (logical) Turn on boostrap resampling for smoothed SNP-based Pi and Fis calculations.
#' Default: \code{bootstrap.pifis = FALSE}.
#' @param bootstrap.fst (logical) Turn on boostrap resampling for smoothed Fst calculations based on pairwise population comparison of SNPs.
#' Default: \code{bootstrap.fst = FALSE}.
#' @param bootstrap.div (logical) Turn on boostrap resampling for smoothed haplotype diveristy and gene diversity calculations based on haplotypes.
#' Default: \code{bootstrap.div = FALSE}.
#' @param bootstrap.phist (logical) Turn on boostrap resampling for smoothed Phi_st calculations based on haplotypes.
#' Default: \code{bootstrap.phist = FALSE}.
#' @param bootstrap.wl (character) Path to a whitelist file. Only use bootstrap loci contained in this whitelist.
#' Default: \code{bootstrap.wl = NULL}.


#' @param ordered.export (logical) If data is reference aligned, exports will be ordered; only a single representative of each overlapping site.
#' Default: \code{ordered.export = FALSE}.
#' @param fasta.loci (logical) Output consensus sequences of all loci, in FASTA format.
#' Default: \code{fasta.loci = FALSE}.
#' @param fasta.samples (logical) Output the sequences of the two haplotypes of each (diploid) sample, for each locus, in FASTA format.
#' Default: \code{fasta.samples = FALSE}.
#' @param fasta.samples_raw (logical) Output all haplotypes observed in each sample, for each locus, in FASTA format.
#' Default: \code{fasta.samples_raw = FALSE}.
#' @param vcf  (logical) Output SNPs in Variant Call Format (VCF).
#' Default: \code{vcf = TRUE}.
#' @param genepop (logical) Output results in GenePop format.
#' Default: \code{genepop = FALSE}.
#' @param structure (logical) Output results in Structure format.
#' Default: \code{structure = FALSE}.
#' @param radpainter (logical) Output results results in fineRADstructure/RADpainter format.
#' Default: \code{radpainter = FALSE}.
#' @param phase (logical) Output genotypes in PHASE format.
#' Default: \code{phase = FALSE}.
#' @param fastphase (logical) Output genotypes in fastPHASE format.
#' Default: \code{fastphase = FALSE}.
# @param beagle (logical) Output genotypes in Beagle format.
# Default: \code{beagle = FALSE}.
# @param beagle.phased (logical) Output haplotypes in Beagle format.
# Default: \code{beagle.phased = FALSE}.
#' @param plink (logical) Output genotypes in PLINK format.
#' Default: \code{plink = FALSE}.
#' @param hzar (logical) Output genotypes in Hybrid Zone Analysis using R (HZAR) format.
#' Default: \code{hzar = FALSE}.
#' @param phylip (logical) Output nucleotides that are fixed-within, and variant among populations in Phylip format for phylogenetic tree construction.
#' Default: \code{phylip = FALSE}.
#' @param phylip.var (logical) Include variable sites in the phylip output encoded using IUPAC notation.
#' Default: \code{phylip.var = FALSE}.
#' @param phylip.var.all (logical) Include all sequence as well as variable sites in the phylip output encoded using IUPAC notation.
#' Default: \code{phylip.var.all = FALSE}.
#' @param treemix (logical) Output SNPs in a format useable for the TreeMix program (Pickrell and Pritchard).
#' Default: \code{treemix = FALSE}.
#' @param no.hap.exports (logical) Omit haplotype outputs.
#' Default: \code{no.hap.exports = FALSE}.



#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @param verbose turn on additional logging.
#' Default: \code{verbose = FALSE}

#' @param v print program version.
#' Default: \code{v = FALSE}
#' @param log.fst.comp (logical) Log components of Fst/Phi_st calculations to a file.
#' Default: \code{log.fst.comp = FALSE}


#' @rdname run_populations
#' @export
#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/populations.php}{populations}
#' returns different output depending on arguments selected.


#' @examples
#' \dontrun{
#' pop <- stackr::run_populations(M = "population.map.tsv")
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/populations.php}{populations}

#' \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS Version 2.0b}


#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @references Guo SW, Thompson EA (1992)
#' Performing the exact test of Hardy-Weinberg proportion for multiple alleles.
#' Biometrics, 48, 361-372.

# populations ------------------------------------------------------------------

run_populations <- function(
  P = "06_ustacks_2_gstacks",
  V = NULL,
  O = "07_populations",
  M,
  parallel.core = parallel::detectCores() - 1,
  batch_size = 10000,
  p = 1,
  r = 0.3,
  R = 0.3,
  H = TRUE,
  min.maf = NULL,
  min.mac = NULL,
  max.obs.het = NULL,
  write.single.snp = FALSE,
  write.random.snp = FALSE,
  B = NULL,
  W = NULL,
  e = NULL,
  merge.sites = FALSE,
  merge.prune.lim = NULL,
  hwe = FALSE,
  fstats = FALSE,
  fst.correction = NULL,
  p.value.cutoff = 0.05,
  k = FALSE,
  smooth.fstats = FALSE,
  smooth.popstats = FALSE,
  sigma = 150000,
  bootstrap = FALSE,
  N = 100,
  bootstrap.pifis = FALSE,
  bootstrap.fst = FALSE,
  bootstrap.div = FALSE,
  bootstrap.phist = FALSE,
  bootstrap.wl = NULL,
  ordered.export = FALSE,
  fasta.samples = FALSE,
  fasta.samples_raw = FALSE,
  fasta.loci = FALSE,
  vcf = TRUE,
  genepop = FALSE,
  structure = FALSE,
  radpainter = FALSE,
  phase = FALSE,
  fastphase = FALSE,
  #beagle = FALSE,
  #beagle.phased = FALSE,
  plink = FALSE,
  hzar = FALSE,
  phylip = FALSE,
  phylip.var = FALSE,
  phylip.var.all = FALSE,
  treemix = FALSE,
  no.hap.exports = FALSE,
  h = FALSE,
  verbose = FALSE,
  v = FALSE,
  log.fst.comp = FALSE
) {

  cat("#######################################################################\n")
  cat("####################### stackr::run_populations #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists("09_log_files")) dir.create("09_log_files")
  if (!dir.exists(O)) dir.create(O)
  if (missing(M)) stop("Population map required")

  # file data and time ---------------------------------------------------------
  file.date.time <- format(Sys.time(), "%Y%m%d@%H%M")

  # logs file ------------------------------------------------------------------
  populations.log.file <- stringi::stri_join("09_log_files/populations_", file.date.time,".log")
  message("For progress, look in the log file:\n", populations.log.file)

  # common options -------------------------------------------------------------
  input.folder <- P
  P <- stringi::stri_join("-P ", shQuote(P))

  # output.folder <- O
  output.folder <- stringi::stri_join("populations_", file.date.time)
  output.folder <- file.path(O, output.folder)
  dir.create(output.folder)
  message("\nOutput files written in:\n", output.folder)
  O <- stringi::stri_join("-O ", shQuote(output.folder))

  if (is.null(V)) {
    V <- ""
  } else {
    V <- stringi::stri_join("-V ", V)
  }

  # pop map
  M <- stringi::stri_join("-M ", M)

  # parallel.core <- t
  t <- stringi::stri_join("-t ", parallel.core)
  batch_size <- stringi::stri_join("--batch_size ", batch_size)

  # Data Filtering -------------------------------------------------------------
  p <- stringi::stri_join("-p ", p)
  r <- stringi::stri_join("-r ", r)
  R <- stringi::stri_join("-R ", R)

  if (H) {
    H <- "-H "
  } else {
    H <- ""
  }

  if (is.null(min.maf)) {
    min.maf <- ""
  } else {
    min.maf <- stringi::stri_join("--min.maf ", min.maf)
  }

  if (is.null(min.mac)) {
    min.mac <- ""
  } else {
    min.mac <- stringi::stri_join("--min.mac ", min.mac)
  }

  if (is.null(max.obs.het)) {
    max.obs.het <- ""
  } else {
    max.obs.het <- stringi::stri_join("--max.obs.het ", max.obs.het)
  }

  if (write.single.snp) {
    write.single.snp <- "--write-single-snp "
  } else {
    write.single.snp <- ""
  }

  if (write.random.snp) {
    write.random.snp <- "--write-random-snp "
  } else {
    write.random.snp <- ""
  }

  if (is.null(B)) {
    B <- ""
  } else {
    B <- stringi::stri_join("-B ", B)
  }

  if (is.null(W)) {
    W <- ""
  } else {
    W <- stringi::stri_join("-W ", W)
  }


  # Merging and Phasing --------------------------------------------------------
  if (is.null(e)) {
    e <- ""
  } else {
    e <- stringi::stri_join("-e ", e)
  }

  if (merge.sites) {
    merge.sites <- "--merge-sites "
  } else {
    merge.sites <- ""
  }

  if (is.null(merge.prune.lim)) {
    merge.prune.lim <- ""
  } else {
    merge.prune.lim <- stringi::stri_join("--merge-prune-lim ", merge.prune.lim)
  }

  # Locus stats ----------------------------------------------------------------
  if (hwe) {
    hwe <- "--hwe "
  } else {
    hwe <- ""
  }

  # Fstats ---------------------------------------------------------------------
  if (fstats) {
    fstats <- "--fstats "
    if (is.null(fst.correction)) {
      fst.correction <- ""
    } else {
      fst.correction <- stringi::stri_join("--fst-correction ", fst.correction)
    }
    p.value.cutoff <- stringi::stri_join("-p-value-cutoff ", p.value.cutoff)
  } else {
    fstats <- ""
    fst.correction <- ""
    p.value.cutoff <- ""
  }

  # Kernel-smoothing algorithm -------------------------------------------------
  if (k) {
    kernel.smoothed <- TRUE
    k <- "-k "
  } else {
    kernel.smoothed <- FALSE
    k <- ""
  }

  if (smooth.fstats) {
    smooth.fstats <- TRUE
    smooth.fstats <- "--smooth-fstats "
  } else {
    smooth.fstats <- FALSE
    smooth.fstats <- ""
  }

  if (smooth.popstats) {
    smooth.popstats <- TRUE
    smooth.popstats <- "--smooth-popstats "
  } else {
    smooth.popstats <- FALSE
    smooth.popstats <- ""
  }

  if (kernel.smoothed) {
    sigma <- stringi::stri_join("--sigma ", sigma)
    if (bootstrap) {
      bootstrap <- "--bootstrap "
    } else {
      bootstrap <- ""
    }
    N <- stringi::stri_join("-N ", N)

    if (bootstrap.pifis) {
      bootstrap.pifis <- "--bootstrap-pifis "
    } else {
      bootstrap.pifis <- ""
    }

    if (bootstrap.fst) {
      bootstrap.fst <- "--bootstrap-fst "
    } else {
      bootstrap.fst <- ""
    }

    if (bootstrap.div) {
      bootstrap.div <- "--bootstrap-div "
    } else {
      bootstrap.div <- ""
    }

    if (bootstrap.phist) {
      bootstrap.phist <- "--bootstrap-phist "
    } else {
      bootstrap.phist <- ""
    }

    if (is.null(bootstrap.wl)) {
      bootstrap.wl <- ""
    } else {
      bootstrap.wl <- stringi::stri_join("--bootstrap-wl ", bootstrap.wl)
    }
  } else {
    sigma <- ""
    bootstrap <- ""
    N <- ""
    bootstrap.pifis <- ""
    bootstrap.fst <- ""
    bootstrap.div <- ""
    bootstrap.phist <- ""
    bootstrap.wl <- ""
  }

  # File output options --------------------------------------------------------
  if (ordered.export) {
    ordered.export <- "--ordered-export "
  } else {
    ordered.export <- ""
  }

  if (fasta.samples) {
    fasta.samples <- "--fasta-samples "
  } else {
    fasta.samples <- ""
  }

  if (fasta.samples_raw) {
    fasta.samples_raw <- "--fasta-samples-raw "
  } else {
    fasta.samples_raw <- ""
  }

  if (fasta.loci) {
    fasta.loci <- "--fasta-loci "
  } else {
    fasta.loci <- ""
  }

  if (vcf) {
    vcf <- "--vcf "
  } else {
    vcf <- ""
  }

  if (genepop) {
    genepop <- "--genepop "
  } else {
    genepop <- ""
  }

  if (structure) {
    structure <- "--structure "
  } else {
    structure <- ""
  }

  if (radpainter) {
    radpainter <- "--radpainter "
  } else {
    radpainter <- ""
  }

  if (phase) {
    phase <- "--phase "
  } else {
    phase <- ""
  }


  if (fastphase) {
    fastphase <- "--fastphase "
  } else {
    fastphase <- ""
  }

  # if (beagle) {
  #   beagle <- "--beagle "
  # } else {
  #   beagle <- ""
  # }
  #
  # if (beagle.phased) {
  #   beagle.phased <- "--beagle-phased "
  # } else {
  #   beagle.phased <- ""
  # }

  if (plink) {
    plink <- "--plink "
  } else {
    plink <- ""
  }


  if (hzar) {
    hzar <- "--hzar "
  } else {
    hzar <- ""
  }


  if (phylip) {
    phylip <- "--phylip "
  } else {
    phylip <- ""
  }

  if (phylip.var) {
    phylip.var <- "--phylip-var "
  } else {
    phylip.var <- ""
  }

  if (phylip.var.all) {
    phylip.var.all <- "--phylip-var-all "
  } else {
    phylip.var.all <- ""
  }

  if (treemix) {
    treemix <- "--treemix "
  } else {
    treemix <- ""
  }

  if (no.hap.exports) {
    no.hap.exports <- "--no-hap-exports "
  } else {
    no.hap.exports <- ""
  }

  # Additional options  --------------------------------------------------------
  if (h) {
    h <- "-h "
  } else {
    h <- ""
  }

  if (verbose) {
    verbose <- "--verbose "
  } else {
    verbose <- ""
  }

  if (v) {
    v <- "-v "
  } else {
    v <- ""
  }

  if (log.fst.comp) {
    log.fst.comp <- "--log-fst-comp "
  } else {
    log.fst.comp <- ""
  }


  # command args ---------------------------------------------------------------
  command.arguments <- paste(
    P, V, O, M, t, batch_size, p, r, min.maf, min.mac, max.obs.het, write.single.snp,
    write.random.snp, B, W, e, merge.sites, merge.prune.lim, hwe, fstats, fst.correction,
    p.value.cutoff, k, smooth.fstats, smooth.popstats, sigma, bootstrap, N, bootstrap.pifis, bootstrap.fst,
    bootstrap.div, bootstrap.phist, bootstrap.wl, ordered.export,
    fasta.samples, fasta.samples_raw, fasta.loci, vcf, genepop, structure,
    radpainter, phase,
    fastphase,
    # beagle, beagle.phased,
    plink, hzar, phylip, phylip.var,
    phylip.var.all, treemix,
    no.hap.exports,
    h, verbose, v, log.fst.comp
  )


  # run command ----------------------------------------------------------------
  system2(command = "populations", args = command.arguments, stderr = populations.log.file)

  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("######################## populations completed ########################\n")
  # return(res)
}# end run_populations
