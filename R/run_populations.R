#' @name run_populations
#' @title Run STACKS Version 2.0Beta8 populations module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/populations.php}{populations}
#' module inside R!

#' @param P (path, character) Path to the directory containing all the STACKS files.
#' Default: \code{P = "06_ustacks_cstacks_sstacks"}.

#' @param V (character) Path to an input VCF file. When this module is used to
#' filter an existing vcf file.
#' Default: \code{V = NULL}.

#' @param O (character) Path to a directory where to write the output files.
#' With default: \code{O = "07_populations"}, the function creates a folder inside
#' \code{07_populations}, with date and time appended to \code{stackr_stacks_populations_}.

#' @param M path to a population map file. The format is a tab-separated file,
#' with first column containing sample name and second column population id.
#' No heather (column name).
#' e.g. \code{M = "02_project_info/population.map.catalog.tsv"}

#' @param t (integer) enable parallel execution with num_threads threads.
#' Default: \code{t = parallel::detectCores() - 1}

#' @param b (integer) Database/batch ID of the input catalog to consider.
#' Advice: don't modify the default.
#' Default: \code{b = "guess"}.

#' @param batch_size (integer) The number of loci (de novo mode) or
#' chromosome (reference mode), to process in a batch.
#' Increase to speed analysis, uses more memory, decrease to save memory).
#' Default in de novo mode (loci/batch): \code{batch_size = 10000}.
#' Default in reference mode (chromosome/batch): \code{batch_size = 1}.

# @param s (logical) Output a file to import results into an SQL database.
# Default: \code{s = FALSE}.


#' @param p (integer) Minimum number of populations a locus must be present in to process a locus.
#' Default: \code{p = 1}.
#' @param r (double) Minimum percentage of individuals in a population required to process a locus for that population.
#' Default: \code{r = 0.3}.
#' @param min_maf (double) Specify a minimum minor allele frequency required to process a nucleotide site at a locus (0 < min_maf < 0.5).
#' Default: \code{min_maf = NULL}. Recommendation: use my other R package RADIATOR to filter data based on MAF.
#' @param max_obs_het (double) Specify a maximum observed heterozygosity required to process a nucleotide site at a locus.
#' Default: \code{max_obs_het = NULL}.  Recommendation: use my other R package RADIATOR to filter data based on heterozygosity.
#' @param m (integer) Specify a minimum stack depth required for individuals at a locus.
#' Default: \code{m = NULL}.
#' @param lnl_lim (integer) Filter loci with log likelihood values below this threshold.
#' Default: \code{lnl_lim = NULL}.
#' @param write_single_snp (logical) Restrict data analysis to only the first SNP per locus.
#' Default: \code{write_single_snp = FALSE}.
#' @param write_random_snp (logical) Restrict data analysis to one random SNP per locus.
#' Default: \code{write_random_snp = FALSE}.
#' @param B (character) Path to a file containing Blacklisted markers to be excluded from the export.
#' Default: \code{B = NULL}.
#' @param W (character) Path to a file containing Whitelisted markers to include in the export.
#' Default: \code{W = NULL}.

#' @param e (character) Restriction enzyme name.
#' Default: \code{e = NULL}.
#' @param merge_sites (logical) Merge loci that were produced from the same restriction enzyme cutsite (requires reference-aligned data).
#' Default: \code{merge_sites = FALSE}.
#' @param merge_prune_lim (integer) When merging adjacent loci, if at least X% samples posses both loci prune the remaining samples out of the analysis.
#' Default: \code{merge_prune_lim = NULL}.

#' @param hwe (logical) Calculate divergence from Hardy-Weinberg equilibrium
#' using the exact test at the SNP level and Guo and Thompson MCMC algorithm at
#' the haplotype level.
#' Default: \code{hwe = FALSE}.

#' @param fstats (logical) Enable SNP and haplotype-based F statistics.
#' Default: \code{fstats = FALSE}.
#' @param fst_correction (character) Specify a correction to be applied to Fst values: 'p_value', 'bonferroni_win', or 'bonferroni_gen'.
#' Default: \code{fst_correction = NULL}.
#' @param p_value_cutoff (double) maximum p-value to keep an Fst measurement.
#' Also used as base for Bonferroni correction.
#' Default: \code{p_value_cutoff = 0.05}.

#' @param k (logical) Enable kernel-smoothed Pi, Fis, Fst, Fst', and Phi_st calculations.
#' Default: \code{k = FALSE}.
#' @param smooth_fstats (logical) Enable kernel-smoothed Fst, Fst', and Phi_st calculations.
#' Default: \code{smooth_fstats = FALSE}.
#' @param smooth_popstats (logical) Enable kernel-smoothed Pi and Fis calculations.
#' Default: \code{smooth_popstats = FALSE}.

#' @param sigma (integer) Standard deviation of the kernel smoothing weight distribution.
#' Default: \code{sigma = 150000} (150kb).
#' @param bootstrap (logical) Turn on boostrap resampling for all smoothed statistics.
#' Default: \code{bootstrap = FALSE}.
#' @param N (integer) Number of bootstrap resamplings to calculate.
#' Default: \code{N = 100}.
#' @param bootstrap_pifis (logical) Turn on boostrap resampling for smoothed SNP-based Pi and Fis calculations.
#' Default: \code{bootstrap_pifis = FALSE}.
#' @param bootstrap_fst (logical) Turn on boostrap resampling for smoothed Fst calculations based on pairwise population comparison of SNPs.
#' Default: \code{bootstrap_fst = FALSE}.
#' @param bootstrap_div (logical) Turn on boostrap resampling for smoothed haplotype diveristy and gene diversity calculations based on haplotypes.
#' Default: \code{bootstrap_div = FALSE}.
#' @param bootstrap_phist (logical) Turn on boostrap resampling for smoothed Phi_st calculations based on haplotypes.
#' Default: \code{bootstrap_phist = FALSE}.
#' @param bootstrap_wl (character) Path to a whitelist file. Only use bootstrap loci contained in this whitelist.
#' Default: \code{bootstrap_wl = NULL}.


#' @param ordered_export (logical) If data is reference aligned, exports will be ordered; only a single representative of each overlapping site.
#' Default: \code{ordered_export = FALSE}.
#' @param genomic (logical) Output each nucleotide position (fixed or polymorphic) in all population members to a file (requires restriction enzyme name \code{e}).
#' Default: \code{genomic = FALSE}.
#' @param fasta_samples (logical) Output the sequences of the two haplotypes of each (diploid) sample, for each locus, in FASTA format.
#' Default: \code{fasta_samples = FALSE}.
#' @param fasta_samples_raw (logical) Output all haplotypes observed in each sample, for each locus, in FASTA format.
#' Default: \code{fasta_samples_raw = FALSE}.
#' @param fasta_loci (logical) Output consensus sequences of all loci, in FASTA format.
#' Default: \code{fasta_loci = FALSE}.
#' @param vcf  (logical) Output SNPs in Variant Call Format (VCF).
#' Default: \code{vcf = TRUE}.
#' @param genepop (logical) Output results in GenePop format.
#' Default: \code{genepop = FALSE}.
#' @param structure (logical) Output results in Structure format.
#' Default: \code{structure = FALSE}.
#' @param finestructure (logical) Output results in FineStructure/FineRADStructure format.
#' Default: \code{finestructure = FALSE}.
#' @param phase (logical) Output genotypes in PHASE format.
#' Default: \code{phase = FALSE}.
#' @param fastphase (logical) Output genotypes in fastPHASE format.
#' Default: \code{fastphase = FALSE}.
#' @param beagle (logical) Output genotypes in Beagle format.
#' Default: \code{beagle = FALSE}.
#' @param beagle_phased (logical) Output haplotypes in Beagle format.
#' Default: \code{beagle_phased = FALSE}.
#' @param plink (logical) Output genotypes in PLINK format.
#' Default: \code{plink = FALSE}.
#' @param hzar (logical) Output genotypes in Hybrid Zone Analysis using R (HZAR) format.
#' Default: \code{hzar = FALSE}.
#' @param phylip (logical) Output nucleotides that are fixed-within, and variant among populations in Phylip format for phylogenetic tree construction.
#' Default: \code{phylip = FALSE}.
#' @param phylip_var (logical) Include variable sites in the phylip output encoded using IUPAC notation.
#' Default: \code{phylip_var = FALSE}.
#' @param phylip_var_all (logical) Include all sequence as well as variable sites in the phylip output encoded using IUPAC notation.
#' Default: \code{phylip_var_all = FALSE}.
#' @param treemix (logical) Output SNPs in a format useable for the TreeMix program (Pickrell and Pritchard).
#' Default: \code{treemix = FALSE}.



#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @param verbose turn on additional logging.
#' Default: \code{verbose = FALSE}

#' @param v print program version.
#' Default: \code{v = FALSE}
#' @param log_fst_comp (logical) Log components of Fst/Phi_st calculations to a file.
#' Default: \code{log_fst_comp = FALSE}


#' @rdname run_populations
#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr mutate filter distinct
#' @importFrom purrr keep walk pwalk pmap
#' @importFrom tibble data_frame

#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/populations.php}{populations}
#' returns different output depending on arguments selected.


#' @examples
#' \dontrun{
#' pop <- stackr::run_populations(M = "population.map.tsv")
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/populations.php}{populations}

#' \href{http://catchenlab.life.illinois.edu/stacks/stacks_v2.php}{STACKS Version 2.0Beta8}


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
  P = "06_ustacks_cstacks_sstacks",
  V = NULL,
  O = "07_populations",
  M,
  t = parallel::detectCores() - 1,
  b = "guess",
  batch_size = 10000,
  # s = FALSE,
  p = 1,
  r = 0.3,
  min_maf = NULL,
  max_obs_het = NULL,
  m = NULL,
  lnl_lim = NULL,
  write_single_snp = FALSE,
  write_random_snp = FALSE,
  B = NULL,
  W = NULL,
  e = NULL,
  merge_sites = FALSE,
  merge_prune_lim = NULL,
  hwe = FALSE,
  fstats = FALSE,
  fst_correction = NULL,
  p_value_cutoff = 0.05,
  k = FALSE,
  smooth_fstats = FALSE,
  smooth_popstats = FALSE,
  sigma = 150000,
  bootstrap = FALSE,
  N = 100,
  bootstrap_pifis = FALSE,
  bootstrap_fst = FALSE,
  bootstrap_div = FALSE,
  bootstrap_phist = FALSE,
  bootstrap_wl = NULL,
  ordered_export = FALSE,
  genomic = FALSE,
  fasta_samples = FALSE,
  fasta_samples_raw = FALSE,
  fasta_loci = FALSE,
  vcf = TRUE,
  genepop = FALSE,
  structure = FALSE,
  finestructure = FALSE,
  phase = FALSE,
  fastphase = FALSE,
  beagle = FALSE,
  beagle_phased = FALSE,
  plink = FALSE,
  hzar = FALSE,
  phylip = FALSE,
  phylip_var = FALSE,
  phylip_var_all = FALSE,
  treemix = FALSE,
  h = FALSE,
  verbose = FALSE,
  v = FALSE,
  log_fst_comp = FALSE
) {

  cat("#######################################################################\n")
  cat("####################### stackr::run_populations #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists("09_log_files")) dir.create("09_log_files")
  if (!dir.exists("07_populations")) dir.create("07_populations")

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
  O <- stringi::stri_join("-O ", shQuote(output.folder))
  message("\nOutput files written in:\n", output.folder)

  V = NULL
  if (is.null(V)) {
    V <- ""
  } else {
    V <- stringi::stri_join("-V ", V)
  }

  if (missing(M)) stop("Population map required")
  M <- stringi::stri_join("-M ", M)


  # parallel.core <- t
  t <- stringi::stri_join("-t ", t)


  if (b == "guess") {
    b <- ""
  } else {
    b <- stringi::stri_join("-b ", b)
  }

  batch_size <- stringi::stri_join("--batch_size ", batch_size)

  # Data Filtering -------------------------------------------------------------
  p <- stringi::stri_join("-p ", p)
  r <- stringi::stri_join("-r ", r)

  if (is.null(min_maf)) {
    min_maf <- ""
  } else {
    min_maf <- stringi::stri_join("--min_maf ", min_maf)
  }

  if (is.null(max_obs_het)) {
    max_obs_het <- ""
  } else {
    max_obs_het <- stringi::stri_join("--max_obs_het ", max_obs_het)
  }

  if (is.null(m)) {
    m <- ""
  } else {
    m <- stringi::stri_join("-m ", m)
  }

  if (is.null(lnl_lim)) {
    lnl_lim <- ""
  } else {
    lnl_lim <- stringi::stri_join("--lnl_lim ", lnl_lim)
  }

  if (write_single_snp) {
    write_single_snp <- "--write_single_snp "
  } else {
    write_single_snp <- ""
  }

  if (write_random_snp) {
    write_random_snp <- "--write_random_snp "
  } else {
    write_random_snp <- ""
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

  if (merge_sites) {
    merge_sites <- "--merge_sites "
  } else {
    merge_sites <- ""
  }

  if (is.null(merge_prune_lim)) {
    merge_prune_lim <- ""
  } else {
    merge_prune_lim <- stringi::stri_join("--merge_prune_lim ", merge_prune_lim)
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
    if (is.null(fst_correction)) {
      fst_correction <- ""
    } else {
      fst_correction <- stringi::stri_join("--fst_correction ", fst_correction)
    }
    p_value_cutoff <- stringi::stri_join("-p_value_cutoff ", p_value_cutoff)
  } else {
    fstats <- ""
    fst_correction <- ""
    p_value_cutoff <- ""
  }





  # Kernel-smoothing algorithm -------------------------------------------------
  if (k) {
    kernel.smoothed <- TRUE
    k <- "-k "
  } else {
    kernel.smoothed <- FALSE
    k <- ""
  }

  if (smooth_fstats) {
    smooth_fstats <- TRUE
    smooth_fstats <- "--smooth_fstats "
  } else {
    smooth_fstats <- FALSE
    smooth_fstats <- ""
  }

  if (smooth_popstats) {
    smooth_popstats <- TRUE
    smooth_popstats <- "--smooth_popstats "
  } else {
    smooth_popstats <- FALSE
    smooth_popstats <- ""
  }

  if (kernel.smoothed) {
    sigma <- stringi::stri_join("--sigma ", sigma)
    if (bootstrap) {
      bootstrap <- "--bootstrap "
    } else {
      bootstrap <- ""
    }
    N <- stringi::stri_join("-N ", N)

    if (bootstrap_pifis) {
      bootstrap_pifis <- "--bootstrap_pifis "
    } else {
      bootstrap_pifis <- ""
    }

    if (bootstrap_fst) {
      bootstrap_fst <- "--bootstrap_fst "
    } else {
      bootstrap_fst <- ""
    }

    if (bootstrap_div) {
      bootstrap_div <- "--bootstrap_div "
    } else {
      bootstrap_div <- ""
    }

    if (bootstrap_phist) {
      bootstrap_phist <- "--bootstrap_phist "
    } else {
      bootstrap_phist <- ""
    }

    if (is.null(bootstrap_wl)) {
      bootstrap_wl <- ""
    } else {
      bootstrap_wl <- stringi::stri_join("--bootstrap_wl ", bootstrap_wl)
    }
  } else {
    sigma <- ""
    bootstrap <- ""
    N <- ""
    bootstrap_pifis <- ""
    bootstrap_fst <- ""
    bootstrap_div <- ""
    bootstrap_phist <- ""
    bootstrap_wl <- ""
  }

  # File output options --------------------------------------------------------
  if (ordered_export) {
    ordered_export <- "--ordered_export "
  } else {
    ordered_export <- ""
  }

  if (genomic) {
    genomic <- "--genomic "
  } else {
    genomic <- ""
  }

  if (fasta_samples) {
    fasta_samples <- "--fasta_samples "
  } else {
    fasta_samples <- ""
  }

  if (fasta_samples_raw) {
    fasta_samples_raw <- "--fasta_samples_raw "
  } else {
    fasta_samples_raw <- ""
  }

  if (fasta_loci) {
    fasta_loci <- "--fasta_loci "
  } else {
    fasta_loci <- ""
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

  if (finestructure) {
    finestructure <- "--finestructure "
  } else {
    finestructure <- ""
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

  if (beagle) {
    beagle <- "--beagle "
  } else {
    beagle <- ""
  }

  if (beagle_phased) {
    beagle_phased <- "--beagle_phased "
  } else {
    beagle_phased <- ""
  }

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

  if (phylip_var) {
    phylip_var <- "--phylip_var "
  } else {
    phylip_var <- ""
  }

  if (phylip_var_all) {
    phylip_var_all <- "--phylip_var_all "
  } else {
    phylip_var_all <- ""
  }

  if (treemix) {
    treemix <- "--treemix "
  } else {
    treemix <- ""
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

  if (log_fst_comp) {
    log_fst_comp <- "--log_fst_comp "
  } else {
    log_fst_comp <- ""
  }


  # command args ---------------------------------------------------------------
  command.arguments <- paste(
    P, V, O, M, t, b, batch_size, p, r, min_maf, max_obs_het, m, lnl_lim, write_single_snp,
    write_random_snp, B, W, e, merge_sites, merge_prune_lim, hwe, fstats, fst_correction,
    p_value_cutoff, k, smooth_fstats, smooth_popstats, sigma, bootstrap, N, bootstrap_pifis, bootstrap_fst,
    bootstrap_div, bootstrap_phist, bootstrap_wl, ordered_export, genomic,
    fasta_samples, fasta_samples_raw, fasta_loci, vcf, genepop, structure,
    finestructure, phase,
    fastphase, beagle, beagle_phased, plink, hzar, phylip, phylip_var,
    phylip_var_all, treemix, h, verbose, v, log_fst_comp
  )


  # run command ----------------------------------------------------------------
  system2(command = "populations", args = command.arguments, stderr = populations.log.file)

  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("######################## populations completed ########################\n")
  # return(res)
}# end run_populations
