# Detect allele problems

#' @name detect_allele_problems

#' @title Detect allele problems

#' @description Detect allele problems.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users. The function computes allele counts and
#' looks for alleles bellow a certain threshold. The summary statistics for the
#' markers with problematic allele is computed based on coverage and
#' genotype likelihood. This function is very fast to highlight: i) bias in
#' representation of allelic copies and unequal coverage and ii) type I error 
#' during genotyping of heterozygote with low coverage data.

#' @inheritParams tidy_genomic_data

#' @param allele.threshold Threshold of allele copies. Below this threshold
#' markers with those allele are blacklisted and summary statistics
#' for read, ref and alt depth and genotype likelihood information is summarized.

#' @return A list with summary, plot and blacklist.


#' @details 
#' \strong{Input files:} see \pkg{stackr} \code{\link[stackr]{tidy_genomic_data}}
#' for detailed information about supported file format.
#' 
#' \strong{Under developement, use with caution}

#' @examples
#' \dontrun{
#' problem <- detect_allele_problems(
#' data = "batch_1.vcf",
#' strata = "strata.salmon.tsv",
#' allele.threshold = 3
#' )
#' # The default with this function:
#' # vcf.metadata = TRUE # VCF metadata (read depth, etc. are imported)
#' # monomorphic.out = TRUE # = discarded
#' # common.markers = TRUE # markers not in common between pop are discarded
#' }

#' @export
#' @rdname detect_allele_problems
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join summarise_each_ funs
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub
#' @importFrom tibble has_name
#' @importFrom tidyr gather
#' @importFrom readr write_tsv

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_allele_problems <- function(
  data,
  allele.threshold = 3,
  vcf.metadata = TRUE,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  blacklist.genotype = NULL,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  max.marker = NULL,
  blacklist.id = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  pop.select = NULL,
  filename = NULL,
  verbose = FALSE
) {
  
  cat("#######################################################################\n")
  cat("################## stackr::detect_allele_problems #####################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  
  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  
  
  # Import data ---------------------------------------------------------------
  input <- stackr::tidy_genomic_data(
    data = data,
    vcf.metadata = TRUE,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = common.markers,
    maf.thresholds = maf.thresholds,
    maf.pop.num.threshold = maf.pop.num.threshold,
    maf.approach = maf.approach,
    maf.operator = maf.operator,
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    filename = NULL
  )  
  
  input$GT <- stringi::stri_replace_all_fixed(
    str = input$GT,
    pattern = c("/", ":", "_", "-", "."),
    replacement = c("", "", "", "", ""),
    vectorize_all = FALSE
  )
  
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input), 
      pattern = "GENOTYPE", 
      replacement = "GT", 
      vectorize_all = FALSE)
  }
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # Check that coverage info is there
  if (!tibble::has_name(input, "READ_DEPTH")) stop("This function requires read depth metadata")
  if (!tibble::has_name(input, "ALLELE_ALT_DEPTH")) stop("This function requires allele depth metadata")
  if (!tibble::has_name(input, "GL")) stop("This function requires genotype likelihood metadata")
  
  
  message("\nComputing allele count...")
  
  # Allele count ---------------------------------------------------------------
  allele.count <- input %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::select(MARKERS, INDIVIDUALS, GT) %>% 
    dplyr::mutate(
      A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
      A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
    ) %>% 
    dplyr::select(-GT) %>% 
    tidyr::gather(
      data = .,
      key = ALLELES,
      value = GT, 
      -c(MARKERS, INDIVIDUALS)
    ) %>% 
    dplyr::arrange(MARKERS, INDIVIDUALS, GT) %>%
    dplyr::group_by(MARKERS, GT) %>% 
    dplyr::tally(.) %>% # count alleles, longest part of the block
    dplyr::ungroup(.) %>% 
    dplyr::group_by(MARKERS) %>% 
    dplyr::filter(n == (min(n))) %>% 
    dplyr::arrange(n)
  
  # Histogram of allele counts--------------------------------------------------
  allele.count.plot <- ggplot2::ggplot(allele.count, aes(n)) +
    geom_bar() +
    labs(x = "Copie(s) of alternate allele") +
    labs(y = "Distribution (number)") +
    # scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1),
    # labels = c("0", "0.05", "0.1", "0.2", "0.5", "1.0")) +
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.text.x = element_text(size = 8, family = "Helvetica")#, angle = 90, hjust = 1, vjust = 0.5) 
    )
  # allele.count.plot
  
  # Summary stats and Blacklist ------------------------------------------------
  problem <- allele.count %>%
    dplyr::rename(ALLELE_SUM = n) %>% 
    dplyr::group_by(ALLELE_SUM) %>% 
    dplyr::tally(.) %>% # count alleles, longest part of the block
    dplyr::ungroup(.) %>% 
    dplyr::filter(ALLELE_SUM <= allele.threshold)
  
  blacklist <- dplyr::filter(.data = allele.count, n <= allele.threshold) %>%
    dplyr::select(MARKERS, ALLELE_COPIES = n)
  
  input.info.needed <- input %>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, READ_DEPTH, ALLELE_REF_DEPTH, ALLELE_ALT_DEPTH, GL) %>%
    dplyr::group_by(MARKERS, CHROM, LOCUS, POS, REF, ALT) %>% 
    dplyr::summarise(
      READ_DEPTH_MEAN = mean(READ_DEPTH, na.rm = TRUE),
      READ_DEPTH_MIN = min(READ_DEPTH, na.rm = TRUE),
      READ_DEPTH_MAX = max(READ_DEPTH, na.rm = TRUE),
      ALLELE_REF_DEPTH_MEAN = mean(ALLELE_REF_DEPTH, na.rm = TRUE),
      ALLELE_REF_DEPTH_MIN = min(ALLELE_REF_DEPTH, na.rm = TRUE),
      ALLELE_REF_DEPTH_MAX = max(ALLELE_REF_DEPTH, na.rm = TRUE),
      ALLELE_ALT_DEPTH_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = TRUE),
      ALLELE_ALT_DEPTH_MIN = min(ALLELE_ALT_DEPTH, na.rm = TRUE),
      ALLELE_ALT_DEPTH_MAX = max(ALLELE_ALT_DEPTH, na.rm = TRUE),
      GL_MEAN = mean(GL, na.rm = TRUE),
      GL_MIN = min(GL, na.rm = TRUE),
      GL_MAX = max(GL, na.rm = TRUE)
    )
  
  blacklist.summary <- dplyr::left_join(blacklist, input.info.needed, by = "MARKERS") %>% 
    dplyr::arrange(ALLELE_COPIES, CHROM, LOCUS, POS)
  
  allele.summary <- dplyr::left_join(blacklist, input, by = "MARKERS") %>% 
    dplyr::filter(GT_VCF %in% c("1/1", "1/0", "0/1")) %>% 
    dplyr::select(-GT_BIN, -GT) %>% 
    dplyr::left_join(
      blacklist.summary %>% 
        dplyr::select(MARKERS, READ_DEPTH_MEAN, READ_DEPTH_MIN, READ_DEPTH_MAX, ALLELE_REF_DEPTH_MEAN, ALLELE_REF_DEPTH_MIN, ALLELE_REF_DEPTH_MAX, ALLELE_ALT_DEPTH_MEAN, ALLELE_ALT_DEPTH_MIN, ALLELE_ALT_DEPTH_MAX, GL_MEAN, GL_MIN, GL_MAX),
      , by = "MARKERS"
    ) %>% 
    dplyr::select(
      ALLELE_COPIES, MARKERS, CHROM, LOCUS, POS, POP_ID, INDIVIDUALS, GT_VCF, REF, ALT, 
      READ_DEPTH, READ_DEPTH_MEAN, READ_DEPTH_MIN, READ_DEPTH_MAX, 
      ALLELE_REF_DEPTH, ALLELE_REF_DEPTH_MEAN, ALLELE_REF_DEPTH_MIN, ALLELE_REF_DEPTH_MAX,
      ALLELE_ALT_DEPTH, ALLELE_ALT_DEPTH_MEAN, ALLELE_ALT_DEPTH_MIN, ALLELE_ALT_DEPTH_MAX,
      GL, GL_MEAN, GL_MIN, GL_MAX)
  
  readr::write_tsv(x = allele.summary, path = "problematic.allele.summary.stats.tsv")
  message("Written to the working directory:\nproblematic.allele.summary.stats.tsv")
  
  blacklist.markers <- dplyr::ungroup(allele.summary) %>% 
    dplyr::distinct(MARKERS, CHROM, LOCUS, POS)
  
  readr::write_tsv(x = blacklist.markers, path = "blacklist.markers.allele.problem.tsv")
  message("Written to the working directory:\nblacklist.markers.allele.problem.tsv")
  
  
  cat("############################### RESULTS ###############################\n")
  res = list(
    allele.summary = allele.summary,
    allele.count.distribution.plot = allele.count.plot,
    number.markers = problem,
    blacklist.markers = blacklist.markers
  )
  message(paste0("Number of markers with allele copies below threshold: ", nrow(blacklist.markers)))
  timing <- proc.time() - timing
  message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
  cat("############################## completed ##############################\n")
  return(res)
}
