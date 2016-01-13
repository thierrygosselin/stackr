# Observed heterozygosity 
#' @title Heterozygosity filter
#' @description Observed heterozygosity filter from a sumstats prepared file
#' or a tidy VCF file.
#' @param data A data frame object or file (using ".tsv")
#' of class sumstats or a tidy VCF summarised.
#' @param approach Character. By \code{"SNP"} or by \code{"haplotype"}. 
#' The function will consider the SNP or haplotype MAF statistics to filter the marker. 
#' Default by \code{"haplotype"}.
#' @param het.threshold Number (0 - 0.5)
#' @param het.diff.threshold Number (0 - 1).
#' @param pop.threshold Fixed number of pop required to keep the locus.
#' @param percent Is the threshold a percentage ? TRUE or FALSE.
#' @param filename Name of the file written to the working directory (optional).
#' @rdname filter_het
#' @export
#' @import stringi
#' @import dplyr
#' @import readr
#' @seealso \link{plot_density_distribution_het}


filter_het <- function(data, approach = "haplotype", het.threshold, het.diff.threshold, pop.threshold, percent, filename) {
  
  POP_ID <- NULL
  LOCUS <- NULL
  HET_DIFF <- NULL
  HET_O <- NULL
  HET_MAX <- NULL
  
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
    message("Using the file in your directory")
  } else {
    data <- data
    message("Using the file from your global environment")
  }
  
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message("Using a proportion threshold...")
    threshold.id <- "of proportion"
  } else if (stri_detect_fixed(percent, "T")) {
    multiplication.number <- 100/pop.number
    message("Using a percentage threshold...")
    threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message("Using a fixed threshold...")
    threshold.id <- "population as a fixed"
  }
  
  if (missing(approach) | approach == "haplotype"){
    message("Approach selected: haplotype")
    het.filter <- data %>%
      select(LOCUS, POS, POP_ID, HET_O) %>%
      group_by (LOCUS, POP_ID) %>%
      summarise(
        HET_DIFF = max(HET_O) - min(HET_O),
        HET_MAX = max(HET_O)
      ) %>%
      filter(HET_DIFF <= het.diff.threshold & HET_MAX <= het.threshold) %>%  
      group_by(LOCUS) %>%
      tally() %>%
      filter((n * multiplication.number) >= pop.threshold) %>%
      select(LOCUS) %>%
      left_join(data, by="LOCUS") %>%
      arrange(LOCUS, POP_ID)
  } else {
    message("Approach selected: SNP")
    het.filter <- data %>%
      select(LOCUS, POS, POP_ID, HET_O) %>%
      group_by(LOCUS, POS, POP_ID) %>%
      summarise(
        HET_DIFF = max(HET_O) - min(HET_O),
        HET_MAX = max(HET_O)
      ) %>%
      filter(HET_DIFF <= het.diff.threshold & HET_MAX <= het.threshold) %>%  
      group_by(LOCUS, POS) %>%
      tally() %>%
      filter((n * multiplication.number) >= pop.threshold) %>%
      select(LOCUS, POS) %>%
      left_join(data, by = c("LOCUS", "POS")) %>%
      arrange(LOCUS, POS, POP_ID)
  }
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(het.filter, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
  }
  
  invisible(cat(sprintf(
    "HET filter:

Markers with an observed heterozygosity <= %s and with a difference of <= %s
between the max and min observed heterozygosity along the reads
(independent of the number of SNP/reads)
in at least %s percent of sites/pop are kept\n
The number of SNP removed by the HET filter = %s SNP
The number of LOCI removed by the HET filter = %s LOCI
The number of SNP before -> after the HET filter: %s -> %s SNP
The number of LOCI before -> after the HET filter: %s -> %s LOCI\n
%s\n
Working directory:
%s",
    het.threshold, het.diff.threshold, pop.threshold,
    n_distinct(data$POS)-n_distinct(het.filter$POS),
    n_distinct(data$LOCUS)-n_distinct(het.filter$LOCUS),
    n_distinct(data$POS), n_distinct(het.filter$POS),
    n_distinct(data$LOCUS), n_distinct(het.filter$LOCUS),
    saving, getwd()
  )))
  return(het.filter)
}







