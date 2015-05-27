# Minor Allele Frequency
#' @title MAF filter
#' @description Minor Allele Frequency filter from a sumstats prepared file
#' or a tidy VCF file.
#' @param data A data frame object or file (using ".tsv")
#' of class sumstats or a tidy VCF summarised.
#' @param local.maf.threshold Number.
#' @param global.maf.threshold Number.
#' @param pop.threshold Fixed number of pop required to keep the locus.
#' @param filename Name of the file written to the working directory (optional).
#' @rdname filter_maf
#' @export
#' @import stringi
#' @import dplyr
#' @import readr
#' @details To help choose a threshold for the local and global MAF, look
#' at the function \link{diagnostic_maf}
#' @seealso \link{summary_stats_vcf_tidy, sumstats_prep}


filter_maf <- function(data, local.maf.threshold, global.maf.threshold, pop.threshold, filename) {
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
    message("Using the file in your directory")
  } else {
    data <- data
    message("Using the file from your global environment")
  }
  
  
  
  maf.filter <- data %>%
    group_by(LOCUS, POP_ID) %>%
    #     summarise(MAF_DIFF = max(FREQ_ALT) - min(FREQ_ALT)) %>% 
    #     filter(MAF_DIFF <= maf.diff.threshold) %>%
    #     group_by(LOCUS) %>%
    #     tally() %>%
    #     filter((n * multiplication.number) >= pop.maf.diff.threshold) %>%
    #     select(LOCUS) %>%
    #     left_join(data, by="LOCUS") %>%
    #     group_by(LOCUS, POP_ID) %>%
    summarise(
      GLOBAL_MAF = mean(GLOBAL_MAF, na.rm = TRUE),  
      LOCAL_MAF = min(FREQ_ALT, na.rm = TRUE)
    ) %>% 
    filter(LOCAL_MAF >= local.maf.threshold | GLOBAL_MAF >= global.maf.threshold) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter(n >= pop.threshold) %>% 
    select(LOCUS) %>%
    left_join(data, by = "LOCUS") %>%
    arrange(LOCUS, POP_ID)
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(maf.filter, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
  }
  
  
  invisible(cat(sprintf(
    "MAF filter:

To pass the filter, markers with a local MAF >= %s or a global MAF >= %s, 
in at least %s sampling sites/pop.
The number of SNP removed by the MAF filter = %s SNP
The number of LOCI removed by the MAF filter = %s LOCI
The number of SNP before -> after the MAF filter: %s -> %s SNP
The number of LOCI before -> after the MAF filter: %s -> %s LOCI\n
%s\n
Working directory:
%s",
    local.maf.threshold,
    global.maf.threshold,
    pop.threshold,
    n_distinct(data$POS)-n_distinct(maf.filter$POS),
    n_distinct(data$LOCUS)-n_distinct(maf.filter$LOCUS),
    n_distinct(data$POS), n_distinct(maf.filter$POS),
    n_distinct(data$LOCUS), n_distinct(maf.filter$LOCUS),
    saving, getwd()
  )))
  
  return(maf.filter)
  
}
