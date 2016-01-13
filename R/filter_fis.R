# Fis 
#' @title Fis filter
#' @description Inbreeding coefficient filter from a sumstats prepared file
#' or a tidy VCF file.
#' @param data A data frame object or file (using ".tsv")
#' of class sumstats or a tidy VCF summarised.
#' @param approach Character. By \code{"SNP"} or by \code{"haplotype"}. 
#' The function will consider the SNP or haplotype MAF statistics to filter the marker. 
#' Default by \code{"haplotype"}.
#' @param fis.min.threshold Number
#' @param fis.max.threshold Number.
#' @param fis.diff.threshold Number (0 - 1)
#' @param pop.threshold Fixed number of pop required to keep the locus.
#' @param percent Is the threshold a percentage ? TRUE or FALSE.
#' @param filename Name of the file written to the working directory (optional).
#' @rdname filter_fis
#' @export
#' @import stringi
#' @import dplyr
#' @import readr


filter_fis <- function(data, approach = "haplotype", fis.min.threshold, fis.max.threshold, fis.diff.threshold, pop.threshold, percent, filename) {
  
  POP_ID <- NULL
  LOCUS <- NULL
  FIS <- NULL
  FIS_MIN <- NULL
  FIS_MAX <- NULL
  FIS_DIFF <- NULL
  
  
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
    fis.filter <- data %>%
      select(LOCUS, POS, POP_ID, FIS) %>%
      group_by (LOCUS, POP_ID) %>%
      summarise(
        FIS_MIN = min(FIS),
        FIS_MAX = max(FIS),
        FIS_DIFF = FIS_MAX-FIS_MIN
      ) %>%
      filter(FIS_MIN >= fis.min.threshold) %>%
      filter(FIS_MAX <= fis.max.threshold) %>%
      filter(FIS_DIFF <= fis.diff.threshold) %>%
      group_by(LOCUS) %>%
      tally() %>%
      filter((n * multiplication.number) >= pop.threshold) %>%
      select(LOCUS) %>%
      left_join(data, by="LOCUS") %>%
      arrange(LOCUS, POP_ID)
  } else {
    message("Approach selected: SNP")
    fis.filter <- data %>%
      select(LOCUS, POS, POP_ID, FIS) %>%
      group_by(LOCUS, POS, POP_ID) %>%
      summarise(
        FIS_MIN = min(FIS),
        FIS_MAX = max(FIS),
        FIS_DIFF = FIS_MAX-FIS_MIN
      ) %>%
      filter(FIS_MIN >= fis.min.threshold) %>%
      filter(FIS_MAX <= fis.max.threshold) %>%
      filter(FIS_DIFF <= fis.diff.threshold) %>%
      group_by(LOCUS, POS) %>%
      tally() %>%
      filter((n * multiplication.number) >= pop.threshold) %>%
      select(LOCUS, POS) %>%
      left_join(data, by = c("LOCUS", "POS")) %>%
      arrange(LOCUS, POS, POP_ID)
}
  
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(fis.filter, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
  }
  
  
  invisible(cat(sprintf(
    "Fis filter:
Fis min >= %s or Fis max <= %s or 
difference along the read/haplotype between the max and min Fis > %s,
all in %s percent of the sampling sites/pop were removed\n
The number of SNP removed by the Fis filter = %s SNP
The number of LOCI removed by the Fis filter = %s LOCI
The number of SNP before -> after the Fis filter: %s -> %s SNP
The number of LOCI before -> after the Fis filter: %s -> %s LOCI\n
%s\n
Working directory:
%s",
    fis.min.threshold, fis.max.threshold, fis.diff.threshold, pop.threshold,
    n_distinct(data$POS)-n_distinct(fis.filter$POS),
    n_distinct(data$LOCUS)-n_distinct(fis.filter$LOCUS),
    n_distinct(data$POS), n_distinct(fis.filter$POS),
    n_distinct(data$LOCUS), n_distinct(fis.filter$LOCUS),
    saving, getwd()
  )))
  return(fis.filter)
}
