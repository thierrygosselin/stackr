# Minor Allele Frequency Diagnostic
#' @title MAF diagnostic
#' @description Minor Allele Frequency diagnostic, help choose a filter threshold.
#' @param data A data frame object or file (using ".tsv")
#' of class sumstats or tidy VCF.
#' @param group.rank The number of group to class the MAF (Number) 
#' @param filename Name of the file written to the working directory (optional).
#' @rdname diagnostic_maf
#' @export
#' @import dplyr
#' @import readr
#' @details Highly recommended to look at the distribution of MAF
#' \link{plot_density_distribution_maf}.
#' @seealso \link{filter_maf}


diagnostic_maf <- function(data, group.rank, filename){
  
  LOCUS <- NULL
  POP_ID <- NULL
  FREQ_ALT <- NULL
  RANK <- NULL
  GLOBAL_MAF <- NULL
  MAF_P <- NULL
  MAF_L <- NULL
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = TRUE)
    message("Using the file in your directory")
    
  } else {
    data <- data
    message("Using the file from your global environment")
  }
  
  # Local 
  test.local <- data %>%
    select(LOCUS, POP_ID, FREQ_ALT) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      MAF_P = min(FREQ_ALT, na.rm = TRUE)
    ) %>% 
    group_by(LOCUS) %>%
    summarise(MAF_L = mean(MAF_P, na.rm =TRUE)) %>%
    group_by(RANK = ntile(MAF_L, group.rank)) %>%
    summarise(
      LOCAL_MAF = mean(MAF_L, na.rm = T),
      n = length(LOCUS)
    ) %>%
    select(-n)
  
  # Global
  
  test.global <- data %>%
    select(LOCUS, POP_ID, GLOBAL_MAF) %>%
    group_by(LOCUS) %>%
    summarise(GLOBAL_MAF = mean(GLOBAL_MAF, na.rm =TRUE)) %>%
    group_by(RANK = ntile(GLOBAL_MAF, group.rank)) %>%
    summarise(
      GLOBAL_MAF = mean(GLOBAL_MAF, na.rm = T),
      n = length(LOCUS)
    ) %>%
    select(GLOBAL_MAF, n)
  
  maf.diagnostic <- bind_cols(test.local, test.global)
  
  if (missing(filename) == "FALSE") {
    message("Saving the table in your working directory...")
    write_tsv(maf.diagnostic, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
  }
  
  
  
  invisible(cat(sprintf(
    "%s\n
Working directory:
%s",
    saving, getwd()
  )))
  
  return(maf.diagnostic)
  
}
