#' @title Individual filter
#' @description Filter individuals number, percentage or fix threshold
#' from a tidy VCF file.
#' @param tidy.vcf A tidy vcf object or file (using ".tsv").
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An character string with your populations ordered.
#' @param pop.labels An optional character string with new populations names.
#' @param ind.threshold The individual threshold, proportion, percentage or 
#' number e.g. 0.70, 70, 15.
#' @param percent Is the threshold a percentage? TRUE or FALSE.
#' This argument is necessary to distinguish percentage from integer individual
#' threshold (e.g. 70 percent or 70 individuals).
#' @param filename Name of the file written to the working directory.
#' @rdname filter_individual
#' @export
#' @import stringi
#' @import dplyr
#' @import readr

filter_individual <- function(tidy.vcf, pop.id.start, pop.id.end, pop.levels, pop.labels, ind.threshold, percent, filename) {
  
  X1 <- NULL
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  N <- NULL
  LOCUS_N <- NULL
  IND <- NULL
  POP <- NULL
  N_TOT <- NULL
  N_SUM <- NULL
  N_MIN <- NULL
  GT <- NULL
  N_VCF <- NULL
  POP_ID <- NULL
  
  if (is.vector(tidy.vcf) == "TRUE") {
    data <- read_tsv(tidy.vcf, col_names = TRUE)
    message("Using the file in your directory")
    
  } else {
    data <- tidy.vcf
    message("Using the tidy.vcf from your global environment")
  }
  
  # tidy VCF 
  if (stri_detect_fixed(percent, "F")) {
    message("Using a fixed threshold..")
    threshold.id <- "individuals as a integer"
    
    ind.filter <- data %>%
      group_by(LOCUS, POS, POP_ID) %>%
      summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
      group_by(LOCUS, POP_ID) %>%
      summarise(N_VCF = ceiling(mean(N_VCF, na.rm = TRUE))) %>%
      group_by(LOCUS) %>%
      summarise (
        POP = length (POP_ID),
        IND = length (N_VCF[N_VCF >= ind.threshold])
      ) %>%
      group_by(LOCUS) %>%
      filter(round(((IND/POP)*100),0) >= 100) %>%
      select(LOCUS) %>%
      left_join(data, by="LOCUS") %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)
      ) %>%
      arrange(POP_ID, LOCUS, POS)
    
  } else if (stri_detect_fixed(ind.threshold, ".") & ind.threshold < 1) {
    
    message("Using a proportion threshold...")
    threshold.id <- " of proportion"
    
    ind.filter <- data %>%
      group_by(POP_ID) %>%
      summarise(N_TOT = as.numeric(n_distinct(INDIVIDUALS))) %>%
      mutate(N_MIN = floor(N_TOT * as.numeric(ind.threshold))) %>%
      full_join(
        data %>%
          group_by(LOCUS, POS, POP_ID) %>%
          summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
          group_by(LOCUS, POP_ID) %>%
          summarise(N_VCF = ceiling(mean(N_VCF, na.rm = TRUE))),
        by = c("POP_ID")
      )%>%
      group_by(LOCUS) %>%
      summarise(
        POP = as.numeric(length (POP_ID)),
        IND = as.numeric(length (N_VCF[N_VCF >= N_MIN]))
      ) %>%
      group_by(LOCUS) %>%
      filter(round(((IND/POP)*100),0) >= 100) %>%
      select(LOCUS) %>%
      left_join(data, by="LOCUS") %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)
      ) %>%
      arrange(POP_ID, LOCUS, POS)
    
  } else {
    
    message("Using a percentage threshold...")
    threshold.id <- "percent"
    
    ind.filter <- data %>%
      group_by(POP_ID) %>%
      summarise(N_TOT = as.numeric(n_distinct(INDIVIDUALS))) %>%
      mutate(N_MIN = floor(N_TOT * as.numeric(ind.threshold / 100))) %>%
      full_join(
        data %>%
          group_by(LOCUS, POS, POP_ID) %>%
          summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
          group_by(LOCUS, POP_ID) %>%
          summarise(N_VCF = ceiling(mean(N_VCF, na.rm = T))),
        by = c("POP_ID")
      )%>%
      group_by(LOCUS) %>%
      summarise(
        POP = as.numeric(length (POP_ID)),
        IND = as.numeric(length (N_VCF[N_VCF >= N_MIN]))
      ) %>%
      group_by(LOCUS) %>%
      filter(round(((IND/POP)*100),0) >= 100) %>%
      select(LOCUS) %>%
      left_join(data, by="LOCUS") %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)
      ) %>%
      arrange(POP_ID, LOCUS, POS)
  }
  
if (missing(filename) == "FALSE") {
  message("Saving the file in your working directory...")
  write_tsv(ind.filter, filename, append = FALSE, col_names = TRUE)
  saving <- paste("Saving was selected, the filename:", filename, sep = " ")
} else {
  saving <- "Saving was not selected"
}

invisible(cat(sprintf(
  "Individual filter: %s %s threshold of genotyped individuals per sampling sites to keep the marker
  The number of SNP removed by the individual filter = %s SNP
  The number of LOCI removed by the individual filter = %s LOCI
  The number of SNP before -> after the individual filter: %s -> %s SNP
  The number of LOCI before -> after the individual filter: %s -> %s LOCI\n
  %s\n
  Working directory:
  %s", 
  ind.threshold,
  threshold.id,
  n_distinct(data$POS)-n_distinct(ind.filter$POS),
  n_distinct(data$LOCUS)-n_distinct(ind.filter$LOCUS),
  n_distinct(data$POS), n_distinct(ind.filter$POS),
  n_distinct(data$LOCUS), n_distinct(ind.filter$LOCUS),
  saving, getwd()
)))
return(ind.filter)
}
