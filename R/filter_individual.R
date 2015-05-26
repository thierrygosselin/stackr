#' @title Individual filter
#' @description Filter individuals number, percentage or fix threshold
#' in a sumstats prepared file or a tidy VCF file.
#' @param data A data frame object or file (using ".tsv")
#' of class sumstats or tidy VCF.
#' @param is.vcf Is the data a tidy vcf file ? TRUE or FALSE.
#' @param population.map The population map or individuals listed in one column.
#' No headers.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An optional character string with your populations ordered.
#' @param ind.threshold The individual threshold, proportion, percentage or 
#' number e.g. 0.70, 70, 15.
#' @param threshold.fixed Is the threshold fixed ? TRUE or FALSE.
#' @param filename Name of the file written to the working directory.
#' @rdname individual_filter
#' @export
#' @import stringi

filter_individual <- function(data, is.vcf, population.map, pop.id.start, pop.id.end, pop.levels, ind.threshold, threshold.fixed, filename) {
  
  if (is.vector(population.map) == "TRUE") {
    
    message("Using the population map in your directory")
    
    population.map <- read_tsv(population.map, col_names = FALSE) %>%
      select(INDIVIDUALS=X1) %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
      )
    
  } else {
    
    message("Using the population map from your global environment")
    
    population.map <- population.map %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =TRUE)
      )
  }
  if (stri_detect_fixed(is.vcf, "F")) {
    
    if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = TRUE)
      message("Using the sumstats file in your directory")
      
    } else {
      data <- data
      message("Using the sumstats from your global environment")
    }
    
    if (stri_detect_fixed(threshold.fixed, "T")) {
      
      message("Using a fixed threshold")
      threshold.id <- "individuals as a fixed"
      
      ind.filter <- data %>%
        group_by (LOCUS, POP_ID) %>%
        summarise(LOCUS_N = ceiling(mean(N))) %>%
        group_by(LOCUS) %>%
        summarise (
          POP = as.numeric(length (POP_ID)),
          IND = as.numeric(length (LOCUS_N[LOCUS_N >= ind.threshold]))
        ) %>%
        group_by(LOCUS) %>%
        filter(round(((IND/POP)*100),0) >= 100) %>%
        select(LOCUS) %>%
        left_join(data, by="LOCUS") %>%
        mutate(
          POP_ID = factor(POP_ID, levels = pop.levels, ordered =TRUE)
        ) %>%
        arrange(LOCUS)
      
    } else if (stri_detect_fixed(ind.threshold, ".") == "TRUE") {
      
      message("Using a proportion threshold")
      threshold.id <- "of proportion"
      
      ind.filter <- population.map %>%
        group_by(POP_ID) %>%
        summarise(
          N_TOT = as.numeric(length(INDIVIDUALS)),
          N_MIN = floor((N_TOT * as.numeric(ind.threshold)))
        ) %>%
        left_join(
          data %>%
            select(LOCUS, POS, POP_ID, N) %>%
            group_by(LOCUS, POP_ID) %>%
            summarise(N_SUM = ceiling(mean(N, na.rm = TRUE))),
          by="POP_ID") %>%
        group_by(LOCUS) %>%
        summarise (
          POP = as.numeric(length (POP_ID)),
          IND = as.numeric(length (N_SUM[N_SUM >= N_MIN]))
        ) %>%
        group_by(LOCUS) %>%
        filter(round(((IND/POP)*100),0) >= 100) %>%
        select(LOCUS) %>%
        left_join(data, by="LOCUS") %>%
        mutate(
          POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)
        ) %>%
        arrange(LOCUS)
      
    } else {
      
      message("Using a percentage threshold")
      threshold.id <- "percent"
      
      ind.filter <- population.map %>%
        group_by(POP_ID) %>%
        summarise(
          N_TOT = as.numeric(length(INDIVIDUALS)),
          N_MIN = floor((N_TOT * as.numeric(ind.threshold)/100))
        ) %>%
        left_join(
          data %>%
            select(LOCUS, POS, POP_ID, N) %>%
            group_by(LOCUS, POP_ID) %>%
            summarise(N_SUM = ceiling(mean(N, na.rm = T))),
          by="POP_ID") %>%
        group_by(LOCUS) %>%
        summarise (
          POP = as.numeric(length (POP_ID)),
          IND = as.numeric(length (N_SUM[N_SUM >= N_MIN]))
        ) %>%
        group_by(LOCUS) %>%
        filter(round(((IND/POP)*100),0) >= 100) %>%
        select(LOCUS) %>%
        left_join(data, by="LOCUS") %>%
        mutate(
          POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)
        ) %>%
        arrange(LOCUS)
    }
    
  } else {
    
    if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      message("Using the tidy vcf file in your directory")
      
    } else {
      data <- data
      message("Using the tidy vcf from your global environment")
    }
    
    if (stri_detect_fixed(threshold.fixed, "T")) {
      message("Using a fixed threshold")
      threshold.id <- "individuals as a fixed"
      
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
      
    } else if (stri_detect_fixed(ind.threshold, ".") == "TRUE") {
      
      message("Using a proportion threshold")
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
      
      message("Using a percentage threshold")
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
  ind.filter
}

