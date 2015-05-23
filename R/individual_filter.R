#' @title Individual filter
#' @description Filter individuals number, percentage or fix threshold
#' in a sumstats prepared file or a tidy VCF file.
#' @param data A data frame object or file (using ".tsv")
#' of class sumstats or tidy VCF.
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

individual_filter <- function(data, population.map, pop.id.start, pop.id.end, pop.levels, ind.threshold, threshold.fixed ) {
  
  if (is.vector(population.map) == "TRUE") {
  
    population.map <- read_tsv(population.map, col_names = F) %>%
      select(INDIVIDUALS=X1) %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
        )

  } else {
  
    population.map <- population.map %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
        )
  }
  
  if (stri_detect_fixed(data, "sumstats")) {
    
    if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
    } else {
      data <- data
    }

    if (stri_detect_fixed(threshold.fixed, "T")) {
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
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
        ) %>%
      arrange(LOCUS)
        
   
    } else if (stri_detect_fixed(ind.threshold, ".") == "TRUE") {
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
          POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
          ) %>%
        arrange(LOCUS)
   
    } else {
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
          POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
          ) %>%
        arrange(LOCUS)
    }
    
  } else {
    
    if (is.vector(data) == "TRUE") {
      stacks.vcf.file <- read_tsv(data, col_names = T)
      } else {
        stacks.vcf.file <- data
      }
    
    if (stri_detect_fixed(threshold.fixed, "T")) {
    
      ind.filter <- stacks.vcf.file %>%
        group_by(LOCUS, POS, POP_ID) %>%
        summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
        group_by(LOCUS, POP_ID) %>%
        summarise(N_VCF = ceiling(mean(N_VCF, na.rm = T))) %>%
        group_by(LOCUS) %>%
        summarise (
          POP = length (POP_ID),
          IND = length (N_VCF[N_VCF >= ind.threshold])
          ) %>%
        group_by(LOCUS) %>%
        filter(round(((IND/POP)*100),0) >= 100) %>%
        select(LOCUS) %>%
        left_join(stacks.vcf.file, by="LOCUS") %>%
        mutate(
          POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
          ) %>%
        arrange(POP_ID, LOCUS, POS)
    
    } else if (stri_detect_fixed(ind.threshold, ".") == "TRUE") {
          
      ind.filter <- stacks.vcf.file %>%
        group_by(POP_ID) %>%
        summarise(N_TOT = as.numeric(n_distinct(INDIVIDUALS))) %>%
        mutate(N_MIN = floor(N_TOT * as.numeric(ind.threshold))) %>%
        full_join(
          stacks.vcf.file %>%
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
        left_join(stacks.vcf.file, by="LOCUS") %>%
        mutate(
          POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
          ) %>%
        arrange(POP_ID, LOCUS, POS)
      
    } else {
      
      ind.filter <- stacks.vcf.file %>%
        group_by(POP_ID) %>%
        summarise(N_TOT = as.numeric(n_distinct(INDIVIDUALS))) %>%
        mutate(N_MIN = floor(N_TOT * as.numeric(ind.threshold / 100))) %>%
        full_join(
          stacks.vcf.file %>%
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
        left_join(stacks.vcf.file, by="LOCUS") %>%
        mutate(
          POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
          ) %>%
        arrange(POP_ID, LOCUS, POS)
    }
  }
  
invisible(
        cat(
         sprintf(
"Individual filter: %s percent or proportion of genotyped individuals per sampling sites to keep the marker
The number of SNP removed by the individual filter = %s SNP
The number of LOCI removed by the individual filter = %s LOCI
The number of SNP after the individual filter = %s SNP
The number of LOCI after the individual filter = %s LOCI", 
ind.threshold,
n_distinct(data$POS)-n_distinct(ind.filter$POS),
n_distinct(data$LOCUS)-n_distinct(ind.filter$LOCUS),
n_distinct(ind.filter$POS),
n_distinct(ind.filter$LOCUS)
        )
        )
      )
ind.filter
}
