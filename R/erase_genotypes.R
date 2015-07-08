#' @title Blacklist erase genotypes
#' @description This function creates a blacklist of loci and individual
#' genotypes useful for the function erase_genotype.
#' @param tidy.vcf.file A data frame object or file (using ".tsv")
#' of class tidy VCF.
#' @param read.depth.threshold Threshold number.
#' @param allele.depth.threshold Threhold number for the min depth 
#' of REF or ALT alleles.
#' @param allele.imbalance.threshold Working on it.
#' @param filename Name of the file written to the working directory.
#' @details 1. Keep the genotypes below the read depth and the 
#' genotype likelihood threshold. 2. Keep the REF and ALT allele below the
#' coverage threshold. 3. Create a blacklist of genotypes to erase.
#' @return The function returns the individuals genotypes, by loci and 
#' individuals, to erase.
#' @rdname blacklist_erase_genotype.
#' @export
#' @import dplyr
#' @import readr


blacklist_erase_genotype <- function(tidy.vcf.file, read.depth.threshold, allele.depth.threshold, allele.imbalance.threshold, filename) {
  
  GT <- NULL
  GL <- NULL
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL
  MIN_REF <- NULL
  MIN_ALT <- NULL
  READ_DEPTH <- NULL
  ALLELE_COVERAGE_RATIO <- NULL
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    tidy.vcf.file <- read_tsv(tidy.vcf.file, col_names = T)
    message("Using the tidy vcf file in your directory")
  } else {
    tidy.vcf.file <- tidy.vcf.file
    message("Using the tidy vcf file from your global environment")
  }
  
  blacklist <- tidy.vcf.file %>%
    filter(GT != "./." & GT != "0/0" & GT != "1/1" ) %>%
    filter(READ_DEPTH == read.depth.threshold) %>%
    filter(ALLELE_REF_DEPTH < allele.depth.threshold | ALLELE_ALT_DEPTH < allele.depth.threshold) %>%
    filter(ALLELE_COVERAGE_RATIO < -allele.imbalance.threshold | ALLELE_COVERAGE_RATIO > allele.imbalance.threshold) %>%
    select(LOCUS, POS, POP_ID, INDIVIDUALS) %>%
    arrange(LOCUS, POS, POP_ID, INDIVIDUALS)
  
  
  
  write.table(blacklist,
              filename,
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F
  )
  
  # interesting stats.
  erased.genotype.number <- length(blacklist$INDIVIDUALS)
  total.genotype.number <- length(tidy.vcf.file$GT[tidy.vcf.file$GT != "./."])
  percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 2), "%", sep = " ")
  
  
  invisible(cat(sprintf(
    "Blacklist erase genotypes:
1. Read depth threshold: %s
2. REF and ALT coverage threshold: %s
3. Allele imbalance ratio threshold: %s
4. %s of genotypes needs erasing with the function erase_genotypes\n
Filename:
%s
Written in the directory:
%s",
    read.depth.threshold,
    allele.depth.threshold,
    allele.imbalance.threshold,
    percent,
    filename, getwd()
  )))
  
  return(blacklist)
}



#' @title Erase genotypes
#' @description This function erase the genotypes of individuals 
#' based on coverage and genotype likelihood thresholds.
#' @param data The 'batch_x.haplotypes.tsv' or a tidy vcf file.
#' @param is.tidy.vcf Using a tidy VCF file: TRUE or FALSE.
#' @param blacklist.genotypes A blacklist of loci and genotypes 
#' containing at least 2 columns header 'LOCUS' and 'INDIVIDUALS'. The file is 
#' in the global environment (myfile) or in the working directory ("myfile.tsv").
#' @param filename The filename saved to the working directory.
#' @details Genotypes below average quality i.e. below threshold for the
#' coverage of REF and/or ALT allele and genotype likelihood
#' are zeroed from the file. The function erase SNP in the VCF file and loci
#' in the haplotypes file.
#' @rdname erase_genotypes
#' @export
#' @import dplyr
#' @import readr

erase_genotypes <- function(data, is.tidy.vcf, blacklist.genotypes, filename) {
  
  GL <- NULL
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  GROUP <- NULL
  VALUE <- NULL
  # ALLELE_P <- NULL
  ALLELE_COVERAGE_RATIO <- NULL
  READ_DEPTH <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL
  CONSENSUS <- NULL
  CONSENSUS_MAX <- NULL
  GT <- NULL
  ALLELE_Q <- NULL
  
  
  if (stri_detect_fixed(is.tidy.vcf, "F")) {
    file.type <- "haplotypes"
    
    
    if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T) %>%
        rename(LOCUS =`Catalog ID`)
      message("Using the haplotypes file in your directory")
    } else {
      data <- data
      message("Using the haplotypes file from your global environment")
    }
    
    
    if (is.vector(blacklist.genotypes) == "TRUE") {
      blacklist <- read_tsv(blacklist.genotypes, col_names = T) %>%
        mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
        select(LOCUS, INDIVIDUALS) %>%
        distinct(LOCUS, INDIVIDUALS)
      message("Using the blacklist file in your directory")
    } else {
      blacklist <- blacklist.genotypes %>%
        mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
        select(LOCUS, INDIVIDUALS) %>%
        distinct(LOCUS, INDIVIDUALS)
      message("Using the blacklist from your global environment")
    }
    
    message("Erasing... Erasing...")
    
    # haplotypes file preparation
    haplo.prep <- data %>%
      gather(INDIVIDUALS, HAPLOTYPES, -c(LOCUS, Cnt)) %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS))
    
    
    # consensus
    consensus.pop <- haplo.prep %>%
      mutate(CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus")) %>%
      group_by(LOCUS) %>%
      summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
      filter(CONSENSUS_MAX > 0) %>%
      select(LOCUS)
    
    
    # interesting stats.
    erased.genotype.number <- length(blacklist$INDIVIDUALS)
    
    haplo.number <- haplo.prep %>%
      filter(HAPLOTYPES != "-") %>%
      select(HAPLOTYPES)
    
    total.genotype.number <- length(haplo.number$HAPLOTYPES)
    percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 2), "%", sep = " ")
    
    
    # Erasing genotype with the blacklist
    erase <- blacklist %>% 
      select(LOCUS, INDIVIDUALS) %>%
      left_join(haplo.prep, by = c("LOCUS", "INDIVIDUALS")) %>%
      mutate(HAPLOTYPES = rep("-", n()))
    
    keep <- haplo.prep %>%
      anti_join(consensus.pop, by = "LOCUS") %>%
      anti_join(
        blacklist %>%
          select(LOCUS, INDIVIDUALS),
        by = c("LOCUS", "INDIVIDUALS")
      )
    
    new.file <- bind_rows(erase, keep) %>%
      arrange(LOCUS, INDIVIDUALS) %>%
      rename(`Catalog ID` = LOCUS) %>%
      spread(INDIVIDUALS, HAPLOTYPES)
    
    
    
  } else {
    file.type <- "tidy vcf"
    
    
    if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      message("Using the tidy vcf file in your directory")
      
    } else {
      data <- data
      message("Using the tidy vcf file from your global environment")
      
    }
    
    if (is.vector(blacklist.genotypes) == "TRUE") {
      blacklist <- read_tsv(blacklist.genotypes, col_names = T)
      message("Using the blacklist file in your directory")
      
    } else {
      blacklist <- blacklist.genotypes
      message("Using the blacklist from your global environment")
    }
    
    message("Erasing... Erasing...")
    
    # interesting stats.
    erased.genotype.number <- length(blacklist$INDIVIDUALS)
    total.genotype.number <- length(data$GT[data$GT != "./."])
    percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 2), "%", sep = " ")
    
    message("Inspecting tidy VCF for coverage problems...")
    
    
    new.file <- data %>%
      mutate(
        GT = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "./.", GT),
        READ_DEPTH = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", READ_DEPTH)
      )
    message("Step ALLELE REF and ALT...")
    
    new.file <- new.file %>% 
      mutate(
        ALLELE_REF_DEPTH = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_REF_DEPTH),
        ALLELE_ALT_DEPTH = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_ALT_DEPTH)
      )
    message("Step GL ...")
    
    new.file <- new.file %>% 
      mutate(
        ALLELE_COVERAGE_RATIO = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_COVERAGE_RATIO),
        GL = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", GL)
      )
    
    message("Step converting to numeric columns...")
    
    new.file <- new.file %>% 
      mutate(
        READ_DEPTH = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", READ_DEPTH),
        ALLELE_REF_DEPTH = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_REF_DEPTH),
        ALLELE_ALT_DEPTH = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_ALT_DEPTH),
        ALLELE_COVERAGE_RATIO = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_COVERAGE_RATIO),
        GL = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", GL)
      ) %>% 
      mutate(
        READ_DEPTH = suppressWarnings(as.numeric(READ_DEPTH)),
        ALLELE_REF_DEPTH = suppressWarnings(as.numeric(ALLELE_REF_DEPTH)),
        ALLELE_ALT_DEPTH = suppressWarnings(as.numeric(ALLELE_ALT_DEPTH)),
        ALLELE_COVERAGE_RATIO = suppressWarnings(as.numeric(ALLELE_COVERAGE_RATIO)),
        GL = suppressWarnings(as.numeric(GL))
      )
    #       mutate(
    #         ALLELE_P = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_P),
    #         ALLELE_Q = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_Q)
    #         )
  }
  
  
  message("Saving the file in your working directory...")
  
  #   write.table(new.file, filename, sep = "\t", row.names = F, col.names = T, quote = F)
  write_tsv(new.file, filename, append = FALSE, col_names = TRUE)
  
  
  invisible(cat(sprintf(
    "Erasing genotypes of individuals in the %s file.
Out of a total of %s genotypes 
Erased genotypes: %s (= %s )\n
Filename:
%s
Written in the directory:
%s",
    file.type, total.genotype.number, erased.genotype.number, percent, filename, getwd()
  )))
  return(new.file)
}



