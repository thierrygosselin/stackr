#' @title Blacklist erase genotypes
#' @description This function creates a blacklist of loci and individual
#' genotypes useful for the function erase_genotype.
#' @param tidy.vcf.file A data frame object or file (using ".tsv")
#' of class tidy VCF.
#' @param allele.depth.threshold Threshold number
#' of min depth of REF or ALT alleles.
#' @param gl.threshold Threshold number of mean genotype likelihood. 
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


blacklist_erase_genotype <- function(tidy.vcf.file, allele.depth.threshold, gl.threshold, filename) {
  GT <- NULL
  GL <- NULL
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL
  MIN_REF <- NULL
  MIN_ALT <- NULL
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    tidy.vcf.file <- read_tsv(tidy.vcf.file, col_names = T)
    message("Using the tidy vcf file in your directory")
  } else {
    tidy.vcf.file <- tidy.vcf.file
    message("Using the tidy vcf file from your global environment")
  }
  
  blacklist <- tidy.vcf.file %>%
    filter(GT != "./.") %>%
    filter(GL < gl.threshold) %>%
    group_by(LOCUS, INDIVIDUALS, POP_ID) %>%
    summarise(
      MIN_REF = min(ALLELE_REF_DEPTH, na.rm = T),
      MIN_ALT = min(ALLELE_ALT_DEPTH, na.rm = T)
    ) %>%
    filter(MIN_REF < allele.depth.threshold | MIN_ALT < allele.depth.threshold) %>%
    ungroup() %>%
    mutate(
      STATUS = rep("erased", n()),
      INDIVIDUALS = as.character(INDIVIDUALS)
    ) %>%
    arrange(LOCUS, INDIVIDUALS)
  
  write.table(blacklist,
              filename,
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F
  )
  
  # interesting stats.
  erased.genotype.number <- length(blacklist$STATUS)
  total.genotype.number <- length(tidy.vcf.file$GT[tidy.vcf.file$GT != "./."])
  percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 2), "%", sep = " ")
  
  
  invisible(cat(sprintf(
    "Blacklist erase genotypes:
1. REF and ALT coverage threshold: %s
2. Genotype likelihood threshold: %s
3. %s of genotypes needs erasing with the function erase_genotypes\n
Filename:
%s
Written in the directory:
%s",
    allele.depth.threshold,
    gl.threshold,
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
  STATUS <- NULL
  GROUP <- NULL
  VALUE <- NULL
  ALLELE_P <- NULL
  ALLELE_COVERAGE_RATIO <- NULL
  READ_DEPTH <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL

  
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
      blacklist.genotypes <- read_tsv(blacklist.genotypes, col_names = T)
      message("Using the blacklist file in your directory")
    } else {
      blacklist.genotypes <- blacklist.genotypes
      message("Using the blacklist from your global environment")
    }
    
    message("Erasing... Erasing...")
    
    # haplotypes file preparation
    haplo.prep <- data %>%
      gather(SAMPLES, HAPLOTYPES, -c(LOCUS, Cnt))
#       melt(id.vars = c("LOCUS", "Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES")
    
    # interesting stats.
    erased.genotype.number <- length(blacklist.genotypes$STATUS)
    
    haplo.number <- haplo.prep %>%
      filter(HAPLOTYPES != "-") %>%
      select(HAPLOTYPES)
    
    total.genotype.number <- length(haplo.number$HAPLOTYPES)
    percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 2), "%", sep = " ")
    
    # Erasing genotype with the blacklist
    new.file <- haplo.prep %>%
      mutate(
        HAPLOTYPES = ifelse(SAMPLES %in% blacklist.genotypes$INDIVIDUALS & LOCUS %in% blacklist.genotypes$LOCUS, "-", HAPLOTYPES)
      ) %>%
      rename(`Catalog ID` = LOCUS) %>%
      spread(SAMPLES, HAPLOTYPES)
    
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
      blacklist.genotypes <- read_tsv(blacklist.genotypes, col_names = T)
      message("Using the blacklist file in your directory")
      
    } else {
      blacklist.genotypes <- blacklist.genotypes
      message("Using the blacklist from your global environment")
    }

    message("Erasing... Erasing...")
    
    # interesting stats.
    erased.genotype.number <- length(blacklist.genotypes$STATUS)
    total.genotype.number <- length(data$GT[data$GT != "./."])
    percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 2), "%", sep = " ")
    
    order.vcf <- names(data)
    
    erase <- blacklist.genotypes %>% 
      select(LOCUS, INDIVIDUALS, STATUS) %>%
      left_join(data, by = c("LOCUS", "INDIVIDUALS")) %>%
      mutate(GT = rep("./.", n())) %>%
      gather(GROUP, VALUE, ALLELE_P:ALLELE_COVERAGE_RATIO) %>%
      mutate(VALUE = rep("NA", n())) %>%
      spread(GROUP, VALUE) %>%
      mutate(
        READ_DEPTH = as.numeric(READ_DEPTH),
        ALLELE_REF_DEPTH = as.numeric(ALLELE_REF_DEPTH),
        ALLELE_ALT_DEPTH = as.numeric(ALLELE_ALT_DEPTH),
        ALLELE_COVERAGE_RATIO = as.numeric(ALLELE_COVERAGE_RATIO),
        GL = as.numeric(GL)
        )
      
      erase <- erase[order.vcf]

    keep <- data %>% 
      anti_join(
        blacklist.genotypes %>%
          select(LOCUS, INDIVIDUALS, STATUS),
        by = c("LOCUS", "INDIVIDUALS")
      )
    
    new.file <- bind_rows(erase, keep) %>%
      arrange(LOCUS, POS, POP_ID, INDIVIDUALS)
    
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
