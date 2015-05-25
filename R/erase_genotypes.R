#' @title Blacklist erase genotypes.
#' @description This function creates a blacklist of loci and individual
#' genotypes useful for the function erase_genotype.
#' @param tidy.vcf.file A data frame object or file (using ".tsv")
#' of class tidy VCF.
#' @param read.depth.threshold Threshold number
#' of maximum read depth.
#' @param allele.depth.threshold Threshold number
#' of min depth of REF or ALT alleles.
#' @param gl.threshold Threshold number of mean genotype likelihood. 
#' @param filename Name of the file written to the working directory.
#' @details 1. Keep the genotypes below the read depth and the 
#' genotype likelihood threshold. 2. Keep the REF and ALT allele below the
#' coverage threshold. 3. Create a blacklist of genotypes to erase.
#' @return The function returns the individuals genotypes, by loci and 
#' individuals, to erase.
#' @rdname blacklist_erase_genotype
#' @export

blacklist_erase_genotype <- function(tidy.vcf.file, read.depth.threshold, allele.depth.threshold, gl.threshold, filename {
  
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, col_names = T)
    message("Using the tidy vcf file in your directory")
  } else {
    data <- tidy.vcf.file
    message("Using the tidy vcf file from your global environment")
  }
  
  
  blacklist <- tidy.vcf.file %>%
    #   filter(GT == "0/0" | GT == "1/1") %>%
    filter(GT != "./.") %>%
    filter(READ_DEPTH < read.depth.threshold) %>%
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
    rename(SAMPLES = INDIVIDUALS) %>%
    arrange(LOCUS, SAMPLES)
  
  write.table(blacklist,
              filename,
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F
  )
  
  invisible(cat(sprintf(
    "Blacklist erase genotypes:
1. Read depth thresold: %
2. REF and ALT coverage threshold: %s
3. Genotype likelihood threshold: %s
Filename:
%s
Written in the directory:
%s",
    read.depth.threshold,
    allele.depth.threshold,
    gl.threshold,
    filename, getwd()
    
  )))
  return(blacklist)
}




#' @title Erase genotypes.
#' @description This function erase the genotype in
#' the 'batch_x.haplotypes.tsv' or a tidy vcf file.
#' @param data The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param is.tidy.vcf Using a tidy VCF file: TRUE or FALSE.
#' @param blacklist.genotypes A blacklist of loci and genotypes 
#' containing at least 2 columns header 'LOCUS' and 'SAMPLES'.
#' @param filename The filename saved to the working directory.
#' @details Genotypes below average quality: below threshold for the
#' read coverage, REF and/or ALT depth coverage and genotype likelihood
#' are zeroed from the file.
#' @rdname erase_genotypes
#' @export

erase_genotypes <- function(data, is.tidy.vcf, blacklist.genotypes, filename) {
  
  if (stri_detect_fixed(is.tidy.vcf, "F")) {
    message("Using the haplotypes file")
    file.type <- "haplotypes"
    
    
    if (is.vector(haplotypes.file) == "TRUE") {
      data <- read_tsv(data, col_names = T) %>%
        rename(LOCUS =`Catalog ID`)
      message("Using the file in your directory")
    } else {
      data <- data
      message("Using the file from your global environment")
    }
    
    
    if (is.vector(blacklist.genotypes) == "TRUE") {
      blacklist.genotypes <- read_tsv(blacklist.genotypes, col_names = T)
      message("Using the blacklist file in your directory")
    } else {
      blacklist.genotypes <- blacklist.genotypes
      message("Using the blacklist from your global environment")
    }
    
    erased.genotype.number <- length(blacklist.genotypes$STATUS)
    
    haplo.number <- data %>%
      melt(id.vars = c("LOCUS", "Cnt"), 
           variable.name = "SAMPLES",
           value.name = "HAPLOTYPES"
      ) %>%
      filter(HAPLOTYPES != "-") %>%
      select(HAPLOTYPES)
    
    total.genotype.number <- length(haplo.number)
    
    # Erasing genotype with the blacklist
    new.file <- data %>%
      melt(id.vars = c("LOCUS", "Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES") %>%
      mutate(
        HAPLOTYPES = ifelse(SAMPLES %in% blacklist.genotypes$SAMPLES & LOCUS %in% blacklist.genotypes$LOCUS, "-", HAPLOTYPES)
      ) %>%
      rename(`Catalog ID` = LOCUS) %>%
      dcast(`Catalog ID`+Cnt~SAMPLES, value.var="HAPLOTYPES")
    
    
    
  } else {
    message("Using the tidy vcf file")
    file.type <- "tidy vcf"
    
    
    if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T) %>%
        rename(LOCUS =`Catalog ID`)
      message("Using the file in your directory")
      
    } else {
      data <- data
      message("Using the file from your global environment")
      
    }
    
    if (is.vector(blacklist.genotypes) == "TRUE") {
      blacklist.genotypes <- read_tsv(blacklist.genotypes, col_names = T)
      message("Using the blacklist file in your directory")
      
    } else {
      blacklist.genotypes <- blacklist.genotypes
      message("Using the blacklist from your global environment")
    }
    
     erased.genotype.number <- length(blacklist.genotypes$STATUS)
      total.genotype.number <- length(vcf.tidy$GT[vcf.tidy$GT != "./."])
      
      new.file <- data %>%
        mutate(
          GT = ifelse(INDIVIDUALS %in% blacklist.genotypes$SAMPLES & LOCUS %in% blacklist.genotypes$LOCUS, "-", INDIVIDUALS)
        )
      
  }
    
    
    write.table(new.file, filename, sep = "\t", row.names = F,
                col.names = T, quote = F)
    
    
    invisible(cat(sprintf(
      "Erasing genotypes of individuals in the %s file.
Erased genotypes: %s.
Out of a total of %s genotypes.
Filename:
%s
Written in the directory:
%s",
      file.type, erased.genotype.number, filename, getwd()
    )))
    return(new.file)
  }
  
