#' @title Read a VCF file produced by STACKS and transform in tidy format
#' @description Import a VCF file created by STACKS and mofify to a tidy format.
#' @param vcf.file The VCF file created by STACKS.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An optional character string with your populations ordered.
#' @param whitelist An optional filter of loci can be applied to the vcf, using
#' a file in the working directory (e.g. "myfile.txt") or an object
#' in the global environment.
#' @param blacklist (optional) Blacklist file with loci to filter out
#' of the vcf, using a file in the working directory (e.g. "myfile.txt")
#' or an object in the global environment.
#' @param filename (optional) The name of the file written in the directory.
#' @rdname read_stacks_vcf
#' @export
#' @import dplyr
#' @import readr


read_stacks_vcf <- function(vcf.file, pop.id.start, pop.id.end, pop.levels, whitelist, blacklist, filename) {
  
  X1 <- NULL
  POP_ID <- NULL
  QUAL <- NULL
  FILTER <- NULL
  FORMAT <- NULL
  ID <- NULL
  `#CHROM` <- NULL
  INFO <- NULL
  N <- NULL
  AF <- NULL
  INDIVIDUALS <- NULL
  REF <- NULL
  ALT <- NULL
  READ_DEPTH <- NULL
  REF_FREQ <- NULL
  ALT_FREQ <- NULL
  ALLELE_DEPTH <- NULL
  GT <- NULL
  GL <- NULL
  ALLELE_P <- NULL
  ALLELE_Q <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL
  
  
  message("Importing the VCF...")
  
  vcf <- read_delim(
    vcf.file, delim = "\t", 
    skip = 9,#n_max = max.read.lines,
    progress = interactive()
  ) %>% 
    select(-c(QUAL, FILTER, FORMAT)) %>%
    rename(LOCUS = ID, CHROM = `#CHROM`)
  
  if (missing(whitelist) == "TRUE") {
    vcf <- vcf
    message("No whitelist to apply to the VCF...")
    
  } else if (is.vector(whitelist) == "TRUE") {
    vcf <- read_tsv(whitelist, col_names = T) %>%
      left_join(vcf, by = "LOCUS")
    message("Filtering the VCF with the whitelist from your directory")
    
  } else {
    vcf <- whitelist %>%
      left_join(vcf, by = "LOCUS")
    message("Filtering the VCF with the whitelist from your global environment")
    
  }
  
  if (missing(blacklist) == "TRUE") {
    vcf <- vcf
    message("No blacklist to apply to the VCF...")
    
  } else if (is.vector(blacklist) == "TRUE") {
    message("Filtering the VCF with the blacklist from your directory")
    
    vcf <- vcf  %>% 
      anti_join(
        read_tsv(blacklist, col_names = T),
        by = "LOCUS"
      )
    
  } else {
    message("Filtering the VCF with the blacklist from your global environment")
    
    vcf <- vcf %>% 
      anti_join(blacklist, by = "LOCUS")
    
  }
  
  
  message("Tidying the VCF...")
  
  vcf <- vcf %>%
    separate(INFO, c("N", "AF"), sep = ";", extra = "error") %>%
    mutate(
      N = as.numeric(stri_replace_all_fixed(N, "NS=", "", vectorize_all=F)),
      AF = stri_replace_all_fixed(AF, "AF=", "", vectorize_all=F)
    ) %>%
    separate(AF, c("REF_FREQ", "ALT_FREQ"), sep = ",", extra = "error") %>%
    mutate(
      REF_FREQ = as.numeric(REF_FREQ),
      ALT_FREQ = as.numeric(ALT_FREQ)
    )
  
  message("Gathering individuals in 1 column...")
  
  vcf <- vcf %>%
    gather(INDIVIDUALS, FORMAT, -c(CHROM:ALT_FREQ)) %>%
    separate(FORMAT, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"), sep = ":",
             extra = "error") %>%
    mutate(
      READ_DEPTH = as.numeric(READ_DEPTH),
      READ_DEPTH = as.numeric(ifelse (READ_DEPTH == "0", "NA", READ_DEPTH))
    ) %>%
    separate(ALLELE_DEPTH, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"), sep = ",",
             extra = "error") %>%
    separate(GT, c("ALLELE_P", "ALLELE_Q"), sep = "/", extra = "error",
             remove = F)
  
  message("Fixing columns...")
  
  vcf <- vcf %>%
    mutate(
      CHROM = as.numeric(stri_replace_all_fixed(CHROM, "un", "1", 
                                                vectorize_all=F)),
      GL = as.numeric(stri_replace_all_fixed(GL, c(".,", ",."), "",
                                             vectorize_all=F)),
      ALLELE_P = ifelse (ALLELE_P == ".", "NA",
                         ifelse(ALLELE_P == "0", REF, ALT)),
      ALLELE_Q = ifelse (ALLELE_Q == ".", "NA",
                         ifelse(ALLELE_Q == "0", REF, ALT)),
      ALLELE_REF_DEPTH = as.numeric(ALLELE_REF_DEPTH),
      ALLELE_ALT_DEPTH = as.numeric(ALLELE_ALT_DEPTH),
      ALLELE_REF_DEPTH = as.numeric(ifelse (ALLELE_REF_DEPTH == "0", "NA", ALLELE_REF_DEPTH)),
      ALLELE_ALT_DEPTH = as.numeric(ifelse (ALLELE_ALT_DEPTH == "0", "NA", ALLELE_ALT_DEPTH)),
      ALLELE_COVERAGE_RATIO = as.numeric(ifelse(GT == "./." | GT == "0/0" | GT == "1/1", "NA",
                                                ((ALLELE_ALT_DEPTH - ALLELE_REF_DEPTH)/(ALLELE_ALT_DEPTH + ALLELE_REF_DEPTH)))),
      POP_ID = factor(str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
                      levels = pop.levels, ordered =T)
    ) %>%
    arrange(LOCUS, POS, POP_ID, INDIVIDUALS)
  
  vcf <- vcf[c("CHROM", "LOCUS", "POS", "N", "REF", "ALT", "REF_FREQ", "ALT_FREQ", "POP_ID", "INDIVIDUALS", "GT", "ALLELE_P", "ALLELE_Q", "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "ALLELE_COVERAGE_RATIO", "GL")]
  
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(vcf, filename, append = FALSE, col_names = TRUE)
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
  vcf
}


