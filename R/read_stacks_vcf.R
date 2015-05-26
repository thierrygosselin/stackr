#' @title Read a VCF file produced by STACKS and transform in tidy format
#' @description Import a VCF file created by STACKS and mofify to a tidy format.
#' @param vcf.file The VCF file created by STACKS.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An optional character string with your populations ordered.
#' @param filter An optional filter of loci can be applied to the vcf, using
#' a file in the working directory (e.g. "myfile.txt") or an object
#' in the global environment.
#' @param filename The name of the file written in the directory.
#' @rdname read_stacks_vcf
#' @export
#' @import plyr
#' @import dplyr
#' @import readr
 
# read_stacks_vcf <- function(vcf.file, skip.line, max.read.lines, pop.id.start, pop.id.end, pop.levels, filter, filename) {

read_stacks_vcf <- function(vcf.file, pop.id.start, pop.id.end, pop.levels, filter, filename) {

  message("Tidying the VCF...")
  
  vcf <- read_delim(
    vcf.file, delim = "\t", 
    skip = 9,#n_max = max.read.lines,
    progress = interactive()
    )
  
  vcf <- vcf %>%
    select(-c(QUAL, FILTER, FORMAT)) %>%
    rename(LOCUS = ID, CHROM = `#CHROM`) %>%
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
    )
  
  vcf <- vcf[c("CHROM", "LOCUS", "POS", "N", "REF", "ALT", "REF_FREQ", "ALT_FREQ", "POP_ID", "INDIVIDUALS", "GT", "ALLELE_P", "ALLELE_Q", "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "ALLELE_COVERAGE_RATIO", "GL")]
  
  if (missing(filter) == "TRUE") {
    vcf <- vcf
  } else if (is.vector(filter) == "TRUE") {
    vcf <- read_tsv(filter, col_names = T) %>%
      left_join(vcf, by = "LOCUS")
  } else {
    vcf <- filter %>%
      left_join(vcf, by = "LOCUS")
  }
  
#   message("Saving the file in your working directory...")
#   write_tsv(vcf, filename, append = FALSE, col_names = TRUE)
# #   write.table(vcf, filename, sep = "\t", row.names = F, col.names = T,
# #               quote = F)

 
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


