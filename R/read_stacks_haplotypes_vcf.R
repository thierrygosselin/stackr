#' @title Make a tidy format of the batch_x.haplotypes.vcf file
#' @description Import and transform in tidy format the batch_x.haplotypes.vcf
#' file produced by STACKS.
#' @param haplotypes.vcf.file The haplotypes VCF file.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An optional character string with your populations ordered.
#' @param filter An optional filter of loci can be applied
#' to the haplotypes VCF, using a file in the working directory
#'  (e.g. "myfile.txt") or an object in the global environment.
#'  Highly recommended to reduce the resulting file size.
#' @param filename The name of the file written in the directory.
#' @rdname read_stacks_haplotypes_vcf
#' @export 
#' @import dplyr
#' @import stringi

read_stacks_haplotypes_vcf <- function(haplotypes.vcf.file, pop.id.start, pop.id.end, pop.levels, filter, filename) {
  
  QUAL <- NULL
  FILTER <- NULL
  FORMAT <- NULL
  ID <- NULL
  `#CHROM` <- NULL
  INFO <- NULL
  N <- NULL
  AF <- NULL
  INDIVIDUALS <- NULL
  ALT <- NULL
  READ_DEPTH <- NULL
  REF_FREQ <- NULL
  . <- NULL
  
  
  message("Tidying the Haplotypes VCF")
  
  vcf <- read_delim(
    haplotypes.vcf.file, delim = "\t", 
    skip = 7,#n_max = max.read.lines,
    progress = interactive()
  )%>%
    select(-c(QUAL, FILTER, FORMAT)) %>%
    rename(LOCUS = ID, CHROM = `#CHROM`) %>%
    tidyr::separate(INFO, c("N", "AF"), sep = ";", extra = "warn") %>%
    mutate(
      N = as.numeric(stri_replace_all_fixed(N, "NS=", "", vectorize_all=F)),
      AF = stri_replace_all_fixed(AF, "AF=", "", vectorize_all=F)
    )
  
  
  
  message("Gathering individuals in 1 column...")
  
  vcf <- vcf %>%
    tidyr::gather(INDIVIDUALS, FORMAT, -c(CHROM:AF)) %>%
    tidyr::separate(FORMAT, c("GT", "READ_DEPTH"), sep = ":",
             extra = "warn") %>%
    tidyr::separate(AF, c("REF_FREQ", stri_join("ALT_FREQ", seq(1, max(stri_count_fixed(.$ALT, pattern = ","))), sep = "_")), sep = ",", extra = "drop") %>%
    tidyr::separate(ALT, stri_join("ALT", seq(1, max(stri_count_fixed(.$ALT, pattern = ","))), sep = "_"), extra = "drop")
  
 
  message("Fixing columns...")
  
  vcf <- vcf %>%
    mutate(
      READ_DEPTH = ifelse (READ_DEPTH == "0", "NA", READ_DEPTH),
      CHROM = stri_replace_all_fixed(CHROM, "un", "1", vectorize_all=F),
      CHROM = as.integer(CHROM),
      REF_FREQ = as.numeric(REF_FREQ),
      READ_DEPTH = as.numeric(READ_DEPTH),
      POP_ID = factor(str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
                      levels = pop.levels, ordered =T)
    )
  
  
  
  if (missing(filter) == "TRUE") {
    vcf <- vcf
  } else if (is.vector(filter) == "TRUE") {
    vcf <- read_tsv(filter, col_names = T) %>%
      left_join(vcf, by = "LOCUS")
  } else {
    vcf <- filter %>%
      left_join(vcf, by = "LOCUS")
  }
  
  
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
