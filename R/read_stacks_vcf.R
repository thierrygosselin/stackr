#' @title read_stacks_vcf.
#' @description Import a VCF file created by STACKS and mofify to a tidy format.
#' @param vcf.file The VCF file created by STACKS.
#' @param skip.line Default is 9 with vcf created in STACKS.
#' @param max.read.lines The number of markers + the header line.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An optional character string with your populations ordered.
#' @param filter To apply a filter to the vcf, use a file or object.
#' @param filename The name of the file written in the directory.

read_stacks_vcf <- function(vcf.file, skip.line, max.read.lines, pop.id.start, pop.id.end, pop.levels, filter, filename) {
  
  vcf <- read_delim(
    vcf.file,
    delim = "\t",
    skip = skip.line,
    n_max = max.read.lines,
    progress = interactive()
    ) %>%
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
    ) %>%
  melt(
    id.vars = c("CHROM", "LOCUS", "POS", "N", "REF", "ALT", "REF_FREQ",
                "ALT_FREQ"),
    variable.name = "INDIVIDUALS", 
    value.name = "FORMAT"
    ) %>%
  separate(FORMAT, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"), sep = ":",
           extra = "error") %>%
  mutate(
    READ_DEPTH = as.numeric(READ_DEPTH),
    READ_DEPTH = as.numeric(ifelse (READ_DEPTH == "0", "NA", READ_DEPTH))
  ) %>%
  separate(ALLELE_DEPTH, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"), sep = ",",
           extra = "error") %>%
  separate(GT, c("ALLELE_P", "ALLELE_Q"), sep = "/", extra = "error",
           remove = F) %>%
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

  if (missing(filter) == "TRUE") {
    vcf <- vcf
  } else if (is.vector(filter) == "TRUE") {
    vcf <- read_tsv(filter, col_names = T) %>%
      left_join(vcf, by = "LOCUS")
  } else {
    vcf <- filter %>%
      left_join(vcf, by = "LOCUS")
  }
  
  write.table(vcf, filename, sep = "\t", row.names = F, col.names = T,
              quote = F)

invisible(cat(sprintf(
"Stacks VCF filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
vcf
}

#' @title tidy_vcf_summary.
#' @description Summarise and prepare the tidy VCF. Summary, by population and markers (SNP), of frequency of the REF and the ALT alleles, the observed and the expected heterozygosity and the inbreeding coefficient. The Global MAF of Loci, with STACKS GBS/RAD loci = read or de novo haplotypes, is included and repeated over SNP.
#' @param data The tidy VCF file created with read_stacks_vcf.

tidy_vcf_summary <- function(data) {

  vcf.summary <- data %>%
    filter(GT != "./.") %>%
    group_by(LOCUS, POS, POP_ID) %>%
    summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT[GT == "0/0"])),
      PQ = as.numeric(length(GT[GT == "0/1" | GT == "1/0"])),
      QQ = as.numeric(length(GT[GT == "1/1"]))
      ) %>%
    mutate(
      FREQ_REF = ((PP*2) + PQ)/(2*N),
      FREQ_ALT = ((QQ*2) + PQ)/(2*N),
      HET_O = PQ/N,
      HET_E = 2 * FREQ_REF * FREQ_ALT,
      FIS = ifelse(HET_O == 0, 0, round (((HET_E - HET_O) / HET_E), 6))
      )
  
  global.maf <- vcf.summary %>%
    group_by(LOCUS, POS) %>%
    summarise_each_(funs(sum), vars = c("N", "PP", "PQ", "QQ")) %>%
    mutate(GLOBAL_MAF = (PQ + (2 * QQ)) / (2*N)) %>%
    select(LOCUS, POS, GLOBAL_MAF)
  
  vcf.prep <- global.maf %>%
    left_join(vcf.summary, by = c("LOCUS", "POS"))
  
  vcf.prep <- vcf.prep[c("LOCUS", "POS", "POP_ID", "N", "PP", "PQ", "QQ", "FREQ_REF", "FREQ_ALT", "GLOBAL_MAF", "HET_O", "HET_E", "FIS")]
  
  return(vcf.prep)
}

