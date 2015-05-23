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

read_stacks_vcf_old <- function(vcf.file, skip.line, max.read.lines, pop.id.start, pop.id.end, pop.levels) {
  read_delim(
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
    READ_DEPTH = as.numeric(stri_replace_all_fixed(READ_DEPTH, "0", "NA",
                                                   vectorize_all = F))
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
    ALLELE_REF_DEPTH = as.numeric(stri_replace_all_fixed(ALLELE_REF_DEPTH,
                                                         ".", "NA",
                                                         vectorize_all=F)),
    ALLELE_ALT_DEPTH = as.numeric(stri_replace_all_fixed(ALLELE_ALT_DEPTH,
                                                         ".", "NA",
                                                         vectorize_all=F)),
    ALLELE_REF_DEPTH = ifelse(GT == "./.", "NA",
                       ifelse(GT == "0/0", READ_DEPTH,
                       ifelse(GT == "1/1", "NA", ALLELE_REF_DEPTH))),
    ALLELE_ALT_DEPTH = ifelse(GT == "./.", "NA",
                       ifelse(GT == "1/1", READ_DEPTH,
                       ifelse(GT == "0/0", "NA", ALLELE_ALT_DEPTH))),
    POP_ID = factor(str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
                    levels = pop.levels, ordered =T)
    )
}

