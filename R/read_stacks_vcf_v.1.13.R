read_stacks_vcf_v.1.13 <- function(vcf.file, skip.line, max.read.lines, pop.id.start, pop.id.end, pop.levels) {
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
  separate(FORMAT, c("GT", "READ_DEPTH", "GL"), sep = ":",
           extra = "error") %>%
  mutate(
    READ_DEPTH = as.numeric(stri_replace_all_fixed(READ_DEPTH, "0", "NA",
                                                   vectorize_all = F))
    ) %>%
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
    POP_ID = factor(str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
                    levels = pop.levels, ordered =T)
    )
}
