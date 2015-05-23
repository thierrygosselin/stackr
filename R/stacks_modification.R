#### STACKS SUMSTATS MODIFICATIONS  ############################################
#' @title STACKS sumstats.tsv file modification
#' @description Import STACKS sumstats.tsv and add useful info (frequency of REF allele, Global MAF) for subsequent use. to a tidy format.
#' @param sumstats The STACKS batch_x.susmtats.tsv file.
#' @param skip.line The number of line without the header to start reading the data
#' @param pop.num The number of populations.
#' @param pop.col.types Integer or Character used in STACKS populations module?
#' @param pop.integer.equi When Integer was used for your population id, give the character equivalence
#' @param pop.levels A character string with your populations in order.
sumstats_prep <- function(sumstats, skip.line, pop.num, pop.col.types, pop.integer.equi, pop.levels) {

  if(pop.col.types == "integer"){

    sumstats.prep <- read_delim(sumstats,
                           delim = "\t",
                           na = "NA",
                          skip = skip.line,
                           progress = interactive(),
                           col_names = c("BATCH","LOCUS","CHROM","POS","COL","POP_ID","ALLELE_P","ALLELE_Q","N","FREQ_ALLELE_P","HET_O","HOM_O","HET_E","HOM_E","PI","SMOOTHED_PI","SMOOTHED_PI_P_VALUE","FIS","SMOOTHED_FIS","SMOOTHED_FIS_P_VALUE","PRIVATE"),
                          col_types = "iiciiiccddddddidddddi"
                          ) %>%
      mutate(
        FREQ_ALLELE_Q = 1 - FREQ_ALLELE_P,
        POP_ID = stri_replace_all_fixed(POP_ID, seq(from = 1, to = pop.num, by = 1), pop.integer.equi, vectorize_all=F),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
        )

  } else if (pop.col.types == "character") {

    sumstats.prep <- read_delim(sumstats,
                           delim = "\t",
                           na = "NA",
                            skip = skip.line,
                            progress = interactive(),
                            col_names = c("BATCH","LOCUS","CHROM","POS","COL","POP_ID","ALLELE_P","ALLELE_Q","N","FREQ_ALLELE_P","HET_O","HOM_O","HET_E","HOM_E","PI","SMOOTHED_PI","SMOOTHED_PI_P_VALUE","FIS","SMOOTHED_FIS","SMOOTHED_FIS_P_VALUE","PRIVATE"),
                            col_types = "iiciicccddddddidddddi"
                            ) %>%
      mutate(
        FREQ_ALLELE_Q = 1 - FREQ_ALLELE_P,
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
        )

  } else {

    col.types = NULL
    
    sumstats.prep <- read_delim(sumstats,
                           delim = "\t",
                          na = "NA",
                           skip = skip.line,
                          progress = interactive(),
                          col_names = c("BATCH","LOCUS","CHROM","POS","COL","POP_ID","ALLELE_P","ALLELE_Q","N","FREQ_ALLELE_P","HET_O","HOM_O","HET_E","HOM_E","PI","SMOOTHED_PI","SMOOTHED_PI_P_VALUE","FIS","SMOOTHED_FIS","SMOOTHED_FIS_P_VALUE","PRIVATE")
                           ) %>%
      mutate(
        FREQ_ALLELE_Q = 1 - FREQ_ALLELE_P,
        POP_ID = stri_replace_all_fixed(POP_ID, seq(from = 1, to = pop.num, by = 1), pop.integer.equi, vectorize_all=F),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
        )
  }
    
  # Global MAF and duplicating columns
  sumstats <- sumstats.prep %>%
    melt(
      id.vars = c("LOCUS", "POS", "POP_ID", "N", "FREQ_ALLELE_P", "FREQ_ALLELE_Q"),
      measure.vars = c("ALLELE_P", "ALLELE_Q"), 
      variable.name = "ALLELE", 
      value.name = "NUCLEOTIDE"
      ) %>%
    filter(NUCLEOTIDE != "-") %>%
    mutate(
      N_IND_POP = ifelse(ALLELE == "ALLELE_P", floor(N * FREQ_ALLELE_P), ceiling(N * FREQ_ALLELE_Q))
      ) %>%
    group_by(LOCUS, POS, ALLELE, NUCLEOTIDE) %>%
    summarise(
      N_SNP = sum(N_IND_POP)
      ) %>%
    group_by(LOCUS, POS) %>%
    summarise(
      N = sum(N_SNP),
      GLOBAL_MAF = N_SNP[ALLELE == "ALLELE_Q"]/N
      ) %>%
    select(-N) %>%
    full_join(sumstats.prep, by = c("LOCUS", "POS")) %>%
    mutate(
      REF = ALLELE_P,
      ALT = ALLELE_Q,
      FREQ_REF = FREQ_ALLELE_P,
      FREQ_ALT = FREQ_ALLELE_Q,
      MAF = FREQ_ALT
      )

  #Change the order of the columns
  sumstats <- sumstats[c("BATCH","LOCUS","CHROM","POS","COL","POP_ID","N","ALLELE_P","ALLELE_Q","FREQ_ALLELE_P","FREQ_ALLELE_Q","HET_O","HOM_O","HET_E","HOM_E","PI","SMOOTHED_PI","SMOOTHED_PI_P_VALUE","FIS","SMOOTHED_FIS","SMOOTHED_FIS_P_VALUE","PRIVATE", "REF", "ALT", "FREQ_REF", "FREQ_ALT", "MAF", "GLOBAL_MAF")]
  
  invisible(cat(sprintf(
  "This is the order of your sampling sites in stacks sumstats file:\n%s\n
  This is the re-order (upstream -> downstream) of your sampling sites now in the stacks sumstats file:\n%s\n
  Inspect Columns classes:
  The column LOCUS is of class %s, and should be integer
  The column POS is of class %s, and should be integer
  The column with the position of the snp on the reads (COL) is of class %s, and should be integer
  The column with the population id (POP_ID) is of class %s, and should be ordered factors
  The column with the frequency of your P allele (FREQ_ALLELE_P) is of class %s, and should be numeric
  The column with the frequency of your Q allele (FREQ_ALLELE_Q) is of class %s, and should be numeric
  The column with the frequency of your REF allele (FREQ_REF) is of class %s, and should be numeric
  The column with the frequency of your ALT allele (FREQ_ALT) is of class %s, and should be numeric
  The column with the frequency of your GLOBAL_MAF is of class %s, and should be numeric
  The column with the Observed Heterosigosity (HET_O)is of class %s, and should be numeric
  The column with the Inbreding coefficient (FIS) is of class %s, and should be numeric
  The column with the number of samples (N) is of class %s, and should be numeric", 
  list(pop.integer.equi),
  list(levels(sumstats$POP_ID)),
  class(sumstats$LOCUS),
  class(sumstats$POS),
  class(sumstats$COL),
  list(class(sumstats$POP_ID)),
  class(sumstats$FREQ_ALLELE_P),
  class(sumstats$FREQ_ALLELE_Q),
  class(sumstats$FREQ_REF),
  class(sumstats$FREQ_ALT),
  class(sumstats$GLOBAL_MAF),
  class(sumstats$HET_O),
  class(sumstats$FIS),
  class(sumstats$N)
  )))

sumstats

}


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
