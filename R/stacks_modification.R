#### STACKS SUMSTATS MODIFICATIONS  ############################################
#' @title Import and prepare for analysis a batch_x.sumstats.tsv file produced by STACKS
#' @description Import STACKS sumstats.tsv and add useful info (frequency of REF allele, Global MAF) for subsequent use. to a tidy format.
#' @param sumstats The STACKS batch_x.susmtats.tsv file.
#' @param skip.line The number of line without the header to start reading the data
#' @param pop.num The number of populations.
#' @param pop.col.types Integer or Character used in STACKS populations module?
#' @param pop.integer.equi When Integer was used for your population id, give the character equivalence
#' @param pop.levels A character string with your populations in order.
#' @param filename The name of the file written in the directory.
#' @export
#' @rdname sumstats_prep
#' @import dplyr
#' @import readr

sumstats_prep <- function(sumstats, skip.line, pop.num, pop.col.types, pop.integer.equi, pop.levels, filename) {
  
  BATCH <- NULL
  LOCUS <- NULL
  CHROM <- NULL
  POS <- NULL
  COL <- NULL
  POP_ID <- NULL
  ALLELE_P <- NULL
  ALLELE_Q <- NULL
  N <- NULL
  FREQ_ALLELE_P <- NULL
  HET_O <- NULL
  HOM_O <- NULL
  HET_E <- NULL
  HOM_E <- NULL
  PI <- NULL
  SMOOTHED_PI <- NULL
  SMOOTHED_PI_P_VALUE <- NULL
  FIS <- NULL
  SMOOTHED_FIS <- NULL
  SMOOTHED_FIS_P_VALUE <- NULL
  PRIVATE <- NULL
  FREQ_ALLELE_P <- NULL
  FREQ_ALLELE_Q <- NULL
  ALLELE <- NULL
  NUCLEOTIDE <- NULL
  N_IND_POP <- NULL
  N_SNP <- NULL
  FREQ_ALT <- NULL
  N_P <- NULL
  NTOT <- NULL
  
  
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
    select(c(LOCUS, POS, POP_ID, N, FREQ_ALLELE_P)) %>%
    group_by(LOCUS, POS, POP_ID) %>%
    mutate(
      N_P = 2 * N * FREQ_ALLELE_P # by locus, snp and pop
    ) %>%
    group_by(LOCUS, POS) %>%
    summarise(
      N_P = sum(N_P), # by locus and snp only
      NTOT = sum(N)
    ) %>%
    mutate(
      GLOBAL_MAF = ((2*NTOT)-N_P) / (2*NTOT)
    ) %>%
    select(-N_P, -NTOT) %>%
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
  
  write.table(sumstats, filename, sep = "\t", row.names = F, col.names = T,
              quote = F)
  
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
  The column with the number of samples (N) is of class %s, and should be numeric\n
  Filename:
  %s
  Written in the directory:
  %s",
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
    class(sumstats$N),
    filename, getwd()
  )))
  
  sumstats
  
}
