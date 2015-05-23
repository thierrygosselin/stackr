#### STACKS SUMSTATS MODIFICATIONS  ############################################
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

# Deprecated, will be discontinued in future version
#### ALLELIC INVERSION SUMSTATS ###################################
# Add a column with the number of individual for each nucleotide
# Summarise by SNP and nucleotide the number of individual
# Remove nucleotide with no individual (in the cases where allele were fixed)
# Define the REF allele
# Remove duplicate (same number of allele...)
# Merge with the total number of individual per SNP
# Calculate the global MAF
# Keep the interesting column and change the name of NUCLEOTIDE to REF
# Switch the frequencies if inversion was found
# FREQ_REF = Frequency of the REFERENCE ALLELE
# FREQ_ALT = Frequency of the ALTERNATE ALLELE = LOCAL MINOR ALLELE FREQUENCY = LOCAL MAF
# the warning message: attributes are not identical across measure variables; they will be dropped is normal

sumstats_inversion <- function (sumstats.modified) {

  inversion <- sumstats %>%
    melt(
    id.vars = c("POS", "POP_ID", "N", "FREQ_ALLELE_P", "FREQ_ALLELE_Q"),
    measure.vars = c("ALLELE_P", "ALLELE_Q"), 
    variable.name = "ALLELE", 
    value.name = "NUCLEOTIDE"
    ) %>%
  mutate(
    NUMBER_INDIVIDUAL = ifelse(ALLELE == "ALLELE_P", (N * FREQ_ALLELE_P), (N * FREQ_ALLELE_Q))
    ) %>%
  group_by(POS, NUCLEOTIDE) %>%
  summarise(
    NUMBER_INDIVIDUAL_TOTAL = sum(NUMBER_INDIVIDUAL)
    ) %>%
  filter(subset =! (NUMBER_INDIVIDUAL_TOTAL == 0)) %>%
  filter(NUMBER_INDIVIDUAL_TOTAL == max(NUMBER_INDIVIDUAL_TOTAL)) %>%
  distinct(POS) %>%
  merge(
    sumstats %>%
      group_by(POS) %>%
      summarise(TOTAL_INDIVIDUAL=sum(N)),
    by="POS"
    ) %>%
  mutate(
    GLOBAL_MAF = (1-(NUMBER_INDIVIDUAL_TOTAL/TOTAL_INDIVIDUAL))
    ) %>%
  select(POS,REF=NUCLEOTIDE,GLOBAL_MAF) %>%
  merge(sumstats, by="POS") %>%
  mutate(
    INVERSION = ifelse(as.character(ALLELE_P) == as.character(REF), "no", "inversion"),
    FREQ_REF = ifelse(INVERSION=="inversion", FREQ_ALLELE_Q, FREQ_ALLELE_P),
    FREQ_ALT = ifelse(INVERSION=="inversion", FREQ_ALLELE_P, FREQ_ALLELE_Q)
    )

  # reorder sumstats.inversion with de novo you don't need all those column
  # inversion <- inversion[c("LOCUS","POS","COL","POP_ID","N", "ALLELE_P", "REF", "INVERSION","ALLELE_Q","FREQ_ALLELE_P", "FREQ_REF", "FREQ_ALLELE_Q", "FREQ_ALT", "GLOBAL_MAF", "HET_O","HOM_O","HET_E","HOM_E","PI","FIS", "PRIVATE")]

  # With a reference genome:
  inversion <- inversion[c("BATCH","LOCUS","CHROM","POS","COL","POP_ID","N", "ALLELE_P", "REF", "INVERSION","ALLELE_Q","FREQ_ALLELE_P", "FREQ_REF", "FREQ_ALLELE_Q", "FREQ_ALT", "GLOBAL_MAF", "HET_O","HOM_O","HET_E","HOM_E","PI","SMOOTHED_PI","SMOOTHED_PI_P_VALUE","FIS","SMOOTHED_FIS","SMOOTHED_FIS_P_VALUE","PRIVATE")]

  #The result might be a bit different than the VCF REF allele, because of the ways you can handle duplicate entries (the 2 alleles at a SNP with the same number of individual )

  # get the number of inversion per pop
  number.inversion <- inversion %>%
    filter (INVERSION == "inversion") %>%
    select (POS, POP_ID)


  number.inversion.pop <- number.inversion %>%
   group_by(POP_ID) %>%
   tally()
  number.inversion.pop

  write.table(number.inversion.pop, "number.inversion.pop.tsv", sep = "\t", row.names = F, 
    col.names = T, quote = F)

  #save the result if you are working a lot on the filter, so skip this step next time 
  write.table(inversion, "sumstats.inversion.tsv", sep = "\t", row.names = F, 
    col.names = T, quote = F)


  invisible(
          cat(
           sprintf(
"A total of %s unique loci were detected and corrected for inversions in the sumstats file
Directory: %s
2 files were written to the directory:
1. number.inversion.pop.tsv 
2. sumstats.inversion.tsv",
n_distinct(number.inversion$POS), getwd()
        )
        )
      )
number.inversion.pop
inversion
}

sumstats_global_maf <- function (sumstats.modified) {

  global.maf <- sumstats %>%
    melt(
    id.vars = c("POS", "POP_ID", "N", "FREQ_ALLELE_P", "FREQ_ALLELE_Q"),
    measure.vars = c("ALLELE_P", "ALLELE_Q"), 
    variable.name = "ALLELE", 
    value.name = "NUCLEOTIDE"
    ) %>%
  mutate(
    NUMBER_INDIVIDUAL = ifelse(ALLELE == "ALLELE_P", (N * FREQ_ALLELE_P), (N * FREQ_ALLELE_Q))
    ) %>%
  group_by(POS, NUCLEOTIDE) %>%
  summarise(
    NUMBER_INDIVIDUAL_TOTAL = sum(NUMBER_INDIVIDUAL)
    ) %>%
  filter(subset =! (NUMBER_INDIVIDUAL_TOTAL == 0)) %>%
  filter(NUMBER_INDIVIDUAL_TOTAL == max(NUMBER_INDIVIDUAL_TOTAL)) %>%
  distinct(POS) %>%
  merge(
    sumstats %>%
      group_by(POS) %>%
      summarise(TOTAL_INDIVIDUAL=sum(N)),
    by="POS"
    ) %>%
  mutate(
    GLOBAL_MAF = (1-(NUMBER_INDIVIDUAL_TOTAL/TOTAL_INDIVIDUAL))
    ) %>%
  select(POS,REF=NUCLEOTIDE,GLOBAL_MAF) %>%
  merge(sumstats, by="POS") %>%
  mutate(
    INVERSION = ifelse(as.character(ALLELE_P) == as.character(REF), "no", "inversion"),
    FREQ_REF = ifelse(INVERSION=="inversion", FREQ_ALLELE_Q, FREQ_ALLELE_P),
    FREQ_ALT = ifelse(INVERSION=="inversion", FREQ_ALLELE_P, FREQ_ALLELE_Q)
    )

  # reorder sumstats.inversion with de novo you don't need all those column
  # inversion <- inversion[c("LOCUS","POS","COL","POP_ID","N", "ALLELE_P", "REF", "INVERSION","ALLELE_Q","FREQ_ALLELE_P", "FREQ_REF", "FREQ_ALLELE_Q", "FREQ_ALT", "GLOBAL_MAF", "HET_O","HOM_O","HET_E","HOM_E","PI","FIS", "PRIVATE")]

  # With a reference genome:
  inversion <- inversion[c("BATCH","LOCUS","CHROM","POS","COL","POP_ID","N", "ALLELE_P", "REF", "INVERSION","ALLELE_Q","FREQ_ALLELE_P", "FREQ_REF", "FREQ_ALLELE_Q", "FREQ_ALT", "GLOBAL_MAF", "HET_O","HOM_O","HET_E","HOM_E","PI","SMOOTHED_PI","SMOOTHED_PI_P_VALUE","FIS","SMOOTHED_FIS","SMOOTHED_FIS_P_VALUE","PRIVATE")]

  #The result might be a bit different than the VCF REF allele, because of the ways you can handle duplicate entries (the 2 alleles at a SNP with the same number of individual )

  # get the number of inversion per pop
  number.inversion <- inversion %>%
    filter (INVERSION == "inversion") %>%
    select (POS, POP_ID)


  number.inversion.pop <- number.inversion %>%
   group_by(POP_ID) %>%
   tally()
  number.inversion.pop

  write.table(number.inversion.pop, "number.inversion.pop.tsv", sep = "\t", row.names = F, 
    col.names = T, quote = F)

  #save the result if you are working a lot on the filter, so skip this step next time 
  write.table(inversion, "sumstats.inversion.tsv", sep = "\t", row.names = F, 
    col.names = T, quote = F)


  invisible(
          cat(
           sprintf(
"A total of %s unique loci were detected and corrected for inversions in the sumstats file
Directory: %s
2 files were written to the directory:
1. number.inversion.pop.tsv 
2. sumstats.inversion.tsv",
n_distinct(number.inversion$POS), getwd()
        )
        )
      )
number.inversion.pop
inversion
}



