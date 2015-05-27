# Write a genepop file from STACKS haplotype file
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE_1", "ALLELE_2", "GENOTYPE", "NUCLEOTIDES"))



#' @title Use the batch_x.haplotypes.tsv file to write a genpop file
#' @description This function can first filter the haplotypes file with a whitelist of loci
#' and a blacklist of individuals (optional).
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci A whitelist of loci and a column header 'LOCUS'.
#' @param blacklist.id A blacklist with individual ID and
#' a column header 'SAMPLES'.
#' @param genepop.filename The name of the file written to the directory.
#' @param genepop.header The first line of the Genepop file 
#' (ex: genepop of 10 pop of Acipenser fulvescens filtered MAF 0.05).
#' @param pop.levels An optional character string with your populations ordered.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_pad

write_genepop <- function(haplotypes.file, 
                          whitelist.loci = NULL, 
                          blacklist.id = NULL, 
                          genepop.filename, genepop.header, 
                          pop.levels, pop.id.start, pop.id.end) {
  
if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == TRUE) {
  
  # No filter
  haplotype.no.filter <- read_tsv(haplotypes.file, col_names = T) %>% 
    rename(Catalog.ID = `Catalog ID`) %>%
    gather(SAMPLES, HAPLOTYPES, -c(Catalog.ID, Cnt)) %>%
    arrange(Catalog.ID)

  
  # Creates a vector containing the loci name
  loci <- unique(haplotype.no.filter$Catalog.ID)
  data <- haplotype.no.filter
  
} else if (is.null(whitelist.loci) == FALSE & is.null(blacklist.id) == TRUE) {
  
  # just whitelist.loci, NO Blacklist of individual
  haplotype.whitelist.loci <- read_tsv(haplotypes.file, col_names = T) %>%
    rename(Catalog.ID = `Catalog ID`) %>%
    gather(SAMPLES, HAPLOTYPES, -c(Catalog.ID, Cnt)) %>%
    right_join(
      read_tsv(whitelist.loci, col_names = T) %>%
        rename(Catalog.ID = LOCUS),
      by = "Catalog.ID") %>%
    arrange(Catalog.ID)
  
  # Creates a vector containing the loci name
  loci <- unique(haplotype.whitelist.loci$Catalog.ID)
  data <- haplotype.whitelist.loci
  
} else if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == FALSE) {
  
  # NO whitelist, JUST Blacklist of individual
  haplotype.blacklist <- read_tsv(haplotypes.file, col_names = T) %>%
    rename(Catalog.ID = `Catalog ID`) %>%
    gather(SAMPLES, HAPLOTYPES, -c(Catalog.ID, Cnt)) %>%
    anti_join(read_tsv(blacklist.id, col_names = T), by = "SAMPLES") %>%
    arrange(Catalog.ID)
  
  
  # Creates a vector containing the loci name
  loci <- unique(haplotype.blacklist$Catalog.ID)
  data <- haplotype.blacklist
  
} else {
  
  # whitelist.loci + Blacklist of individual
  haplotype.whitelist.blacklist <- read_tsv(haplotypes.file, col_names = T) %>% 
    rename(Catalog.ID = `Catalog ID`) %>%
    gather(SAMPLES, HAPLOTYPES, -c(Catalog.ID, Cnt)) %>%
    right_join(
      read_tsv(whitelist.loci, col_names = T) %>%
        rename(Catalog.ID = LOCUS),
      by = "Catalog.ID") %>% 
    anti_join(read_tsv(blacklist.id, col_names = T),by = "SAMPLES") %>%
    arrange(Catalog.ID)
  
  # Creates a vector containing the loci name
  loci <- unique(haplotype.whitelist.blacklist$Catalog.ID)
  data <- haplotype.whitelist.blacklist
}
  # genepop factory...
  # Work on the haplotype file to create a genepop dataframe
  # 1. Starts with a filtered batch_1.haplotypes.tsv stacks file
  # 2. Replace the missing '-' for 'NA'
  # 3. separate haplotype in allele and drop loci > 2 alleles
  # 4. if only 1 allèle is found, copy it
  # 5. shorcut to transform the dataframe in genepop '001' format
  # 6. merge the 2 alleles
  # 7. pivot the table 
  # 8. add a ',' at the end of the sample id

  genepop <- data %>%
    mutate(
      HAPLOTYPES = stri_replace_all_fixed(HAPLOTYPES, "-", "NA", 
                                          vectorize_all=F)
      ) %>%
    separate(
      col = HAPLOTYPES, into = c("ALLELE_1", "ALLELE_2"), 
      sep = "/", extra = "drop", remove = T
      ) %>%
    mutate(
      ALLELE_1 = stri_replace_all_fixed(ALLELE_1, "NA", "0", vectorize_all=F),
      ALLELE_2 = stri_replace_na(ALLELE_2, replacement = "no_allele"),
      ALLELE_2 = ifelse(ALLELE_2 == "no_allele", ALLELE_1, ALLELE_2)
      ) %>%
    melt(id.vars = c("Catalog.ID","Cnt", "SAMPLES"),
         measure.vars = c("ALLELE_1", "ALLELE_2"), 
         variable.name = "ALLELE", 
         value.name = "NUCLEOTIDES"
         ) %>%
    group_by(Catalog.ID) %>%
    mutate(
      NUCLEOTIDES = factor(NUCLEOTIDES),
      NUCLEOTIDES = as.numeric(NUCLEOTIDES),
      NUCLEOTIDES = NUCLEOTIDES-1, # this removes 1 levels to enable missing values = 0
      NUCLEOTIDES = str_pad(NUCLEOTIDES, 3, side = "left", pad = "0")
      ) %>%
    select(-Cnt) %>%
    dcast(Catalog.ID+SAMPLES~ALLELE, value.var="NUCLEOTIDES") %>%
    unite(GENOTYPE, ALLELE_1:ALLELE_2, sep = "") %>%
    dcast(SAMPLES~Catalog.ID, value.var="GENOTYPE") %>%
    mutate(SAMPLES = paste(SAMPLES, ",", sep = ""))

  # Create a vector with the population ordered by levels
    if (missing(pop.levels) == "TRUE") {
    pop <- substr(genepop$SAMPLES, pop.id.start, pop.id.end)
    } else {
    pop <- factor(substr(genepop$SAMPLES, 
                         pop.id.start, pop.id.end),
                  levels = pop.levels, ordered = T
                  )  
    }

  # split that genepop dataframe by populations (with the population vector)
  genepop.split <- split(genepop, pop)
  
# Write the file in genepop format 
  newfile <- file(genepop.filename, "write")
  cat(genepop.header, "\n", file = newfile, append = TRUE)
  cat(loci, sep = "\n", file = newfile, append = TRUE)
  
  for (i in 1:length(genepop.split)) {
    cat("pop\n", file = newfile, append = TRUE)
    write.table(
      genepop.split[[i]],
      file = newfile, 
      append = TRUE, 
      col.names = FALSE,
      row.names = FALSE, 
      sep = " ", 
      quote = FALSE
      )
  }
  close(newfile)

  invisible(cat(sprintf("Genepop file name:\n%s\n\nWritten in the directory:\n%s", 
  genepop.filename, getwd())))
}
