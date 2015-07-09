# Write a genepop file from STACKS haplotype file
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE_1", "ALLELE_2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX"))


#' @name haplo2genepop
#' @title Use the batch_x.haplotypes.tsv file to write a genpop file
#' @description This function can first filter the haplotypes file 
#' with a whitelist of loci
#' and a blacklist of individuals (optional).
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param genepop.filename The name of the file written to the directory.
#' @param genepop.header The first line of the Genepop file 
#' (ex: genepop of 10 pop of Acipenser fulvescens filtered MAF 0.05).
#' @param pop.levels An optional character string with your populations ordered.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' 
#' 
#' @export
#' @rdname haplo2genepop
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_pad
#' @references Rousset F (2008) genepop’007: a complete re-implementation
#' of the genepop software for Windows and Linux.
#' Molecular Ecology Resources, 8, 103-106.
#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}



haplo2genepop <- function(haplotypes.file, 
                          whitelist.loci = NULL, 
                          blacklist.id = NULL, 
                          genepop.filename, genepop.header, 
                          pop.levels, pop.id.start, pop.id.end) {
  
  # Whitelist
  if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "TRUE") {
    message("Using the whitelist from the directory")
    whitelist <- read_tsv(whitelist.loci, col_names = T) %>%
      rename(Catalog.ID = LOCUS)
  } else if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "FALSE") {
    message("Using whitelist from your global environment")
    whitelist <- whitelist.loci %>%
      rename(Catalog.ID = LOCUS)
  } else {
    message("No whitelist")
    whitelist <- NULL
  }
  
  
  # Blacklist
  if (missing(blacklist.id) == "FALSE" & is.vector(blacklist.id) == "TRUE") {
    message("Using the blacklist of id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)    
  } else if (missing(blacklist.id) == "FALSE" & is.vector(blacklist.id) == "FALSE") {
    message("Using blacklist of id from your global environment")
    blacklist.id <- blacklist.id
    
  } else {
    message("No blacklist of id")
    blacklist.id <- NULL
  }
  
  
  # Haplotype file
  haplotype <- read_tsv(haplotypes.file, col_names = T) %>% 
    rename(Catalog.ID = `Catalog ID`) %>%
    gather(INDIVIDUALS, HAPLOTYPES, -c(Catalog.ID, Cnt)) %>%
    arrange(Catalog.ID)
  
  
  
  if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == TRUE) {
    
    message("Combination 1: No whitelist and No blacklist")
    
    # No filter
    haplotype.no.filter <- haplotype    
    
    # Creates a vector containing the loci name
    loci <- unique(haplotype.no.filter$Catalog.ID)
    data <- haplotype.no.filter
    
  } else if (is.null(whitelist.loci) == FALSE & is.null(blacklist.id) == TRUE) {
    
    message("Combination 2: Using whitelist, but No blacklist")
    
    # just whitelist.loci, NO Blacklist of individual
    haplotype.whitelist.loci <- haplotype %>%
      semi_join(whitelist, by = "Catalog.ID") %>% # instead of right_join
#       might result in less loci then the whitelist, but this make sure
#   we don't end up with unwanted or un-caracterized locus, e.g. comparing
#   adults and juveniles samples... the whitelist is developed with the adults,
#   but the haplotype file comes from the juveniles... 
      arrange(Catalog.ID)
    
    # Creates a vector containing the loci name
    loci <- unique(haplotype.whitelist.loci$Catalog.ID)
    data <- haplotype.whitelist.loci
    
  } else if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == FALSE) {
    message("Combination 3: Using a blacklist of id, but No whitelist")
    
    # NO whitelist, JUST Blacklist of individual
    haplotype.blacklist <- haplotype %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(Catalog.ID)
    
    
    # Creates a vector containing the loci name
    loci <- unique(haplotype.blacklist$Catalog.ID)
    data <- haplotype.blacklist
    
  } else {
    message("Combination 4: Using a whitelist and blacklist")
    
    # whitelist.loci + Blacklist of individual
    haplotype.whitelist.blacklist <- haplotype %>%
      semi_join(whitelist, by = "Catalog.ID") %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
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
  
  
  # Paralogs..
  message("Looking for paralogs...")
  
  paralogs <- data %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(Catalog.ID) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    group_by(Catalog.ID) %>%
    select(Catalog.ID) %>%
    distinct(Catalog.ID)
  
  nparalogs <- stri_join("Found and removed", n_distinct(paralogs$Catalog.ID), "paralogs", sep = " ")
  message(nparalogs)
  
  message("Haplotypes into Genepop factory ...")
  
  genepop <- data %>%
    anti_join(paralogs, by = "Catalog.ID") %>%
    filter(HAPLOTYPES != "consensus") %>%    
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
    melt(id.vars = c("Catalog.ID","Cnt", "INDIVIDUALS"),
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
    dcast(Catalog.ID + INDIVIDUALS ~ ALLELE, value.var = "NUCLEOTIDES") %>%
    unite(GENOTYPE, ALLELE_1:ALLELE_2, sep = "") %>%
    dcast(INDIVIDUALS ~ Catalog.ID, value.var = "GENOTYPE") %>%
    mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))
  
  # Create a vector with the population ordered by levels
  if (missing(pop.levels) == "TRUE") {
    pop <- substr(genepop$INDIVIDUALS, pop.id.start, pop.id.end)
  } else {
    pop <- factor(substr(genepop$INDIVIDUALS, 
                         pop.id.start, pop.id.end),
                  levels = pop.levels, ordered = T)
    
  }
  
  # split that genepop dataframe by populations (with the population vector)
  genepop.split <- split(genepop, pop)
  
  # Write the file in genepop format 
  message("Saving the file in your working directory...")
  
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
