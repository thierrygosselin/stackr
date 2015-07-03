# Write a strataG gtypes file from STACKS haplotypes file

# to get rid of notes in build check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE_1", "ALLELE_2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX", "contains", "other"))


#' @name haplo2gtypes
#' @title Convert between batch_x.haplotypes.tsv and \code{gtypes} objects
#' @description This function can first filter the haplotypes file 
#' with a whitelist of loci
#' and a blacklist of individuals (optional). Then it will convert the file
#' to a strataG \code{gtypes} objects.
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param filename The name of the file written to the directory (optional).
#' @param pop.levels An optional character string with your populations ordered.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param description a string naming or describing this dataset.
#' @return a \code{gtypes} objects.
#' @export
#' @rdname haplo2gtypes
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_pad
# @importFrom strataG.devel df2gtypes
#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @seealso strataG is available on CRAN and https://github.com/EricArcher/
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

haplo2gtypes <- function(haplotypes.file,
                         whitelist.loci = NULL, 
                         blacklist.id = NULL, 
                         filename, 
                         pop.levels,
                         pop.id.start, 
                         pop.id.end, 
                         description) {
  
  
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
  haplotype <- read_tsv(file = "batch_1.haplotypes.tsv", col_names = T) %>%
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
  
  message("Haplotypes into conversion to strataG gtypes factory ...")
  
  # Haplo prep
  # Remove paralogs
  # Remove consensus loci
  
  haplo.prep <- data %>%
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
      NUCLEOTIDES = as.numeric(factor(NUCLEOTIDES)),
      NUCLEOTIDES = NUCLEOTIDES-1, # this removes 1 levels to enable missing values = 0
      NUCLEOTIDES = str_pad(NUCLEOTIDES, 3, side = "left", pad = "0")
    ) %>%
    select(-Cnt) %>%
    dcast(INDIVIDUALS ~ Catalog.ID + ALLELE, value.var = "NUCLEOTIDES") %>%
    mutate(
      POP_ID = factor(substr(INDIVIDUALS, 
                             pop.id.start, pop.id.end),
                      levels = pop.levels, ordered = T)
    ) %>%
    select(INDIVIDUALS, POP_ID, contains("ALLELE"))  
  
  
  # Write file to directory (optional)
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(haplo.prep, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
  }
  
  # convert to gtypes
  res <- strataG.devel::df2gtypes(x = haplo.prep,
                                  ploidy = 2,
                                  id.col = "INDIVIDUALS",
                                  strata.col = "POP_ID",
                                  loc.col = 3, 
                                  sequences = NULL,
                                  description = description
  )
  
# Message at the end
  invisible(cat(sprintf(
    "%s\n
Working directory:
%s",
    saving, getwd()
  )))
  return(res)
}
