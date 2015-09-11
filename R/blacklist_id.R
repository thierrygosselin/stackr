#' @name missing_genotypes
#' @title Missing genotypes information and blacklist of individuals based 
#' on a threshold of missing genotype.
#' @description Missing genotypes information summary per individuals and 
#' per population. Create a blacklist of individuals based 
#' on a threshold of missing genotype. 
#' This function accept a whitelist of loci to create a blacklist 
#' of individuals before or after filtering of loci. 
#' Paralogs are automatically removed from STACKS haplotype file.
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci with a column header 
#' 'LOCUS'. If the whitelist is written in the directory 
#' \code{whitelist.loci = "whitelist.txt"}. If the whitelist is in 
#' the global environment \code{whitelist.loci = whitelist.1000loci}
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param missing.geno.threshold (integer) Percentage of missing genotype 
#' allowed per individuals. e.g. for a maximum of 30% of missing genotype 
#' per individuals \code{missing.geno.threshold = 30}.
#' @return a list with 4 data frames: $missing.genotypes.ind,
#' $missing.genotypes.pop, $blacklisted.id, $plot.missing.
#' @details For the plot, to see the information with the popualtion in 
#' different facet, use \code{+facet_wrap(~POP_ID, nrow=2,ncol=5)} 
#' after the object of the plot, e.g. \code{fig <- missing$plot.missing}, to have
#' facet by pop \code{fig +facet_wrap(~POP_ID, nrow=2,ncol=5)}
#' where \code{nrow} and \code{ncol} in this exeample would spread 
#' the 10 populations on 2 rows and 5 columns.
#' @export
#' @rdname missing_genotypes
#' @import reshape2
#' @import dplyr
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

missing_genotypes <- function(haplotypes.file,
                              whitelist.loci = NULL,
                              pop.id.start, pop.id.end,
                              missing.geno.threshold
) {
  
  CONSENSUS <- NULL
  CONSENSUS_MAX <- NULL
  MISSING <- NULL
  TOTAL <- NULL
  GENOTYPED <- NULL
  PERC_MISSING <- NULL
  ..scaled.. <- NULL
  MISSING_GROUP <- NULL
  N_IND <- NULL
  
  # Haplotype file--------------------------------------------------------------
  haplotype <- read_tsv(file = haplotypes.file, col_names = T) %>%
    select(-Cnt) %>% 
    rename(LOCUS = `Catalog ID`) %>%
    tidyr::gather(INDIVIDUALS, HAPLOTYPES, -LOCUS)
  
  
  # Whitelist-------------------------------------------------------------------
  if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "TRUE") {
    message("Using the whitelist from the directory")
    whitelist <- read_tsv(whitelist.loci, col_names = T)
  } else if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "FALSE") {
    message("Using whitelist from your global environment")
    whitelist <- whitelist.loci
  } else {
    message("No whitelist")
    whitelist <- NULL
  }
  
  
  if (is.null(whitelist.loci) == TRUE) {
    
    # Combination 1: No whitelist ----------------------------------------------
    haplotype <- haplotype    
    
  } else if (is.null(whitelist.loci) == FALSE) {
    
    # Combination 2: Using whitelist -------------------------------------------
    
    # just whitelist.loci
    haplotype <- haplotype %>% 
      semi_join(whitelist, by = "LOCUS") %>% 
      arrange(LOCUS)
  }    
  
  # dump unused object
  whitelist <- NULL
  
  # Paralogs-------------------------------------------------------------------
  blacklist.loci.paralogs <- haplotype %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(LOCUS) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    group_by(LOCUS) %>%
    select(LOCUS) %>%
    distinct(LOCUS) %>% 
    arrange(LOCUS)
  
  nparalogs <- stri_join("Found and/or removed", n_distinct(blacklist.loci.paralogs$LOCUS), "paralogs", sep = " ")
  message(nparalogs)
  
  # Consensus-------------------------------------------------------------------
  
  blacklist.loci.consensus <- haplotype %>%
    mutate(CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus")) %>%
    group_by(LOCUS) %>%
    summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
    filter(CONSENSUS_MAX > 0) %>%
    group_by(LOCUS) %>%
    select (LOCUS) %>%
    distinct(LOCUS) %>%
    arrange(LOCUS)
  
  # Missing -------------------------------------------------------------------
  
  missing.genotypes.ind <- haplotype %>%
    filter(subset =! LOCUS %in% blacklist.loci.consensus$LOCUS) %>%
    filter(subset =! LOCUS %in% blacklist.loci.paralogs$LOCUS) %>%
    group_by(INDIVIDUALS) %>%
    summarise(
      TOTAL = n(),
      GENOTYPED = length(INDIVIDUALS[HAPLOTYPES != "-"]),
      MISSING = length(INDIVIDUALS[HAPLOTYPES == "-"]),
      PERC_MISSING = round((MISSING/TOTAL)*100, 2)
    ) %>% 
    mutate(POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end)) %>% 
    select(POP_ID, INDIVIDUALS, TOTAL, GENOTYPED, MISSING, PERC_MISSING)
  
  write.table(missing.genotypes.ind, 
              "missing.genotypes.ind.tsv", 
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F)
  
  # After applying threshold
  # missing.geno.threshold<-10
  blacklisted.id <- missing.genotypes.ind %>% 
    filter(PERC_MISSING > missing.geno.threshold) %>% 
    select(INDIVIDUALS) %>% 
    distinct(INDIVIDUALS)
  
  filename <- paste("blacklisted.id", missing.geno.threshold, "txt", sep = ".")
  write.table(blacklisted.id, 
              filename, 
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F)
  
  # Figure : Density of missing genotype before and after threshold
  plot.missing <- ggplot(data = missing.genotypes.ind, aes(x = PERC_MISSING, na.rm = T))+
    geom_line(aes(y = ..scaled.., color = POP_ID), stat = "density", size = 0.7, adjust = 1)+
    #     scale_colour_manual(name = "Filters", values = c("black", "blue"))+
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))+
    labs(x = "Missing genotype (Percent)")+
    labs(y = "Samples (scaled density)")+
    expand_limits(x = 0)+
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
    )
  
  
  # Create a table of missing individual by pop 
  
  missing.genotypes.pop <- missing.genotypes.ind %>% 
    mutate(
      MISSING_GROUP = cut(PERC_MISSING, 
                          breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), 
                          labels = c("0-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100"))
    ) %>% 
    group_by(POP_ID) %>% 
    mutate(N_IND = length(INDIVIDUALS)) %>% 
    group_by(POP_ID, MISSING_GROUP) %>%
    summarise(
      N = mean(N_IND),
      INDIVIDUALS = n()
    ) %>% 
    dcast(POP_ID + N ~ MISSING_GROUP, value.var = "INDIVIDUALS")
  
  write.table(missing.genotypes.pop, 
              "missing.genotypes.pop.tsv", 
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F)
  
  
  res <- list()
  res$missing.genotypes.ind <- missing.genotypes.ind
  res$missing.genotypes.pop <- missing.genotypes.pop
  res$blacklisted.id <- blacklisted.id
  res$plot.missing <- plot.missing
  
  invisible(cat(sprintf(
    "
Missing genotypes information per ind filename: missing.genotypes.ind.tsv
Missing genotypes information per pop filename: missing.genotypes.pop.tsv
The missing genotype percentage threshold: >= %s
The number of individuals blacklisted = %s
Blacklist individual filename:  %s
A plot with missing genotypes information was created: $plot.missing

Written in the directory:
  %s", 
    missing.geno.threshold,
    n_distinct(blacklisted.id$INDIVIDUALS),
    filename,
    getwd()
  )))
  
  return(res)
}
