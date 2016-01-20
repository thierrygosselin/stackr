#' @name missing_genotypes
#' @title Missing genotypes information and blacklist of individuals based 
#' on a threshold of missing genotype.
#' @description Use this function to have a summary of missing genotypes 
#' per individuals and per population. Density distribution plot of misssing
#' information is produced. This function also create a blacklist
#' of individuals based on the desired threshold of missing genotype. 
#' When a whitelist of loci is provided, \code{_filtered} will be appended 
#' to the filename. Locus with more than 2 alleles by individual (paralogs and/or sequencing errors) can be deleted or 
#' genotypes with more than 2 alleles by individual can be erased.
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci with a column header 
#' 'LOCUS'. If the whitelist is written in the directory 
#' \code{whitelist.loci = "whitelist.txt"}. If the whitelist is in 
#' the global environment \code{whitelist.loci = whitelist.1000loci}
#' @param erase (character) Default when missing argument = \code{"loci"}. Loci with more than 2 alleles 
#' (paralogs and/or sequencing errors) will be removed. \code{"genotypes"}: genotypes 
#' with more than 2 alleles by individual (paralogs and/or sequencing errors) will be erased.
#' Keeping the loci for other individuals.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels A character string with your populations ordered.
#' @param pop.labels (optional character string) If you need to rename 
#' pop/sampling sites in \code{pop.levels} or combined sites/pop
#' into a different names, here is the place (e.g. 3 sites: 
#' c("ONT", "QUE", "SAS") and you want to combine "ONT" and "SAS" into the same 
#' grouping "CAN" use : c("CAN", "QUE", "CAN").
#' @param missing.geno.threshold (integer) Percentage of missing genotype 
#' allowed per individuals. e.g. for a maximum of 30% of missing genotype 
#' per individuals \code{missing.geno.threshold = 30}.
#' @return a list with 4 data frames: $missing.genotypes.ind,
#' $missing.genotypes.pop, $blacklisted.id, $plot.missing. A 5th data frame is
#' added when \code{erase = genotypes} is selected: $blacklist.genotypes.
#' @details For the plot, to see the information with the popualtion in 
#' different facet, use \code{+facet_wrap(~POP_ID, nrow=2,ncol=5)} 
#' after the object of the plot, e.g. \code{fig <- missing$plot.missing}, to have
#' facet by pop \code{fig +facet_wrap(~POP_ID, nrow=2,ncol=5)}
#' where \code{nrow} and \code{ncol} in this exeample would spread 
#' the 10 populations on 2 rows and 5 columns.

#' @examples
#' \dontrun{
#' Missing genotypes before filtering
#' I want a blacklist of individuals with more than 30% missing data.
#' missing.info.before.filter <- missing_genotypes(
#' haplotypes.file = "batch_1.haplotypes.tsv",
#' pop.id.start = 5,
#' pop.id.end = 7, 
#' pop.levels = c("QUE", "ONT", "MAN"),
#' erase = "loci", 
#' missing.geno.threshold = 30
#' )
#' 
#' To update the blacklist with the filtered dataset:
#' missing.info.after.filter <- missing_genotypes(
#' haplotypes.file = "batch_1.haplotypes.tsv", 
#' pop.id.start = 5, 
#' pop.id.end = 7, 
#' pop.levels = c("QUE", "ONT", "MAN"),
#' whitelist.loci = "whitelist.loci.txt", 
#' erase = "loci", 
#' missing.geno.threshold = 30
#' )
#' }

#' @export
#' @rdname missing_genotypes
#' @import reshape2
#' @import dplyr
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

missing_genotypes <- function(haplotypes.file,
                              whitelist.loci = NULL,
                              erase,
                              pop.id.start, pop.id.end, pop.levels, pop.labels, 
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
  
  # Checking for missing and/or default arguments ******************************
  if (missing(haplotypes.file)) stop("haplotypes.file required")
  if (missing(whitelist.loci)) whitelist.loci <- NULL # no Whitelist
  if (missing(erase)) erase <- NULL # no genotype to erase
  if (missing(pop.levels)) stop("pop.levels required")
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels
  if (missing(pop.id.start)) stop("pop.id.start required")
  if (missing(pop.id.end)) stop("pop.id.end required")
  if (missing(missing.geno.threshold)) stop("missing.geno.threshold required")
  
  
  
  # Haplotype file--------------------------------------------------------------
  haplotype <- suppressWarnings(
    read_tsv(file = haplotypes.file, col_names = T, na = "-") %>%
      select(-Cnt) %>% 
      rename(LOCUS = `Catalog ID`) %>%
      tidyr::gather(INDIVIDUALS, HAPLOTYPES, -LOCUS) %>% 
      mutate(HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = "-"))
  )
  
  
  # Whitelist-------------------------------------------------------------------
  if (is.null(whitelist.loci)) {
    message("No whitelist to apply to the haplotypes file")
    whitelist.loci <- NULL
    haplotype <- haplotype
  } else if (is.vector(whitelist.loci) == "TRUE") {
    message("Using the whitelist from the directory")
    haplotype <- haplotype %>% 
      semi_join(
        read_tsv(whitelist.loci, col_names = TRUE), 
        by = "LOCUS") %>% 
      arrange(LOCUS)
  } else {
    message("Using whitelist from your global environment")
    haplotype <- haplotype %>% 
      semi_join(whitelist.loci, by = "LOCUS") %>% 
      arrange(LOCUS)
  }
  
  # Consensus-------------------------------------------------------------------
  
  blacklist.loci.consensus <- haplotype %>%
    filter(HAPLOTYPES == "consensus") %>% 
    select (LOCUS) %>%
    distinct(LOCUS) %>%
    arrange(LOCUS)
  
  haplotype <- haplotype %>% 
    filter(subset =! LOCUS %in% blacklist.loci.consensus$LOCUS)

  # Paralogs-------------------------------------------------------------------
  message("Looking for paralogs...")
  
  if(erase == "loci"){
    message("Loci with more than 2 alleles will be removed")
    
    paralogs <- haplotype %>%
      mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
      group_by(LOCUS) %>%
      summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
      filter(POLYMORPHISM_MAX > 1) %>%
      group_by(LOCUS) %>%
      select(LOCUS) %>%
      distinct(LOCUS)
    
    haplotype <- suppressWarnings(
      haplotype %>%
        filter(subset =! LOCUS %in% paralogs$LOCUS)
    )
    message(stri_join("Found and/or removed", n_distinct(paralogs$LOCUS), "paralogs", sep = " "))
    
  } else {
    message("Erasing genotypes with more than 2 alleles")
    
    # get the number of genotypes...
    haplo.number <- haplotype %>%
      filter(HAPLOTYPES != "-") %>%
      select(HAPLOTYPES)
    
    haplotype <- haplotype %>% 
      mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/"))
    
    erased.genotype.number <- length(haplotype$INDIVIDUALS[haplotype$POLYMORPHISM > 1])
    total.genotype.number.haplo <- length(haplo.number$HAPLOTYPES)
    percent.haplo <- paste(round(((erased.genotype.number/total.genotype.number.haplo)*100), 2), "%", sep = " ")
    message(stri_paste("Out of a total of ", total.genotype.number.haplo, " genotypes, ", percent.haplo, " (", erased.genotype.number, ")"," will be erased"))
    
    # keep track of blacklisted genotypes
    blacklist.genotypes <- haplotype %>% 
      filter(POLYMORPHISM > 1) %>% 
      select(INDIVIDUALS, LOCUS)
    
    if (is.null(whitelist.loci)) {
      blacklist.genotypes.filename <- paste("blacklist.genotypes", missing.geno.threshold, "txt", sep = ".")
      write.table(blacklist.genotypes, blacklist.genotypes.filename, 
                  sep = "\t", row.names = F, col.names = T, quote = F)
    } else {
      blacklist.genotypes.filename <- paste("blacklist.genotypes.filtered", missing.geno.threshold, "txt", sep = ".")
      write.table(blacklist.genotypes, blacklist.genotypes.filename, 
                  sep = "\t", row.names = F, col.names = T, quote = F)
    }
    message("A blacklist of genotypes was written in your working directory")
    
    # Erasing genotype with the blacklist
    message("Erasing... Erasing...")
    haplotype <- suppressWarnings(
      haplotype %>%
        mutate(HAPLOTYPES = ifelse(POLYMORPHISM > 1, "-", HAPLOTYPES))
    )
  }
  
  
  # Missing -------------------------------------------------------------------
  missing.genotypes.ind <- suppressWarnings(
    haplotype %>%
      group_by(INDIVIDUALS) %>%
      summarise(
        TOTAL = n(),
        GENOTYPED = length(INDIVIDUALS[HAPLOTYPES != "-"]),
        MISSING = length(INDIVIDUALS[HAPLOTYPES == "-"]),
        PERC_MISSING = round((MISSING/TOTAL)*100, 2)
      ) %>% 
      mutate(
        POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = F), levels = pop.labels, ordered =T),
        POP_ID = droplevels(POP_ID)
      ) %>% 
      select(POP_ID, INDIVIDUALS, TOTAL, GENOTYPED, MISSING, PERC_MISSING)
  )
  
  if (is.null(whitelist.loci)) {
    write.table(missing.genotypes.ind, "missing.genotypes.ind.tsv", 
                sep = "\t", row.names = F, col.names = T, quote = F)
    mis.geno.ind.filename <- "missing.genotypes.ind.tsv"
  } else {
    write.table(missing.genotypes.ind, "missing.genotypes.ind.filtered.tsv", 
                sep = "\t", row.names = F, col.names = T, quote = F)
    mis.geno.ind.filename <- "missing.genotypes.ind.filtered.tsv"
  }
  
  # After applying threshold
  # missing.geno.threshold<-10
  blacklisted.id <- missing.genotypes.ind %>% 
    filter(PERC_MISSING > missing.geno.threshold) %>% 
    select(INDIVIDUALS) %>% 
    distinct(INDIVIDUALS)
  
  if (is.null(whitelist.loci)) {
    blacklist.id.filename <- paste("blacklisted.id", missing.geno.threshold, "txt", sep = ".")
    write.table(blacklisted.id, blacklist.id.filename, 
                sep = "\t", row.names = F, col.names = T, quote = F)
  } else {
    blacklist.id.filename <- paste("blacklisted.id.filtered", missing.geno.threshold, "txt", sep = ".")
    write.table(blacklisted.id, blacklist.id.filename, 
                sep = "\t", row.names = F, col.names = T, quote = F)
  }
  
  
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
    dcast(POP_ID + N ~ MISSING_GROUP, value.var = "INDIVIDUALS") %>% 
    arrange(POP_ID)
  
  if (is.null(whitelist.loci)) {
    write.table(missing.genotypes.pop, "missing.genotypes.pop.filtered.tsv", 
                sep = "\t", row.names = F, col.names = T, quote = F)
    mis.geno.pop.filename <- "missing.genotypes.pop.filtered.tsv"
  } else {
    write.table(missing.genotypes.pop, "missing.genotypes.pop.tsv", 
                sep = "\t", row.names = F, col.names = T, quote = F)
    mis.geno.pop.filename <- "missing.genotypes.pop.tsv"
  }
  
  
  res <- list()
  res$missing.genotypes.ind <- missing.genotypes.ind
  res$missing.genotypes.pop <- missing.genotypes.pop
  res$blacklisted.id <- blacklisted.id
  res$plot.missing <- plot.missing
  
  if(erase == "genotypes"){
    res$blacklist.genotypes <- blacklist.genotypes
  }
    
  invisible(cat(sprintf(
    "
Missing genotypes information per ind filename: %s
Missing genotypes information per pop filename: %s
The missing genotype percentage threshold: >= %s
The number of individuals blacklisted = %s
Blacklist individual filename:  %s
A plot with missing genotypes information was created: $plot.missing

Written in the directory:
  %s",
    mis.geno.ind.filename,
    mis.geno.pop.filename,
    missing.geno.threshold,
    n_distinct(blacklisted.id$INDIVIDUALS),
    blacklist.id.filename,
    getwd()
  )))
  
  return(res)
}
