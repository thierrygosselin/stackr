### HAPLOTYPE FILE SUMMARY : consensus, mono- poly-morphic ###

haplotype_file_summary <- function(haplotype.file, pop.id.start, pop.id.end, pop.levels) {

  # Import haplotype file
  batch_1.haplotypes <- read_tsv(haplotype.file, col_names = T) %>%
    rename(LOCUS =`Catalog ID`) %>%
    melt(id.vars = c("LOCUS","Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES") %>%
    mutate(
        POP_ID = str_sub(SAMPLES, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
    )

  
  # Locus with > 2 alleles
  # Create a blacklist of catalog loci with paralogs

  paralogs <- batch_1.haplotypes %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    mutate(PARALOGS = rep("paralogs", times = n())) %>%
    select(LOCUS, POP_ID, PARALOGS)
  
  
 
  blacklist.loci.paralogs <- paralogs %>%
    group_by(LOCUS) %>%
    select (LOCUS) %>%
    distinct(LOCUS) %>%
    arrange(LOCUS)
  write.table(blacklist.loci.paralogs, "blacklist.loci.paralogs.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

  paralogs
  
  # Locus with concensus alleles
  consensus <- batch_1.haplotypes %>%
    mutate(CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus")) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
    filter(CONSENSUS_MAX > 0) %>%
    mutate(CONSENSUS = rep("consensus", times = n())) %>%
    select(LOCUS, POP_ID, CONSENSUS)
  

# Number of loci in the catalog with consensus alleles
  
  blacklist.loci.consensus <- consensus %>%
    group_by(LOCUS) %>%
    select (LOCUS) %>%
    distinct(LOCUS) %>%
    arrange(LOCUS)
  write.table(blacklist.loci.consensus, "blacklist.loci.consensus.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

  consensus  
  
  # Summary dataframe
  summary <- batch_1.haplotypes %>%
    filter(subset =! LOCUS %in% consensus$LOCUS) %>%
    mutate(
      ALLELES_COUNT = stri_count_fixed(HAPLOTYPES, "/")
      ) %>%
    arrange(LOCUS) %>%
    group_by(POP_ID, LOCUS) %>%
    summarise(
      ALLELES_COUNT_SUM = sum(ALLELES_COUNT)
      ) %>%
    mutate(
      POP_LEVEL_POLYMORPHISM = ifelse(ALLELES_COUNT_SUM == 0, "mono", "poly")
      ) %>%
    summarise(
      MONOMORPHIC = length(POP_LEVEL_POLYMORPHISM[POP_LEVEL_POLYMORPHISM == "mono"]),
      POLYMORPHIC = length(POP_LEVEL_POLYMORPHISM[POP_LEVEL_POLYMORPHISM == "poly"])
      ) %>%
    full_join(
      consensus %>%
        group_by(POP_ID) %>%
        summarise(CONSENSUS = n_distinct(LOCUS)),
      by = "POP_ID"
      ) %>%
    group_by(POP_ID) %>%
    mutate(
      TOTAL = MONOMORPHIC + POLYMORPHIC + CONSENSUS
      ) %>%
    full_join(
      paralogs %>%
        group_by(POP_ID) %>%
        summarise(PARALOGS = n_distinct(LOCUS)),
      by = "POP_ID"
      ) %>%
    mutate(
      MONOMORPHIC_PROP = round(MONOMORPHIC/TOTAL, 4),
      POLYMORPHIC_PROP = round(POLYMORPHIC/TOTAL, 4),
      PARALOG_PROP = round(PARALOGS/TOTAL, 4)
      )
  write.table(summary, "haplotype.catalog.loci.summary.pop.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
  
  invisible(
        cat(
         sprintf(
"The number of loci in the catalog = %s LOCI
The number of putative paralogs loci in the catalog (> 2 alleles) = %s LOCI
The number of loci in the catalog with consensus alleles = %s LOCI
3 files were written in this directory: %s
1. blacklist.loci.paralogs.tsv
2. blacklist.loci.consensus.tsv
3. haplotype.catalog.loci.summary.pop.tsv", 
            n_distinct(batch_1.haplotypes$LOCUS),
            n_distinct(paralogs$LOCUS),
            n_distinct(consensus$LOCUS),
            getwd()
        )
        )
      )
  summary
}


paralogs <- function(haplotype.file, pop.id.start, pop.id.end, pop.levels) {
  
  paralogs <- read_tsv(haplotype.file, col_names = T) %>%
    rename(LOCUS =`Catalog ID`) %>%
    melt(id.vars = c("LOCUS","Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES") %>%
    mutate(
        POP_ID = str_sub(SAMPLES, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
        ) %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(LOCUS) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    mutate(PARALOGS = rep("paralogs", times = n())) %>%
    select(LOCUS, PARALOGS)
    
}


# Function: Haplotypes statistics summary 'batch_x.hapstats.tsv'

hapstats_summary <- function(data, pop.num, pop.col.types, pop.integer.equi, pop.levels) {
  
  skip.lines <- pop.num + 1
  
  if(pop.col.types == "integer"){
    col.types = "iiciiiiddddddc"
  } 
  if(pop.col.types == "character") {
    col.types = "iiciciiddddddc"
  } else {
    col.types = NULL
  }
  hapstats <- read_tsv(data,
                       na = "NA",
                       skip = skip.lines,
                       progress = interactive(),
                       col_names = c("BATCH_ID", "LOCUS", "CHR", "BP", "POP_ID", "N", "HAPLOTYPE_CNT", "GENE_DIVERSITY", "SMOOTHED_GENE_DIVERSITY", "SMOOTHED_GENE_DIVERSITY_PVALUE", "HAPLOTYPE_DIVERSITY", "SMOOTHED_HAPLOTYPE_DIVERSITY", "SMOOTHED_HAPLOTYPE_DIVERSITY_PVALUE", "HAPLOTYPES"),
                       col_types = col.types) %>%
  mutate (
    POP_ID = stri_replace_all_fixed(POP_ID, seq(from = 1, to = pop.num, by = 1), pop.integer.equi, vectorize_all=F),
    POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
    ) %>%
    arrange(LOCUS, POP_ID)
#   separate(HAPLOTYPES, c("ALLELE_P", "ALLELE_Q"), sep = "/", extra = "error", remove = F) %>%
}

# Figure: Density distribution of diversity (Gene and Haplotypes)
figure_distribution_diversity <- function(data, aes.x, aes.colour,  x.title, y.title) {
  ggplot(hapstats.summary, aes.x)+
  geom_line(aes.colour, stat = "density", adjust = 0.8)+
#   scale_colour_manual(name = "Populations", values = colour_palette_sites.pink, breaks = c("BUR", "GRA", "GUL", "LLI", "ANG", "WEI", "HAY", "GOD"))+
#   geom_density(aes(fill=POP_ID, color=NA), alpha=0.4)+
#   scale_fill_manual(name="Populations", values=colour_palette_sites.pink, breaks = c("BUR", "GRA", "GUL", "LLI", "ANG", "WEI", "HAY", "GOD"))+
#   scale_color_manual(name="Populations", values=colour_palette_sites.pink, breaks = c("BUR", "GRA", "GUL", "LLI", "ANG", "WEI", "HAY", "GOD"))+
  labs(x = x.title)+
  labs(y = y.title)+
  expand_limits(y=0)+
  theme(
    legend.position = "none",
    axis.title.x=element_text(size=12, family="Helvetica",face="bold"),
    axis.title.y=element_text(size=12, family="Helvetica",face="bold"),
    legend.title=element_text(size=12, family="Helvetica",face="bold"),
    legend.text=element_text(size=12, family="Helvetica",face="bold"),
    strip.text.y=element_text(angle=0,size=12, family="Helvetica",face="bold"),
    strip.text.x=element_text(size=12, family="Helvetica",face="bold")
    )
}

# Figure: Box Plot of diversity (Gene and Haplotypes)
figure_box_plot_diversity <- function(data, aes.x.y, y.title) {
  ggplot(data, aes.x.y)+
  geom_violin(trim = F)+
  geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
  stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
  labs(x = "Sampling sites")+
  labs(y = y.title)+
  theme(
     legend.position = "none",
     axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
     axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"),
     legend.title = element_text(size = 12, family = "Helvetica", face = "bold"),
     legend.text = element_text(size = 12, family = "Helvetica", face = "bold"),
     strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
     )
}

# Phistats
phistats_summary <- function(data, skip.lines) {
  
  phistats <- read_tsv(data,
    na = "NA",
    skip = skip.lines,
    progress = interactive(),
    col_types = "iiciiddddddddddddddd",
    col_names = c("BATCH_ID", "LOCUS", "CHR", "BP", "POP_COUNT", "PHI_ST", "SMOOTHED_PHI_ST", "SMOOTHED_PHI_ST_P_VALUE", "PHI_CT", "SMOOTHED_PHI_CT", "SMOOTHED_PHI_CT_P_VALUE", "PHI_SC", "SMOOTHED_PHI_SC", "SMOOTHED_PHI_SC_P_VALUE", "FST_PRIME", "SMOOTHED_FST_PRIME", "SMOOTHED_FST_PRIME_P_VALUE", "D_EST", "SMOOTHED_D_EST", "SMOOTHED_D_EST_P_VALUE")
  ) %>%
    select(-c(BATCH_ID, CHR, SMOOTHED_PHI_ST, SMOOTHED_PHI_ST_P_VALUE, SMOOTHED_PHI_CT, SMOOTHED_PHI_CT_P_VALUE, SMOOTHED_PHI_SC, SMOOTHED_PHI_SC_P_VALUE, SMOOTHED_FST_PRIME, SMOOTHED_FST_PRIME_P_VALUE, SMOOTHED_D_EST, SMOOTHED_D_EST_P_VALUE)) %>%
    melt(
      id.vars = c("LOCUS","BP","POP_COUNT"),
      variable.name = c("PHI_ST", "FST_PRIME", "D_EST"),
      value.name = "VALUE"
      )
}

#### Haplotype file erase genotype ####

#' Erase genotypes that didn't pass min coverage threshold
#' This function modify the 'batch_x.haplotypes.tsv' file
#' with a blacklist of loci and individuals.
#' @param haplotype.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param blacklist.genotypes A blacklist of loci and genotypes containing at least 2 columns header 'LOCUS' and 'SAMPLES'.

erase_genotypes <- function(haplotype.file, blacklist.genotypes, filename) {
  
  if (is.vector(haplotype.file) == "TRUE") {
      batch_1.haplotypes <- read_tsv(haplotype.file, col_names = T) %>%
        rename(LOCUS =`Catalog ID`)
    } else {
      batch_1.haplotypes <- haplotype.file
  }

  if (is.vector(blacklist.genotypes) == "TRUE") {
      blacklist.genotypes <- read_tsv(blacklist.genotypes, col_names = T)
    } else {
      blacklist.genotypes <- blacklist.genotypes
  }
    

  # filtering haplotypes file with the whitelist
  haplotypes.modified <- batch_1.haplotypes %>%
    melt(id.vars = c("LOCUS", "Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES") %>%
    mutate(
      HAPLOTYPES = ifelse(SAMPLES %in% blacklist.genotypes$SAMPLES & LOCUS %in% blacklist.genotypes$LOCUS, "-", HAPLOTYPES)
      ) %>%
    rename(`Catalog ID` = LOCUS) %>%
    dcast(`Catalog ID`+Cnt~SAMPLES, value.var="HAPLOTYPES")
  

  write.table(erase_genotypes, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)

  invisible(cat(sprintf(
  "Modified haplotypes filename:
  %s\n
  Written in the directory:
  %s", 
  filename, getwd()
  )))
haplotypes.modified
}
