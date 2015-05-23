#### IBM analysis: Identity-by-Missingness #####################################

# IBM table preparation function
ibm_table_prep <- function(plink.mds, pop.map) {
  
  read_table(plink.mds, 
                   col_names = T,
                   col_types = "ccidddd"
                   ) %>%
  select(-c(IID), INDIVIDUALS=FID) %>%
  left_join(
    read_tsv(pop.map, 
             col_names = T, 
             col_types = "ccic"), 
    by = "INDIVIDUALS") %>%
  mutate(
    POP_ID = as.factor(POP_ID),
    LANES = as.factor(LANES),
    SEQUENCER = as.factor(SEQUENCER)
    )
}

# IBM figure function
figure_ibm <- function(data, aes.colour, figure.title) {
  graph <- ggplot(data, aes(x = C1, y = C2))+
  geom_point(aes.colour)+
  labs(title = figure.title)+
  labs(x = "PC1")+
  labs(y = "PC2")+
  theme(axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"),
        legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
        legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
        strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
        )
}
#### MISSING GENOTYPES BEFORE FILTERS ##########################################

# Missing genotype before filter function
missing_genotype <- function(haplotype.file, pop.levels, pop.id.start, pop.id.end) {
  
  batch_1.haplotypes <- read_tsv(haplotype.file, col_names = T)
  
  NUMBER_CATALOG_LOCI <- n_distinct(batch_1.haplotypes$`Catalog ID`)

  missing.genotypes.before.filter <- batch_1.haplotypes %>%
    rename(LOCUS =`Catalog ID`) %>%
    melt(id.vars = c("LOCUS", "Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES") %>%
    filter(HAPLOTYPES != "consensus" & HAPLOTYPES == "-") %>%
    group_by(SAMPLES) %>%
    tally() %>%
    mutate(
      MISSING_LOCI_PERC = round(((n/NUMBER_CATALOG_LOCI)*100),2),
      POP_ID = factor(str_sub(SAMPLES, pop.id.start, pop.id.end), levels = pop.levels, ordered = T)
      )
  missing.genotypes.before.filter
}

# Figure of missing genotype before filter function
figure_missing_genotype <- function(data, aes.colour){
  ggplot(data, aes(x = MISSING_LOCI_PERC, na.rm = T))+
    geom_line(aes.colour, stat = "density", size = 0.7, adjust = 1)+
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))+
    labs(x = "Missing genotype (Percent)")+
    labs(y = "Samples (scaled density)")+
    expand_limits(x = 0)+
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold"))
}

# Blacklist individual function
# Create a blacklist of individuals based on their missing info
blacklist_individuals <- function(missing.file, missing.loci.perc.threshold) {
  blacklisted.id <- missing.file %>%
    filter(MISSING_LOCI_PERC >= missing.loci.perc.threshold) %>%
    select(SAMPLES) %>%
    distinct(SAMPLES) %>%
    arrange(SAMPLES)
  
  invisible(
        cat(
         sprintf(
"The missing percentage threshold: >= %s
The number of individuals blacklisted = %s", 
        missing.loci.perc.threshold,
        n_distinct(blacklisted.id$SAMPLES)
        )
        )
      )
  blacklisted.id
}
#### MISSING GENOTYPES AFTER FILTERS #########################################

#' Missing genotype AFTER filter function.
#' This function filters the 'batch_x.haplotypes.tsv' file
#' with a whitelist of loci.
#' @param haplotype.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist A whitelist of loci containing a column header 'LOCUS'.
#' @param pop.levels A character string with your populations ordered.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.

missing_genotype_filtered <- function(haplotype.file, whitelist, pop.levels, 
                                      pop.id.start, pop.id.end, blacklist.genotypes) {
  
  if (is.vector(haplotype.file) == "TRUE") {
      batch_1.haplotypes <- read_tsv(haplotype.file, col_names = T) %>%
        rename(LOCUS =`Catalog ID`)
    } else {
      batch_1.haplotypes <- haplotype.file
  }

  NUMBER_CATALOG_LOCI <- n_distinct(batch_1.haplotypes$LOCUS)

  if (is.vector(whitelist) == "TRUE") {
      whitelist.file <- read_tsv(whitelist, col_names = T)
    } else {
      whitelist.file <- whitelist
  }
  
  if (is.vector(blacklist.genotypes) == "TRUE") {
      blacklist.genotypes <- read_tsv(blacklist.genotypes, col_names = T)
    } else {
      blacklist.genotypes <- blacklist.genotypes
  }
    
  NUMBER_LOCI_FILTERED <- n_distinct(whitelist.file$LOCUS)


  # filtering haplotypes file with the whitelist
  missing.genotypes.after.filter <- batch_1.haplotypes %>%
    melt(id.vars = c("LOCUS", "Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES") %>%
    mutate(
      HAPLOTYPES = ifelse(SAMPLES %in% blacklist.genotypes$SAMPLES & LOCUS %in% blacklist.genotypes$LOCUS, "-", HAPLOTYPES )
    ) %>%
    filter(HAPLOTYPES != "consensus" & HAPLOTYPES == "-") %>%
    group_by(SAMPLES) %>%
    tally() %>%
    mutate(
      MISSING_LOCI_PERC = round(((n / NUMBER_CATALOG_LOCI) * 100), 2),
      POP_ID = factor(str_sub(SAMPLES, pop.id.start, pop.id.end), levels = pop.levels, ordered = T),
      GROUP = rep("pre-filters", n())
      ) %>%
    rbind(
      batch_1.haplotypes %>%
        right_join(whitelist.file, 
                   by = "LOCUS"
                   ) %>%
        melt(id.vars = c("LOCUS", "Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES") %>%
        filter(HAPLOTYPES != "consensus" & HAPLOTYPES == "-") %>%
        group_by(SAMPLES) %>%
        tally() %>%
        mutate(
          MISSING_LOCI_PERC = round(((n / NUMBER_LOCI_FILTERED) * 100), 2),
          POP_ID = factor(str_sub(SAMPLES, pop.id.start, pop.id.end), levels = pop.levels, ordered = T),
          GROUP = rep("post-filters", n())
          )
      ) %>%
    mutate(GROUP = factor(GROUP, levels = c("pre-filters", "post-filters"), ordered = T))

missing.genotypes.after.filter
}


# Figure function: Density of missing genotype before and after filter 1 big group
figure_missing_genotype_before_after <- function(data, aes.colour) {

  ggplot(data = data, aes(x = MISSING_LOCI_PERC, na.rm = T))+
    geom_line(aes.colour, stat = "density", size = 0.7, adjust = 1)+
    scale_colour_manual(name = "Filters", values = c("black", "blue"))+
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
}

# Create a table of missing individual by pop 
missing_genotypes_individuals_after_filter_pop_table <- function(missing.file) {

  total.individuals.pop <- missing.file %>%
    filter(GROUP == "post-filters") %>%
    group_by(POP_ID) %>%
    summarise(N = length(SAMPLES))

  individuals.blacklisted <- missing.file %>% filter(GROUP == "post-filters") %>% 
  mutate(MISSING_GROUP = cut(MISSING_LOCI_PERC, breaks = c(0, 10, 20, 
    30, 40, 50, 60, 70, 80, 90, 100), labels = c("0-10", "11-20", "21-30", 
    "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100"))) %>% 
  group_by(POP_ID, MISSING_GROUP) %>% summarise(COUNT = n()) %>% dcast(POP_ID ~ 
  MISSING_GROUP, value.var = "COUNT") %>% full_join(total.individuals.pop, 
  by = "POP_ID")
  
  individuals.blacklisted
}

blacklist_individuals_after_filters <- function(missing.file, 
                                               missing.loci.perc.threshold,
                                               filename) {

  if (is.vector(missing.file) == "TRUE") {
      missing.file <- read_tsv(missing.file, col_names = T)
    } else {
      missing.file <- missing.file
  }
  
  blacklisted.id <- missing.file %>%
    filter(GROUP == "post-filters") %>%
    filter(MISSING_LOCI_PERC >= missing.loci.perc.threshold) %>%
    select(SAMPLES) %>%
    distinct(SAMPLES) %>%
    arrange(SAMPLES)

  write.table(blacklisted.id, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)

  invisible(cat(sprintf(
  "The missing percentage threshold: >= %s
  The number of individuals blacklisted = %s\n
  Blacklist id filename:
  %s\n
  Written in the directory:
  %s", 
  missing.loci.perc.threshold,
  n_distinct(blacklisted.id$SAMPLES),
  filename, getwd()
  )))
  
  blacklisted.id
}

