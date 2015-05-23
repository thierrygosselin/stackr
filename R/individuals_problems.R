# In R
rm(list=ls())
ls()

# Load or install necessary librairies
# install.packages(c("devtools", "roxygen2", "testthat", "knitr", "reshape2", "ggplot2", "stringr", "stringi", "plyr", "dplyr", "tidyr", "ggvis", "ggdendro", "readr"))
devtools::install_github("hadley/readr")
# devtools::install_github("hadley/dplyr")
# install.packages("grid")
# install.packages("gridExtra")
# install.packages("gplots")
# install_github("thibautjombart/adegenet")
# install.packages("ape")
# install.packages("ade4")
# install.packages("qvalue")
# install.packages("HardyWeinberg")
# install.packages("diveRsity")
# install_github("jgx65/hierfstat", auth_token=Sys.getenv("GITHUB_TOKEN"))
# install_github("jgx65/hierfstat")
# install.packages("microbenchmark")
library(devtools)
library(reshape2)
library(ggplot2)
library(stringr)
library(stringi)
library(plyr)
library(dplyr) # load this package after plyr to work properly
library(tidyr)
library(readr)
library(ggvis)
library(grid)
library(ggdendro)
library(ggtree)
library(gridExtra)
library(gplots)
library(adegenet)
library(ape)
library(ade4)
# library(qvalue)
library(HardyWeinberg)
library(diveRsity)
library(hierfstat)
library(microbenchmark)
library(formatR)
library(roxygen2)
library(grid)
library(pegas)
formatR::tidy_app()

setwd("/Users/thierry/Dropbox/brook_charr_pop/01_stacks_populations/populations_10pop")


# TEST INDIVIDUALS PROBLEMS

# import the files
# HAPLOTYPES 
haplotype.file <- "batch_1.haplotypes.tsv"
pop.id.start <- 5
pop.id.end <- 7
pop.levels <- SITES_LEVELS
pop.map <- "population.map.10pop.adu.txt"

haplotypes <- read_tsv(haplotype.file, col_names = T) %>%
    rename(LOCUS =`Catalog ID`) %>%
    melt(id.vars = c("LOCUS","Cnt"), variable.name = "SAMPLES", value.name = "HAPLOTYPES") %>%
    mutate(
        POP_ID = str_sub(SAMPLES, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
    )

# comparaison du nombre ID dans la pop map et dans le fichier haplotype

ind.haplo <- haplotypes %>%
  group_by(POP_ID) %>%
  summarise(N_HAPLO = n_distinct(SAMPLES))

# POP MAP
ind.pop.map <- read_tsv(pop.map, col_names = F) %>%
          select(INDIVIDUALS=X1) %>%
          mutate(
            INDIVIDUALS = as.character(INDIVIDUALS),
            POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
            POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
            ) %>%
          group_by(POP_ID) %>%
          summarise(
            N_POPMAP = length(INDIVIDUALS)
            )

# VCF
source("/Users/thierry/Dropbox/stackr/R/read_stacks_vcf.R")
SITES_LEVELS

stacks.vcf.before.filters <- read_stacks_vcf(
  vcf.file = "batch_1.vcf",
  skip.line = 9,
  max.read.lines = 18903,
  pop.id.start = 5,
  pop.id.end = 7,
  pop.levels = SITES_LEVELS,
  filename = "stacks.vcf.before.filters.tsv"
  )

ind.vcf <- stacks.vcf.before.filters %>%
  group_by(POP_ID) %>%
  summarise(N_VCF = n_distinct(INDIVIDUALS))

# Comparaison du nombre d'individus génotypés
comp.ind <- ind.pop.map %>%
  bind_cols(ind.haplo, ind.vcf)

# sumstats
source("/Users/thierry/Dropbox/stackr/R/stacks_modification.R")
sumstats <- sumstats_prep(sumstats = "batch_1.sumstats.tsv", 
                          skip.line = 11,
                          pop.num = 10,
                          pop.col.types = "character", 
                          pop.integer.equi = c("CNC", "FRE", "GOO", "LIM", "MOO", "NIN", "SKY", "SUN", "TWE", "WEI"), 
                          pop.levels = c("SKY", "LIM", "TWE" , "NIN", "CNC", "MOO", "SUN", "GOO", "WEI", "FRE")
                          )


# LOCI and SNP

# sumstats
loci.sumstats <- sumstats %>%
  select (LOCUS) %>%
  distinct(LOCUS) %>%
  arrange(LOCUS)

snp.sumstats <- sumstats %>%
  select (POS) %>%
  distinct(POS) %>%
  arrange(POS)

write.table(snp.sumstats, "snp.sumstats", sep = "\t", row.names = F, col.names = T, quote = F)

# VCF
loci.vcf <- stacks.vcf.before.filters %>%
  select (LOCUS) %>%
  distinct(LOCUS) %>%
  arrange(LOCUS)

snp.vcf <- stacks.vcf.before.filters %>%
  select (POS) %>%
  distinct(POS) %>%
  arrange(POS)
write.table(snp.vcf, "snp.vcf", sep = "\t", row.names = F, col.names = T, quote = F)

# Haplotypes
loci.haplo <- haplotypes %>%
  mutate(
    CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus"),
    POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")
    ) %>%
  group_by(LOCUS) %>%
  summarise(
    CONSENSUS_MAX = max(CONSENSUS),
    POLYMORPHISM_MAX = max(POLYMORPHISM)
    ) %>%
  filter(CONSENSUS_MAX == 0 & POLYMORPHISM_MAX <= 1) %>%
  select(LOCUS) %>%
  distinct(LOCUS) %>%
  arrange(LOCUS)

paralogs <- haplotypes %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    group_by(LOCUS) %>%
    select(LOCUS) %>%
    distinct(LOCUS) %>%
    arrange(LOCUS)

consensus <- haplotypes %>%
    mutate(CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus")) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
    filter(CONSENSUS_MAX > 0) %>%
    group_by(LOCUS) %>%
    select(LOCUS) %>%
    distinct(LOCUS) %>%
    arrange(LOCUS)


  
  
  
# quick look at the caracteristics of the loci in the haplotypes files that are 
# not used in the sumstats...

loci.haplo.summary <- haplotypes %>%
  group_by(LOCUS) %>%
  select(LOCUS) %>%
  distinct(LOCUS) %>%
  arrange(LOCUS) %>%
  mutate(
    GROUP = ifelse(LOCUS %in% paralogs$LOCUS, "paralog",
            ifelse(LOCUS %in% consensus$LOCUS, "consensus", "pass")),
    SUM = ifelse(LOCUS %in% loci.sumstats$LOCUS, "sum", "absent"),
    VCF = ifelse(LOCUS %in% loci.vcf$LOCUS, "vcf", "absent")
    )

loci.sum.absent <- loci.haplo.summary %>%
  filter(SUM == "absent") %>%
  group_by(GROUP) %>%
  summarise(N = n())
  
 loci.sum.passed <- loci.haplo.summary %>%
  filter(SUM == "sum") %>%
  group_by(GROUP) %>%
  summarise(N = n()) 

### conclusion, better to use the paralog filter, 
### accepting the whitelist form a filtered haplotypes file is no guarantee...
### for the brook charr 5 loci pass the haplotype filter but are not present in
### the sumstats or vcf file. They show very low polymorphism and might be 
### filtered by stacks populations somehow before the sumstats or vcf are created.

# Apply the whitelist of haplotypes filtered loci to the VCF and sumstats and look at the proportion of missing genotype per pop...

names(sumstats)
# SNP
number.id.sumstats <- sumstats %>%
  select(LOCUS, POS, POP_ID, N) %>%
  rename(N_SUM = N) %>%
  arrange(LOCUS, POS, POP_ID)

# LOCUS
number.id.sumstats <- sumstats %>%
  select(LOCUS, POS, POP_ID, N) %>%
  group_by(LOCUS, POP_ID) %>%
  summarise(N_SUM = ceiling(mean(N, na.rm = T))) %>%
  semi_join(loci.haplo, by = "LOCUS") %>%
  arrange(LOCUS,POP_ID)

write.table(number.id.sumstats, "number.id.sumstats", sep = "\t", row.names = F, col.names = T, quote = F)
  
#   group_by(LOCUS, POS, POP_ID) %>%
#   summarise(N = mean(N))

# SNP
number.id.vcf <- stacks.vcf.before.filters %>%
  group_by(LOCUS, POS, POP_ID) %>%
  summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
  select(LOCUS, POS, POP_ID, N_VCF) %>%
  arrange(LOCUS, POS, POP_ID)

# LOCUS *** Important to use ceiling in the mean calculation, otherwise the numbers will differ between stacks format (vcf, sumstats, haplotypes) that use SNP or LOCI internally...
# the next one gives identical results with the sumstats..
number.id.vcf <- stacks.vcf.before.filters %>%
  group_by(LOCUS, POS, POP_ID) %>%
  summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
  group_by(LOCUS, POP_ID) %>%
  summarise(N_VCF = ceiling(mean(N_VCF, na.rm = T))) %>%
  select(LOCUS, POP_ID, N_VCF) %>%
  semi_join(loci.haplo, by = "LOCUS") %>%
  arrange(LOCUS, POP_ID)
write.table(number.id.vcf, "number.id.vcf", sep = "\t", row.names = F, col.names = T, quote = F)

# not this one
# what might be a difference is that with VCF you could also do this...
number.id.vcf2 <- stacks.vcf.before.filters %>%
  group_by(LOCUS, INDIVIDUALS, POP_ID) %>%
  distinct() %>%
  group_by(LOCUS, POP_ID) %>%
  summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
  select(LOCUS, POP_ID, N_VCF) %>%
  arrange(LOCUS, POP_ID)
write.table(number.id.vcf2, "number.id.vcf2", sep = "\t", row.names = F, col.names = T, quote = F)


number.id.haplo <- haplotypes %>%
  mutate(
    CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus"),
    POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")
    ) %>%
  group_by(LOCUS) %>%
  summarise(
    CONSENSUS_MAX = max(CONSENSUS),
    POLYMORPHISM_MAX = max(POLYMORPHISM)
    ) %>%
  filter(CONSENSUS_MAX == 0 & POLYMORPHISM_MAX <= 1) %>%
  left_join(haplotypes, by = "LOCUS") %>%
  group_by(LOCUS, POP_ID) %>%
  filter(HAPLOTYPES != "-") %>%
#   distinct(SAMPLES) %>%
  summarise(N_HAPLO = length(SAMPLES)) %>%
  select(LOCUS, POP_ID, N_HAPLO) %>%
  filter(LOCUS %in% number.id.sumstats$LOCUS) %>%
  arrange(LOCUS, POP_ID)

write.table(number.id.haplo, "number.id.haplo", sep = "\t", row.names = F, col.names = T, quote = F)


# compare difference between pop number
number.pop.sumstats <- number.id.sumstats %>%
  group_by(LOCUS) %>%
  summarise(POP_SUM = length(POP_ID)) %>%
  arrange(POP_SUM)
range(number.pop.sumstats$POP_SUM)
hist(number.pop.sumstats$POP_SUM)

write.table(number.pop.sumstats, "number.pop.sumstats", sep = "\t", row.names = F, col.names = T, quote = F)

number.pop.vcf <- number.id.vcf %>%
  group_by(LOCUS) %>%
  summarise(POP_VCF = length(POP_ID)) %>%
  arrange(POP_VCF)

write.table(number.pop.vcf, "number.pop.vcf", sep = "\t", row.names = F, col.names = T, quote = F)

number.pop.haplo <- number.id.haplo %>%
  group_by(LOCUS) %>%
  summarise(POP_HAPLO = length(POP_ID))
write.table(number.pop.haplo, "number.pop.haplo", sep = "\t", row.names = F, col.names = T, quote = F)

