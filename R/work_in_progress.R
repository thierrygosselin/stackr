



min_coverage_filter <- function (data, population.map, pop.id.start, pop.id.end, pop.levels, min.depth.threshold, ind.depth.threshold, percent) {

  if (is.vector(data) == "TRUE") {
      stacks.vcf.file <- read_tsv(data, col_names = T)
    } else {
      stacks.vcf.file <- data
    }
  
  if (is.vector(population.map) == "TRUE") {
    population.map <- read_tsv(population.map, col_names = F) %>%
      select(INDIVIDUALS=X1) %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
        ) %>%
      arrange(INDIVIDUALS, POP_ID)
    } else {
    population.map <- population.map
    }  
  
  
 if(stri_detect_fixed(ind.depth.threshold, ".") & ind.depth.threshold < 1) {
    multiplication.number <- 1
    } else {
      multiplication.number <- 100
    }  
  
  number.individuals.per.pop <- population.map %>%
  group_by(POP_ID) %>%
  tally() %>%
  rename(N = n)

  filter <- stacks.vcf.file %>%
    filter(GT == "0/0" | GT == "1/1") %>%
    group_by(LOCUS, INDIVIDUALS, POP_ID) %>%
    summarise(
      MIN_REF = min(ALLELE_REF_DEPTH, na.rm = T),
      MIN_ALT = min(ALLELE_ALT_DEPTH, na.rm = T)
      ) %>%
    filter(MIN_REF < min.depth.threshold | MIN_ALT < min.depth.threshold) %>%
    group_by(POP_ID, LOCUS) %>%
    summarise(PROBLEM = n()) %>%
    full_join(number.individuals.per.pop, by = "POP_ID") %>%  
    filter(PROBLEM / N * multiplication.number > ind.depth.threshold) %>%
    group_by(LOCUS) %>%
    select(LOCUS) %>%
    distinct(LOCUS)


# filter
invisible(
  cat(
    sprintf(
"Minimum coverage filter:
1. Min depth for REF and ALT alleles: %s
2. Proportion or Percentage of individuals by pop threshold = %s \n
The number of markers before the min coverage filter: LOCI = %s
The number of markers removed by the min coverage filter: LOCI = %s
The number of markers after the min coverage filter: LOCI = %s",
min.depth.threshold,
ind.depth.threshold,
n_distinct(stacks.vcf.file$LOCUS), n_distinct(filter$LOCUS), n_distinct(stacks.vcf.file$LOCUS)-n_distinct(filter$LOCUS)
)))
filter
}


blacklist_loci_gl_vcf_venn <- function (data, ref.allele.coverage.min.threshold, read.mean.coverage.max.threshold, gl.mean.threshold, gl.min.threshold, pop.threshold, percent, filename) {

  if (is.vector(data) == "TRUE") {
      stacks.vcf.file <- read_tsv(data, col_names = T)
    } else {
      stacks.vcf.file <- data
  }
  
  pop.number <- n_distinct(stacks.vcf.file$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    } else if (stri_detect_fixed(percent, "T")) {
      multiplication.number <- 100/pop.number
      } else {
        multiplication.number <- 1
        }
  
  loci.snp.individuals <- stacks.vcf.file %>%
    group_by(LOCUS, INDIVIDUALS) %>%
    summarise(SNP_TOT = n_distinct(POS))


  filter <- stacks.vcf.file %>%
    filter(ALLELE_REF_DEPTH >= ref.allele.coverage.min.threshold) %>% # min coverage REF allele
    filter(READ_DEPTH <= read.mean.coverage.max.threshold) %>% # max read coverage to include secondary reads
    group_by(LOCUS,  POP_ID) %>%
    summarise(
      GL_MEAN = mean(GL, na.rm = T),
      GL_MIN = min(GL, na.rm = T),
      GL_MAX = max(GL, na.rm = T),
      GL_DIFF = GL_MAX - GL_MIN,
      COVERAGE_MAX = max(READ_DEPTH, na.rm = T)
      ) %>%
#   filter(GL_MEAN >= gl.mean.threshold) %>%  # mean threshold only
    filter(GL_MEAN >= gl.mean.threshold & GL_MIN >= gl.min.threshold & COVERAGE_MAX < coverage.max.threshold) %>%  # GL
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.threshold) %>%
    select(LOCUS) %>%
    left_join(stacks.vcf.file, by = "LOCUS") %>%
    arrange(LOCUS, POS, POP_ID)
    
  blacklist <- anti_join(stacks.vcf.file, filter, by = "LOCUS") %>%
    select(LOCUS, POS)
  
  blacklist.loci <- blacklist %>%
    ungroup() %>%
    select(LOCUS) %>%
    distinct(LOCUS)

  write.table(blacklist.loci, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)

  invisible(cat(sprintf(
  "Blacklisted loci by the GL filter from the original dataset, necessary for Venn diagram
  1. Min individuals REF allele coverage threshold: %s
  2. Maximum read coverage threshold: %s
  3. Mean and min loci GL threshold: %s mean GL, %s min GL with 
  population threshold = %s percent\n
  The number of SNP blacklisted = %s SNP,
  The number of LOCI blacklisted = %s LOCI\n,
  Blacklist filename:
  %s\n
  Written in the directory:
  %s",
  ref.allele.coverage.min.threshold,
  read.mean.coverage.max.threshold,
  gl.mean.threshold, gl.min.threshold,
  pop.threshold,
  n_distinct(blacklist$POS),
  n_distinct(blacklist$LOCUS),
  filename, getwd())))
blacklist.loci
}




#### HAPLOTYPES FILES ERASE TO ZERO #####
haplotype.erase <- read_tsv(haplotypes.file, col_names = T) %>% 
  rename(Catalog.ID = `Catalog ID`) %>%
  melt(
    id.vars = c("Catalog.ID", "Cnt"),
    variable.name = "SAMPLES",
    value.name = "HAPLOTYPES",
    factorAsStrings = F
    ) %>%
  mutate(
    HAPLOTYPES = ifelse(SAMPLES %in% blacklist$INDIVIDUALS & Catalog.ID %in% blacklist$LOCUS, "-", HAPLOTYPES )
    ) %>%
  dcast(Catalog.ID+Cnt~SAMPLES, value.var="LOCUS")

write.table(haplotype.erase,"haplotype.erase.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

  
  separate(
      col = HAPLOTYPES, into = c("ALLELE_1", "ALLELE_2"), 
      sep = "/", extra = "drop", remove = T
      ) %>%
    mutate(
      ALLELE_1 = stri_replace_all_fixed(ALLELE_1, "NA", "0", vectorize_all=F),
      ALLELE_2 = stri_replace_na(ALLELE_2, replacement = "no_allele"),
      ALLELE_2 = ifelse(ALLELE_2 == "no_allele", ALLELE_1, ALLELE_2)
      ) %>%
    melt(
      id.vars = c("Catalog.ID","Cnt", "SAMPLES"),
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

####  Haplotypes vcf stacks file 'batch_1.haplotypes.vcf' within R #######################################
setwd("/Users/thierry/Dropbox/esturgeon_dropbox/stacks_populations_2015/01_stacks_populations/populations_8pop")

batch_1.haplotypes.modified.vcf <- read_tsv ("batch_1.haplotypes.vcf", col_names = T, skip = 7) %>%
  rename(LOCUS = ID) %>%
  right_join(whitelist.stacks.pop, by = "LOCUS")


# Now write the output 
newfile <- file("batch_1.haplotypes.modified.vcf", "write")
cat(paste('##fileformat=VCFv4.2',
          '##fileDate=20150410',
          '##source="Stacks v1.29 modified in R"',
          '##INFO=<ID=NS,Number=1,Type=Integer, Description="Number of Samples With Data">',
          '##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">',
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
          sep="\n"), sep = "\n", file = newfile, append = TRUE)
write.table(batch_1.haplotypes.modified.vcf, file = newfile, append = TRUE, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
close(newfile)

#### COUNT HAPLOTYPE ALLELE BY LOCUS
test <- read_tsv("batch_1.haplotypes.tsv", col_names = T) %>% 
    rename(Catalog.ID = `Catalog ID`) %>%
    melt(
      id.vars = c("Catalog.ID", "Cnt"),
      variable.name = "SAMPLES",
      value.name = "HAPLOTYPES",
      factorAsStrings = F
      ) %>%
    arrange(Catalog.ID) %>%
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
    melt(
      id.vars = c("Catalog.ID","Cnt", "SAMPLES"),
      measure.vars = c("ALLELE_1", "ALLELE_2"), 
      variable.name = "ALLELE", 
      value.name = "NUCLEOTIDES"
      ) %>%
    group_by(Catalog.ID) %>%
    mutate(
      NUCLEOTIDES = factor(NUCLEOTIDES),
      NUCLEOTIDES = as.numeric(NUCLEOTIDES),
      NUCLEOTIDES = NUCLEOTIDES-1
      )  %>%  # this removes 1 levels to enable missing values = 0
    summarise(HAPLO = max(NUCLEOTIDES))




##### CoVERAGE IMBALANCE
names(stacks.vcf.before.filters)

test.imbalance <- stacks.vcf.before.filters %>%
  filter(READ_DEPTH <= 8) %>%
  filter(ALLELE_COVERAGE_RATIO != "NA")

class(stacks.vcf.before.filters$ALLELE_COVERAGE_RATIO)
range(test.imbalance$ALLELE_COVERAGE_RATIO)

hist(test.imbalance$ALLELE_COVERAGE_RATIO)

imbalance.coverage <- stacks.vcf.before.filters %>%
  mutate(
    GROUP_COVERAGE = ifelse(READ_DEPTH <= 8, "READ DEPTH <= 8", "READ DEPTH > 8"),
    GROUP_GL = ifelse(GL < 15, "POOR GL", "GOOD GL")
    ) %>%
  filter(GROUP_COVERAGE != "NA")

figure_imbalance_coverage <- function(data, pop.levels, aes.colour, adjust.bin) {

  if (is.vector(before.filter.data) == "TRUE") {
  data <- read_tsv(data, col_names = T, col_types = "diidccddccccdddddc") %>%
    mutate(INDIVIDUALS = factor(INDIVIDUALS))
  } else {
  data <- data
  }
  
  if (missing(pop.levels) == "TRUE") {
  data <- data
  } else {
  data <- data %>%
    mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
  }
  
  
  imbalance.coverage <- data %>%
    mutate(
      GROUP_COVERAGE = ifelse(READ_DEPTH <= 8, "READ DEPTH <= 8", "READ DEPTH > 8"),
      GROUP_GL = ifelse(GL < 15, "POOR GL", "GOOD GL")
      )
aes.colour <- aes(y = ..count.., colour = GROUP_COVERAGE)
adjust.bin <- 0.5
adjust.bin <- 0.9

graph <- ggplot(imbalance.coverage, aes(x = ALLELE_COVERAGE_RATIO))+
#   geom_bar()+
  geom_line(aes.colour, stat = "density", adjust = adjust.bin)+
  labs(x="Coverage imbalance between allele")+
  labs(y="Distribution (number)")+
  facet_grid(~GROUP_COVERAGE)+
  theme(axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
        legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
        legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
        strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold"))

  graph
}

