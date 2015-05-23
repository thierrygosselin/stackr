setwd("/Users/thierry/Dropbox/brook_charr_pop/01_stacks_populations/populations_10pop")
getwd()
# data(nancycats)
# x <- as.loci(nancycats)
haplotypes.file <- "batch_1.haplotypes.tsv"
pop.id.start <- 5
pop.id.end <- 7
SITES_LEVELS <- c("SKY", "LIM", "TWE" , "NIN", "CNC", "MOO", "SUN", "GOO", "WEI", "FRE") # 10 pop
pop.levels <- SITES_LEVELS


  haplo <- read_tsv(haplotypes.file, col_names = T) %>% 
    rename(LOCUS = `Catalog ID`) %>%
    melt(
      id.vars = c("LOCUS", "Cnt"),
      variable.name = "SAMPLES",
      value.name = "HAPLOTYPES",
      factorAsStrings = F
      ) %>%
    filter(HAPLOTYPES != "consensus") %>%
    arrange(LOCUS) %>%
    mutate(
      population = str_sub(SAMPLES, pop.id.start, pop.id.end),
      population = factor(population, levels = pop.levels, ordered = T),
      POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")
      ) %>%
    filter(POLYMORPHISM <= 1) %>%
    select(LOCUS, SAMPLES, population, HAPLOTYPES) %>%
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
      id.vars = c("LOCUS", "SAMPLES", "population"),
      measure.vars = c("ALLELE_1", "ALLELE_2"), 
      variable.name = "ALLELE", 
      value.name = "NUCLEOTIDES"
      ) %>%
    group_by(LOCUS) %>%
    mutate(
      NUCLEOTIDES = factor(NUCLEOTIDES),
      NUCLEOTIDES = as.numeric(NUCLEOTIDES),
      NUCLEOTIDES = NUCLEOTIDES-1, # this removes 1 levels to enable missing values = 0
      NUCLEOTIDES = str_pad(NUCLEOTIDES, 3, side = "right", pad = "0"),
      NUCLEOTIDES = as.numeric(NUCLEOTIDES)
      ) %>%
    dcast(LOCUS + SAMPLES + population ~ ALLELE, value.var = "NUCLEOTIDES") %>%
    unite(GENOTYPE, ALLELE_1:ALLELE_2, sep = "/") %>%
    mutate(
      GENOTYPE = stri_replace_all_fixed(GENOTYPE, "0/0", "NA", vectorize_all=F)
    ) %>%
    dcast(SAMPLES + population ~ LOCUS, value.var = "GENOTYPE")

class(haplo$SAMPLES)
id <- as.character(haplo$SAMPLES)
pop <- haplo$population
class(pop)

geno <- haplo %>%
  select(-SAMPLES, -population)

test.genind <- df2genind(X = geno, sep = "/", ind.names = id, loc.names = loci, pop = pop, NA.char = "NA", ploidy = 2)
class(test.genind)



class(haplo$NUCLEOTIDES) 
class(haplo$GENOTYPE) 

test <- as.loci(x = haplo, allele.sep = "/", col.pop = haplo$population)
class(test)
hw <- hw.test(x, B = 0)
hw <- hw.test(test, B = 0)
x <- haplo



# mean heterozygosity of a locus by population
het <- read_tsv(haplotypes.file, col_names = T) %>%
    select(-Cnt) %>%
    rename(LOCUS = `Catalog ID`) %>%
    melt(
      id.vars = c("LOCUS"),
      variable.name = "SAMPLES",
      value.name = "HAPLOTYPES",
      factorAsStrings = F
      ) %>%
    filter(HAPLOTYPES != "consensus") %>%
    filter(HAPLOTYPES != "-") %>%
    arrange(LOCUS) %>%
    mutate(
      POP_ID = str_sub(SAMPLES, pop.id.start, pop.id.end),
      POP_ID = factor(POP_ID, levels = pop.levels, ordered = T),
      POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")
      ) %>%
    filter(POLYMORPHISM <= 1) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      HOM = length(POLYMORPHISM[POLYMORPHISM == 0]),
      HET = length(POLYMORPHISM[POLYMORPHISM == 1])
      ) %>%
    mutate(
      HOM_O = HOM / (HOM + HET),
      HET_O = HET / (HOM + HET)
      ) %>%
    select(-HOM, -HET) %>%
  mutate(
    PP = round(N * FREQ_REF * FREQ_REF, 0),
    PQ2 = round(2 * N * FREQ_REF * FREQ_ALT, 0),
    QQ = round(N * FREQ_ALT * FREQ_ALT, 0)
    )


hist(het$HET_O)
class(het$HET)




#### Use with adegenet
haplo <- read_tsv(haplotypes.file, col_names = T) %>% 
    rename(LOCUS = `Catalog ID`) %>%
    melt(
      id.vars = c("LOCUS", "Cnt"),
      variable.name = "SAMPLES",
      value.name = "HAPLOTYPES",
      factorAsStrings = F
      ) %>%
    filter(HAPLOTYPES != "consensus") %>%
    arrange(LOCUS) %>%
    mutate(
      population = str_sub(SAMPLES, pop.id.start, pop.id.end),
      population = factor(population, levels = pop.levels, ordered = T),
      POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")
      ) %>%
    filter(POLYMORPHISM <= 1) %>%
    select(LOCUS, SAMPLES, population, HAPLOTYPES) %>%
    mutate(
      HAPLOTYPES = stri_replace_all_fixed(HAPLOTYPES, "-", "NA", 
                                          vectorize_all=F)
      ) %>%
    dcast(SAMPLES + population ~ LOCUS, value.var = "HAPLOTYPES") %>%
  arrange(population)

class(haplo$SAMPLES)
id <- as.character(haplo$SAMPLES)
pop <- haplo$population
class(pop)

geno <- haplo %>%
  select(-SAMPLES, -population)

loci <- colnames(geno)
class(loci)

test.genind <- df2genind(X = geno, sep = "/", ind.names = id, loc.names = loci, pop = pop, NA.char = "NA")
class(test.genind)


#### testing HardyWeinberg package on the vcf file

source("/Users/thierry/Dropbox/stackr/R/read_stacks_vcf.R")
SITES_LEVELS

stacks.vcf.before.filters <- read_stacks_vcf(
  vcf.file = "batch_1.vcf",
  skip.line = 9,
  pop.id.start = 5,
  pop.id.end = 7,
  pop.levels = SITES_LEVELS,
  filter = "whitelist.loci.IND.POP.txt",
  filename = "stacks.vcf.before.filters.tsv"
  )
View(stacks.vcf.before.filters)
names(stacks.vcf.before.filters)

sample.number <- stacks.vcf.before.filters %>%
  filter(GT != "./.") %>%
  group_by(LOCUS, POS, POP_ID) %>%
  summarise(N = n())
  

prep.hwe <- stacks.vcf.before.filters %>%
  filter(GT != "./.") %>%
  group_by(LOCUS, POS, POP_ID) %>%
  summarise(
    N = n(),
    PP = as.numeric(length(GT[GT == "0/0"])),
    PQ = as.numeric(length(GT[GT == "0/1" | GT == "1/0"])),
    QQ = as.numeric(length(GT[GT == "1/1"]))
    ) %>%
  mutate(
    FREQ_P = ((PP*2) + PQ)/(2*N),
    FREQ_Q = ((QQ*2) + PQ)/(2*N),
    HET_O = PQ/N,
    HET_E = 2 * FREQ_P * FREQ_Q,
    FIS = (HET_E - HET_O)/HET_E
  )


hwe.fis <- prep.hwe %>%
  group_by(LOCUS, POS, POP_ID) %>%
  do(
  FIS=HWf(c(.$PP,.$PQ, .$QQ))
  )

hwe.exact <- prep.hwe %>%
  group_by(LOCUS, POS, POP_ID) %>%
  do(
    EXACT = HWExact(c(.$PP,.$PQ, .$QQ), alternative = "two.sided", pvaluetype = "dost")
  ) %>%
  mutate(
    PVAL = EXACT$pval,
    PSAMPLE = EXACT$pofthesample
    ) %>%
  select(-EXACT) %>%
  full_join(prep.hwe, by = c("LOCUS", "POS", "POP_ID"))

HWQqplot(X = hwe.matrix, nsim = 100, fit = "curve", logplot = T, main = "Q-Q plot for HWE", pvaluetype = "selome")




hwe.perm <- prep.hwe %>%
  group_by(LOCUS, POS, POP_ID) %>%
  do(
    PERM = HWPerm(c(.$PP,.$PQ, .$QQ), nperm = 100, verbose = T, function(y) 1 - HWExact(y)$pval)
  ) %>%
  mutate(
    STAT = PERM$stat,
    PVAL = PERM$pval
    ) %>%
  select(-PERM)


hwe <- prep.hwe %>%
  ungroup() %>%
  select(PP, PQ, QQ)

hwe.matrix <- as.matrix(hwe)
class(hwe.matrix)

results.hwe <- HWChisqMat(X = hwe.matrix)
output <- prep.hwe %>%
  cbind(results.hwe$chisqvec, results.hwe$pvalvec, results.hwe$Dvec)

n <- sum(hwe.matrix)
nN <- mac(X = hwe.matrix)

pw4 <- HWPower(n, nN, alpha = 0.05, test = "exact", theta = 4, pvaluetype = "selome")
print(pw4)
pw4 <- HWPower(n, nN, alpha = 0.05, test = "chisq", theta = 4, pvaluetype = "selome")
print(pw4)


pw8 <- HWPower(n, nN, alpha = 0.05, test = "exact", theta = 8, pvaluetype = "selome")
pw8 <- HWPower(n, nN, alpha = 0.05, test = "chisq", theta = 8, pvaluetype = "selome")

print(pw8)


#### Filter prep ###

haplotypes.file <- "batch_1.haplotypes.tsv"
pop.id.start <- 5
pop.id.end <- 7
SITES_LEVELS <- c("SKY", "LIM", "TWE" , "NIN", "CNC", "MOO", "SUN", "GOO", "WEI", "FRE") # 10 pop
pop.levels <- SITES_LEVELS

data <- stacks.vcf.before.filters
names(data)
class(data$INDIVIDUALS)



  filter(GT != "./.") %>%
  group_by(LOCUS, POS, POP_ID) %>%
  summarise(
    N = n(),
    PP = as.numeric(length(GT[GT == "0/0"])),
    PQ = as.numeric(length(GT[GT == "0/1" | GT == "1/0"])),
    QQ = as.numeric(length(GT[GT == "1/1"]))
    ) %>%
  mutate(
    FREQ_P = ((PP*2) + PQ)/(2*N),
    FREQ_Q = ((QQ*2) + PQ)/(2*N),
    HET_O = PQ/N,
    HET_E = 2 * FREQ_P * FREQ_Q,
    FIS = (HET_E - HET_O)/HET_E
  )

data <- stacks.vcf.before.filters
ind.threshold <- 10

individual_vcf_filter <- function(data, population.map, pop.id.start, pop.id.end, pop.levels, ind.threshold, threshold.fixed ) {
  
  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      } else {
        data <- data
        }

  summary <- data %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(N_TOT = n_distinct(INDIVIDUALS)) %>%
    full_join(
      data %>%
        group_by(LOCUS, POP_ID) %>%
        summarise (IND = length (GT[GT != "./."])),
      by = c("LOCUS", "POP_ID")
      ) %>%
    mutate(
      PERC = (IND / N_TOT) * 100,
      PROP = IND / N_TOT
      )
  
  if (stri_detect_fixed(threshold.fixed, "T")) {
    
    ind.filter <- summary %>%
        group_by(LOCUS) %>%
        summarise (
          POP = length (POP_ID),
          IND = length (IND[IND >= ind.threshold])
          ) %>%
        group_by(LOCUS) %>%
        filter(round(((IND/POP)*100),0) >= 100) %>%
        select(LOCUS) %>%
        left_join(data, by="LOCUS") %>%
        mutate(
          POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
          ) %>%
        arrange(POP_ID, LOCUS, POS)
  
  } else if (stri_detect_fixed(ind.threshold, ".") == "TRUE") {
    
    ind.filter <- summary %>%
      group_by(LOCUS) %>%
      summarise(
        POP = length (POP_ID),
        IND = length ([PROP >= ind.threshold])
        ) %>%
      group_by(LOCUS) %>%
      filter(round(((IND/POP)*100),0) >= 100) %>%
      select(LOCUS) %>%
      left_join(data, by="LOCUS") %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
        ) %>%
        arrange(POP_ID, LOCUS, POS)
  
  } else {
    ind.filter <- summary %>%
      group_by(LOCUS) %>%
      summarise(
        POP = length (POP_ID),
        IND = length ([PERC >= ind.threshold])
        ) %>%
      group_by(LOCUS) %>%
      filter(round(((IND/POP)*100),0) >= 100) %>%
      select(LOCUS) %>%
      left_join(data, by="LOCUS") %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =T)
        ) %>%
        arrange(POP_ID, LOCUS, POS)

}
