

# min coverage filter
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

# Genotype likelihood
gl_vcf_filter <- function (data, read.depth.max.threshold, gl.mean.threshold, gl.min.threshold, pop.threshold, percent) {

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
    
#   loci.snp.individuals <- stacks.vcf.file %>%
#     group_by(LOCUS, INDIVIDUALS) %>%
#     summarise(SNP_TOT = n_distinct(POS))

   filter <- stacks.vcf.file %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      READ_DEPTH_MAX = max(READ_DEPTH, na.rm = T),
      GL_MEAN = mean(GL, na.rm = T),
      GL_MIN = min(GL, na.rm = T),
      GL_MAX = max(GL, na.rm = T),
      GL_DIFF = GL_MAX - GL_MIN
      ) %>%
    filter(GL_MEAN >= gl.mean.threshold & GL_MIN >= gl.min.threshold) %>%  # GL
    filter(READ_DEPTH_MAX <= read.depth.max.threshold) %>% # read coverage
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.threshold) %>%
    select(LOCUS) %>%
    left_join(stacks.vcf.file, by = "LOCUS") %>%
    arrange(LOCUS, POS, POP_ID)
  

# filter
invisible(
  cat(
    sprintf(
"Genotype likelihood (GL) filter:
1. Max read depth threshold: %s
2. Mean and min loci GL threshold: %s mean GL, %s min GL with 
population threshold = %s percent\n
The number of markers before the GL filter: SNP = %s, LOCI = %s
The number of markers removed by the GL filter: SNP = %s, LOCI = %s
The number of markers after the GL filter: SNP = %s, LOCI = %s",
read.depth.max.threshold,
gl.mean.threshold, gl.min.threshold,
pop.threshold,
n_distinct(stacks.vcf.file$POS), n_distinct(stacks.vcf.file$LOCUS), (n_distinct(stacks.vcf.file$POS)-n_distinct(filter$POS)), (n_distinct(stacks.vcf.file$LOCUS)-n_distinct(filter$LOCUS)), n_distinct(filter$POS), n_distinct(filter$LOCUS)
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


# Minor Allele Frequency
maf_filter <- function(data, maf.diff.threshold, pop.maf.diff.threshold, percent, local.maf.threshold, global.maf.threshold, pop.maf.threshold) {

  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
    } else {
      data <- data
  }
  
  pop.number <- n_distinct(data$POP_ID)

  if (stri_detect_fixed(pop.maf.diff.threshold, ".") & pop.maf.diff.threshold < 
    1) {
    multiplication.number <- 1/pop.number
  } else if (stri_detect_fixed(percent, "T")) {
    multiplication.number <- 100/pop.number
  } else {
    multiplication.number <- 1
  }

  maf.filter <- data %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(MAF_DIFF = max(FREQ_ALT) - min(FREQ_ALT)) %>% 
    filter(MAF_DIFF <= maf.diff.threshold) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.maf.diff.threshold) %>%
    select(LOCUS) %>%
    left_join(data, by="LOCUS") %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      GLOBAL_MAF = mean(GLOBAL_MAF),  
      LOCAL_MAF = min(FREQ_ALT)
      ) %>% 
    filter(LOCAL_MAF >= local.maf.threshold | GLOBAL_MAF >= global.maf.threshold) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter(n >= pop.maf.threshold) %>% 
    select(LOCUS) %>%
    left_join(data, by = "LOCUS") %>%
    arrange(LOCUS, POP_ID)
  
  invisible(
        cat(
         sprintf(
"MAF filter v.1.0:
1. To pass the filter, at least %s percent of sites/pop need to 
have a difference between the max and min MAF along the read/haplotype <= %s 
2. To pass the filter, markers with a local MAF >= %s or a global MAF >= %s, 
in at least %s sampling sites/pop.
The number of SNP removed by the MAF filter = %s SNP
The number of LOCI removed by the MAF filter = %s LOCI
The number of SNP after the MAF filter = %s SNP
The number of LOCI after the MAF filter = %s LOCI", 
pop.maf.diff.threshold,
maf.diff.threshold,
local.maf.threshold,
global.maf.threshold,
pop.maf.threshold,
n_distinct(data$POS)-n_distinct(maf.filter$POS),
n_distinct(data$LOCUS)-n_distinct(maf.filter$LOCUS),
n_distinct(maf.filter$POS),
n_distinct(maf.filter$LOCUS)
        )
        )
      )
maf.filter
}

blacklist_loci_maf_venn <- function(data, maf.diff.threshold, pop.maf.diff.threshold, percent, local.maf.threshold, global.maf.threshold, pop.maf.threshold, filename) {

  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
    } else {
      data <- data
  }
  
  pop.number <- n_distinct(data$POP_ID)

  if (stri_detect_fixed(pop.maf.diff.threshold, ".") & pop.maf.diff.threshold < 
    1) {
    multiplication.number <- 1/pop.number
  } else if (stri_detect_fixed(percent, "T")) {
    multiplication.number <- 100/pop.number
  } else {
    multiplication.number <- 1
  }
  
  maf.filter <- data %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(MAF_DIFF = max(FREQ_ALT) - min(FREQ_ALT)) %>% 
    filter(MAF_DIFF <= maf.diff.threshold) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.maf.diff.threshold) %>%
    select(LOCUS) %>%
    left_join(data, by="LOCUS") %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      GLOBAL_MAF = mean(GLOBAL_MAF),  
      LOCAL_MAF = min(FREQ_ALT)
      ) %>% 
    filter(LOCAL_MAF >= local.maf.threshold | GLOBAL_MAF >= global.maf.threshold) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter(n >= pop.maf.threshold) %>% 
    ungroup() %>%
    select(LOCUS) %>%
    distinct(LOCUS)
    
  blacklist <- anti_join(data, maf.filter, by = "LOCUS") %>%
  select(LOCUS, POS)
  
  blacklist.loci <- blacklist %>%
    ungroup() %>%
    select(LOCUS) %>%
    distinct(LOCUS)
  
  write.table(blacklist.loci, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)
  
  invisible(
        cat(
         sprintf(
"Blacklisted loci by the MAF filter from the original dataset, necessary for Venn diagram
MAF filter v.1.0:
1. To pass the filter, at least %s percent of sites/pop need to 
have a difference between the max and min MAF along the read/haplotype <= %s 
2. To pass the filter, markers with a local MAF >= %s or a global MAF >= %s, 
in at least %s sampling sites/pop.
The number of SNP blacklisted = %s SNP
The number of LOCI blacklisted = %s LOCI\n
Blacklist filename:
%s\n
Written in the directory:
%s",
pop.maf.diff.threshold,
maf.diff.threshold,
local.maf.threshold,
global.maf.threshold,
pop.maf.threshold,
n_distinct(blacklist$POS),
n_distinct(blacklist$LOCUS),
filename, getwd()
        )
        )
      )
blacklist.loci
}

# Observed heterozygosity 
het_filter <- function(data, het.threshold, het.diff.threshold, pop.threshold, percent) {
  
  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      } else {
        data <- data
        }
  
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    } else if (stri_detect_fixed(percent, "T")) {
      multiplication.number <- 100/pop.number
      } else {
        multiplication.number <- 1
        }

  het.filter <- data %>%
    group_by (LOCUS, POP_ID) %>%
    summarise(
      HET_DIFF = max(HET_O) - min(HET_O),
      HET_MAX = max(HET_O)
      ) %>%
    filter(HET_DIFF <= het.diff.threshold & HET_MAX <= het.threshold) %>%  
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.threshold) %>%
    select(LOCUS) %>%
    left_join(data, by="LOCUS") %>%
    arrange(LOCUS, POP_ID)
  
  invisible(
        cat(
         sprintf(
"HET filter:
Markers with an observed heterozygosity <= %s and with a difference of <= %s
between the max and min observed heterozygosity along the reads
(independent of the number of SNP/reads)
in at least %s percent of sites/pop are kept\n
The number of SNP removed by the HET filter = %s SNP
The number of LOCI removed by the HET filter = %s LOCI
The number of SNP after the HET filter = %s SNP
The number of LOCI after the HET filter = %s LOCI", 
het.threshold, het.diff.threshold, pop.threshold,
n_distinct(data$POS)-n_distinct(het.filter$POS),
n_distinct(data$LOCUS)-n_distinct(het.filter$LOCUS),
n_distinct(het.filter$POS),
n_distinct(het.filter$LOCUS)
        )
        )
      )
het.filter
}

blacklist_loci_het_venn <- function(data, het.threshold, het.diff.threshold, pop.threshold, percent, filename) {
  
  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      } else {
        data <- data
        }
  
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    } else if (stri_detect_fixed(percent, "T")) {
      multiplication.number <- 100/pop.number
      } else {
        multiplication.number <- 1
        }

  het.filter <- data %>%
    group_by (LOCUS, POP_ID) %>%
    summarise(
      HET_DIFF = max(HET_O) - min(HET_O),
      HET_MAX = max(HET_O)
      ) %>%
    filter(HET_DIFF <= het.diff.threshold & HET_MAX <= het.threshold) %>%  
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.threshold) %>%
    select(LOCUS) %>%
    distinct(LOCUS)
  
  blacklist <- anti_join(data, het.filter, by = "LOCUS") %>%
    select(LOCUS, POS)
  
  blacklist.loci <- blacklist %>%
    ungroup() %>%
    select(LOCUS) %>%
    distinct(LOCUS)

  write.table(blacklist.loci, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)
  
  invisible(
        cat(
         sprintf(
"Blacklisted loci by the HET filter from the original dataset, necessary for Venn diagram:
Markers with an observed heterozygosity <= %s and with a difference of <= %s
between the max and min observed heterozygosity along the reads
(independent of the number of SNP/reads)
in at least %s percent of sites/pop are kept\n
The number of SNP blacklisted = %s SNP
The number of LOCI blacklisted = %s LOCI\n
Blacklist filename:
%s\n
Written in the directory:
%s", 
het.threshold, het.diff.threshold, pop.threshold,
n_distinct(blacklist$POS),
n_distinct(blacklist$LOCUS),
filename, getwd()
        )
        )
      )
blacklist.loci
}

# Fis
fis_filter <- function(data, fis.min.threshold = -1, fis.max.threshold = 1, fis.diff.threshold = 1, pop.threshold, percent) {
  
  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      } else {
        data <- data
        }
  
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    } else if (stri_detect_fixed(percent, "T")) {
      multiplication.number <- 100/pop.number
      } else {
        multiplication.number <- 1
        }

  fis.filter <- data %>%
    group_by (LOCUS, POP_ID) %>%
    summarise(
      FIS_MEAN = mean(FIS),
      FIS_MIN = min(FIS),
      FIS_MAX = max(FIS),
      FIS_DIFF = FIS_MAX-FIS_MIN
      ) %>%
    filter(FIS_MIN >= fis.min.threshold) %>%
    filter(FIS_MAX <= fis.max.threshold) %>%
    filter(FIS_DIFF <= fis.diff.threshold) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.threshold) %>%
    select(LOCUS) %>%
    left_join(data, by="LOCUS") %>%
    arrange(LOCUS, POP_ID)

#   invisible(
#         cat(
#          sprintf(
# "Fis filter:
# markers along the read/haplotype with difference
# between their max and min values of FIS > %s in %s percent 
# of the sampling sites/pop were removed\n
# The number of SNP removed by the Fis filter = %s SNP
# The number of LOCI removed by the Fis filter = %s LOCI
# The number of SNP after the Fis filter = %s SNP
# The number of LOCI after the Fis filter = %s LOCI", 
# fis.max.threshold, pop.threshold,
# n_distinct(data$POS)-n_distinct(fis.filter$POS),
# n_distinct(data$LOCUS)-n_distinct(fis.filter$LOCUS),
# n_distinct(fis.filter$POS),
# n_distinct(fis.filter$LOCUS)
#         )
#         )
#       )
  invisible(
        cat(
         sprintf(
"Fis filter:
Fis min >= %s or Fis max <= %s or 
difference along the read/haplotype between the max and min Fis > %s,
all in %s percent of the sampling sites/pop were removed\n
The number of SNP removed by the Fis filter = %s SNP
The number of LOCI removed by the Fis filter = %s LOCI
The number of SNP after the Fis filter = %s SNP
The number of LOCI after the Fis filter = %s LOCI", 
fis.min.threshold, fis.max.threshold, fis.diff.threshold, pop.threshold,
n_distinct(data$POS)-n_distinct(fis.filter$POS),
n_distinct(data$LOCUS)-n_distinct(fis.filter$LOCUS),
n_distinct(fis.filter$POS),
n_distinct(fis.filter$LOCUS)
        )
        )
      )
fis.filter

}

blacklist_loci_fis_venn <- function(data, fis.min.threshold = -1, fis.max.threshold = 1, fis.diff.threshold = 1, pop.threshold, percent, filename) {
  
  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      } else {
        data <- data
        }
  
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    } else if (stri_detect_fixed(percent, "T")) {
      multiplication.number <- 100/pop.number
      } else {
        multiplication.number <- 1
        }

  fis.filter <- data %>%
    group_by (LOCUS, POP_ID) %>%
    summarise(
      FIS_MEAN = mean(FIS),
      FIS_MIN = min(FIS),
      FIS_MAX = max(FIS),
      FIS_DIFF = FIS_MAX-FIS_MIN
      ) %>%
    filter(FIS_DIFF <= fis.diff.threshold) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter(round(((n / pop.number) * 100), 0) >= pop.threshold) %>%
    select(LOCUS) %>%
    distinct(LOCUS)
  
  blacklist <- anti_join(data, fis.filter, by = "LOCUS") %>%
    select(LOCUS, POS)
  
  blacklist.loci <- blacklist %>%
    ungroup() %>%
    select(LOCUS) %>%
    distinct(LOCUS)

  write.table(blacklist.loci, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)

  
  invisible(
        cat(
         sprintf(
"Blacklisted loci by the Fis filter from the original dataset, necessary for Venn diagram:
Fis min >= $s or Fis max <= $s or 
difference along the read/haplotype between the max and min Fis > %s,
all in %s percent of the sampling sites/pop were removed\n
The number of SNP blacklisted = %s SNP
The number of LOCI blacklisted = %s LOCI\n
Blacklist filename:
%s\n
Written in the directory:
%s", 
fis.min.threshold, fis.max.threshold, fis.diff.threshold, pop.threshold,
n_distinct(blacklist$POS),
n_distinct(blacklist$LOCUS),
filename, getwd()
        )
        )
      )
blacklist.loci
}


# SNP number per haplotype
snp_number_filter <- function(data, max.snp.number, pop.threshold, percent) {
  
  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      } else {
        data <- data
        }
  
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    } else if (stri_detect_fixed(percent, "T")) {
      multiplication.number <- 100/pop.number
      } else {
        multiplication.number <- 1
        }

  snp.filter <- data %>%
    group_by (LOCUS, POP_ID) %>%
    summarise(SNP_N = n_distinct(POS)) %>%
    filter(SNP_N <= max.snp.number) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.threshold) %>%
    select(LOCUS) %>%
    left_join(data, by="LOCUS") %>%
    arrange(LOCUS, POP_ID)

  invisible(
        cat(
         sprintf(
"SNP number per haplotype filter:
max number of SNP allowed = %s in %s percent 
of the sampling sites/pop\n
The number of SNP removed by the SNP number filter = %s SNP
The number of LOCI removed by the SNP number filter = %s LOCI
The number of SNP after the SNP number filter = %s SNP
The number of LOCI after the SNP number filter = %s LOCI", 
max.snp.number, pop.threshold,
n_distinct(data$POS)-n_distinct(snp.filter$POS),
n_distinct(data$LOCUS)-n_distinct(snp.filter$LOCUS),
n_distinct(snp.filter$POS),
n_distinct(snp.filter$LOCUS)
        )
        )
      )
snp.filter

}

blacklist_loci_snp_number_venn <- function(data, max.snp.number, pop.threshold, percent, filename) {
  
  if (is.vector(data) == "TRUE") {
      data <- read_tsv(data, col_names = T)
      } else {
        data <- data
        }
  
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    } else if (stri_detect_fixed(percent, "T")) {
      multiplication.number <- 100/pop.number
      } else {
        multiplication.number <- 1
        }

  snp.filter <- data %>%
    group_by (LOCUS, POP_ID) %>%
    summarise(SNP_N = n_distinct(POS)) %>%
    filter(SNP_N <= max.snp.number) %>%
    group_by(LOCUS) %>%
    tally() %>%
    filter((n * multiplication.number) >= pop.threshold) %>%
    select(LOCUS) %>%
    distinct(LOCUS)
  
  blacklist <- anti_join(data, snp.filter, by = "LOCUS") %>%
    select(LOCUS, POS)
  
  blacklist.loci <- blacklist %>%
    ungroup() %>%
    select(LOCUS) %>%
    distinct(LOCUS)

  write.table(blacklist.loci, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)

  
  invisible(
        cat(
         sprintf(
"Blacklisted loci by the SNP numer filter from the original dataset, necessary for Venn diagram:
max number of SNP allowed = %s in %s percent 
of the sampling sites/pop\n
The number of SNP blacklisted = %s SNP
The number of LOCI blacklisted = %s LOCI\n
Blacklist filename:
%s\n
Written in the directory:
%s", 
max.snp.number, pop.threshold,
n_distinct(blacklist$POS),
n_distinct(blacklist$LOCUS),
filename, getwd()
        )
        )
      )
blacklist.loci
}

#' @title blacklist_loci_filter
#' @description This function use a blacklist of loci to filter a dataset of class sumstats.
#' @param data A data frame object or file (using ".tsv") of class sumstats.
#' @param blacklist The blacklist is a data frame or file (using ".tsv") with one column named `LOCUS`.
#' @param type A description of the type of blacklist filtering used.

blacklist_loci_filter <- function (data, blacklist, type) {
  
  if (is.vector(blacklist) == "TRUE") {
      blacklist.loci <- read_tsv(blacklist, col_names = T)
    } else {
      blacklist.loci <- blacklist
      }
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
    } else {
      data <- data
      }
  
  blacklist.filter <- data %>%
  anti_join(blacklist.loci, by = "LOCUS") %>%
  arrange(LOCUS, POS, POP_ID)
  
invisible(
  cat(
    sprintf(
"The number of markers before the blacklist %s filter: SNP = %s, LOCI = %s
The number of markers removed by the blacklist %s filter: SNP = %s, LOCI = %s
The number of markers after the blacklist %s filter: SNP = %s, LOCI = %s",
type, n_distinct(data$POS), n_distinct(data$LOCUS), type, (n_distinct(data$POS)-n_distinct(blacklist.filter$POS)), (n_distinct(data$LOCUS)-n_distinct(blacklist.filter$LOCUS)), type, n_distinct(blacklist.filter$POS), n_distinct(blacklist.filter$LOCUS)
)))
blacklist.filter
}


# Filter with a whitelist of individuals
#' @title whitelist_id_filter
#' @description This function use a whitelist of individuals to filter a dataset.
#' @param data A data frame object or file (using ".tsv") with an 'INDIVIDUALS'
#' column.
#' @param whitelist.id The whitelist of individuals with no header column,
#' the first column is assumed to be the indivuals ID.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels A character string with your populations ordered.
#' @param filename The name of the file written to the directory.

whitelist_id_filter <- function (data, whitelist.id, pop.id.start, pop.id.end, pop.levels, filename) {
if (is.vector(whitelist.id) == "TRUE") {
      whitelist.id <- read_tsv(whitelist.id, col_names = F) %>%
        select(INDIVIDUALS=X1) %>%
        mutate(INDIVIDUALS = factor(INDIVIDUALS))
      
    } else {
      whitelist.id <- whitelist.id %>%
        select(INDIVIDUALS) %>%
        mutate(INDIVIDUALS = factor(INDIVIDUALS)) %>%
        ungroup()
        
      }
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
#     data <- read_tsv(data, col_names = T) %>%
#       mutate(INDIVIDUALS = as.character(INDIVIDUALS))
    
    } else {
      data <- data %>%
#         mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
        ungroup()
    }
  
  whitelist.id.filter <- whitelist.id %>%
  left_join(data, by = "INDIVIDUALS") %>%
  mutate(
    POP_ID = factor(str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
                    levels = pop.levels, ordered =T)
    ) %>%
  arrange(LOCUS)
  
  write.table(whitelist.id.filter, filename, sep = "\t", row.names = F,
              col.names = T, quote = F)

invisible(cat(sprintf(
"Whitelist individuals filter 
Filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
whitelist.id.filter
}



