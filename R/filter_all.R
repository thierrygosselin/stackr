#' @title Combination of all filters
#' @description You need to get a number of markers rapidly, and you tested
#' individually each filter before, you know what you're doing... this is the
#' all-in-one solution filter.
#' @param haplotypes Working on it.
#' @param vcf Working on it.
#' @param pop.id.start Working on it.
#' @param pop.id.end Working on it. 
#' @param pop.levels Working on it.
#' @param read.depth.threshold Working on it.
#' @param allele.depth.threshold Working on it.
#' @param allele.imbalance.threshold Working on it.
#' @param allele.min.depth.threshold Working on it.
#' @param read.depth.max.threshold Working on it.
#' @param gl.mean.threshold Working on it.
#' @param gl.min.threshold Working on it.
#' @param gl.diff.threshold Working on it.
#' @param gl.pop.threshold Working on it.
#' @param gl.percent Working on it.
#' @param ind.threshold Working on it.
#' @param threshold.fixed Working on it.
#' @param pop.threshold Working on it.
#' @param pop.percent Working on it.
#' @param local.maf.threshold Working on it.
#' @param global.maf.threshold Working on it.
#' @param maf.pop.threshold Working on it.
#' @param het.threshold Working on it.
#' @param het.diff.threshold Working on it.
#' @param het.pop.threshold Working on it.
#' @param het.percent Working on it.
#' @param fis.min.threshold Working on it.
#' @param fis.max.threshold Working on it.
#' @param fis.diff.threshold Working on it.
#' @param fis.pop.threshold Working on it.
#' @param fis.percent Working on it.
#' @param max.snp.number Working on it.
#' @rdname filter_all
#' @export
#' @import stringi
#' @import dplyr
#' @import readr

filter_all <- function (haplotypes, vcf,
                        pop.id.start, pop.id.end, pop.levels,
                        read.depth.threshold, allele.depth.threshold, 
                        allele.imbalance.threshold,
                        allele.min.depth.threshold, read.depth.max.threshold, 
                        gl.mean.threshold, gl.min.threshold, gl.diff.threshold, 
                        gl.pop.threshold, gl.percent,
                        ind.threshold, threshold.fixed,
                        pop.threshold, pop.percent,
                        local.maf.threshold, global.maf.threshold, maf.pop.threshold,
                        het.threshold, het.diff.threshold, het.pop.threshold, het.percent,
                        fis.min.threshold, fis.max.threshold, fis.diff.threshold, fis.pop.threshold, fis.percent,
                        max.snp.number) {
  
  `#CHROM` <- NULL
  AF <- NULL
  ALLELE_ALT_DEPTH <- NULL
  ALLELE_COVERAGE_RATIO <- NULL
  ALLELE_DEPTH <- NULL
  # ALLELE_P <- NULL
  # ALLELE_Q <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALT <- NULL
  ALT_FREQ <- NULL
  FILTER <- NULL
  FIS <- NULL
  FIS_MAX <- NULL
  FIS_MIN <- NULL
  FIS_DIFF <- NULL
  FORMAT <- NULL
  FREQ_ALT <- NULL
  FREQ_REF <- NULL
  GL <- NULL
  GL_DIFF <- NULL
  GL_MAX <- NULL
  GL_MEAN <- NULL
  GL_MIN <- NULL
  GLOBAL_MAF <- NULL
  GROUP <- NULL
  GT <- NULL
  HET <- NULL
  HET_DIFF <- NULL
  HET_E <- NULL
  HET_MAX <- NULL
  HET_O <- NULL
  HOM_E <- NULL
  HOM_O <- NULL
  ID <- NULL
  IND <- NULL
  INDIVIDUALS <- NULL
  INFO <- NULL
  MAF <- NULL
  LOCAL_MAF <- NULL
  LOCUS_N <- NULL
  MIN_ALT <- NULL
  MIN_REF <- NULL
  N <- NULL
  N_MIN <- NULL
  N_SUM <- NULL
  N_TOT <- NULL
  N_VCF <- NULL
  POP <- NULL
  POP_ID <- NULL
  PP <- NULL
  PQ <- NULL
  QUAL <- NULL
  READ_DEPTH <- NULL
  READ_DEPTH_MAX <- NULL
  REF <- NULL
  REF_FREQ <- NULL
  STATUS <- NULL
  VALUE <- NULL
  X1 <- NULL
  POLYMORPHISM <- NULL
  POLYMORPHISM_MAX <- NULL
  QQ <- NULL
  SNP_N <- NULL
  
  
  # Haplotypes file
  message("Importing the haplotypes file...")
  
  haplo <- read_tsv(haplotypes, col_names = T) %>%
    rename(LOCUS =`Catalog ID`) %>%
    tidyr::gather(SAMPLES, HAPLOTYPES, -c(LOCUS, Cnt)) %>%
    mutate(
      POP_ID = str_sub(SAMPLES, pop.id.start, pop.id.end),
      POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
    )
  
  
  # Locus with > 2 alleles by pop
  # Create a blacklist of catalog loci with paralogs
  
  message("Looking for paralogs...")
  
  paralogs <- haplo %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    group_by(LOCUS) %>%
    select (LOCUS) %>%
    distinct(LOCUS)
  
  nparalogs <- stri_join("Found", n_distinct(paralogs$LOCUS), "paralogs", sep = " ")
  message(nparalogs)
  
  message("Importing the vcf file...")

  vcf.before.filters <- read_delim(
    vcf, delim = "\t", 
    skip = 9,
    progress = interactive()
  ) %>% 
    select(-c(QUAL, FILTER, FORMAT)) %>%
    rename(LOCUS = ID, CHROM = `#CHROM`)
  
  message("Removing paralogs from the VCF")
  
  vcf.paralogs <- vcf.before.filters %>% 
    anti_join(paralogs, by = "LOCUS")
  
  
  message("Tidying the VCF...")
  
 vcf <- vcf.paralogs %>%
    tidyr::separate(INFO, c("N", "AF"), sep = ";", extra = "warn") %>%
    mutate(
      N = as.numeric(stri_replace_all_fixed(N, "NS=", "", vectorize_all=F)),
      AF = stri_replace_all_fixed(AF, "AF=", "", vectorize_all=F)
    ) %>%
    tidyr::separate(AF, c("REF_FREQ", "ALT_FREQ"), sep = ",", extra = "warn") %>%
    mutate(
      REF_FREQ = as.numeric(REF_FREQ),
      ALT_FREQ = as.numeric(ALT_FREQ)
    )
  # Gather individuals in 1 colummn --------------------------------------------
  vcf <- tidyr::gather(vcf, INDIVIDUALS, FORMAT, -c(CHROM:ALT_FREQ))
  
  message("Gathering individuals in 1 column")
  
  # Separate FORMAT and COVERAGE columns ---------------------------------------
  message("Tidying the VCF...")
  
  vcf <- vcf %>%
    tidyr::separate(FORMAT, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"),
             sep = ":", extra = "warn") %>%
    tidyr::separate(ALLELE_DEPTH, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
             sep = ",", extra = "warn")
  
  # Work with Mutate on CHROM and GL -------------------------------------------
  message("Fixing columns...")
  
  vcf <- vcf %>%
    mutate(
      CHROM = suppressWarnings(as.numeric(stri_replace_all_fixed(CHROM, "un", "1", vectorize_all=F))),
      GL = suppressWarnings(as.numeric(stri_replace_all_fixed(GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all=F)))
    ) %>%
    # Work with Mutate on ALLELE_P and ALLELE_Q
    # vcf <- vcf %>% 
    #   mutate(
    #     ALLELE_P = ifelse (ALLELE_P == ".", "NA",
    #                        ifelse(ALLELE_P == "0", REF, ALT)),
    #     ALLELE_Q = ifelse (ALLELE_Q == ".", "NA",
    #                        ifelse(ALLELE_Q == "0", REF, ALT))
    #   )
    # Mutate read and alleles REF/ALT depth
    mutate(READ_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(READ_DEPTH, "^0$", "NA", vectorize_all=F))),
           ALLELE_REF_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(ALLELE_REF_DEPTH, "^0$", "NA", vectorize_all = TRUE))),
           ALLELE_ALT_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(ALLELE_ALT_DEPTH, "^0$", "NA", vectorize_all = TRUE))),
           # Mutate coverage ration for allelic imbalance
           ALLELE_COVERAGE_RATIO = suppressWarnings(
             as.numeric(ifelse(GT == "./." | GT == "0/0" | GT == "1/1", "NA",
                               ((ALLELE_ALT_DEPTH - ALLELE_REF_DEPTH)/(ALLELE_ALT_DEPTH + ALLELE_REF_DEPTH)))))
    )
  # Populations levels ---------------------------------------------------------
  message("Adding a population column ...")
  
  vcf.tidy <- mutate(
    vcf, 
    POP_ID = factor(str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
                    levels = pop.levels, ordered =T)
  ) %>%
    arrange(LOCUS, POS, POP_ID, INDIVIDUALS)  
    
    
    
  message("Inspecting tidy VCF for coverage problems...")
  
  
  blacklist <- vcf.tidy %>%
    filter(GT != "./." & GT != "0/0" & GT != "1/1" ) %>%
    filter(READ_DEPTH == read.depth.threshold) %>%
    filter(ALLELE_REF_DEPTH < allele.depth.threshold | ALLELE_ALT_DEPTH < allele.depth.threshold) %>%
    filter(ALLELE_COVERAGE_RATIO < -allele.imbalance.threshold | ALLELE_COVERAGE_RATIO > allele.imbalance.threshold) %>%
    select(LOCUS, POS, POP_ID, INDIVIDUALS)
  
  # interesting stats.
  erased.genotype.number <- length(blacklist$INDIVIDUALS)
  total.genotype.number <- length(vcf.tidy$GT[vcf.tidy$GT != "./."])
  
  
  erasing <- stri_join("Erasing", round(((erased.genotype.number/total.genotype.number)*100), 2), "% of", total.genotype.number, "genotypes...", sep = " ")
  message(erasing)
  
  
  vcf.erased <- suppressWarnings(vcf.tidy %>%
    mutate(
      GT = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "./.", GT),
      READ_DEPTH = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", READ_DEPTH)),
      ALLELE_REF_DEPTH = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_REF_DEPTH)),
      ALLELE_ALT_DEPTH = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_ALT_DEPTH)),
      ALLELE_COVERAGE_RATIO = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_COVERAGE_RATIO)),
      GL = as.numeric(ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", GL))
      # ALLELE_P = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_P),
      # ALLELE_Q = ifelse(LOCUS %in% blacklist$LOCUS & POS %in% blacklist$POS & INDIVIDUALS %in% blacklist$INDIVIDUALS, "NA", ALLELE_Q)
    )
  )
  
  
  message("Filter: Genotype likelihood ...")
  
  
  pop.number <- n_distinct(vcf.tidy$POP_ID)
  
  if(stri_detect_fixed(gl.pop.threshold, ".") & gl.pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message("Using a proportion threshold...")
    threshold.id <- "of proportion"
  } else if (stri_detect_fixed(gl.percent, "T")) {
    multiplication.number <- 100/pop.number
    message("Using a percentage threshold...")
    threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message("Using a fixed threshold...")
    threshold.id <- "population as a fixed"
    
  }
  
  
  gl.filter <- vcf.erased %>%
    group_by(LOCUS, POP_ID) %>% # at the population level
    summarise(
      MIN_REF = min(ALLELE_REF_DEPTH, na.rm = T),
      MIN_ALT = min(ALLELE_ALT_DEPTH, na.rm = T),
      READ_DEPTH_MAX = max(READ_DEPTH, na.rm = T),
      GL_MEAN = mean(GL, na.rm = T),
      GL_MIN = min(GL, na.rm = T),
      GL_MAX = max(GL, na.rm = T),
      GL_DIFF = GL_MAX - GL_MIN
    ) %>%
    filter(MIN_REF >= allele.min.depth.threshold | MIN_ALT >= allele.min.depth.threshold) %>%
    filter(READ_DEPTH_MAX <= read.depth.max.threshold) %>% # read coverage
    filter(GL_MEAN >= gl.mean.threshold & GL_MIN >= gl.min.threshold) %>%  # GL
    filter(GL_DIFF <= gl.diff.threshold) %>%
    group_by(LOCUS) %>%
    tally() %>% # Globally accross loci
    filter((n * multiplication.number) >= gl.pop.threshold) %>%
    select(LOCUS) %>%
    left_join(vcf.erased, by = "LOCUS")
  
  message("Filter: individuals...")
  
  if (stri_detect_fixed(threshold.fixed, "T")) {
    message("Using a fixed threshold")
    threshold.id <- "individuals as a fixed"
    
    ind.filter <- gl.filter %>%
      group_by(LOCUS, POS, POP_ID) %>%
      summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
      group_by(LOCUS, POP_ID) %>%
      summarise(N_VCF = ceiling(mean(N_VCF, na.rm = TRUE))) %>%
      group_by(LOCUS) %>%
      summarise (
        POP = length (POP_ID),
        IND = length (N_VCF[N_VCF >= ind.threshold])
      ) %>%
      group_by(LOCUS) %>%
      filter(round(((IND/POP)*100),0) >= 100) %>%
      select(LOCUS) %>%
      left_join(gl.filter, by="LOCUS") %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)
      )
    
  } else if (stri_detect_fixed(ind.threshold, ".") == "TRUE") {
    
    message("Using a proportion threshold")
    threshold.id <- " of proportion"
    
    ind.filter <- gl.filter %>%
      group_by(POP_ID) %>%
      summarise(N_TOT = as.numeric(n_distinct(INDIVIDUALS))) %>%
      mutate(N_MIN = floor(N_TOT * as.numeric(ind.threshold))) %>%
      full_join(
        data %>%
          group_by(LOCUS, POS, POP_ID) %>%
          summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
          group_by(LOCUS, POP_ID) %>%
          summarise(N_VCF = ceiling(mean(N_VCF, na.rm = TRUE))),
        by = c("POP_ID")
      )%>%
      group_by(LOCUS) %>%
      summarise(
        POP = as.numeric(length (POP_ID)),
        IND = as.numeric(length (N_VCF[N_VCF >= N_MIN]))
      ) %>%
      group_by(LOCUS) %>%
      filter(round(((IND/POP)*100),0) >= 100) %>%
      select(LOCUS) %>%
      left_join(gl.filter, by="LOCUS") %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)
      )
    
  } else {
    
    message("Using a percentage threshold")
    threshold.id <- "percent"
    
    
    ind.filter <- gl.filter %>%
      group_by(POP_ID) %>%
      summarise(N_TOT = as.numeric(n_distinct(INDIVIDUALS))) %>%
      mutate(N_MIN = floor(N_TOT * as.numeric(ind.threshold / 100))) %>%
      full_join(
        gl.filter %>%
          group_by(LOCUS, POS, POP_ID) %>%
          summarise(N_VCF = length(INDIVIDUALS[GT != "./."])) %>%
          group_by(LOCUS, POP_ID) %>%
          summarise(N_VCF = ceiling(mean(N_VCF, na.rm = T))),
        by = c("POP_ID")
      )%>%
      group_by(LOCUS) %>%
      summarise(
        POP = as.numeric(length (POP_ID)),
        IND = as.numeric(length (N_VCF[N_VCF >= N_MIN]))
      ) %>%
      group_by(LOCUS) %>%
      filter(round(((IND/POP)*100),0) >= 100) %>%
      select(LOCUS) %>%
      left_join(gl.filter, by="LOCUS") %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)
      )  
    
  }
  message("Filter: populations...")
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    message("Using a proportion threshold")
    threshold.id <- "of proportion"
    
    multiplication.number <- 1/pop.number
    
  } else if (stri_detect_fixed(pop.percent, "T")) {
    message("Using a percentage threshold")
    threshold.id <- "percent"
    multiplication.number <- 100/pop.number
    
  } else {
    message("Using a fixed threshold")
    threshold.id <- "population(s) as a fixed"
    multiplication.number <- 1
    
  }
  
  pop.filter <- ind.filter %>%
    group_by (LOCUS) %>%
    filter((n_distinct(POP_ID) * multiplication.number) >= pop.threshold)  
  
  message("Calculating summary statistics...")
  
  vcf.summary <- pop.filter %>%
    filter(GT != "./.") %>%
    group_by(LOCUS, POS, POP_ID) %>%
    summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT[GT == "0/0"])),
      PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
      QQ = as.numeric(length(GT[GT == "1/1"]))
    ) %>%
    mutate(
      FREQ_REF = ((PP*2) + PQ)/(2*N),
      FREQ_ALT = ((QQ*2) + PQ)/(2*N),
      HET_O = PQ/N,
      HET_E = 2 * FREQ_REF * FREQ_ALT,
      FIS = ifelse(HET_O == 0, 0, round (((HET_E - HET_O) / HET_E), 6))
    )
  
  global.maf <- vcf.summary %>%
    group_by(LOCUS, POS) %>%
    summarise_each_(funs(sum), vars = c("N", "PP", "PQ", "QQ")) %>%
    mutate(GLOBAL_MAF = (PQ + (2 * QQ)) / (2*N)) %>%
    select(LOCUS, POS, GLOBAL_MAF)
  
  vcf.prep <- global.maf %>%
    left_join(vcf.summary, by = c("LOCUS", "POS"))
  
  
  
  message("Filtering: MAF, HET and Fis")
  
  if(stri_detect_fixed(het.pop.threshold, ".") & het.pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message("Using a proportion for the HET threshold...")
    het.threshold.id <- "of proportion"
  } else if (stri_detect_fixed(het.percent, "T")) {
    multiplication.number <- 100/pop.number
    message("Using a percentage for the HET threshold...")
    het.threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message("Using a fixed HET threshold...")
    het.threshold.id <- "population as a fixed"
  }
  
  if(stri_detect_fixed(fis.pop.threshold, ".") & fis.pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message("Using a proportion for the Fis threshold...")
    fis.threshold.id <- "of proportion"
  } else if (stri_detect_fixed(fis.percent, "T")) {
    multiplication.number <- 100/pop.number
    message("Using a percentage for the Fis threshold...")
    fis.threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message("Using a fixed Fis threshold...")
    fis.threshold.id <- "population as a fixed"
  }
  
  
  all.filters <- vcf.prep %>%
    select(LOCUS, POS, POP_ID, GLOBAL_MAF, FREQ_ALT, HET_O, FIS) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      GLOBAL_MAF = mean(GLOBAL_MAF, na.rm = TRUE),  
      LOCAL_MAF = min(FREQ_ALT, na.rm = TRUE),
      HET_DIFF = max(HET_O) - min(HET_O),
      HET_MAX = max(HET_O),
      FIS_MIN = min(FIS),
      FIS_MAX = max(FIS),
      FIS_DIFF = FIS_MAX-FIS_MIN,
      SNP_N = n_distinct(POS)
    ) %>%
    group_by(LOCUS) %>%
    summarise(
      MAF = length(POP_ID[LOCAL_MAF >= local.maf.threshold | GLOBAL_MAF >= global.maf.threshold]),
      HET = length(POP_ID[HET_DIFF <= het.diff.threshold & HET_MAX <= het.threshold]),
      FIS = length(POP_ID[FIS_MIN >= fis.min.threshold & FIS_MAX <= fis.max.threshold & FIS_DIFF <= fis.diff.threshold]),
      SNP_N = max(SNP_N)
    ) %>%
    filter(MAF >= maf.pop.threshold) %>% 
    filter((HET * multiplication.number) >= het.pop.threshold) %>%
    filter((FIS * multiplication.number) >= fis.pop.threshold) %>%
    filter(SNP_N <= max.snp.number)
  
}
