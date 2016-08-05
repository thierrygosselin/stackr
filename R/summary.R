#' @title Import and summarise the batch_x.hapstats.tsv file
#' @description Import and summarise the batch_x.hapstats.tsv file.
#' Necessary preparation for density distribution and box plot figures.
#' @param data The 'batch_x.hapstats.tsv' created by STACKS.
#' @param pop.num The number of populations analysed.
#' @param pop.col.types \code{"integer"} or \code{"character"} used in STACKS populations module?
#' @param pop.integer.equi When Integer was used for your population id,
#' give the character equivalence
#' @param pop.levels A character string with your populations in order.
#' @rdname summary_hapstats
#' @export 

summary_hapstats <- function(data, pop.num, pop.col.types, pop.integer.equi, pop.levels) {
  
  POP_ID <- NULL
  
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
}




## VCF
#' @title Summary statistics of a tidy VCF by population and markers
#' @description Summarise and prepare the tidy VCF. 
#' Summary, by population and markers (SNP), of frequency of the REF 
#' and the ALT alleles, the observed and the expected heterozygosity 
#' and the inbreeding coefficient. The Global MAF of Loci, 
#' with STACKS GBS/RAD loci = read or de novo haplotypes, 
#' is included and repeated over SNP.
#' @param filename (optional) Name of the file written to the working directory.
#' @param data The tidy VCF file created with \link{tidy_genomic_data}.
#' @rdname summary_stats_vcf_tidy
#' @export

summary_stats_vcf_tidy <- function(data, filename) {
  
  
  GT <- NULL
  GL <- NULL
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  N <- NULL
  HET_O <- NULL
  HOM_O <- NULL
  HET_E <- NULL
  HOM_E <- NULL
  FREQ_ALT <- NULL
  FREQ_REF <- NULL
  GLOBAL_MAF <- NULL
  PP <- NULL
  PQ <- NULL
  QQ <- NULL
  
  
  
  
  vcf.summary <- data %>%
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
  
  vcf.prep <- vcf.prep[c("LOCUS", "POS", "POP_ID", "N", "PP", "PQ", "QQ", "FREQ_REF", "FREQ_ALT", "GLOBAL_MAF", "HET_O", "HET_E", "FIS")]
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(vcf.prep, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected"
  }
  
  return(vcf.prep)
}

#' @title Summary statistics of a tidy VCF by population
#' @description Summarise the tidy VCF. 
#' The populations summary on :  frequency of the REF 
#' and the ALT alleles, the observed and the expected heterozygosity 
#' and the inbreeding coefficient. The Global MAF of Loci, 
#' with STACKS GBS/RAD loci = read or de novo haplotypes, 
#' is included and repeated over SNP.
#' @param filename (optional) Name of the file written to the working directory.
#' @param data The tidy VCF file created with tidy_genomic_data.
#' @rdname summary_stats_pop
#' @export

summary_stats_pop <- function(data, filename) {
  
  
  POP_ID <- NULL
  N <- NULL
  HET_O <- NULL
  HET_E <- NULL
  FREQ_REF <- NULL
  FIS <- NULL
  SNP <- NULL
  LOCUS <- NULL
  
  
  vcf.summary <- data %>%
    group_by(POP_ID) %>%
    summarise(
      SNP = length(unique(POS)),
      LOCUS = length(unique(LOCUS)),
      N = max(N, na.rm = TRUE),
      FREQ_REF = mean(FREQ_REF, na.rm = TRUE),
      HET_O = mean(HET_O, na.rm = TRUE),
      HET_E = mean(HET_E, na.rm = TRUE),
      FIS = mean(FIS, na.rm = TRUE)
    ) %>%
    select(POP_ID, N, SNP, LOCUS, FREQ_REF, HET_O, HET_E, FIS)  
  
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(vcf.summary, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected"
  }
  
  return(vcf.summary)
}



## Coverage
#' @title Coverage summary
#' @description This function create a table summary of the important
#' coverage statistics from the tidy vcf created with tidy_genomic_data.
#' @param data A tidy vcf object or file (using ".tsv"), 
#' created with tidy_genomic_data.
#' @param pop.levels Character string defining your ordered populations.
#' @param filename Name of the file saved to the working directory.
#' @details The tables contains summary statistics (mean, median, min, max)
#' of read, ref and alt allele coverage. To access 
#' the two tables, use $. The table that summarise by populations was created
#' using average nested: loci -> individuals -> populations.
#' The long format is used for creating figures.
#' @return A list with 2 tables: the long format of loci and populations
#' coverage statistics and the short format by populations.
#' The short-format is more user-friendly and
#' is written to the working directory.
#' @rdname summary_coverage
#' @importFrom stats median
#' @export

summary_coverage <- function (data, pop.levels, filename) {
  
  POP_ID <- NULL
  READ_DEPTH <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL
  INDIVIDUALS <- NULL
  
  
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T, col_types = "iiiiccddcdccddddc")
    message("Using the file in your directory")
    
  } else {
    data = data
    message("Using the file from your global environment")
    
  }
  
  coverage.sum.loci <- data %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      READ_MEAN = mean(READ_DEPTH, na.rm = TRUE),
      READ_MEDIAN = stats::median(READ_DEPTH, na.rm = TRUE),
      READ_MIN = min(READ_DEPTH, na.rm = TRUE),
      READ_MAX = max(READ_DEPTH, na.rm = TRUE),
      REF_MEAN = mean(ALLELE_REF_DEPTH, na.rm = TRUE),
      REF_MEDIAN = stats::median(ALLELE_REF_DEPTH, na.rm = TRUE),
      REF_MIN = min(ALLELE_REF_DEPTH, na.rm = TRUE),
      REF_MAX = max(ALLELE_REF_DEPTH, na.rm = TRUE),
      ALT_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = TRUE),
      ALT_MEDIAN = stats::median(ALLELE_ALT_DEPTH, na.rm = TRUE),
      ALT_MIN = min(ALLELE_ALT_DEPTH, na.rm = TRUE),
      ALT_MAX = max(ALLELE_ALT_DEPTH, na.rm = TRUE)
    ) %>%
    melt(
      id.vars = c("LOCUS", "POP_ID"),
      #     measure.vars = c(), # if left blank will use all the non id.vars
      variable.name = "COVERAGE_GROUP", 
      value.name = "VALUE"
    )
  if (missing(pop.levels) == "TRUE") {
    coverage <- coverage.sum.loci
  } else {
    coverage <- coverage.sum.loci %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
  }
  
  # by pop
  coverage.sum.pop <- data %>%
    group_by(POP_ID, INDIVIDUALS) %>%
    summarise(
      READ_DEPTH_MEAN = mean(READ_DEPTH, na.rm = TRUE),
      READ_DEPTH_MEDIAN = stats::median(READ_DEPTH, na.rm = TRUE),
      READ_DEPTH_MIN = min(READ_DEPTH, na.rm = TRUE),
      READ_DEPTH_MAX = max(READ_DEPTH, na.rm = TRUE),
      ALLELE_REF_DEPTH_MEAN = mean(ALLELE_REF_DEPTH, na.rm = TRUE),
      ALLELE_REF_DEPTH_MEDIAN = stats::median(ALLELE_REF_DEPTH, na.rm = TRUE),
      ALLELE_REF_DEPTH_MIN = min(ALLELE_REF_DEPTH, na.rm = TRUE),
      ALLELE_REF_DEPTH_MAX = max(ALLELE_REF_DEPTH, na.rm = TRUE),
      ALLELE_ALT_DEPTH_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = TRUE),
      ALLELE_ALT_DEPTH_MEDIAN = stats::median(ALLELE_ALT_DEPTH, na.rm = TRUE),
      ALLELE_ALT_DEPTH_MIN = min(ALLELE_ALT_DEPTH, na.rm = TRUE),
      ALLELE_ALT_DEPTH_MAX = max(ALLELE_ALT_DEPTH, na.rm = TRUE)
    ) %>%
    group_by(POP_ID) %>%
    summarise_each_(funs(mean), vars = c("READ_DEPTH_MEAN", "READ_DEPTH_MEDIAN", "READ_DEPTH_MIN", "READ_DEPTH_MAX", "ALLELE_REF_DEPTH_MEAN", "ALLELE_REF_DEPTH_MEDIAN", "ALLELE_REF_DEPTH_MIN", "ALLELE_REF_DEPTH_MAX", "ALLELE_ALT_DEPTH_MEAN", "ALLELE_ALT_DEPTH_MEDIAN", "ALLELE_ALT_DEPTH_MIN", "ALLELE_ALT_DEPTH_MAX")) %>%
    melt(
      id.vars = c("POP_ID"),
      variable.name = "GENOTYPE_LIKELIHOOD_GROUP", 
      value.name = "VALUE"
    )
  
  coverage.summary.total <- coverage.sum.pop %>%
    summarise_each(funs(mean))
  
  coverage.summary.total[1,1] <- "TOTAL"
  
  if (missing(pop.levels) == "TRUE") {
    coverage.summary.pop.total <- coverage.sum.pop %>%
      rbind(coverage.summary.total)
  } else {
    coverage.summary.pop.total <- coverage.sum.pop %>%
      rbind(coverage.summary.total) %>% 
      mutate(POP_ID = factor(POP_ID, levels = c(pop.levels, "TOTAL"), ordered = T)) %>%
      arrange(POP_ID)
  }
  coverage.summary.pop.total
  
  write.table(coverage.summary.pop.total, filename, sep = "\t", row.names = F, col.names = T, quote = F)
  
  invisible(cat(sprintf(
    "Filename: %s
    Written in this directory: %s", 
    filename,
    getwd()
  )))
  
  # results
  results <- list()
  results$coverage.summary.long <- coverage.sum.loci
  results$coverage.summary.pop <- coverage.summary.pop.total
  
  return(results)
  
}




#' @title Table of low coverage genotypes
#' @description This function create a table summary of the genotypes
#' below a user-define threshold.
#' coverage statistics by populations.
#' @param tidy.vcf A tidy vcf object or file (using ".tsv"), 
#' created with tidy_genomic_data.
#' @param pop.levels Character string defining your ordered populations.
#' @param read.depth.threshold The read depth threshold to evaluate.
#' @param filename.low.coverage Filename of the low coverage table written
#' in the working directory.
#' @param filename.low.coverage.imbalance Filename of ...
#' @return a list of 2 tables (accessed with $). The values in the tables
#' represent percentage of samples.
#' @details work in progress....
#' Table 1: low coverage summary $low.coverage.summary (homo- and
#' hetero- zygotes genotypes).
#' Table 2: summary of coverage imbalance between alleles in the heterozygotes.
#' 0/0 : homozygote REF allele.
#' 1/1 : homozygote ALT allele.
#' 0/1 : heterozygote with coverage REF > ALT allele.
#' 1/0 : heterozygote with coverage REF < ALT allele.
#' @rdname table_low_coverage_summary
#' @export


table_low_coverage_summary <- function(tidy.vcf,
                                       pop.levels, 
                                       read.depth.threshold,
                                       filename.low.coverage,
                                       filename.low.coverage.imbalance) {
  
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  GT <- NULL
  READ_DEPTH <- NULL
  ALLELE_COVERAGE_RATIO <- NULL
  SAMPLES_NUMBER <- NULL
  TOTAL_NUMBER <- NULL
  IMBALANCE_NUMBER <- NULL
  
  if (is.vector(tidy.vcf) == "TRUE") {
    data <- read_tsv(tidy.vcf, 
                     col_names = T, 
                     col_types = "diidccddccccdddddc") %>%
      mutate(INDIVIDUALS = factor(INDIVIDUALS))
    message("Using the file in your directory")
    
  } else {
    data <- tidy.vcf
    message("Using the file from your global environment")
    
  }
  
  if (missing(pop.levels) == "TRUE") {
    data <- tidy.vcf
    
  } else {
    data <- tidy.vcf %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
  }
  
  low.coverage.summary <- data %>%
    filter(GT != "./.") %>%
    select(GT, POP_ID) %>%
    group_by(GT, POP_ID) %>%
    summarise(
      TOTAL_NUMBER = n()
    ) %>%
    full_join(
      data %>%
        filter(READ_DEPTH < read.depth.threshold & GT != "./.") %>%
        group_by(GT, POP_ID) %>%
        summarise(
          SAMPLES_NUMBER = n()
        ),
      by = c("GT", "POP_ID")
    ) %>%
    full_join(
      data %>%
        filter(READ_DEPTH < read.depth.threshold & GT != "./.") %>%
        filter(ALLELE_COVERAGE_RATIO != "NA" & ALLELE_COVERAGE_RATIO != 0 ) %>%
        group_by(GT, POP_ID) %>%
        summarise(
          IMBALANCE_NUMBER = n()
        ),
      by = c("GT", "POP_ID")
    ) %>%
    mutate(
      LOW_COVERAGE_PERCENT = round(SAMPLES_NUMBER / TOTAL_NUMBER * 100, 2),
      IMBALANCE_PERCENT = round(IMBALANCE_NUMBER / TOTAL_NUMBER * 100, 2)
    )
  
  low.coverage.summary.table <- low.coverage.summary %>%
    dcast(POP_ID ~ GT, value.var = "LOW_COVERAGE_PERCENT")
  
  write.table(low.coverage.summary.table, filename.low.coverage, sep = "\t", row.names = F, col.names = T, quote = F)
  
  low.coverage.imbalance.summary.table <- low.coverage.summary %>%
    filter(GT != "0/0" & GT != "1/1") %>%
    dcast(POP_ID ~ GT, value.var = "IMBALANCE_PERCENT")
  
  write.table(low.coverage.imbalance.summary.table, filename.low.coverage.imbalance, sep = "\t", row.names = F, col.names = T, quote = F)
  
  invisible(
    cat(
      sprintf(
        "2 files:
%s
%s\n
Written in the directory:
%s",
        filename.low.coverage, filename.low.coverage.imbalance, getwd()
      )))
  
  res <- list()
  res$low.coverage.summary <- low.coverage.summary.table
  res$heterozygote.imbalance <- low.coverage.imbalance.summary.table
  
  return(res)
}




## Genotype likelihood ###

#' @title Genotype likelihood summary
#' @description This function create 3 summary tables with 
#' genotype likelihood statistics. Individuals, populations and markers information is provided.
#' The input data is a tidy vcf created.
#' with \link{tidy_genomic_data}.
#' @param tidy.vcf A tidy vcf object or file (using ".tsv"), 
#' created with tidy_genomic_data.
#' @param pop.levels (optional) Character string defining your ordered populations.
#' @param approach Character. By \code{"SNP"} or by \code{"haplotype"}. 
#' The function will consider the SNP or haplotype GL statistics to filter the marker. 
#' Default: \code{approach = "haplotype"}.
#' @param folder (optional) Path to folder where results will be saved. 
#' Default: \code{folder = NULL}, the working directory is used.

#' @return A list with 3 object (data frames): the statistics at the 
#' individual level ($gl.individuals, 
#' saved as: "genotype.likelihood.individual.summary.tsv"), a long format with information by 
#' marker and population ($gl.summary.marker.pop) and the overall summary by 
#' population ($gl.summary.pop). Those 3 data frames are also written in your
#' working directory.
#' @details The table contains summary statistics: mean, median, min, max and 
#' diff (max-min), of genotype likelihood for individuals, markers and populations. To access 
#' the two tables, use $. The table that summarise by populations was created
#' using average nested: loci -> individuals -> populations.
#' Using the haplotype approach: The summary statistics are averaged
#' and nested SNP -> individuals -> population -> loci. e.g. the mean GL is the average
#' genotype likelihood for all individuals of pop x for loci x.

#' @importFrom stats median
#' @rdname summary_genotype_likelihood
#' @export

summary_genotype_likelihood <- function(tidy.vcf, pop.levels, approach, folder){
  
  POP_ID <- NULL
  GL <- NULL
  GL_MAX <- NULL
  GL_MIN <- NULL
  INDIVIDUALS <- NULL
  GENOTYPE_LIKELIHOOD_GROUP <- NULL
  
  
  if (missing(pop.levels)) pop.levels <- NULL
  if (missing(approach)) approach <- "haplotype"
  
  if (approach == "haplotype"){
    message("Approach: haplotype")
  } else {
    message("Approach: SNP")
  }
  
  if (missing(folder)) folder <- NULL
  
  # set the working directory to save results
  if (is.null(folder)) {
    folder <- getwd()
  }
  
  # import data
  message("Importing data")
  if (is.vector(tidy.vcf) == "TRUE") {
    data <- data.table::fread(
      input = tidy.vcf,
      sep = "\t",
      stringsAsFactors = FALSE, 
      header = TRUE,
      showProgress = TRUE,
      verbose = FALSE
    ) %>% 
      as_data_frame()
    message("Using the file in your directory")
  } else {
    data <- tidy.vcf
    message("Using the file from your global environment")
    
  }
  
  # make sure it's a tidy vcf with Allele Depth info
  columns.names.tidy.vcf <- names(data) 
  if("ALLELE_REF_DEPTH" %in% columns.names.tidy.vcf){
    message("Looking for genotype likelihood information")
  }else{
    stop("This is not a tidy vcf with Allele Depth information,
         you need to run STACKS with a version > 1.29 for this to work.")
  }
  
  if (!is.null(pop.levels)) {
    data <- data %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)) %>%
      arrange(POP_ID)
  }
  
  # summary individuals
  message("Summarising GL info: individual")
  gl.individuals <- data %>% 
    group_by(POP_ID, INDIVIDUALS) %>% 
    summarise(
      GL_MEAN = mean(GL, na.rm = TRUE),
      GL_MEDIAN = median(GL, na.rm = TRUE),
      GL_MIN = min(GL, na.rm = TRUE),
      GL_MAX = max(GL, na.rm = TRUE),
      GL_DIFF = GL_MAX - GL_MIN
    ) %>%  
    arrange(POP_ID, INDIVIDUALS)
  write_tsv(x = gl.individuals, path = stri_paste(folder,"/genotype.likelihood.individual.summary.tsv"), col_names = TRUE)
  
  # summary by locus and pop
  message("Summarising GL info: marker and pop")
  if (approach == "haplotype"){ # by haplotype
    gl.summary.marker.pop <- data %>%
      ungroup() %>% 
      group_by(LOCUS, POP_ID) %>%
      summarise(
        GL_MEAN = mean(GL, na.rm = T),
        GL_MEDIAN = median(GL, na.rm = T),
        GL_MIN = min(GL, na.rm = T),
        GL_MAX = max(GL, na.rm = T),
        GL_DIFF = GL_MAX - GL_MIN
      ) %>%
      melt(
        id.vars = c("LOCUS", "POP_ID"),
        variable.name = "GENOTYPE_LIKELIHOOD_GROUP", 
        value.name = "VALUE"
      ) %>% 
      arrange(LOCUS, POP_ID, GENOTYPE_LIKELIHOOD_GROUP)
  } else { # by SNP
    gl.summary.marker.pop <- data %>%
      ungroup() %>% 
      group_by(LOCUS, POS, POP_ID) %>%
      summarise(
        GL_MEAN = mean(GL, na.rm = T),
        GL_MEDIAN = median(GL, na.rm = T),
        GL_MIN = min(GL, na.rm = T),
        GL_MAX = max(GL, na.rm = T),
        GL_DIFF = GL_MAX - GL_MIN
      ) %>%
      melt(
        id.vars = c("LOCUS", "POS", "POP_ID"),
        variable.name = "GENOTYPE_LIKELIHOOD_GROUP", 
        value.name = "VALUE"
      ) %>% 
      arrange(LOCUS, POS, POP_ID, GENOTYPE_LIKELIHOOD_GROUP)
  }
  write_tsv(x = gl.summary.marker.pop, path = stri_paste(folder,"/genotype.likelihood.summary.marker.pop.tsv"), col_names = TRUE)
  
  # Summary by pop overall locus
  message("Summarising GL info: pop")
  gl.summary.pop <- gl.individuals %>%
    group_by(POP_ID) %>%
    summarise_each_(funs(mean), vars = c("GL_MEAN", "GL_MEDIAN", "GL_MIN", "GL_MAX", "GL_DIFF"))
  write_tsv(x = gl.summary.pop, path = stri_paste(folder,"/genotype.likelihood.summary.pop.tsv"), col_names = TRUE)
  
  
  # results
  results <- list()
  results$gl.individuals <- gl.individuals
  results$gl.summary.marker.pop <- gl.summary.marker.pop
  results$gl.summary.pop <- gl.summary.pop
  return(results)
}




#' @title Import and summarise the batch_x.phistats.tsv file
#' @description Import and summarise the batch_x.phistats.tsv file.
#' Necessary preparation for density distribution and box plot figures.
#' @param data The 'batch_x.phistats.tsv' created by STACKS.
#' @param skip.lines The number of line without the header 
#' to start reading the data.
#' @rdname summary_phistats
#' @export

summary_phistats <- function(data, skip.lines) {
  
  BATCH_ID <- NULL
  LOCUS <- NULL
  CHR <- NULL
  BP <- NULL
  POP_COUNT <- NULL
  PHI_ST <- NULL
  SMOOTHED_PHI_ST <- NULL
  SMOOTHED_PHI_ST_P_VALUE <- NULL
  PHI_CT <- NULL
  SMOOTHED_PHI_CT <- NULL
  SMOOTHED_PHI_CT_P_VALUE <- NULL
  PHI_SC <- NULL
  SMOOTHED_PHI_SC <- NULL
  SMOOTHED_PHI_SC_P_VALUE <- NULL
  FST_PRIME <- NULL
  SMOOTHED_FST_PRIME <- NULL
  SMOOTHED_FST_PRIME_P_VALUE <- NULL
  D_EST <- NULL
  SMOOTHED_D_EST <- NULL
  SMOOTHED_D_EST_P_VALUE <- NULL
  
  
  
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

