## Summary and tables


## VCF
#' @title Summary statistics of a tidy VCF by population and markers.
#' @description Summarise and prepare the tidy VCF. 
#' Summary, by population and markers (SNP), of frequency of the REF 
#' and the ALT alleles, the observed and the expected heterozygosity 
#' and the inbreeding coefficient. The Global MAF of Loci, 
#' with STACKS GBS/RAD loci = read or de novo haplotypes, 
#' is included and repeated over SNP.
#' @param data The tidy VCF file created with read_stacks_vcf.
#' @rdname summary_vcf_tidy
#' @export

summary_vcf_tidy <- function(data) {

  vcf.summary <- data %>%
    filter(GT != "./.") %>%
    group_by(LOCUS, POS, POP_ID) %>%
    summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT[GT == "0/0"])),
      PQ = as.numeric(length(GT[GT == "0/1" | GT == "1/0"])),
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
  
  return(vcf.prep)
}



## Coverage
#' @title Coverage summary.
#' @description This function create a table summary of the important
#' coverage statistics from the tidy vcf created with read_stacks_vcf.
#' @param tidy.vcf.file The tidy VCF file created with read_stacks_vcf.
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
#' @export

summary_coverage <- function (tidy.vcf.file, pop.levels, filename) {
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, col_names = T, col_types = "iiiiccddcdccddddc")
    message("Using the file in your directory")
    
  } else {
    data = tidy.vcf.file
    message("Using the file from your global environment")
    
  }
  
  coverage.sum.loci <- data %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      READ_MEAN = mean(READ_DEPTH, na.rm = T),
      READ_MEDIAN = median(READ_DEPTH, na.rm = T),
      READ_MIN = min(READ_DEPTH, na.rm = T),
      READ_MAX = max(READ_DEPTH, na.rm = T),
      REF_MEAN = mean(ALLELE_REF_DEPTH, na.rm = T),
      REF_MEDIAN = median(ALLELE_REF_DEPTH, na.rm = T),
      REF_MIN = min(ALLELE_REF_DEPTH, na.rm = T),
      REF_MAX = max(ALLELE_REF_DEPTH, na.rm = T),
      ALT_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = T),
      ALT_MEDIAN = median(ALLELE_ALT_DEPTH, na.rm = T),
      ALT_MIN = min(ALLELE_ALT_DEPTH, na.rm = T),
      ALT_MAX = max(ALLELE_ALT_DEPTH, na.rm = T)
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
      READ_DEPTH_MEAN = mean(READ_DEPTH, na.rm = T),
      READ_DEPTH_MEDIAN = median(READ_DEPTH, na.rm = T),
      READ_DEPTH_MIN = min(READ_DEPTH, na.rm = T),
      READ_DEPTH_MAX = max(READ_DEPTH, na.rm = T),
      ALLELE_REF_DEPTH_MEAN = mean(ALLELE_REF_DEPTH, na.rm = T),
      ALLELE_REF_DEPTH_MEDIAN = median(ALLELE_REF_DEPTH, na.rm = T),
      ALLELE_REF_DEPTH_MIN = min(ALLELE_REF_DEPTH, na.rm = T),
      ALLELE_REF_DEPTH_MAX = max(ALLELE_REF_DEPTH, na.rm = T),
      ALLELE_ALT_DEPTH_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = T),
      ALLELE_ALT_DEPTH_MEDIAN = median(ALLELE_ALT_DEPTH, na.rm = T),
      ALLELE_ALT_DEPTH_MIN = min(ALLELE_ALT_DEPTH, na.rm = T),
      ALLELE_ALT_DEPTH_MAX = max(ALLELE_ALT_DEPTH, na.rm = T)
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




#' @title Table of low coverage genotypes.
#' @description This function create a table summary of the genotypes
#' below a user-define threshold.
#' coverage statistics by populations.
#' @param tidy.vcf.file The tidy VCF file created with read_stacks_vcf.
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


table_low_coverage_summary <- function(tidy.vcf.file,
                                       pop.levels, 
                                       read.depth.threshold,
                                       filename.low.coverage,
                                       filename.low.coverage.imbalance) {
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, 
                     col_names = T, 
                     col_types = "diidccddccccdddddc") %>%
      mutate(INDIVIDUALS = factor(INDIVIDUALS))
    message("Using the file in your directory")
    
  } else {
    data <- tidy.vcf.file
    message("Using the file from your global environment")
    
  }
  
  if (missing(pop.levels) == "TRUE") {
    data <- tidy.vcf.file
    
  } else {
    data <- tidy.vcf.file %>%
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

#' @title Genotype likelihood summary.
#' @description This function create 2 tables summary of the important
#' genotype likelihood statistics from the tidy vcf created with read_stacks_vcf.
#' @param tidy.vcf.file The tidy VCF file created with read_stacks_vcf.
#' @param pop.levels Character string defining your ordered populations.
#' @return A list with 2 tables: the long format of loci and populations
#' genotype likelihood statistics and the short format by populations.
#' The short-format is more user-friendly and
#' is written to the working directory.
#' @details The table contains summary statistics: mean, median, min, max and 
#' diff (max-min), of genotype likelihood by locus and populations. To access 
#' the two tables, use $. The table that summarise by populations was created
#' using average nested: loci -> individuals -> populations.

summary_genotype_likelihood <- function(tidy.vcf.file, pop.levels, filename){
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, col_names = T, col_types = "iiiiccddcdccddddc")
    message("Using the file in your directory")
  } else {
    data <- tidy.vcf.file
    message("Using the file from your global environment")
    
  }
  
  GL.loci.pop <- data %>%
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
    )
  
  if (missing(pop.levels) == "TRUE") {
    GL.loci.pop.summary <- GL.loci.pop
  } else {
    GL.loci.pop.summary <- GL.loci.pop %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)) %>%
      arrange(POP_ID)
  }
  
  write.table(GL.loci.pop.summary, 
              filename,
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F
  )
  
  GL.pop <- data %>%
    group_by(POP_ID, INDIVIDUALS) %>%
    summarise(
      GL_MEAN = mean(GL, na.rm = T),
      GL_MEDIAN = median(GL, na.rm = T),
      GL_MIN = min(GL, na.rm = T),
      GL_MAX = max(GL, na.rm = T),
      GL_DIFF = GL_MAX - GL_MIN
    ) %>%
    group_by(POP_ID) %>%
    summarise_each_(funs(mean), vars = c("GL_MEAN", "GL_MEDIAN", "GL_MIN", "GL_MAX", "GL_DIFF")) %>%
    melt(
      id.vars = c("POP_ID"),
      variable.name = "GENOTYPE_LIKELIHOOD_GROUP", 
      value.name = "VALUE"
    )
  
  if (missing(pop.levels) == "TRUE") {
    GL.pop.summary <- GL.pop
  } else {
    GL.pop.summary <- GL.pop %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)) %>%
      arrange(POP_ID)
  }
  
  GL.pop.summary.table <- GL.pop.summary %>%
    dcast(POP_ID ~ GENOTYPE_LIKELIHOOD_GROUP, value.var = "VALUE")
  
  invisible(cat(sprintf(
    "Filename:
%s
Written in the directory:
%s",
    filename, getwd()
  )))
  
  # results
  results <- list()
  results$gl.summary.long <- GL.loci.pop.summary
  results$gl.summary.pop <- GL.pop.summary.table
  
  return(results)
  
}



