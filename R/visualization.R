#### Visualization

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

coverage_summary <- function (tidy.vcf.file, pop.levels, filename) {
  
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




#' @title Figure density distribution of coverage summary statistics.
#' @description Create density distribution of coverage summary statistics.
#' Use the coverage summary file created with coverage_summary function.
#' @param data Coverage summary file.
#' @param aes.colour GGPLOT2 aesthetics colour, 
#' e.g. aes(y = ..scaled.., color = COVERAGE_GROUP).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).

figure_density_distribution_coverage <- function(data, aes.colour, adjust.bin) {
  ggplot(data, aes(x = VALUE, na.rm = T))+
    geom_line(aes.colour, stat = "density", size = 0.5, adjust = adjust.bin)+
    labs(x = "Depth of coverage(read number)")+
    labs(y = "Loci (scaled density)")+
    expand_limits(x = 0)+
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold"))
}




#' @title Figure box plot of coverage summary statistics.
#' @description Create box plots of coverage summary statistics.
#' Use the coverage summary file created with coverage_summary function.
#' @param data Coverage summary file.

figure_box_plot_coverage <- function(data) {
  
  ggplot(data, aes(x = factor(POP_ID), y = VALUE, na.rm = T))+
    geom_violin(trim = F)+
    geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
    stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
    labs(x = "Sampling sites")+
    labs(y = "Read depth coverage")+
    facet_wrap(facets = ~COVERAGE_GROUP, scales = "free")+
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
}




#' @title Visual diagnostic of coverage imbalance.
#' @description GBS data and STACKS pipeline sometimes output REF and ALT
#' alleles that have coverage imbalance, i.e. the coverage is not equal 
#' and skewed towards the REF or ALT alleles. 
#' Thw density distribution figure of coverage imbalance between REF and ALT
#' alleles will highlight the problem in you data.
#' @param tidy.vcf.file The tidy VCF file created with read_stacks_vcf.
#' @param pop.levels Character string defining your ordered populations.
#' @param read.depth.threshold Define the threshold you wish to analyse.
#' @param aes.colour GGPLOT2 aesthetics, 
#' e.g. aes(y = ..count..).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @return 4-plots highlighting the different combination of under
#'  or over the coverage threshold and mean genotype likelihood. 
#'  Y- axis show the distribution of genotypes.
#'  X- axis show the coverage imbalance the 
#'  Negative ratio (left x axis) : REF > ALT.
#'  Positive ratio (right x axis) : ALT > REF.
#' @details The figures shows the results of the of coverage threshold
#' selected and mean genotype likelihood. 
#' You can test different threshold to inspect your data.
#' Ideally the lower left pannel of the 4-plot should be empty. If it is, this
#' shows that setting the threshold of the genotype likelihood filter
#' to the mean or close to it take care of the allelic coverage imbalance.
#' #' e.g. fig <- figure_coverage_imbalance_diagnostic(
#' tidy.vcf.file, pop.levels, read.depth.threshold, aes.colour, adjust.bin)
#' Use ( fig + facet_grid(GROUP_GL ~ GROUP_COVERAGE)).


figure_coverage_imbalance_diagnostic <- function(tidy.vcf.file, pop.levels, read.depth.threshold, aes.colour, adjust.bin) {
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, col_names = T, col_types = "diidccddccccdddddc") %>%
      mutate(INDIVIDUALS = factor(INDIVIDUALS))
  } else {
    data <- tidy.vcf.file
  }
  
  if (missing(pop.levels) == "TRUE") {
    data <- data
  } else {
    data <- data %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
  }
  
  below.threshold <- stri_join("READ DEPTH <=", read.depth.threshold, sep = " ")
  over.threshold <- stri_join("READ DEPTH >", read.depth.threshold, sep = " ")
  
  imbalance.coverage <- data %>%
    mutate(
      #       GROUP_COVERAGE = ifelse(READ_DEPTH <= read.depth.threshold, "READ DEPTH <= 8", "READ DEPTH > 8"),
      GROUP_COVERAGE = ifelse(READ_DEPTH <= read.depth.threshold, below.threshold, over.threshold),
      GROUP_GL = ifelse(GL < (mean(data$GL, na.rm = T)), "< mean GL", ">= mean GL"),
      GROUP_COVERAGE = factor(GROUP_COVERAGE, levels = c(below.threshold, over.threshold), ordered = T),
      #       GROUP_COVERAGE = factor(GROUP_COVERAGE, levels = c("READ DEPTH <= 8", "READ DEPTH > 8"), ordered = T),
      GROUP_GL = factor(GROUP_GL, levels = c("< mean GL", ">= mean GL"), ordered = T)
    ) %>%
    filter(GROUP_COVERAGE != "NA" & GROUP_GL != "NA")
  
  ggplot(imbalance.coverage, aes(x = ALLELE_COVERAGE_RATIO, na.rm = T))+
    #   geom_bar()+
    geom_line(aes.colour, stat = "density", adjust = adjust.bin)+
    labs(x="Coverage imbalance between REF and ALT alleles (ratio)")+
    labs(y="Distribution of genotypes")+
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = element_text(size = 12, family = "Helvetica", face = "bold")
    )
  
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




##### Genotype likelihood ###

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

genotype_likelihood_summary <- function(tidy.vcf.file, pop.levels, filename){
  
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




#' @title Figure density distribution of genotype likelihood summary statistics.
#' @description Create density distribution of genotype likelihood 
#' summary statistics.
#' Use the coverage summary file created with 
#' genotype_likelihood_summary function.
#' @param data Genotype likelihood summary file.
#' @param aes.colour GGPLOT2 aesthetics colour, 
#' e.g. aes(y = ..scaled.., color = GENOTYPE_LIKELIHOOD_GROUP).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).

figure_density_distribution_genotype_likelihood <- function(data, aes.colour, 
                                                            adjust.bin) {
  # BREAKS <- seq(0, 150, by = 20)
  ggplot(data, aes(x = VALUE, na.rm = T))+
    geom_line(aes.colour, stat = "density", size = 0.5, adjust = adjust.bin)+
    #  scale_colour_manual(name = "Filters", values = c("black", "blue"))+
    #  scale_x_continuous(breaks = BREAKS)+
    labs(x = "Loci Genotype Likelihood")+
    labs(y = "Loci (scaled density)")+
    expand_limits(x = 0)+
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold"))
}




# Figure: Box Plot of genotype likelihood summary of loci
#' @title Figure box plot of genotype likelihood summary statistics.
#' @description Create box plots of genotype likelihood summary statistics.
#' Use the genotype likelihood summary file created
#' with genotype_likelihood_summary function.
#' @param data genotype likelihood summary file.

figure_box_plot_genotype_likelihood <- function(data) {
  ggplot(data, aes(x = factor(POP_ID), y = VALUE, na.rm = T))+
    geom_violin(trim = F)+
    geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
    stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
    labs(x = "Sampling sites")+
    labs(y = "Genotype likelihood of loci")+
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
}




#' @title Figure density distribution of minor allele frequency (MAF) 
#' summary statistics.
#' @description Create density distribution of MAF summary statistics.
#' @param data sumstats or tidy vcf files.
#' @param maf.group The GGPLOT2 aes (e.g. aes(x = FREQ_ALT, na.rm = F)).
#' @param aes.colour GGPLOT2 aesthetics colour, 
#' e.g. aes(y = ..scaled.., color = GROUP).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @param x.title Title of the x-axis.

figure_density_distribution_maf <- function(data, maf.group, aes.colour, adjust.bin, x.title) {
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
  } else {
    data <- data
  }
  #   font.group <- 
  graph <- ggplot(data, maf.group)+
    geom_line(aes.colour, stat = "density", adjust = adjust.bin)+ # pop colored
    #   scale_colour_manual(name ="Sampling sites", values = colour_palette_sites.pink)+
    scale_x_continuous(breaks=c(0, 0.05, 0.1, 0.2, 0.5, 1),
                       labels = c("0", "0.05", "0.1", "0.2", "0.5", "1.0"))+
    labs(x = x.title)+
    labs(y = "Density of SNP (scaled)")+
    expand_limits(y = 0)+
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
    )
}





# HET Figure
#' @title Figure density distribution of the observed heterozygosity
#' summary statistics.
#' @description Create density distribution of the observed heterozygosity
#' summary statistics.
#' @param data sumstats or tidy vcf files.
#' @param pop.levels Character string defining your ordered populations.
#' @param het.group = aes(x = HET_MAX, na.rm = F)
#' @param aes.colour GGPLOT2 aesthetics colour, 
#' e.g. aes(y = ..scaled.., color = GROUP).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @param x.title Title of the x-axis.
#' 
figure_density_distribution_het <- function(data, pop.levels, het.group, aes.colour, adjust.bin, x.title){
  
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
  } else {
    data = data
  }
  
  data.summary <- data %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      HET_MEAN = mean(HET_O),
      HET_MAX = max(HET_O),
      HET_MIN = min(HET_O),
      HET_DIFF = HET_MAX - HET_MIN
    ) %>%
    melt(
      id.vars = c("LOCUS", "POP_ID"),
      variable.name = "HET_GROUP", 
      value.name = "VALUE"
    )
  
  if (missing(pop.levels) == "TRUE") {
    data.summary <- data.summary
  } else {
    data.summary <- data.summary %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)) %>%
      arrange(POP_ID)
  }
  
  graph <- ggplot(data.summary, aes(x = VALUE, na.rm = F))+
    geom_line(aes.colour, stat = "density", adjust = adjust.bin)+ # pop colored
    #   scale_colour_manual(name ="Sampling sites", values = colour_palette_sites.pink)+
    #     scale_x_continuous(breaks=c(0, 0.05, 0.1, 0.2, 0.5, 1),
    #                        labels = c("0", "0.05", "0.1", "0.2", "0.5", "1.0"))+
    labs(x = x.title)+
    labs(y = "Density of SNP (scaled)")+
    expand_limits(y = 0)+
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
    )
  
}



#' @title Figure of the distribution of SNP per locus before and after filters.
#' @param before.filter.data Data set before filter.
#' @param after.filter.data Data set after filter.
figure_snp_number_loci <- function(before.filter.data, after.filter.data) {
  
  if (is.vector(before.filter.data) == "TRUE") {
    before.filter.data <- read_tsv(before.filter.data, col_names = T)
  } else {
    before.filter.data <- before.filter.data
  }
  
  if (is.vector(after.filter.data) == "TRUE") {
    after.filter.data <- read_tsv(after.filter.data, col_names = T)
  } else {
    after.filter.data <- after.filter.data
  }
  
  number.snp.loci <- before.filter.data %>% # Before
    group_by (LOCUS) %>%
    summarise(SNP_N = n_distinct(POS)) %>%
    mutate(GROUP = rep("pre-filters", n())) %>%
    rbind(
      after.filter.data %>% # After
        group_by (LOCUS) %>%
        summarise(SNP_N = n_distinct(POS)) %>%
        mutate(GROUP=rep("post-filters", n()))) %>%
    mutate (
      GROUP = factor(GROUP, levels = c("pre-filters", "post-filters"),
                     ordered = T)
    )
  
  graph <- ggplot(number.snp.loci, aes(factor(SNP_N)))+
    geom_bar()+
    labs(x="Number of SNP per haplotypes")+
    labs(y="Distribution (number)")+
    facet_wrap(~GROUP, nrow = 1, ncol = 2)+
    theme(axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
          legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
          legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
          strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold"))
  
  graph
}

#' @title Figure of the distribution of SNP nucleotide position alond the read.
#' @param before.filter.data Data set before filter.
#' @param gl.blacklist Blacklist to show the impact of the filter on marker
#' density distribution.

#  Distribution of SNP nucleotide position along the read 
nucleotide_number_position <- function(before.filter.data, 
                                       gl.blacklist, 
                                       maf.blacklist, het.blacklist, fis.blacklist, all.filters.blacklist) {
  
  if (is.vector(before.filter.data) == "TRUE") {
    data <- read_tsv(before.filter.data, col_names = T)
  } else {
    data <- before.filter.data
  }
  
  if (is.vector(gl.blacklist) == "TRUE") {
    gl.bl <- read_tsv(gl.blacklist, col_names = T)
  } else {
    gl.bl <- gl.blacklist
  }
  
  if (is.vector(maf.blacklist) == "TRUE") {
    maf.bl <- read_tsv(maf.blacklist, col_names = T)
  } else {
    maf.bl <- maf.blacklist
  }
  
  if (is.vector(het.blacklist) == "TRUE") {
    het.bl <- read_tsv(gl.blacklist, col_names = T)
  } else {
    het.bl <- het.blacklist
  }
  
  if (is.vector(fis.blacklist) == "TRUE") {
    fis.bl <- read_tsv(fis.blacklist, col_names = T)
  } else {
    fis.bl <- fis.blacklist
  }
  
  if (is.vector(all.filters.blacklist) == "TRUE") {
    all.filters <- read_tsv(all.filters.blacklist, col_names = T)
  } else {
    all.filters <- all.filters.blacklist
  }
  
  GL  <-   gl.bl$LOCUS
  MAF <-  maf.bl$LOCUS
  HET <-  het.bl$LOCUS
  FIS <-  fis.bl$LOCUS
  ALL <-  all.filters$LOCUS
  
  data <- data %>%
    select(POP_ID, LOCUS, POS, COL) %>%
    mutate(
      PRE_FILTERS = rep("pre-filters", n()),
      COL = as.numeric(COL),
      GL_FILTER = ifelse(LOCUS %in% GL, "blacklist", "gl.whitelist"),
      MAF_FILTER = ifelse(LOCUS %in% MAF, "blacklist", "maf.whitelist"),
      HET_FILTER = ifelse(LOCUS %in% HET, "blacklist", "het.whitelist"),
      FIS_FILTER = ifelse(LOCUS %in% FIS, "blacklist", "fis.whitelist"),
      ALL_FILTERS = ifelse(LOCUS %in% ALL, "blacklist", "all.whitelist")
    ) %>%
    melt(
      id.vars = c("POP_ID", "LOCUS", "POS", "COL"),
      variable.name = "FILTERS", 
      value.name = "VALUE"
    ) %>%
    filter(VALUE != "blacklist") %>%
    select(-VALUE) %>%
    mutate(
      FILTERS = factor(FILTERS,
                       levels = c("PRE_FILTERS", "GL_FILTER","MAF_FILTER", "HET_FILTER", "FIS_FILTER", "ALL_FILTERS")),
      GL = ifelse(FILTERS == "PRE_FILTERS" | FILTERS == "GL_FILTER", "GL", "delete"),
      MAF = ifelse(FILTERS == "PRE_FILTERS" | FILTERS == "MAF_FILTER", "MAF", "delete"),
      HET = ifelse(FILTERS == "PRE_FILTERS" | FILTERS == "HET_FILTER", "HET", "delete"),
      FIS = ifelse(FILTERS == "PRE_FILTERS" | FILTERS == "FIS_FILTER", "FIS", "delete"),
      ALL = ifelse(FILTERS == "PRE_FILTERS" | FILTERS == "ALL_FILTERS", "ALL", "delete")
    ) %>%
    melt(
      id.vars = c("POP_ID", "LOCUS", "POS", "COL", "FILTERS"),
      measure.vars = c("GL", "MAF", "HET", "FIS", "ALL"),
      variable.name = "GROUP", 
      value.name = "VALUE"
    ) %>%
    filter(VALUE != "delete") %>%
    select(-VALUE)
}

figure_nucleotide_number_position <- function(data, aes.colour, y.title) {
  
  ggplot(data, aes(x = COL, na.rm = T))+
    geom_line((aes.colour), stat = "density", size = 1, adjust = 0.7)+
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))+
    labs(x = "Nucleotide position on the read")+
    labs(y = y.title)+
    expand_limits(x = 0)+
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
    )
}

