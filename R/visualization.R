#### Visualization

#' @title Figure density distribution of coverage summary statistics
#' @description Create density distribution of coverage summary statistics.
#' Use the coverage summary file created with coverage_summary function.
#' @param data Coverage summary file.
#' @param aes.colour GGPLOT2 aesthetics colour, 
#' e.g. aes(y = ..scaled.., color = COVERAGE_GROUP).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @export
#' @rdname figure_density_distribution_coverage
#' @import ggplot2
#' @import dplyr
#' @import readr

figure_density_distribution_coverage <- function(data, aes.colour, adjust.bin) {
  
  VALUE <- NULL

  
  
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




#' @title Figure box plot of coverage summary statistics
#' @description Create box plots of coverage summary statistics.
#' Use the coverage summary file created with coverage_summary function.
#' @param data Coverage summary file.
#' @export
#' @rdname figure_box_plot_coverage

figure_box_plot_coverage <- function(data) {
  
  POP_ID <- NULL
  VALUE <- NULL
  POP_ID <- NULL
  VALUE <- NULL
  

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




#' @title Visual diagnostic of coverage imbalance
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
#' @export
#' @rdname figure_coverage_imbalance_diagnostic


figure_coverage_imbalance_diagnostic <- function(tidy.vcf.file, pop.levels, read.depth.threshold, aes.colour, adjust.bin) {
  
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  READ_DEPTH <- NULL
  GL <- NULL
  GROUP_COVERAGE <- NULL
  GROUP_GL <- NULL
  ALLELE_COVERAGE_RATIO <- NULL
  
  
  
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


#' @title Figure density distribution of genotype likelihood summary statistics
#' @description Create density distribution of genotype likelihood 
#' summary statistics.
#' Use the long version of coverage summary file created with 
#' genotype_likelihood_summary function ($gl.summary.long).
#' @param data Genotype likelihood summary file.
#' @param aes.colour GGPLOT2 aesthetics colour, 
#' e.g. aes(y = ..scaled.., color = GENOTYPE_LIKELIHOOD_GROUP).
#' @param adjust.bin Adjust GGPLOT2 bin size (0 to 1).
#' @export
#' @rdname figure_density_distribution_genotype_likelihood
#' @seealso \link{summary_genotype_likelihood}

figure_density_distribution_genotype_likelihood <- function(data, aes.colour, 
                                                            adjust.bin) {
  
  VALUE <- NULL

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
#' @title Figure box plot of genotype likelihood summary statistics
#' @description Create box plots of genotype likelihood summary statistics.
#' Use the genotype likelihood summary file created
#' with genotype_likelihood_summary function.
#' @param data genotype likelihood summary file.
#' @export
#' @rdname figure_box_plot_genotype_likelihood

figure_box_plot_genotype_likelihood <- function(data) {
  
  POP_ID <- NULL
  VALUE <- NULL
  
  
  
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
#' @export
#' @rdname figure_density_distribution_maf

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
#' @export
#' @rdname figure_density_distribution_het

figure_density_distribution_het <- function(data, pop.levels, het.group, aes.colour, adjust.bin, x.title){
  
  POP_ID <- NULL
  HET_O <- NULL
  HET_MAX <- NULL
  HET_MIN <- NULL
  VALUE <- NULL
  
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



#' @title Figure of the distribution of SNP per locus before and after filters
#' @param before.filter.data Data set before filter.
#' @param after.filter.data Data set after filter.
#' @export
#' @rdname figure_snp_number_loci

figure_snp_number_loci <- function(before.filter.data, after.filter.data) {
  
  GROUP <- NULL
  SNP_N <- NULL
  POP_ID <- NULL
  
  
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




#' @title Figure of the distribution of SNP nucleotide position alond the read
#' @param data Data for the figure.
#' @param aes.colour GGPLOT2 aesthetic.
#' @param y.title Title of the Y-axis.
#' @export
#' @rdname figure_nucleotide_number_position



figure_nucleotide_number_position <- function(data, aes.colour, y.title) {
  
  COL <- NULL
  
  
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




#' @title Density distribution of diversity (Gene and Haplotypes)
#' @description GGPLOT2 Density distribution of diversity (Gene and Haplotypes).
#' @param data The hapstats summary file or object.
#' @param aes.x GGPLOT2 aesthetics, 
#' e.g. aes.x = aes(x = GENE_DIVERSITY, na.rm = T).
#' @param aes.colour GGPLOT2 aesthetics colour, 
#' e.g. aes.colour = aes(y = ..scaled.., colour = POP_ID).
#' @param x.title Title of the x-axis.
#' @param y.title Title of the y-axis.
#' @export
#' @rdname figure_distribution_diversity

figure_distribution_diversity <- function(data, aes.x, aes.colour, x.title, y.title) {
  
  hapstats.summary <- NULL

  
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




#' @title Box plot of the diversity (Gene and Haplotypes)
#' @description GGPLOT2 Box plot of the diversity (Gene and Haplotypes).
#' @param data The hapstats summary file or object.
#' @param aes.x.y The GGPLOT2 aesthetics, 
#' e.g. aes.x.y = aes(x = factor(POP_ID), y = GENE_DIVERSITY, na.rm = T). 
#' @param y.title Title of the y-axis.
#' @export
#' @rdname figure_box_plot_diversity
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

