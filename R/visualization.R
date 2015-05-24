#### Visualization

#BEFORE FILTERS

##### Coverage ###
# Coverage summary function
coverage_summary <- function (stacks.vcf.file, pop.levels) {
    
  if (is.vector(stacks.vcf.file) == "TRUE") {
      data <- read_tsv(stacks.vcf.file, col_names = T, col_types = "iiiiccddcdccddddc")
    } else {
      data = stacks.vcf.file
  }
  
  coverage.sum <- data %>%
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
    coverage <- coverage.sum
  } else {
    coverage.levels <- coverage.sum %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
    coverage <- coverage.levels
  }
  coverage
}

# Figure coverage summary density
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

# Figure coverage summary Box Plot
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

# Coverage summary by pop before filters TABLE
coverage_summary_table_pop <- function(stacks.vcf.file, pop.levels){
      
  if (is.vector(stacks.vcf.file) == "TRUE") {
      data <- read_tsv(stacks.vcf.file, col_names = T, col_types = "iiiiccddcdccddddc")
    } else {
      data = stacks.vcf.file
  }
  
  coverage.sum <- data %>%
    group_by(POP_ID) %>%
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
      )
  coverage.summary.total <- coverage.sum %>%
    summarise_each(funs(mean))
  
  coverage.summary.total[1,1] <- "TOTAL"
  
  if (missing(pop.levels) == "TRUE") {
    coverage.summary.pop.total <- coverage.sum %>%
      rbind(coverage.summary.total)
  } else {
    coverage.summary.pop.total <- coverage.sum %>%
      rbind(coverage.summary.total) %>% 
      mutate(POP_ID = factor(POP_ID, levels = c(pop.levels, "TOTAL"), ordered = T)) %>%
      arrange(POP_ID)
  }
  coverage.summary.pop.total
}

# Figure coverage imbalance allele
figure_coverage_imbalance_alleles <- function(data, pop.levels, aes.colour, adjust.bin) {

  if (is.vector(data) == "TRUE") {
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
      GROUP_GL = ifelse(GL < (mean(data$GL, na.rm = T)), "< mean GL", ">= mean GL"),
      GROUP_COVERAGE = factor(GROUP_COVERAGE, levels = c("READ DEPTH <= 8", "READ DEPTH > 8"), ordered = T),
      GROUP_GL = factor(GROUP_GL, levels = c("< mean GL", ">= mean GL"), ordered = T)
      ) %>%
    filter(GROUP_COVERAGE != "NA" & GROUP_GL != "NA")

  ggplot(imbalance.coverage, aes(x = ALLELE_COVERAGE_RATIO, na.rm = T))+
#   geom_bar()+
    geom_line(aes.colour, stat = "density", adjust = adjust.bin)+
    labs(x="Coverage imbalance between REF and ALT alleles")+
    labs(y="Distribution")+
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = element_text(size = 12, family = "Helvetica", face = "bold")
      )

}

# Table summary of genotype with read coverage < 8
table_low_coverage_summary <- function(data, pop.levels, read.depth.threshold, filename.low.coverage, filename.low.coverage.imbalance) {
  
  if (is.vector(data) == "TRUE") {
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
    LOW_COVERAGE_PERCENT = SAMPLES_NUMBER / TOTAL_NUMBER * 100,
    IMBALANCE_PERCENT = IMBALANCE_NUMBER / TOTAL_NUMBER * 100
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

# Genotype likelihood summary by loci function
genotype_likelihood_summary_loci <- function(stacks.vcf.file, pop.levels){
    
  if (is.vector(stacks.vcf.file) == "TRUE") {
      data <- read_tsv(stacks.vcf.file, col_names = T, col_types = "iiiiccddcdccddddc")
    } else {
      data = stacks.vcf.file
  }
  
  GL <- data %>%
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
#     measure.vars = c(), # if left blank will use all the non id.vars
      variable.name = "GENOTYPE_LIKELIHOOD_GROUP", 
      value.name = "VALUE"
      )

  if (missing(pop.levels) == "TRUE") {
    GL.summary <- GL
  } else {
    GL.summary <- GL %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)) %>%
      arrange(POP_ID)
  }
GL.summary
}

# Figure: Density distribution of genotype likelihood summary of loci function
figure_density_distribution_genotype_likelihood <- function(data, aes.colour) {
# BREAKS <- seq(0, 150, by = 20)
  ggplot(data, aes(x = VALUE, na.rm = T))+
    geom_line(aes.colour, stat = "density", size = 0.5, adjust = 1)+
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

# MAF Figure
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



# AFTER FILTERS
# Distribution number of SNP per LOCUS BEFORE and AFTER filters
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

# DAPC ADMIXTURE
# figure_dapc_admixture <- function(data) {
  figure_dapc_admixture <- function(data, colour.palette, pop.levels) {

  ggplot(data, aes(x = INDIVIDUALS, y = PERC, fill = ANCESTRY))+
  geom_bar(stat = "identity",  position = "fill")+
#   scale_fill_manual(name = "Ancestry K", limits = SITES_LEVELS)+
  scale_fill_manual(name = "Ancestry", values = colour.palette, breaks = pop.levels)+
  scale_y_continuous(expand  =  c(0,  0))+ #in order for the line not to expand beyond the graph!
  labs(y = "Ancestry")+
  labs(x = "Individuals")+
  theme_bw()+
  theme(
    panel.grid.major.x  =  element_blank(), 
    panel.grid.minor.x  =  element_blank(), 
    panel.grid.major.y  =  element_line(colour = "grey60",  linetype = "dashed"), 
    axis.title.x = element_text(size = 10,  family = "Helvetica", face = "bold"), 
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size = 10,  family = "Helvetica", face = "bold"), 
    axis.text.y = element_text(size = 8, family = "Helvetica"),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle=0),
    panel.margin.y = unit(0.5, "lines"),
    panel.margin.x = unit(0.1, "lines"),
    legend.position = "none") +
  guides(fill = guide_legend(
    label.position = "bottom",
    title.position = "bottom",
    title.hjust = 0.5,
    label.hjust = 0.5)
    )
}
