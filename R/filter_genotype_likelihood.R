#' @title Genotype likelihood filter
#' @description Filter markers based coverage and/or genotype likelihood
#' using a tidy VCF file.
#' @param tidy.vcf A tidy vcf object or file (using ".tsv").
#' @param approach Character. By \code{"SNP"} or by \code{"haplotype"}. 
#' The function will consider the SNP or haplotype MAF statistics to filter the marker. 
#' Default by \code{"haplotype"}.

#' @param allele.min.depth.threshold (integer) The minimum depth of coverage for
#' the alleles (REF and ALT).
#' @param read.depth.max.threshold (integer) The maximum depth of coverage for the
#' read.
#' @param gl.mean.threshold Threshold number of mean genotype likelihood. 
#' @param gl.min.threshold Threshold number of min genotype likelihood. 
#' @param gl.diff.threshold Threshold number of diff genotype likelihood,
#' the difference between the max and min GL over a loci/read/haplotype.

#' @param pop.threshold A threshold number: proportion, percentage
#' or fixed number e.g. 0.50, 50 or 5.
#' @param percent Is the threshold a percentage? TRUE or FALSE.
#' This argument is necessary to distinguish percentage from integer population threshold.
#' (e.g. 5 percent or 5 populations).

#' @param filename Name of the file written to the working directory.
#' @details Using the haplotype approach: The summary statistics are averaged
#' and nested SNP -> individuals -> population -> loci. e.g. the mean GL is the average
#' genotype likelihood for all individuals of pop x for loci x.
#' Here the filter focus on global trend accross pop for individual loci.
#' The gl.diff.threshold can spot big difference in GL for a loci, e.g. when
#' 3 SNP have similar GL the difference will be close to 0, higher value, higher
#' the difference between the SNP inside the loci.
#' @rdname filter_genotype_likelihood
#' @export
#' @import dplyr
#' @import readr

#' @examples
#' \dontrun{
#' If I have 10 populations, the codes below will give the same numbers:
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(tidy.vcf = beluga.vcf.tidy, 
#' approach = "haplotype", allele.min.depth.threshold = 10, 
#' read.depth.max.threshold = 100, gl.mean.threshold = 20, 
#' gl.min.threshold = 5, gl.diff.threshold = 100, 
#' pop.threshold = 50, percent = TRUE)
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(tidy.vcf = beluga.vcf.tidy, 
#' approach = "haplotype", allele.min.depth.threshold = 10, 
#' read.depth.max.threshold = 100, gl.mean.threshold = 20, 
#' gl.min.threshold = 5, gl.diff.threshold = 100, 
#' pop.threshold = 0.5, percent = FALSE)
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(tidy.vcf = beluga.vcf.tidy, 
#' approach = "haplotype", allele.min.depth.threshold = 10, 
#' read.depth.max.threshold = 100, gl.mean.threshold = 20, 
#' gl.min.threshold = 5, gl.diff.threshold = 100, 
#' pop.threshold = 5, percent = FALSE)
#' }

filter_genotype_likelihood <- function (tidy.vcf, approach = "haplotype", allele.min.depth.threshold, read.depth.max.threshold, gl.mean.threshold, gl.min.threshold, gl.diff.threshold, pop.threshold, percent, filename) {
  
  
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL
  READ_DEPTH <- NULL
  GL <- NULL
  GL_MAX <- NULL
  GL_MIN <- NULL
  MIN_REF <- NULL
  MIN_ALT <- NULL
  READ_DEPTH_MAX <- NULL
  GL_MEAN <- NULL
  GL_DIFF <- NULL
  POP_ID <- NULL
  
  
  if (is.vector(tidy.vcf) == "TRUE") {
    data <- read_tsv(tidy.vcf, col_names = TRUE)
    message("Using the tidy vcf file in your directory")
  } else {
    data <- tidy.vcf
    message("Using the tidy vcf file from your global environment")
  }
  
  pop.number <- n_distinct(data$POP_ID) # get the number of population
  
  # make sure it's a tidy vcf with Allele Depth info
  columns.names.tidy.vcf <- names(data) 
  if("ALLELE_REF_DEPTH" %in% columns.names.tidy.vcf){
    message("Looking for genotype likelihood and coverage information")
  }else{
    stop("This is not a tidy vcf with Allele Depth information,
         you need to run STACKS with a version > 1.29 for this to work.")
  }
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message("Using a proportion threshold...")
    threshold.id <- "of proportion"
  } else if (stri_detect_fixed(percent, "T")) {
    multiplication.number <- 100/pop.number
    message("Using a percentage threshold...")
    threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message("Using a fixed threshold...")
    threshold.id <- "population as a fixed"
  }
  
  if (missing(approach) | approach == "haplotype"){
    message("Approach selected: haplotype")
    filter <- data %>%
      group_by(LOCUS, POP_ID) %>% # at the population level
      summarise(
        MIN_REF = min(ALLELE_REF_DEPTH, na.rm = T),
        MIN_ALT = min(ALLELE_ALT_DEPTH, na.rm = T),
        READ_DEPTH_MAX = max(READ_DEPTH, na.rm = T),
        GL_MEAN = mean(GL, na.rm = T),
        GL_MEDIAN = stats::median(GL, na.rm = T),
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
      filter((n * multiplication.number) >= pop.threshold) %>%
      select(LOCUS) %>%
      left_join(data, by = "LOCUS") %>%
      arrange(LOCUS, POS, POP_ID)
  } else {
    message("Approach selected: SNP")
    filter <- data %>%
      group_by(LOCUS, POS, POP_ID) %>% # at the population level
      summarise(
        MIN_REF = min(ALLELE_REF_DEPTH, na.rm = T),
        MIN_ALT = min(ALLELE_ALT_DEPTH, na.rm = T),
        READ_DEPTH_MAX = max(READ_DEPTH, na.rm = T),
        GL_MEAN = mean(GL, na.rm = T),
        GL_MEDIAN = stats::median(GL, na.rm = T),
        GL_MIN = min(GL, na.rm = T),
        GL_MAX = max(GL, na.rm = T),
        GL_DIFF = GL_MAX - GL_MIN
      ) %>%
      filter(MIN_REF >= allele.min.depth.threshold | MIN_ALT >= allele.min.depth.threshold) %>%
      filter(READ_DEPTH_MAX <= read.depth.max.threshold) %>% # read coverage
      filter(GL_MEAN >= gl.mean.threshold & GL_MIN >= gl.min.threshold) %>%  # GL
      filter(GL_DIFF <= gl.diff.threshold) %>%
      group_by(LOCUS, POS) %>%
      tally() %>% # Globally accross loci
      filter((n * multiplication.number) >= pop.threshold) %>%
      select(LOCUS, POS) %>%
      left_join(data, by = c("LOCUS", "POS")) %>%
      arrange(LOCUS, POS, POP_ID)
  }
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(filter, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected"
  }
  
  invisible(cat(sprintf(
    "Genotype likelihood (GL) filter:
%s %s threshold of populations to keep the marker
1. Min REF or ALT coverage threshold: %s
2. Max read depth threshold: %s
3. Mean, min and diff (max-min) loci GL threshold: %s mean GL, %s min GL, %s diff GL
The number of markers removed by the GL filter: SNP = %s, LOCI = %s
The number of markers before -> after the GL filter: %s -> %s SNP, %s -> %s LOCI\n
%s\n
Working directory:
%s",
    pop.threshold,
    threshold.id,
    allele.min.depth.threshold,
    read.depth.max.threshold,
    gl.mean.threshold, gl.min.threshold, gl.diff.threshold,
    (n_distinct(data$POS)-n_distinct(filter$POS)), (n_distinct(data$LOCUS)-n_distinct(filter$LOCUS)),
    n_distinct(data$POS), n_distinct(filter$POS), n_distinct(data$LOCUS), n_distinct(filter$LOCUS),
    saving, getwd()
    
  )))
  filter
}

