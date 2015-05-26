#' @title Genotype likelihood filter
#' @description Filter markers based coverage and/or genotype likelihood
#' using a tidy VCF file.
#' @param tidy.vcf.file A data frame object or file (using ".tsv")
#' of class tidy VCF.
#' @param allele.min.depth.threshold Threshold number
#' of min depth of REF or ALT alleles.
#' @param read.depth.max.threshold Threshold number
#' of maximum read depth.
#' @param gl.mean.threshold Threshold number of mean genotype likelihood. 
#' @param gl.min.threshold Threshold number of min genotype likelihood. 
#' @param gl.diff.threshold Threshold number of diff genotype likelihood,
#' the difference between the max and min GL over a loci/read/haplotype. 
#' @param pop.threshold A threshold number: proportion, percentage
#' or fixed number e.g. 0.70, 70 or 15.
#' @param percent Is the threshold a percentage ? TRUE or FALSE.
#' @param filename Name of the file written to the working directory.
#' @details The summary statistics are averaged
#' and nested SNP -> individuals -> population -> loci. e.g. the mean GL is the average
#' genotype likelihood for all individuals of pop x for loci x.
#' Here the filter focus on global trend accross pop for individual loci.
#' The gl.diff.threshold can spot big difference in GL for a loci, e.g. when
#' 3 SNP have similar GL the difference will be close to 0, higher value, higher
#' the difference between the SNP inside the loci.
#' @rdname filter_genotype_likelihood
#' @export

filter_genotype_likelihood <- function (tidy.vcf.file, allele.min.depth.threshold, read.depth.max.threshold, gl.mean.threshold, gl.min.threshold, gl.diff.threshold, pop.threshold, percent, filename) {
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, col_names = T)
    message("Using the tidy vcf file in your directory")
  } else {
    data <- tidy.vcf.file
    message("Using the tidy vcf file from your global environment")
  }
  
  pop.number <- n_distinct(data$POP_ID)
  
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
  
  
  filter <- data %>%
    group_by(LOCUS, POP_ID) %>% # at the population level
    summarise(
      MIN_REF = min(ALLELE_REF_DEPTH, na.rm = T),
      MIN_ALT = min(ALLELE_ALT_DEPTH, na.rm = T),
      READ_DEPTH_MAX = max(READ_DEPTH, na.rm = T),
      GL_MEAN = mean(GL, na.rm = T),
      GL_MEDIAN = median(GL, na.rm = T),
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
