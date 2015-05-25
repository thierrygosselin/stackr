#' @title Population filter
#' @description Filter loci based on population threshold with loci: number, 
#' percentage or fix threshold.
#' @param data A data frame object or file (using ".tsv").
#' The data is a sumstats prepared file,
#' a tidy VCF file or any of the last two file format modified by other filters.
#' @param pop.threshold A threshold number: proportion, percentage
#' or fixed number e.g. 0.70, 70 or 15.
#' @param percent Is the threshold a percentage ? TRUE or FALSE.
#' @param filename Name of the file written to the working directory.
#' @rdname filter_population
#' @export

# Population
filter_population <- function(data, pop.threshold, percent, filename) {

  if (is.vector(data) == "TRUE") {
    message("Using the file in your directory")
    data <- read_tsv(data, col_names = T)
  } else {
    message("Using the file from your global environment")
    data <- data
  }
   
  pop.number <- n_distinct(data$POP_ID)
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    message("Using a proportion threshold")
    threshold.id <- "of proportion"
    
    multiplication.number <- 1/pop.number
    
  } else if (stri_detect_fixed(percent, "T")) {
    message("Using a percentage threshold")
    threshold.id <- "percent"
    multiplication.number <- 100/pop.number
    
  } else {
    message("Using a fixed threshold")
    threshold.id <- "population(s) as a fixed"
    multiplication.number <- 1
    
  }
  
  pop.filter <- data %>%
    group_by (LOCUS) %>%
    filter((n_distinct(POP_ID) * multiplication.number) >= pop.threshold) %>%
    arrange (LOCUS, POP_ID)
  
  write.table(pop.filter,
              filename,
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F
              )
  

  invisible(cat(sprintf(
  "Population filter: %s %s threshold as sampling sites/pop required to keep the marker
  The number of SNP removed by the population filter = %s SNP
  The number of LOCI removed by the population filter = %s LOCI
  The number of SNP after the population filter = %s SNP
  The number of LOCI after the population filter = %s LOCI
  Filename:
  %s
  Written in the directory:
  %s", 
  pop.threshold,
  threshold.id,
  n_distinct(data$POS)-n_distinct(pop.filter$POS),
  n_distinct(data$LOCUS)-n_distinct(pop.filter$LOCUS),
  n_distinct(pop.filter$POS),
  n_distinct(pop.filter$LOCUS),
  filename, getwd()
  )))
  
  pop.filter
}
