# Create a whitelist from an object with LOCUS column

#' @title Whitelist loci
#' @description This function creates a whitelist of loci, used after applying a filter
#' to keep track of the loci kept by a filter.
#' @param data The data frame after the filter. Object
#' or file (using ".tsv") of class sumstats. 
#' @param filename The name of the file written in the directory.
#' @param col.header TRUE and the loci will have a column header 'LOCUS'.
#' @rdname whitelist_loci
#' @export

whitelist_loci <- function(data, filename, col.header) {
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
  } else {
    data <- data
  }
  
  whitelist <- data %>% select(LOCUS) %>% distinct(LOCUS) %>% arrange(LOCUS)
  
    write.table(whitelist, filename, sep = "\t", row.names = F, col.names = col.header,
              quote = F)
invisible(
  cat(
    sprintf(
"Whitelist filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
  whitelist
}




#' @title Whitelist loci and snp.
#' @description This function creates a whitelist of loci and snp,
#' useful in the populations module of STACKS.
#' @param data The data frame after the filter. Object
#' or file (using ".tsv") of class sumstats. 
#' @param filename The name of the file written in the directory.
#' @param col.header TRUE and the loci will have a column header 'LOCUS'.
#' @rdname whitelist_loci_snp
#' @export

whitelist_loci_snp <- function(data, filename, col.header) {
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
  } else {
    data <- data
  }
  
  whitelist <- data %>%
    group_by(LOCUS, POS) %>%
    select(LOCUS, POS) %>%
    distinct(POS) %>% 
    arrange(LOCUS, POS)
  
  write.table(whitelist, filename, sep = "\t", row.names = F, col.names = col.header,
              quote = F)
  invisible(
  cat(
    sprintf(
"Whitelist filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
  whitelist
}




#' @title Whitelist loci for VCF tools
#' @description This function creates a whitelist of loci for VCF tools.
#' With 2 columns (CHROM and ID).
#' @param data The data frame after the filter. Object
#' or file (using ".tsv") of class sumstats. 
#' @param filename The name of the file written in the directory.
#' @rdname whitelist_loci_vcf
#' @export

whitelist_loci_vcf <- function(data, filename) {
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
  } else {
    data <- data
  }
  
  whitelist <- data %>% group_by(POS) %>% select(POS) %>% distinct(POS) %>% 
    mutate(CHROM = rep("1", n())) %>% 
    arrange(CHROM, POS) %>% 
    group_by(CHROM) %>%
    select(CHROM, POS)
  
  write.table(whitelist, filename, sep = "\t", row.names = F, col.names = T,
              quote = F)
  invisible(
  cat(
    sprintf(
"VCF whitelist filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
  whitelist
}



#' @title Blacklist loci
#' @description This function creates a blacklist of loci, used after applying a filter
#' to keep track of the loci removed by a filter.
#' @param before.filter.data The data frame before the filter. Object
#' or file (using ".tsv") of class sumstats. 
#' @param after.filter.data The data frame after the filter. Object
#' or file (using ".tsv") of class sumstats. 
#' @param filename The name of the file written in the directory.
#' @param col.header TRUE and the loci will have a column header 'LOCUS'.
#' @rdname blacklist_loci
#' @export

blacklist_loci <- function(before.filter.data, after.filter.data, filename, 
                           col.header) {
  
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
  
  blacklist <- before.filter.data %>%
    group_by(LOCUS) %>% 
    select(LOCUS) %>% 
    distinct(LOCUS) %>% 
    anti_join(after.filter.data %>% 
                group_by(LOCUS) %>% 
                select(LOCUS) %>%
                distinct(LOCUS),
              by = "LOCUS") %>%
    arrange(LOCUS)
  
  write.table(blacklist, filename, sep = "\t", row.names = F,
              col.names = col.header, quote = F)
  
  invisible(
  cat(
    sprintf(
"Blacklist filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
  
  blacklist
}


