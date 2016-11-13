#' @name snp_ld
#' @title GBS/RADseq short distance linkage disequilibrium pruning
#' @description \strong{For VCF file only}. 
#' SNP short distance linkage disequilibrium pruning. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 3 options (see argument below).
#' 
#' For long distance linkage disequilibrium pruning, see details below.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.

#' @param data A data set with ID (LOCUS) and POS (SNP) information.
#' Usually a VCF file in the working directory or 
#' object in the global environment in wide or long (tidy) formats. 

#' @param snp.ld (character) \strong{For VCF file only}. 3 options: 
#' \code{snp.ld = "random"} selection, \code{snp.ld = "first"} or
#' \code{snp.ld = "last"} for SNP on the same short read/haplotype.
#' Default: \code{snp.ld = "first"}.

#' @export
#' @rdname snp_ld


#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name
#' @importFrom dplyr select distinct group_by sample_n summarise semi_join n_distinct

snp_ld <- function(data, snp.ld = "first") {
  snp.ld <- match.arg(snp.ld, c("first", "random", "last"))
  
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  
  # Import data ---------------------------------------------------------------
  
  if (is.vector(data)) {
    input <- stackr::read_long_tidy_wide(data = data, import.metadata = FALSE)
  } else {
    input <- data
  }
  
  # check genotype column naming ---------------------------------------------
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input), 
      pattern = "GENOTYPE", 
      replacement = "GT", 
      vectorize_all = FALSE)
  }
  
  # Check that fiel format as ID and POS -------------------------------------
  if (!tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "POS")) {
    stop("snp.ld is only available for VCF file and/or files with ID and POS info")
  }
  
  message("Minimizing LD...")
  snp.locus <- dplyr::select(.data = input, LOCUS, POS) %>% dplyr::distinct(POS, .keep_all = TRUE)
  
  locus.stats <- dplyr::group_by(.data = snp.locus, LOCUS) %>% 
    dplyr::tally(.) %>% dplyr::rename(SNP_N = n) %>% 
    dplyr::group_by(SNP_N) %>%
    dplyr::tally(.) #%>% dplyr::filter(SNP_N != 1)
  
  if (nrow(locus.stats) > 1) {
    range.number.snp.locus <- range(locus.stats$SNP_N, na.rm = TRUE)
    message(stringi::stri_join("The range in the number of SNP/reads or locus is: ", range.number.snp.locus))
  } else {
    message("There is no variation in the number of SNP/reads or locus across the data")
  }
  
  # Random selection ---------------------------------------------------------
  if (snp.ld == "random") {
    snp.select <- snp.locus %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::sample_n(tbl = ., size = 1, replace = FALSE)
    message(stringi::stri_join("Number of original SNP = ", dplyr::n_distinct(snp.locus$POS), "\n", "Number of SNP randomly selected to keep 1 SNP per read/haplotype = ", dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS)))
  }
  
  # Fist SNP on the read -----------------------------------------------------
  if (snp.ld == "first") {
    snp.select <- snp.locus %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::summarise(POS = min(POS))
    message(stringi::stri_join("Number of original SNP = ", dplyr::n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS)))
  }
  
  # Last SNP on the read -----------------------------------------------------
  if (snp.ld == "last") {
    snp.select <- snp.locus %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::summarise(POS = max(POS))
    message(stringi::stri_join("Number of original SNP = ", dplyr::n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS)))
  }
  
  # filtering the VCF to minimize LD -----------------------------------------
  input <- dplyr::semi_join(input, snp.select, by = c("LOCUS", "POS"))
  message("Filtering the tidy VCF to minimize LD by keeping only 1 SNP per short read/haplotype")
  return(input)
}#End snp_ld
