# Detect if markers are biallelic

#' @name allele_frequencies

#' @title Compute allele frequencies per markers and populations
#' @description Compute allele frequencies per markers and populations.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{stackr} \code{\link{tidy_genomic_data}}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty 
#' during execution.
#' Default: \code{verbose = TRUE}.

#' @return A list with allele frequencies in a data frame in long and wide format,
#' and a matrix.

#' @export
#' @rdname allele_frequencies
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join summarise_each_ funs
#' @importFrom stringi stri_replace_all_fixed stri_sub stri_join
#' @importFrom tibble has_name
#' @importFrom tidyr gather spread complete nesting

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

allele_frequencies <- function(data, verbose = TRUE) {
  
  if (verbose) {
    cat("#######################################################################\n")
    cat("##################### stackr::allele_frequencies ######################\n")
    cat("#######################################################################\n")
    timing <- proc.time()
    message("Calculating allele frequencies...")
  }
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  
  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- stackr::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }
  
  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input), 
      pattern = "GENOTYPE", 
      replacement = "GT", 
      vectorize_all = FALSE)
  }
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  if (tibble::has_name(input, "CHROM")) {
    metadata.markers <- dplyr::distinct(input, MARKERS, CHROM, LOCUS, POS)
  }
  
  input <- input %>%
    dplyr::filter(GT != "000000") %>% 
    dplyr::mutate(
      A1 = stringi::stri_sub(GT, 1, 3),
      A2 = stringi::stri_sub(GT, 4,6)
    ) %>% 
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>% 
    tidyr::gather(key = ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, POP_ID, MARKERS)) %>%
    dplyr::group_by(MARKERS, ALLELES, POP_ID) %>% 
    dplyr::tally(.) %>%
    dplyr::ungroup(.) %>%
    tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, ALLELES), fill = list(n = 0)) %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    dplyr::mutate(
      NAPL = sum(n),
      FREQ_APL = n / NAPL # Frequency of alleles per pop and locus
    ) %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::mutate(FREQ = sum(n) / sum(NAPL)) %>% #Frequency of alleles per locus
    dplyr::select(MARKERS, POP_ID, ALLELES, FREQ) %>% 
    dplyr::arrange(MARKERS, POP_ID, ALLELES)
  
  if (tibble::has_name(input, "CHROM")) {
    input <- dplyr::full_join(input, metadata.markers, by = "MARKERS") %>% 
      dplyr::select(MARKERS, CHROM, LOCUS, POS, POP_ID, ALLELES, FREQ)
  }
  
  
  freq.wide <- dplyr::ungroup(input) %>% 
    dplyr::select(MARKERS, ALLELES, POP_ID, FREQ) %>%
    dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
    dplyr::select(-MARKERS, -ALLELES) %>% 
    dplyr::arrange(MARKERS_ALLELES, POP_ID) %>%
    dplyr::group_by(POP_ID) %>% 
    tidyr::spread(data = ., key = MARKERS_ALLELES, value = FREQ)
  
  freq.mat <- suppressWarnings(
    freq.wide %>%
      tibble::remove_rownames(df = .) %>% 
      tibble::column_to_rownames(df = ., var = "POP_ID") %>% 
      as.matrix(.)
  )  
  if (verbose) {
    message(stringi::stri_join("Computation time: ", round((proc.time() - timing)[[3]]), " sec"))
    cat("############################## completed ##############################\n")
  }
  res <- list(
    freq.long = input,
    freq.wide = freq.wide,
    mat = freq.mat
  )
  return(res)
}
