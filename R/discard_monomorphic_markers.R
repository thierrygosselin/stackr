# Discard monomorphic markers

#' @name discard_monomorphic_markers

#' @title Discard monomorphic markers

#' @description Discard monomorphic markers.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.
#' 
#' @param data A tidy genomic data set in the working directory or 
#' object in the global environment in wide or long (tidy) formats. 
#' See details for more info. 

#' @return A list with the filtered input file and the blacklist of markers removed.


#' @details \strong{Input data:}
#'  
#' To discriminate the long from the wide format, 
#' the function \pkg{stackr} \code{\link[stackr]{read_long_tidy_wide}} searches 
#' for \code{MARKERS or LOCUS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns: 
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals), 
#' the remaining columns are the markers in separate columns storing genotypes.
#' 
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info. 
#' (e.g. from a VCF see \pkg{stackr} \code{\link{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID}, 
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#' 
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 of 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 of 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#' 
#' \emph{How to get a tidy data frame ?}
#' \pkg{stackr} \code{\link{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame.

#' @export
#' @rdname discard_monomorphic_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter semi_join n_distinct
#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

discard_monomorphic_markers <- function(data) {
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
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
  markers.df <- dplyr::distinct(.data = input, MARKERS, CHROM, LOCUS, POS)
  }
  
  message("Scanning for monomorphic markers...")
  
  mono.markers <- input %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT) %>%
    dplyr::mutate(
      A1 = stringi::stri_sub(GT, 1, 3),
      A2 = stringi::stri_sub(GT, 4,6)
    ) %>% 
    dplyr::select(-GT) %>% 
    tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
    dplyr::group_by(MARKERS, ALLELES) %>% 
    dplyr::tally(.) %>%
    dplyr::ungroup() %>% 
    dplyr::select(MARKERS) %>% 
    dplyr::group_by(MARKERS) %>% 
    dplyr::tally(.) %>% 
    dplyr::filter(n == 1) %>% 
    dplyr::select(MARKERS)
  
  # Remove the markers from the dataset
  message(paste0("Number of monomorphic markers removed = ", dplyr::n_distinct(mono.markers$MARKERS)))
  
  if (length(mono.markers$MARKERS) > 0) {
    input <- dplyr::anti_join(input, mono.markers, by = "MARKERS")
  }
  
  if (tibble::has_name(input, "CHROM")) {
    mono.markers <- dplyr::left_join(mono.markers, markers.df, by = "MARKERS")
  }
  
  res <- list(input = input, blacklist.momorphic.markers = mono.markers)
  return(res)
} # end discar mono markers

