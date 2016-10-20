# Detect if markers are biallelic

#' @name detect_biallelic_markers

#' @title detect if markers in tidy data set are biallelic

#' @description Detect if markers in tidy data set are biallelic.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.
#' 
#' @param data A tidy genomic data set in the working directory or 
#' object in the global environment in wide or long (tidy) formats. 
#' See details for more info. 

#' @return A logical character string (TRUE/FALSE). That answer the question if
#' the data set is biallelic or not.

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
#' @rdname detect_biallelic_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter
#' @importFrom stringi stri_replace_all_fixed stri_sub
#' @importFrom tibble has_name
#' @importFrom tidyr gather
#' @importFrom purrr flatten_chr

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_biallelic_markers <- function(data) {
  
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
  
  # Detecting biallelic markers-------------------------------------------------
  input <- input %>%
    dplyr::filter(GT != "000000") %>% 
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
    dplyr::mutate(
      A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
      A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
    ) %>% 
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
    tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% 
    dplyr::group_by(MARKERS, ALLELES) %>% 
    dplyr::tally(.) %>%
    dplyr::ungroup() %>% 
    dplyr::select(MARKERS) %>% 
    dplyr::group_by(MARKERS) %>% 
    dplyr::tally(.) %>% 
    dplyr::summarise(BIALLELIC = max(n, na.rm = TRUE)) %>%
    purrr::flatten_chr(.x = .)
  
  if (input != 2) {
    biallelic <- FALSE
  } else {
    biallelic <- TRUE
  }
  return(biallelic)
} # End detect_biallelic_markers
