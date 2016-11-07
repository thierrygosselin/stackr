# Keep markers in common between all populations

#' @name keep_common_markers

#' @title Keep markers in common between all populations

#' @description Keep markers in common between all populations.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.
#' 
#' @param data A tidy genomic data set in the working directory or 
#' object in the global environment in wide or long (tidy) formats. 
#' See details for more info. 

#' @return The filtered input data set with markers found in all the populations
#' present in the data set.


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
#' @rdname keep_common_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter semi_join n_distinct
#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

keep_common_markers <- function(data) {
  
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
  
  message("Using markers common in all populations:")
  pop.number <- dplyr::n_distinct(input$POP_ID)
  
  pop.filter <- dplyr::filter(input, GT != "000000")
  
  pop.filter <- pop.filter %>% 
    dplyr::group_by(MARKERS) %>%
    dplyr::filter(dplyr::n_distinct(POP_ID) == pop.number) %>%
    dplyr::arrange(MARKERS) %>%
    dplyr::distinct(MARKERS)
  
  markers.input <- dplyr::n_distinct(input$MARKERS)
  markers.in.common <- dplyr::n_distinct(pop.filter$MARKERS)
  blacklist.markers.common <- markers.input - markers.in.common
  
  message(
    stringi::stri_join(
      "Number of original markers = ", markers.input,
      "\n", "Number of markers removed = ", blacklist.markers.common, "\n",
      "Number of markers present in all the populations = ",
      markers.in.common
    )
  )
  
  if (blacklist.markers.common > 0) {
    input <- suppressWarnings(input %>% dplyr::semi_join(pop.filter, by = "MARKERS"))
  }
  return(input)
}
