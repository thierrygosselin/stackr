# Keep markers in common between all populations

#' @name keep_common_markers

#' @title Keep markers in common between all populations

#' @description Keep markers in common between all populations.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.


#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{stackr} \code{\link{tidy_genomic_data}}.


#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty 
#' during execution. 
#' Default: \code{verbose = FALSE}.

#' @return The filtered input data set with markers found in all the populations
#' present in the data set.



#' @export
#' @rdname keep_common_markers
#' @importFrom dplyr select mutate group_by ungroup rename tally filter semi_join n_distinct
#' @importFrom stringi stri_replace_all_fixed stri_join
#' @importFrom tibble has_name

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

keep_common_markers <- function(data, verbose = FALSE) {
  
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
  
  if (verbose) message("Using markers common in all populations:")
  pop.filter <- dplyr::select(.data = input, MARKERS, POP_ID, GT) %>% 
    dplyr::filter(GT != "000000") %>%
    dplyr::distinct(MARKERS, POP_ID) %>%
    dplyr::count(x = ., MARKERS) %>% 
    dplyr::filter(n == dplyr::n_distinct(input$POP_ID)) %>%
    dplyr::distinct(MARKERS) %>% 
    dplyr::arrange(MARKERS)
  
  markers.input <- dplyr::n_distinct(input$MARKERS)
  markers.in.common <- nrow(pop.filter)
  blacklist.markers.common <- markers.input - markers.in.common
  
  if (verbose) message(stringi::stri_join("    Number of markers before = ", markers.input))
  if (verbose) message(stringi::stri_join("    Number of markers removed = ", blacklist.markers.common))
  if (verbose) message(stringi::stri_join("    Number of markers after (common between populations) = ", markers.in.common))

  
  if (blacklist.markers.common > 0) {
    input <- suppressWarnings(dplyr::semi_join(input, pop.filter, by = "MARKERS"))
  }
  return(input)
}
