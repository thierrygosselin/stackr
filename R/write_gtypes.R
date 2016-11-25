# write a strataG gtypes object from a tidy data frame

#' @name write_gtypes
#' @title Tidy genomic data to \code{\link[strataG]{gtypes}} 
#' @description Write a \code{\link[strataG]{gtypes}} object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.
#' 
#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info. 

#' @return An object of the class \code{\link[strataG]{gtypes}} is returned.

#' @details
#' \strong{Details for Input data:}
#' 
#' To discriminate the long from the wide format, 
#' the function \pkg{stackr} \code{\link[stackr]{read_long_tidy_wide}} searches 
#' for "MARKERS" in column names (TRUE = long format).
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
#' \code{MARKERS} and \code{GENOTYPE or GT}. The rest are considered metata info.
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
#' @rdname write_gtypes
#' @import strataG
#' @importFrom tidyr gather unite spread 
#' @importFrom data.table fread
#' @importFrom methods new
#' @importFrom stringi stri_replace_all_fixed stri_sub
#' @importFrom dplyr select arrange rename mutate


#' @seealso \code{strataG.devel} is available on github \url{https://github.com/EricArcher/}

#' @references Eric Archer, Paula Adams and Brita Schneiders (2016). 
#' strataG: Summaries and Population Structure Analyses of
#' Genetic Data. R package version 1.0.5. https://CRAN.R-project.org/package=strataG

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_gtypes <- function(data) {
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file necessary to write the hierfstat file is missing")
  
  # Import data ---------------------------------------------------------------
  input <- data
  
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input), 
    pattern = "GENOTYPE", 
    replacement = "GT", 
    vectorize_all = FALSE)
  
  # Switch colnames LOCUS to MARKERS if found
  if ("LOCUS" %in% colnames(input)) input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  
  input <- input %>% 
    dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS) %>% 
    dplyr::mutate(
      GT = replace(GT, which(GT == "000000"), NA),
      POP_ID = as.character(POP_ID),
      `1` = stringi::stri_sub(str = GT, from = 1, to = 3), # most of the time: faster than tidyr::separate
      `2` = stringi::stri_sub(str = GT, from = 4, to = 6)
    ) %>%
    dplyr::select(-GT) %>% 
    tidyr::gather(
      data = ., 
      key = ALLELES, 
      value = GT, 
      -c(MARKERS, INDIVIDUALS, POP_ID)
    ) %>% 
    tidyr::unite(data = ., MARKERS_ALLELES, MARKERS, ALLELES, sep = ".") %>% 
    tidyr::spread(data = ., key = MARKERS_ALLELES, value = GT)

  res <- suppressWarnings(
    methods::new("gtypes", gen.data = input[, -(1:2)], ploidy = 2, ind.names = input$INDIVIDUALS, 
                 strata = input$POP_ID, schemes = NULL, sequences = NULL, 
                 description = NULL, other = NULL)
  )      
  return(res)
}# End write_gtypes
