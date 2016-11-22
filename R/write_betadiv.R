# write a betadiv from a tidy data frame

#' @name write_betadiv

#' @title Used internally in stackr to write a betadiv file from a tidy 
#' data frame

#' @description Write a betadiv file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.
#' 
#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info. 

#' @return A betadiv object is returned.

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
#' @rdname write_betadiv
#' @import reshape2
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @references Lamy T, Legendre P, Chancerelle Y, Siu G, Claudet J (2015) 
#' Understanding the Spatio-Temporal Response of Coral Reef Fish Communities to 
#' Natural Disturbances: Insights from Beta-Diversity Decomposition. 
#' PLoS ONE, 10, e0138696.

#' @seealso \code{beta.div} is available on Pierre Legendre web site \url{http://adn.biol.umontreal.ca/~numericalecology/Rcode/} \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}

#' @author Laura Benestan \email{laura.benestan@@icloud.com} and
#' Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_betadiv <- function(data) {
  
  input <- stackr::read_long_tidy_wide(data = data)
  
  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                              pattern = "GENOTYPE", 
                                              replacement = "GT", 
                                              vectorize_all = FALSE)
  }
  
  # Switch colnames LOCUS to MARKERS if found
  if ("LOCUS" %in% colnames(input)) input <- rename(.data = input, MARKERS = LOCUS)
  
  
  # Compute count and Minor Allele Frequency -----------------------------------
  # We split the alleles here to prep for MAF
  # need to compute REF/ALT allele for non VCF file
  if (!tibble::has_name(input, "GT_VCF")) {
    ref.alt.alleles.change <- ref_alt_alleles(data = input)
    input <- left_join(input, ref.alt.alleles.change$input, by = c("MARKERS", "INDIVIDUALS"))
  }
  
  # MAF
  betadiv <- input %>% 
    group_by(MARKERS, POP_ID) %>%
    summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT_VCF[GT_VCF == "0/0"])),
      PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT == "0/1"])),
      QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
    ) %>%
    mutate(MAF = ((QQ*2) + PQ)/(2*N)) %>%
    select(POP_ID, MARKERS, MAF) %>% 
    group_by(POP_ID) %>% 
    tidyr::spread(data = ., key = MARKERS, value = MAF) %>%
    ungroup() %>% 
    mutate(POP_ID = as.integer(POP_ID))
  return(betadiv)
} # end write_betadiv

