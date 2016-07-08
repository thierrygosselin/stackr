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
  
  if (is.vector(data)) {
    # input <- structure.prep # test
    input <- stackr::read_long_tidy_wide(data = data)
  } else {
    # data <- structure.prep # test
    input <- data
  }
  
  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                              pattern = "GENOTYPE", 
                                              replacement = "GT", 
                                              vectorize_all = FALSE)
  }
  
  # Compute count and Minor Allele Frequency -----------------------------------
  # We split the alleles here to prep for MAF
  # need to compute REF/ALT allele for non VCF file
  if (!tibble::has_name(input, "GT_VCF")) {
    input.select <- input %>%
      select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
      #faster than: tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
      mutate(
        A1 = stri_sub(str = GT, from = 1, to = 3),
        A2 = stri_sub(str = GT, from = 4, to = 6)
      ) %>%
      select(-GT)
    
    new.ref.alt.alleles <- input.select %>% 
      tidyr::gather(
        data = ., key = ALLELES, 
        value = GT, 
        -c(MARKERS, INDIVIDUALS, POP_ID)
      ) %>% # just miliseconds longer than data.table.melt so keeping this one for simplicity
      filter(GT != "000") %>%  # remove missing "000"
      group_by(MARKERS, GT) %>%
      tally %>% 
      ungroup() %>% 
      mutate(ALLELE_NUMBER = rep(c("A1", "A2"), each = 1, times = n()/2)) %>% 
      group_by(MARKERS) %>%
      mutate(
        PROBLEM = ifelse(n[ALLELE_NUMBER == "A1"] == n[ALLELE_NUMBER == "A2"], "equal_number", "ok"),
        GROUP = ifelse(n == max(n), "REF", "ALT")
      ) %>% 
      group_by(MARKERS, GT) %>% 
      mutate(
        ALLELE = ifelse(PROBLEM == "equal_number" & ALLELE_NUMBER == "A1", "REF", 
                        ifelse(PROBLEM == "equal_number" & ALLELE_NUMBER == "A2", "ALT", 
                               GROUP)
        )
      ) %>% 
      select(MARKERS, GT, ALLELE) %>%
      group_by(MARKERS) %>% 
      tidyr::spread(data = ., ALLELE, GT) %>% 
      select(MARKERS, REF, ALT)
    
    ref.alt.alleles.change <- full_join(input.select, new.ref.alt.alleles, by = "MARKERS") %>% 
      mutate(
        A1 = replace(A1, which(A1 == "000"), NA),
        A2 = replace(A2, which(A2 == "000"), NA),
        GT_VCF_A1 = if_else(A1 == REF, "0", "1", missing = "."),
        GT_VCF_A2 = if_else(A2 == REF, "0", "1", missing = "."),
        GT_VCF = stri_paste(GT_VCF_A1, GT_VCF_A2, sep = "/"),
        GT_BIN = stri_replace_all_fixed(
          str = GT_VCF, 
          pattern = c("0/0", "1/1", "0/1", "1/0", "./."), 
          replacement = c("0", "2", "1", "1", NA), 
          vectorize_all = FALSE
        )
      ) %>% 
      select(MARKERS, INDIVIDUALS, REF, ALT, GT_VCF, GT_BIN)
    
    input <- left_join(input, ref.alt.alleles.change, by = c("MARKERS", "INDIVIDUALS"))
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

