# Compute the ref and alt alleles of a tidy dataset

#' @name ref_alt_alleles

#' @title Used internally in stackr to compute the REF and ALT alleles from a 
#' genomic tidy data frame

#' @description compute the REF and ALT alleles from a 
#' genomic tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.
#' 
#' @param data A genomic data set in the working directory or 
#' object in the global environment in wide or long (tidy) formats. 
#' See details for more info. 

#' @return A tidy data frame with 6 columns: 
#' \code{MARKERS, INDIVIDUALS, REF, ALT, GT_VCF, GT_BIN}.
#' \code{GT_VCF}: the genotype in VCF format
#' \code{GT_BIN}: coding used internally to easily convert to genlight, 
#' the coding \code{0, 1, 2, NA} stands for the number of ALT allele in the 
#' genotype and \code{NA} for missing genotype.

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
#' @rdname ref_alt_alleles
#' @import reshape2
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


ref_alt_alleles <- function (data) {
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                              pattern = "GENOTYPE", 
                                              replacement = "GT", 
                                              vectorize_all = FALSE)
  }
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (!tibble::has_name(input, "MARKERS") && tibble::has_name(input, "LOCUS")) {
    input <- rename(.data = input, MARKERS = LOCUS)
  }
  
  input.select <- input %>%
    select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
    #faster than: tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
    mutate(
      A1 = stri_sub(str = GT, from = 1, to = 3),
      A2 = stri_sub(str = GT, from = 4, to = 6)
    ) %>%
    select(-GT)
  
  new.ref.alt.alleles <- tidyr::gather(
    data = input.select, 
    key = ALLELES, 
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
  return(ref.alt.alleles.change)
}
