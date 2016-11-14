# write a genind file from a tidy data frame

#' @name write_genind
#' @title Used internally in stackr to write a genind object from a tidy data frame
#' @description Write a genind object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.
#' 
#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info. 

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
#' @rdname write_genind
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom data.table fread as.data.table dcast.data.table
#' @importFrom tidyr spread gather separate complete
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_replace_na 
#' @importFrom tibble has_name as_data_frame

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1: 
#' new tools for the analysis of genome-wide SNP data. 
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genind <- function(data) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the genepop file is missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input), 
    pattern = "GENOTYPE", 
    replacement = "GT", 
    vectorize_all = FALSE
  )
  
  strata.genind <- dplyr::distinct(.data = input, INDIVIDUALS, POP_ID)
  
  # When VCF data available
  if (tibble::has_name(input, "GT_VCF")) {
    genind.prep <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT_VCF) %>%
      dplyr::mutate(
        A1_A2 = stringi::stri_replace_all_fixed(
          str = GT_VCF, 
          pattern = c("0/0", "1/1", "0/1", "1/0", "./."), 
          replacement = c("2_0", "0_2", "1_1", "1_1", NA),
          vectorize_all = FALSE
        )
      ) %>% 
      dplyr::mutate(POP_ID = factor(as.character(POP_ID))) %>%# xvalDapc doesn't accept pop as ordered factor
      dplyr::mutate(
        A1 = stringi::stri_sub(str = A1_A2, from = 1, to = 1),
        A2 = stringi::stri_sub(str = A1_A2, from = 3, to = 3)
      ) %>% 
      dplyr::select(-GT_VCF, -A1_A2) %>% 
      tidyr::gather(data = ., key = ALLELES, value = n, -c(INDIVIDUALS, POP_ID, MARKERS)) %>% 
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>% 
      dplyr::select(-MARKERS, -ALLELES)

    genind.prep <- data.table::dcast.data.table(
      data = data.table::as.data.table(genind.prep), 
      formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES, 
      value.var = "n") %>% 
      tibble::as_data_frame()
  } else {
    missing.geno <- dplyr::ungroup(input) %>%
      dplyr::filter(GT == "000000") %>%
      dplyr::select(MARKERS, INDIVIDUALS) %>% 
      dplyr::mutate(MISSING = rep("blacklist", n()))
    
    genind.prep <- input %>%
      dplyr::filter(GT != "000000") %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>% 
      dplyr::mutate(
        A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
        A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
      ) %>% 
      dplyr::select(-GT) %>% 
      tidyr::gather(
        data = .,
        key = ALLELES,
        value = GT, 
        -c(MARKERS, INDIVIDUALS)
      ) %>% 
      dplyr::arrange(MARKERS, INDIVIDUALS, GT) %>%
      dplyr::group_by(INDIVIDUALS, MARKERS, GT) %>% 
      dplyr::tally(.) %>% # count alleles, longest part of the block
      dplyr::ungroup(.) %>%
      tidyr::complete(data = ., INDIVIDUALS, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
      dplyr::anti_join(missing.geno, by = c("MARKERS", "INDIVIDUALS")) %>% 
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, GT, sep = ".")) %>%
      dplyr::select(-MARKERS, -GT) %>%
      dplyr::right_join(strata.genind, by = "INDIVIDUALS") %>%#include strata
      dplyr::mutate(POP_ID = factor(as.character(POP_ID))) %>%# xvalDapc doesn't accept pop as ordered factor
      dplyr::arrange(MARKERS_ALLELES, INDIVIDUALS)
    
    missing.geno <- NULL
    
    genind.prep <- data.table::dcast.data.table(
      data = data.table::as.data.table(genind.prep), 
      formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES, 
      value.var = "n") %>% 
      tibble::as_data_frame()
  }
  
  # genind arguments common to all data.type
  ind <- genind.prep$INDIVIDUALS
  pop <- genind.prep$POP_ID
  genind.df <-  dplyr::ungroup(genind.prep) %>% dplyr::select(-c(INDIVIDUALS, POP_ID))
  suppressWarnings(rownames(genind.df) <- ind)
  loc.names <- colnames(genind.df)
  # strata <- dplyr::ungroup(genind.prep) %>% dplyr::distinct(INDIVIDUALS, POP_ID)
  
  # genind constructor
  prevcall <- match.call()
  res <- adegenet::genind(
    tab = genind.df,
    pop = pop,
    prevcall = prevcall,
    ploidy = 2,
    type = "codom",
    strata = strata.genind,
    hierarchy = NULL
  )
  
  return(res)
} # End write_genind
