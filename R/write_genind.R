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
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join summarise_each_ funs
#' @importFrom data.table fread as.data.table
#' @importFrom tidyr spread gather separate
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_replace_na 


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
  # Switch colnames LOCUS to MARKERS if found
  # if ("LOCUS" %in% colnames(input)) input <- rename(.data = input, MARKERS = LOCUS)
  
  genind.prep <- input %>% 
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
    #faster than: tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
    dplyr::mutate(
      A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
      A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
    ) %>% 
    dplyr::select(-GT) %>% 
    tidyr::gather(
      data = .,
      key = ALLELES,
      value = GT, 
      -c(MARKERS, INDIVIDUALS, POP_ID)
    ) %>% # just miliseconds longer than data.table.melt so keeping this one for simplicity
    dplyr::filter(GT != "000") # remove missing "000"
  
  # this reintroduce the missing, but with NA
  genind.prep <- data.table::dcast.data.table(
    data = data.table::as.data.table(genind.prep), 
    formula = POP_ID + INDIVIDUALS + ALLELES ~ MARKERS, 
    value.var = "GT") %>% 
    tibble::as_data_frame() %>%
    dplyr::mutate_all(.tbl = ., .funs = factor, exclude = NA) %>% 
    dplyr::mutate(INDIVIDUALS = as.character(INDIVIDUALS))
  
  # The next part is longer than it used to be with VCF file only, 
  # but it as the advantage of working and simplifying the use for other file type.
  genind.prep <- suppressWarnings(dplyr::mutate_each(tbl = genind.prep, dplyr::funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)))
  
  genind.prep <- tidyr::gather(
    data = genind.prep, 
    key = MARKERS, 
    value = GT, 
    -c(INDIVIDUALS, POP_ID, ALLELES)
  ) %>% # faster than data.table.melt...
    dplyr::mutate(GT = stringi::stri_replace_na(str = GT, replacement = "000")) %>%
    dplyr::filter(GT != "000") %>%
    dplyr::select(-ALLELES) %>%
    dplyr::group_by(POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
    dplyr::tally(.) %>% # count alleles, longest part of the block
    dplyr::ungroup(.)
  
  genind.prep <- genind.prep %>%
    dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, GT, sep = ":")) %>%  # faster then: tidyr::unite(MARKERS_ALLELES, MARKERS, GT, sep = ":", remove = TRUE)
    dplyr::select(-GT, -MARKERS) %>% 
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES)
  
  genind.prep <- data.table::dcast.data.table(
    data = data.table::as.data.table(genind.prep), 
    formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES, 
    value.var = "n") %>% 
    tibble::as_data_frame()
  
  genind.prep <- tidyr::gather(data = genind.prep, key = MARKERS_ALLELES, value = COUNT, -c(INDIVIDUALS, POP_ID)) %>% 
    tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = ":", remove = TRUE) %>% 
    dplyr::mutate(COUNT = as.numeric(stringi::stri_replace_na(str = COUNT, replacement = "0"))) %>% 
    dplyr::group_by(INDIVIDUALS, MARKERS) %>%
    dplyr::mutate(MAX_COUNT_MARKERS = max(COUNT, na.rm = TRUE)) %>%
    dplyr::ungroup(.) %>% 
    dplyr::mutate(COUNT = ifelse(MAX_COUNT_MARKERS == 0, "erase", COUNT)) %>%
    dplyr::select(-MAX_COUNT_MARKERS) %>% 
    dplyr::mutate(COUNT = replace(COUNT, which(COUNT == "erase"), NA)) %>% 
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES)
  
  genind.prep <- genind.prep %>%
    dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%  # faster then: tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE)
    dplyr::select(-MARKERS, -ALLELES) %>% 
    dplyr::mutate(
      POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
      POP_ID = factor(POP_ID) # xvalDapc doesn't accept pop as ordered factor
    )
  
  genind.prep <- data.table::dcast.data.table(
    data = data.table::as.data.table(genind.prep), 
    formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES, 
    value.var = "COUNT") %>% 
    tibble::as_data_frame() %>%
    dplyr::arrange(POP_ID, INDIVIDUALS)
  
  # genind arguments common to all data.type
  ind <- genind.prep$INDIVIDUALS
  pop <- genind.prep$POP_ID
  genind.df <-  dplyr::ungroup(genind.prep) %>% dplyr::select(-c(INDIVIDUALS, POP_ID))
  suppressWarnings(rownames(genind.df) <- ind)
  loc.names <- colnames(genind.df)
  strata <- dplyr::ungroup(genind.prep) %>% dplyr::distinct(INDIVIDUALS, POP_ID)
  
  # genind constructor
  prevcall <- match.call()
  res <- adegenet::genind(
    tab = genind.df,
    pop = pop,
    prevcall = prevcall,
    ploidy = 2,
    type = "codom",
    strata = strata,
    hierarchy = NULL
  )
  
  return(res)
} # End write_genind
