# write a genind file from a tidy data frame

#' @name write_genlight
#' @title Used internally in stackr to write a genlight object from a tidy data frame
#' @description Write a genlight object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
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
#' @rdname write_genlight
#' @import reshape2
#' @import dplyr
#' @import stringi
#' @importFrom methods new
#' @importFrom data.table fread

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1: 
#' new tools for the analysis of genome-wide SNP data. 
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genlight <- function(data) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the genepop file is missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
  colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                            pattern = "GENOTYPE", 
                                            replacement = "GT", 
                                            vectorize_all = FALSE)
  
  # data = input
  marker.meta <- distinct(.data = input, MARKERS, CHROM, LOCUS, POS)
  
  if (!tibble::has_name(input, "GT_BIN")) {
    input$GT_BIN <- stri_replace_all_fixed(
      str = input$GT_VCF, 
      pattern = c("0/0", "1/1", "0/1", "1/0", "./."), 
      replacement = c("0", "2", "1", "1", NA), 
      vectorize_all = FALSE
      )
  }
  
  
  geno.df <- select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT_BIN) %>% 
    mutate(GT_BIN = as.integer(GT_BIN)) %>% 
    arrange(MARKERS)
  
  geno.df <- data.table::dcast.data.table(
    data = as.data.table(geno.df), 
    formula = INDIVIDUALS + POP_ID ~ MARKERS, 
    value.var = "GT_BIN") %>% 
    as_data_frame() %>% 
    arrange(POP_ID, INDIVIDUALS)
  
  # Generate genlight
  genlight.object <- methods::new("genlight", geno.df[,-(1:2)], parallel=FALSE)
  adegenet::indNames(genlight.object) <- geno.df$INDIVIDUALS
  adegenet::pop(genlight.object) <- geno.df$POP_ID
  adegenet::chromosome(genlight.object) <- marker.meta$CHROM
  adegenet::locNames(genlight.object)   <- marker.meta$LOCUS
  adegenet::position(genlight.object)   <- marker.meta$POS
  
  
  # test
  # genlight.object@n.loc
  # genlight.object@ind.names
  # genlight.object@chromosome
  # genlight.object@position
  # genlight.object@loc.names
  # genlight.object@pop
  # genlight.object@strata
  # adegenet::nLoc(genlight.object)
  # adegenet::popNames(genlight.object)
  # adegenet::indNames(genlight.object)
  # adegenet::nPop(genlight.object)
  # adegenet::NA.posi(genlight.object)
  
  
  return(genlight.object)
} # End write_genlight
