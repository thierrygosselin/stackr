# write a plinkfile from a tidy data frame

#' @name write_plink
#' @title Used internally in stackr to write a plink tped/tfam file from a tidy 
#' data frame
#' @description Write a plink file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.
#' 
#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info. 

#' @param filename (optional) The file name prefix for tped/tfam files 
#' written to the working directory. With default: \code{filename = NULL}, 
#' the date and time is appended to \code{stackr_plink_}.


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
#' @rdname write_plink
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, et al. 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. 
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_plink <- function(data, filename = NULL) {
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
  colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                            pattern = "GENOTYPE", 
                                            replacement = "GT", 
                                            vectorize_all = FALSE)
  
  
  
  tped <- input %>% 
    arrange(INDIVIDUALS) %>% 
    mutate(
      COL1 = rep("0", n()),
      COL3 = rep("0", n()),
      COL4 = rep("0", n())
    ) %>% 
    select(COL1, MARKERS, COL3, COL4, INDIVIDUALS, GT) %>% 
    mutate(
      A1 = stri_sub(str = GT, from = 1, to = 3),
      A2 = stri_sub(str = GT, from = 4, to = 6)
    ) %>% 
    select(-GT) %>% 
    tidyr::gather(ALLELES, GENOTYPE, -c(COL1, MARKERS, COL3, COL4, INDIVIDUALS)) %>%
    mutate(
      GENOTYPE = as.character(as.numeric(GENOTYPE)),
      GENOTYPE = stri_pad_left(GENOTYPE, width = 2, pad = "0")
    ) %>%  
    arrange(INDIVIDUALS, ALLELES) %>% 
    tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_") %>%
    group_by(COL1, MARKERS, COL3, COL4) %>% 
    tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GENOTYPE) %>% 
    arrange(MARKERS)
  
  tfam <- input %>%
    distinct(POP_ID, INDIVIDUALS) %>% 
    arrange(INDIVIDUALS) %>% 
    mutate(
      COL3 = rep("0",n()),
      COL4 = rep("0",n()),
      COL5 = rep("0",n()),
      COL6 = rep("-9",n())
    )
  
  # Create a filename to save the output files ********************************
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename.tped <- stri_paste("stackr_plink_", file.date, ".tped")
    filename.tfam <- stri_paste("stackr_plink_", file.date, ".tfam")
  } else {
    filename.tped <- stri_paste(filename, ".tped")
    filename.tfam <- stri_paste(filename, ".tfam")
  }
  write_delim(x = tped, path = filename.tped, col_names = FALSE, delim = " ")
  write_delim(x = tfam, path = filename.tfam, col_names = FALSE, delim = " ")
} # end write_plink
