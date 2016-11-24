# write a hierfstat file from a tidy data frame

#' @name write_hierfstat
#' @title Used internally in stackr to write a hierfstat file from a tidy data frame
#' @description Write a hierfstat file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.
#' 
#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info. 

#' @param filename (optional) The file name prefix for the hierfstat file 
#' written to the working directory. With default: \code{filename = NULL}, 
#' the date and time is appended to \code{stackr_hierfstat_}.

#' @return A hierfstat file is saved to the working directory. 

#' @details \strong{Details for the sep argument}
#' This character is directly used in regular expressions using strigi. 
#' Some characters need to be preceeded by double backslashes \code{\\}. 
#' For instance, "/" works but "|" must be coded as "\\|".
#' 
#' 
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
#' @rdname write_hierfstat
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @references Goudet, J. (1995) FSTAT (Version 1.2): A computer program to 
#' calculate F- statistics. Journal of Heredity, 86, 485-486.
#' @references Goudet, J. (2005) hierfstat, a package for r to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5, 184-186.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_hierfstat <- function(data, filename = NULL) {
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file necessary to write the hierfstat file is missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data)
  
  colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                            pattern = "GENOTYPE", 
                                            replacement = "GT", 
                                            vectorize_all = FALSE)
  
  # Switch colnames LOCUS to MARKERS if found
  if ("LOCUS" %in% colnames(input)) input <- rename(.data = input, MARKERS = LOCUS)
  
  input <- select(.data = input, POP_ID, INDIVIDUALS, MARKERS, GT) %>% arrange(MARKERS, POP_ID, INDIVIDUALS)
  
  # Create a marker vector  ----------------------------------------------------
  markers <- unique(input$MARKERS)
  
  # Get the number of sample (pop) for hierfstat -------------------------------
  np <- nlevels(droplevels(input$POP_ID))
  np.message <- stri_paste("    * Number of sample pop, np = ", np, sep = "")
  message(np.message)
  
  # Get the number of loci -----------------------------------------------------
  nl <- length(markers)
  nl.message <- stri_paste("    * Number of markers, nl = ", nl, sep = "")
  message(nl.message)
  
  input <- input %>%
    select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
    mutate(
      A1 = as.numeric(stri_sub(str = GT, from = 1, to = 3)),
      A2 = as.numeric(stri_sub(str = GT, from = 4, to = 6))
    ) %>%
    select(-GT)
  
  # Get the highest number used to label an allele -----------------------------
  nu <- max(c(unique(input$A1), unique(input$A2)), na.rm = TRUE)
  nu.message <- stri_paste("    * The highest number used to label an allele, nu = ", 
                           nu, sep = "")
  message(nu.message)
  
  # prep the data  -------------------------------------------------------------
  input <- suppressWarnings(
    input %>% 
      tidyr::unite(GT, A1, A2, sep = "") %>%
      group_by(POP_ID, INDIVIDUALS) %>% 
      tidyr::spread(data = ., MARKERS, GT) %>%
      ungroup %>% 
      arrange(POP_ID, INDIVIDUALS) %>% 
      mutate(POP_ID = as.integer(POP_ID)) %>% 
      select(-INDIVIDUALS)
  )
  
  # allele coding --------------------------------------------------------------
  allele.coding <- 1
  message("    * The alleles are encoded with one digit number")
  
  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename <- stri_paste("stackr_hierfstat_", file.date, ".dat")
  } else {
    filename <- stri_paste(filename, ".dat")
  }
  
  
  # FSTAT: write the first line ------------------------------------------------
  fstat.first.line <- stri_paste(np, nl, nu, allele.coding, sep = " ")
  fstat.first.line <- as.data.frame(fstat.first.line)
  write_delim(x = fstat.first.line, path = filename, delim = "\n", append = FALSE, 
              col_names = FALSE)
  
  # FSTAT: write the locus name to the file
  loci.table <- as.data.frame(markers)
  write_delim(x = loci.table, path = filename, delim = "\n", append = TRUE, 
              col_names = FALSE)
  
  # FSTAT: write the pop and genotypes
  write_delim(x = input, path = filename, delim = "\t", append = TRUE, 
              col_names = FALSE)
  return(input)
}# End write_hierfstat
