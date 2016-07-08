# write a structure file from a tidy data frame

#' @name write_structure
#' @title Used internally in stackr to write a structure file from a tidy data frame
#' @description Write a structure file from a tidy data frame
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info.

#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}.

#' @param markers.line (optional, logical) In the genepop and structure
#' file, you can write the markers on a single line separated by 
#' commas \code{markers.line = TRUE}, 
#' or have markers on a separate line, i.e. in one column, for the genepop file
#' (not very useful with thousands of markers) and not printed at all for the
#' structure file.
#' Default: \code{markers.line = TRUE}.


#' @param filename (optional) The file name prefix for the genepop file 
#' written to the working directory. With default: \code{filename = NULL}, 
#' the date and time is appended to \code{stackr_structure_}.

#' @param ... other parameters passed to the function.

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

#' @return A structure file is saved to the working directory. 

#' @export
#' @rdname write_structure
#' @import reshape2
#' @import dplyr
#' @import stringi

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# to get rid of notes in build check
# if(getRversion() >= "2.15.1") {
#   utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", 
#                            "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", 
#                            "SAMPLES", "ALLELES", 'A1', 'A2', 'COUNT', 
#                            "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", 
#                            "POLYMORPHISM", "POLYMORPHISM_MAX", "other", 
#                            "strata", "hierarchy", "GROUP", ".", 'MARKERS', 
#                            'MARKERS_ALLELES', 'STRATA'
#   )
#   )
# }


write_structure <- function(
  data,
  pop.levels = NULL, 
  markers.line = TRUE, 
  filename = NULL,
  ...
  ) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the structure file is missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- stackr::read_long_tidy_wide(data = data)
  } else {
    input <- data
  }
  
  colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                            pattern = "GENOTYPE", 
                                            replacement = "GT", 
                                            vectorize_all = FALSE)
  
  # Switch colnames LOCUS to MARKERS if found
  if ("LOCUS" %in% colnames(input)) input <- rename(.data = input, MARKERS = LOCUS)
  
  input <- input %>% 
    select(POP_ID, INDIVIDUALS, MARKERS, GT)
  
  # pop.levels -----------------------------------------------------------------
  if (!is.null(pop.levels)) {
    input <- input %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =TRUE),
        POP_ID = droplevels(POP_ID)
        ) %>% 
      arrange(POP_ID, INDIVIDUALS, MARKERS)
  } else {
    input <- input %>% 
      mutate(POP_ID = factor(POP_ID)) %>% 
      arrange(POP_ID, INDIVIDUALS, MARKERS)
  }
  
  # Create a marker vector  ------------------------------------------------
  markers <- input %>% distinct(MARKERS) %>% arrange(MARKERS)
  markers <- markers$MARKERS
  
  # Structure format ----------------------------------------------------------------
  input <- input %>%
    tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3, extra = "drop", remove = TRUE) %>%
    tidyr::gather(data = ., key = ALLELES, value = GT, -c(POP_ID, INDIVIDUALS, MARKERS)) %>% 
    mutate(
      GT = stri_replace_all_fixed(str = GT, pattern = "000", replacement = "-9", vectorize_all = FALSE),
      GT = as.numeric(GT)
    ) %>%
    select(INDIVIDUALS, POP_ID, MARKERS, ALLELES, GT) %>% 
    tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
    mutate(POP_ID = as.numeric(POP_ID)) %>% 
    select(-ALLELES) %>% 
    arrange(POP_ID, INDIVIDUALS)

  # Write the file in structure format -----------------------------------------
  
  # Filename
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename <- stri_paste("stackr_structure_", file.date, ".gen")
  } else {
    filename <- stri_paste(filename, ".str")
  }
  
  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = stri_paste(markers, sep = "\t", collapse = "\t"), con = filename.connection, sep = "\n") 
  close(filename.connection) # close the connection
  write_tsv(x = input, path = filename, append = TRUE, col_names = FALSE)
} # end write_structure
