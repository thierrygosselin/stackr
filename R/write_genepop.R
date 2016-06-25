# write a genepop file from a tidy data frame

#' @name write_genepop
#' @title Used internally in stackr to write a genepop file from a tidy data frame
#' @description Write a genepop file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.
#' 
#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info. 
#' @param sep (optional) A character string separating alleles. 
#' Default: \code{sep = NULL}.
#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}.
#' @param genepop.header The first line of the Genepop file.
#' Default: \code{genepop.header = "my firt genepop"}.
#' @param markers.line.format (optional, character) You can write the markers 
#' on a single line separated by commas \code{markers.line.format = "line"}, 
#' or have markers on a separate line, i.e. in one column, 
#' \code{markers.line.format = "column"} 
#' (not very useful with thousands of markers). 
#' Default: \code{markers.line.format = "line"}.
#' @param filename The name of the file written to the working directory.
#' Use the extension ".gen" at the end. 
#' Default: \code{filename = "stackr_genepop.gen"}.
#' @param ... other parameters passed to the function.
#' @return A genepop file is saved to the working directory. 
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
#' @rdname write_genepop
#' @import reshape2
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# to get rid of notes in build check
# if(getRversion() >= "2.15.1") {
#   utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", 
#                            "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", 
#                            "SAMPLES", "ALLELES", 'A1', 'A2', 'COUNT', 
#                            "GT", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", 
#                            "POLYMORPHISM", "POLYMORPHISM_MAX", "other", 
#                            "strata", "hierarchy", "GROUP", ".", 'MARKERS', 
#                            'MARKERS_ALLELES', 'STRATA'
#   )
#   )
# }


write_genepop <- function(
  data,
  sep = NULL, 
  pop.levels = NULL, 
  genepop.header = "my firt genepop", 
  markers.line.format = "line", 
  filename = "stackr_genepop.gen",
  ...
  ) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the genepop file is missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- stackr::read_long_tidy_wide(data = data)
    colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                              pattern = "GENOTYPE", 
                                              replacement = "GT", 
                                              vectorize_all = FALSE)
  } else {
    input <- data
    colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                              pattern = "GENOTYPE", 
                                              replacement = "GT", 
                                              vectorize_all = FALSE)
  }
  
  # Switch colnames LOCUS to MARKERS if found
  if ("LOCUS" %in% colnames(input)) input <- rename(.data = input, MARKERS = LOCUS)
  
  input <- input %>% 
    select(POP_ID, INDIVIDUALS, MARKERS, GT)
  
  if (!is.null(sep)) {
    input <- input %>% 
      mutate(GT = stri_replace_all_fixed(str = GT, pattern = sep, replacement = "", vectorize_all = FALSE))
  }
  
  # pop.levels -----------------------------------------------------------------
  if (!is.null(pop.levels)) {
    input <- input %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =TRUE),
        POP_ID = droplevels(POP_ID)
        ) %>% 
      arrange(POP_ID, INDIVIDUALS, MARKERS)
  }
  
  # Create a marker vector  ------------------------------------------------
  markers <- input %>% select(MARKERS) %>% distinct(MARKERS) %>% arrange(MARKERS)
  markers <- markers$MARKERS
  
  # Wide format ----------------------------------------------------------------
  input <- input %>%
    arrange(MARKERS) %>% 
    group_by(POP_ID, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
    arrange(POP_ID, INDIVIDUALS) %>%
    ungroup() %>% 
    mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))
  
  # Write the file in genepop format -------------------------------------------
  pop <- input$POP_ID # Create a population vector
  input <- split(select(.data = input, -POP_ID), pop) # split genepop by populations
  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = genepop.header, con = filename.connection, sep = "\n") # write the genepop header
  if (markers.line.format == "line") { # write the markers on a single line
    writeLines(text = stri_paste(markers, sep = ",", collapse = ", "), con = filename.connection, sep = "\n") 
  } else {# write the markers on a single column (separate lines)
    writeLines(text = stri_paste(markers, sep = "\n"), con = filename.connection, sep = "\n")
  }
  close(filename.connection) # close the connection
  for (i in 1:length(input)) {
    write_delim(x = as.data.frame("pop"), path = filename, delim = "\n", append = TRUE, col_names = FALSE)
    write_delim(x = input[[i]], path = filename, delim = " ", append = TRUE, col_names = FALSE)
  }
}# End write_genepop
