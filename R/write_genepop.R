# write a genepop file from a tidy data frame

#' @name write_genepop

#' @title Used internally in stackr to write a genepop file from a tidy data frame

#' @description Write a genepop file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info. 

#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}.

#' @param genepop.header The first line of the Genepop file.
#' Default: \code{genepop.header = NULL} will use "stackr genepop with date".

#' @param markers.line (optional, logical) In the genepop and structure
#' file, you can write the markers on a single line separated by 
#' commas \code{markers.line = TRUE}, 
#' or have markers on a separate line, i.e. in one column, for the genepop file
#' (not very useful with thousands of markers) and not printed at all for the
#' structure file.
#' Default: \code{markers.line = TRUE}.

#' @param filename (optional) The file name prefix for the genepop file 
#' written to the working directory. With default: \code{filename = NULL}, 
#' the date and time is appended to \code{stackr_genepop_}.

#' @param ... other parameters passed to the function.

#' @return A genepop file is saved to the working directory. 

#' @details
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
#' If a separator is present, it is automatically removed. 
#' 
#' 
#' 
#' \emph{How to get a tidy data frame ?}
#' \pkg{stackr} \code{\link{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame.
#' @export
#' @rdname write_genepop
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_genepop <- function(
  data,
  pop.levels = NULL, 
  genepop.header = NULL, 
  markers.line = TRUE, 
  filename = NULL,
  ...
) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the genepop file is missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data)
  
  colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                            pattern = "GENOTYPE", 
                                            replacement = "GT", 
                                            vectorize_all = FALSE)
  
  # Switch colnames LOCUS to MARKERS if found
  if ("LOCUS" %in% colnames(input)) input <- rename(.data = input, MARKERS = LOCUS)
  
  input <- input %>% 
    select(POP_ID, INDIVIDUALS, MARKERS, GT)
  
  input <- input %>% 
    mutate(
      GT = stri_replace_all_fixed(
        str = as.character(GT), 
        pattern = c("/", ":", "_", "-", "."), 
        replacement = "", 
        vectorize_all = FALSE),
      GT = stri_pad_left(str = as.character(GT), pad = "0", width = 6)
    )
  
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
  markers <- input %>% distinct(MARKERS) %>% arrange(MARKERS)
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
  
  # Filename ------------------------------------------------------------------
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename <- stri_paste("stackr_genepop_", file.date, ".gen")
  } else {
    filename <- stri_paste(filename, ".gen")
  }
  
  
  # genepop header  ------------------------------------------------------------
  if (is.null(genepop.header)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    genepop.header <- stri_paste("stackr genepop ", file.date)
  }
  
  # genepop construction
  pop <- input$POP_ID # Create a population vector
  input <- split(select(.data = input, -POP_ID), pop) # split genepop by populations
  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = genepop.header, con = filename.connection, sep = "\n") # write the genepop header
  if (markers.line) { # write the markers on a single line
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
