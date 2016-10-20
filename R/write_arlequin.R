# write an Arlequin file from a tidy data frame

#' @name write_arlequin
#' @title Used internally in stackr to write a arlequin file from a tidy data frame
#' @description Write a arlequin file from a tidy data frame
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info.

#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}.

#' @param filename (optional) The file name prefix for the arlequin file 
#' written to the working directory. With default: \code{filename = NULL}, 
#' the date and time is appended to \code{stackr_arlequin_}.

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

#' @return An arlequin file is saved to the working directory. 

#' @references Excoffier, L.G. Laval, and S. Schneider (2005) 
#' Arlequin ver. 3.0: An integrated software package for population genetics 
#' data analysis. Evolutionary Bioinformatics Online 1:47-50.

#' @export
#' @rdname write_arlequin
#' @import reshape2
#' @import dplyr
#' @import stringi

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_arlequin <- function(
  data,
  pop.levels = NULL, 
  filename = NULL,
  ...
) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the arlequin file is missing")
  
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
  
  # pop.levels -----------------------------------------------------------------
  if (!is.null(pop.levels)) {
    input <- input %>%
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE),
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
  npop <- length(unique(input$POP_ID))
  
  # arlequin format ----------------------------------------------------------------
  input <- input %>%
    tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3, extra = "drop", remove = TRUE) %>%
    tidyr::gather(data = ., key = ALLELES, value = GT, -c(POP_ID, INDIVIDUALS, MARKERS)) %>% 
    mutate(
      GT = stri_replace_all_fixed(str = GT, pattern = "000", replacement = "-9", vectorize_all = FALSE)
    ) %>%
    select(INDIVIDUALS, POP_ID, MARKERS, ALLELES, GT) %>% 
    tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
    select(-ALLELES) %>% 
    arrange(POP_ID, INDIVIDUALS)
  
  # Write the file in arlequin format -----------------------------------------
  
  # Filename
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename <- stri_paste("stackr_arlequin_", file.date, ".arp")
  } else {
    filename <- stri_paste(filename, ".arp")
  }
  
  filename.connection <- file(filename, "w") # open the connection to the file
  # Profile section
  write("[Profile]", file = filename.connection) 
  write(paste('Title = "stackr_arlequin_export', " ", file.date, '"'), file = filename.connection, append = TRUE) 
  write(paste("NbSamples = ", npop), file = filename.connection, append = TRUE) # number of pop
  write(paste("GenotypicData = 1"), file = filename.connection, append = TRUE)
  write(paste("LocusSeparator = WHITESPACE"), file = filename.connection, append = TRUE)
  write(paste("GameticPhase = 0"), file = filename.connection, append = TRUE) # 0 = unknown gametic phase
  write(paste("MissingData = '?'"), file = filename.connection, append = TRUE) 
  write(paste("DataType = STANDARD"), file = filename.connection, append = TRUE) 
  write(paste("[Data]"), file = filename.connection, append = TRUE)
  write(paste("[[Samples]]"), file = filename.connection, append = TRUE)
  
  pop <- input$POP_ID # Create a population vector
  input.split <- split(input, pop) # split genepop by populations
  for (i in 1:length(input.split)) {
    # i <- 1
    pop.data <- input.split[[i]]
    pop.name <- unique(pop.data$POP_ID)
    n.ind <- n_distinct(pop.data$INDIVIDUALS)
    write(paste("SampleName = ", pop.name), file = filename.connection, append = TRUE)
    write(paste("SampleSize = ", n.ind), file = filename.connection, append = TRUE)
    write(paste("SampleData = {"), file = filename.connection, append = TRUE)
    pop.data$INDIVIDUALS[seq(from = 2, to = n.ind * 2, by = 2)] <- ""
    # test <- data.frame(input.split[[i]])
    # close(filename.connection) # close the connection
    pop.data <- select(.data = pop.data, -POP_ID)
    ncol.data <- ncol(pop.data)
    pop.data <- as.matrix(pop.data)
    write(x = t(pop.data), ncolumns = ncol.data, sep = "\t", append = TRUE, file = filename.connection)
    write("}", file = filename.connection, append = TRUE)
  }
  
  write(paste("[[Structure]]"), file = filename.connection, append = TRUE)
  write(paste("StructureName = ", "\"", "One cluster", "\""), file = filename.connection, append = TRUE)
  write(paste("NbGroups = 1"), file = filename.connection, append = TRUE)
  write(paste("Group = {"), file = filename.connection, append = TRUE)
  for (i in 1:length(input.split)) {
    pop.data <- input.split[[i]]
    pop.name <- unique(pop.data$POP_ID)
    write(paste("\"", pop.name, "\"", sep = ""), file = filename.connection, append = TRUE)
  }
  write(paste("}"), file = filename.connection, append = TRUE)
} # end write_arlequin
