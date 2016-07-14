# write a vcf file from a tidy data frame

#' @name write_vcf
#' @title Used internally in stackr to write a vcf file from a tidy 
#' data frame
#' @description Write a vcf file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.
#' 
#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. See details for more info. 

#' @param pop.info (optional, logical) Should the population information be 
#' included in the FORMAT field (along the GT info for each samples ?). To make
#' the VCF population-ready use \code{pop.info = TRUE}. The populatio information
#' must be included in the \code{POP_ID} column of the tidy dataset.
#' Default: \code{pop.info = FALSE}.

#' @param filename (optional) The file name prefix for the vcf file 
#' written to the working directory. With default: \code{filename = NULL}, 
#' the date and time is appended to \code{stackr_vcf_file_}.

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
#' @rdname write_vcf
#' @import reshape2
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_vcf <- function(data, pop.info = FALSE, filename = NULL) {
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
  colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                            pattern = "GENOTYPE", 
                                            replacement = "GT", 
                                            vectorize_all = FALSE)
  
  
  # REF/ALT Alleles and VCF genotype format ------------------------------------
  if (!tibble::has_name(input, "GT_VCF")) {
    ref.alt.alleles.change <- ref_alt_alleles(data = input)
    input <- left_join(input, ref.alt.alleles.change, by = c("MARKERS", "INDIVIDUALS"))
  }
  
  # Include CHROM, LOCUS, POS --------------------------------------------------
  if (!tibble::has_name(input, "CHROM")) {
    input <- mutate(
      .data = input, 
      CHROM = rep("1", n()),
      LOCUS = MARKERS,
      POS = MARKERS
    )
  }
  
  # Remove the POP_ID column ---------------------------------------------------
  if (tibble::has_name(input, "POP_ID") & (!pop.info)) {
    input <- select(.data = input, -POP_ID)
  }
  
  # Info field -----------------------------------------------------------------
  info.field <- suppressWarnings(
    input %>% 
      group_by(MARKERS) %>%
      filter(GT_VCF != "./.") %>% 
      tally %>% 
      mutate(INFO = stri_paste("NS=", n, sep = "")) %>% 
      select(-n)
  )
  
  # VCF body  ------------------------------------------------------------------
  GT_VCF_POP_ID <- NULL
  if (pop.info) {
    output <- suppressWarnings(
      left_join(input, info.field, by = "MARKERS") %>% 
        select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INFO, INDIVIDUALS, GT_VCF, POP_ID) %>%
        mutate(GT_VCF_POP_ID = stri_paste(GT_VCF, POP_ID, sep = ":")) %>%
        select(-c(GT_VCF, POP_ID)) %>% 
        group_by(MARKERS, CHROM, LOCUS, POS, INFO, REF, ALT) %>% 
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT_VCF_POP_ID) %>%
        ungroup() %>% 
        mutate(
          QUAL = rep(".", n()),
          FILTER = rep("PASS", n()),
          FORMAT = rep("GT:POP", n())
        )
    )
    
  } else {
    output <- suppressWarnings(
      left_join(input, info.field, by = "MARKERS") %>% 
        select(MARKERS, CHROM, LOCUS, POS, REF, ALT, INDIVIDUALS, GT_VCF, INFO) %>% 
        group_by(MARKERS, CHROM, LOCUS, POS, INFO, REF, ALT) %>% 
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT_VCF) %>%
        ungroup() %>% 
        mutate(
          QUAL = rep(".", n()),
          FILTER = rep("PASS", n()),
          FORMAT = rep("GT", n())
        )
    )
  }
  
  output <- output %>% 
    arrange(CHROM, LOCUS, POS) %>%
    ungroup() %>%
    select(-MARKERS) %>%
    select('#CHROM' = CHROM, POS, ID = LOCUS, REF, ALT, QUAL, FILTER, INFO, FORMAT, everything())
  
  
  # Filename ------------------------------------------------------------------
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename <- stri_paste("stackr_vcf_file_", file.date, ".vcf")
  } else {
    filename <- stri_paste(filename, ".vcf")
  }
  
  # File format ----------------------------------------------------------------
  write_delim(x = data_frame("##fileformat=VCFv4.2"), path = filename, delim = " ", append = FALSE, col_names = FALSE)
  
  # File date ------------------------------------------------------------------
  file.date <- stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")
  file.date <- stri_paste("##fileDate=", file.date, sep = "")
  write_delim(x = data_frame(file.date), path = filename, delim = " ", append = TRUE, col_names = FALSE)
  
  # Source ---------------------------------------------------------------------
  write_delim(x = data_frame(stri_paste("##source=stackr_v.", utils::packageVersion("stackr"))), path = filename, delim = " ", append = TRUE, col_names = FALSE)

  # Info field 1 ---------------------------------------------------------------
  info1 <- as.data.frame('##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">')
  utils::write.table(x = info1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)


  # Format field 1 -------------------------------------------------------------
  format1 <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  format1 <- as.data.frame(format1)
  utils::write.table(x = format1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Format field 2 ---------------------------------------------------------------
  if (pop.info) {
    format2 <- as.data.frame('##FORMAT=<ID=POP_ID,Number=1,Type=Character,Description="Population identification of Sample">')
    utils::write.table(x = format2, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  # Write the prunned vcf to the file ------------------------------------------
  suppressWarnings(write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE))
} # end write_vcf
