# write a SNPRelate object from a tidy data frame

#' @name write_snprelate
#' @title Write a SNPRelate object from a tidy data frame
#' @description Write a \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}
#' object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{stackr} \code{\link{tidy_genomic_data}}.
#' \strong{The genotypes are biallelic.}


# @param filename (optional) The file name of Genomic Data Structure (GDS file).
# The file must finish with \code{.gds}. If not provided a filename is provided
# with current date and time. If filename chosen is already present in the
# working directory, the default is chosen.
# Default: \code{filename = NULL}.

#' @export
#' @rdname write_snprelate

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom tibble has_name
#' @importFrom tidyr spread
#' @importFrom SNPRelate snpgdsOpen snpgdsClose snpgdsCreateGeno

#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_snprelate <- function(data) {
  
  # Check that snprelate is installed
  if (!"SNPRelate" %in% utils::installed.packages()[,"Package"]) {
    stop("Please install SNPRelate:\n
         github::zhengxwen/SNPRelate")
  }
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  
  
  # if (!is.null(filename)) {
  #   # filename <- "stackr.gds"# test
  #   scan.gds.file <- list.files(path = getwd(), pattern = filename, full.names = FALSE)[1]
  #   if (!is.na(scan.gds.file)) filename <- NULL
  # }
  
  # if (is.null(filename)) {
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stringi::stri_replace_all_fixed(
      str = file.date,
      pattern = c("-", " ", ":"),
      replacement = c("", "_", ""), 
      vectorize_all = FALSE
    )
    filename <- stringi::stri_join("stackr_", file.date, ".gds", sep = "")
  # }

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- stackr::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }
  # check genotype column naming
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input), 
    pattern = "GENOTYPE", 
    replacement = "GT", 
    vectorize_all = FALSE
  )
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # SNPRelate prep -------------------------------------------------------------
  biallelic <- stackr::detect_biallelic_markers(data = input)
  if (!biallelic) stop("SNPRelate requires biallelic genotypes")
  verbose <- FALSE
  
  if (verbose) message("Generating GDS format...")
  # keep markers in common
  # gds.genotypes <- suppressMessages(stackr::keep_common_markers(data = input))
  
  strata.df <- dplyr::distinct(input, POP_ID, INDIVIDUALS) %>%
    dplyr::mutate(POP_ID = factor(POP_ID))
  
  snp.id <- dplyr::distinct(.data = input, MARKERS) %>%
    dplyr::arrange(MARKERS) %>%
    purrr::flatten_chr(.)
  
  input <- suppressWarnings(
    dplyr::select(.data = input, MARKERS, INDIVIDUALS, GT_BIN) %>%
      dplyr::group_by(MARKERS) %>% 
      tidyr::spread(data = ., key = INDIVIDUALS, value = GT_BIN) %>% 
      dplyr::arrange(MARKERS) %>%
      tibble::column_to_rownames(df = ., var = "MARKERS") %>% 
      data.matrix(.)
  )
  
  # Generate GDS format --------------------------------------------------------
  
  SNPRelate::snpgdsCreateGeno(
    gds.fn = filename,
    genmat = input,
    sample.id = strata.df$INDIVIDUALS,
    snp.id = snp.id,
    snp.rs.id = NULL,
    snp.chromosome = NULL,
    snp.position = NULL,
    snp.allele = NULL,
    snpfirstdim = TRUE,
    compress.annotation = "ZIP_RA.max",
    compress.geno = "",
    other.vars = NULL
  )
  
  gds.file.connection <- SNPRelate::snpgdsOpen(filename)

  message("\nNote: \n
GDS filename: ", filename)

  message("\nTo close the connection use SNPRelate::snpgdsClose(OBJECT_NAME)")
  return(gds.file.connection)
} # End write_snprelate
