# Detect detect_genomic_format

#' @name detect_genomic_format
#' @title Used internally in stackr to detect the file format
#' @description Detect file format of genomic data set.
#' @param data A genomic data set in the global environment
#' @rdname detect_genomic_format
#' @importFrom stringi stri_detect_fixed stri_replace_all_fixed
#' @importFrom adegenet is.genind
#' @importFrom tibble has_name
#' @import strataG
# @keywords internal
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_genomic_format <- function(data){
  
  if (!is.vector(data)) {
    if (tibble::has_name(data, "POP_ID") & tibble::has_name(data, "INDIVIDUALS") & tibble::has_name(data, "MARKERS")) {
      data.type <- "df.file"
    } else {
      if (adegenet::is.genind(data)) {
        data.type <- "genind.file"
        # message("File type: genind object")
      } else if (strataG::is.gtypes(data)) {
        data.type <- "gtypes"
        } else {
        stop("Input file not recognised")
      }
    }
  } else {
    data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
    
    if (identical(data.type, "##fileformat=VCF") | stringi::stri_detect_fixed(str = data, pattern = ".vcf")) {
      data.type <- "vcf.file"
      # message("File type: VCF")
    }
    
    if (stringi::stri_detect_fixed(str = data, pattern = ".tped")) {
      data.type <- "plink.file"
      # message("File type: PLINK")
      if (!file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        stop("Missing tfam file with the same prefix as your tped")
      }
    }
    
    if (stringi::stri_detect_fixed(str = data.type, pattern = "POP_ID") | stringi::stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stringi::stri_detect_fixed(str = data.type, pattern = "MARKERS") | stringi::stri_detect_fixed(str = data.type, pattern = "LOCUS")) {
      data.type <- "df.file"
      # message("File type: data frame of genotypes")
    }
    if (stringi::stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      data.type <- "haplo.file"
    }
    if (stringi::stri_detect_fixed(str = data, pattern = ".gen")) {
      data.type <- "genepop.file"
    }
  } # end file type detection
  return(data.type)
} # End detect_genomic_format
