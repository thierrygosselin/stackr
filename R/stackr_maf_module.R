# MAF module

#' @name stackr_maf_module

#' @title Used internally in stackr to compile and/or filter maf data
#' from a genomic tidy data frame

#' @description compute the minor allele frequency (local and global) using 
#' a genomic tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and might be of interest for users.
#' 
#' @param data A biallelic genomic data set in the working directory or 
#' object in the global environment in wide or long (tidy) formats. 
#' See details for more info. 

#' @inheritParams tidy_genomic_data

#' @return A tidy data frame with 6 columns: 
#' \code{MARKERS, INDIVIDUALS, REF, ALT, GT_VCF, GT_BIN}.
#' \code{GT_VCF}: the genotype in VCF format
#' \code{GT_BIN}: coding used internally to easily convert to genlight, 
#' the coding \code{0, 1, 2, NA} stands for the number of ALT allele in the 
#' genotype and \code{NA} for missing genotype.

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
#' @rdname stackr_maf_module
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs summarise_at
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub
#' @importFrom tibble has_name
#' @importFrom tidyr gather complete separate
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

stackr_maf_module <- function(
  data,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR"
) {
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
  # check genotype column naming
  if (tibble::has_name(input, "GENOTYPE")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input), 
      pattern = "GENOTYPE", 
      replacement = "GT", 
      vectorize_all = FALSE)
  }
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # Minor Allele Frequency filter ----------------------------------------------
  # maf.thresholds <- c(0.05, 0.1) # test
  if (!is.null(maf.thresholds)) { # with MAF
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
    message("MAF filter: yes")
    
    if (tibble::has_name(input, "GT_VCF")) {
      maf.local <- input %>%
        dplyr::filter(GT_VCF != "./.") %>%
        dplyr::group_by(MARKERS, POP_ID, REF, ALT) %>%
        dplyr::summarise(
          N = as.numeric(n()),
          PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
          QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
        ) %>%
        dplyr::mutate(MAF_LOCAL = ((QQ * 2) + PQ) / (2 * N))
      
      maf.global <- maf.local %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::summarise_at(.tbl = ., .cols = c("N", "PQ", "QQ"), dplyr::funs(sum)) %>%
        dplyr::mutate(MAF_GLOBAL = ((QQ * 2) + PQ) / (2 * N)) %>%
        dplyr::select(MARKERS, MAF_GLOBAL)
      
      maf.data <- maf.global %>%
        dplyr::left_join(maf.local, by = c("MARKERS")) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
      
      maf.local <- maf.global <- NULL
    } else {
      # We split the alleles here to prep for MAF
      maf.data <- input %>%
        dplyr::filter(GT != "000000") %>% 
        dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6)
        ) %>% 
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
        dplyr::group_by(MARKERS, GT, POP_ID) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup() %>%
        tidyr::complete(data = ., POP_ID, nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::rename(n.al.pop = n) %>% 
        dplyr::arrange(MARKERS, GT) %>% 
        dplyr::group_by(MARKERS, GT) %>%
        dplyr::mutate(n.al.tot = sum(n.al.pop)) %>% 
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(MAF_GLOBAL = min(n.al.tot)/sum(n.al.pop)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(MAF_LOCAL = n.al.pop/sum(n.al.pop)) %>% 
        dplyr::arrange(MARKERS, POP_ID, GT) %>% 
        dplyr::group_by(MARKERS, POP_ID) %>% 
        dplyr::filter(n.al.pop == min(n.al.pop)) %>% 
        dplyr::distinct(MARKERS, POP_ID, .keep_all = TRUE) %>% 
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
    }# end maf calculations
    
    write_tsv(
      x = maf.data, 
      path = "maf.data.tsv",
      col_names = TRUE, 
      append = FALSE
    )
    message("The MAF table was written in your folder")
    
    if (maf.approach == "haplotype") {
      if (!tibble::has_name(input, "CHROM")) {
        stop("The haplotype approach for maf needs locus and snp info from vcf")
      }
      vcf.maf <- tidyr::separate(
        data = maf.data, 
        col = MARKERS, 
        into = c("CHROM", "LOCUS", "POS"), 
        sep = "__", 
        remove = FALSE, 
        extra = "warn"
      )
      if (maf.operator == "OR") {
        vcf.maf <- vcf.maf %>%
          dplyr::group_by(LOCUS, POP_ID) %>%
          dplyr::summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          dplyr::filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
          dplyr::group_by(LOCUS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n >= maf.pop.num.threshold) %>%
          dplyr::select(LOCUS) %>%
          dplyr::left_join(input, by = "LOCUS") %>%
          dplyr::arrange(LOCUS, POP_ID)
      } else {# AND operator between local and global maf
        vcf.maf <- vcf.maf %>%
          dplyr::group_by(LOCUS, POP_ID) %>%
          dplyr::summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          dplyr::filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
          dplyr::group_by(LOCUS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n >= maf.pop.num.threshold) %>%
          dplyr::select(LOCUS) %>%
          dplyr::left_join(input, by = "LOCUS") %>%
          dplyr::arrange(LOCUS, POP_ID)
      }
      vcf.maf <- dplyr::select(vcf.maf, -c(CHROM, LOCUS, POS))
    } # end maf haplotype approach
    
    if (maf.approach == "SNP") { # SNP approach
      if (maf.operator == "OR") {
        vcf.maf <- maf.data %>%
          dplyr::group_by(MARKERS, POP_ID) %>%
          dplyr::summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          dplyr::filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n >= maf.pop.num.threshold) %>%
          dplyr::select(MARKERS) %>%
          dplyr::left_join(input, by = "MARKERS") %>%
          dplyr::arrange(MARKERS, POP_ID)
      } else {# AND operator between local and global maf
        vcf.maf <- maf.data %>%
          dplyr::group_by(MARKERS, POP_ID) %>%
          dplyr::summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          dplyr::filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
          dplyr::group_by(MARKERS) %>%
          dplyr::tally(.) %>%
          dplyr::filter(n >= maf.pop.num.threshold) %>%
          dplyr::select(MARKERS) %>%
          dplyr::left_join(input, by = "MARKERS") %>%
          dplyr::arrange(MARKERS, POP_ID)
      }
    } # end maf snp approach
    
    
    message(
      stringi::stri_join(
        "The number of MARKERS removed by the MAF filters = ", 
        dplyr::n_distinct(input$MARKERS) - dplyr::n_distinct(vcf.maf$MARKERS), "\n", 
        "The number of MARKERS before -> after the MAF filters: ", 
        dplyr::n_distinct(input$MARKERS)," -> ", dplyr::n_distinct(vcf.maf$MARKERS), 
        " MARKERS"
      )
    )
    res <- list(input = vcf.maf, maf.data = maf.data)
    return(res)
  } # End of MAF filters
}# End function
