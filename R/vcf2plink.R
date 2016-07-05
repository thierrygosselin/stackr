# VCF data imputation using Random Forest

#' @name vcf2plink
#' @title VCF to plink with filters and data imputation
#' @description This function can first filter the vcf file 
#' with a whitelist of loci and a blacklist of individuals (optional). 
#' Then it will convert the file
#' to a PLINK tped and tfam file.
#' Map-independent imputation using Random Forest or the most frequent category
#' is also available as an option.
#' @param data The VCF file created by STACKS.
#' @param whitelist.markers (optional) A whitelist containing CHROM (character or integer) and/or LOCUS (integer) and/or 
#' POS (integer) columns header. To filter by CHROM and/or locus and/or by snp.
#' The whitelist is in the directory (e.g. "whitelist.txt"). de novo CHROM column with 'un' need to be changed to 1.
#' @param blacklist.genotype (optional) Useful to erase genotype with below 
#' average quality, e.g. genotype with more than 2 alleles in diploid likely 
#' sequencing errors or genotypes with poor genotype likelihood or coverage. 
#' The blacklist as a minimum of 2 column headers (markers and individuals). 
#' Markers can be 1 column (CHROM or LOCUS or POS), 
#' a combination of 2 (e.g. CHROM and POS or CHROM and LOCUS or LOCUS and POS) or 
#' all 3 (CHROM, LOCUS, POS) The markers columns must be designated: CHROM (character
#' or integer) and/or LOCUS (integer) and/or POS (integer). The id column designated
#' INDIVIDUALS (character) columns header. The blacklist must be in the working 
#' directory (e.g. "blacklist.genotype.txt"). For de novo VCF, CHROM column 
#' with 'un' need to be changed to 1. Default \code{NULL} for no blacklist of 
#' genotypes to erase.
#' @param snp.ld (optional) Minimize linkage disequilibrium (LD) by choosing 
#' among these 3 options: \code{"random"} selection, \code{"first"} or 
#' \code{"last"} SNP on the same read/haplotype. Default = \code{NULL}.
#' @param common.markers (optional) Logical. Default = \code{FALSE}. 
#' With \code{TRUE}, will keep markers present in all the populations.
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param pop.levels (required) A character string with your populations ordered.
#' @param pop.labels (optional) A character string of your populations labels.
#' If you need to rename sampling sites in \code{pop.levels} or combined sites/pop
#' into a different names, here is the place.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param strata (optional) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. Default: \code{strata = NULL}. 
#' This is especially useful if you have a hierarchical or factorial sampling 
#' design. See also \code{\link[adegenet]{genind}} for details.
#' @param pop.select (string) Conduct the conversion on a
#' selected list of populations. Default = \code{NULL} for no selection and keep
#' all population.
#' e.g. \code{pop.select = "QUE"} to select QUE population samples.
#' \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population samples.

#' @param imputation.method (character, optional) 
#' Methods available for map-independent imputations of missing genotype: 
#' (1) \code{"max"} to use the most frequent category for imputations.
#' (2) \code{"rf"} using Random Forest algorithm. 
#' Default: no imputation \code{imputation.method = NULL}.

#' @param impute (character, optional) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
#' Default: \code{"genotype"}.

#' @param imputations.group (character, optional) \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default = \code{"populations"}.

#' @param num.tree (integer, optional) The number of trees to grow in Random Forest. 
#' Default: \code{num.tree = 100}.

#' @param iteration.rf (integer, optional) The number of iterations of missing data algorithm
#' in Random Forest. 
#' Default: \code{iteration.rf = 10}.

#' @param split.number (integer, optional) Non-negative integer value used to specify
#' random splitting in Random Forest. 
#' Default: \code{split.number = 100}.

#' @param verbose (logical, optional) Should trace output be enabled on each iteration
#' in Random Forest ? 
#' Default: \code{verbose = FALSE}.

#' @param parallel.core (optional) The number of core for OpenMP shared-memory parallel
#' programming of Random Forest imputations. For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.


#' @param filename (optional) The prefix for the PLINK file written to the directory.
#' "tped" and "tfam" will be appended by the function. 
#' Default \code{data_date@time.tped, data_date@time.tped}.

#' @details The imputations using Random Forest requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set 
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals 
#' will require 15 min.

#' @return When no imputation is selected the PLINK files are saved to the 
#' working directory. When imputation is selected 4 files are saved to
#' the working directory. The PLINK file returned is in the \code{tped/tfam} 
#' format. The \code{strata} or the population arguments above are used 
#' for first 2 columns of the \code{tfam} file.


#' @export
#' @rdname vcf2plink
#' @import reshape2
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, et al. 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. 
#' American Journal of Human Genetics. 2007; 81: 559–575. doi:10.1086/519795
#' @seealso \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# to get rid of notes in build check
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", 
                           "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", 
                           "SAMPLES", "ALLELES", 'A1', 'A2', 'COUNT', 
                           "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", 
                           "POLYMORPHISM", "POLYMORPHISM_MAX", "other", 
                           "strata", "hierarchy", "GROUP", ".", 'MARKERS', 
                           'MARKERS_ALLELES', 'STRATA'
  )
  )
}


vcf2plink <- function(
  data,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  max.marker = NULL,
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  imputation.method = NULL,
  impute = "genotype",
  imputations.group = "populations",
  num.tree = 100,
  iteration.rf = 10,
  split.number = 100,
  verbose = FALSE,
  parallel.core = detectCores()-1,
  filename = NULL
  
  
) {
  cat("#######################################################################\n")
  cat("######################### stackr::vcf2plink ###########################\n")
  cat("#######################################################################\n")
  
  
  # Checking for missing and/or default arguments ##############################
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  if (is.null(imputation.method)) {
    message("vcf2plink: no imputation...")
  } else {
    message("vcf2plink: with imputations...")
  }
  
  
  # File type detection#########################################################
  if(adegenet::is.genind(data)){
    data.type <- "genind.file"
    message("File type: genind object")
  } else {
    data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
    if (identical(data.type, "##fileformat=VCF") | stri_detect_fixed(str = data, pattern = ".vcf")) {
      data.type <- "vcf.file"
      # message("File type: VCF")
    }
    if (stri_detect_fixed(str = data, pattern = ".tped")) {
      data.type <- "plink.file"
      # message("File type: PLINK")
      if (!file.exists(stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        stop("Missing tfam file with the same prefix as your tped")
      }
    } 
    if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS")) {
      data.type <- "df.file"
      # message("File type: data frame of genotypes")
    }
    if (stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      # data.type <- "haplo.file"
      message("File type: haplotypes from stacks")
      if (is.null(blacklist.genotype)) {
        stop("blacklist.genotype file missing. 
             Use stackr's missing_genotypes function to create this blacklist")
      }
    }
    if (stri_detect_fixed(str = data, pattern = ".gen")) {
      # data.type <- "genepop.file"
      message("File type: genepop")
      message("Multilocus genepop file won't work, only for biallelic markers")
    } 
    
  } # end file type detection
  
  # Strata argument required for VCF and haplotypes files#######################
  if (data.type == "vcf.file" & is.null(strata)) stop("strata argument is required")
  if (data.type == "haplo.file") stop("This function is for biallelic dataset only")
  
  # Import #####################################################################
  input <- stackr::tidy_genomic_data(
    data = data, 
    vcf.metadata = FALSE,
    blacklist.id = blacklist.id, 
    blacklist.genotype = blacklist.genotype, 
    whitelist.markers = whitelist.markers, 
    monomorphic.out = monomorphic.out, 
    max.marker = max.marker,
    snp.ld = snp.ld, 
    common.markers = common.markers, 
    maf.thresholds = maf.thresholds, 
    maf.pop.num.threshold = maf.pop.num.threshold, 
    maf.approach = maf.approach, 
    maf.operator = maf.operator,
    strata = strata, 
    pop.levels = pop.levels, 
    pop.labels = pop.labels, 
    pop.select = pop.select,
    filename = NULL
  )
  
  # create a strata.df
  strata.df <- input %>% 
    distinct(INDIVIDUALS, POP_ID)
  strata <- strata.df
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  # Conversion into PLINK -----------------------------------------------------
  # Tidy VCF
  input <- arrange(.data = input, MARKERS, POP_ID, INDIVIDUALS)
  
  # results no imputation--------------------------------------------------------------------
  write_plink <- function(x, filename) {
    tped <- x %>% 
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
    
    tfam <- x %>%
      distinct(POP_ID, INDIVIDUALS) %>% 
      arrange(INDIVIDUALS) %>% 
      mutate(
        COL3 = rep("0",n()),
        COL4 = rep("0",n()),
        COL5 = rep("0",n()),
        COL6 = rep("-9",n())
      )
    
    # Create a filename to save the output files ********************************
    # Get date and time to have unique filenaming
    if (is.null(filename)) {
      file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
      file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
      file.date <- stri_sub(file.date, from = 1, to = 13)
      filename.tped <- stri_paste("data_", file.date, ".tped")
      filename.tfam <- stri_paste("data_", file.date, ".tfam")
      if (is.null(imputation.method)) {
        filename.tped <- stri_paste("data_imputed_", file.date, ".tped")
        filename.tfam <- stri_paste("data_imputed_", file.date, ".tfam")
      }
    } else {
      filename.tped <- stri_paste(filename, ".tped")
      filename.tfam <- stri_paste(filename, ".tfam")
      if (imputation.method != FALSE) {
        filename.tped <- stri_paste(filename, "_imputed", ".tped")
        filename.tfam <- stri_paste(filename, "_imputed", ".tfam")
      }
    }
    
    write_delim(x = tped, path = filename.tped, col_names = FALSE, delim = " ")
    write_delim(x = tfam, path = filename.tfam, col_names = FALSE, delim = " ")
  } # end write_plink
  
  message("Generating the PLINK tped and tfam files")
  write_plink (x = input, filename = filename)
  
  # Imputations-----------------------------------------------------------------
  if (!is.null(imputation.method)) {
    
    input.imp <- stackr::stackr_imputations_module(
      data = input, 
      imputation.method = imputation.method, 
      impute = impute, 
      imputations.group = imputations.group, 
      num.tree = num.tree, 
      iteration.rf = iteration.rf, 
      split.number = split.number, 
      verbose = verbose, 
      parallel.core = parallel.core, 
      filename = "dataset.imputed.tsv"
    )
    input <- NULL # no longer needed
    
    message("Generating the PLINK tped and tfam files: with imputations")
    write_plink (x = input.imp, filename = filename)
    
  } # End imputations
}
