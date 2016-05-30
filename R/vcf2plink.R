# VCF data imputation using Random Forest

#' @name vcf2plink
#' @title VCF to genepop with filters and data imputation
#' @description This function can first filter the vcf file 
#' with a whitelist of loci and a blacklist of individuals (optional). 
#' Then it will convert the file
#' to a PLINK tped and tmap file.
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

#' @param genotype.format The genotype are in \code{character} (e.g. \code{A/T}), 
#' or \code{integer} (e.g. \code{01/04}) format. 
#' Where \code{A = 01, C = 02, G = 03 and T = 04}. 
#' Default: \code{genotype.format = "character"}.

#' @param imputation.method Should a map-independent imputations of markers be
#' computed. Available choices are: (1) \code{FALSE} for no imputation.
#' (2) \code{"max"} to use the most frequent category for imputations.
#' (3) \code{"rf"} using Random Forest algorithm. Default = \code{FALSE}.
#' @param impute (character) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
#' @param imputations.group \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default = \code{"populations"}.
#' @param num.tree The number of trees to grow in Random Forest. Default is 100.
#' @param iteration.rf The number of iterations of missing data algorithm 
#' in Random Forest. Default is 10.
#' @param split.number Non-negative integer value used to specify 
#' random splitting in Random Forest. Default is 100.
#' @param verbose Logical. Should trace output be enabled on each iteration 
#' in Random Forest ? Default is \code{FALSE}.
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
#' the working directory. The PLINK file returned is in the \code{tped/tmap} 
#' format. The \code{strata} or the population arguments above are used 
#' for first 2 columns of the \code{tmap} file.


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


vcf2plink <- function(data, 
                      whitelist.markers,
                      blacklist.genotype,
                      snp.ld,
                      common.markers,
                      blacklist.id,
                      pop.id.start,
                      pop.id.end,
                      pop.levels,
                      pop.labels,
                      strata,
                      pop.select,
                      genotype.format,
                      imputation.method,
                      impute,
                      imputations.group,
                      num.tree,
                      iteration.rf,
                      split.number,
                      verbose,
                      parallel.core,
                      filename
) {
  if (missing(data)) stop("Input file missing")
  if (missing(whitelist.markers)) whitelist.markers <- NULL # no Whitelist
  if (missing(blacklist.id)) blacklist.id <- NULL # No blacklist of ID
  if (missing(blacklist.genotype)) blacklist.genotype <- NULL # no genotype to erase
  if (missing(pop.id.start)) pop.id.start <- NULL
  if (missing(pop.id.end)) pop.id.end <- NULL
  if (missing(strata)) strata <- NULL
  if (missing(pop.levels)) stop("pop.levels required")
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels
  if (missing(pop.select)) pop.select <- NULL
  if (missing(snp.ld)) snp.ld <- NULL
  if (missing(common.markers)) common.markers <- FALSE
  if (missing(genotype.format)) genotype.format <- "character"
  if (missing(imputation.method)) imputation.method <- FALSE
  if (imputation.method != FALSE & missing(impute)) stop("impute argument is necessary")
  if (imputation.method == FALSE & missing(impute)) impute <- NULL
  if (missing(imputations.group)) imputations.group <- "populations"
  if (missing(num.tree)) num.tree <- 100
  if (missing(iteration.rf)) iteration.rf <- 10
  if (missing(split.number)) split.number <- 100
  if (missing(verbose)) verbose <- FALSE
  if (missing(parallel.core)) parallel.core <- detectCores()-1
  if (missing(filename)) filename <- NULL
  
  if (imputation.method == "FALSE") {
    message("vcf2plink: without imputation...")
  } else {
    message("vcf2plink: with/without imputations...")
  }
  
  # Create a filename to save the output files ********************************
  # Get date and time to have unique filenaming
  if (is.null(filename)) {
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename.tped <- stri_paste("data_", file.date, ".tped")
    filename.tmap <- stri_paste("data_", file.date, ".tmap")
    if (imputation.method != FALSE) {
      filename.tped.imp <- stri_paste("data_imputed_", file.date, ".tped")
      filename.tmap.imp <- stri_paste("data_imputed_", file.date, ".tmap")
    }
  } else {
    filename.tped <- stri_paste(filename, ".tped")
    filename.tmap <- stri_paste(filename, ".tmap")
    if (imputation.method != FALSE) {
      filename.tped.imp <- stri_paste(filename, "_imputed", ".tped")
      filename.tmap.imp <- stri_paste(filename, "_imputed", ".tmap")
    }
  }
  
  
  # Import whitelist of markers ************************************************
  if (is.null(whitelist.markers)) { # no Whitelist
    message("Whitelist of markers: no")
  } else { # with Whitelist of markers
    message("Whitelist of markers: yes")
    whitelist.markers <- read_tsv(whitelist.markers, col_names = TRUE)
    columns.names.whitelist <- colnames(whitelist.markers)
    if ("CHROM" %in% columns.names.whitelist) {
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
  }
  
  # Import blacklist id ********************************************************
  if (is.null(blacklist.id)) { # No blacklist of ID
    message("Blacklisted individuals: no")
  } else { # With blacklist of ID
    message("Blacklisted individuals: yes")
    blacklist.id <- read_tsv(blacklist.id, col_names = TRUE)
  }
  # Import/read VCF ------------------------------------------------------------- 
  message("Importing the VCF...")
  
  if (is.null(strata) & is.null(pop.id.start) & is.null(pop.id.end)) {
    stop("pop.id.start and pop.id.end or strata arguments are required to 
         identify your populations with a VCF file")
  }
  
  input <- data.table::fread(
    input = data,
    sep = "\t",
    stringsAsFactors = FALSE, 
    header = TRUE,
    # Automatically filter with blacklist of id
    drop = c("QUAL", "FILTER", "INFO", blacklist.id$INDIVIDUALS),
    skip = "CHROM",
    showProgress = TRUE,
    verbose = FALSE
  ) %>% 
    as_data_frame() %>% 
    rename(LOCUS = ID, CHROM = `#CHROM`) %>%
    mutate(
      CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
    )
  
  # Filter with whitelist of markers
  if (!is.null(whitelist.markers)) {
    input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
  }
  
  # Detect STACKS version
  if (stri_detect_fixed(input$FORMAT[1], "AD")) {
    stacks.version <- "new"
  } else {
    stacks.version <- "old"
  }
  input <- input %>% select(-FORMAT)
  
  # Tidying the VCF to make it easy to work on the data for conversion
  message("Making the VCF population wise")
  input <- input %>%
    tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, LOCUS, POS, REF, ALT)) # Gather individuals in 1 colummn
  
  if (is.null(strata)){
    input <- input %>%
      mutate( # Make population ready
        POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered = TRUE),
        INDIVIDUALS =  as.character(INDIVIDUALS)
      )
  } else { # Make population ready with the strata provided
    strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
      rename(POP_ID = STRATA)
    
    input <- input %>%
      mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
      left_join(strata.df, by = "INDIVIDUALS") %>% 
      mutate(POP_ID = factor(POP_ID, levels = unique(pop.labels), ordered =TRUE))
  }
  
  # Pop select
  if (!is.null(pop.select)) {
    message(stri_join(length(pop.select), "population(s) selected", sep = " "))
    input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
  }
  
  # Conversion into PLINK -----------------------------------------------------
  # Tidy VCF
  message("Tidy vcf into PLINK factory...")
  
  if (stacks.version == "new") { # with new version of stacks > v.1.29
    input <- input %>%
      tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"),
                      sep = ":", extra = "warn") %>%
      select(-c(READ_DEPTH, ALLELE_DEPTH, GL))
  } else { # stacks version prior to v.1.29 had no Allele Depth field...
    input <- input %>%
      tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "GL"),
                      sep = ":", extra = "warn") %>%
      select(-c(READ_DEPTH, GL))
  }
  
  # Blacklist genotypes ********************************************************
  if (is.null(blacklist.genotype)) { # no Whitelist
    message("Erasing genotype: no")
  } else {
    message("Erasing genotype: yes")
    blacklist.genotype <- read_tsv(blacklist.genotype, col_names = TRUE)
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    if ("CHROM" %in% columns.names.blacklist.genotype) {
      columns.names.blacklist.genotype$CHROM <- as.character(columns.names.blacklist.genotype$CHROM)
    }
    
    # control check to keep only individuals in pop.select
    if (!is.null(pop.select)) {
      message("Control check to keep only individuals present in pop.select")
      # updating the blacklist.genotype
      if (is.null(strata)){
        blacklist.genotype <- suppressWarnings(
          blacklist.genotype  %>% 
            mutate( # Make population ready
              POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
              POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = F), levels = unique(pop.labels), ordered =T),
              INDIVIDUALS =  as.character(INDIVIDUALS) 
            ) %>% 
            filter(POP_ID %in% pop.select) %>% 
            select(-POP_ID)
        )
      } else {
        blacklist.genotype <- suppressWarnings(
          blacklist.genotype %>%
            mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
            left_join(strata.df, by = "INDIVIDUALS") %>% 
            filter(POP_ID %in% pop.select) %>% 
            select(-POP_ID)
        )
      }
    }
    
    # control check to keep only whitelisted markers from the blacklist of genotypes
    if (!is.null(whitelist.markers)) {
      blacklist.genotype <- blacklist.genotype
      message("Control check to keep only whitelisted markers present in the blacklist of genotypes to erase.")
      # updating the whitelist of markers to have all columns that id markers
      whitelist.markers.ind <- input %>% select(CHROM, LOCUS, POS, INDIVIDUALS) %>% distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      
      
      # updating the blacklist.genotype
      blacklist.genotype <- suppressWarnings(semi_join(whitelist.markers.ind, blacklist.genotype, by = columns.names.blacklist.genotype))
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # control check to remove blacklisted individuals from the blacklist of genotypes
    if (!is.null(blacklist.id)) {
      message("Control check to remove blacklisted individuals present in the blacklist of genotypes to erase.")
      blacklist.genotype <- suppressWarnings(anti_join(blacklist.genotype, blacklist.id, by = "INDIVIDUALS"))
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # Add one column that will allow to include the blacklist in the dataset 
    # by x column(s) of markers
    blacklist.genotype <- mutate(.data = blacklist.genotype, ERASE = rep("erase", n()))
    
    input <- suppressWarnings(
      input %>%
        full_join(blacklist.genotype, by = columns.names.blacklist.genotype) %>%
        mutate(
          ERASE = stri_replace_na(str = ERASE, replacement = "ok"),
          GT = ifelse(ERASE == "erase", "./.", GT)
        ) %>% 
        select(-ERASE)
    )
  }# End erase genotypes
  
  
  # dump unused object
  blacklist.genotype <- NULL
  blacklist.id <- NULL
  whitelist.markers <- NULL
  
  # LD control... keep only 1 SNP per haplotypes/reads (optional) ************
  if (!is.null(snp.ld)) {
    message("Minimizing LD...")
    snp.locus <- input %>% select(LOCUS, POS) %>% distinct(POS)
    # Random selection
    if (snp.ld == "random") {
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>%
        sample_n(size = 1, replace = FALSE)
      message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP randomly selected to keep 1 SNP per read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # Fist SNP on the read
    if (snp.ld == "first") {
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>%
        summarise(POS = min(POS))
      message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # Last SNP on the read
    if (snp.ld == "last") {
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>%
        summarise(POS = max(POS))
      message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # filtering the VCF to minimize LD
    input <- input %>% semi_join(snp.select, by = c("LOCUS", "POS"))
    message("Filtering the tidy VCF to minimize LD by keeping only 1 SNP per short read/haplotype")
  } # End of snp.ld control
  
  # Unique markers id: combine CHROM, LOCUS and POS into MARKERS *************
  input <- input %>%
    mutate(
      POS = stri_pad_left(str = POS, width = 8, pad = "0"),
      LOCUS = stri_pad_left(str = LOCUS, width = 8, pad = "0")
    ) %>%
    arrange(CHROM, LOCUS, POS) %>%
    tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
  
  # Markers in common between all populations (optional) *********************
  if (common.markers == TRUE) { # keep only markers present in all pop
    message("Using markers common in all populations")
    pop.number <- n_distinct(input$POP_ID)
    
    pop.filter <- input %>%
      filter(GT != "./.") %>%
      group_by(MARKERS) %>%
      filter(n_distinct(POP_ID) == pop.number) %>%
      arrange(MARKERS) %>%
      select(MARKERS) %>%
      distinct(MARKERS)
    
    
    message(stri_join("Number of original markers = ", n_distinct(input$MARKERS), 
                      "\n", "Number of markers present in all the populations = ", 
                      n_distinct(pop.filter$MARKERS), "\n", 
                      "Number of markers removed = ", 
                      n_distinct(input$MARKERS) - n_distinct(pop.filter$MARKERS))
    )
    input <- suppressWarnings(input %>% semi_join(pop.filter, by = "MARKERS"))
    pop.filter <- NULL # ununsed object
  } # End common markers
  
  # Change the genotype coding  **********************************************
  message("Recoding genotypes for PLINK")
  
  if (genotype.format == "integer") {
    input <- input %>%
      mutate(
        REF= stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE), # replace nucleotide with numbers
        ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE),# replace nucleotide with numbers
        REF = stri_pad_left(str = REF, width = 2, pad = "0"),
        ALT = stri_pad_left(str = ALT, width = 2, pad = "0"),
        GT = ifelse(GT == "0/0", stri_join(REF, REF, sep = "_"),
                    ifelse(GT == "1/1",  stri_join(ALT, ALT, sep = "_"),
                           ifelse(GT == "0/1", stri_join(REF, ALT, sep = "_"),
                                  ifelse(GT == "1/0", stri_join(ALT, REF, sep = "_"), "00_00")
                           )
                    )
        )
      ) %>%
      arrange(MARKERS, POP_ID) %>%
      select(-c(REF, ALT))
  } else {
    input <- input %>%
      mutate(
        GT = ifelse(GT == "0/0", stri_join(REF, REF, sep = "_"),
                    ifelse(GT == "1/1",  stri_join(ALT, ALT, sep = "_"),
                           ifelse(GT == "0/1", stri_join(REF, ALT, sep = "_"),
                                  ifelse(GT == "1/0", stri_join(ALT, REF, sep = "_"), "00_00")
                           )
                    )
        )
      ) %>%
      arrange(MARKERS, POP_ID) %>%
      select(-c(REF, ALT))
  }
  
  # results no imputation--------------------------------------------------------------------
  message("Output: No imputation")
  
  # to create a PLINK tped and tmap
  tped <- input %>% 
    mutate(GT = ifelse(GT == "00_00", "0_0", GT)) %>% 
    arrange(INDIVIDUALS) %>% 
    mutate(
      COL1 = rep("0", n()),
      COL3 = rep("0", n()),
      COL4 = rep("0", n())
    ) %>% 
    select(COL1, MARKERS, COL3, COL4, INDIVIDUALS, GT) %>% 
    tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>% 
    tidyr::gather(ALLELES, GT, -c(COL1, MARKERS, COL3, COL4, INDIVIDUALS)) %>% 
    arrange(INDIVIDUALS, ALLELES) %>% 
    tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_") %>%
    group_by(COL1, MARKERS, COL3, COL4) %>% 
    tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GT) %>% 
    arrange(MARKERS)
  
  write_delim(x = tped, path = filename.tped, col_names = FALSE, delim = " ")
  
  # Create a tmap file
  tfam <- input %>% 
    select(POP_ID, INDIVIDUALS) %>% 
    distinct(POP_ID, INDIVIDUALS) %>% 
    arrange(POP_ID, INDIVIDUALS) %>% 
    mutate(
      COL3 = rep("0",n()),
      COL4 = rep("0",n()),
      COL5 = rep("0",n()),
      COL6 = rep("-9",n())
    )
  
  # the tfam must have the same name as the tped 
  write_delim(x = tfam, path = filename.tmap, col_names = FALSE, delim = " ")

  # Imputations ***************************************************************
  if (imputation.method != "FALSE") {
    
    if (impute == "genotype") {
      input.prep <- input %>%
        mutate(
          GT = stri_replace_all_fixed(GT, pattern = "00_00", replacement = "NA", vectorize_all = FALSE),
          GT = replace(GT, which(GT == "NA"), NA)
        ) %>%
        group_by(INDIVIDUALS, POP_ID) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>%
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
    }
    
    if (impute == "allele") {
      input.prep <- input %>%
        tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
        tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
        mutate(
          GT = stri_replace_all_fixed(GT, pattern = "00", replacement = "NA", vectorize_all = FALSE),
          GT = replace(GT, which(GT == "NA"), NA)
        ) %>% 
        # tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ":") %>%
        group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
    }
    input <- NULL # unused object
    
    # Imputation with Random Forest
    if (imputation.method == "rf") {
      # Parallel computations options
      options(rf.cores = parallel.core, mc.cores = parallel.core)
      
      # imputations using Random Forest with the package randomForestSRC
      impute_genotype_rf <- function(x) {
        randomForestSRC::impute.rfsrc(data = x,
                                      ntree = num.tree,
                                      nodesize = 1,
                                      nsplit = split.number,
                                      nimpute = iteration.rf,
                                      do.trace = verbose)
      } # End of imputation function
      
      # Random Forest by pop
      if (imputations.group == "populations") {
        message("Imputations computed by populations, take a break...")
        df.split.pop <- split(x = input.prep, f = input.prep$POP_ID) # slip data frame by population
        pop.list <- names(df.split.pop) # list the pop
        imputed.dataset <-list() # create empty list
        
        # Function to go through the populations
        impute_rf_pop <- function(pop.list, ...){
          sep.pop <- df.split.pop[[pop.list]]
          sep.pop <- suppressWarnings(
            plyr::colwise(factor, exclude = NA)(sep.pop)
          )
          # message of progress for imputations by population
          message(paste("Completed imputations for pop ", pop.list, sep = ""))
          imputed.dataset <- impute_genotype_rf(sep.pop)
          return(imputed.dataset)
        } # End impute_rf_pop
        
        input.imp <- list()
        input.imp <- parallel::mclapply(
          X = pop.list, 
          FUN = impute_rf_pop, 
          mc.preschedule = FALSE, 
          mc.silent = FALSE, 
          mc.cores = parallel.core
        )
        
        # Compiling the results
        message("Compiling imputations results")
        input.imp <- suppressWarnings(bind_rows(input.imp))
        
        # Second round of imputations (globally) to remove introduced NA 
        # In case that some pop don't have the markers
        input.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(input.imp)) # Make the columns factor
        input.imp <- impute_genotype_rf(input.imp) # impute globally
        
        # dump unused objects
        df.split.pop <- NULL
        pop.list <- NULL
        sep.pop <- NULL
        imputed.dataset <- NULL
        input.prep <- NULL
        
      } # End imputation RF populations
      # Random Forest global
      if (imputations.group == "global") { # Globally (not by pop_id)
        message("Imputations computed globally, take a break...")
        input.prep <- plyr::colwise(factor, exclude = NA)(input.prep)
        input.imp <- impute_genotype_rf(input.prep)
        
        input.prep <- NULL # remove unused object
      } # End imputation RF global
    } # End imputation RF
    
    # Imputation using the most common genotype
    if (imputation.method == "max") { # End imputation max
      if (imputations.group == "populations") {
        message("Imputations computed by populations")
        
        if (impute == "genotype"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
              group_by(MARKERS, POP_ID) %>%
              mutate(
                GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                GT = replace(GT, which(GT == "NA"), NA)
              ) %>%
              # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
              # will take the global observed values by markers for those cases.
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              group_by(INDIVIDUALS, POP_ID) %>% 
              tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
        }
        
        if (impute == "allele"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
              group_by(MARKERS, POP_ID) %>%
              mutate(
                GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                GT = replace(GT, which(GT == "NA"), NA)
              ) %>%
              # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
              # will take the global observed values by markers for those cases.
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
              tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
        }
        input.prep <- NULL # remove unused object
        
      } # End imputation max populations 
      if (imputations.group == "global") {
        # Globally (not by pop_id)
        message("Imputations computed globally")
        
        if (impute == "genotype"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              group_by(INDIVIDUALS, POP_ID) %>% 
              tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
        }
        
        if (impute == "allele"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
              tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
        }
        input.prep <- NULL # remove unused object
      } # End imputation max global 
    } # End imputations max
    
    # transform the imputed dataset into genind object ------------------------
    message("Imputed vcf into PLINK factory...")
    if (impute == "genotype") {
      plink.prep.imp <- suppressWarnings(
        input.imp %>%
          tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% # make tidy
          tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
          tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
          tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_", remove = FALSE) %>%
          arrange(MARKERS, POP_ID, INDIVIDUALS_ALLELES)
      )
    }
    
    if (impute == "allele") {
      plink.prep.imp <- suppressWarnings(
        input.imp %>%
          tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
          tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_", remove = FALSE) %>%
          arrange(MARKERS, POP_ID, INDIVIDUALS_ALLELES)
      ) 
    }
    
    # to create a PLINK tped and tfam
    tped.imp <- plink.prep.imp %>% 
      select(-INDIVIDUALS, -ALLELES) %>% 
      mutate(
        COL1 = rep("0", n()),
        COL3 = rep("0", n()),
        COL4 = rep("0", n())
      ) %>% 
      arrange(POP_ID, INDIVIDUALS_ALLELES) %>% 
      select(COL1, MARKERS, COL3, COL4, INDIVIDUALS_ALLELES, GT) %>%
      group_by(COL1, MARKERS, COL3, COL4) %>% 
      tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GT) %>% 
      arrange(MARKERS)
    
    write_delim(x = tped.imp, path = filename.tped.imp, col_names = FALSE, delim = " ")
    
    tfam.imp <- plink.prep.imp %>% 
      select(POP_ID, INDIVIDUALS) %>% 
      distinct(POP_ID, INDIVIDUALS) %>% 
      arrange(POP_ID, INDIVIDUALS) %>% 
      mutate(
        COL3 = rep("0",n()),
        COL4 = rep("0",n()),
        COL5 = rep("0",n()),
        COL6 = rep("-9",n())
      )
    
    # the tfam must have the same name as the tped 
    write_delim(x = tfam.imp, path = filename.tmap.imp, col_names = FALSE, delim = " ")
    
    # removed objects
    input <- NULL
    input.prep <- NULL
    input.imp <- NULL # remove unused object
    plink.prep <- NULL
    plink.prep.imp <- NULL
  } # End imputation
}
