# Import and prune data file

#' @name data_pruning
#' @title Data set pruning
#' @description This function enable to prune the data set with different 
#' argument (see below) in order to be prep for subsequent analysis.
#' Various input files is offered. Individuals, populations and
#' markers can be filtered and/or selected in several ways using blacklist,
#' whitelist and other arguments.

#' @param data Options include the VCF (1) or an haplotype files (2) created in STACKS 
#' (\code{data = "batch_1.vcf"} and \code{data = "batch_1.haplotypes.tsv"}, 
#' respectively) or a data frame (3) with tab separating the 
#' genotypes in columns (\code{data = "data.tsv"}). 
#' The 1st column is the \code{POP_ID}, 2nd colum 
#' the \code{INDIVIDUALS} and the remaining columns are the markers IDs
#' containing genotypes in the format: 3 digits per allele
#' and no space between alleles (e.g. 235240 : allele1 = 235 and allele2 = 240).
#' Missing genotypes are coded \code{0} or \code{000000}. 
#' Note that the \code{POP_ID} column can be any hierarchical grouping. 
#' See the argument \code{strata} for other means of controlling grouping.
#' The last option for data input is a PLINK file in 
#' \code{tped/tfam} format (e.g. \code{data =  "data.assignment.tped"}). 
#' The first 2 columns of the \code{tfam} file will be used for the 
#' \code{strata} argument below, unless a new one is provided. 
#' Columns 1, 3 and 4 of the \code{tped} are discarded. The remaining columns 
#' correspond to the genotype in the format \code{01/04} 
#' where \code{A = 01, C = 02, G = 03 and T = 04}. For \code{A/T} format, use 
#' PLINK or bash to convert.
#' Use VCFTOOLS \url{http://vcftools.sourceforge.net/} with \code{--plink-tped} 
#' to convert very large VCF file. For \code{.ped} file conversion to 
#' \code{.tped} use PLINK \url{http://pngu.mgh.harvard.edu/~purcell/plink/} 
#' with \code{--recode transpose}.
#' 
#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by chromosome and/or locus and/or by snp.
#' The whitelist is in the working directory (e.g. "whitelist.txt").
#' de novo CHROM column with 'un' need to be changed to 1. 
#' Default \code{NULL} for no whitelist of markers. In the VCF, the column ID is
#' the LOCUS identification.

#' @param monomorphic.out (optional) For PLINK file, should the monomorphic 
#' markers present in the dataset be filtered out ? 
#' Default: \code{monomorphic.out = TRUE}.

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

#' @param snp.ld (optional) For VCF file only. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 3 options: \code{"random"} selection, \code{"first"} or
#' \code{"last"} SNP on the same short read/haplotype. 
#' Default: \code{snp.ld = NULL}.
#' Note that for other file type, use stackr package for haplotype file and 
#' create a whitelist, for plink and data frames, use PLINK linkage 
#' disequilibrium based SNP pruning option.
#' @param common.markers (optional) Logical. Default: \code{common.markers = TRUE}, 
#' will only keep markers in common (genotyped) between all the populations.


#' @param maf.thresholds (string, double, optional) String with 
#' local/populations and global/overall maf thresholds, respectively.
#' Default: \code{maf.thresholds = NULL}. 
#' e.g. \code{maf.thresholds = c(0.05, 0.1)} for a local maf threshold 
#' of 0.05 and a global threshold of 0.1. Available for VCF, PLINK and data frame 
#' files. Use stackr for haplotypes files.
#' @param maf.pop.num.threshold (integer, optional) When maf thresholds are used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. Default: \code{maf.pop.num.threshold = 1}
#' @param maf.approach (character, optional). By \code{maf.approach = "SNP"} or 
#' by \code{maf.approach = "haplotype"}.
#' The function will consider the SNP or ID/LOCUS/haplotype/read MAF statistics 
#' to filter the markers.
#' Default is \code{maf.approach = "SNP"}. The \code{haplotype} approach is 
#' restricted to VCF file.
#' @param maf.operator (character, optional) \code{maf.operator = "AND"} or default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").

#' @param filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. Default \code{assignment_data.txt}.
#' The number of markers used will be appended to the name of the file.

#' @param pop.levels (optional) A character string with your populations ordered.
#' @param pop.labels (optional) A character string for your populations labels.
#' If you need to rename sampling sites in \code{pop.levels} or combined sites/pop
#' into a different names, here is the place.
#' @param pop.id.start The start of your population id
#' in the name of your individual sample. Your individuals are identified 
#' in this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
#' then, \code{pop.id.start} = 5. If you didn't name your individuals
#' with the pop id in it, use the \code{strata} argument. 
#' @param pop.id.end The end of your population id
#' in the name of your individual sample. Your individuals are identified 
#' in this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
#' then, \code{pop.id.end} = 7. If you didn't name your individuals
#' with the pop id in it, use the \code{strata} argument.
#' @param strata (optional) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. Default: \code{strata = NULL}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping.

#' @param pop.select (string) Conduct the assignment analysis on a
#' selected list of populations. Default = \code{NULL} for no selection and keep
#' all population.
#' e.g. \code{pop.select = "QUE"} to select QUE population samples.
#' \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population samples.

#' @details You need to have either the \code{pop.id.start} and \code{pop.id.end}
#' or the \code{strata} argument, to identify your populations.

#' @return Depending on arguments selected, several files are written to the your
#' working directory or \code{folder}

#' @export
#' @rdname data_pruning
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @examples
#' \dontrun{
#' data_pruning(data = "plink.tped",
#' pop.levels = c("pop_1", "pop_2", "pop_3", "pop_4", "pop_5", "pop_6", "pop_7"),
#' monomorphic.out <- TRUE,
#' maf.thresholds <- c(0.05, 0.1),
#' maf.pop.num.threshold <- 1,
#' maf.approach <- "SNP",
#' maf.operator <- "OR",
#' common.markers <- TRUE,
#' pop.select <- c("pop_1", "pop_2", "pop_3")
#' ) 
#' }


#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, et al. 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. 
#' American Journal of Human Genetics. 2007; 81: 559–575. doi:10.1086/519795
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("ID", "#CHROM", "CHROM", "FORMAT", "INDIVIDUALS", "FORMAT_ID", "LOCUS",
      "POS", "REF", "ALT", "POP_ID", "READ_DEPTH", "ALLELE_DEPTH", "GL",
      "ERASE", "GT", "MARKERS", "QQ", "PQ", "N", "MAF_GLOBAL", "MAF_LOCAL",
      "ALLELES", "POP_ID", "GT", "INDIVIDUALS", "MARKERS", "POP_ID", "nal",
      "ALLELES_GROUP", "ALLELES", "N_IND_GENE", "P", "N", "nal_sq",
      "nal_sq_sum", "nal_sq_sum_nt", "npl", "het", "mho", "mhom", "dum",
      "dum1", "SSG", "ntal", "SSP", "ntalb", "SSi", "MSI", "sigw", "MSP",
      "siga", "sigb", "lsiga", "lsigb", "lsigw", "FST", "MARKERS",
      "MARKERS_ALLELES", "ALLELES", "POP_ID", "INDIVIDUALS", "filename",
      "ID", "KEEPER", "ASSIGN", "OTHERS", "CURRENT", "INFERRED",
      "SECOND_BEST_POP", "SCORE", "SECOND_BEST_SCORE", "NUMBER", "INDIVIDUALS_ALLELES",
      "MARKER_NUMBER", "MISSING_DATA", "TOTAL", "ASSIGNMENT_PERC",
      "MARKERS", "CURRENT", "INFERRED", "MISSING_DATA",
      "ITERATIONS", "METHOD", "TOTAL", "MEAN_i", "MEAN", "ASSIGNMENT_PERC",
      "SE", "MEDIAN", "MIN", "MAX", "QUANTILE25", "QUANTILE75", "SE_MIN",
      "SE_MAX", ".", "QUAL", "FILTER", "INFO", "pb", "SUBSAMPLE", "STRATA", 
      "sum.pop", "A1", "A2", "INDIVIDUALS_2", "Cnt", "Catalog ID", "GROUP",
      "COUNT", "MAX_COUNT_MARKERS", "hierarchy", "COL1", "COL3", "COL4"
    )
  )
}

data_pruning <- function(data,
                         whitelist.markers,
                         monomorphic.out,
                         blacklist.genotype,
                         blacklist.id,
                         pop.id.start, 
                         pop.id.end,
                         strata,
                         pop.levels,
                         pop.labels,
                         common.markers,
                         pop.select,
                         snp.ld,
                         maf.thresholds,
                         maf.pop.num.threshold,
                         maf.approach,
                         maf.operator,
                         filename
                         ) {
  
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (missing(whitelist.markers)) whitelist.markers <- NULL # no Whitelist
  if (missing(monomorphic.out)) monomorphic.out <- TRUE # remove monomorphic
  if (missing(blacklist.genotype)) blacklist.genotype <- NULL # no genotype to erase
  if (missing(snp.ld)) snp.ld <- NULL
  if (missing(maf.thresholds)) maf.thresholds <- NULL
  if (missing(maf.pop.num.threshold)) maf.pop.num.threshold <- 1
  if (missing(maf.approach)) maf.approach <- "SNP"
  if (missing(maf.operator)) maf.operator <- "OR"
  if (missing(blacklist.id)) blacklist.id <- NULL # No blacklist of ID
  
  if (missing(pop.levels)) pop.levels <- NULL
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels
  if (missing(pop.id.start)) pop.id.start <- NULL
  if (missing(pop.id.end)) pop.id.end <- NULL
  if (missing(strata)) strata <- NULL
  
  if (missing(pop.select)) pop.select <- NULL
  if (missing(common.markers)) common.markers <- TRUE
  
  if (missing(filename)) filename <- NULL

  # pop.information ---------------------------------------------------------
  
  if (is.null(strata) & is.null(pop.id.start) & is.null(pop.id.end)) {
    pop.info <- FALSE
  } else {
    pop.info <- TRUE
  }
  
  if (pop.info == FALSE & !is.null(pop.select)) {
    stop("For pop.select to work, you need population information: read arguments strata, pop.id.start and pop.id.end")
  }
  
  if (pop.info == FALSE & !is.null(common.markers)) {
    stop("For common.markers to work, you need population information: read arguments strata, pop.id.start and pop.id.end")
  }
  
  # File type detection ********************************************************
  data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
  
  if (identical(data.type, "##fileformat=VCF") | stri_detect_fixed(str = data, pattern = ".vcf")) {
    data.type <- "vcf.file"
    message("File type: VCF")
  }
  
  if (stri_detect_fixed(str = data, pattern = ".tped")) {
    data.type <- "plink.file"
    message("File type: PLINK")
    if (!file.exists(stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
      stop("Missing tfam file with the same prefix as your tped")
    }
  } 
  
  if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS")) {
    data.type <- "df.file"
    message("File type: data frame of genotypes")
  }
  
  if (stri_detect_fixed(str = data.type, pattern = "Catalog")) {
    data.type <- "haplo.file"
    message("File type: haplotypes from stacks")
    if (is.null(blacklist.genotype)) {
      stop("blacklist.genotype file missing. 
           Use stackr's missing_genotypes function to create this blacklist")
    }
  }
  
  # Create a filename to save the output files ********************************
  # Get date and time to have unique filenaming
  if (is.null(filename)) {
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    if (data.type == "vcf.file") {
      filename <- stri_replace_all_fixed(data, pattern = ".vcf", replacement = (stri_paste("_pruning_", file.date, ".vcf")), vectorize_all = FALSE)
    }
    
    if (data.type == "haplo.file") {
      filename <- stri_replace_all_fixed(data, pattern = ".tsv", replacement = (stri_paste("_pruning_", file.date, ".tsv")), vectorize_all = FALSE)
    }
    
    if (data.type == "df.file") {
      filename <- stri_paste(data, "_pruning_", file.date, ".tsv")
    }
    
    if (data.type == "plink.file") {
      filename.tped <- stri_replace_all_fixed(data, pattern = ".tped", replacement = (stri_paste("_pruning_", file.date, ".tped")), vectorize_all = FALSE)
      filename.tmap <- stri_replace_all_fixed(data, pattern = ".tped", replacement = (stri_paste("_pruning_", file.date, ".tmap")), vectorize_all = FALSE)
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
    if ("LOCUS" %in% columns.names.whitelist) {
      whitelist.markers$LOCUS <- as.character(whitelist.markers$LOCUS)
    }
    if ("POS" %in% columns.names.whitelist) {
      whitelist.markers$POS <- as.character(whitelist.markers$POS)
    }
  }
  
  if (data.type == "haplo.file") {
    whitelist.markers <- select(.data = whitelist.markers, LOCUS)
    columns.names.whitelist <- colnames(whitelist.markers)
  }
  
  # Import blacklist id ********************************************************
  if (is.null(blacklist.id)) { # No blacklist of ID
    message("Blacklisted individuals: no")
  } else { # With blacklist of ID
    message("Blacklisted individuals: yes")
    blacklist.id <- read_tsv(blacklist.id, col_names = TRUE)
  }
  
  # Import data ****************************************************************
  if (data.type == "vcf.file") { # VCF
    message("Importing the VCF...")
    
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
        CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1"),
        POS = as.character(POS),
        LOCUS = as.character(LOCUS)
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
    
    if (pop.info == TRUE) {
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
    } else {
      input <- input %>% mutate(INDIVIDUALS =  as.character(INDIVIDUALS))
    }
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # Tidy VCF
    message("Tidying the vcf...")
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
  } # End import VCF
  
  if (data.type == "plink.file") { # PLINK
    message("Importing the PLINK files...")
    strata.df <- data.table::fread(
      input = stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE),
      sep = " ", 
      header = FALSE, 
      stringsAsFactors = FALSE,
      verbose = FALSE,
      select = c(1,2),
      colClasses=list(character = c(1,2)),
      col.names = c("POP_ID", "INDIVIDUALS"),
      showProgress = TRUE, 
      data.table = FALSE)
    
    # remove "_" in individual name and replace with "-"
    strata.df$INDIVIDUALS <- stri_replace_all_fixed(str = strata.df$INDIVIDUALS, pattern = "_", replacement = "-", vectorize_all = TRUE)
    
    tped.header.prep <- strata.df %>% 
      select(INDIVIDUALS) %>%
      mutate(NUMBER = seq(1,n())) %>%
      mutate(ALLELE1 = rep("A1", n()), ALLELE2 = rep("A2", n())) %>%
      tidyr::gather(ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, NUMBER)) %>%
      arrange(NUMBER) %>% 
      select(-ALLELES_GROUP) %>% 
      tidyr::unite(INDIVIDUALS_ALLELES, c(INDIVIDUALS, ALLELES), sep = "_", remove = FALSE) %>% 
      arrange(NUMBER) %>% 
      mutate(NUMBER = seq(from = (1 + 4), to = n() + 4)) %>% 
      select(-ALLELES)
    
    tped.header.names <- c("LOCUS", tped.header.prep$INDIVIDUALS_ALLELES)
    tped.header.integer <- c(2, tped.header.prep$NUMBER)
    
    if (!is.null(blacklist.id)) { # using the blacklist of individuals
      whitelist.id <- tped.header.prep %>% 
        anti_join(blacklist.id, by = "INDIVIDUALS") %>% 
        arrange(NUMBER)
      tped.header.names <- c("LOCUS", whitelist.id$INDIVIDUALS_ALLELES)
      tped.header.integer <- c(2, whitelist.id$NUMBER)
    }
    
    input <- data.table::fread( # import PLINK
      input = data, 
      sep = " ", 
      header = FALSE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      select = tped.header.integer,
      col.names = tped.header.names,
      showProgress = TRUE,
      data.table = FALSE) %>% 
      as_data_frame() %>% 
      mutate(LOCUS = as.character(LOCUS))
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Using the argument strata if provided to replace the current one
    if (!is.null(strata)) {
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA)
    }
    
    # Make tidy
    message("Tidying the PLINK file and integrating the tfam/strata file, for large dataset this may take several minutes...")
    input <- input %>% 
      tidyr::gather(key = INDIVIDUALS_ALLELES, value = GT, -LOCUS) %>%
      mutate(INDIVIDUALS = stri_replace_all_fixed(str = INDIVIDUALS_ALLELES, pattern = c("_A1", "_A2"), replacement = "", vectorize_all = FALSE)) %>% 
      left_join(strata.df, by = "INDIVIDUALS") %>% 
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =TRUE),
        GT = stri_pad_left(str = GT, width = 3, pad = "0")
      )
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # removing untyped markers across all-pop
    remove.missing.gt <- input %>%
      select(LOCUS, GT) %>%
      filter(GT != "000")
    
    untyped.markers <- n_distinct(input$LOCUS) - n_distinct(remove.missing.gt$LOCUS)
    if (untyped.markers > 0) {
      message(paste0("Number of marker with 100 % missing genotypes: ", untyped.markers))
      input <- suppressWarnings(
        semi_join(input, 
                  remove.missing.gt %>% 
                    select(LOCUS) %>% 
                    distinct(LOCUS), 
                  by = "LOCUS")
      )
    }
    
    # Removing monomorphic markers
    if (monomorphic.out == TRUE) {
      message("Removing monomorphic markers...")
      mono.markers <- remove.missing.gt %>%
        group_by(LOCUS, GT) %>%
        distinct %>% 
        group_by(LOCUS) %>%
        tally %>% 
        filter(n == 1) %>% 
        select(LOCUS) %>% 
        arrange(LOCUS)
      
      # Remove the markers from the dataset
      input <- anti_join(input, mono.markers, by = "LOCUS")
      message(paste0("Number of monomorphic markers removed: ", n_distinct(mono.markers$LOCUS)))
    }
    # Unused objects
    tped.header.prep <- NULL
    tped.header.integer <- NULL
    tped.header.names <- NULL
    remove.missing.gt <- NULL
    mono.markers <- NULL
  } # End import PLINK
  
  if (data.type == "df.file") { # DATA FRAME OF GENOTYPES
    message("Importing the data frame")
    input <- data.table::fread(
      input = data, 
      sep = "\t", 
      header = TRUE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      showProgress = TRUE,
      data.table = FALSE
    ) %>% 
      as_data_frame() %>% 
      tidyr::gather(key = LOCUS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% 
      mutate(GT = stri_pad_left(str = GT, width = 6, pad = "0"))
    
    if (!is.null(pop.levels)) {
      input <- input %>% 
        mutate(POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered =T))
    }
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      input <- suppressWarnings(anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    # Using the argument strata if provided to replace the current one
    if (is.null(strata)) {
      strata.df <- input %>% select(INDIVIDUALS, POP_ID) %>% distinct(INDIVIDUALS)
    } else {
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA)
      input <- left_join(input, strata.df, by = "INDIVIDUALS")
    }
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
  } # End import data frame of genotypes
  
  if (data.type == "haplo.file") { # Haplotype file
    message("Importing the stacks haplotype file")
    input <- data.table::fread(
      input = data, 
      sep = "\t", 
      header = TRUE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      showProgress = TRUE,
      data.table = FALSE
    ) %>% 
      as_data_frame() %>% 
      select(-Cnt) %>% 
      rename(LOCUS = `Catalog ID`) %>%
      tidyr::gather(INDIVIDUALS, GT, -LOCUS) %>% 
      mutate(LOCUS = as.character(LOCUS))
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      input <- suppressWarnings(anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    if (pop.info == TRUE) {
      if (is.null(strata)){
        input <- input %>%
          mutate( # Make population ready
            POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
            INDIVIDUALS =  as.character(INDIVIDUALS)
          )
        if (!is.null(pop.levels)) {
          input <- input %>% 
            mutate(
              POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered = TRUE)
            )
        } else {
          input <- input %>% 
            mutate(POP_ID = factor(POP_ID))
        }
        
        
      } else { # Make population ready with the strata provided
        strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
          rename(POP_ID = STRATA)
        
        input <- input %>%
          mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
          left_join(strata.df, by = "INDIVIDUALS")
        
        if (!is.null(pop.levels)) {
          input <- input %>% 
            mutate(POP_ID = factor(POP_ID, levels = unique(pop.labels), ordered =TRUE))
          
        } else {
          input <- input %>% 
            mutate(POP_ID = factor(POP_ID))
        }
      }
    }
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
  } # End import haplotypes file
  
  # Blacklist genotypes ********************************************************
  if (is.null(blacklist.genotype)) { # no Whitelist
    message("Erasing genotype: no")
  } else {
    message("Erasing genotype: yes")
    blacklist.genotype <- read_tsv(blacklist.genotype, col_names = TRUE)
    blacklist.genotype <- suppressWarnings(
      plyr::colwise(as.character, exclude = NA)(blacklist.genotype)
    )
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    
    if ("CHROM" %in% columns.names.blacklist.genotype) {
      columns.names.blacklist.genotype$CHROM <- as.character(columns.names.blacklist.genotype$CHROM)
    }
    
    if (data.type == "haplo.file") {
      blacklist.genotype <- select(.data = blacklist.genotype, INDIVIDUALS, LOCUS)
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
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
      if (data.type == "vcf.file"){
        whitelist.markers.ind <- input %>% select(CHROM, LOCUS, POS, INDIVIDUALS) %>% distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      } else {
        whitelist.markers.ind <- input %>% select(LOCUS, INDIVIDUALS) %>% distinct(LOCUS, INDIVIDUALS)
      }
      
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
        mutate(ERASE = stri_replace_na(str = ERASE, replacement = "ok"))
    )
    
    if (data.type == "vcf.file") {
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "./.", GT)) %>% 
        select(-ERASE)
    }
    
    if (data.type == "plink.file") {
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "000", GT)) %>% 
        select(-ERASE)
    }
    
    if (data.type == "df.file") {
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "000000", GT)) %>% 
        select(-ERASE)
    }
    
    if (data.type == "haplo.file") {
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "-", GT)) %>% 
        select(-ERASE)
    }
    
  } # End erase genotypes
  
  # dump unused object
  blacklist.id <- NULL
  whitelist.markers <- NULL
  whitelist.markers.ind <- NULL
  blacklist.genotype <- NULL
  
  # LD control... keep only 1 SNP per haplotypes/reads (optional) ************
  if (!is.null(snp.ld)) {
    if (data.type != "vcf.file") {
      stop("snp.ld is only available for VCF file, use stackr package for 
             haplotype file and create a whitelist, for other file type, use 
             PLINK linkage disequilibrium based SNP pruning option")
    }
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
  
  # Unique markers id: for VCF combine CHROM, LOCUS and POS into MARKERS *****
  # For other type of file change LOCUS to MARKERS
  if (data.type != "vcf.file") {
    input <- input %>% rename(MARKERS = LOCUS)
  } else {
    input <- input %>%
      arrange(CHROM, LOCUS, POS) %>%
      tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
  } # End Unique markers id
  
  # Markers in common between all populations (optional) *********************
  if (common.markers == TRUE) { # keep only markers present in all pop
    message("Using markers common in all populations:")
    pop.number <- n_distinct(input$POP_ID)
    
    if (data.type == "vcf.file") pop.filter <- input %>% filter(GT != "./.")
    if (data.type == "plink.file") pop.filter <- input %>% filter(GT != "000")
    if (data.type == "df.file") pop.filter <- input %>% filter(GT != "000000")
    if (data.type == "haplo.file") pop.filter <- input %>% filter(GT != "-")
    
    pop.filter <- pop.filter %>% 
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
  
  # Minor Allele Frequency filter ********************************************
  if (!is.null(maf.thresholds)) { # with MAF
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
    message("MAF filter: yes")
    
    if (data.type == "vcf.file") {
      maf.local <- input %>%
        filter(GT != "./.") %>%
        group_by(MARKERS, POP_ID, REF, ALT) %>%
        summarise(
          N = as.numeric(n()),
          PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
          QQ = as.numeric(length(GT[GT == "1/1"]))
        ) %>%
        mutate(MAF_LOCAL = ((QQ * 2) + PQ) / (2 * N))
      
      maf.global <- maf.local %>%
        group_by(MARKERS) %>%
        summarise_each_(funs(sum), vars = c("N", "PQ", "QQ")) %>%
        mutate(MAF_GLOBAL = ((QQ * 2) + PQ) / (2 * N)) %>%
        select(MARKERS, MAF_GLOBAL)
      
      maf.data <- maf.global %>%
        left_join(maf.local, by = c("MARKERS")) %>%
        select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
      
      maf.local <- NULL
      maf.global <- NULL
    } # end maf calculations with vcf
    
    if (data.type == "plink.file" | data.type == "df.file") {
      message("Calculating global and local MAF, this may take some time on large data set")
      
      # For data frame we split the alleles here to prep for MAF
      if (data.type == "df.file") { # for data frame of genotypes
        maf.data <- input %>% 
          tidyr::separate(data = ., col = GT, into = .(A1, A2), sep = 3, remove = TRUE) %>% 
          tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
          select(MARKERS, GT, POP_ID) %>% 
          filter(GT != "000")
      }
      
      if (data.type == "plink.file") { # For PLINK and common code below
        maf.data <- input %>%
          select(MARKERS, GT, POP_ID) %>% 
          filter(GT != "000")
      }
      
      maf.data <- maf.data %>%
        group_by(MARKERS, GT, POP_ID) %>%
        tally %>%
        arrange(MARKERS, GT) %>% 
        group_by(MARKERS, GT) %>%
        mutate(sum.pop = sum(n)) %>% 
        group_by(MARKERS) %>%
        mutate(MAF_GLOBAL = min(sum.pop)/sum(n)) %>% 
        group_by(MARKERS, POP_ID) %>%
        mutate(MAF_LOCAL = n/sum(n)) %>% 
        arrange(MARKERS, POP_ID, GT) %>% 
        group_by(MARKERS, POP_ID) %>% 
        filter(n == min(n)) %>% 
        distinct(MARKERS, POP_ID) %>% 
        select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
    }# end maf calculations with PLINK or data frame of genotypes
    
    if (data.type == "haplo.file") {
      stop("MAF filtering is only available for haplotype file, use stackr
             package and update your whitelist")
    }
    
    write_tsv(x = maf.data, 
              path = "maf.data.tsv", 
              col_names = TRUE, 
              append = FALSE
    )
    message("The MAF table was written in your folder")
    
    # # update the vcf with the maf info
    # input <- full_join(input, maf.data, by = c("MARKERS", "POP_ID"))
    if (maf.approach == "haplotype") {
      if (data.type != "vcf.file") {
        stop("The haplotype approach during MAF filtering is for VCF files only")
      }
      vcf.maf <- tidyr::separate(data = maf.data, 
                                 col = MARKERS, 
                                 into = c("CHROM", "LOCUS", "POS"), 
                                 sep = "_", 
                                 remove = FALSE, 
                                 extra = "warn"
      )
      
      if (maf.operator == "OR") {
        vcf.maf <- maf.data %>%
          group_by(LOCUS, POP_ID) %>%
          summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
          group_by(LOCUS) %>%
          tally() %>%
          filter(n >= maf.pop.num.threshold) %>%
          select(LOCUS) %>%
          left_join(input, by = "LOCUS") %>%
          arrange(LOCUS, POP_ID)
      } else { # AND operator between local and global maf
        vcf.maf <- maf.data %>%
          group_by(LOCUS, POP_ID) %>%
          summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
          group_by(LOCUS) %>%
          tally() %>%
          filter(n >= maf.pop.num.threshold) %>%
          select(LOCUS) %>%
          left_join(input, by = "LOCUS") %>%
          arrange(LOCUS, POP_ID)
      }
      vcf.maf <- vcf.maf %>% select(-c(CHROM, LOCUS, POS))
    } # end maf haplotype approach
    
    if (maf.approach == "SNP") { # SNP approach
      if (maf.operator == "OR") {
        vcf.maf <- maf.data %>%
          group_by(MARKERS, POP_ID) %>%
          summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
          group_by(MARKERS) %>%
          tally() %>%
          filter(n >= maf.pop.num.threshold) %>%
          select(MARKERS) %>%
          left_join(input, by = "MARKERS") %>%
          arrange(MARKERS, POP_ID)
      } else { # AND operator between local and global maf
        vcf.maf <- maf.data %>%
          group_by(MARKERS, POP_ID) %>%
          summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
          group_by(MARKERS) %>%
          tally() %>%
          filter(n >= maf.pop.num.threshold) %>%
          select(MARKERS) %>%
          left_join(input, by = "MARKERS") %>%
          arrange(MARKERS, POP_ID)
      }
    } # end maf snp approach
    
    
    message(stri_join("The number of MARKERS removed by the MAF filters = ", 
                      n_distinct(input$MARKERS)-n_distinct(vcf.maf$MARKERS), "\n", 
                      "The number of MARKERS before -> after the MAF filters: ", 
                      n_distinct(input$MARKERS)," -> ", n_distinct(vcf.maf$MARKERS), 
                      " MARKERS"))
    
    input <- vcf.maf
    
    # unused object
    vcf.maf <- NULL 
    maf.data <- NULL
  } # End of MAF filters
  
  # back to original format ****************************************************
  if (data.type == "vcf.file") {
    info.field <- suppressWarnings(
      input %>% 
        group_by(MARKERS) %>%
        filter(GT != "./.") %>% 
        tally %>% 
        mutate(INFO = stri_paste("NS=", n, sep = "")) %>% 
        select(-n)
    )
    
     output <- suppressWarnings(
       left_join(input, info.field, by = "MARKERS") %>% 
        group_by(MARKERS, INFO) %>% 
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT) %>%
        ungroup() %>% 
        tidyr::separate(MARKERS, c("CHROM", "ID", "POS"), sep = "_", extra = "warn") %>%
        mutate(
          ID = as.numeric(ID),
          POS = as.numeric(POS),
          QUAL = rep(".", n()),
          FILTER = rep("PASS", n()),
          FORMAT = rep("GT", n())
        ) %>% 
        arrange(CHROM, ID, POS) %>% 
        select('#CHROM' = CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, everything())
    )
     
     # File format
     file.format <- "##fileformat=VCFv4.0"
     file.format <- as.data.frame(file.format)
     write.table(x = file.format, file = filename, sep = " ", append = FALSE, col.names = FALSE, row.names = FALSE, quote = FALSE)
     
     # File date
     file.date <- stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")
     file.date <- stri_paste("##fileDate=", file.date, sep = "")
     write.table(x = file.date, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
     
     # Source
     file.source <- as.data.frame(stri_paste("##source=stackr v.", packageVersion("stackr"), sep = ""))
     write.table(x = file.source, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
     
     # Info field 1
     info1 <- '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">'
     info1 <- as.data.frame(info1)
     write.table(x = info1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
     
     # Format field 1
     format1 <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
     format1 <- as.data.frame(format1)
     write.table(x = format1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
     
     # Write the prunned vcf to the file
     suppressWarnings(
       write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE)
     )
  }
  if (data.type == "haplo.file") {
    cnt.field <- suppressWarnings(
      input %>% 
        group_by(MARKERS) %>%
        filter(GT != "-") %>% 
        tally %>%
        rename(Cnt = n) %>% 
        arrange(as.integer(MARKERS))
    )
    
    output <- suppressWarnings(
      left_join(input, cnt.field, by = "MARKERS") %>% 
        group_by(MARKERS) %>% 
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT) %>%
        ungroup() %>% 
        arrange(as.integer(MARKERS)) %>% 
        rename(`Catalog ID` = MARKERS)
    )

    # Write the prunned vcf to the file
    suppressWarnings(
      write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE)
    )
  } # end haplotypes file
  
  if (data.type == "plink.file") {
    # to create a PLINK tped and tfam
    tped <- input %>% 
      mutate(
        COL1 = rep("0", n()),
        COL3 = rep("0", n()),
        COL4 = rep("0", n())
      ) %>% 
      arrange(POP_ID, INDIVIDUALS, INDIVIDUALS_ALLELES) %>% 
      select(COL1, MARKERS, COL3, COL4, INDIVIDUALS_ALLELES, GT) %>%
      group_by(COL1, MARKERS, COL3, COL4) %>% 
      tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GT) %>% 
      arrange(MARKERS)
    
    write_delim(x = tped, path = filename.tped, col_names = FALSE, delim = " ")
    
    # Create a tfam file
    tfam <- strata.df %>% 
      select(POP_ID, INDIVIDUALS) %>% 
      arrange(POP_ID, INDIVIDUALS) %>% 
      mutate(
        COL3 = rep("0",n()),
        COL4 = rep("0",n()),
        COL5 = rep("0",n()),
        COL6 = rep("-9",n())
      )
    
    # the tfam must have the same name as the tped 
    write_delim(x = tfam, path = filename.tmap, col_names = FALSE, delim = " ")
  } # end plink
  if (data.type == "df.file") {
    output <- input %>% 
      select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
      arrange(MARKERS, POP_ID, INDIVIDUALS) %>% 
      group_by(POP_ID, INDIVIDUALS) %>% 
      tidyr::spread(data = ., key = MARKERS, value = GT)
      
    # Write the prunned vcf to the file
    suppressWarnings(
      write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE)
    )
      }# end df

} # End data_pruning



