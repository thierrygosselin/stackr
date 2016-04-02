#' @name vcf2dadi
#' @title Create a \code{dadi} SNP input file from a STACKS vcf file.
#' @description This function will create a \code{dadi} SNP input file and a 
#' joint allele frequency spectrum file (folded) using STACKS \code{batch_x.vcf} file
#' or an imputed vcf file see \code{\link[stackr]{vcf_imputation}}. If your VCF 
#' is not filtered, you can supply the function a whitelist of loci and a 
#' blacklist of individuals.

#' @param data The VCF file created by STACKS or imputed VCF file created
#' by \code{\link[stackr]{vcf_imputation}}.

#' @param whitelist.markers (optional) A whitelist containing CHROM (character or integer) and/or LOCUS (integer) and/or 
#' POS (integer) columns header. To filter by CHROM and/or locus and/or by snp.
#' The whitelist is in the directory (e.g. "whitelist.txt"). de novo CHROM column with 'un' need to be changed to 1.

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

#' @param pop.select (string) Conduct the assignment analysis on a
#' selected list of populations. Default = \code{NULL} for no selection and keep
#' all population.
#' e.g. \code{pop.select = "QUE"} to select QUE population samples.
#' \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population samples.

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

#' @param fasta.outgroup (optional) The fasta output file from STACKS. This file is 
#' required to use an outgroup. Default: \code{fasta.outgroup = NULL}.

#' @param fasta.ingroup (optional) The fasta output file from STACKS. This file is 
#' required to use with an outgroup. Leave empty if no outgroup. 
#' Default: \code{fasta.ingroup = NULL}.

#' @param sumstats.outgroup (optional) The sumstats output file from STACKS 
#' when running STACKS for the outgroup fasta file. This file is 
#' required to use an outgroup. Default: \code{sumstats.outgroup = NULL}.

#' @param sumstats.ingroup (optional) The sumstats output file from STACKS
#' when running STACKS for the ingroup fasta file.This file is 
#' required to use with an outgroup. Leave empty if no outgroup. 
#' Default: \code{sumstats.ingroup = NULL}.


#' @param dadi.input.filename (optional) Name of the \code{dadi} SNP input file 
#' written to the working directory. e.g. \code{dadi.file.txt}. 
#' Default use date and time to make the file. If used, the file extension
#' need to finish with \code{.txt}.

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
#' @details The imputations using Random Forest requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set 
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals 
#' will require 15 min.

#' @export
#' @rdname vcf2dadi
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @examples
#' \dontrun{
#' vcf2dadi(
#' data = "batch_1.vcf", 
#' whitelist.markers = "whitelist.loci.txt",
#' pop.id.start = 5,
#' pop.id.end = 7, 
#' pop.levels = c("PAN", "COS"),
#' common.markers = TRUE, 
#' fasta.ingroup = "batch_1.ingroup.fa", 
#' fasta.outgroup = "batch_1.outgroup.fa", 
#' sumstats.ingroup = "batch_1.sumstats.ingroup.tsv", 
#' sumstats.outgroup = "batch_1.sumstats.outgroup.tsv"
#' )
#' 
#' With Imputations:
#' vcf2dadi(
#' data = "batch_1.vcf", 
#' whitelist.markers = "whitelist.loci.txt",
#' pop.id.start = 5,
#' pop.id.end = 7, 
#' pop.levels = c("PAN", "COS"),
#' common.markers = TRUE, 
#' fasta.ingroup = "batch_1.ingroup.fa", 
#' fasta.outgroup = "batch_1.outgroup.fa", 
#' sumstats.ingroup = "batch_1.sumstats.ingroup.tsv", 
#' sumstats.outgroup = "batch_1.sumstats.outgroup.tsv",
#' imputation.method = "max", 
#' impute = "allele", 
#' imputations.group = "populations", 
#' num.tree = 100, 
#' iteration.rf = 10, 
#' split.number = 100, 
#' verbose = FALSE, 
#' parallel.core = 8
#' )
#' }

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @references Gutenkunst RN, Hernandez RD, Williamson SH, Bustamante CD (2009)
#' Inferring the Joint Demographic History of Multiple Populations 
#' from Multidimensional SNP Frequency Data (G McVean, Ed,). 
#' PLoS genetics, 5, e1000695.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and
#' Anne-Laure Ferchaud \email{annelaureferchaud@@gmail.com} 

# remove NOTE about no visible binding for global variable during Build ------
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("#CHROM","QUAL", "FILTER", "INFO", "FORMAT", "FORMAT_ID",
                           "ID", "REF", "ALT", "READ_DEPTH", 'ALLELE_DEPTH', 'GL', 'MARKERS', 
                           "GT", "PP", "QQ", "PQ", "MAF", "fasta", "ALLELE_GROUP",
                           "Allele1", "Allele2", "POP", "IN_GROUP", "OUT_GROUP",
                           "ID.FILTER", "ANCESTRAL", "SEQUENCES", "GARBAGE", 
                           "Chr", "Locus", "BP", "Col", "SNP_READ_POS", "FASTA_REF",
                           "contains", "Locus ID", "i", "m"
  )
  )
}

vcf2dadi <- function(
  data,
  whitelist.markers = NULL, 
  blacklist.id = NULL,
  pop.id.start, 
  pop.id.end,
  pop.levels,
  pop.labels,
  strata = NULL,
  pop.select = NULL,
  blacklist.genotype = NULL,
  snp.ld = NULL,
  common.markers = FALSE,
  fasta.ingroup = NULL,
  fasta.outgroup = NULL,
  sumstats.ingroup = NULL,
  sumstats.outgroup = NULL,
  dadi.input.filename = NULL,
  imputation.method = FALSE,
  impute,
  imputations.group = "populations",
  num.tree = 100,
  iteration.rf = 10,
  split.number = 100,
  verbose = FALSE,
  parallel.core
){
  
  # dadi unicode character: \u2202
  
  if (missing(data)) stop("Input file missing")
  if (missing(whitelist.markers)) whitelist.markers <- NULL # no Whitelist
  if (missing(blacklist.id)) blacklist.id <- NULL # No blacklist of ID
  if (missing(strata)) strata <- NULL
  if (missing(pop.select)) pop.select <- NULL
  if (missing(blacklist.genotype)) blacklist.genotype <- NULL # no genotype to erase
  if (missing(pop.levels)) stop("pop.levels required")
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels
  if (missing(common.markers)) common.markers <- FALSE
  if (missing(snp.ld)) snp.ld <- NULL
  if (missing(fasta.outgroup)) fasta.outgroup <- NULL
  if (missing(fasta.ingroup)) fasta.ingroup <- NULL
  if (missing(sumstats.outgroup)) sumstats.outgroup <- NULL
  if (missing(sumstats.ingroup)) sumstats.ingroup <- NULL
  
  if (missing(dadi.input.filename)) dadi.input.filename <- NULL
  
  
  if (missing(imputation.method)) imputation.method <- FALSE
  if (imputation.method != FALSE & missing(impute)) stop("impute argument is necessary")
  if (imputation.method == FALSE & missing(impute)) impute <- NULL
  if (missing(imputations.group)) imputations.group <- "populations"
  if (missing(num.tree)) num.tree <- 100
  if (missing(iteration.rf)) iteration.rf <- 10
  if (missing(split.number)) split.number <- 100
  if (missing(verbose)) verbose <- FALSE
  if (missing(parallel.core)) parallel.core <- detectCores()-1
  
  if (imputation.method == "FALSE") {
    message("vcf2dadi: no imputation...")
  } else {
    message("vcf2dadi: with imputations...")
  }
  # File type detection ********************************************************
  data.type <- read_lines(file = data, skip = 2, n_max = 1)
  
  if (stri_detect_fixed(str = data.type, pattern = "##source=stackr")) {
    data.type <- "vcf.imputed"
    message("File type: imputed VCF")
  } else {
    data.type <- "vcf.raw"
    message("File type: VCF non imputed")
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
  
  # no imputation
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
  
  # Conversion -----------------------------------------------------
  # Tidy VCF
  message("Tidying the vcf ...")
  
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
  
  # Compute count and Minor Allele Frequency
  # MAF = the ALT allele in the VCF
  message("Computing the Allele Frequency Spectrum")
  
  input.count <- input %>% 
    group_by(MARKERS, POP_ID, REF, ALT) %>%
    summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT[GT == "0/0"])),
      PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
      QQ = as.numeric(length(GT[GT == "1/1"]))
    ) %>%
    mutate(MAF = ((QQ*2) + PQ)/(2*N)) 
  
  # Function to make dadi input  data format ***********************************
  message("Preparing \u2202a\u2202i input SNP data format")
  write_dadi <- function(input, write.imputation, ...){
    # input <- input.count.imp # testing
    
    if (is.null(fasta.ingroup)) {
      dadi.input <- suppressWarnings(
        input %>%
          group_by(MARKERS, POP_ID, REF, ALT) %>%
          mutate(
            A1 = (2 * PP) + PQ,
            A2 = (2 * QQ) + PQ
          ) %>%
          ungroup() %>% 
          select(POP_ID, Allele1 = REF, A1, Allele2 = ALT, A2, MARKERS) %>%
          tidyr::gather(ALLELE_GROUP, COUNT, -c(POP_ID, MARKERS, Allele1, Allele2)) %>% 
          tidyr::unite(POP, POP_ID, ALLELE_GROUP, sep = "_") %>% 
          group_by(MARKERS, Allele1, Allele2) %>% 
          tidyr::spread(data = ., key = POP, value = COUNT) %>%
          mutate(
            IN_GROUP = rep("---", n()), #version 2
            OUT_GROUP = rep("---", n())
          ) %>% 
          select(IN_GROUP, OUT_GROUP, Allele1, contains("A1"), Allele2, contains("A2"), MARKERS) %>% 
          arrange(MARKERS)
      )
    } # no fasta, no outgroup
    if (!is.null(fasta.ingroup)) {# With outgroup and ingroup fasta file
      message("using outgroup info")
      # Get the list of ref. allele in the vcf of the ingroup
      ref.allele.vcf.ingroup <- input %>% 
        ungroup() %>% 
        select(MARKERS, REF) %>%
        distinct(MARKERS, REF) %>%
        arrange(MARKERS) %>% 
        tidyr::separate(MARKERS, c("CHROM", "LOCUS", "POS"), sep = "_") %>%
        distinct(CHROM, LOCUS, POS, REF) %>%
        mutate(
          CHROM = as.character(stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>%
        arrange(CHROM, LOCUS, POS)
      # View(ref.allele.vcf.ingroup)
      
      # keep the list of markers
      markers <- ref.allele.vcf.ingroup %>% ungroup %>% select(-REF)
      
      # import out- and in- group fasta files ********************************
      fasta.outgroup <- data.table::fread(
        input = fasta.outgroup,
        sep = "\t",
        stringsAsFactors = FALSE,
        showProgress = FALSE,
        verbose = FALSE, 
        header = FALSE,
        col.names = "LOCUS"
      ) %>%
        mutate(
          ID.FILTER = rep(c("LOCUS", "SEQ"), times = n()/2),
          ANCESTRAL = rep("outgroup", times = n())
        )
      
      fasta.ingroup <- data.table::fread(
        input = fasta.ingroup,
        sep = "\t",
        stringsAsFactors = FALSE,
        showProgress = FALSE,
        verbose = FALSE, 
        header = FALSE,
        col.names = "LOCUS"
      ) %>%
        mutate(
          ID.FILTER = rep(c("LOCUS", "SEQ"), times = n()/2),
          ANCESTRAL = rep("ingroup", times = n())
        )
      
      fasta.data <- bind_rows(fasta.outgroup, fasta.ingroup)
      
      fasta.ingroup <- NULL
      fasta.outgroup <- NULL
      
      message("Preparing fasta sequences")
      fasta.prep <- fasta.data %>% 
        filter(ID.FILTER == "LOCUS") %>% 
        select(-ID.FILTER) %>% 
        bind_cols(fasta.data %>% 
                    filter(ID.FILTER == "SEQ") %>% 
                    select(-ID.FILTER, -ANCESTRAL, SEQUENCES = LOCUS)
        ) %>% 
        tidyr::separate(LOCUS, c("LOCUS", "ALLELE"), sep = "_Sample_", extra = "warn" ) %>%
        tidyr::separate(ALLELE, c("GARBAGE", "ALLELE"), sep = "_Allele_", extra = "warn" ) %>%
        tidyr::separate(ALLELE, c("ALLELE", "INDIVIDUALS"), sep = " ", extra = "warn" ) %>% 
        mutate(LOCUS = as.integer(stri_replace_all_fixed(LOCUS, pattern = ">CLocus_", replacement = "", vectorize_all = FALSE))) %>%
        select(-GARBAGE, -INDIVIDUALS) %>% 
        distinct(LOCUS, ALLELE, ANCESTRAL, SEQUENCES) 
      
      fasta.data <- NULL
      
      # import out- and in- group sumstats files*****************************
      
      # Note: only one sumstats might be necessary to have the info needed, need checking
      sumstats.outgroup <- data.table::fread(
        input = sumstats.outgroup,
        sep = "\t",
        stringsAsFactors = FALSE, 
        skip = "Batch ID",
        select = c("Chr", "Locus ID", "BP", "Col"),
        showProgress = TRUE,
        verbose = FALSE
      ) %>% 
        as_data_frame() %>%
        select(CHROM = Chr, LOCUS = `Locus ID`, POS = BP, SNP_READ_POS = Col) %>%
        mutate(
          CHROM = as.character(stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>% 
        distinct(CHROM, LOCUS, SNP_READ_POS) %>% 
        semi_join(markers, by = c("CHROM", "LOCUS", "POS")) %>% 
        arrange(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        mutate(ANCESTRAL = rep("outgroup", times = n())) %>%
        # the SNP position on the read is +1 for the fasta file
        mutate(SNP_READ_POS = SNP_READ_POS + 1)
      
      sumstats.ingroup <- data.table::fread(
        input = sumstats.ingroup,
        sep = "\t",
        stringsAsFactors = FALSE, 
        skip = "Batch ID",
        select = c("Chr", "Locus ID", "BP", "Col"),
        showProgress = TRUE,
        verbose = FALSE
      ) %>% 
        as_data_frame() %>%
        select(CHROM = Chr, LOCUS = `Locus ID`, POS = BP, SNP_READ_POS = Col) %>%
        mutate(
          CHROM = as.character(stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>% 
        distinct(CHROM, LOCUS, SNP_READ_POS) %>% 
        semi_join(markers, by = c("CHROM", "LOCUS", "POS")) %>% 
        arrange(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        mutate(ANCESTRAL = rep("ingroup", times = n())) %>%
        # the SNP position on the read is +1 for the fasta file
        mutate(SNP_READ_POS = SNP_READ_POS + 1)
      
      markers <- NULL # remove unused object
      
      # When REF and ALT allele are equal in number, this can be problematic between ANCESTRAL group of alleles.
      # For those loci, we will set them based on both group (out- and in-group), so that there is no differences
      
      # THINKING: when ingroup REF and ALT equal, set to REF in VCF, and EQUAL in outgroup
      # when outgroup REF and ALT equal, set the REF based on the ingroup VCF
      
      max.length.read <- stri_length(str = fasta.prep[1,4])
      
      message("combining information from the fasta and sumstats file")
      ref.allele.vcf.ingroup.fasta <- ref.allele.vcf.ingroup %>% 
        # The sumstats is used to get the position of the markers along the sequence read
        left_join(
          sumstats.ingroup %>% 
            select(-ANCESTRAL)
          , by = c("CHROM", "LOCUS", "POS")
        ) %>% 
        mutate(SNP_READ_POS = stri_replace_na(SNP_READ_POS, replacement = "NA")) %>%
        filter(SNP_READ_POS != "NA") %>%
        mutate(SNP_READ_POS = as.integer(SNP_READ_POS)) %>% 
        arrange(CHROM, LOCUS, POS) %>%
        ungroup() %>% 
        # get the fasta sequence for those LOCUS...
        left_join(
          fasta.prep %>% 
            filter(ANCESTRAL == "ingroup") %>% 
            select (-ALLELE, -ANCESTRAL)
          , by = "LOCUS"
        ) %>%
        mutate(SNP_READ_POS = stri_replace_na(SNP_READ_POS, replacement = "NA")) %>%
        filter(SNP_READ_POS != "NA") %>%
        # get the flanking bases, left and right, of the SNP
        mutate(
          SNP_READ_POS = as.integer(SNP_READ_POS),
          FASTA_REF = stri_sub(SEQUENCES, from = SNP_READ_POS, to = SNP_READ_POS),
          IN_GROUP = stri_sub(SEQUENCES, from = SNP_READ_POS-1, to = SNP_READ_POS+1)
        ) %>%
        # remove lines with no match between all the alleles in the fasta file and the REF in the VCF
        filter(FASTA_REF == REF) %>%
        group_by(CHROM, LOCUS, POS, REF) %>%
        distinct(CHROM, LOCUS, POS, REF) %>% 
        mutate(
          IN_GROUP = ifelse(SNP_READ_POS == max.length.read, stri_pad(IN_GROUP, width = 3, side = "right", pad = "-"), IN_GROUP),
          IN_GROUP = ifelse(SNP_READ_POS == 1, stri_pad(IN_GROUP, width = 3, side = "left", pad = "-"), IN_GROUP)
        )
      
      
      ingroup <- ref.allele.vcf.ingroup.fasta %>% 
        ungroup() %>% 
        select(CHROM, LOCUS, POS, IN_GROUP) %>% 
        mutate(
          POS = stri_pad_left(str = POS, width = 8, pad = "0"),
          LOCUS = stri_pad_left(str = LOCUS, width = 8, pad = "0")
        ) %>%
        arrange(CHROM, LOCUS, POS) %>%
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
      
      markers.id <- ref.allele.vcf.ingroup.fasta %>% 
        ungroup %>% 
        select(CHROM, LOCUS, POS, SNP_READ_POS)
      
      markers.id.ingroup.nuc <- ref.allele.vcf.ingroup.fasta %>% 
        ungroup %>% 
        select(CHROM, LOCUS, POS, IN_GROUP)
      
      # remove unused objects
      ref.allele.vcf.ingroup <- NULL
      ref.allele.vcf.ingroup.fasta <- NULL
      sumstats.ingroup <- NULL
      sumstats.outgroup <- NULL
      
      # Same thing but with outgroup
      ref.allele.outgroup.fasta <- markers.id %>% 
        left_join(
          fasta.prep %>% 
            filter(ANCESTRAL == "outgroup") %>% 
            select (-ALLELE, -ANCESTRAL)
          , by = "LOCUS"
        ) %>%
        mutate(SEQUENCES = stri_replace_na(SEQUENCES, replacement = "NA")) %>%
        filter(SEQUENCES != "NA") %>%
        # get the flanking bases, left and right, of the SNP
        mutate(OUT_GROUP = stri_sub(SEQUENCES, from = SNP_READ_POS-1, to = SNP_READ_POS+1)) %>%
        select(CHROM, LOCUS, POS, OUT_GROUP, SNP_READ_POS) %>%
        group_by(CHROM, LOCUS, POS, OUT_GROUP, SNP_READ_POS) %>%
        tally %>% 
        group_by(CHROM, LOCUS, POS) %>%
        filter(n == max(n))
      
      fasta.prep <- NULL # remove unused object
      
      # The problem is that some markers in the outgroup will have number of REF and ALT alleles equals...
      # Making the ancestral allele call ambiguous (50/50 chance of differences...)
      ambiguous.ancestral.allele <- ref.allele.outgroup.fasta %>%
        ungroup() %>%
        select(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        group_by(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        tally %>% 
        filter(n > 1) %>% 
        select(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        left_join(markers.id.ingroup.nuc, by = c("CHROM", "LOCUS", "POS")) %>% 
        rename(OUT_GROUP = IN_GROUP) %>% 
        ungroup() %>% 
        arrange(CHROM, LOCUS, POS)
      
      markers.id.ingroup.nuc <- NULL # remove unused object
      
      outgroup <- ref.allele.outgroup.fasta %>% 
        select(-n) %>% 
        anti_join(ambiguous.ancestral.allele, by = c("CHROM", "LOCUS", "POS")) %>%
        ungroup() %>% 
        arrange(CHROM, LOCUS, POS) %>% 
        bind_rows(ambiguous.ancestral.allele) %>%
        ungroup() %>% 
        arrange(CHROM, LOCUS, POS) %>%
        mutate(
          OUT_GROUP = ifelse(SNP_READ_POS == max.length.read, stri_pad(OUT_GROUP, width = 3, side = "right", pad = "-"), OUT_GROUP),
          OUT_GROUP = ifelse(SNP_READ_POS == 1, stri_pad(OUT_GROUP, width = 3, side = "left", pad = "-"), OUT_GROUP)
        ) %>% 
        select(CHROM, LOCUS, POS, OUT_GROUP) %>% 
        mutate(
          POS = stri_pad_left(str = POS, width = 8, pad = "0"),
          LOCUS = stri_pad_left(str = LOCUS, width = 8, pad = "0")
        ) %>%
        arrange(CHROM, LOCUS, POS) %>%
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
      
      ambiguous.ancestral.allele <- NULL # remove unused object
      ref.allele.outgroup.fasta <- NULL # remove unused object
      
      # common markers between ingroup and outgroup
      markers.ingroup <- ingroup %>% select(MARKERS)
      markers.outgroup <- outgroup %>% select(MARKERS)
      
      markers.ingroup.outgroup.common <- bind_rows(markers.ingroup, markers.outgroup) %>% 
        group_by(MARKERS) %>%
        tally %>% 
        filter(n == 2) %>%
        arrange(MARKERS) %>%
        select(MARKERS) %>%
        distinct(MARKERS)
      
      markers.ingroup <- NULL
      markers.outgroup <- NULL
      
      message(stri_join("Number of markers common between in- and out- group = ", n_distinct(markers.ingroup.outgroup.common$MARKERS)))
      
      dadi.input <- suppressWarnings(
        input %>%
          group_by(MARKERS, POP_ID, REF, ALT) %>%
          mutate(
            A1 = (2 * PP) + PQ,
            A2 = (2 * QQ) + PQ
          ) %>%
          ungroup() %>%
          select(POP_ID, Allele1 = REF, A1, Allele2 = ALT, A2, MARKERS) %>%
          tidyr::gather(ALLELE_GROUP, COUNT, -c(POP_ID, MARKERS, Allele1, Allele2)) %>%
          tidyr::unite(POP, POP_ID, ALLELE_GROUP, sep = "_") %>%
          group_by(MARKERS, Allele1, Allele2) %>%
          tidyr::spread(data = ., key = POP, value = COUNT) %>% 
          select(Allele1, contains("A1"), Allele2, contains("A2"), MARKERS) %>%
          arrange(MARKERS) %>% 
          ungroup %>% 
          semi_join(markers.ingroup.outgroup.common, by = "MARKERS") %>% 
          inner_join(ingroup, by = "MARKERS") %>% 
          inner_join(outgroup, by = "MARKERS") %>% 
          select(IN_GROUP, OUT_GROUP, Allele1, contains("A1"), Allele2, contains("A2"), MARKERS) %>% 
          arrange(MARKERS)
      )
      
      # remove unused objects
      markers.id <- NULL
      markers.ingroup.outgroup.common <- NULL
      ingroup <- NULL
      outgroup <- NULL
    }
    
    # We need to modify the header to remove "_A1" and "_A2" that were necessary for spreading the info accross columns
    header.dadi <- colnames(dadi.input)
    colnames(dadi.input) <- stri_replace_all_fixed(str = header.dadi, pattern = c("_A1", "_A2"), replacement = "", vectorize_all = FALSE)
    
    if(is.null(dadi.input.filename)){
      # when filename is not provided will save the 'dadi.input' with date and time
      file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
      file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "_", ""), vectorize_all = FALSE)
      if (write.imputation == FALSE) {
        dadi.input.filename <- stri_paste("dadi.input", file.date, "txt", sep = ".")
      }
      
      if (write.imputation == TRUE) {
        dadi.input.filename.imp <- stri_paste("dadi.input.imputed", file.date, "txt", sep = ".")
      }
    } else {
      if (write.imputation == FALSE) {
        dadi.input.filename <- dadi.input.filename
      }
      
      if (write.imputation == TRUE) {
        dadi.input.filename.imp <- stri_replace_all_fixed(dadi.input.filename,
                                                          pattern = "txt",
                                                          replacement = stri_join(
                                                            i, m, 
                                                            "imputed.txt", sep = "_"
                                                          )
        )
      }
    }
    file.version <- suppressWarnings(stri_paste("#\u2202a\u2202i SNP input file generated with stackr v.", packageVersion("stackr"), sep = ""))
    file.date <- suppressWarnings(stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = ""))
    file.header.line <- suppressWarnings(as.data.frame(stri_paste(file.version, file.date, sep = " ")))
    
    if (write.imputation == FALSE) {
      write_tsv(x = file.header.line, path = dadi.input.filename, append = FALSE, col_names = FALSE) # write the header line
      write_tsv(x = dadi.input, path = dadi.input.filename, append = TRUE, col_names = TRUE) # write the data frame
      message(stri_paste("\u2202a\u2202i input file name is: ", dadi.input.filename, "\n", "Saved here: ", getwd()))
    }
    
    if (write.imputation == TRUE) {
      write_tsv(x = file.header.line, path = dadi.input.filename.imp, append = FALSE, col_names = FALSE) # write the header line
      write_tsv(x = dadi.input, path = dadi.input.filename.imp, append = TRUE, col_names = TRUE) # write the data frame
      message(stri_paste("\u2202a\u2202i input file name is: ", dadi.input.filename.imp, "\n", "Saved here: ", getwd()))
    }
    
    
  } # End function write dadi
  
  # without imputations (automatic)
  write_dadi(input = input.count, write.imputation = FALSE)
  
  input.count <- NULL # remove unused object
  
  # Imputations **************************************************************
  if (imputation.method != "FALSE") {
    
    if (imputation.method == "max"){
      message("Calculating map-independent imputations using the most frequent allele.")
    } 
    if (imputation.method == "rf"){
      message("Calculating map-independent imputations using random forest")
    }
    
    # common part for imputation and vcf2genind
    if (impute == "genotype") {
      input.prep <- input %>%
        mutate(GT = stri_replace_all_fixed(str = GT, pattern = "/", replacement = ":", vectorize_all = FALSE)) %>%
        mutate(GT = stri_replace_all_fixed(str = GT, pattern = c("0:0", "1:1", "0:1", "1:0", ".:."), replacement = c("2_0", "0_2", "1_1", "1_1", "NA_NA"), vectorize_all = FALSE)) %>%
        select(-REF, -ALT) %>% 
        arrange(MARKERS, POP_ID) %>% 
        mutate(
          GT = stri_replace_all_fixed(GT, pattern = "NA_NA", replacement = "NA", vectorize_all = FALSE),
          GT = replace(GT, which(GT == "NA"), NA)
        ) %>%
        group_by(INDIVIDUALS, POP_ID) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>%
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
    }
    if (impute == "allele") {
      input.prep <- input %>%
        mutate(GT = stri_replace_all_fixed(str = GT, pattern = "/", replacement = ":", vectorize_all = FALSE)) %>% 
        mutate(GT = stri_replace_all_fixed(str = GT, pattern = c("0:0", "1:1", "0:1", "1:0", ".:."), replacement = c("REF_REF", "ALT_ALT", "REF_ALT", "ALT_REF", "NA_NA"), vectorize_all = FALSE)) %>% 
        select(-REF, -ALT) %>% 
        arrange(MARKERS, POP_ID) %>% 
        tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
        tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% 
        mutate(GT = replace(GT, which(GT == "NA"), NA)) %>%
        group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
    }
    
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
          # imputed.dataset[[i]] <- impute_markers_rf(sep.pop) # test with foreach
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
    
    # transform the imputed dataset  ------------------------
    if (impute == "genotype") {
      # Compute count and Minor Allele Frequency
      message("Computing the Allele Frequency Spectrum for the imputed data")
      
      input.count.imp <- input.imp %>%
        tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% # make tidy
        mutate(GT = stri_replace_all_fixed(str = GT, pattern = c("2_0", "0_2", "1_1"), replacement = c("0/0", "1/1", "0/1"), vectorize_all = FALSE)) %>% 
        group_by(MARKERS, POP_ID) %>%
        summarise(
          N = as.numeric(n()),
          PP = as.numeric(length(GT[GT == "0/0"])),
          PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
          QQ = as.numeric(length(GT[GT == "1/1"]))
        ) %>%
        mutate(MAF = ((QQ*2) + PQ)/(2*N)) %>% 
        arrange(MARKERS, POP_ID) %>% 
        inner_join(
          input %>% 
            select(MARKERS, POP_ID, REF, ALT) %>% 
            distinct(MARKERS, POP_ID, REF, ALT)
          , by = c("MARKERS", "POP_ID")
        ) %>% 
        select(MARKERS, POP_ID, REF, ALT, N, PP, PQ, QQ, MAF)
    }
    
    if (impute == "allele") {
      # Compute count and Minor Allele Frequency
      message("Computing the Allele Frequency Spectrum for the imputed data")
      
      input.count.imp <- input.imp %>%
        tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
        tidyr::spread(data = ., key = ALLELES, value = GT) %>%
        tidyr::unite(GT, A1, A2, sep = "_") %>% 
        mutate(GT = stri_replace_all_fixed(str = GT, pattern = c("REF_REF", "ALT_ALT", "REF_ALT", "ALT_REF"), replacement = c("0/0", "1/1", "0/1", "1/0"), vectorize_all = FALSE)) %>% 
        group_by(MARKERS, POP_ID) %>%
        summarise(
          N = as.numeric(n()),
          PP = as.numeric(length(GT[GT == "0/0"])),
          PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
          QQ = as.numeric(length(GT[GT == "1/1"]))
        ) %>%
        mutate(MAF = ((QQ*2) + PQ)/(2*N)) %>% 
        arrange(MARKERS, POP_ID) %>% 
        inner_join(
          input %>% 
            select(MARKERS, POP_ID, REF, ALT) %>% 
            distinct(MARKERS, POP_ID, REF, ALT)
          , by = c("MARKERS", "POP_ID")
        ) %>% 
        select(MARKERS, POP_ID, REF, ALT, N, PP, PQ, QQ, MAF)
    }
    input.imp <- NULL # remove unused object

    # output in dadi
    write_dadi(input = input.count.imp, write.imputation = TRUE)
    } # end imputations
} # End vcf2dadi


