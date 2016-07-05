# VCF data imputation using Random Forest

#' @name vcf2genepop
#' @title VCF to genepop with filters and data imputation
#' @description This function can first filter the vcf file 
#' with a whitelist of loci and a blacklist of individuals (optional). 
#' Then it will convert the file
#' to a genepop file.
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
#' @param hierarchy (optional) A formula that explicitely defines hierarchical levels 
#' in your strata. See \code{\link[adegenet]{genind}} for details.


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

#' @param genepop.filename The name of the file written to the directory.
#' Use the extension ".gen" at the end. 
#' Default: \code{genepop.filename = "genepop.gen"}.
#' @param genepop.header The first line of the Genepop file.
#' Default: \code{genepop.header = "my firt genepop"}.


#' @details The imputations using Random Forest requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set 
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals 
#' will require 15 min.

#' @return When no imputation is selected a genepop file is saved to the 
#' working directory. When imputation is selected 2 genepop files are saved to
#' the working directory.

#' @export
#' @rdname vcf2genepop
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
#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1: 
#' new tools for the analysis of genome-wide SNP data. 
#' Bioinformatics, 27, 3070-3071.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.
#' @seealso \code{adegenet} is available on CRAN \url{http://cran.r-project.org/web/packages/adegenet/} and github \url{https://github.com/thibautjombart/}
#' \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}
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


vcf2genepop <- function(data, 
                        whitelist.markers = NULL,
                        blacklist.genotype = NULL,
                        snp.ld = NULL,
                        common.markers = FALSE,
                        blacklist.id = NULL,
                        pop.id.start, pop.id.end,
                        pop.levels,
                        pop.labels,
                        strata = NULL,
                        pop.select = NULL,
                        hierarchy = NULL,
                        imputation.method = NULL,
                        impute = "genotype",
                        imputations.group = "populations",
                        num.tree = 100,
                        iteration.rf = 10,
                        split.number = 100,
                        verbose = FALSE,
                        parallel.core = detectCores()-1,
                        genepop.filename = "genepop.gen",
                        genepop.header = "my firt genepop"
) {
  if (missing(data)) stop("Input file missing")
  if (missing(pop.id.start)) pop.id.start <- NULL
  if (missing(pop.id.end)) pop.id.end <- NULL
  if (missing(pop.levels)) stop("pop.levels required")
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels

  if (is.null(imputation.method)) {
    message("vcf2genepop: without imputation...")
  } else {
    message("vcf2genepop: with/without imputations...")
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
  
  # Conversion into genepop -----------------------------------------------------
  # Tidy VCF
  message("Tidy vcf into genepop factory...")
  
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
      whitelist.markers.ind <- input %>% distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      
      
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
    snp.locus <- input %>% distinct(LOCUS, POS)
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
  message("Recoding genotypes for genepop")
  input <- input %>%
    mutate(
      REF= stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE), # replace nucleotide with numbers
      ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE),# replace nucleotide with numbers
      REF = stri_pad_left(str = REF, width = 3, pad = "0"),
      ALT = stri_pad_left(str = ALT, width = 3, pad = "0"),
      GT = ifelse(GT == "0/0", stri_join(REF, REF, sep = "_"),
                  ifelse(GT == "1/1",  stri_join(ALT, ALT, sep = "_"),
                         ifelse(GT == "0/1", stri_join(REF, ALT, sep = "_"),
                                ifelse(GT == "1/0", stri_join(ALT, REF, sep = "_"), "000_000")
                         )
                  )
      )
    ) %>%
    arrange(MARKERS, POP_ID) %>%
    select(-c(REF, ALT))
  
  loci <- input %>% distinct(MARKERS) %>% arrange(MARKERS)
  
  # results no imputation--------------------------------------------------------------------
  genepop.prep <- input %>%
    mutate(
      GT= stri_replace_all_fixed(str = GT, pattern = "_", replacement = "", vectorize_all = FALSE)
    ) %>% 
    arrange(MARKERS) %>% 
    select(POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
    group_by(POP_ID, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
    arrange(POP_ID, INDIVIDUALS) %>%
    ungroup() %>% 
    mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))
  
  # Create a vector with the population ordered by levels
  pop <- genepop.prep$POP_ID
  # levels(pop)
  
  # split genepop by populations
  genepop.prep <- select(.data = genepop.prep, -POP_ID)
  genepop.split <- split(genepop.prep, pop)
  
  # Write the file in genepop format
  message("Output: No imputation")
  
  genepop.header <- as.data.frame(genepop.header)
  write.table(genepop.header, file = genepop.filename, col.names = FALSE, row.names = FALSE, quote = FALSE)
  write_delim(x = loci, path = genepop.filename, delim = "\n", append = TRUE, col_names = FALSE)
  pop.sep <- as.data.frame("pop")
  for (i in 1:length(genepop.split)) {
    write_delim(x = pop.sep, path = genepop.filename, delim = "\n", append = TRUE, col_names = FALSE)
    write_delim(x = genepop.split[[i]], path = genepop.filename, delim = " ", append = TRUE, col_names = FALSE)
  }
  message(paste0("Genepop file: ", genepop.filename))
  message(paste0("Written to the working directory: \n", getwd()))
  
  genepop.prep <- NULL
  genepop.split <- NULL
  
  # Imputations: genind with imputed haplotypes using Random Forest*************
  if (!is.null(imputation.method)) {
    
    if (impute == "genotype") {
      input.prep <- input %>%
        mutate(
          GT = stri_replace_all_fixed(GT, pattern = "000_000", replacement = "NA", vectorize_all = FALSE),
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
          GT = stri_replace_all_fixed(GT, pattern = "000", replacement = "NA", vectorize_all = FALSE),
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
    message("Imputed vcf into genepop factory...")
    if (impute == "genotype") {
      genepop.prep.imp <- suppressWarnings(
        input.imp %>%
          tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% # make tidy
          mutate(
            GT= stri_replace_all_fixed(str = GT, pattern = "_", replacement = "", vectorize_all = FALSE)
          ) %>%
          arrange(MARKERS, POP_ID, INDIVIDUALS) %>%
          group_by(POP_ID, INDIVIDUALS) %>%
          tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
          arrange(POP_ID, INDIVIDUALS)
      )
    }
    
    if (impute == "allele") {
      genepop.prep.imp <- suppressWarnings(
        input.imp %>%
          tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
          mutate(
            GT= stri_replace_all_fixed(str = GT, pattern = "_", replacement = "", vectorize_all = FALSE)
          ) %>% 
          tidyr::spread(data = ., key = ALLELES, value = GT) %>%
          tidyr::unite(GT, A1, A2, sep = "") %>%
          arrange(MARKERS, POP_ID, INDIVIDUALS) %>% 
          group_by(POP_ID, INDIVIDUALS) %>%
          tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
          arrange(POP_ID, INDIVIDUALS)
      )
    }
    input.imp <- NULL # remove unused object
    
    genepop.prep.imp <- genepop.prep.imp %>%
      arrange(POP_ID, INDIVIDUALS) %>%
      ungroup() %>% 
      mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = "")) %>% 
      select(-POP_ID)
    
    # split genepop by populations
    genepop.imp.split <- split(genepop.prep.imp, pop)
    
    # Write the file in genepop format 
    message("Output: with imputations")
    
    # Add "_imputed" to the filename
    genepop.filename.imp <- stri_replace_all_fixed(genepop.filename,
                                                   pattern = ".gen",
                                                   replacement = "_imputed.gen")
    
    
    genepop.header <- genepop.header$genepop.header
    genepop.header <- stri_join(genepop.header, "with imputation", sep = " ")
    genepop.header <- as.data.frame(genepop.header)
    write.table(genepop.header, file = genepop.filename.imp, col.names = FALSE, row.names = FALSE, quote = FALSE)
    loci.table <- as.data.frame(loci)
    write_delim(x = loci.table, path = genepop.filename.imp, delim = "\n", append = TRUE, col_names = FALSE)
    pop.sep <- as.data.frame("pop")
    for (i in 1:length(genepop.imp.split)) {
      write_delim(x = pop.sep, path = genepop.filename.imp, delim = "\n", append = TRUE, col_names = FALSE)
      write_delim(x = genepop.imp.split[[i]], path = genepop.filename.imp, delim = " ", append = TRUE, col_names = FALSE)
    }    
    
    
    message(paste0("Genepop file: ", genepop.filename.imp))
    message(paste0("Written to the working directory: \n", getwd()))
    
    genepop.prep <- NULL
    genepop.split <- NULL
  } # End imputation
}
