# Write a betadiv object from STACKS VCF file

# to get rid of notes in build check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE", "ALLELE1", "ALLELE2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX", "other", "strata", "hierarchy", "GROUP", "."))


#' @name vcf2betadiv
#' @title Create a \code{betadiv} object from a \code{STACKS} vcf file
#' @description This function can first filter the vcf file 
#' with a whitelist of loci
#' and a blacklist of individuals (optional). Then it will convert the file
#' to a \code{betadiv} object.
#' Map-independent imputation using Random Forest or the most frequent category
#' is also available as an option.
#' @param vcf.file The VCF file created by STACKS.
#' @param whitelist.markers (optional) A whitelist containing CHROM (character or integer) and/or LOCUS (integer) and/or 
#' POS (integer) columns header. To filter by CHROM and/or locus and/or by snp.
#' The whitelist is in the directory (e.g. "whitelist.txt"). de novo CHROM column with 'un' need to be changed to 1.
#' @param snp.LD (optional) Minimize linkage disequilibrium (LD) by choosing 
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
#' @param imputations Should a map-independent imputations of markers be
#' computed. Available choices are: (1) \code{FALSE} for no imputation.
#' (2) \code{"max"} to use the most frequent category for imputations.
#' (3) \code{"rf"} using Random Forest algorithm. Default = \code{FALSE}.
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
#' @return When no imputation is selected a betadiv object is returned.
#' When imputation is selected a list with 2 objects is returned
#' and accessed with \code{$no.imputation} or \code{$imputed}.
#' @export
#' @rdname vcf2betadiv
#' @import reshape2
#' @import dplyr
#' @references Lamy T, Legendre P, Chancerelle Y, Siu G, Claudet J (2015) 
#' Understanding the Spatio-Temporal Response of Coral Reef Fish Communities to 
#' Natural Disturbances: Insights from Beta-Diversity Decomposition. 
#' PLoS ONE, 10, e0138696.
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
#' @seealso \code{beta.div} is available on Pierre Legendre web site \url{http://adn.biol.umontreal.ca/~numericalecology/Rcode/} \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Laura Benestan \email{laura.benestan@@icloud.com} and
#' Thierry Gosselin \email{thierrygosselin@@icloud.com}

vcf2betadiv <- function(vcf.file, 
                       whitelist.markers = NULL,
                       snp.LD = NULL,
                       common.markers = FALSE,
                       blacklist.id = NULL, 
                       pop.id.start, pop.id.end,
                       pop.levels,
                       pop.labels,
                       imputations = FALSE,
                       imputations.group = "populations",
                       num.tree = 100,
                       iteration.rf = 10,
                       split.number = 100,
                       verbose = FALSE,
                       parallel.core = 2
) {
  
  # remove NOTE about no visible binding for global variable during Build 
  QUAL <- NULL
  FILTER <- NULL
  FORMAT <- NULL
  FORMAT_ID <- NULL
  ID <- NULL
  '#CHROM' <- NULL
  INFO <- NULL
  REF <- NULL
  ALT <- NULL
  READ_DEPTH <- NULL
  ALLELE_DEPTH <- NULL
  GT <- NULL
  GL <- NULL
  MARKERS <- NULL
  MARKERS_ALLELES <- NULL
  ALLELES <- NULL
  COUNT <- NULL
  i <- NULL
  QQ <- NULL
  PQ <- NULL
  MAF <- NULL
  
  
  if (imputations == "FALSE") {
    message("vcf2betadiv: without imputation...")
  } else {
    message("vcf2betadiv: with imputations...")
  }
  
  # Import/read VCF ------------------------------------------------------------- 
  message("Importing the VCF...")
  
  vcf <- read_delim(
    vcf.file, 
    delim = "\t", 
    comment = "##",
    progress = interactive()
  ) %>% 
    select(-c(QUAL, FILTER, INFO, REF, ALT)) %>% 
    rename(LOCUS = ID, CHROM = `#CHROM`) %>% 
    mutate(CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1"))
  
  
  # Detect STACKS version used to create the vcf
  if(stri_detect_fixed(vcf$FORMAT[1], "AD")) {
    stacks.version <- "new"
  } else{
    stacks.version <- "old"
  }
  vcf <- vcf %>% select(-FORMAT)
  
  # Whitelist of markers ------------------------------------------------------------
  
  if (is.null(whitelist.markers) | missing(whitelist.markers)) { # no Whitelist
    message("No whitelist to apply to the VCF")
    whitelist.markers <- NULL
    vcf <- vcf
  } else { # with Whitelist of markers
    message("Filtering the VCF with the whitelist from your directory")
    whitelist.markers <- read_tsv(whitelist.markers, col_names = TRUE)
    columns.names.whitelist <- colnames(whitelist.markers)
    if("CHROM" %in% columns.names.whitelist){
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
    vcf <- vcf %>% semi_join(whitelist.markers, by = columns.names.whitelist)
  }
  
  # Tidying the VCF to make it easy to work on the data for conversion----------
  # Preping the pop.labels
  if(missing(pop.labels)){
    pop.labels <- pop.levels
  } else {
    pop.labels <- pop.labels
  }
  
  message("Making the VCF population wise")
  vcf <- suppressWarnings(
    vcf %>%
      tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, LOCUS, POS)) %>% # Gather individuals in 1 colummn
      mutate( # Make population ready
        POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
        POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = F), levels = pop.labels, ordered =T),
        POP_ID = droplevels(POP_ID),
        INDIVIDUALS =  as.character(INDIVIDUALS)
      )
  )
  
  # Blacklist id -----------------------------------------------------------------
  if (is.null(blacklist.id) | missing(blacklist.id)) { # No blacklist of ID
    message("No individual blacklisted")
    blacklist.id <- NULL
    vcf <- vcf
  } else { # With blacklist of ID
    message("Using the blacklisted id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)
    vcf <- suppressWarnings(
      vcf %>% 
        anti_join(blacklist.id, by = "INDIVIDUALS") %>% 
        mutate(POP_ID = droplevels(POP_ID))
    )
  }
  
  # dump unused object
  blacklist.id <- NULL
  whitelist.markers <- NULL
  
  # Conversion into betadiv -----------------------------------------------------
  message("Tidy vcf into factory for conversion into betadiv ...")
  
  if(stacks.version == "new"){ # with new version of stacks > v.1.29
    vcf <- suppressWarnings(
      vcf %>%
        tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"), # no imputation
                        sep = ":", extra = "warn") %>%  # no imputation
        select(-c(READ_DEPTH, ALLELE_DEPTH, GL)) # no imputation
    )
  } else { # stacks version prior to v.1.29 had no Allele Depth field...
    vcf <- suppressWarnings(
      vcf %>%
        tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "GL"), # no imputation
                        sep = ":", extra = "warn") %>%  # no imputation
        select(-c(READ_DEPTH, GL))  # no imputation
    )
  }
  message("step 1/2: completed")
  
  # Keep only 1 SNP per haplotypes/reads
  if(missing(snp.LD) | is.null(snp.LD)){
    vcf <- vcf
  } else{
    snp.locus <- vcf %>% distinct(LOCUS, POS)
    # Random selection
    if(snp.LD == "random"){
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>% 
        sample_n(size = 1, replace = FALSE)
      message(stri_c("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP randomly selected to keep 1 SNP per read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # Fist SNP on the read
    if(snp.LD == "first"){
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>% 
        summarise(POS = min(POS))
      message(stri_c("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # Last SNP on the read
    if(snp.LD == "last"){
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>% 
        summarise(POS = max(POS))
      message(stri_c("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # filtering the VCF to minimize LD
    vcf <- vcf %>% semi_join(snp.select, by = c("LOCUS", "POS"))
    message("Filtering the tidy VCF to minimize LD by keeping only 1 SNP per short read/haplotype")
  }
  
  # combine CHROM, LOCUS and POS into MARKERS
  vcf <- vcf %>%
    mutate(
      POS = stri_pad_left(str = POS, width = 7, pad = "0"),
      LOCUS = stri_pad_left(str = LOCUS, width = 7, pad = "0")
    ) %>%
    arrange(CHROM, LOCUS, POS) %>% 
    tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
  
  # keeping or not markers in common in all the populations
  if(missing(common.markers) | is.null(common.markers) | common.markers == FALSE) {
    vcf <- vcf
  }
  
  if(common.markers == TRUE){ # keep only markers present in all pop
    pop.number <- n_distinct(vcf$POP_ID)
    
    pop.filter <- vcf %>%
      filter(GT != "./.") %>%
      group_by(MARKERS) %>%
      filter(n_distinct(POP_ID) == pop.number) %>%
      arrange(MARKERS) %>% 
      distinct(MARKERS)
    
    message(stri_c("Number of original markers = ", n_distinct(vcf$MARKERS), "\n", "Number of markers present in all the populations = ", n_distinct(pop.filter$MARKERS), "\n", "Number of markers removed = ", n_distinct(vcf$MARKERS) - n_distinct(pop.filter$MARKERS)))
    vcf <- vcf %>% semi_join(pop.filter, by = "MARKERS")
  }
  
  
  message("step 2/2: completed")
  
  # results no imputation--------------------------------------------------------------------
  res <- list()
  res$no.imputation <- vcf %>% 
    group_by(MARKERS, POP_ID) %>%
    summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT[GT == "0/0"])),
      PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
      QQ = as.numeric(length(GT[GT == "1/1"]))
    ) %>%
    mutate(MAF = ((QQ*2) + PQ)/(2*N)) %>%
    select(POP_ID, MARKERS, MAF) %>% 
    group_by(POP_ID) %>% 
    tidyr::spread(key = MARKERS, value = MAF) %>% 
    mutate(POP_ID = as.integer(POP_ID))
  
  
  if (imputations == "FALSE") {
    message("A large 'betadiv' object (no imputation) was created in your Environment")
  } else if (imputations == "max"){
    message("Calculating map-independent imputations using the most frequent allele.")
  } else {
    message("Calculating map-independent imputations using random forest")
  }
  
    # Imputations: betadiv with imputed haplotypes using Random Forest ------------------
  if (imputations != "FALSE"){
    
    vcf.prep <- vcf %>%
      mutate(GT = replace(GT, which(GT == "./."), NA)) %>% 
      dcast(INDIVIDUALS + POP_ID ~ MARKERS, value.var = "GT") %>% 
      arrange(POP_ID, INDIVIDUALS)
    
    vcf <- NULL # remove unused object
    
    if (imputations == "rf") {
      # Parallel computations options
      if (missing(parallel.core) == "TRUE"){
        # Automatically select all the core -1 
        options(rf.cores=detectCores()-1, mc.cores=detectCores()-1)
        # Start cluster registration backend using n - 1 CPU
        cl <- makeCluster(detectCores() - 1)
        registerDoParallel(cl, cores = detectCores() - 1)
      } else {
        options(rf.cores = parallel.core, mc.cores = parallel.core)
        # Start cluster registration backend using n - 1 CPU
        cl <- makeCluster(parallel.core)
        registerDoParallel(cl, cores = parallel.core)
      }
      
      # imputations using Random Forest with the package randomForestSRC
      
      impute_markers_rf <- function(x){
        randomForestSRC::impute.rfsrc(data = x, 
                                      ntree = num.tree, 
                                      nodesize = 1, 
                                      nsplit = split.number, 
                                      nimpute = iteration.rf, 
                                      do.trace = verbose)
      }
      
      # imputations by populations (default) or globally -------------------------
      
      # default by pop
      if (missing(imputations.group) == "TRUE" | imputations.group == "populations"){
        message("Imputations computed by populations, take a break...")
        df.split.pop <- split(x = vcf.prep, f = vcf.prep$POP_ID) # slip data frame by population
        pop.list <- names(df.split.pop) # list the pop
        imputed.dataset <-list() # create empty list 
        imputed.dataset <- foreach(i=pop.list, .packages = c("magrittr", "plyr", "dplyr", "tidyr", "stringi", "readr", "randomForestSRC", "reshape2")) %dopar% {
          sep.pop <- df.split.pop[[i]]
          sep.pop <- suppressWarnings(
            plyr::colwise(factor, exclude = NA)(sep.pop)
          )
          # message of progress for imputations by population
          message(paste("Completed imputations for pop ", i, sep = ""))
          imputed.dataset[[i]] <- impute_markers_rf(sep.pop)
        }
        # close parallel connection settings
        stopCluster(cl)
        message("Almost finished with the imputations...")
        vcf.imp <- suppressWarnings(as.data.frame(bind_rows(imputed.dataset)))
        
        # Second round of imputations: remove introduced NA if some pop don't have the markers by using
        # RF globally
        vcf.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(vcf.imp)) # Make the columns factor
        vcf.imp <- impute_markers_rf(vcf.imp) # impute globally
        
        # dump unused objects
        df.split.pop <- NULL
        pop.list <- NULL
        sep.pop <- NULL
        imputed.dataset <- NULL
        vcf.prep <- NULL
        
      } else if (imputations.group == "global"){
        # Globally (not by pop_id)
        message("Imputations computed globally, take a break...")
        vcf.prep <- plyr::colwise(factor, exclude = NA)(vcf.prep)
        vcf.imp <- impute_markers_rf(vcf.prep)
        
        vcf.prep <- NULL # remove unused object
        
      } 
      
    } else if (imputations == "max") {
      
      if (missing(imputations.group) == "TRUE" | imputations.group == "populations"){
        message("Imputations computed by populations")
        
        vcf.imp <- suppressWarnings(
          vcf.prep %>%
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
            dcast(INDIVIDUALS + POP_ID ~ MARKERS, value.var = "GT")
        )
        
        vcf.prep <- NULL # remove unused object
        
      } else if (imputations.group == "global"){
        # Globally (not by pop_id)
        message("Imputations computed globally")
        
        vcf.imp <- suppressWarnings(
          vcf.prep %>%
            tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
            group_by(MARKERS) %>% 
            mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>% 
            dcast(INDIVIDUALS + POP_ID ~ MARKERS, value.var = "GT")
        )
        
        vcf.prep <- NULL # remove unused object
      }
    }
    
    # transform the imputed dataset into betadiv object ------------------------
    
    message("Imputed vcf into factory for conversion into betadiv...")
    
    res$imputed <- suppressWarnings(
      vcf.imp %>%
        tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% # make tidy
        group_by(MARKERS, POP_ID) %>%
        summarise(
          N = as.numeric(n()),
          PP = as.numeric(length(GT[GT == "0/0"])),
          PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
          QQ = as.numeric(length(GT[GT == "1/1"]))
        ) %>%
        mutate(MAF = ((QQ*2) + PQ)/(2*N)) %>%
        select(POP_ID, MARKERS, MAF) %>% 
        group_by(POP_ID) %>% 
        tidyr::spread(key = MARKERS, value = MAF) %>% 
        mutate(POP_ID = as.integer(POP_ID))
    )
    vcf.imp <- NULL # remove unused object
    message("Imputed betadiv data set saved into the Global Environment")
    }
  # outout results -------------------------------------------------------------
  return(res)
}
