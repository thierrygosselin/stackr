# Write a hierfstat object from STACKS VCF file

# to get rid of notes in build check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE", "ALLELE1", "ALLELE2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX", "other", "strata", "hierarchy", "GROUP", "."))


#' @name vcf2hierfstat
#' @title Create a \code{hierfstat} object from a \code{STACKS} vcf file
#' @description This function can first filter the vcf file 
#' with a whitelist of loci
#' and a blacklist of individuals (optional). Then it will convert the file
#' to a \code{hierfstat} object.
#' Map-independent imputation using Random Forest or the most frequent category
#' is also available as an option.
#' @param vcf.file The VCF file created by STACKS.
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
#' @param fstat.filename The name of the file written to the directory.
#' Use the extension '.dat' at the end. Default \code{fstat_gbs.dat}.
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
#' @return When no imputation is selected an hierfstat object is returned.
#' When imputation is selected a list with 2 hierfstat dataframe objects is returned
#' and accessed with \code{$no.imputation} or \code{$imputed}.
#' @export
#' @rdname vcf2hierfstat
#' @import reshape2
#' @import dplyr
#' @importFrom stringr str_pad
#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Goudet, J. (1995) FSTAT (Version 1.2): A computer program to 
#' calculate F- statistics. Journal of Heredity, 86, 485-486.
#' @references Goudet, J. (2005) hierfstat, a package for r to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5, 184-186.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.
#' @seealso \code{hierfstat} is available on CRAN \url{http://cran.r-project.org/web/packages/hierfstat/} and github \url{https://github.com/jgx65/hierfstat/}
#' \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

vcf2hierfstat <- function(vcf.file, 
                          whitelist.markers = NULL, 
                          blacklist.id = NULL, 
                          pop.id.start, 
                          pop.id.end,
                          pop.levels,
                          pop.labels,
                          fstat.filename = "fstat_gbs.dat",
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
  
  
  if (imputations == "FALSE") {
    message("vcf2hierfstat: without imputation...")
  } else {
    message("vcf2hierfstat: with imputations...")
  }
  
  # Import/read VCF ------------------------------------------------------------- 
  message("Importing the VCF...")
  
  vcf <- read_delim(
    vcf.file, 
    delim = "\t", 
    comment = "##",
    progress = interactive()
  ) %>% 
    # select(-c(QUAL, FILTER, INFO, REF, ALT)) %>% 
    select(-c(QUAL, FILTER, INFO)) %>% 
    rename(LOCUS = ID, CHROM = `#CHROM`) %>% 
    mutate(CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1"))
  
  # Detect STACKS version
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
      tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_") %>% # group markers info
      tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(MARKERS, REF, ALT)) %>% # Gather individuals in 1 colummn
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
  
  # Get the number of sample (pop) for hierfstat
  np <- nlevels(droplevels(vcf$POP_ID))
  np.message <- stri_paste("Number of sample pop, np = ", np, sep = "")
  message(np.message)
  
  # get the list of loci after filter
  loci <- unique(vcf$MARKERS)
  
  # Get the number of loci
  nl <- length(loci)
  nl.message <- stri_paste("Number of markers, nl = ", nl, sep = "")
  message(nl.message)
  
  # dump unused object
  blacklist.id <- NULL
  whitelist.markers <- NULL
  
  # Conversion into hierfstat -----------------------------------------------------
  message("Tidy vcf into factory for conversion into hierfstat ...")
  
  if(stacks.version == "new"){ # with new version of stacks > v.1.29
    vcf <- vcf %>%
      tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"), # no imputation
                      sep = ":", extra = "warn") %>%  # no imputation
      select(-c(READ_DEPTH, ALLELE_DEPTH, GL)) # no imputation
  } else { # stacks version prior to v.1.29 had no Allele Depth field...
    vcf <- vcf %>%
      tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "GL"), # no imputation
                      sep = ":", extra = "warn") %>%  # no imputation
      select(-c(READ_DEPTH, GL))  # no imputation
  }
  message("step 1/3: completed")
  
  # testing
  # vcf.bk <- vcf
  # vcf <- vcf.bk
  
  # Get the highest number used to label an allele
  nu <- 4
  
  # allele coding
  allele.coding <- 1
  message("The alleles are encoded with one digit number")
  
  # Change the genotype coding for easier integration in downstream conversion to hierfstat
  vcf <- vcf %>% 
    mutate(
      REF= stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE), # replace nucleotide with numbers
      ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE),# replace nucleotide with numbers
      GT = ifelse(GT == "0/0", stri_c(REF, REF, sep = "_"),
                  ifelse(GT == "1/1",  stri_c(ALT, ALT, sep = "_"),
                         ifelse(GT == "0/1", stri_c(REF, ALT, sep = "_"),
                                ifelse(GT == "1/0", stri_c(ALT, REF, sep = ""), "0_0")
                         )
                  )
      )
    ) %>% 
    arrange(MARKERS, POP_ID) %>% 
    select(-c(REF, ALT))
  
  message("step 2/3: completed")
  
  # results no imputation--------------------------------------------------------------------
  # convert to hierfstat
  
  hierfstat.prep <- vcf %>%
    mutate(
      GT = stri_replace_all_fixed(str = GT, pattern = "_", replacement = "", vectorize_all = FALSE),
      GT = stri_replace_all_fixed(GT, pattern = "00", replacement = "NA", vectorize_all = FALSE),
      GT = replace(GT, which(GT == "NA"), NA),
      GT = as.integer(GT)
    ) %>% 
    dcast(INDIVIDUALS + POP_ID ~ MARKERS, value.var = "GT") %>% 
    arrange(POP_ID, INDIVIDUALS)
  
  rownames(hierfstat.prep) <- hierfstat.prep$INDIVIDUALS
  
  hierfstat.prep <- select(.data = hierfstat.prep, -INDIVIDUALS)
  
  res <- hierfstat.prep
  message("step 3/3: completed")
  hierfstat.prep$POP_ID <- as.integer(hierfstat.prep$POP_ID) # Change pop id to integer
  
  # fstat filename
  if (missing(fstat.filename) == "TRUE") {
    fstat.filename <- "fstat_gbs.dat"
  } else {
    fstat.filename <- fstat.filename
  }
  
  # FSTAT: write the first line
  fstat.first.line <- stri_paste(np, nl, nu, allele.coding, sep = " ")
  fstat.first.line <- as.data.frame(fstat.first.line)
  write_delim(x = fstat.first.line, path = fstat.filename, delim = "\n", append = FALSE, 
              col_names = FALSE)
  
  # FSTAT: write the locus name to the file
  loci.table <- as.data.frame(loci)
  write_delim(x = loci.table, path = fstat.filename, delim = "\n", append = TRUE, 
              col_names = FALSE)
  
  # FSTAT: write the pop and genotypes
  write_delim(x = hierfstat.prep, path = fstat.filename, delim = "\t", append = TRUE, 
              col_names = FALSE)
  savind.message <- stri_c("The hierfstat file (no imputation) was saved in your working directory: \n", getwd())
  message(savind.message)
  
  if (imputations == "max"){
    message("Calculating map-independent imputations using the most frequent allele.")
  } else {
    message("Calculating map-independent imputations using random forest")
  }
  
  # dump unused objects
  # hierfstat.prep <- NULL
  message("A hierfstat dataframe, without imputation, was created in your Environment")
  
  # Imputations: hierfstat with imputed haplotypes using Random Forest ------------------
  if (imputations != "FALSE"){
    
    vcf.prep <- vcf %>%
      mutate(
        # GT = stri_replace_all_fixed(GT, pattern = "000/000", replacement = "NA", vectorize_all = FALSE),
        GT = stri_replace_all_fixed(GT, pattern = "0_0", replacement = "NA", vectorize_all = FALSE),
        GT = replace(GT, which(GT == "NA"), NA)
      ) %>% 
      dcast(INDIVIDUALS + POP_ID ~ MARKERS, value.var = "GT") %>% 
      arrange(POP_ID, INDIVIDUALS)
    
    vcf <- NULL # remove unused object
    
    if (imputations == "rf") {
      # Parallel computations options
      if (missing(parallel.core) == "TRUE"){
        # Automatically select all the core -1 
        options(rf.cores=detectCores()-1, mc.cores=detectCores()-1)
      } else {
        options(rf.cores = parallel.core, mc.cores = parallel.core)
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
        for (i in pop.list) {
          sep.pop <- df.split.pop[[i]]
          sep.pop <- suppressWarnings(
            plyr::colwise(factor, exclude = NA)(sep.pop)
          )
          imputed.dataset[[i]] <- impute_markers_rf(sep.pop)
          # message of progress for imputations by population
          pop.imputed <- paste("Completed imputations for pop ", i, sep = "")
          message(pop.imputed)
        }
        vcf.imp <- suppressWarnings(as.data.frame(bind_rows(imputed.dataset)))
        
        # Second round of imputations: remove introduced NA if some pop don't have the markers by using
        # RF globally
        vcf.imp <- impute_markers_rf(vcf.imp)
        
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
    
    # transform the imputed dataset into hierfstat object ------------------------
    
    message("Imputed haplotypes into factory for conversion into hierfstat...")
    hierfstat.prep.imp <- suppressWarnings(
      vcf.imp %>%
        tidyr::gather(key = MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>% # make tidy
        mutate(
          GT = stri_replace_all_fixed(str = GT, pattern = "_", replacement = "", vectorize_all = FALSE),
          GT = as.integer(GT)
        ) %>%
        dcast(INDIVIDUALS + POP_ID ~ MARKERS, value.var = "GT") %>% 
        arrange(POP_ID, INDIVIDUALS)
    )
    
    rownames(hierfstat.prep.imp) <- hierfstat.prep.imp$INDIVIDUALS
    hierfstat.prep.imp <- select(.data = hierfstat.prep.imp, -INDIVIDUALS)
    
    no.imputation <- res
    res <- list()
    res$no.imputation <- no.imputation
    no.imputation <- NULL # drop unused object
    res$hierfstat.imputed <- hierfstat.prep.imp
    
    hierfstat.prep.imp$POP_ID <- as.integer(hierfstat.prep.imp$POP_ID) # Change pop id to integer

    vcf.imp <- NULL # remove unused object
    
    # results ------------------------------------------------------------------
    # fstat filename
    if (missing(fstat.filename) == "TRUE") {
      fstat.filename <- "fstat_gbs.dat"
    } else {
      fstat.filename <- fstat.filename
    }
    # Add '_imputed' to the filename
    fstat.filename.imp <- stri_replace_all_fixed(fstat.filename, pattern = ".dat", 
                                                 replacement = "_imputed.dat")
    
    # FSTAT: write the first line
    fstat.first.line <- stri_paste(np, nl, nu, allele.coding, sep = " ")
    fstat.first.line <- as.data.frame(fstat.first.line)
    write_delim(x = fstat.first.line, path = fstat.filename.imp, delim = "\n", 
                append = FALSE, col_names = FALSE)
    
    # FSTAT: write the locus name to the file
    loci.table <- as.data.frame(loci)
    write_delim(x = loci.table, path = fstat.filename.imp, delim = "\n", append = TRUE, 
                col_names = FALSE)
    
    # FSTAT: write the pop and genotypes
    write_delim(x = hierfstat.prep.imp, path = fstat.filename.imp, delim = "\t", append = TRUE, 
                col_names = FALSE)
    
    message("A large list with hierfstat dataframe was created in your Environment (with and without imputations)")
  }
  return(res)
}

