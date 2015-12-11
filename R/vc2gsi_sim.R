# Write a gsi_sim file from STACKS VCF file
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID",
                                                        "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt",
                                                        "HAPLOTYPES", "SAMPLES", "ALLELE", "ALLELE1", "ALLELE2", "GENOTYPE",
                                                        "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX",
                                                        "colwise", "detectCores", "mc.cores", "."))

#' @name vcf2gsi_sim
#' @title Use the STACKS VCF file to run assignment analysis in gsi_sim.
#' @description \code{gsi_sim} is a tool for doing and simulating genetic stock identification. 
#' The \code{vcf2gsi_sim} function can first filter the VCF file 
#' with a whitelist of markers and a blacklist of individuals. 
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
#' @param marker.number (Integer or string of number or "all") Calculations with 
#' fixed or subsample of your markers. Default= \code{"all"}.
#' e.g. To test 500, 1000, 2000 and all  the markers: 
#' \code{marker.number = c(500, 1000, 2000, "all"}. 
#' To use only 500 makers \code{marker.number = 500}.
#' @param sampling.method (character) Should the markers be randomly selected 
#' \code{"random"} for a classic Leave-One-Out (LOO) assignment or 
#' chosen based on ranked Fst \code{"ranked"}, used in a 
#' Traing-Holdout-Leave One Out (THL) assignment ?
#' @param THL (integer, proportion) For the \code{ranked} sampling method. 
#' Default \code{1}: 1 individual is used as holdout. This individual is not 
#' participating in the markers ranking. For each marker marker number, 
#' the analysis will be repeated with all the indiviuals in the data set 
#' (e.g. 500 individuals, 500 times 500, 1000, 2000 markers). 
#' If a proportion is used e.g. \code{0.15},= 15% of individuals in each 
#' populations are chosen randomly as holdout individuals. 
#' You can create different holdout individuals lists with the \code{iterations} 
#' argument.
#' @param iterations With random marker selection the iterations argument = 
#' the number of iterations to repeat marker resampling, default is \code{10} 
#' With \code{marker.number = c(500, 1000)} and default iterations setting, 
#' 500 markers will be randomly chosen 10 times and 1000 markers will be randomly
#' chosen 10 times. For the ranked method, using \code{THL = 1}, the analysis 
#' will be repeated for each individuals in the data set for every 
#' \code{marker.number} selected. With a proportion argument \code{THL = 0.15},
#' 15% of individuals in each populations are chosen randomly as holdout 
#' individuals and this process is reapeated the number of times chosen by the 
#' \code{iterations} value.
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param subsample (Integer or Proportion) Default = \code{subsample = NULL}. 
#' With a proportion argument \code{subsample = 0.15}, 15 percent of individuals
#' in each populations are chosen randomly to represent the dataset. 
#' With \code{subsample = 36}, 36 individuals in each populations are chosen 
#' randomly to represent the dataset.

#' @param gsi_sim.filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. Default \code{gsi_sim_data.txt}.
#' The number of markers used will be appended to the name of the file.
#' @param keep.gsi.files (Boolean) Default \code{FALSE} The input and output gsi_sim files
#' will be deleted from the directory when finished processing. 
#' With \code{TRUE}, remember to allocate a large chunk of the disk space for the analysis.
#' @param pop.levels (required) A character string with your populations ordered.
#' @param pop.labels (optional) A character string of your populations labels.
#' If you need to rename sampling sites in \code{pop.levels} or combined sites/pop
#' into a different names, here is the place.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.

#' @param baseline (optional) A character string with your baseline id. 
#' From the \code{pop.id.start} and \code{pop.id.end} you isolate 
#' the baseline and mixture group. Here you need to give the id for 
#' baseline e.g. \code{c("QUE-ADU", "ONT-ADU")}.
#' @param mixture (optional) But required if bseline was selected. A character 
#' string with your mixture id. e.g. \code{c("QUE-JUV", "ONT-JUV")}.

#' @param imputations Should a map-independent imputations of markers be
#' computed. Available choices are: (1) \code{FALSE} for no imputation.
#' (2) \code{"max"} to use the most frequent category for imputations.
#'  (3) \code{"rf"} using Random Forest algorithm. Default = \code{FALSE}.
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
#' @details The imputations using Random Forest requires more time to compute 
#' and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set 
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals 
#' will require 15 min.
#' The Fst is based on Weir and Cockerham 1984 equations.
#' @return When no imputation is selected a gsi_sim file is saved to the 
#' working directory. When imputation is selected 2 gsi_sim files are saved to
#' the working directory.
#' @export
#' @rdname vcf2gsi_sim
#' @import reshape2
#' @import dplyr
#' @import foreach
#' @import parallel
#' @import doParallel
#' @importFrom stringr str_pad
#' @references Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008) 
#' An improved method for predicting the accuracy of genetic stock identification.
#' Canadian Journal of Fisheries and Aquatic Sciences 65, 7:1475-1486.
#' @references Anderson, E. C. (2010) Assessing the power of informative subsets of 
#' loci for population assignment: standard methods are upwardly biased.
#' Molecular ecology resources 10, 4:701-710.
#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @references Weir BS, Cockerham CC (1984) Estimating F-Statistics for the 
#' Analysis of Population Structure. Evolution, 38, 1358–1370.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


vcf2gsi_sim <- function(vcf.file, 
                        whitelist.markers = NULL,
                        snp.LD = NULL,
                        common.markers = FALSE,
                        marker.number = "all",
                        sampling.method,
                        THL,
                        iterations = 10,
                        blacklist.id = NULL,
                        subsample = NULL,
                        gsi_sim.filename = "gsi_sim_data.txt",
                        keep.gsi.files,
                        pop.levels, pop.labels, pop.id.start, pop.id.end,
                        baseline = NULL, 
                        mixture = NULL,
                        imputations = FALSE,
                        imputations.group = "populations",
                        num.tree = 100,
                        iteration.rf = 10,
                        split.number = 100,
                        verbose = FALSE,
                        parallel.core = 2) {
  
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
  nal <- NULL
  ALLELES_GROUP <- NULL
  N_IND_GENE <- NULL
  P <- NULL
  pb <- NULL
  nal_sq <- NULL
  nal_sq_sum <- NULL
  nal_sq_sum_nt <- NULL
  npl <- NULL
  het <- NULL
  mho <- NULL
  mhom <- NULL
  dum <- NULL
  dum1 <- NULL
  SSG <- NULL
  ntal <- NULL
  SSP <- NULL
  ntalb <- NULL
  SSi <- NULL
  MSI <- NULL
  sigw <- NULL
  MSP <- NULL
  siga <- NULL
  sigb <- NULL
  lsiga <- NULL
  lsigb <- NULL
  lsigw <- NULL
  FST <- NULL
  KEEPER <- NULL
  ASSIGN <- NULL
  OTHERS <- NULL
  CURRENT <- NULL
  INFERRED <- NULL
  SECOND_BEST_POP <- NULL
  SCORE <- NULL
  SECOND_BEST_SCORE <- NULL
  MARKER_NUMBER <- NULL
  TOTAL <- NULL
  ASSIGNMENT_PERC <- NULL
  MISSING_DATA <- NULL
  MEAN <- NULL
  METHOD <- NULL
  SE <- NULL
  
  
  
  # Control of arguments
  # Create a folder based on filename to save the output files ------------------
  if (imputations == "FALSE") {
    message("vcf2gsi_sim: without imputation.")
    directory <- stri_paste(getwd(),"/","method_", sampling.method, "_no_imputations","/", sep = "")
    dir.create(file.path(directory))
    
  } else {
    message("vcf2gsi_sim: with imputations.")
    directory <- stri_paste(getwd(),"/","method_", sampling.method, "_imputations_", imputations,"_", imputations.group,"/", sep = "")
    dir.create(file.path(directory))
  }
  
  if(missing(subsample) | is.null(subsample)) {
    subsample <- FALSE
  }
  
  if (missing(gsi_sim.filename)){
    gsi_sim.filename <- "gsi_sim_data.txt"
  }
  
  if(missing(keep.gsi.files)){
    keep.gsi.files <- FALSE
  }
  
  
  # To work inside foreach ...
  if(missing(mixture) | is.null(mixture)){
    # if(missing(mixture)){
    mixture <- FALSE
  } else{
    mixture <- mixture
  }
  
  if(missing(baseline)| is.null(baseline)){
    # if(missing(baseline)){
    baseline <- FALSE
  } else{
    baseline <- baseline
  }
  
  if(missing(THL)){
    THL <- 1
  }
  
  # Import/read VCF ------------------------------------------------------------- 
  message("Importing the VCF...")
  
  vcf <- read_delim(
    vcf.file, 
    delim = "\t", 
    comment = "##",
    progress = interactive()
  ) %>% 
    select(-c(QUAL, FILTER, INFO)) %>% 
    rename(LOCUS = ID, CHROM = `#CHROM`) %>% 
    mutate(
      CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
    )
  
  # Detect STACKS version
  if(stri_detect_fixed(vcf$FORMAT[1], "AD")) {
    stacks.version <- "new"
  } else{
    stacks.version <- "old"
  }
  vcf <- vcf %>% select(-FORMAT)
  
  # Whitelist of markers -------------------------------------------------------
  
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
    vcf <- vcf %>% 
      semi_join(whitelist.markers, by = columns.names.whitelist)
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
      tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, LOCUS, POS, REF, ALT)) %>% # Gather individuals in 1 colummn
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
  
  # subsample ------------------------------------------------------------------
  if (subsample == FALSE){
    vcf <- vcf 
  } else{
    ind.pop.df <- vcf %>% select(POP_ID, INDIVIDUALS) %>% distinct(POP_ID, INDIVIDUALS)
    if (subsample > 1){
      message(stri_c("subsampling your dataset with ", subsample, " individuals per population", sep = ""))
      subsample.select <- ind.pop.df %>%
        group_by(POP_ID) %>% 
        sample_n(subsample, replace = FALSE) %>% # sampling individuals for each pop
        arrange(POP_ID, INDIVIDUALS)
    } 
    if (subsample < 1){ # proportion
      message(stri_c("subsampling your dataset with ", subsample * 100, " % of your individuals per population", sep = ""))
      subsample.select <- ind.pop.df %>%
        group_by(POP_ID) %>% 
        sample_frac(subsample, replace = FALSE) %>% # sampling individuals for each pop
        arrange(POP_ID, INDIVIDUALS)
    }
    vcf <- vcf %>% 
      semi_join(subsample.select, by = c("POP_ID", "INDIVIDUALS"))
  }
  
  # dump unused object
  blacklist.id <- NULL
  whitelist.markers <- NULL
  ind.pop.df <- NULL
  subsample.select <- NULL
  
  # Tidy VCF ------------ -----------------------------------------------------
  message("Tidy vcf into factory for conversion into gsi_sim ...")
  
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
  
  # LD control... keep only 1 SNP per haplotypes/reads (optional) --------------
  if(missing(snp.LD) | is.null(snp.LD)){
    vcf <- vcf
  } else{
    snp.locus <- vcf %>% select(LOCUS, POS) %>% distinct(POS)
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
  
  # Unique markers id: combine CHROM, LOCUS and POS into MARKERS ---------------
  vcf <- vcf %>%
    mutate(
      POS = stri_pad_left(str = POS, width = 8, pad = "0"),
      LOCUS = stri_pad_left(str = LOCUS, width = 8, pad = "0")
    ) %>%
    arrange(CHROM, LOCUS, POS) %>% 
    tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
  
  # Markers in common between all populations (optional)------------------------
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
      select(MARKERS) %>% 
      distinct(MARKERS)
    
    message(stri_c("Number of original markers = ", n_distinct(vcf$MARKERS), "\n", "Number of markers present in all the populations = ", n_distinct(pop.filter$MARKERS), "\n", "Number of markers removed = ", n_distinct(vcf$MARKERS) - n_distinct(pop.filter$MARKERS)))
    vcf <- vcf %>% semi_join(pop.filter, by = "MARKERS")
  }
  
  # Change the genotype coding for easier integration in downstream conversion to gsi_sim
  # this step will be use for the imputation part, so we stop the object assignment to vcf with this step
  vcf <- vcf %>% 
    mutate(
      REF= stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE), # replace nucleotide with numbers
      ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE),# replace nucleotide with numbers
      GT = ifelse(GT == "0/0", stri_c(REF, REF, sep = "_"),
                  ifelse(GT == "1/1",  stri_c(ALT, ALT, sep = "_"),
                         ifelse(GT == "0/1", stri_c(REF, ALT, sep = "_"),
                                ifelse(GT == "1/0", stri_c(ALT, REF, sep = "_"), "0_0")
                         )
                  )
      )
    ) %>% 
    arrange(MARKERS, POP_ID) %>% 
    select(-c(REF, ALT))
  message("step 2/3: completed")
  
  # more prep for the no imputation part
  gsim.prep <- vcf %>%
    tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
    tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID))
  message("step 3/3: completed")
  
  # save.image("assignment.lobster.RData")  
  # load("assignment.lobster.RData")
  
  # Imputations: gsi_sim with imputed genotypes using Random Forest or the most frequent allele-----------
  if (imputations != "FALSE"){
    
    vcf.prep <- vcf %>%
      mutate(
        GT = stri_replace_all_fixed(GT, pattern = "0_0", replacement = "NA", vectorize_all = FALSE),
        GT = replace(GT, which(GT == "NA"), NA)
      ) %>% 
      dcast(INDIVIDUALS + POP_ID ~ MARKERS, value.var = "GT") %>% 
      arrange(POP_ID, INDIVIDUALS)
    
    # vcf <- NULL # remove unused object
    
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
      
      # imputations by populations (default) or globally------------------------
      # default by pop
      if (missing(imputations.group) == "TRUE" | imputations.group == "populations"){
        message("Imputations computed by populations, take a break...")
        df.split.pop <- split(x = vcf.prep, f = vcf.prep$POP_ID) # slip data frame by population
        pop.list <- names(df.split.pop) # list the pop
        imputed.dataset <-list() # create empty list 
        # for (i in pop.list) {
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
    
    # transform the imputed dataset into gsi_sim  ------------------------
    
    message("Imputed VCF into factory for conversion into gsi_sim...")
    gsi.prep.imp <- suppressWarnings(
      vcf.imp %>%
        tidyr::gather(key = MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>% # make tidy
        tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
        tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) # make tidy
    )
  } # End imputations
  
  # Sampling of markers -------------------------------------------------------
  
  # get the unique list of markers for argument 'marker.number'
  # if "all" is present in the list, change to the maximum number of markers
  unique.markers <- gsim.prep %>% select(MARKERS) %>% distinct(MARKERS) %>% arrange(MARKERS)
  marker.number <- stri_replace_all_fixed(str = marker.number, pattern = "all", replacement = nrow(unique.markers), vectorize_all = TRUE)
  
  # Random method
  if (sampling.method == "random"){
    # keep track of random seed number
    filename.seed <- "random_seed_numbers_gsi_sim.tsv"
    random.seed.table <- data_frame(MARKER_NUMBER = integer(0), ITERATIONS = integer(0), RANDOM_NUMBER = integer(0)) #create an empty dataframe
    write_tsv(x = random.seed.table, path = paste0(directory,filename.seed), col_names = TRUE, append = FALSE) #create an empty file
    
    # iterations with random markers
    # Number of times to repeat the sampling of markers
    if (missing(iterations)){
      iterations.list <- 1:10
    } else {
      iterations.list <- iterations
    }
    # Marker number loop
    for (i in marker.number) {
      
      message(paste("Marker number: ", i))
      # i <- 100
      i <- as.numeric(i)
      
      # Start cluster registration backend using n - 1 CPU
      cl <- makeCluster(detectCores() - 1)
      registerDoParallel(cl, cores = detectCores() - 1)
      
      
      j <- NULL
      foreach(j=iterations.list, .packages = c("magrittr", "plyr", "dplyr", "tidyr", "stringi", "readr", "randomForestSRC", "reshape2")) %dopar% {
        # j <-100
        
        # Sampling markers ----------------------------------------------------
        if (sampling.method == "random"){
          # Set seed for random sampling
          # random.seed <- sample(x = 1:10000, size = 1)
          random.seed <- i + j
          random.seed.table <- data_frame(MARKER_NUMBER = i, ITERATIONS = j, RANDOM_NUMBER = random.seed)
          write_tsv(x = random.seed.table, path = paste0(directory, filename.seed), col_names = FALSE, append = TRUE)
          
          # RANDOM sample markers
          select.markers <- sample_n(tbl = unique.markers, size = i, replace = FALSE) %>% 
            arrange(MARKERS)
        } # end random sampling
        
        # No imputation ----------------------------------------------------------
        gsim.prep.no.imputation.random <- gsim.prep %>% 
          semi_join(select.markers, by = "MARKERS") %>% 
          arrange(MARKERS) %>%  # make tidy
          tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>% 
          arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
          dcast(INDIVIDUALS + POP_ID ~ MARKERS_ALLELES, value.var = "GT") %>% 
          arrange(POP_ID, INDIVIDUALS)
        
        # get the list of loci after filter  
        loci <- unique(select.markers$MARKERS)
        n.markers <- length(loci)
        
        if(mixture == FALSE){
          message("No baseline or mixture data")
          gsim.prep.no.imputation.random <- gsim.prep.no.imputation.random
        } else {
          baseline.data <- suppressWarnings(
            gsim.prep.no.imputation.random %>%
              filter(POP_ID %in% baseline) %>% 
              arrange(POP_ID) %>% 
              mutate(POP_ID = droplevels(POP_ID))
          )
          
          mixture.data <- suppressWarnings(
            gsim.prep.no.imputation.random %>%
              filter(POP_ID %in% mixture) %>% 
              arrange(POP_ID) %>% 
              mutate(POP_ID = droplevels(POP_ID))
          )
        }
        # results no imputation-------------------------------------------------------
        if(mixture == FALSE){
          n.individuals <- n_distinct(gsim.prep.no.imputation.random$INDIVIDUALS)
          
          # Create a vector with the population ordered by levels
          pop <- gsim.prep.no.imputation.random$POP_ID
          gsim.prep.no.imputation.random <- suppressWarnings(gsim.prep.no.imputation.random %>% select(-POP_ID))
          
          # split gsi_sim by populations
          gsi_sim.split <- split(gsim.prep.no.imputation.random, pop)
          
          # Write the file in gsi_sim format
          message("Output...No imputation")
          
          # gsi_sim filename
          filename <- gsi_sim.filename
          marker.number.in.filename <- stri_paste(i, j,"txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = marker.number.in.filename)
          
          # Add "_no.imputation" to the filename
          filename <- stri_replace_all_fixed(filename,
                                             pattern = ".txt",
                                             replacement = "_no.imputation.txt")
          
          # Line 1: number of individuals and the number of markers
          line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
          write.table(line1_gsi_sim, file = paste0(directory, filename), col.names = FALSE, row.names = FALSE, quote = FALSE)
          
          # Markers names
          loci.table <- as.data.frame(loci)
          write_delim(x = loci.table, path = paste0(directory, filename), delim = "\n", append = TRUE, col_names = FALSE)
          
          # remaining lines, individuals and genotypes
          for (k in levels(pop)) {
            pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
            write_delim(x = pop.line, path = paste0(directory, filename), delim = "\n", append = TRUE, col_names = FALSE)
            write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename), delim = " ", append = TRUE, col_names = FALSE)
          }
          message(stri_paste("Data file (no imputation):", filename, "\nWritten to the working directory:", directory, sep = " "))
          
        } else { # baseline and mixture output
          # 1. Baseline
          n.individuals <- n_distinct(baseline.data$INDIVIDUALS)
          pop <- baseline.data$POP_ID
          baseline.data <- baseline.data %>% select(-POP_ID)
          
          # split gsi_sim by populations
          gsi_sim.split <- split(baseline.data, pop)
          
          # gsi_sim filename
          filename <- gsi_sim.filename
          marker.number.in.filename <- stri_paste("baseline", i, j, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = marker.number.in.filename)
          
          # Add "_no.imputation" to the filename
          filename <- stri_replace_all_fixed(filename,
                                             pattern = ".txt",
                                             replacement = "_no.imputation.txt")
          
          
          # Line 1: number of individuals and the number of markers
          line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
          write.table(line1_gsi_sim, file = paste0(directory, filename), col.names = FALSE, row.names = FALSE, quote = FALSE)
          
          # Markers names
          loci.table <- as.data.frame(loci)
          write_delim(x = loci.table, path = paste0(directory, filename), delim = "\n", append = TRUE, col_names = FALSE)
          
          # remaining lines, individuals and genotypes
          for (k in levels(pop)) {
            pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
            write_delim(x = pop.line, path = paste0(directory, filename), delim = "\n", append = TRUE, col_names = FALSE)
            write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename), delim = " ", append = TRUE, col_names = FALSE)
          }
          message(stri_paste("Baseline data file:", filename, "\nWritten to the working directory:", directory, sep = " "))
          
          # 2. Mixture
          n.individuals <- n_distinct(mixture.data$INDIVIDUALS)
          pop <- mixture.data$POP_ID
          mixture.data <- mixture.data %>% select(-POP_ID)
          
          # split gsi_sim by populations
          gsi_sim.split <- split(mixture.data, pop)
          
          # gsi_sim filename
          filename <- gsi_sim.filename
          marker.number.in.filename <- stri_paste("mixture", i, j, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = marker.number.in.filename)
          
          # Add "_no.imputation" to the filename
          filename <- stri_replace_all_fixed(filename,
                                             pattern = ".txt",
                                             replacement = "_no.imputation.txt")
          
          # Line 1: number of individuals and the number of markers
          line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
          write.table(line1_gsi_sim, file = paste0(directory, filename), col.names = FALSE, row.names = FALSE, quote = FALSE)
          
          # Markers names
          loci.table <- as.data.frame(loci)
          write_delim(x = loci.table, path = paste0(directory, filename), delim = "\n", append = TRUE, col_names = FALSE)
          
          # remaining lines, individuals and genotypes
          for (k in levels(pop)) {
            pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
            write_delim(x = pop.line, path = paste0(directory, filename), delim = "\n", append = TRUE, col_names = FALSE)
            write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename), delim = " ", append = TRUE, col_names = FALSE)
          }
          message(stri_paste("Mixture data file:", filename, "\nWritten to the working directory:", directory, sep = " "))
        }
        # with imputed dataset -------------------------------------------------
        gsim.prep.imputation.random <- suppressWarnings(
          gsi.prep.imp %>% 
            semi_join(select.markers, by = "MARKERS") %>% 
            arrange(MARKERS) %>%  # make tidy
            tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>% 
            arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
            dcast(INDIVIDUALS + POP_ID ~ MARKERS_ALLELES, value.var = "GT") %>%
            mutate(POP_ID = factor(POP_ID, levels = pop.labels, ordered = TRUE), POP_ID = droplevels(POP_ID)) %>% 
            arrange(POP_ID, INDIVIDUALS)
        )
        if(mixture == FALSE){
          message("No baseline or mixture data")
          gsim.prep.imputation.random <- gsim.prep.imputation.random
        } else {
          # Baseline dataset
          baseline.data <- suppressWarnings(
            gsim.prep.imputation.random %>%
              filter(POP_ID %in% baseline) %>% 
              arrange(POP_ID) %>% 
              mutate(POP_ID = droplevels(POP_ID))
          )
          # Mixture dataset
          mixture.data <- suppressWarnings(
            gsim.prep.imputation.random %>%
              filter(POP_ID %in% mixture) %>% 
              arrange(POP_ID) %>% 
              mutate(POP_ID = droplevels(POP_ID))
          )
        }
        # results imputation----------------------------------------------------
        if(mixture == FALSE){
          n.individuals <- n_distinct(gsim.prep.imputation.random$INDIVIDUALS)
          
          # Create a vector with the population ordered by levels
          pop <- gsim.prep.imputation.random$POP_ID
          gsim.prep.imputation.random <- suppressWarnings(gsim.prep.imputation.random %>% select(-POP_ID))
          
          # split gsi_sim by populations
          gsi_sim.split <- split(gsim.prep.imputation.random, pop)
          
          # Write the file in gsi_sim format 
          message("Output of imputed data...")
          
          # gsi_sim filename
          filename <- gsi_sim.filename
          marker.number.in.filename <- stri_paste(i, j, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = marker.number.in.filename)
          
          # Add "_imputed" to the filename
          filename.imp <- stri_replace_all_fixed(filename,
                                                 pattern = ".txt",
                                                 replacement = "_imputed.txt")
          
          # Line 1: number of individuals and the number of markers
          line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
          write.table(line1_gsi_sim, file = paste0(directory, filename.imp), col.names = FALSE, row.names = FALSE, quote = FALSE)
          
          # Markers names
          loci.table <- as.data.frame(loci)
          write_delim(x = loci.table, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
          
          # remaining lines, individuals and genotypes
          for (k in levels(pop)) {
            pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
            write_delim(x = pop.line, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
            write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename.imp), delim = " ", append = TRUE, col_names = FALSE)
          }
          message(stri_paste("Data file (with imputations):", filename.imp, "\nWritten to the working directory:", directory, sep = " "))
          
        } else { # baseline and mixture output
          # 1. Baseline
          n.individuals <- n_distinct(baseline.data$INDIVIDUALS)
          pop <- baseline.data$POP_ID
          baseline.data <- baseline.data %>% select(-POP_ID)
          
          # split gsi_sim by populations
          gsi_sim.split <- split(baseline.data, pop)
          
          # gsi_sim filename
          filename <- gsi_sim.filename
          marker.number.in.filename <- stri_paste("baseline", i, j, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = marker.number.in.filename)
          
          # Add "_imputed" to the filename
          filename.imp <- stri_replace_all_fixed(filename,
                                                 pattern = ".txt",
                                                 replacement = "_imputed.txt")
          
          # Line 1: number of individuals and the number of markers
          line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
          write.table(line1_gsi_sim, file = paste0(directory, filename.imp), col.names = FALSE, row.names = FALSE, quote = FALSE)
          
          # Markers names
          loci.table <- as.data.frame(loci)
          write_delim(x = loci.table, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
          
          # remaining lines, individuals and genotypes
          for (k in levels(pop)) {
            pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
            write_delim(x = pop.line, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
            write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename.imp), delim = " ", append = TRUE, col_names = FALSE)
          }
          message(stri_paste("Baseline data file (with imputations):", filename, "\nWritten to the working directory:", directory, sep = " "))
          
          # 2. Mixture
          n.individuals <- n_distinct(mixture.data$INDIVIDUALS)
          pop <- mixture.data$POP_ID
          mixture.data <- mixture.data %>% select(-POP_ID)
          
          # split gsi_sim by populations
          gsi_sim.split <- split(mixture.data, pop)
          
          # gsi_sim filename
          filename <- gsi_sim.filename
          marker.number.in.filename <- stri_paste("mixture", i, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = marker.number.in.filename)
          
          # Add "_imputed" to the filename
          filename.imp <- stri_replace_all_fixed(filename,
                                                 pattern = ".txt",
                                                 replacement = "_imputed.txt")
          
          # Line 1: number of individuals and the number of markers
          line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
          write.table(line1_gsi_sim, file = paste0(directory, filename.imp), col.names = FALSE, row.names = FALSE, quote = FALSE)
          
          # Markers names
          loci.table <- as.data.frame(loci)
          write_delim(x = loci.table, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
          
          # remaining lines, individuals and genotypes
          for (k in levels(pop)) {
            pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
            write_delim(x = pop.line, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
            write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename.imp), delim = " ", append = TRUE, col_names = FALSE)
          }
          message(stri_paste("Mixture data file (with imputations):", filename, "\nWritten to the working directory:", directory, sep = " "))
        }
      } # End of iterations
      # close parallel connection settings
      stopCluster(cl)
    } # End marker number
  } # end method random
  
  # Ranked method
  if (sampling.method == "ranked"){
    message("Using THL method, ranking Fst with training samples...")
    message("Holdout samples saved in the directory")
    ind.pop.df<- vcf %>% select(POP_ID, INDIVIDUALS) %>% distinct(POP_ID, INDIVIDUALS)# List of all individuals
    
    # THL selection-------------------------------------------------------------
    if(THL == 1){
      # Will go through the individuals in the list one by one.
      iterations.list <- ind.pop.df$INDIVIDUALS 
      
      # Keep track of holdout individuals
      holdout.individuals <- ind.pop.df %>% 
        mutate(ITERATIONS = stri_c("HOLDOUT", seq(1:n()), sep = "_"))
      write_tsv(x = holdout.individuals, path = paste0(directory,"holdout.individuals.tsv"), col_names = TRUE, append = FALSE)
    } else {
      # Create x (iterations) list of y (THL) proportion of individuals per pop.
      if(stri_detect_fixed(THL, ".") & THL < 1) {
        # iterations <- 30 # test
        # THL <- 0.15 # test
        iterations.list <- list()
        for (x in 1:iterations){
          holdout.individuals <- ind.pop.df %>% 
            group_by(POP_ID) %>%
            sample_frac(THL, replace = FALSE) %>%  # sampling fraction for each pop
            arrange(POP_ID, INDIVIDUALS) %>% 
            ungroup() %>% 
            select(INDIVIDUALS) %>% 
            mutate(ITERATIONS = rep(x, n()))
          iterations.list[[x]] <- holdout.individuals
        }
        holdout.individuals <- NULL
        
        # Keep track of holdout individuals
        holdout.individuals <- as.data.frame(bind_rows(iterations.list))
        write_tsv(x = holdout.individuals, 
                  path = paste0(directory,"holdout.individuals.tsv"), 
                  col_names = TRUE, 
                  append = FALSE
        )
      }
      
      # Create x (iterations) list of y (THL) individuals per pop.
      if(THL > 1) {
        iterations.list <- list()
        for (x in 1:iterations){
          holdout.individuals <- ind.pop.df %>% 
            group_by(POP_ID) %>%
            sample_n(THL, replace = FALSE) %>% # sampling individuals for each pop
            arrange(POP_ID, INDIVIDUALS) %>% 
            ungroup() %>% 
            select(INDIVIDUALS) %>% 
            mutate(ITERATIONS = rep(x, n()))
          iterations.list[[x]] <- holdout.individuals
        }
        holdout.individuals <- NULL
        holdout.individuals <- as.data.frame(bind_rows(iterations.list))
        write_tsv(x = holdout.individuals, 
                  path = paste0(directory,"holdout.individuals.tsv"), 
                  col_names = TRUE, 
                  append = FALSE
        )
      }
    } # end tracking holdout individuals
    
    # Fst: Weir & Cockerham 1984 -----------------------------------------------
    # Currently the Fst is computed on the non imputed data
    # Would be easy to computed on imputed data
    
    # Fst Function
    fst_WC84 <- function(data, holdout.individuals){
      pop.number <- n_distinct(data$POP_ID)
      
      data.genotyped <- data %>%
        filter(GT != "0_0") %>% # remove missing genotypes
        filter(!INDIVIDUALS %in% holdout.individuals) # remove supplementary individual before ranking markers with Fst
      
      # Nombre de population pour chaque locus
      n.pop.locus <- data.genotyped %>%
        select(MARKERS, POP_ID) %>%
        group_by(MARKERS) %>% 
        distinct(POP_ID) %>% 
        tally %>% 
        rename(npl = n)
      
      # Fréquence de chaque allele pour chaque locus et pop
      ind.count.locus <- data.genotyped %>% 
        group_by(MARKERS) %>% 
        tally
      
      ind.count.locus.pop <- data.genotyped %>% 
        group_by(POP_ID, MARKERS) %>% 
        tally %>%
        rename(nal = n) %>% 
        mutate(
          nal_sq = nal^2,
          N_IND_GENE = nal*2
        )
      
      #common
      freq.al.locus <- data.genotyped %>%
        mutate(
          A1 = stri_sub(GT, 1, 1),
          A2 = stri_sub(GT, 3, 3)
        ) %>% 
        select(-GT) %>% 
        tidyr::gather(key = ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, POP_ID, MARKERS))
      
      #pop
      freq.al.locus.pop <- freq.al.locus %>% 
        group_by(POP_ID, MARKERS, ALLELES) %>% 
        tally %>% 
        full_join(ind.count.locus.pop, by = c("POP_ID", "MARKERS")) %>% 
        mutate(P = n/N_IND_GENE) %>% # Freq. Allele per pop
        select(POP_ID, MARKERS, ALLELES, P) %>% 
        group_by(MARKERS, ALLELES) %>% 
        tidyr::spread(POP_ID, P) %>% 
        tidyr::gather(key = POP_ID, value = P, -c(MARKERS, ALLELES)) %>% 
        mutate(P = as.numeric(stri_replace_na(str = P, replacement = 0))) %>% 
        full_join(ind.count.locus.pop, by = c("POP_ID", "MARKERS"))
      
      # global
      freq.al.locus.global <- freq.al.locus %>% 
        group_by(MARKERS, ALLELES) %>% 
        tally %>% 
        full_join(ind.count.locus%>% rename(N = n), by = "MARKERS") %>% 
        mutate(pb = n/(2*N)) %>% # Global Freq. Allele
        select(MARKERS, ALLELES, pb)
      
      # correction
      # nombre total d'individus corrigé par locus
      mean.n.pop.corrected.per.locus <- ind.count.locus.pop %>% 
        group_by(MARKERS) %>%
        summarise(nal_sq_sum = sum(nal_sq, na.rm = TRUE)) %>% 
        full_join(ind.count.locus, by = "MARKERS") %>% 
        mutate(nal_sq_sum_nt = (n - nal_sq_sum/n)) %>% 
        full_join(n.pop.locus, by = "MARKERS") %>% 
        mutate(ncal = nal_sq_sum_nt/(npl-1)) %>% 
        select(MARKERS, ncal)
      
      # Nombre d'individus par allele et locus global
      # Nombre corrigé d'individus par locus et allele global
      
      ncal <- freq.al.locus %>% 
        select(MARKERS, ALLELES) %>% 
        group_by(MARKERS, ALLELES) %>% 
        distinct(MARKERS, ALLELES) %>% 
        full_join(ind.count.locus, by = "MARKERS") %>% 
        rename(ntal = n) %>% 
        full_join(mean.n.pop.corrected.per.locus, by = "MARKERS") 
      
      fst.ranked <- data.genotyped %>% 
        mutate(het = ifelse(stri_sub(GT, 1, 1) != stri_sub(GT, 3, 3), 1, 0)) %>% 
        group_by(MARKERS, POP_ID) %>%
        summarise(mho = sum(het, na.rm = TRUE)) %>%  # = the number of heterozygote individuals per pop and markers
        group_by(MARKERS) %>% 
        tidyr::spread(POP_ID, mho) %>% 
        tidyr::gather(key = POP_ID, value = mho, -MARKERS) %>%
        mutate(mho = as.numeric(stri_replace_na(str = mho, replacement = 0))) %>%
        full_join(freq.al.locus.pop, by = c("POP_ID", "MARKERS")) %>%
        mutate(
          mhom = round(((2 * nal * P - mho)/2), 0),
          dum = nal * (P - 2 * P^2) + mhom
        ) %>% 
        group_by(MARKERS, ALLELES) %>% 
        full_join(freq.al.locus.global, by = c("MARKERS", "ALLELES")) %>%
        mutate(
          SSi = sum(dum, na.rm = TRUE),
          dum1 = nal * (P - pb)^2
        ) %>% 
        group_by(MARKERS, ALLELES) %>% 
        mutate(SSP = 2 * sum(dum1, na.rm = TRUE)) %>%
        group_by(MARKERS, POP_ID) %>% 
        mutate(SSG = nal * P - mhom) %>%
        group_by(MARKERS, ALLELES) %>% 
        full_join(ncal, by = c("MARKERS", "ALLELES")) %>%
        full_join(n.pop.locus, by = "MARKERS") %>%
        rename(ntalb = npl) %>%
        mutate(
          sigw = round(sum(SSG, na.rm = TRUE), 2)/ntal,
          MSP = SSP/(ntalb - 1),
          MSI = SSi/(ntal - ntalb),
          sigb = 0.5 * (MSI - sigw),
          siga = 1/2/ncal * (MSP - MSI)
        ) %>% 
        group_by(MARKERS) %>% 
        summarise(
          lsiga = sum(siga, na.rm = TRUE),
          lsigb = sum(sigb, na.rm = TRUE),
          lsigw = sum(sigw, na.rm = TRUE),
          FST = round(lsiga/(lsiga + lsigb + lsigw), 6),
          FIS = round(lsigb/(lsigb + lsigw), 6)
        ) %>% 
        arrange(desc(FST)) %>% 
        select(MARKERS, FST) %>% 
        mutate(RANKING = seq(from = 1, to = n()))
      
      # select(MARKERS, FIS, FST)
      return(fst.ranked)
    } # end Fst function
    
    # Function needed to write GSI_SIM file ------------------------------------
    write_gsi <- function (data, imputations, filename){
      
      # get some info
      n.individuals <- n_distinct(data$INDIVIDUALS) # number of individuals
      pop <- data$POP_ID # Create a vector with the population ordered by levels
      
      # prep
      data <- suppressWarnings(data %>% select(-POP_ID)) # remove pop id
      gsi_sim.split <- split(data, pop) # split gsi_sim by populations
      
      # gsi_sim filename
      filename <- filename
      marker.number.in.filename <- stri_c(iteration.id, j,"txt", sep = ".")
      filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                         replacement = marker.number.in.filename)
      
      
      # filename modification based with or without imputations
      if(imputations == FALSE) {
        message("Output...No imputation")
        # Add "_no.imputation" to the filename
        filename <- stri_replace_all_fixed(filename,
                                           pattern = ".txt",
                                           replacement = "_no.imputation.txt")
      }else{
        message("Output...With imputations")
        filename <- stri_replace_all_fixed(filename,
                                           pattern = ".txt",
                                           replacement = "_imputed.txt")
      }
      
      
      # directory <- getwd() # test
      # Line 1: number of individuals and the number of markers
      line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
      write.table(line1_gsi_sim, file = paste0(directory, filename), col.names = FALSE, row.names = FALSE, quote = FALSE)
      
      # Markers names
      loci.table <- as.data.frame(loci)
      write_delim(x = loci.table, path = paste0(directory, filename), delim = "\n", append = TRUE, col_names = FALSE)
      
      # remaining lines, individuals and genotypes
      for (k in levels(pop)) {
        pop.line <- as.data.frame(stri_c("pop", k, sep = " "))
        write_delim(x = pop.line, path = paste0(directory, filename), delim = "\n", append = TRUE, col_names = FALSE)
        write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename), delim = " ", append = TRUE, col_names = FALSE)
      }
      message(stri_c("Data file (no imputation):", filename, "\nWritten to the working directory:", directory, sep = " "))
      return(filename)
    } # end write gsi
    
    message("Starting parallel computations, progress bar not yet implemented, so keep track in the directory for activity")
    # Start cluster registration backend
    if (missing(parallel.core) == "TRUE"){
      # Automatically select all the core -1 
      cl <- makeCluster(detectCores() - 1)
      registerDoParallel(cl, cores = detectCores() - 1)
    } else {
      cl <- makeCluster(parallel.core)
      registerDoParallel(cl, cores = parallel.core)
    }
    
    # Going through the loop of holdout individuals-----------------------------
    i <- NULL
    assignment.res <- list()
    assignment.res <-foreach(i=iterations.list, .packages = c("magrittr", "plyr", "dplyr", "tidyr", "stringi", "readr", "reshape2", "purrr")) %dopar% {
      
      holdout <- data.frame(i)
      # holdout <- data.frame(iterations.list[2])
      
      # Ranked the Fst without the holdout individuals
      fst.ranked <- fst_WC84(data = vcf, holdout.individuals = holdout$INDIVIDUALS)
      
      # Saving Fst
      if(THL == 1){
        colnames(holdout) <- "INDIVIDUALS"
        iteration.id <-i
        # iteration.id <- "BON_19" #test
        fst.ranked.filename <- stri_paste("fst.ranked.", "holdout_", iteration.id, ".tsv", sep = "") 
        write_tsv(x = fst.ranked, path = paste0(directory, fst.ranked.filename), col_names = TRUE, append = FALSE)
      } else { # for THL != 1 (numbers and proportions)
        iteration.id <- unique(holdout$ITERATIONS)
        fst.ranked.filename <- stri_paste("fst.ranked.", "iteration_", iteration.id, ".tsv", sep = "") 
        write_tsv(x = fst.ranked, path = paste0(directory, fst.ranked.filename), col_names = TRUE, append = FALSE)
      }
      
      # Markers numbers loop -----------------------------------------------------
      # Create empty lists to feed the results
      for (j in marker.number) {
        message("Marker number: ", j)
        # j <- 200 # test
        j <- as.numeric(j)
        
        # Select i markers
        select.markers.fst <- filter(fst.ranked, row_number() <= j)
        select.markers <- select.markers.fst %>% select(MARKERS)
        
        if(imputations == FALSE){ # No imputation: prep ------------------------
          gsim.prep.no.imputation.ranked <- gsim.prep %>% 
            semi_join(select.markers, by = "MARKERS") %>% 
            arrange(MARKERS) %>%  # make tidy
            tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>% 
            arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
            dcast(INDIVIDUALS + POP_ID ~ MARKERS_ALLELES, value.var = "GT") %>% 
            arrange(POP_ID, INDIVIDUALS)
          
          # get the list of loci after filter  
          loci <- unique(select.markers$MARKERS)
          n.markers <- length(loci)
          
          if(mixture == FALSE){
            message("No baseline or mixture data")
            
            input <- write_gsi(data = gsim.prep.no.imputation.ranked, imputations = imputations, filename = gsi_sim.filename)
            
          } else {
            # Baseline
            baseline.data <- suppressWarnings(
              gsim.prep.no.imputation.ranked %>%
                filter(POP_ID %in% baseline) %>% 
                arrange(POP_ID) %>% 
                mutate(POP_ID = droplevels(POP_ID))
            )
            
            # gsi_sim baseline filename
            baseline.filename <- gsi_sim.filename
            marker.number.in.filename <- stri_paste("baseline", iteration.id, j, "txt", sep = ".")
            baseline.filename <- stri_replace_all_fixed(baseline.filename, pattern = "txt",
                                                        replacement = marker.number.in.filename)
            
            # save file
            baseline.input <- write_gsi(data = baseline.data, imputations = imputations, filename = baseline.filename)
            message(stri_paste("Baseline data file:", filename, "\nWritten to the working directory:", directory, sep = " "))
            
            # Mixture
            mixture.data <- suppressWarnings(
              gsim.prep.no.imputation.ranked %>%
                filter(POP_ID %in% mixture) %>% 
                arrange(POP_ID) %>% 
                mutate(POP_ID = droplevels(POP_ID))
            )
            
            # gsi_sim mixture filename
            mixture.filename <- gsi_sim.filename
            marker.number.in.filename <- stri_paste("baseline", iteration.id, j, "txt", sep = ".")
            mixture.filename <- stri_replace_all_fixed(mixture.filename, pattern = "txt",
                                                       replacement = marker.number.in.filename)
            
            # save file
            mixture.input <- write_gsi(data = mixture.data, imputations = imputations, filename = mixture.filename)
            message(stri_paste("Mixture data file:", filename, "\nWritten to the working directory:", directory, sep = " "))
          } # end writing gsi files to disk
          
          # Run GSI_SIM ------------------------------------------------------------
          if(mixture == FALSE){
            input.gsi <- stri_c(directory,input)
            output.gsi <- stri_replace_all_fixed(input.gsi, pattern = "txt", replacement = "output.txt")
            system(paste("gsisim -b", input.gsi, "--self-assign > ", output.gsi))
          } else{
            message("need to finish the function for this option :)")
          }
          # Option remove the input file from directory to save space
          if(keep.gsi.files == FALSE){
            file.remove(input.gsi)
          }
          
          # Get Assignment results -------------------------------------------------
          
          # Keep track of the holdout individual
          
          holdout.id <- holdout$INDIVIDUALS
          
          # Number of markers
          n.locus <- j
          
          assignment <- suppressWarnings(
            read_delim(output.gsi, col_names = "ID", delim = "\t") %>%
              tidyr::separate(ID, c("KEEPER", "ASSIGN"), sep = ":/", extra = "warn") %>% 
              filter(KEEPER == "SELF_ASSIGN_A_LA_GC_CSV") %>%
              tidyr::separate(ASSIGN, c("INDIVIDUALS", "ASSIGN"), sep = ";", extra = "merge") %>% 
              tidyr::separate(ASSIGN, c("INFERRED", "OTHERS"), sep = ";", convert = TRUE, numerals = "no.loss", extra = "merge") %>% 
              tidyr::separate(OTHERS, c("SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss", extra = "merge") %>% 
              tidyr::separate(OTHERS, c("SECOND_BEST_POP", "OTHERS"), sep = ";", convert = TRUE, numerals = "no.loss", extra = "merge") %>% 
              tidyr::separate(OTHERS, c("SECOND_BEST_SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss") %>%
              mutate(
                CURRENT = factor(stri_sub(INDIVIDUALS, pop.id.start, pop.id.end), levels = pop.levels, labels = pop.labels, ordered = T),
                CURRENT = droplevels(CURRENT),
                INFERRED = factor(INFERRED, levels = pop.labels, ordered = T),
                INFERRED = droplevels(INFERRED),
                SECOND_BEST_POP = factor(SECOND_BEST_POP, levels = pop.labels, ordered = T),
                SECOND_BEST_POP = droplevels(SECOND_BEST_POP),
                SCORE = round(SCORE, 2),
                SECOND_BEST_SCORE = round(SECOND_BEST_SCORE, 2),
                MARKER_NUMBER = as.numeric(rep(n.locus, n()))
              ) %>% 
              select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER) %>% 
              arrange(CURRENT) %>% 
              filter(INDIVIDUALS %in% holdout.id) %>% 
              select(INDIVIDUALS, CURRENT, INFERRED, SCORE, MARKER_NUMBER)
          )
          # write_tsv(x = assignment, path = "assignment.tsv", col_names = TRUE)
          if(keep.gsi.files == FALSE){
            file.remove(output.gsi)
          }
          
          if(THL != 1){
            assignment <- assignment %>% 
              group_by(CURRENT, INFERRED, MARKER_NUMBER) %>% 
              tally %>% 
              group_by(CURRENT, MARKER_NUMBER) %>% 
              mutate(TOTAL = sum(n)) %>% 
              ungroup() %>% 
              mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>%
              select(CURRENT, INFERRED, MARKER_NUMBER, ASSIGNMENT_PERC) %>%
              group_by(CURRENT, MARKER_NUMBER) %>% 
              tidyr::spread(INFERRED, ASSIGNMENT_PERC) %>% 
              tidyr::gather(INFERRED, ASSIGNMENT_PERC, -c(CURRENT, MARKER_NUMBER)) %>% 
              mutate(ASSIGNMENT_PERC = as.numeric(stri_replace_na(ASSIGNMENT_PERC, replacement = 0))) %>%
              filter(as.character(CURRENT) == as.character(INFERRED)) %>% 
              select(CURRENT, INFERRED, ASSIGNMENT_PERC, MARKER_NUMBER) %>% 
              mutate(ITERATIONS = rep(unique(holdout$ITERATIONS), n()))
          }
          
          #compile assignment results each marker number for the iteration
          j <- as.character(j)
          assignment.res[[j]] <- assignment
        } # end no imputation in marker number loop
        
        if(imputations != FALSE){ # imputed dataset: prep-----------------------
          
          # require some additional work 
          # iteration.id to replace i during filenaming...
          
          gsim.prep.imputation.ranked <- suppressWarnings(
            gsi.prep.imp %>% 
              semi_join(select.markers, by = "MARKERS") %>% 
              arrange(MARKERS) %>%  # make tidy
              tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>% 
              arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
              dcast(INDIVIDUALS + POP_ID ~ MARKERS_ALLELES, value.var = "GT") %>%
              mutate(POP_ID = factor(POP_ID, levels = pop.labels, ordered = TRUE), POP_ID = droplevels(POP_ID)) %>% 
              arrange(POP_ID, INDIVIDUALS)
          )
          if(mixture == FALSE){
            message("No baseline or mixture data")
            gsim.prep.imputation.random <- gsim.prep.imputation.random
          } else {
            # Baseline dataset
            baseline.data <- suppressWarnings(
              gsim.prep.imputation.random %>%
                filter(POP_ID %in% baseline) %>% 
                arrange(POP_ID) %>% 
                mutate(POP_ID = droplevels(POP_ID))
            )
            # Mixture dataset
            mixture.data <- suppressWarnings(
              gsim.prep.imputation.random %>%
                filter(POP_ID %in% mixture) %>% 
                arrange(POP_ID) %>% 
                mutate(POP_ID = droplevels(POP_ID))
            )
          }
          
          if(mixture == FALSE){# imputed datasets : output----------------------
            n.individuals <- n_distinct(gsim.prep.imputation.random$INDIVIDUALS)
            
            # Create a vector with the population ordered by levels
            pop <- gsim.prep.imputation.random$POP_ID
            gsim.prep.imputation.random <- suppressWarnings(gsim.prep.imputation.random %>% select(-POP_ID))
            
            # split gsi_sim by populations
            gsi_sim.split <- split(gsim.prep.imputation.random, pop)
            
            # Write the file in gsi_sim format 
            message("Output of imputed data...")
            
            # gsi_sim filename
            filename <- gsi_sim.filename
            marker.number.in.filename <- stri_paste(i, j, "txt", sep = ".")
            filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                               replacement = marker.number.in.filename)
            # Add "_imputed" to the filename
            filename.imp <- stri_replace_all_fixed(filename,
                                                   pattern = ".txt",
                                                   replacement = "_imputed.txt")
            
            # Line 1: number of individuals and the number of markers
            line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
            write.table(line1_gsi_sim, file = paste0(directory, filename.imp), col.names = FALSE, row.names = FALSE, quote = FALSE)
            
            # Markers names
            loci.table <- as.data.frame(loci)
            write_delim(x = loci.table, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
            
            # remaining lines, individuals and genotypes
            for (k in levels(pop)) {
              pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
              write_delim(x = pop.line, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
              write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename.imp), delim = " ", append = TRUE, col_names = FALSE)
            }
            message(stri_paste("Data file (with imputations):", filename.imp, "\nWritten to the working directory:", directory, sep = " "))
            
          } else { # baseline and mixture output
            # 1. Baseline
            n.individuals <- n_distinct(baseline.data$INDIVIDUALS)
            pop <- baseline.data$POP_ID
            baseline.data <- baseline.data %>% select(-POP_ID)
            
            # split gsi_sim by populations
            gsi_sim.split <- split(baseline.data, pop)
            
            # gsi_sim filename
            filename <- gsi_sim.filename
            marker.number.in.filename <- stri_paste("baseline", i, j, "txt", sep = ".")
            filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                               replacement = marker.number.in.filename)
            
            # Add "_imputed" to the filename
            filename.imp <- stri_replace_all_fixed(filename,
                                                   pattern = ".txt",
                                                   replacement = "_imputed.txt")
            
            # Line 1: number of individuals and the number of markers
            line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
            write.table(line1_gsi_sim, file = paste0(directory, filename.imp), col.names = FALSE, row.names = FALSE, quote = FALSE)
            
            # Markers names
            loci.table <- as.data.frame(loci)
            write_delim(x = loci.table, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
            
            # remaining lines, individuals and genotypes
            for (k in levels(pop)) {
              pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
              write_delim(x = pop.line, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
              write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename.imp), delim = " ", append = TRUE, col_names = FALSE)
            }
            message(stri_paste("Baseline data file (with imputations):", filename, "\nWritten to the working directory:", directory, sep = " "))
            
            # 2. Mixture
            n.individuals <- n_distinct(mixture.data$INDIVIDUALS)
            pop <- mixture.data$POP_ID
            mixture.data <- mixture.data %>% select(-POP_ID)
            
            # split gsi_sim by populations
            gsi_sim.split <- split(mixture.data, pop)
            
            # gsi_sim filename
            filename <- gsi_sim.filename
            marker.number.in.filename <- stri_paste("mixture", i, "txt", sep = ".")
            filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                               replacement = marker.number.in.filename)
            
            # Add "_imputed" to the filename
            filename.imp <- stri_replace_all_fixed(filename,
                                                   pattern = ".txt",
                                                   replacement = "_imputed.txt")
            
            # Line 1: number of individuals and the number of markers
            line1_gsi_sim <- as.data.frame(stri_paste(n.individuals, n.markers, sep = " "))
            write.table(line1_gsi_sim, file = paste0(directory, filename.imp), col.names = FALSE, row.names = FALSE, quote = FALSE)
            
            # Markers names
            loci.table <- as.data.frame(loci)
            write_delim(x = loci.table, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
            
            # remaining lines, individuals and genotypes
            for (k in levels(pop)) {
              pop.line <- as.data.frame(stri_paste("pop", k, sep = " "))
              write_delim(x = pop.line, path = paste0(directory, filename.imp), delim = "\n", append = TRUE, col_names = FALSE)
              write_delim(x = gsi_sim.split[[k]], path = paste0(directory, filename.imp), delim = " ", append = TRUE, col_names = FALSE)
            }
            message(stri_paste("Mixture data file (with imputations):", filename, "\nWritten to the working directory:", directory, sep = " "))
          }
        } # end imputation in marker number loop
        # when imputation method is finished, here you should merge the no.imputation and imputed data
        # res[[i]] <-assignment
      } # End marker number loop for both with and without imputations
      return(assignment.res)
    } # End holdout individuals loop
    stopCluster(cl) # close parallel connection settings
    
    # Missing data info
    if(imputations == FALSE) {
      missing.data <- "no.imputation"
    } else if (imputations == "rf"){
      if(imputations.group == "populations"){
        missing.data <- "imputed RF populations"
      } else{
        missing.data <- "imputed RF global"
      }
    }else{
      if(imputations.group == "populations"){
        missing.data <- "imputed max populations"
      } else{
        missing.data <- "imputed max global"
      }
    }
    # summary of results
    assignment.res.summary <- suppressWarnings(
      as_data_frame(bind_rows(purrr::flatten(assignment.res))) %>% 
        mutate(
          MISSING_DATA = rep(missing.data, n()),
          METHOD = rep("THL", n())
        ))
    
    if(THL == 1){
      assignment.stats.pop <- suppressWarnings(
        assignment.res.summary %>% 
          group_by(CURRENT, INFERRED, MARKER_NUMBER, MISSING_DATA, METHOD) %>% 
          tally %>% 
          group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>% 
          mutate(TOTAL = sum(n)) %>% 
          ungroup() %>% 
          mutate(MEAN = round(n/TOTAL*100, 0)) %>% 
          filter(as.character(CURRENT) == as.character(INFERRED)) %>% 
          select(CURRENT, MEAN, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
          mutate(CURRENT = factor(CURRENT, levels = pop.labels, ordered = T)) %>% 
          group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>% 
          mutate(# Below those stats are dummy, so that I can put the stats for the overall within the same dataframe...
            SE = round(sqrt(var(MEAN)/length(MEAN)), 2), #not use
            MIN = round(min(MEAN), 2),#not use
            MAX = round(max(MEAN), 2),#not use
            MEDIAN = round(median(MEAN), 2),#not use
            QUANTILE25 = round(quantile(MEAN, 0.25), 2),#not use
            QUANTILE75 = round(quantile(MEAN, 0.75), 2)#not use
          ) %>%
          arrange(CURRENT, MARKER_NUMBER)
      )
      
      pop.levels.assignment.stats.overall <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
      
      assignment.stats.overall <- assignment.stats.pop %>% 
        group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
        rename(ASSIGNMENT_PERC = MEAN) %>%
        summarise(
          MEAN = round(mean(ASSIGNMENT_PERC), 2),
          SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
          MIN = round(min(ASSIGNMENT_PERC), 2),
          MAX = round(max(ASSIGNMENT_PERC), 2),
          MEDIAN = round(median(ASSIGNMENT_PERC), 2),
          QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
          QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2)
        ) %>% 
        mutate(CURRENT = rep("OVERALL", n())) %>% 
        arrange(CURRENT, MARKER_NUMBER)
      
      assignment.summary.stats <- bind_rows(assignment.stats.pop, assignment.stats.overall) %>% 
        mutate(CURRENT = factor(CURRENT, levels = pop.levels.assignment.stats.overall, ordered = TRUE)) %>% 
        arrange(CURRENT, MARKER_NUMBER) %>% 
        mutate(
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE
        )
    } # end THL == 1
    
    if(THL != 1){
      # summary stats
      assignment.stats.pop <- assignment.res.summary %>% 
        group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
        summarise(
          MEAN = round(mean(ASSIGNMENT_PERC), 2),
          SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
          MIN = round(min(ASSIGNMENT_PERC), 2),
          MAX = round(max(ASSIGNMENT_PERC), 2),
          MEDIAN = round(median(ASSIGNMENT_PERC), 2),
          QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
          QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2),
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE
        ) %>%
        arrange(CURRENT, MARKER_NUMBER)
      
      pop.levels.assignment.stats.overall <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
      
      assignment.stats.overall <- assignment.stats.pop %>% 
        group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
        rename(ASSIGNMENT_PERC = MEAN) %>% 
        summarise(
          MEAN = round(mean(ASSIGNMENT_PERC), 2),
          SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
          MIN = round(min(ASSIGNMENT_PERC), 2),
          MAX = round(max(ASSIGNMENT_PERC), 2),
          MEDIAN = round(median(ASSIGNMENT_PERC), 2),
          QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
          QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2),
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE
        ) %>% 
        mutate(CURRENT = rep("OVERALL", n())) %>% 
        arrange(CURRENT, MARKER_NUMBER)
      
      assignment.summary.stats <- suppressWarnings(
        bind_rows(assignment.stats.pop, assignment.stats.overall) %>% 
          mutate(CURRENT = factor(CURRENT, levels = pop.levels.assignment.stats.overall, ordered = TRUE)) %>% 
          arrange(CURRENT, MARKER_NUMBER)
      )
    } # end THL != 1
    
    # Write the tables to directory
    # assignment results
    filename.assignment.res <- stri_paste("assignment.res", missing.data, sampling.method, "tsv", sep = ".")
    write_tsv(x = assignment.res.summary, path = paste0(directory,filename.assignment.res), col_names = TRUE, append = FALSE)
    
    # assignment summary stats
    filename.assignment.sum <- stri_paste("assignment.summary.stats", missing.data, sampling.method, "tsv", sep = ".")
    write_tsv(x = assignment.summary.stats, path = paste0(directory,filename.assignment.sum), col_names = TRUE, append = FALSE)
  } # end of ranked THL method
  return(assignment.summary.stats)
} # end function




