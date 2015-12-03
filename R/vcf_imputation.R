# VCF data imputation using Random Forest


#' @name vcf_imputation
#' @title VCF data imputation
#' @description Map-independent imputation of a VCF using Random Forest or the 
#' most frequent category.
#' @param vcf.file The VCF file created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param snp.LD (optional) Minimize linkage disequilibrium (LD) by choosing 
#' among these 3 options: \code{"random"} selection, \code{"first"} or 
#' \code{"last"} SNP on the same read/haplotype. Default = \code{NULL}.
#' @param common.markers (optional) Logical. Default = \code{FALSE}. 
#' With \code{TRUE}, will keep markers present in all the populations.
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param blacklist.genotypes A blacklist of genotypes 
#' containing 3 columns header 'LOCUS', 'POS' and 'INDIVIDUALS'.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param filename The name of the file written to the directory.
#' Use the extension ".vcf" at the end. Default \code{batch_1.vcf.imputed.vcf}.
#' @param imputations For map-independent imputations available options are:
#' (1) \code{"max"} to use the most frequent category for imputations.
#'  (2) \code{"rf"} using Random Forest algorithm. Default = \code{"rf"}.
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
#' @param parallel.core (optional) The number of core for OpenMP 
#' shared-memory parallel programming of Random Forest imputations. 
#' For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.
#' @details The imputations using Random Forest requires more time to compute 
#' and can take several minutes and hours depending on the size of 
#' the dataset and polymorphism of the species used. e.g. with a 
#' low polymorphic taxa, and a data set containing 30\% missing data, 
#' 5 000 haplotypes loci and 500 individuals will require 15 min.
#' @return When no imputation is selected ....
#' @export
#' @rdname vcf_imputation
#' @import reshape2
#' @import dplyr
#' @import lazyeval
#' @importFrom stringr str_pad
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
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @seealso \code{randomForestSRC} is available on 
#' CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} 
#' and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

vcf_imputation <- function(vcf.file,
                           whitelist.loci = NULL,
                           snp.LD = NULL,
                           common.markers = FALSE,
                           blacklist.id = NULL,
                           blacklist.genotypes = NULL,
                           pop.id.start, 
                           pop.id.end, 
                           filename = "batch_1.vcf.imputed.vcf",
                           imputations = "rf",
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
  everything <- NULL
  i <- NULL
  
  # Import/read VCF ------------------------------------------------------------- 
  message("Importing the VCF...")
  
  vcf <- read_delim(
    vcf.file, 
    delim = "\t", 
    comment = "##",
    progress = interactive()
  ) %>%
    select(-c(QUAL, FILTER, INFO)) %>%
    rename(LOCUS = ID, CHROM = `#CHROM`)
  
  # detect stacks version
  if(stri_detect_fixed(vcf$FORMAT[1], "AD")) {
    stacks.version <- "new"
  } else{
    stacks.version <- "old"
  }
  
  vcf <- vcf %>% select(-FORMAT)
  
  # Whitelist-------------------------------------------------------------------
  if (missing(whitelist.loci) == "TRUE") {
    vcf <- vcf
    message("No whitelist to apply to the VCF")
  } else if (is.vector(whitelist) == "TRUE") {
    message("Filtering the VCF with the whitelist from your directory")
    whitelist.markers <- read_tsv(whitelist.markers, col_names = TRUE)
    columns.names.whitelist <- colnames(whitelist.markers)
    if("CHROM" %in% columns.names.whitelist){
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
    vcf <- vcf %>% semi_join(whitelist.markers, by = columns.names.whitelist)
  } else {
    message("Filtering the VCF with the whitelist from your global environment")
    columns.names.whitelist <- colnames(whitelist.markers)
    if("CHROM" %in% columns.names.whitelist){
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
    vcf <- vcf %>% semi_join(whitelist.markers, by = columns.names.whitelist)
  }  
  
  # Gather individuals in 1 colummn --------------------------------------------
  vcf <- tidyr::gather(vcf, INDIVIDUALS, FORMAT, -c(CHROM, LOCUS, POS, REF, ALT))
  message("Gathering individuals in 1 column")
  
  
  # Blacklist id----------------------------------------------------------------
  
  if (missing(blacklist.id) == "TRUE") {
    message("No individual blacklisted")
    blacklist.id <- NULL
    vcf <- vcf
  } else if (is.vector(blacklist.id) == "TRUE") {
    message("Using the blacklisted id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)
    vcf <- suppressWarnings(
      vcf %>% 
        anti_join(blacklist.id, by = "INDIVIDUALS")
    )
  } else {
    message("Using the blacklisted id from your global environment")
    blacklist.id <- blacklist.id
    vcf <- suppressWarnings(
      vcf %>% 
        anti_join(blacklist.id, by = "INDIVIDUALS")
    )
  }
  
  # Separate FORMAT and COVERAGE columns ---------------------------------------
  message("Tidying the VCF...")
  
  if(stacks.version == "new"){ # with new version of stacks > v.1.29
    vcf <- vcf %>%
      tidyr::separate(FORMAT, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"),
                      sep = ":", extra = "warn") %>% 
      select(-c(ALLELE_DEPTH, READ_DEPTH, GL))
    
  } else { # stacks version prior to v.1.29 had no Allele Depth field...
    message("Hum....")
    message("you are using an older version of STACKS...")
    message("It's not too late to use the last STACKS version, see STACKS change log for problems associated with older vcf files, for more details see: http://catchenlab.life.illinois.edu/stacks/")
    message("Continuing to work on tidying your VCF, for this time ;)")
    
    vcf <- vcf %>%
      tidyr::separate(FORMAT, c("GT", "READ_DEPTH", "GL"),
                      sep = ":", extra = "warn") %>% 
      select(-c(READ_DEPTH, GL))
  }
  
  # dump unused object
  whitelist <- NULL
  blacklist.id <- NULL
  
  # Keep only 1 SNP per haplotypes/reads
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
  
  # Erasing genotypes with poor coverage and/or genotype likelihood-------------
  if (is.null(blacklist.genotypes) == FALSE) {
    message("Erasing genotypes in the VCF using the list of blacklisted genotypes...")
    
    bad.geno <- read_tsv(blacklist.genotypes, col_names = T) %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      select(-POP_ID)
    
    # interesting stats.
    erased.genotype.number <- length(bad.geno$INDIVIDUALS)
    total.genotype.number <- length(vcf$GT[vcf$GT != "./."])
    genotype.percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 4), "%", sep = " ")
    
    vcf <- vcf %>%
      mutate(
        GT = ifelse(LOCUS %in% bad.geno$LOCUS & POS %in% bad.geno$POS & INDIVIDUALS %in% bad.geno$INDIVIDUALS, "./.", GT)
      )
    
    # Get the number of imputed genotypes and percentage
    imputed.genotype.number <- length(vcf$GT[vcf$GT == "./."])
    imputed.genotype.percent <- paste(round((imputed.genotype.number/(imputed.genotype.number+total.genotype.number))*100, 4), "%", sep = " ")
    
  } else {
    vcf <- vcf
    
    # interesting stats
    erased.genotype.number <- 0
    total.genotype.number <- length(vcf$GT[vcf$GT != "./."])
    genotype.percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 4), "%", sep = " ")
    imputed.genotype.number <- length(vcf$GT[vcf$GT == "./."])
    imputed.genotype.percent <- paste(round((imputed.genotype.number/(imputed.genotype.number+total.genotype.number))*100, 4), "%", sep = " ")
    
  }
  # dump unused objects
  bad.geno <- NULL

  # combine CHROM, LOCUS and POS into MARKERS
  vcf <- vcf %>%
    arrange(CHROM, LOCUS, POS) %>% 
    tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
  
  # keeping or not markers in common in all the populations---------------------
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
  
  # create a keeper list of MARKERS, CHROM, ID, POS, REF and ALT
  vcf.keeper <- vcf %>%
    select(MARKERS, REF, ALT) %>% 
    distinct(MARKERS)
  
  # last step before imputation ------------------------------------------------
  
  vcf <- vcf %>%
    select(-REF, -ALT) %>% 
    mutate(
      GT = stri_replace_all_fixed(GT, "./.", "NA", vectorize_all=F),
      GT = replace(GT, which(GT == "NA"), NA),
      POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end)
      ) %>%
    dcast(INDIVIDUALS + POP_ID ~ MARKERS, value.var = "GT")
  
  if (imputations == "max"){
    message("Calculating map-independent imputations using the most frequent allele.")
  } else {
    message("Calculating map-independent imputations using random forest")
  }
  
  # Imputations: VCF with imputed haplotypes using Random Forest ------------
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
      
      df.split.pop <- split(x = vcf, f = vcf$POP_ID) # slip data frame by population
      pop.list <- names(df.split.pop) # list the pop
      imputed.dataset <-list() # create empty list
      imputed.dataset <- foreach(i=pop.list, .packages = c("magrittr", "plyr", "dplyr", "tidyr", "stringi", "readr", "randomForestSRC", "reshape2")) %dopar% {
      # for (i in pop.list) {
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
      
    } else if (imputations.group == "global"){
      # Globally (not by pop_id)
      message("Imputations computed globally, take a break...")
      vcf <- plyr::colwise(factor, exclude = NA)(vcf)
      vcf.imp <- impute_markers_rf(vcf)
    } 
    
  } else if (imputations == "max") {
    
    if (missing(imputations.group) == "TRUE" | imputations.group == "populations"){
      message("Imputations computed by populations")
      
      vcf.imp <- suppressWarnings(
        vcf %>%
          tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
          group_by(MARKERS, POP_ID) %>% 
          mutate(
            GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
            GT = replace(GT, which(GT == "NA"), NA)
          ) %>%
          # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
          # will take the global observed values by markers for those cases.
          group_by(MARKERS) %>%
          mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)))
      )
      
    } else if (imputations.group == "global"){
      # Globally (not by pop_id)
      message("Imputations computed globally")
      
      vcf.imp <- suppressWarnings(
        vcf %>%
          tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
          group_by(MARKERS) %>% 
          mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)))
      )
    }
  }
  
  # transform the imputed dataset back into a VCF format------------------------
  
  vcf.imp <- suppressWarnings(
    vcf.imp %>% 
      tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
      arrange(POP_ID, INDIVIDUALS) %>% 
      group_by(MARKERS) %>% 
      mutate(INFO = stri_paste("NS=", n(), sep = "")) %>% 
      select(-POP_ID) %>% 
      dcast(MARKERS + INFO ~ INDIVIDUALS, value.var = "GT")
  )
  
  vcf.imp <- suppressWarnings(
    full_join(vcf.keeper, vcf.imp, by = "MARKERS") %>%
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
  
  # Write VCF header to the file ---------------------------------------------
  
  # filename
  if (missing(filename) == "TRUE"){
    filename <- "batch_1.vcf.imputed.vcf"
  } else {
    filename <- filename
  }
  
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
  
  #     # Info field 2
  #     info2 <- '##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">'
  #     info2 <- as.data.frame(info2)
  #     write.table(x = info2, file = filename, sep = " ", append = FALSE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Format field 1
  format1 <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  format1 <- as.data.frame(format1)
  write.table(x = format1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Write the imputed vcf to the file
  suppressWarnings(
    write.table(x = vcf.imp, file = filename, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE, quote = FALSE)
  )
  
  invisible(cat(sprintf(
    "VCF imputation
Total of %s genotypes 
Erased genotypes: %s (= %s )
Imputation of %s genotypes (%s)\n
Filename:
%s
Written in the directory:
%s",
    total.genotype.number, erased.genotype.number, genotype.percent, imputed.genotype.number, imputed.genotype.percent, filename, getwd()
  )))
}
