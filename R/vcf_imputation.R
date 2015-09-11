# VCF data imputation using Random Forest


#' @name vcf_imputation
#' @title VCF data imputation
#' @description Map-independent imputation of a VCF using Random Forest or the 
#' most frequent category.
#' @param vcf.file The VCF file created by STACKS.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param blacklist.genotypes A blacklist of genotypes 
#' containing 3 columns header 'LOCUS', 'POS' and 'INDIVIDUALS'.
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
                           pop.id.start, 
                           pop.id.end, 
                           blacklist.id = NULL,
                           whitelist.loci = NULL,
                           blacklist.genotypes = NULL,
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
  MARKER <- NULL
  everything <- NULL
  
  # Import/read VCF ------------------------------------------------------------- 
  message("Importing the VCF...")
  
  vcf <- read_delim(
    vcf.file, delim = "\t", 
    skip = 9,
    progress = interactive()
  ) %>%
    rename(CHROM = `#CHROM`)
  
  
  # VCF prep
  vcf <- vcf %>% 
    select(-c(QUAL, FILTER, INFO, FORMAT)) %>% 
    tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, ID, POS, REF, ALT)) %>%
    tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"),
             sep = ":", extra = "error") %>% 
    select(-c(READ_DEPTH, ALLELE_DEPTH, GL))
  
  
  # Whitelist-------------------------------------------------------------------
  if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "TRUE") {
    message("Using the whitelist from the directory")
    whitelist <- read_tsv(whitelist.loci, col_names = T) %>%
      rename(ID = LOCUS)
  } else if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "FALSE") {
    message("Using whitelist from your global environment")
    whitelist <- whitelist.loci %>%
      rename(ID = LOCUS)
  } else {
    message("No whitelist")
    whitelist <- NULL
  }
  
  # Blacklist id----------------------------------------------------------------
  if (missing(blacklist.id) == "FALSE" & is.vector(blacklist.id) == "TRUE") {
    message("Using the blacklisted id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)    
  } else if (missing(blacklist.id) == "FALSE" & is.vector(blacklist.id) == "FALSE") {
    message("Using the blacklisted id from your global environment")
    blacklist.id <- blacklist.id
  } else {
    message("No individual blacklisted")
    blacklist.id <- NULL
  }
  
  if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == TRUE) {
    # Combination 1: No whitelist and No blacklist -----------------------------
    vcf <- vcf
  } else if (is.null(whitelist.loci) == FALSE & is.null(blacklist.id) == TRUE) {
    
    # Combination 2: Using whitelist, but No blacklist -------------------------
    
    # just whitelist.loci, NO Blacklist of individual
    vcf <- vcf %>% 
      semi_join(whitelist, by = "ID") %>% 
      arrange(ID)
    
  } else if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == FALSE) {
    
    # Combination 3: Using a blacklist of id, but No whitelist -----------------
    # NO whitelist, JUST Blacklist of individual
    vcf <- vcf %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(ID)
    
  } else {
    # Combination 4: Using a whitelist and blacklist---------------------------
    
    # whitelist.loci + Blacklist of individual
    vcf <- vcf %>% 
      semi_join(whitelist, by = "ID") %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(ID)
  }
  
  # dump unused object
  whitelist <- NULL
  blacklist.id <- NULL
  
  # create a keeper list of MARKER, CHROM, ID, POS, REF and ALT
  vcf.keeper <- vcf %>%
    select(CHROM, ID, POS, REF, ALT) %>% 
    tidyr::unite(MARKER, CHROM, ID, POS, sep = "_", remove = TRUE) %>% 
    distinct(MARKER)
  
  
  # Erasing genotypes with poor coverage and/or genotype likelihood-------------
  if (is.null(blacklist.genotypes) == FALSE) {
    message("Erasing genotypes in the VCF using the list of blacklisted genotypes...")
    
    bad.geno <- read_tsv(blacklist.genotypes, col_names = T) %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      rename(ID = LOCUS) %>% 
      select(-POP_ID)
    
    # interesting stats.
    erased.genotype.number <- length(bad.geno$INDIVIDUALS)
    total.genotype.number <- length(vcf$GT[vcf$GT != "./."])
    genotype.percent <- paste(round(((erased.genotype.number/total.genotype.number)*100), 4), "%", sep = " ")
    
    vcf <- vcf %>%
      mutate(
        GT = ifelse(ID %in% bad.geno$ID & POS %in% bad.geno$POS & INDIVIDUALS %in% bad.geno$INDIVIDUALS, "./.", GT)
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
  
  # last step before imputation ------------------------------------------------
  
  vcf <- vcf %>%
    select(-REF, -ALT) %>% 
    tidyr::unite(MARKER, CHROM, ID, POS, sep = "_", remove = TRUE) %>% 
    mutate(GT = stri_replace_all_fixed(GT, "./.", "NA", vectorize_all=F)) %>% 
    mutate(POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end)) %>%
    dcast(INDIVIDUALS + POP_ID ~ MARKER, value.var = "GT")
  
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
      
      df.split.pop <- split(x = vcf, f = vcf$POP_ID) # slip data frame by population
      pop.list <- names(df.split.pop) # list the pop
      imputed.dataset <-list() # create empty list 
      for (i in pop.list) {
        sep.pop <- df.split.pop[[i]]
        sep.pop <- suppressWarnings(
          plyr::colwise(factor, exclude = "NA")(sep.pop)
        )
        imputed.dataset[[i]] <- impute_markers_rf(sep.pop)
        # message of progress for imputations by population
        pop.imputed <- paste("Completed imputations for pop ", i, sep = "")
        message(pop.imputed)
      }
      vcf.imp <- suppressWarnings(as.data.frame(bind_rows(imputed.dataset)))
      
      # dump unused objects
      df.split.pop <- NULL
      pop.list <- NULL
      sep.pop <- NULL
      imputed.dataset <- NULL
      
    } else if (imputations.group == "global"){
      # Globally (not by pop_id)
      message("Imputations computed globally, take a break...")
      vcf <- plyr::colwise(factor, exclude = "NA")(vcf)
      vcf.imp <- impute_markers_rf(vcf)
    } 
    
  } else if (imputations == "max") {
    
    if (missing(imputations.group) == "TRUE" | imputations.group == "populations"){
      message("Imputations computed by populations")
      
      vcf.imp <- vcf %>%
        tidyr::gather(MARKER, GT, -c(INDIVIDUALS, POP_ID)) %>%
        mutate(GT = replace(GT, which(GT=="NA"), NA)) %>%
        group_by(MARKER, POP_ID) %>% 
        mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)))
      
    } else if (imputations.group == "global"){
      # Globally (not by pop_id)
      message("Imputations computed globally")
      
      vcf.imp <- vcf %>%
        tidyr::gather(MARKER, GT, -c(INDIVIDUALS, POP_ID)) %>%
        mutate(GT = replace(GT, which(GT=="NA"), NA)) %>%
        group_by(MARKER) %>% 
        mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)))
    }
  }
  
  # transform the imputed dataset back into a VCF format------------------------
  
  vcf.imp <- suppressWarnings(
    vcf.imp %>% 
      tidyr::gather(MARKER, GT, -c(INDIVIDUALS, POP_ID)) %>%
      arrange(POP_ID, INDIVIDUALS) %>% 
      group_by(MARKER) %>% 
      mutate(INFO = stri_paste("NS=", n(), sep = "")) %>% 
      select(-POP_ID) %>% 
      dcast(MARKER + INFO ~ INDIVIDUALS, value.var = "GT")
  )
  
  vcf.imp <- suppressWarnings(
    full_join(vcf.keeper, vcf.imp, by = "MARKER") %>%
      tidyr::separate(MARKER, c("CHROM", "ID", "POS"), sep = "_", extra = "error") %>%
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
