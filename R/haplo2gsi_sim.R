# Write a gsi_sim file from STACKS haplotype file
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID",
                                                        "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt",
                                                        "HAPLOTYPES", "SAMPLES", "ALLELE", "ALLELE1", "ALLELE2", "GENOTYPE",
                                                        "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX",
                                                        "colwise", "detectCores", "mc.cores", "."))

#' @name haplo2gsi_sim
#' @title Use the batch_x.haplotypes.tsv file to write a gsi_sim input file.
#' @description \code{gsi_sim} is a tool for doing and simulating genetic stock identification. 
#' The \code{haplo2gsi_sim} function can first filter the haplotypes file 
#' with a whitelist of loci and a blacklist of individuals. 
#' Map-independent imputation using Random Forest or the most frequent category
#' is also available as an option.
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param erase.paralogs (logical) \code{FALSE} Loci with more than 2 alleles 
#' (paralogs and/or sequencing errors) will be removed. \code{TRUE} genotypes 
#' with more than 2 alleles (paralogs and/or sequencing errors) will be erased.
#' @param sample.markers (number) Should the output contain a random 
#' subsample of your markers. Default= \code{"all"}.
#' e.g. \code{sample.markers = 500} to use only 500 randomly chosen markers.
#' @param iterations The number of iterations for marker resampling.
#' Default is 10.
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param gsi_sim.filename The name of the file written to the directory.
#' Use the extension ".txt" at the end. Default \code{gsi_sim_data.txt}.
#' The number of markers used will be appended to the name of the file.
#' @param pop.levels A character string with your populations ordered.
#' @param pop.labels An optional character string with new populations names.
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
#' @return When no imputation is selected a gsi_sim file is saved to the 
#' working directory. When imputation is selected 2 gsi_sim files are saved to
#' the working directory.
#' @export
#' @rdname haplo2gsi_sim
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
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


haplo2gsi_sim <- function(haplotypes.file, 
                          whitelist.loci = NULL,
                          erase.paralogs = NULL,
                          sample.markers = "all",
                          iterations = 10,
                          blacklist.id = NULL, 
                          gsi_sim.filename = "gsi_sim_data.txt",
                          pop.levels, pop.labels, pop.id.start, pop.id.end,
                          baseline, mixture,
                          imputations = FALSE,
                          imputations.group = "populations",
                          num.tree = 100,
                          iteration.rf = 10,
                          split.number = 100,
                          verbose = FALSE,
                          parallel.core = 2) {
  
  #Create a folder based on filename to save the output files
  directory <- stri_paste(getwd(),"/",imputations,"_",imputations.group,"/", sep = "")
  dir.create(file.path(directory))
  
  if (imputations == "FALSE") {
    message("haplo2gsi_sim: without imputation...")
  } else {
    message("haplo2gsi_sim: with imputations...")
  }
  
  # To work inside foreach ...
  if(missing(mixture)){
    mixture <- FALSE
  } else{
    mixture <- mixture
  }
  
  if(missing(baseline)){
    baseline <- FALSE
  } else{
    baseline <- baseline
  }
  
  # Haplotype file--------------------------------------------------------------
  haplotype <- read_tsv(file = haplotypes.file, col_types = cols(.default = col_character()), col_names = T, na = "-") %>% 
    select(-Cnt) %>% 
    rename(Catalog.ID = `Catalog ID`) %>%
    tidyr::gather(INDIVIDUALS, HAPLOTYPES, -Catalog.ID) %>% 
    mutate(
      HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = "-"),
      Catalog.ID = as.integer(Catalog.ID)
    )
  
  
  # Whitelist-------------------------------------------------------------------
  if (is.null(whitelist.loci) | missing(whitelist.loci)) {
    message("No whitelist")
    haplotype.whitelist <- haplotype
    
  } else if (is.vector(whitelist.loci)) {
    message("Using the whitelist from the directory")
    whitelist <- read_tsv(whitelist.loci, col_names = T) %>%
      rename(Catalog.ID = LOCUS)
    
    haplotype.whitelist <- haplotype %>% 
      semi_join(whitelist, by = "Catalog.ID") %>% 
      arrange(Catalog.ID)
    
  } else {
    message("Using whitelist from your global environment")
    whitelist <- whitelist.loci %>%
      rename(Catalog.ID = LOCUS)
    
    haplotype.whitelist <- haplotype %>% 
      semi_join(whitelist, by = "Catalog.ID") %>% 
      arrange(Catalog.ID)
  }
  
  # Blacklist-------------------------------------------------------------------
  if (is.null(blacklist.id) | missing(blacklist.id)) {
    message("No individual blacklisted")
    haplotype.whitelist.blacklist.id <- haplotype.whitelist
    
  } else if (is.vector(blacklist.id)) {
    message("Using the blacklisted id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)
    
    haplotype.whitelist.blacklist.id <- haplotype.whitelist %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(Catalog.ID)
    
  } else {
    message("Using the blacklisted id from your global environment")
    
    haplotype.whitelist.blacklist.id <- haplotype.whitelist %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(Catalog.ID)
  }
  
  # dump unused object
  haplotype <- NULL
  whitelist <- NULL
  blacklist.id <- NULL
  haplotype.whitelist <- NULL
  
  # Paralogs-------------------------------------------------------------------
  message("Looking for paralogs...")
  
  if(missing(erase.paralogs) | erase.paralogs == FALSE){
    message("Loci with more than 2 alleles will be removed")
    
    paralogs <- haplotype.whitelist.blacklist.id %>%
      mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
      group_by(Catalog.ID) %>%
      summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
      filter(POLYMORPHISM_MAX > 1) %>%
      group_by(Catalog.ID) %>%
      select(Catalog.ID) %>%
      distinct(Catalog.ID)
    
    haplotype.whitelist.blacklist.id.paralogs <- suppressWarnings(
      haplotype.whitelist.blacklist.id %>%
        anti_join(paralogs, by = "Catalog.ID") %>%
        filter(HAPLOTYPES != "consensus")
    )
    message(stri_join("Found and/or removed", n_distinct(paralogs$Catalog.ID), "paralogs", sep = " "))
    
  } else {
    message("Erasing genotypes with more than 2 alleles")
    
    # get the number of genotypes...
    haplo.number <- haplotype.whitelist.blacklist.id %>%
      filter(HAPLOTYPES != "-") %>%
      select(HAPLOTYPES)
    
    haplo.poly <- haplotype.whitelist.blacklist.id %>% 
      mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/"))
    
    erased.genotype.number <- length(haplo.poly$INDIVIDUALS[haplo.poly$POLYMORPHISM > 1])
    total.genotype.number.haplo <- length(haplo.number$HAPLOTYPES)
    percent.haplo <- paste(round(((erased.genotype.number/total.genotype.number.haplo)*100), 2), "%", sep = " ")
    message(stri_paste("Out of a total of ", total.genotype.number.haplo, " genotypes, ", percent.haplo, " (", erased.genotype.number, ")"," will be erased"))
    
    message("Erasing... Erasing...")
    
    # Erasing genotype with the blacklist
    
    haplotype.whitelist.blacklist.id.paralogs <- suppressWarnings(
      haplo.poly %>%
        filter(HAPLOTYPES != "consensus") %>% 
        mutate(HAPLOTYPES = ifelse(POLYMORPHISM > 1, "-", HAPLOTYPES))
    )
  }
  # dump unused object
  paralogs <- NULL
  haplo.number <- NULL
  haplo.poly <- NULL
  erased.genotype.number <- NULL
  total.genotype.number.haplo <- NULL
  percent.haplo <- NULL
  haplotype.whitelist.blacklist.id <- NULL
  
  # Conversion into gsi_sim -----------------------------------------------------
  
  message("Haplotypes into gsi_sim factory ...")
  
  # Haplo prep
  # Remove paralogs
  # Remove consensus loci
  
  if(missing(pop.labels)){
    pop.labels <- pop.levels
  } else {
    pop.labels <- pop.labels
  }
  
  # This data set will be used with and without imputations
  haplo <- suppressWarnings(
    haplotype.whitelist.blacklist.id.paralogs %>%
      mutate(
        HAPLOTYPES = stri_replace_all_fixed(HAPLOTYPES, "-", "NA", 
                                            vectorize_all=F),
        HAPLOTYPES = replace(HAPLOTYPES, which(HAPLOTYPES == "NA"), NA),
        POP_ID = factor(substr(INDIVIDUALS, pop.id.start, pop.id.end), 
                        levels = pop.levels, labels = pop.labels, ordered = T),
        POP_ID = droplevels(POP_ID)
      ) %>% 
      arrange(Catalog.ID) %>%
      dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES") %>% 
      arrange(POP_ID, INDIVIDUALS)
  )
  
  # dump unused objects
  haplotype.whitelist.blacklist.id.paralogs <- NULL
  
  # No imputation --------------------------------------------------------------
  # Further work on the data
  haplo.no.imputation <- suppressWarnings(
    haplo %>% 
      melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
        sep = "/", extra = "drop", remove = T
      ) %>%
      mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>% 
      melt(
        id.vars = c("Catalog.ID", "INDIVIDUALS", "POP_ID"),
        measure.vars = c("ALLELE1", "ALLELE2"), 
        variable.name = "ALLELE", 
        value.name = "NUCLEOTIDES"
      ) %>% 
      dcast(INDIVIDUALS + POP_ID + ALLELE ~ Catalog.ID, value.var = "NUCLEOTIDES")
  )
  message("step 1/3: completed")
  
  haplo.no.imputation <- suppressWarnings(
    haplo.no.imputation %>% 
      plyr::colwise(factor, exclude = NA)(.)
  )
  
  message("step 2/3: completed")
  
  haplo.no.imputation <- suppressWarnings(
    haplo.no.imputation %>%
      mutate(GROUP = rep(1, times = nrow(.))) %>% 
      group_by(GROUP) %>% 
      mutate_each(funs(as.integer), -c(ALLELE, INDIVIDUALS, POP_ID, GROUP)) %>%
      ungroup() %>% 
      select(-GROUP) %>% 
      melt(id.vars = c("INDIVIDUALS", "POP_ID", "ALLELE"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>%
      mutate(HAPLOTYPES = stri_replace_na(str = HAPLOTYPES, replacement = "0")) %>% 
      mutate(POP_ID = droplevels(POP_ID)))
  
  message("step 3/3: completed")
  
  
  
  # Imputations ----------------------------------------------------------------
  if (imputations == "max"){
    message("Calculating map-independent imputations using the most frequent allele.")
  } else if (imputations == "rf"){
    message("Calculating map-independent imputations using random forest")
  }
  
  # Imputations: gsi_sim with imputed haplotypes using Random Forest or the most frequent allele-----------
  if (imputations != "FALSE"){
    
    if (imputations == "rf") {
      # A different format is required for the imputations 
      # Transformed columns into factor excluding the "NA"
      haplo.imputation <- suppressWarnings(
        plyr::colwise(factor, exclude = NA)(haplo)
      )
      
      # Parallel computations options
      if (missing(parallel.core) == "TRUE"){
        # Automatically select all the core -1 
        options(rf.cores=detectCores()-1, mc.cores=detectCores()-1)
      } else {
        options(rf.cores = parallel.core, mc.cores = parallel.core)
      }
      
      # imputations function using Random Forest with the package randomForestSRC
      impute_markers_rf <- function(x){
        randomForestSRC::impute.rfsrc(data = x, 
                                      ntree = num.tree, 
                                      nodesize = 1, 
                                      nsplit = split.number, 
                                      nimpute = iteration.rf, 
                                      do.trace = verbose)
      }
      
      # imputations by populations (default) or globally -------------------------
      if (missing(imputations.group) == "TRUE" | imputations.group == "populations"){
        message("Imputations computed by populations, take a break...")
        
        # By pop using for loop to imputed with message when completed
        df.split.pop <- split(x = haplo.imputation, f = haplo.imputation$POP_ID) # slip data frame by population
        pop.list <- names(df.split.pop) # list the pop
        imputed.dataset <-list() # create empty list
        
        # loop trough the pop
        for (i in pop.list) {
          sep.pop <- df.split.pop[[i]]
          imputed.dataset[[i]] <- impute_markers_rf(sep.pop)
          
          # message of progress for imputations by population
          message(paste("Completed imputations for pop ", i, sep = ""))
        }
        # bind rows and transform into a data frame
        haplo.imputation <- as.data.frame(bind_rows(imputed.dataset))
        
        # remove introduced NA if some pop don't have the markers by using
        # RF globally
        haplo.imputation <- impute_markers_rf(haplo.imputation)
        
        # dump unused objects
        df.split.pop <- NULL
        pop.list <- NULL
        imputed.dataset <- NULL
        sep.pop <- NULL
        
      } else if (imputations.group == "global")
        # RF Globally
        message("Imputations computed globally, take a break...")
      haplo.imputation <- impute_markers_rf(haplo.imputation)
      
    } else if (imputations == "max") { # most frequent allele by pop
      
      if (missing(imputations.group) == "TRUE" | imputations.group == "populations"){
        message("Imputations computed by populations")
        
        haplo.imputation <- suppressWarnings(
          haplo %>%
            melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
            mutate(HAPLOTYPES = replace(HAPLOTYPES, which(HAPLOTYPES == "NA"), NA)) %>%
            group_by(Catalog.ID, POP_ID) %>%
            mutate(
              HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = max(HAPLOTYPES, na.rm = TRUE)),
              HAPLOTYPES = replace(HAPLOTYPES, which(HAPLOTYPES == "NA"), NA)
            ) %>%
            # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
            # will take the global observed values by markers for those cases.
            group_by(Catalog.ID) %>%
            mutate(HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = max(HAPLOTYPES, na.rm = TRUE))) %>% 
            dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES") %>% 
            arrange(POP_ID, INDIVIDUALS)
        )
        
      } else if (imputations.group == "global"){ # global
        # Globally (not by pop_id)
        message("Imputations computed globally")
        
        haplo.imputation <- haplo %>%
          melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
          mutate(HAPLOTYPES = replace(HAPLOTYPES, which(HAPLOTYPES == "NA"), NA)) %>%
          group_by(Catalog.ID) %>% 
          mutate(HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = max(HAPLOTYPES, na.rm = TRUE))) %>% 
          dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES") %>% 
          arrange(POP_ID, INDIVIDUALS)
      }
    }
    
    # transform the imputed dataset into gsi_sim object ------------------------
    message("Imputed haplotypes into gsi_sim factory ...")
    
    haplo.imputation <- suppressWarnings(
      haplo.imputation %>% 
        melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>%
        mutate(Catalog.ID = as.integer(as.character(Catalog.ID))) %>% 
        arrange(Catalog.ID) %>% 
        tidyr::separate(
          col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
          sep = "/", extra = "drop", remove = TRUE
        ) %>%
        mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>% 
        melt(
          id.vars = c("Catalog.ID", "INDIVIDUALS", "POP_ID"),
          measure.vars = c("ALLELE1", "ALLELE2"), 
          variable.name = "ALLELE", 
          value.name = "NUCLEOTIDES"
        ) %>% 
        dcast(INDIVIDUALS + POP_ID + ALLELE ~ Catalog.ID, value.var = "NUCLEOTIDES")
    )
    message("step 1/3: completed")
    
    haplo.imputation <- suppressWarnings(
      haplo.imputation %>% 
        plyr::colwise(factor, exclude = NA)(.)
    )
    message("step 2/3: completed")
    
    haplo.imputation <- suppressWarnings(
      haplo.imputation %>%
        mutate(GROUP = rep(1, times = nrow(.))) %>% 
        group_by(GROUP) %>% 
        mutate_each(funs(as.integer), -c(ALLELE, INDIVIDUALS, POP_ID, GROUP)) %>%
        ungroup() %>% 
        select(-GROUP) %>% 
        melt(id.vars = c("INDIVIDUALS", "POP_ID", "ALLELE"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>%
        mutate(HAPLOTYPES = stri_replace_na(str = HAPLOTYPES, replacement = "0")) %>% 
        mutate(POP_ID = droplevels(POP_ID)))
    
    message("step 3/3: completed")
  }  
  
  # Random sampling of markers -------------------------------------------------
  
  # get the unique list of markers for sample.markers that might have "all"
  # in it, and we need to change this to the maximum number of markers 
  unique.markers <- haplo.no.imputation %>% select(Catalog.ID) %>% distinct(Catalog.ID) %>% arrange(Catalog.ID)
  sample.markers <- stri_replace_all_fixed(str = sample.markers, pattern = "all", replacement = nrow(unique.markers), vectorize_all = TRUE)
  
  # create a random seed number file
  filename.seed <- "random_seed_numbers_gsi_sim.tsv"
  random.seed.table <- data_frame(MARKER_NUMBER = integer(0), ITERATIONS = integer(0), RANDOM_NUMBER = integer(0))
  write_tsv(x = random.seed.table, path = paste0(directory,filename.seed), col_names = TRUE, append = FALSE)
  
  for (i in sample.markers) {
    
    message(paste("Marker number: ", i))
    # i <- 500
    i <- as.numeric(i)
    
    # Start cluster registration backend using n - 1 CPU
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
    
    # Number of times to repeat the sampling of markers
    if (missing(iterations)){
      iterations <- 10
    } else {
      iterations <- iterations
    }
    
    j <- NULL
    foreach(j=1:iterations, .packages = c("magrittr", "plyr", "dplyr", "tidyr", "stringi", "readr", "randomForestSRC", "reshape2")) %dopar% {
      # for(j in 1:iterations) {
      # j <-100
      # Subsampling markers --------------------------------------------------------
      # Set seed for random sampling
      random.seed <- sample(x = 1:10000, size = 1)
      random.seed.table <- data_frame(MARKER_NUMBER = i, ITERATIONS = j, RANDOM_NUMBER = random.seed)
      write_tsv(x = random.seed.table, path = paste0(directory, filename.seed), col_names = FALSE, append = TRUE)
      
      # sample markers
      random.markers <- sample_n(tbl = unique.markers, size = i, replace = FALSE) %>% 
        arrange(Catalog.ID)
      
      # No imputation ----------------------------------------------------------
      haplo.no.imputation.random <- haplo.no.imputation %>% 
        semi_join(random.markers, by = "Catalog.ID") %>% 
        arrange(Catalog.ID) %>% 
        dcast(INDIVIDUALS + POP_ID ~ Catalog.ID + ALLELE, value.var = "HAPLOTYPES") %>% 
        arrange(POP_ID)
      
      # get the list of loci after filter  
      loci <- unique(random.markers$Catalog.ID)
      n.markers <- length(loci)
      
      if(mixture == FALSE){
        message("No baseline or mixture data")
        haplo.no.imputation.random <- haplo.no.imputation.random
      } else {
        baseline.data <- suppressWarnings(
          haplo.no.imputation.random %>%
            filter(POP_ID %in% baseline) %>% 
            arrange(POP_ID) %>% 
            mutate(POP_ID = droplevels(POP_ID))
        )
        
        mixture.data <- suppressWarnings(
          haplo.no.imputation.random %>%
            filter(POP_ID %in% mixture) %>% 
            arrange(POP_ID) %>% 
            mutate(POP_ID = droplevels(POP_ID))
        )
      }
      # results no imputation-------------------------------------------------------
      if(mixture == FALSE){
        n.individuals <- n_distinct(haplo.no.imputation.random$INDIVIDUALS)
        
        # Create a vector with the population ordered by levels
        pop <- haplo.no.imputation.random$POP_ID
        haplo.no.imputation.random <- suppressWarnings(haplo.no.imputation.random %>% select(-POP_ID))
        
        # split gsi_sim by populations
        gsi_sim.split <- split(haplo.no.imputation.random, pop)
        
        # Write the file in gsi_sim format
        message("Output...No imputation")
        
        # gsi_sim filename
        if (missing(gsi_sim.filename) == "TRUE"){
          filename <- stri_paste("gsi_sim.data", i, j, "txt", sep = ".")
        } else {
          filename <- gsi_sim.filename
          sample.markers.in.filename <- stri_paste(i, j,"txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = sample.markers.in.filename)
        }
        
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
        if (missing(gsi_sim.filename) == "TRUE"){
          filename <- stri_paste("gsi_sim.baseline.data", i, j, "txt", sep = ".")
        } else {
          filename <- gsi_sim.filename
          sample.markers.in.filename <- stri_paste("baseline", i, j, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = sample.markers.in.filename)
        }
        
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
        if (missing(gsi_sim.filename) == "TRUE"){
          filename <- stri_paste("gsi_sim.mixture.data", i, j, "txt", sep = ".")
        } else {
          filename <- gsi_sim.filename
          sample.markers.in.filename <- stri_paste("mixture", i, j, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = sample.markers.in.filename)
        }
        
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
      haplo.imputation.random <- haplo.imputation %>% 
        semi_join(random.markers, by = "Catalog.ID") %>% 
        arrange(Catalog.ID) %>% 
        dcast(INDIVIDUALS + POP_ID ~ Catalog.ID + ALLELE, value.var = "HAPLOTYPES") %>% 
        arrange(POP_ID)
      
      if(mixture == FALSE){
        message("No baseline or mixture data")
        haplo.imputation.random <- haplo.imputation.random
      } else {
        # Baseline dataset
        baseline.data <- suppressWarnings(
          haplo.imputation.random %>%
            filter(POP_ID %in% baseline) %>% 
            arrange(POP_ID) %>% 
            mutate(POP_ID = droplevels(POP_ID))
        )
        # Mixture dataset
        mixture.data <- suppressWarnings(
          haplo.imputation.random %>%
            filter(POP_ID %in% mixture) %>% 
            arrange(POP_ID) %>% 
            mutate(POP_ID = droplevels(POP_ID))
        )
      }
      # results imputation----------------------------------------------------
      if(mixture == FALSE){
        n.individuals <- n_distinct(haplo.imputation.random$INDIVIDUALS)
        
        # Create a vector with the population ordered by levels
        pop <- haplo.imputation.random$POP_ID
        haplo.imputation.random <- suppressWarnings(haplo.imputation.random %>% select(-POP_ID))
        
        # split gsi_sim by populations
        gsi_sim.split <- split(haplo.imputation.random, pop)
        
        # Write the file in gsi_sim format 
        message("Output of imputed data...")
        
        # gsi_sim filename
        if (missing(gsi_sim.filename) == "TRUE"){
          filename <- stri_paste("gsi_sim.data", i, j, "txt", sep = ".")
        } else {
          filename <- gsi_sim.filename
          sample.markers.in.filename <- stri_paste(i, j, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = sample.markers.in.filename)
        }
        
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
        if (missing(gsi_sim.filename) == "TRUE"){
          filename <- stri_paste("gsi_sim.baseline.data", i, j, "txt", sep = ".")
        } else {
          filename <- gsi_sim.filename
          sample.markers.in.filename <- stri_paste("baseline", i, j, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = sample.markers.in.filename)
        }
        
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
        if (missing(gsi_sim.filename) == "TRUE"){
          filename <- stri_paste("gsi_sim.mixture.data", i,"txt", sep = ".")
        } else {
          filename <- gsi_sim.filename
          sample.markers.in.filename <- stri_paste("mixture", i, "txt", sep = ".")
          filename <- stri_replace_all_fixed(filename, pattern = "txt",
                                             replacement = sample.markers.in.filename)
        }
        
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
    }
    # close parallel connection settings
    stopCluster(cl)
  }
}
