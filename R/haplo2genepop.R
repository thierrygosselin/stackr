# Write a genepop file from STACKS haplotype file
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE1", "ALLELE2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX", "colwise", "detectCores", "mc.cores"))

#' @name haplo2genepop
#' @title Use the batch_x.haplotypes.tsv file to write a genpop file.
#' @description This function can first filter the haplotypes file 
#' with a whitelist of loci and a blacklist of individuals. 
#' Map-independent imputation using Random Forest 
#' is also available as an option.
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param genepop.filename The name of the file written to the directory.
#' Use the extension ".gen" at the end. Default \code{genepop.gen}.
#' @param genepop.header The first line of the Genepop file.
#' Default \code{my firt genepop}.
#' @param pop.levels An optional character string with your populations ordered.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
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
#' @details The imputations using Random Forest requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set 
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals 
#' will require 15 min.
#' @return When no imputation is selected a genepop file is saved to the 
#' working directory. When imputation is selected 2 genepop files are saved to
#' the working directory.
#' @export
#' @rdname haplo2genepop
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_pad
#' @references Rousset F (2008) genepop’007: a complete re-implementation
#' of the genepop software for Windows and Linux.
#' Molecular Ecology Resources, 8, 103-106.
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


haplo2genepop <- function(haplotypes.file, 
                          whitelist.loci = NULL, 
                          blacklist.id = NULL, 
                          genepop.filename = "genepop.gen",
                          genepop.header = "my first genepop", 
                          pop.levels, pop.id.start, pop.id.end,
                          imputations = FALSE,
                          imputations.group = "populations",
                          num.tree = 100,
                          iteration.rf = 10,
                          split.number = 100,
                          verbose = FALSE,
                          parallel.core = 2
) {  
  
  
  if (imputations == "FALSE") {
    message("haplo2genind: without imputation...")
  } else {
    message("haplo2genind: with imputations...")
  }
  
  
  # Haplotype file--------------------------------------------------------------
  haplotype <- read_tsv(file = "batch_1.haplotypes.tsv", col_names = T) %>%
    select(-Cnt) %>% 
    rename(Catalog.ID = `Catalog ID`) %>%
    melt(id.vars = "Catalog.ID", variable.name = "INDIVIDUALS", value.name = "HAPLOTYPES")
  # gather(INDIVIDUALS, HAPLOTYPES, -Catalog.ID)
  
  
  # Whitelist-------------------------------------------------------------------
  if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "TRUE") {
    message("Using the whitelist from the directory")
    whitelist <- read_tsv(whitelist.loci, col_names = T) %>%
      rename(Catalog.ID = LOCUS)
  } else if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "FALSE") {
    message("Using whitelist from your global environment")
    whitelist <- whitelist.loci %>%
      rename(Catalog.ID = LOCUS)
  } else {
    message("No whitelist")
    whitelist <- NULL
  }
  
  
  # Blacklist-------------------------------------------------------------------
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
    haplotype <- haplotype    
    
  } else if (is.null(whitelist.loci) == FALSE & is.null(blacklist.id) == TRUE) {
    
    # Combination 2: Using whitelist, but No blacklist -------------------------
    
    # just whitelist.loci, NO Blacklist of individual
    haplotype <- haplotype %>% 
      group_by(Catalog.ID) %>% 
      semi_join(whitelist, by = "Catalog.ID") %>% 
      arrange(Catalog.ID)
    
  } else if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == FALSE) {
    
    # Combination 3: Using a blacklist of id, but No whitelist -----------------
    
    # NO whitelist, JUST Blacklist of individual
    haplotype <- haplotype %>%
      group_by(Catalog.ID) %>% 
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(Catalog.ID)
    
  } else {
    # Combination 4: Using a whitelist and blacklist---------------------------
    
    # whitelist.loci + Blacklist of individual
    haplotype <- haplotype %>%
      group_by(Catalog.ID) %>% 
      semi_join(whitelist, by = "Catalog.ID") %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(Catalog.ID)
  }
  
  # dump unused object
  whitelist <- NULL
  blacklist.id <- NULL
  
  # Paralogs-------------------------------------------------------------------
  message("Looking for paralogs...")
  
  paralogs <- haplotype %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(Catalog.ID) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    group_by(Catalog.ID) %>%
    select(Catalog.ID) %>%
    distinct(Catalog.ID)
  
  nparalogs <- stri_join("Found and/or removed", n_distinct(paralogs$Catalog.ID), "paralogs", sep = " ")
  message(nparalogs)
  
  # Conversion into genepop -----------------------------------------------------
  
  message("Haplotypes into genepop factory ...")
  
  # Haplo prep
  # Remove paralogs
  # Remove consensus loci
  
  haplo.filtered <- suppressWarnings(
    haplotype %>%
      group_by(Catalog.ID) %>% 
      anti_join(paralogs, by = "Catalog.ID") %>%
      filter(HAPLOTYPES != "consensus") %>%    
      mutate(
        HAPLOTYPES = stri_replace_all_fixed(HAPLOTYPES, "-", "NA", 
                                            vectorize_all=F),
        POP_ID = factor(substr(INDIVIDUALS, pop.id.start, pop.id.end), 
                        levels = pop.levels, ordered = T)
      ) %>% 
      arrange(Catalog.ID)
  )
  # get the list of loci after filter  
  loci <- unique(haplo.filtered$Catalog.ID)
  
  haplo.filtered <- haplo.filtered %>%
    dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES")
  
  message("step 1/5: completed")
  
  # dump unused objects
  paralogs <- NULL
  haplotype <- NULL
  
  
  # No imputation --------------------------------------------------------------
  # Further work on the data
  haplo.prep <- suppressWarnings(
    haplo.filtered %>% 
      # gather(Catalog.ID, HAPLOTYPES, -c(INDIVIDUALS, POP_ID)) %>% 
      melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
      separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
        sep = "/", extra = "drop", remove = T
      ) %>%
      mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2))
  )
  
  message("step 2/5: completed")
  
  haplo.prep <- haplo.prep %>% 
    melt(
      id.vars = c("Catalog.ID", "INDIVIDUALS", "POP_ID"),
      measure.vars = c("ALLELE1", "ALLELE2"), 
      variable.name = "ALLELE", 
      value.name = "NUCLEOTIDES"
    ) %>% 
    dcast(INDIVIDUALS + POP_ID + ALLELE ~ Catalog.ID, value.var = "NUCLEOTIDES")
  
  message("step 3/5: completed")
  
  haplo.prep <- suppressWarnings(
    haplo.prep %>% 
      colwise(factor, exclude = "NA")(.)
  )
  message("step 4/5: completed")
  
  # This part is different than gtypes and genepop...
  haplo.prep <- suppressWarnings(
    haplo.prep %>%
      mutate(GROUP = rep(1, times = nrow(.))) %>% 
      group_by(GROUP) %>% 
      mutate_each(funs(as.integer), -c(ALLELE, INDIVIDUALS, POP_ID, GROUP)) %>%
      ungroup() %>% 
      select(-GROUP, -POP_ID) %>% 
      melt(id.vars = c("INDIVIDUALS", "ALLELE"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>%
      mutate(HAPLOTYPES = as.character(HAPLOTYPES)) %>% 
      mutate(HAPLOTYPES = stri_pad_left(str = HAPLOTYPES, width = 3, pad = "0")) %>% 
      mutate(HAPLOTYPES = stri_replace_na(str = HAPLOTYPES, replacement = "000")) %>% 
      dcast(Catalog.ID + INDIVIDUALS ~ ALLELE, value.var = "HAPLOTYPES") %>%
      unite(GENOTYPE, ALLELE1:ALLELE2, sep = "") %>%
      arrange(Catalog.ID) %>% 
      dcast(INDIVIDUALS ~ Catalog.ID, value.var = "GENOTYPE") %>%
      mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))
  )
  
  message("step 5/5: completed")
  
  # results no imputation-------------------------------------------------------
  # convert to genepop     
  
  # Create a vector with the population ordered by levels
  if (missing(pop.levels) == "TRUE") {
    pop <- substr(haplo.prep$INDIVIDUALS, pop.id.start, pop.id.end)
  } else {
    pop <- factor(substr(haplo.prep$INDIVIDUALS, 
                         pop.id.start, pop.id.end),
                  levels = pop.levels, ordered = T)
  }
  
  # split genepop by populations
  genepop.split <- split(haplo.prep, pop)
  
  # Write the file in genepop format
  message("Output...No imputation")
  
  # genepop filename
  if (missing(genepop.filename) == "TRUE"){
    genepop.filename <- "genepop.gen"
  } else {
    genepop.filename <- genepop.filename
  }
  
  # genepop header
  if (missing(genepop.header) == "TRUE"){
    genepop.header <- "my first genepop"
  } else {
    genepop.header <- genepop.header
  }
  
  file.no.imputation <- file(genepop.filename, "write")
  cat(genepop.header, "\n", file = file.no.imputation, append = TRUE)
  cat(loci, sep = "\n", file = file.no.imputation, append = TRUE)
  
  for (i in 1:length(genepop.split)) {
    cat("pop\n", file = file.no.imputation, append = TRUE)
    write.table(
      genepop.split[[i]],
      file = file.no.imputation, 
      append = TRUE, 
      col.names = FALSE,
      row.names = FALSE, 
      sep = " ", 
      quote = FALSE
    )
  }
  close(file.no.imputation)
  
  invisible(cat(sprintf("Genepop (no imputation) file name:\n%s\n\nWritten in the directory:\n%s", 
                        genepop.filename, getwd())))
  
  if (imputations == "max"){
    message("Calculating map-independent imputations using the most frequent allele.")
  } else if (imputations == "rf"){
    message("Calculating map-independent imputations using random forest")
  }
  
  # dump unused objects
  haplo.prep <- NULL
  
  # Imputations: genepop with imputed haplotypes using Random Forest ------------------
  if (imputations != "FALSE"){
    
    if (imputations == "rf") {
      # A different format is required for the imputations 
      # Transformed columns into factor excluding the "NA"
      haplo.prep <- suppressWarnings(
        plyr::colwise(factor, exclude = "NA")(haplo.filtered)
      )
      
      # Parallel computations options
      if (missing(parallel.core) == "TRUE"){
        # Automatically select all the core -1 
        options(rf.cores=detectCores()-1, mc.cores=detectCores()-1)
      } else {
        options(rf.cores = parallel.core, mc.cores = mc.cores)
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
        
        # By pop using dplyr
        #     INDIVIDUALS <- imp.prep %>% select(INDIVIDUALS)
        #     imp.test <- bind_cols(INDIVIDUALS,
        #                           imp.prep %>%
        #                             select(-INDIVIDUALS) %>% 
        #                             group_by(POP_ID) %>% 
        #                             colwise(factor, exclude = "NA")(.)%>% 
        #                             do(impute_markers_rf(.))
        #     )
        
        # By pop using for loop to imputed with message when completed
        df.split.pop <- split(x = haplo.prep, f = haplo.prep$POP_ID) # slip data frame by population
        pop.list <- names(df.split.pop) # list the pop
        imputed.dataset <-list() # create empty list 
        for (i in pop.list) {
          sep.pop <- df.split.pop[[i]]
          imputed.dataset[[i]] <- impute_markers_rf(sep.pop)
          # message of progress for imputations by population
          pop.imputed <- paste("Completed imputations for pop ", i, sep = "")
          message(pop.imputed)
        }
        haplo.imp <- as.data.frame(bind_rows(imputed.dataset))
        
        # dump unused objects
        haplo.filtered <- NULL
        haplo.prep <- NULL
        df.split.pop <- NULL
        pop.list <- NULL
        sep.pop <- NULL
        imputed.dataset <- NULL
        
      } else if (imputations.group == "global"){
        # Globally (not by pop_id)
        message("Imputations computed globally, take a break...")
        
        #         INDIVIDUALS <- haplo.prep %>% select(INDIVIDUALS)
        #         haplo.prep <- suppressWarnings(
        #           haplo.prep %>%
        #             select(-INDIVIDUALS) %>% 
        #             colwise(factor, exclude = "NA")(.)%>% 
        #             do(impute_markers_rf(.))
        #         )
        #         
        #         haplo.imp <- bind_cols(INDIVIDUALS, haplo.prep)
        
        
        haplo.imp <- impute_markers_rf(haplo.prep)
        
        # dump unused objects
        haplo.prep <- NULL
      } 
      
    } else if (imputations == "max") {
      
      if (missing(imputations.group) == "TRUE" | imputations.group == "populations"){
        message("Imputations computed by populations")
        
        haplo.imp <- haplo.filtered %>%
          melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
          mutate(HAPLOTYPES = replace(HAPLOTYPES, which(HAPLOTYPES=="NA"), NA)) %>%
          group_by(Catalog.ID, POP_ID) %>% 
          mutate(HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = max(HAPLOTYPES, na.rm = TRUE))) %>% 
          dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES")
        
        
      } else if (imputations.group == "global"){
        # Globally (not by pop_id)
        message("Imputations computed globally")
        
        haplo.imp <- haplo.filtered %>%
          melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
          mutate(HAPLOTYPES = replace(HAPLOTYPES, which(HAPLOTYPES=="NA"), NA)) %>%
          group_by(Catalog.ID) %>% 
          mutate(HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = max(HAPLOTYPES, na.rm = TRUE))) %>% 
          dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES")
        
      }
    }
    
    # transform the imputed dataset into genepop object ------------------------
    message("Imputed haplotypes into genepop factory ...")
    
    haplo.imp <- suppressWarnings(
      haplo.imp %>% 
        # gather(Catalog.ID, HAPLOTYPES, -c(INDIVIDUALS, POP_ID)) %>% 
        melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
        separate(
          col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
          sep = "/", extra = "drop", remove = T
        ) %>%
        mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2))
    )
    
    message("step 1/4: completed")
    
    haplo.imp <- haplo.imp %>% 
      melt(
        id.vars = c("Catalog.ID", "INDIVIDUALS", "POP_ID"),
        measure.vars = c("ALLELE1", "ALLELE2"), 
        variable.name = "ALLELE", 
        value.name = "NUCLEOTIDES"
      ) %>% 
      dcast(INDIVIDUALS + POP_ID + ALLELE ~ Catalog.ID, value.var = "NUCLEOTIDES")
    
    message("step 2/4: completed")
    
    haplo.imp <- suppressWarnings(
      haplo.imp %>% 
        colwise(factor, exclude = "NA")(.)
    )
    message("step 3/4: completed")
    # This part is different than gtypes...
    haplo.imp <- suppressWarnings(
      haplo.imp %>%
        mutate(GROUP = rep(1, times = nrow(.))) %>% 
        group_by(GROUP) %>% 
        mutate_each(funs(as.integer), -c(ALLELE, INDIVIDUALS, POP_ID, GROUP)) %>%
        ungroup() %>% 
        select(-GROUP, -POP_ID) %>% 
        melt(id.vars = c("INDIVIDUALS", "ALLELE"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>%
        mutate(HAPLOTYPES = as.character(HAPLOTYPES)) %>% 
        mutate(HAPLOTYPES = stri_pad_left(str = HAPLOTYPES, width = 3, pad = "0")) %>% 
        mutate(HAPLOTYPES = stri_replace_na(str = HAPLOTYPES, replacement = "000")) %>% 
        dcast(Catalog.ID + INDIVIDUALS ~ ALLELE, value.var = "HAPLOTYPES") %>%
        unite(GENOTYPE, ALLELE1:ALLELE2, sep = "") %>%
        arrange(Catalog.ID) %>% 
        dcast(INDIVIDUALS ~ Catalog.ID, value.var = "GENOTYPE") %>% 
        mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))
    )
    
    message("step 4/4: completed")
    
    # Create a vector with the population ordered by levels
    if (missing(pop.levels) == "TRUE") {
      pop <- substr(haplo.imp$INDIVIDUALS, pop.id.start, pop.id.end)
    } else {
      pop <- factor(substr(haplo.imp$INDIVIDUALS, 
                           pop.id.start, pop.id.end),
                    levels = pop.levels, ordered = T)
      
    }
    
    # split genepop by populations
    genepop.imp.split <- split(haplo.imp, pop)
    
    
    # Write the file in genepop format 
    message("Output of imputed data...")
    
    
    # genepop filename
    if (missing(genepop.filename) == "TRUE"){
      genepop.filename <- "genepop.gen"
    } else {
      genepop.filename <- genepop.filename
    }
    # Add "_imputed" to the filenam
    genepop.filename.imp <- stri_replace_all_fixed(genepop.filename,
                                                   pattern = ".gen",
                                                   replacement = "_imputed.gen")
    
    
    # genepop header default
    if (missing(genepop.header) == "TRUE"){
      genepop.header <- "my first genepop"
      genepop.header <- stri_join(genepop.header, "with imputation", sep = " ")
    } else {
      genepop.header <- genepop.header
      genepop.header <- stri_join(genepop.header, "with imputation", sep = " ")
    }
    
    
    file.imputation <- file(genepop.filename.imp, "write")
    cat(genepop.header, "\n", file = file.imputation, append = TRUE)
    cat(loci, sep = "\n", file = file.imputation, append = TRUE)
    
    for (i in 1:length(genepop.imp.split)) {
      cat("pop\n", file = file.imputation, append = TRUE)
      write.table(
        genepop.imp.split[[i]],
        file = file.imputation, 
        append = TRUE, 
        col.names = FALSE,
        row.names = FALSE, 
        sep = " ", 
        quote = FALSE
      )
    }
    close(file.imputation)
    
    invisible(cat(sprintf("Genepop (with imputation) file name:\n%s\n\nWritten in the directory:\n%s", genepop.filename.imp, getwd())))
  }
}
