# Write a strataG gtypes file from STACKS haplotypes file

# to get rid of notes in build check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE1", "ALLELE2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX", "contains", "other", ".", "ALLELE"))


#' @name haplo2gtypes
#' @title Convert between batch_x.haplotypes.tsv and \code{strataG.dev} \code{gtypes} object
#' @description This function can first filter the haplotypes file 
#' with a whitelist of loci
#' and a blacklist of individuals (optional). Then it will convert the file
#' to a \code{strataG.dev} \code{gtypes} object.
#' Map-independent imputation using Random Forest is also available
#' as an option.
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param filename The name of the file written to the directory (optional).
#' @param pop.levels An optional character string with your populations ordered.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param description a string naming or describing this dataset.
#' @param imputation.rf Logical. Should a map-independent imputation of markers 
#' using Random Forest be enabled. This will write to the directory 2 files, 
#' a non-imputed and an imputed genepop files.
#' @param imputation.group \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. 
#' Default = \code{"populations"}.
#' @param num.tree The number of trees to grow. Default is 100.
#' @param iteration.rf Number of iterations of missing data algorithm.
#' Default is 10.
#' @param split.number Non-negative integer value used to specify 
#' random splitting. Default is 100.
#' @param verbose Logical. Should trace output be enabled on each iteration?
#' Default is \code{FALSE}.
#' @param parallel.core (optional) The number of core for OpenMP shared-memory parallel
#' programming of Random Forest imputation. For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.
#' @details The imputation requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set 
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals 
#' will require 15 min.
#' @return When no imputation is selected an object of the 
#' class \code{\link[strataG.devel]{gtypes}} is returned.
#' When imputation is selected a list with 2 objects is returned
#' and accessed with \code{$gtypes.no.imputation} or \code{$gtypes.imputed}
#' @export
#' @rdname haplo2gtypes
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_pad
# @importFrom strataG.devel df2gtypes
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
#' @seealso \code{strataG.devel} is available on github \url{https://github.com/EricArcher/}
#' \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

haplo2gtypes <- function(haplotypes.file,
                         whitelist.loci = NULL, 
                         blacklist.id = NULL, 
                         filename, 
                         pop.levels,
                         pop.id.start, 
                         pop.id.end, 
                         description,
                         imputation.rf = FALSE,
                         imputation.group = "populations",
                         num.tree = 100,
                         iteration.rf = 10,
                         split.number = 100,
                         verbose = FALSE,
                         parallel.core = 2
) {
  
  if (imputation.rf == "FALSE") {
    message("haplo2gtypes: without imputation...")
  } else {
    message("haplo2gtypes: with imputation...")
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
    message("Using the blacklist of id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)    
  } else if (missing(blacklist.id) == "FALSE" & is.vector(blacklist.id) == "FALSE") {
    message("Using blacklist of id from your global environment")
    blacklist.id <- blacklist.id
    
  } else {
    message("No blacklist of id")
    blacklist.id <- NULL
  }
  
  
  if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == TRUE) {
    
    # Combination 1: No whitelist and No blacklist -----------------------------
    haplotype <- haplotype    
    
  } else if (is.null(whitelist.loci) == FALSE & is.null(blacklist.id) == TRUE) {
    
    # Combination 2: Using whitelist, but No blacklist -------------------------
    
    # just whitelist.loci, NO Blacklist of individual
    haplotype <- haplotype %>%
      semi_join(whitelist, by = "Catalog.ID") %>% 
      arrange(Catalog.ID)
    
  } else if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == FALSE) {
    
    # Combination 3: Using a blacklist of id, but No whitelist -----------------
    
    # NO whitelist, JUST Blacklist of individual
    haplotype <- haplotype %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(Catalog.ID)
    
  } else {
    # Combination 4: Using a whitelist and blacklist---------------------------
    
    # whitelist.loci + Blacklist of individual
    haplotype <- haplotype %>%
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
  
  # Conversion into gtypes -----------------------------------------------------
  
  message("Haplotypes into strataG gtypes factory ...")
  
  # Haplo prep
  # Remove paralogs
  # Remove consensus loci
  
  haplo.prep <- suppressWarnings(
    haplotype %>%
      anti_join(paralogs, by = "Catalog.ID") %>%
      filter(HAPLOTYPES != "consensus") %>%    
      mutate(
        HAPLOTYPES = stri_replace_all_fixed(HAPLOTYPES, "-", "NA", 
                                            vectorize_all=F),
        POP_ID = factor(substr(INDIVIDUALS, pop.id.start, pop.id.end), 
                        levels = pop.levels, ordered = T)
      ) %>%
      dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES") #%>% 
  )
  
  # dump unused objects
  paralogs <- NULL
  haplotype <- NULL
  
  
  # No imputation --------------------------------------------------------------
  # Further work on the data
  gtypes.prep <- suppressWarnings(
    haplo.prep %>% 
      # gather(Catalog.ID, HAPLOTYPES, -c(INDIVIDUALS, POP_ID)) %>% 
      melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
      separate(
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
      dcast(INDIVIDUALS + POP_ID + ALLELE ~ Catalog.ID, value.var = "NUCLEOTIDES"))
  
  gtypes.prep <- suppressWarnings(
    gtypes.prep %>% 
      colwise(factor, exclude = "NA")(.)
  )
  
  gtypes.prep <- suppressWarnings(
    gtypes.prep %>%
      mutate_each(funs(as.integer), -c(ALLELE, INDIVIDUALS, POP_ID)) %>%
      melt(id.vars = c("INDIVIDUALS", "POP_ID", "ALLELE"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>% 
      mutate(HAPLOTYPES = as.character(HAPLOTYPES)) %>% 
      mutate(HAPLOTYPES = stri_pad_left(str = HAPLOTYPES, width = 3, pad = "0")) %>% 
      mutate(HAPLOTYPES = stri_replace_na(str = HAPLOTYPES, replacement = "000")) %>% 
      dcast(INDIVIDUALS + POP_ID ~ Catalog.ID + ALLELE, value.var = "HAPLOTYPES")
  )
  
  # results no imputation--------------------------------------------------------------------
  # Write file to directory (optional)
  if (missing(filename) == "FALSE") {
    message("Saving the file without imputation in your working directory...")
    write_tsv(gtypes.prep, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
  }
  
  # convert to gtypes
  res <- strataG.devel::df2gtypes(x = gtypes.prep,
                                  ploidy = 2,
                                  id.col = "INDIVIDUALS",
                                  strata.col = "POP_ID",
                                  loc.col = 3, 
                                  sequences = NULL,
                                  description = description
  )
  
  # dump unused objects
  gtypes.prep <- NULL
  
  # Imputation: gtypes with imputed data using Random Forest ------------------
  
  if (imputation.rf == "TRUE") {
    message("Imputation was selected\nCalculating map-independent imputation using Random Forest")
    
    # A different format is required for the imputation 
    
    # Transformed columns into factor excluding the "NA"
    haplo.prep <- suppressWarnings(
      plyr::colwise(factor, exclude = "NA")(haplo.prep)
    )
    
    # Parallel computations options
    if (missing(parallel.core) == "TRUE"){
      # Automatically select all the core -1 
      options(rf.cores=detectCores()-1, mc.cores=detectCores()-1)
    } else {
      options(rf.cores = parallel.core, mc.cores = mc.cores)
    }
    
    # Imputation using Random Forest with the package randomForestSRC
    
    impute_markers_rf <- function(x){
      randomForestSRC::impute.rfsrc(data = x, 
                                    ntree = num.tree, 
                                    nodesize = 1, 
                                    nsplit = split.number, 
                                    nimpute = iteration.rf, 
                                    do.trace = verbose)
    }
    
    # Imputation by populations (default) or globally -------------------------
    
    # default by pop
    if (missing(imputation.group) == "TRUE" | imputation.group == "populations"){
      message("Imputation computed by populations, take a break...")
      
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
        # message of progress for imputation by population
        pop.imputed <- paste("Completed imputation for pop ", i, sep = "")
        message(pop.imputed)
      }
      haplo.imp <- as.data.frame(bind_rows(imputed.dataset))
      
      # dump unused objects
      haplo.prep <- NULL
      df.split.pop <- NULL
      pop.list <- NULL
      sep.pop <- NULL
      imputed.dataset <- NULL
      
    } else if (imputation.group == "global"){
      # Globally (not by pop_id)
      message("Imputation computed globally, take a break...")
      
      INDIVIDUALS <- haplo.prep %>% select(INDIVIDUALS)
      haplo.imp <- bind_cols(INDIVIDUALS,
                             haplo.prep %>%
                               select(-INDIVIDUALS) %>% 
                               colwise(factor, exclude = "NA")(.)%>% 
                               do(impute_markers_rf(.))
      )
      # dump unused objects
      haplo.prep <- NULL
    }
    
    # transform the imputed dataset into gtypes object ------------------------
    message("Imputed data into gtypes factory...")
    haplo.imp <- suppressWarnings(
      haplo.imp %>% 
        # gather(Catalog.ID, HAPLOTYPES, -c(INDIVIDUALS, POP_ID)) %>% 
        melt(id.vars = c("INDIVIDUALS", "POP_ID"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES")
    ) %>% 
      separate(
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
    
    haplo.imp <- suppressWarnings(
      haplo.imp %>%
        colwise(factor, exclude = "NA")(.)
    )
    
    haplo.imp <- haplo.imp %>%
      mutate_each(funs(as.integer), -c(ALLELE, INDIVIDUALS, POP_ID)) %>%
      melt(id.vars = c("INDIVIDUALS", "POP_ID", "ALLELE"), variable.name = "Catalog.ID", value.name = "HAPLOTYPES") %>%
      mutate(HAPLOTYPES = as.character(HAPLOTYPES)) %>%
      mutate(HAPLOTYPES = stri_pad_left(str = HAPLOTYPES, width = 3, pad = "0")) %>%
      mutate(HAPLOTYPES = stri_replace_na(str = HAPLOTYPES, replacement = "000")) %>%
      dcast(INDIVIDUALS + POP_ID ~ Catalog.ID + ALLELE, value.var = "HAPLOTYPES")
    
    # results WITH imputation---------------------------------------------------
    # 1) the gtypes without imputation is modified and put in a new list
    no.imputation <- res
    res <- list()
    res$no.imputation <- no.imputation
    
    # 2) the gtypes with imputation
    res$imputed <- strataG.devel::df2gtypes(x = haplo.imp,
                                            ploidy = 2,
                                            id.col = "INDIVIDUALS",
                                            strata.col = "POP_ID",
                                            loc.col = 3, 
                                            sequences = NULL,
                                            description = description
    )
    
  }
  
  # outout results -------------------------------------------------------------
  # Message at the end
  invisible(cat(sprintf(
    "%s\n
Working directory:
%s",
    saving, getwd()
  )))
  message("A large 'gtypes' object was created in your Environment")
  
  return(res)
}
