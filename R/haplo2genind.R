# Write a adegenet genind object from STACKS haplotypes file

# to get rid of notes in build check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE1", "ALLELE2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX", "other", "strata", "hierarchy"))


#' @name haplo2genind
#' @title Convert between batch_x.haplotypes.tsv and \code{adegenet} 
#' \code{\link[adegenet]{genind}} object
#' @description This function can first filter the haplotypes file 
#' with a whitelist of loci
#' and a blacklist of individuals (optional). Then it will convert the file
#' to a \code{adegenet} \code{\link[adegenet]{genind}} object.
#' Map-independent imputation using Random Forest is also available
#' as an option.
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param pop.levels (optional) A character string with your populations ordered.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param strata (optional) A data frame that defines population stratifications
#'   for your samples. This is especially useful if you have a hierarchical or
#'   factorial sampling design. See \code{\link[adegenet]{genind}} for details.
#' @param hierarchy (optional) A formula that explicitely defines hierarchical levels 
#' in your strata. See \code{\link[adegenet]{genind}} for details.
#' @param imputation.rf Logical. Should a map-independent imputation of markers 
#' using Random Forest be enabled. This will write to the directory 2 files, 
#' a non-imputed and an imputed genepop files.
#' @param imputation.group Should the imputations be computed globally or by populations. 
#' \code{Default = "populations"}.
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
#' @details The imputation requires more time and can take up to an hour
#' depending on the size of the dataset. A data set with 30\% of missing data,
#' 5 000 haplotypes loci and 500 individuals will require 15-20 min.
#' @return When no imputation is selected an object of the 
#' class \code{\link[adegenet]{genind}} is returned.
#' When imputation is selected a list with 2 objects is returned
#' and accessed with \code{$genind.no.imputation} or \code{$genind.imputed}
#' @export
#' @rdname haplo2genind
# @importFrom adegenet df2genind
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_pad
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
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @seealso \code{adegenet} is available on CRAN \url{http://cran.r-project.org/web/packages/adegenet/} and github \url{https://github.com/thibautjombart/}
#' \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

haplo2genind <- function(haplotypes.file, 
                         whitelist.loci = NULL, 
                         blacklist.id = NULL, 
                         pop.levels, pop.id.start, pop.id.end,
                         strata = NULL, hierarchy = NULL,
                         imputation.rf = FALSE,
                         imputation.group = "populations",
                         num.tree = 100,
                         iteration.rf = 10,
                         split.number = 100,
                         verbose = FALSE,
                         parallel.core = 2) {
  
  if (imputation.rf == "FALSE") {
    message("haplo2genind: without imputation...")
  } else {
    message("haplo2genind: with imputation...")
  }
  
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
  
  
  # Haplotype file--------------------------------------------------------------
  haplotype <- read_tsv(file = haplotypes.file, col_names = T) %>%
    rename(Catalog.ID = `Catalog ID`) %>%
    gather(INDIVIDUALS, HAPLOTYPES, -c(Catalog.ID, Cnt)) %>%
    arrange(Catalog.ID)
  
  
  
  if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == TRUE) {
    
    # Combination 1: No whitelist and No blacklist----------------------------------
    
    # No filter
    haplotype.no.filter <- haplotype    
    
    # Creates a vector containing the loci name
    loci <- unique(haplotype.no.filter$Catalog.ID)
    data <- haplotype.no.filter
    
  } else if (is.null(whitelist.loci) == FALSE & is.null(blacklist.id) == TRUE) {
    
    # Combination 2: Using whitelist, but No blacklist--------------------------
    
    # just whitelist.loci, NO Blacklist of individual
    haplotype.whitelist.loci <- haplotype %>%
      semi_join(whitelist, by = "Catalog.ID") %>% # instead of right_join
      #       might result in less loci then the whitelist, but this make sure
      #   we don't end up with unwanted or un-caracterized locus, e.g. comparing
      #   adults and juveniles samples... the whitelist is developed with the adults,
      #   but the haplotype file comes from the juveniles... 
      arrange(Catalog.ID)
    
    # Creates a vector containing the loci name
    loci <- unique(haplotype.whitelist.loci$Catalog.ID)
    data <- haplotype.whitelist.loci
    
  } else if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == FALSE) {
    
    # Combination 3: Using a blacklist of id, but No whitelist------------------
    
    # NO whitelist, JUST Blacklist of individual
    haplotype.blacklist <- haplotype %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(Catalog.ID)
    
    
    # Creates a vector containing the loci name
    loci <- unique(haplotype.blacklist$Catalog.ID)
    data <- haplotype.blacklist
    
  } else {
    # Combination 4: Using a whitelist and blacklist---------------------------
    
    # whitelist.loci + Blacklist of individual
    haplotype.whitelist.blacklist <- haplotype %>%
      semi_join(whitelist, by = "Catalog.ID") %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(Catalog.ID)
    
    # Creates a vector containing the loci name
    loci <- unique(haplotype.whitelist.blacklist$Catalog.ID)
    data <- haplotype.whitelist.blacklist
  }
  
  # Paralogs-------------------------------------------------------------------
  message("Looking for paralogs...")
  
  paralogs <- data %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(Catalog.ID) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    group_by(Catalog.ID) %>%
    select(Catalog.ID) %>%
    distinct(Catalog.ID)
  
  nparalogs <- stri_join("Found and/or removed", n_distinct(paralogs$Catalog.ID), "paralogs", sep = " ")
  message(nparalogs)
  
  # Conversion into genind -----------------------------------------------------
  message("Haplotypes into genind factory for conversion...")
  
  # genind prep
  # removes paralogs
  # removes consensus loci
  
  genind.prep <- data %>%
    anti_join(paralogs, by = "Catalog.ID") %>%
    filter(HAPLOTYPES != "consensus") %>%    
    mutate(
      HAPLOTYPES = stri_replace_all_fixed(HAPLOTYPES, "-", "NA", 
                                          vectorize_all=F)
    ) %>%
    separate(
      col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
      sep = "/", extra = "drop", remove = T
    ) %>%
    mutate(
      ALLELE1 = stri_replace_all_fixed(ALLELE1, "NA", "0", vectorize_all=F),
      ALLELE2 = stri_replace_na(ALLELE2, replacement = "no_allele"),
      ALLELE2 = ifelse(ALLELE2 == "no_allele", ALLELE1, ALLELE2)
    ) %>%
    melt(id.vars = c("Catalog.ID","Cnt", "INDIVIDUALS"),
         measure.vars = c("ALLELE1", "ALLELE2"), 
         variable.name = "ALLELE", 
         value.name = "NUCLEOTIDES"
    ) %>%
    group_by(Catalog.ID) %>%
    mutate(
      NUCLEOTIDES = as.numeric(factor(NUCLEOTIDES)),
      NUCLEOTIDES = NUCLEOTIDES-1 # this removes 1 levels to enable missing values = 0
    )
  
  
  # No imputation --------------------------------------------------------------
  
  # Further work on the data, pad the genotypes for 3 digits output
  genind.data <- genind.prep %>% 
    mutate(
      NUCLEOTIDES = str_pad(NUCLEOTIDES, 3, side = "left", pad = "0")
    ) %>%
    select(-Cnt) %>%
    dcast(Catalog.ID + INDIVIDUALS ~ ALLELE, value.var = "NUCLEOTIDES") %>%
    unite(GENOTYPE, ALLELE1:ALLELE2, sep = "/") %>%
    dcast(INDIVIDUALS ~ Catalog.ID, value.var = "GENOTYPE") %>% 
    mutate(
      POP_ID = factor(substr(INDIVIDUALS, 
                             pop.id.start, pop.id.end),
                      levels = pop.levels, ordered = T)
    )
  # results no imputation--------------------------------------------------------------------
  
  ind <- genind.data$INDIVIDUALS
  pop <- genind.data$POP_ID
  genind.df <- genind.data %>%
    select(-c(INDIVIDUALS, POP_ID))
  res <- adegenet::df2genind(X = genind.df, sep = "/", ind.names = ind, pop = pop, ploidy = 2, strata = strata, hierarchy = hierarchy)
  
  
  # Imputation: genind with imputed data using Random Forest ------------------
  
  if (imputation.rf == "TRUE") {
    message("Imputation was selected\nCalculating map-independent imputation using Random Forest")
    
    # A different format is required for the imputation 
    imp <- suppressWarnings(
      genind.prep %>% 
        mutate(
          NUCLEOTIDES = stri_replace_all_fixed(NUCLEOTIDES, "0", "NA", 
                                               vectorize_all=F),
          POP_ID = factor(substr(INDIVIDUALS, pop.id.start, pop.id.end), 
                          levels = pop.levels, ordered = T)
        ) %>%
        arrange(Catalog.ID) %>% 
        group_by(Catalog.ID) %>%
        select(-Cnt) %>%
        dcast(INDIVIDUALS+POP_ID ~ Catalog.ID + ALLELE, value.var = "NUCLEOTIDES")
    )
    
    # Transformed columns into factor excluding the "NA"
    imp <- suppressWarnings(
      plyr::colwise(factor, exclude = "NA")(imp)
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
    
    df.split.pop <- split(x = imp, f = imp$POP_ID) # slip data frame by population
    pop.list <- names(df.split.pop) # list the pop
    imputed.dataset <-list() # create empty list 
    genind.imp <- list() # create empty list 
    for (i in pop.list) {
      # group <- paste(i)
      sep.pop <- df.split.pop[[i]]
      imputed.dataset[[i]] <- impute_markers_rf(sep.pop)
      # message of progress for imputation by population
      pop.imputed <- paste("Completed imputation for pop ", i, sep = "")
      message(pop.imputed)
    }
    genind.imp <- as.data.frame(bind_rows(imputed.dataset))
    
    # transform the imputed dataset into genind object
    
    message("Imputed data into genind factory...")
    
    genind.imp <- suppressWarnings(
      genind.imp %>%
        gather(Catalog.ID, NUCLEOTIDES, -c(INDIVIDUALS, POP_ID)) %>%
        separate(Catalog.ID, c("Catalog.ID", "ALLELE"), sep = "_", extra = "error") %>%
        mutate(Catalog.ID = as.integer(Catalog.ID)) %>% 
        arrange(Catalog.ID) %>% 
        mutate(NUCLEOTIDES = str_pad(NUCLEOTIDES, 3, side = "left", pad = "0")) %>%
        dcast(Catalog.ID + INDIVIDUALS + POP_ID ~ ALLELE, value.var = "NUCLEOTIDES") %>%
        unite(GENOTYPE, ALLELE1:ALLELE2, sep = "") %>%
        dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "GENOTYPE") %>%
        mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
    )
    # results WITH imputation---------------------------------------------------
    # 1) the genind without imputation is modified and put in a new list
    genind.no.imputation <- res
    res <- list()
    res$genind.no.imputation <- genind.no.imputation
    
    ind <- genind.imp$INDIVIDUALS
    pop <- genind.imp$POP_ID
    genind.df <- genind.imp %>%
      select(-c(INDIVIDUALS, POP_ID))
    
    res$genind.imputed <- adegenet::df2genind(X = genind.df, sep = "/",
                                              ind.names = ind,
                                              pop = pop,
                                              ploidy = 2,
                                              strata = strata,
                                              hierarchy = hierarchy
    )
  }
  # outout results -------------------------------------------------------------
  message("A large genind object was created in your Environment")
  return(res)
}

