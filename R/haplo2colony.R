# From STACKS haplotypes file write COLONY input files to the working directory

# to get rid of notes in build check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE", "ALLELE1", "ALLELE2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX", "other", "strata", "hierarchy", "GROUP", "N", "."))


#' @name haplo2colony
#' @title Conduct parentage analysis in COLONY easily by converting 
#' a \code{batch_x.haplotypes.tsv} STACKS file into a \code{COLONY} input files 
#' @description This function can first filter the haplotypes file 
#' with a whitelist of loci and a blacklist of individuals (optional). 
#' Then it will convert the file to the required \code{COLONY} input files.
#' Map-independent imputation using Random Forest or the most frequent category
#' is also available as an option.
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.txt").
#' @param whitelist.loci (optional) A whitelist of loci and 
#' a column header 'LOCUS'.
#' The whitelist is in the directory (e.g. "whitelist.txt").
#' @param sample.markers (number) Should COLONY output contain a random 
#' subsample of your markers. 
#' Default= \code{"all"}.
#' e.g. \code{sample.markers = 500} to use only 500 randomly chosen markers.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.select (string) Should the COLONY output be on all populations
#' or a select group. Default = \code{"all"}. 
#' e.g. \code{pop.select = "QUE"} to select QUE population sample. 
#' \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population sample.
#' @param allele.freq (optional string) Population string chosen to estimate
#' the allele frequency for each locus. Default = \code{FALSE}. 
#' e.g. \code{allele.freq = "QUE"} or \code{allele.freq = c("QUE", "ONT")} 
#' or \code{allele.freq = "overall"}.
#' @param inbreeding (boolean) 0/1 no inbreeding/inbreeding. Default = \code{0}
#' @param mating.sys.males (boolean) Mating system in males.
#' 0/1 polygyny/monogyny. Default = \code{0}.
#' @param mating.sys.females (boolean) Mating system in females.
#' 0/1 polygyny/monogyny. Default = \code{0}.
#' @param clone (boolean) Should clones and duplicated individuals be inferred.
#' 0/1, yes/no. Default = \code{0}.
#' @param run.length (integer) Length of run. 1 (short), 2 (medium), 3 (long),
#' 4 (very long). Default = \code{2}.
#' @param analysis (integer) Analysis method. 
#' 0 (Pairwise-Likelihood Score), 1 (Full Likelihood),
#' 2 (combined Pairwise-Likelihood Score and Full Likelihood).
#' Default = \code{1}.
#' @param allelic.dropout Locus allelic dropout rate. Default = \code{0} See Colony manual.
#' @param error.rate Locus error rate. Default = \code{0.02}, see COLONY manual.
#' @param print.all.colony.opt (logical) Should all COLONY options be printed in the file.
#' This require manual curation, for the file to work directly with COLONY. 
#' Default = \code{FALSE}.
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
#' @param parallel.core (optional) The number of core for OpenMP 
#' shared-memory parallel programming of Random Forest imputations. 
#' For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.
#' @param filename Name of the acronym for filenaming in the working directory.
#' @details It is highly recommended to read the user guide distributed with
#' COLONY to find out the details for input and output of the software. 
#' Not all options are provided here. 
#' However, all required options to run COLONY will be printed 
#' in the file written in your working directory. Change the values accordingly.
#' The imputations using Random Forest requires more time to compute 
#' and can take several minutes and hours depending on the size of 
#' the dataset and polymorphism of the species used. e.g. with a 
#' low polymorphic taxa, and a data set containing 30\% missing data, 
#' 5 000 haplotypes loci and 500 individuals will require 15 min.
#' @return When no imputation is selected ....
#' When imputation is selected a list with 2 objects is returned
#' and accessed with \code{$no.imputation} or \code{$imputed}.
#' @export
#' @rdname haplo2colony
#' @import reshape2
#' @import dplyr
#' @import lazyeval
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
#' @references Jones OR, Wang J (2010) COLONY: a program for parentage and 
#' sibship inference from multilocus genotype data. 
#' Molecular Ecology Resources, 10, 551–555.
#' @references Wang J (2012) Computationally Efficient Sibship and 
#' Parentage Assignment from Multilocus Marker Data. Genetics, 191, 183–194.
#' @seealso COLONY is available on Jinliang Wang www 
#' \url{http://www.zsl.org/science/software/colony}
#' \code{randomForestSRC} is available on 
#' CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} 
#' and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

haplo2colony <- function(haplotypes.file,
                         blacklist.id = NULL,
                         whitelist.loci = NULL,
                         sample.markers = "all",
                         pop.id.start, pop.id.end,
                         pop.select = "all",
                         allele.freq = FALSE,
                         inbreeding = 0,
                         mating.sys.males = 0,
                         mating.sys.females = 0,
                         clone = 0,
                         run.length =2,
                         analysis = 1,
                         allelic.dropout = 0,
                         error.rate = 0.02,
                         print.all.colony.opt = FALSE,
                         imputations = FALSE,
                         imputations.group = "populations",
                         num.tree = 100,
                         iteration.rf = 10,
                         split.number = 100,
                         verbose = FALSE,
                         parallel.core = 2,
                         filename = "colony"
) {
  
  if (imputations == "FALSE") {
    message("haplo2colony: without imputation...")
  } else {
    message("haplo2colony: with imputations...")
  }
  
  
  # Haplotype file--------------------------------------------------------------
  haplotype <- read_tsv(file = haplotypes.file, col_types = cols(.default = col_character()), col_names = T, na = "-") %>% 
    select(-Cnt) %>% 
    rename(Catalog.ID = `Catalog ID`) %>%
    melt(id.vars = "Catalog.ID", variable.name = "INDIVIDUALS", value.name = "HAPLOTYPES") %>% 
    mutate(
      HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = "-"),
      Catalog.ID = as.integer(Catalog.ID),
      POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end)
      )
  
  # Pop select -----------------------------------------------------------------
  if(pop.select == "all" | missing(pop.select) == TRUE){
    haplotype <- haplotype
  } else {
    target.pop <- pop.select
    haplotype <- suppressWarnings(
      haplotype %>% 
        filter(POP_ID %in% target.pop)
    )
  }
  
  
  # Whitelist-------------------------------------------------------------------
  if (is.null(whitelist.loci) | missing(whitelist.loci)) {
    message("No whitelist")
    whitelist <- NULL
  } else if (is.vector(whitelist.loci)) {
    message("Using the whitelist from the directory")
    whitelist <- read_tsv(whitelist.loci, col_names = T) %>%
      rename(Catalog.ID = LOCUS)
  } else {
    message("Using whitelist from your global environment")
    whitelist <- whitelist.loci %>%
      rename(Catalog.ID = LOCUS)
  }
  
  
  
  # Blacklist-------------------------------------------------------------------
  if (is.null(blacklist.id) | missing(blacklist.id)) {
    message("No individual blacklisted")
    blacklist.id <- NULL
  } else if (is.vector(blacklist.id)) {
    message("Using the blacklisted id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)    
  } else {
    message("Using the blacklisted id from your global environment")
    blacklist.id <- blacklist.id
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
    distinct(Catalog.ID)
  
  nparalogs <- stri_join("Found and/or removed", n_distinct(paralogs$Catalog.ID), "paralogs", sep = " ")
  message(nparalogs)
  
  
  # Subsampling markers --------------------------------------------------------
  if(sample.markers == "all" | missing(sample.markers) == TRUE){
    
    message("Haplotypes into factory for conversion into COLONY ...")
    
    # Haplo prep
    # Remove paralogs
    # Remove consensus loci
    
    haplo.filtered <- suppressWarnings(
      haplotype %>%
        anti_join(paralogs, by = "Catalog.ID") %>%
        filter(HAPLOTYPES != "consensus") %>%    
        mutate(
          HAPLOTYPES = stri_replace_all_fixed(HAPLOTYPES, "-", "NA", 
                                              vectorize_all=F)
        )%>%
        dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES") #%>% 
    )
  } else{
    sampling.markers.message <- paste("Randomly subsampling ", sample.markers, " markers...")
    message(sampling.markers.message)
    
    haplo.filtered <- suppressWarnings(
      haplotype %>%
        anti_join(paralogs, by = "Catalog.ID") %>%
        filter(HAPLOTYPES != "consensus") %>%
        select(-POP_ID) %>% 
        dcast(Catalog.ID ~ INDIVIDUALS, value.var = "HAPLOTYPES") %>% 
        # group_by(Catalog.ID) %>% 
        do(sample_n(., sample.markers)) %>% 
        arrange(Catalog.ID) %>%
        melt(id.vars = "Catalog.ID", variable.name = "INDIVIDUALS", value.name = "HAPLOTYPES") %>%    
        mutate(
          HAPLOTYPES = stri_replace_all_fixed(HAPLOTYPES, "-", "NA", 
                                              vectorize_all=F),
          POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end)
        ) %>% 
        dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES")
    )
    message("Haplotypes into factory for conversion into COLONY ...")
  }
  message("step 1/5: completed")
  # dump unused objects
  paralogs <- NULL
  haplotype <- NULL
  # No imputation --------------------------------------------------------------
  # Further work on the data
  haplo.prep <- suppressWarnings(
    haplo.filtered %>% 
      melt(
        id.vars = c("INDIVIDUALS", "POP_ID"), 
        variable.name = "Catalog.ID", 
        value.name = "HAPLOTYPES"
      ) %>% 
      tidyr::separate(
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
  # Allele frequency per locus
  if (allele.freq == "overall"){
    
    allele.per.locus <- haplo.prep %>% select(-INDIVIDUALS, -POP_ID, -ALLELE) %>% 
      colwise(nlevels)(.)
    
    frequency.markers <- suppressWarnings(
      haplo.prep %>%
        tidyr::gather(Catalog.ID, NUCLEOTIDES, -c(INDIVIDUALS, POP_ID, ALLELE)) %>% 
        select(-INDIVIDUALS, -ALLELE, -POP_ID) %>% 
        group_by(Catalog.ID) %>%
        filter(NUCLEOTIDES != "NA") %>% 
        mutate(N = length(NUCLEOTIDES)) %>% 
        group_by(Catalog.ID, NUCLEOTIDES) %>%
        mutate(n = n()) %>% 
        distinct(Catalog.ID, NUCLEOTIDES, .keep_all = TRUE) %>% 
        mutate(FREQ = n/N) %>% 
        select(-N, -n) %>% 
        group_by(Catalog.ID) %>% 
        mutate(NUCLEOTIDES = as.numeric(factor(NUCLEOTIDES))) %>% 
        arrange(Catalog.ID, NUCLEOTIDES) %>% 
        dcast(Catalog.ID ~ NUCLEOTIDES, value.var = "FREQ", drop = FALSE) %>% 
        select(-Catalog.ID) %>% 
        mutate_each(funs(stri_replace_na(str = ., replacement = "")))
    )    
    message("step 4/5: completed")
    
  } else if (missing(allele.freq) != TRUE & allele.freq != FALSE){
    
    allele.per.locus <- haplo.prep %>% 
      filter(POP_ID %in% allele.freq) %>%
      select(-INDIVIDUALS, -POP_ID, -ALLELE) %>% 
      colwise(nlevels)(.)
    
    frequency.markers <- suppressWarnings(
      haplo.prep %>%
        filter(POP_ID %in% allele.freq) %>% 
        tidyr::gather(Catalog.ID, NUCLEOTIDES, -c(INDIVIDUALS, POP_ID, ALLELE)) %>% 
        select(-INDIVIDUALS, -ALLELE, -POP_ID) %>% 
        group_by(Catalog.ID) %>%
        filter(NUCLEOTIDES != "NA") %>% 
        mutate(N = length(NUCLEOTIDES)) %>% 
        group_by(Catalog.ID, NUCLEOTIDES) %>%
        mutate(n = n()) %>% 
        distinct(Catalog.ID, NUCLEOTIDES, .keep_all = TRUE) %>% 
        mutate(FREQ = n/N) %>% 
        select(-N, -n) %>% 
        group_by(Catalog.ID) %>% 
        mutate(NUCLEOTIDES = as.numeric(factor(NUCLEOTIDES))) %>% 
        arrange(Catalog.ID, NUCLEOTIDES) %>% 
        dcast(Catalog.ID ~ NUCLEOTIDES, value.var = "FREQ", drop = FALSE) %>% 
        select(-Catalog.ID) %>% 
        mutate_each(funs(stri_replace_na(str = ., replacement = "")))
    )
    message("step 4/5: completed")
    
  } else {
    message("step 4/5: completed")
  }
  haplo.prep <- suppressWarnings(
    haplo.prep %>%
      mutate(GROUP = rep(1, times = nrow(.))) %>% 
      group_by(GROUP) %>% 
      mutate_each(funs(as.integer), -c(ALLELE, INDIVIDUALS, POP_ID, GROUP)) %>%
      ungroup() %>% 
      select(-GROUP) %>% 
      melt(
        id.vars = c("INDIVIDUALS", "POP_ID", "ALLELE"), 
        variable.name = "Catalog.ID", 
        value.name = "HAPLOTYPES"
      ) %>%
      mutate(HAPLOTYPES = as.character(HAPLOTYPES)) %>% 
      mutate(HAPLOTYPES = stri_pad_left(str = HAPLOTYPES, width = 3, pad = "1")) %>% 
      mutate(HAPLOTYPES = stri_replace_na(str = HAPLOTYPES, replacement = "0")) %>%
      group_by(Catalog.ID) %>% 
      arrange(Catalog.ID)
  )
  markers.name <- haplo.prep %>% distinct(Catalog.ID)
  marker.num <- nrow(markers.name)
  markers.name <- t(markers.name)
  
  haplo.prep <- haplo.prep %>% 
    dcast(INDIVIDUALS + POP_ID ~ Catalog.ID + ALLELE, value.var = "HAPLOTYPES")
  
  message("step 5/5: completed")
  
  res <- haplo.prep %>% select(-POP_ID)
  
  # results no imputation-------------------------------------------------------
  # convert to colony
  # Line 1 = Dataset name
  dataset.opt <- "`My first COLONY run`                ! Dataset name"
  dataset.opt <- as.data.frame(dataset.opt)
  write.table(x = dataset.opt, file = filename, sep = " ", append = FALSE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 2 = Output filename
  # out.name.opt <- paste("`", filename, "`", "                         ! Output file name", sep = "")
  colony.output.filename <- stri_replace_all_fixed(filename,
                                                   pattern = ".dat",
                                                   replacement = "")
  
  out.name.opt <- paste(colony.output.filename, "                         ! Output file name", sep = "") 
  out.name.opt <- as.data.frame(out.name.opt)
  write.table(x = out.name.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 3 = Offspring number
  off.num.opt <- nrow(res)
  off.num.opt <- paste(off.num.opt, "                                  ! Number of offspring in the sample", sep = "")
  off.num.opt <- as.data.frame(off.num.opt)
  write.table(x = off.num.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 4 = Number of loci
  marker.number <- (ncol(haplo.prep)-2)/2
  marker.num.opt <- paste(marker.number, "                                 ! Number of loci", sep = "")
  marker.num.opt <- as.data.frame(marker.num.opt)
  write.table(x = marker.num.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 5 = Seed random number generator  
  seed.opt <- "1234                                 ! Seed for random number generator"
  seed.opt <- as.data.frame(seed.opt)
  write.table(x = seed.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 6 = Updating allele frequency
  update.allele.freq.opt <- "0                                    ! 0/1=Not updating/updating allele frequency"
  update.allele.freq.opt <- as.data.frame(update.allele.freq.opt)
  write.table(x = update.allele.freq.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 7 = Dioecious/Monoecious species
  dioecious.opt <- "2                                    ! 2/1=Dioecious/Monoecious species"
  dioecious.opt <- as.data.frame(dioecious.opt)
  write.table(x = dioecious.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 8 = inbreeding
  inbreeding.opt <- paste(inbreeding, "                                    ! 0/1=No inbreeding/inbreeding", sep = "") 
  inbreeding.opt <- as.data.frame(inbreeding.opt)
  write.table(x = inbreeding.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 9 = Ploidy
  ploidy.opt <- "0                                    ! 0/1=Diploid species/HaploDiploid species"
  ploidy.opt <- as.data.frame(ploidy.opt)
  write.table(x = ploidy.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 10 = Mating system (polygamous: 0, monogamous: 1).
  mating.opt <- paste(mating.sys.males,"  ", mating.sys.females, "                                 ! 0/1=Polygamy/Monogamy for males & females", sep = "") 
  mating.opt <- as.data.frame(mating.opt)
  write.table(x = mating.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 11 = Clone inference
  clone.opt <- paste(clone, "                                    ! 0/1=Clone inference =No/Yes", sep = "") 
  clone.opt <- as.data.frame(clone.opt)
  write.table(x = clone.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 12 = Sibship size scaling
  sib.size.scal.opt <- "1                                    ! 0/1=Full sibship size scaling =No/Yes"
  sib.size.scal.opt <- as.data.frame(sib.size.scal.opt)
  write.table(x = sib.size.scal.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 13 = Sibship prior indicator (Integer), average paternal sibship size (Real, optional), average maternal sibship size (Real, optional)
  sib.prior.opt <- "0 0 0                                ! 0, 1, 2, 3 = No, weak, medium, strong sibship size prior; mean paternal & maternal sibship size"
  sib.prior.opt <- as.data.frame(sib.prior.opt)
  write.table(x = sib.prior.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Line 14 = Population allele frequency indicator
  if (missing(allele.freq) | allele.freq == FALSE){
    allele.freq.ind.opt <- "0                                    ! 0/1=Unknown/Known population allele frequency"
    allele.freq.ind.opt <- as.data.frame(allele.freq.ind.opt)
    write.table(x = allele.freq.ind.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  } else{
    allele.freq.ind.opt <- "1                                    ! 0/1=Unknown/Known population allele frequency"
    allele.freq.ind.opt <- as.data.frame(allele.freq.ind.opt)
    write.table(x = allele.freq.ind.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Line 15 and more = Numbers of alleles per locus (Integer, optional). Required when the Population allele frequency indicator is set to 1
  if (allele.freq != FALSE){
    write_delim(x = allele.per.locus, path = filename, delim = " ", append = TRUE, col_names = FALSE)
  }
  
  # Alleles and their frequencies per locus (Integer, Real, optional)
  
  if (allele.freq != FALSE){
    write_delim(x = frequency.markers, path = filename, delim = " ", append = TRUE, col_names = FALSE)
  }
  
  # Number of runs
  num.run.opt <- "\n1                                    ! Number of runs"
  num.run.opt <- as.data.frame(num.run.opt)
  write.table(x = num.run.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Length of run, give a value of 1, 2, 3, 4 to indicate short, medium, long, very long run.
  run.length.opt <- paste(run.length, "                                    ! Length of run", sep = "") 
  run.length.opt <- as.data.frame(run.length.opt)
  write.table(x = run.length.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Monitor method (Time in second)
  monitor.met.opt <- "0                                    ! 0/1=Monitor method by Iterate"
  monitor.met.opt <- as.data.frame(monitor.met.opt)
  write.table(x = monitor.met.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
  # Monitor interval (Time in second)
  monitor.int.opt <- "10000                                ! Monitor interval in Iterate"
  monitor.int.opt <- as.data.frame(monitor.int.opt)
  write.table(x = monitor.int.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # WindowsGUI/DOS, 0 when running Colony in DOS mode or on other platforms
  windows.gui.opt <- "0                                    ! non-Windows version"
  windows.gui.opt <- as.data.frame(windows.gui.opt)
  write.table(x = windows.gui.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Analysis method : 0, 1 or 2 for Pairwise-Likelihood score (PLS), full likelihood method (FL), or the FL and PLS combined method (FPLS). More on these methods are explained above in the Windows GUI data input section.
  analysis.opt <- paste(analysis, "                                    ! Analysis 0 (Pairwise-Likelihood Score), 1 (Full Likelihood), 2 (combined Pairwise-Likelihood Score and Full Likelihood)", sep = "") 
  analysis.opt <- as.data.frame(analysis.opt)
  write.table(x = analysis.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Precision
  precision.opt <- "3                                    ! 1/2/3=low/medium/high Precision for Full likelihood\n"
  precision.opt <- as.data.frame(precision.opt)
  write.table(x = precision.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Marker IDs/Names
  # snp.id <- seq(from = 1, to = marker.number, by = 1)
  markers.name.opt <- c(markers.name, "        ! Marker IDs")
  write.table(x = t(as.data.frame(markers.name.opt)), file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Marker types 
  marker.type.opt <- c(rep(0, marker.number),"         ! Marker types, 0/1=Codominant/Dominant") # marker.type (codominant/dominant)
  write.table(x = t(as.data.frame(marker.type.opt)), file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Allelic dropout rates
  dropout <- c(rep(allelic.dropout, marker.number), "        ! Allelic dropout rate at each locus") # allelic dropout rate
  write.table(x = t(as.data.frame(dropout)), file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Error rates
  error <- c(rep(error.rate, marker.number), "        ! False allele rate\n\n") # genotyping error rate
  write.table(x = t(as.data.frame(error)), file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Offspring IDs and genotype
  # write_delim(x = res, path = filename, delim = " ", append = TRUE, col_names = FALSE)
  write.table(x = res, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  message("Including genotypes of individuals in COLONY format")
  
  # Probabilities that the father and mother of an offspring are included in the candidate males and females. The two numbers must be provided even if there are no candidate males or/and females.
  prob.opt <- "\n\n0  0                                 ! Prob. of dad/mum included in the candidates"
  prob.opt <- as.data.frame(prob.opt)
  write.table(x = prob.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Numbers of candidate males and females
  candidate.opt <- "0 0                                  ! Numbers of candidate males & females"
  candidate.opt <- as.data.frame(candidate.opt)
  write.table(x = candidate.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Candidate male IDs/names and genotypes
  # cat(male.genotype, sep = "\n", file = colony.filename, append = TRUE)
  if (print.all.colony.opt != FALSE){
    candidate.male.id.geno.opt <- "!Candidate male ID and genotypes"
    candidate.male.id.geno.opt <- as.data.frame(candidate.male.id.geno.opt)
    write.table(x = candidate.male.id.geno.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Candidate female IDs/names and genotypes 
  # cat(female.genotype, sep = "\n", file = colony.filename, append = TRUE)
  if (print.all.colony.opt != FALSE){
    candidate.female.id.geno.opt <- "!Candidate female ID and genotypes"
    candidate.female.id.geno.opt <- as.data.frame(candidate.female.id.geno.opt)
    write.table(x = candidate.female.id.geno.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Number of offspring with known paternity
  known.paternity.opt <- "0                                    ! Number of offspring with known father"
  known.paternity.opt <- as.data.frame(known.paternity.opt)
  write.table(x = known.paternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Known offspring-father dyad
  if (print.all.colony.opt != FALSE){
    known.father.dyad.opt <- "! Offspring ID and known father ID (Known offspring-father dyad)"
    known.father.dyad.opt <- as.data.frame(known.father.dyad.opt)
    write.table(x = known.father.dyad.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Number of offspring with known maternity
  known.maternity.opt <- "0                                    ! Number of offspring with known mother"
  known.maternity.opt <- as.data.frame(known.maternity.opt)
  write.table(x = known.maternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Known offspring-mother dyad
  if (print.all.colony.opt != FALSE){
    known.mother.dyad.opt <- "! Offspring ID and known mother ID (Known offspring-mother dyad)"
    known.mother.dyad.opt <- as.data.frame(known.mother.dyad.opt)
    write.table(x = known.mother.dyad.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Number of known paternal sibships (Integer)
  known.paternal.sibships.opt <- "0                                    ! Number of known paternal sibships"
  known.paternal.sibships.opt <- as.data.frame(known.paternal.sibships.opt)
  write.table(x = known.paternal.sibships.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
  # Paternal sibship size and members (Integer, String, optional).
  if (print.all.colony.opt != FALSE){
    known.paternal.sibships.size.opt <- "! Paternal sibship size and members"
    known.paternal.sibships.size.opt <- as.data.frame(known.paternal.sibships.size.opt)
    write.table(x = known.paternal.sibships.size.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Number of known maternal sibships (Integer)
  known.maternal.sibships.opt <- "0                                    ! Number of known maternal sibships"
  known.maternal.sibships.opt <- as.data.frame(known.maternal.sibships.opt)
  write.table(x = known.maternal.sibships.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
  # Maternal sibship size and members (Integer, String, optional).
  if (print.all.colony.opt != FALSE){
    known.maternal.sibships.size.opt <- "! Maternal sibship size and members"
    known.maternal.sibships.size.opt <- as.data.frame(known.maternal.sibships.size.opt)
    write.table(x = known.maternal.sibships.size.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Number of offspring with known excluded paternity (Integer). 
  offspring.known.excl.paternity.opt <- "0                                    ! Number of offspring with known excluded fathers"
  offspring.known.excl.paternity.opt <- as.data.frame(offspring.known.excl.paternity.opt)
  write.table(x = offspring.known.excl.paternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
  # Excluded paternity 
  if (print.all.colony.opt != FALSE){
    excl.paternity.opt <- "! Offspring ID, number of excluded fathers, and excluded father IDs"
    excl.paternity.opt <- as.data.frame(excl.paternity.opt)
    write.table(x = excl.paternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Number of offspring with known excluded maternity (Integer). 
  offspring.known.excl.maternity.opt <- "0                                    ! Number of offspring with known excluded mothers"
  offspring.known.excl.maternity.opt <- as.data.frame(offspring.known.excl.maternity.opt)
  write.table(x = offspring.known.excl.maternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
  # Excluded maternity
  if (print.all.colony.opt != FALSE){
    excl.maternity.opt <- "! Offspring ID, number of excluded mothers, and excluded father IDs"
    excl.maternity.opt <- as.data.frame(excl.maternity.opt)
    write.table(x = excl.maternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Number of offspring with known excluded paternal sibships
  offspring.known.excl.paternal.sibships.opt <- "0                                    ! Number of offspring with known excluded paternal sibships"
  offspring.known.excl.paternal.sibships.opt <- as.data.frame(offspring.known.excl.paternal.sibships.opt)
  write.table(x = offspring.known.excl.paternal.sibships.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Excluded paternal siblings
  if (print.all.colony.opt != FALSE){
    excluded.paternal.siblings.opt <- "! Excluded paternal siblings"
    excluded.paternal.siblings.opt <- as.data.frame(excluded.paternal.siblings.opt)
    write.table(x = excluded.paternal.siblings.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  
  # Number of offspring with known excluded maternal sibships
  offspring.known.excl.maternal.sibships.opt <- "0                                    ! Number of offspring with known excluded maternal sibships"
  offspring.known.excl.maternal.sibships.opt <- as.data.frame(offspring.known.excl.maternal.sibships.opt)
  write.table(x = offspring.known.excl.maternal.sibships.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Excluded maternal siblings
  if (print.all.colony.opt != FALSE){
    excluded.maternal.siblings.opt <- "! Excluded maternal siblings"
    excluded.maternal.siblings.opt <- as.data.frame(excluded.maternal.siblings.opt)
    write.table(x = excluded.maternal.siblings.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  message.write.colony <- paste("A COLONY input file project was created and saved in your working directory", filename, sep = " : ")
  
  if (imputations == "FALSE") {
    message(message.write.colony)
  } else if (imputations == "max"){
    message("Calculating map-independent imputations using the most frequent allele.")
  } else {
    message("Calculating map-independent imputations using random forest")
  }
  
  # dump unused objects
  haplo.prep <- NULL
  
  # Imputations: colony with imputed haplotypes using Random Forest ------------
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
        
      } else if (imputations.group == "global"){
        # Globally (not by pop_id)
        message("Imputations computed globally")
        
        haplo.imp <- haplo.filtered %>%
          melt(
            id.vars = c("INDIVIDUALS", "POP_ID"),
            variable.name = "Catalog.ID", 
            value.name = "HAPLOTYPES"
          ) %>% 
          mutate(HAPLOTYPES = replace(HAPLOTYPES, which(HAPLOTYPES=="NA"), NA)) %>%
          group_by(Catalog.ID) %>% 
          mutate(HAPLOTYPES = stri_replace_na(HAPLOTYPES, replacement = max(HAPLOTYPES, na.rm = TRUE))) %>% 
          dcast(INDIVIDUALS + POP_ID ~ Catalog.ID, value.var = "HAPLOTYPES")
        
      }
    }
    
    # transform the imputed dataset into colony  ------------------------
    message("Imputed haplotypes into factory for conversion into COLONY...")
    
    haplo.imp <- suppressWarnings(
      haplo.imp %>% 
        melt(
          id.vars = c("INDIVIDUALS", "POP_ID"), 
          variable.name = "Catalog.ID", 
          value.name = "HAPLOTYPES"
        ) %>% 
        tidyr::separate(
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
    
    # Allele frequency per locus
    
    if (allele.freq == "overall"){
      
      allele.per.locus <- haplo.imp %>% select(-INDIVIDUALS, -POP_ID, -ALLELE) %>% 
        colwise(nlevels)(.)
      
      frequency.markers <- suppressWarnings(
        haplo.imp %>%
          tidyr::gather(Catalog.ID, NUCLEOTIDES, -c(INDIVIDUALS, POP_ID, ALLELE)) %>% 
          select(-INDIVIDUALS, -ALLELE, -POP_ID) %>% 
          group_by(Catalog.ID) %>%
          filter(NUCLEOTIDES != "NA") %>% 
          mutate(N = length(NUCLEOTIDES)) %>% 
          group_by(Catalog.ID, NUCLEOTIDES) %>%
          mutate(n = n()) %>% 
          distinct(Catalog.ID, NUCLEOTIDES) %>% 
          mutate(FREQ = n/N) %>% 
          select(-N, -n) %>% 
          group_by(Catalog.ID) %>% 
          mutate(NUCLEOTIDES = as.numeric(factor(NUCLEOTIDES))) %>% 
          arrange(Catalog.ID, NUCLEOTIDES) %>% 
          dcast(Catalog.ID ~ NUCLEOTIDES, value.var = "FREQ", drop = FALSE) %>% 
          select(-Catalog.ID) %>% 
          mutate_each(funs(stri_replace_na(str = ., replacement = "")))
      )    
      message("step 3/4: completed")
      
    } else if (missing(allele.freq) != TRUE & allele.freq != FALSE){
      
      allele.per.locus <- haplo.imp %>% 
        filter(POP_ID %in% allele.freq) %>%
        select(-INDIVIDUALS, -POP_ID, -ALLELE) %>% 
        colwise(nlevels)(.)
      
      frequency.markers <- suppressWarnings(
        haplo.imp %>%
          filter(POP_ID %in% allele.freq) %>% 
          tidyr::gather(Catalog.ID, NUCLEOTIDES, -c(INDIVIDUALS, POP_ID, ALLELE)) %>% 
          select(-INDIVIDUALS, -ALLELE, -POP_ID) %>% 
          group_by(Catalog.ID) %>%
          filter(NUCLEOTIDES != "NA") %>% 
          mutate(N = length(NUCLEOTIDES)) %>% 
          group_by(Catalog.ID, NUCLEOTIDES) %>%
          mutate(n = n()) %>% 
          distinct(Catalog.ID, NUCLEOTIDES) %>% 
          mutate(FREQ = n/N) %>% 
          select(-N, -n) %>% 
          group_by(Catalog.ID) %>% 
          mutate(NUCLEOTIDES = as.numeric(factor(NUCLEOTIDES))) %>% 
          arrange(Catalog.ID, NUCLEOTIDES) %>% 
          dcast(Catalog.ID ~ NUCLEOTIDES, value.var = "FREQ", drop = FALSE) %>% 
          select(-Catalog.ID) %>% 
          mutate_each(funs(stri_replace_na(str = ., replacement = "")))
      )
      message("step 3/4: completed")
      
    } else {
      message("step 3/4: completed")
    }
    # This part is different than gtypes...
    haplo.imp <- suppressWarnings(
      haplo.imp %>%
        mutate(GROUP = rep(1, times = nrow(.))) %>% 
        group_by(GROUP) %>% 
        mutate_each(funs(as.integer), -c(ALLELE, INDIVIDUALS, POP_ID, GROUP)) %>%
        ungroup() %>% 
        select(-GROUP) %>% 
        melt(
          id.vars = c("INDIVIDUALS", "POP_ID", "ALLELE"), 
          variable.name = "Catalog.ID", 
          value.name = "HAPLOTYPES"
        ) %>%
        mutate(HAPLOTYPES = as.character(HAPLOTYPES)) %>% 
        mutate(HAPLOTYPES = stri_pad_left(str = HAPLOTYPES, width = 3, pad = "1")) %>% 
        mutate(HAPLOTYPES = stri_replace_na(str = HAPLOTYPES, replacement = "0")) %>%
        group_by(Catalog.ID) %>% 
        arrange(Catalog.ID)
    )
    markers.name <- haplo.imp %>% distinct(Catalog.ID)
    marker.num <- nrow(markers.name)
    markers.name <- t(markers.name)
    
    haplo.imp <- haplo.imp %>% 
      dcast(INDIVIDUALS + POP_ID ~ Catalog.ID + ALLELE, value.var = "HAPLOTYPES")
    
    message("step 4/4: completed")
    
    # results WITH imputations ------------------------------------------------
    res <- haplo.imp %>% select(-POP_ID)
    
    # convert to colony
    # Add "_imputed" to the filename
    filename <- stri_replace_all_fixed(filename,
                                       pattern = ".dat",
                                       replacement = "_imputed.dat")
    
    # Line 1 = Dataset name
    dataset.opt <- "`My first COLONY run`                ! Dataset name"
    dataset.opt <- as.data.frame(dataset.opt)
    write.table(x = dataset.opt, file = filename, sep = " ", append = FALSE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 2 = Output filename
    colony.output.filename <- stri_replace_all_fixed(filename,
                                                     pattern = ".dat",
                                                     replacement = "")
    
    out.name.opt <- paste(colony.output.filename, "                         ! Output file name", sep = "")
    
    out.name.opt <- as.data.frame(out.name.opt)
    write.table(x = out.name.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 3 = Offspring number
    off.num.opt <- nrow(res)
    off.num.opt <- paste(off.num.opt, "                                  ! Number of offspring in the sample", sep = "")
    off.num.opt <- as.data.frame(off.num.opt)
    write.table(x = off.num.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 4 = Number of loci
    marker.number <- (ncol(haplo.imp)-2)/2
    marker.num.opt <- paste(marker.number, "                                 ! Number of loci", sep = "")
    marker.num.opt <- as.data.frame(marker.num.opt)
    write.table(x = marker.num.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 5 = Seed random number generator  
    seed.opt <- "1234                                 ! Seed for random number generator"
    seed.opt <- as.data.frame(seed.opt)
    write.table(x = seed.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 6 = Updating allele frequency
    update.allele.freq.opt <- "0                                    ! 0/1=Not updating/updating allele frequency"
    update.allele.freq.opt <- as.data.frame(update.allele.freq.opt)
    write.table(x = update.allele.freq.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 7 = Dioecious/Monoecious species
    dioecious.opt <- "2                                    ! 2/1=Dioecious/Monoecious species"
    dioecious.opt <- as.data.frame(dioecious.opt)
    write.table(x = dioecious.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 8 = inbreeding
    inbreeding.opt <- paste(inbreeding, "                                    ! 0/1=No inbreeding/inbreeding", sep = "") 
    inbreeding.opt <- as.data.frame(inbreeding.opt)
    write.table(x = inbreeding.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 9 = Ploidy
    ploidy.opt <- "0                                    ! 0/1=Diploid species/HaploDiploid species"
    ploidy.opt <- as.data.frame(ploidy.opt)
    write.table(x = ploidy.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 10 = Mating system (polygamous: 0, monogamous: 1).
    mating.opt <- paste(mating.sys.males,"  ", mating.sys.females, "                                 ! 0/1=Polygamy/Monogamy for males & females", sep = "") 
    mating.opt <- as.data.frame(mating.opt)
    write.table(x = mating.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 11 = Clone inference
    clone.opt <- paste(clone, "                                    ! 0/1=Clone inference =No/Yes", sep = "") 
    clone.opt <- as.data.frame(clone.opt)
    write.table(x = clone.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 12 = Sibship size scaling
    sib.size.scal.opt <- "1                                    ! 0/1=Full sibship size scaling =No/Yes"
    sib.size.scal.opt <- as.data.frame(sib.size.scal.opt)
    write.table(x = sib.size.scal.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 13 = Sibship prior indicator (Integer), average paternal sibship size (Real, optional), average maternal sibship size (Real, optional)
    sib.prior.opt <- "0 0 0                                ! 0, 1, 2, 3 = No, weak, medium, strong sibship size prior; mean paternal & maternal sibship size"
    sib.prior.opt <- as.data.frame(sib.prior.opt)
    write.table(x = sib.prior.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Line 14 = Population allele frequency indicator
    if (missing(allele.freq) | allele.freq == FALSE){
      allele.freq.ind.opt <- "0                                    ! 0/1=Unknown/Known population allele frequency"
      allele.freq.ind.opt <- as.data.frame(allele.freq.ind.opt)
      write.table(x = allele.freq.ind.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    } else{
      allele.freq.ind.opt <- "1                                    ! 0/1=Unknown/Known population allele frequency"
      allele.freq.ind.opt <- as.data.frame(allele.freq.ind.opt)
      write.table(x = allele.freq.ind.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Line 15 and more = Numbers of alleles per locus (Integer, optional). Required when the Population allele frequency indicator is set to 1
    if (allele.freq != FALSE){
      write_delim(x = allele.per.locus, path = filename, delim = " ", append = TRUE, col_names = FALSE)
    }
    
    # Alleles and their frequencies per locus (Integer, Real, optional)    
    if (allele.freq != FALSE){
      write_delim(x = frequency.markers, path = filename, delim = " ", append = TRUE, col_names = FALSE)
    }
    
    # Number of runs
    num.run.opt <- "\n1                                    ! Number of runs"
    num.run.opt <- as.data.frame(num.run.opt)
    write.table(x = num.run.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Length of run, give a value of 1, 2, 3, 4 to indicate short, medium, long, very long run.
    run.length.opt <- paste(run.length, "                                    ! Length of run", sep = "") 
    run.length.opt <- as.data.frame(run.length.opt)
    write.table(x = run.length.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Monitor method (Time in second)
    monitor.met.opt <- "0                                    ! 0/1=Monitor method by Iterate"
    monitor.met.opt <- as.data.frame(monitor.met.opt)
    write.table(x = monitor.met.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    
    # Monitor interval (Time in second)
    monitor.int.opt <- "10000                                ! Monitor interval in Iterate"
    monitor.int.opt <- as.data.frame(monitor.int.opt)
    write.table(x = monitor.int.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # WindowsGUI/DOS, 0 when running Colony in DOS mode or on other platforms
    windows.gui.opt <- "0                                    ! non-Windows version"
    windows.gui.opt <- as.data.frame(windows.gui.opt)
    write.table(x = windows.gui.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Analysis method : 0, 1 or 2 for Pairwise-Likelihood score (PLS), full likelihood method (FL), or the FL and PLS combined method (FPLS). More on these methods are explained above in the Windows GUI data input section.
    analysis.opt <- paste(analysis, "                                    ! Analysis 0 (Pairwise-Likelihood Score), 1 (Full Likelihood), 2 (combined Pairwise-Likelihood Score and Full Likelihood)", sep = "") 
    analysis.opt <- as.data.frame(analysis.opt)
    write.table(x = analysis.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Precision
    precision.opt <- "3                                    ! 1/2/3=low/medium/high Precision for Full likelihood\n"
    precision.opt <- as.data.frame(precision.opt)
    write.table(x = precision.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    
    # insert a space
    # cat("", sep = "\n", fill = TRUE, file = filename, append = TRUE)
    
    # Marker IDs/Names
    # snp.id <- seq(from = 1, to = marker.number, by = 1)
    markers.name.opt <- c(markers.name, "        ! Marker IDs")
    write.table(x = t(as.data.frame(markers.name.opt)), file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Marker types 
    marker.type.opt <- c(rep(0, marker.number),"         ! Marker types, 0/1=Codominant/Dominant") # marker.type (codominant/dominant)
    write.table(x = t(as.data.frame(marker.type.opt)), file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Allelic dropout rates
    dropout <- c(rep(allelic.dropout, marker.number), "        ! Allelic dropout rate at each locus") # allelic dropout rate
    write.table(x = t(as.data.frame(dropout)), file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Error rates
    error <- c(rep(error.rate, marker.number), "        ! False allele rate\n\n") # genotyping error rate
    write.table(x = t(as.data.frame(error)), file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Offspring IDs and genotype
    # write_delim(x = res, path = filename, delim = " ", append = TRUE, col_names = FALSE)
    write.table(x = res, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    message("Including genotypes of individuals in COLONY format")
    
    
    # Probabilities that the father and mother of an offspring are included in the candidate males and females. The two numbers must be provided even if there are no candidate males or/and females.
    prob.opt <- "\n\n0  0                                 ! Prob. of dad/mum included in the candidates"
    prob.opt <- as.data.frame(prob.opt)
    write.table(x = prob.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Numbers of candidate males and females
    candidate.opt <- "0 0                                  ! Numbers of candidate males & females"
    candidate.opt <- as.data.frame(candidate.opt)
    write.table(x = candidate.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Candidate male IDs/names and genotypes
    # cat(male.genotype, sep = "\n", file = colony.filename, append = TRUE)
    if (print.all.colony.opt != FALSE){
      candidate.male.id.geno.opt <- "!Candidate male ID and genotypes"
      candidate.male.id.geno.opt <- as.data.frame(candidate.male.id.geno.opt)
      write.table(x = candidate.male.id.geno.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Candidate female IDs/names and genotypes 
    # cat(female.genotype, sep = "\n", file = colony.filename, append = TRUE)
    if (print.all.colony.opt != FALSE){
      candidate.female.id.geno.opt <- "!Candidate female ID and genotypes"
      candidate.female.id.geno.opt <- as.data.frame(candidate.female.id.geno.opt)
      write.table(x = candidate.female.id.geno.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Number of offspring with known paternity
    known.paternity.opt <- "0                                    ! Number of offspring with known father"
    known.paternity.opt <- as.data.frame(known.paternity.opt)
    write.table(x = known.paternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Known offspring-father dyad
    if (print.all.colony.opt != FALSE){
      known.father.dyad.opt <- "! Offspring ID and known father ID (Known offspring-father dyad)"
      known.father.dyad.opt <- as.data.frame(known.father.dyad.opt)
      write.table(x = known.father.dyad.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Number of offspring with known maternity
    known.maternity.opt <- "0                                    ! Number of offspring with known mother"
    known.maternity.opt <- as.data.frame(known.maternity.opt)
    write.table(x = known.maternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Known offspring-mother dyad
    if (print.all.colony.opt != FALSE){
      known.mother.dyad.opt <- "! Offspring ID and known mother ID (Known offspring-mother dyad)"
      known.mother.dyad.opt <- as.data.frame(known.mother.dyad.opt)
      write.table(x = known.mother.dyad.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Number of known paternal sibships (Integer)
    known.paternal.sibships.opt <- "0                                    ! Number of known paternal sibships"
    known.paternal.sibships.opt <- as.data.frame(known.paternal.sibships.opt)
    write.table(x = known.paternal.sibships.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    
    # Paternal sibship size and members (Integer, String, optional).
    if (print.all.colony.opt != FALSE){
      known.paternal.sibships.size.opt <- "! Paternal sibship size and members"
      known.paternal.sibships.size.opt <- as.data.frame(known.paternal.sibships.size.opt)
      write.table(x = known.paternal.sibships.size.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Number of known maternal sibships (Integer)
    known.maternal.sibships.opt <- "0                                    ! Number of known maternal sibships"
    known.maternal.sibships.opt <- as.data.frame(known.maternal.sibships.opt)
    write.table(x = known.maternal.sibships.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    
    # Maternal sibship size and members (Integer, String, optional).
    if (print.all.colony.opt != FALSE){
      known.maternal.sibships.size.opt <- "! Maternal sibship size and members"
      known.maternal.sibships.size.opt <- as.data.frame(known.maternal.sibships.size.opt)
      write.table(x = known.maternal.sibships.size.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Number of offspring with known excluded paternity (Integer). 
    offspring.known.excl.paternity.opt <- "0                                    ! Number of offspring with known excluded fathers"
    offspring.known.excl.paternity.opt <- as.data.frame(offspring.known.excl.paternity.opt)
    write.table(x = offspring.known.excl.paternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    
    # Excluded paternity 
    if (print.all.colony.opt != FALSE){
      excl.paternity.opt <- "! Offspring ID, number of excluded fathers, and excluded father IDs"
      excl.paternity.opt <- as.data.frame(excl.paternity.opt)
      write.table(x = excl.paternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    
    # Number of offspring with known excluded maternity (Integer). 
    offspring.known.excl.maternity.opt <- "0                                    ! Number of offspring with known excluded mothers"
    offspring.known.excl.maternity.opt <- as.data.frame(offspring.known.excl.maternity.opt)
    write.table(x = offspring.known.excl.maternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    
    # Excluded maternity 
    if (print.all.colony.opt != FALSE){
      excl.maternity.opt <- "! Offspring ID, number of excluded mothers, and excluded father IDs"
      excl.maternity.opt <- as.data.frame(excl.maternity.opt)
      write.table(x = excl.maternity.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Number of offspring with known excluded paternal sibships
    offspring.known.excl.paternal.sibships.opt <- "0                                    ! Number of offspring with known excluded paternal sibships"
    offspring.known.excl.paternal.sibships.opt <- as.data.frame(offspring.known.excl.paternal.sibships.opt)
    write.table(x = offspring.known.excl.paternal.sibships.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Excluded paternal siblings
    if (print.all.colony.opt != FALSE){
      excluded.paternal.siblings.opt <- "! Excluded paternal siblings"
      excluded.paternal.siblings.opt <- as.data.frame(excluded.paternal.siblings.opt)
      write.table(x = excluded.paternal.siblings.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Number of offspring with known excluded maternal sibships
    offspring.known.excl.maternal.sibships.opt <- "0                                    ! Number of offspring with known excluded maternal sibships"
    offspring.known.excl.maternal.sibships.opt <- as.data.frame(offspring.known.excl.maternal.sibships.opt)
    write.table(x = offspring.known.excl.maternal.sibships.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Excluded maternal siblings
    if (print.all.colony.opt != FALSE){
      excluded.maternal.siblings.opt <- "! Excluded maternal siblings"
      excluded.maternal.siblings.opt <- as.data.frame(excluded.maternal.siblings.opt)
      write.table(x = excluded.maternal.siblings.opt, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }   
    message.write.colony <- paste("A COLONY input file project was created and saved in your working directory", filename, sep = " : ")
  }
}



