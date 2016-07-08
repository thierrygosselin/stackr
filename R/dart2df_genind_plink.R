# Import, filter and transform a dart output file to different formats

#' @name dart2df_genind_plink

#' @title swiss army knife tool to prepare [DArT](http://www.diversityarrays.com) output file for population genetics analysis.

#' @description Import, filter and transform a DArT output file to different formats: 
#' data frame of genotypes, genind object and/or PLINK \code{tped/tfam} format.

#' @param data DArT output file in wide format or binary format tipically 
#' used by CSIRO genomic projects.

#' @param strata A tab delimited file with columns header:
#' \code{INDIVIDUALS} and \code{POP_ID}. 
#' Note: the column \code{POP_ID} refers to any grouping of individuals. If a third column
#' named \code{NEW_ID} is used, this column will be used to replace the
#' \code{INDIVIDUALS} column in the main data file.

#' @param pop.levels (optional string) A character string with your populations ordered.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

#' @param pop.select (optional string) Keep specific populations. 
#' Default = \code{NULL} for no selection and keep all population.
#' e.g. \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population samples.
#' 
#' @param filter.monomorphic (optional, logical) Should the monomorphic markers present in 
#' the dataset be filtered out ? Default: \code{filter.monomorphic = TRUE}.

#' @param common.markers (optional, logical) Default: \code{common.markers = TRUE}, 
#' will only keep markers in common (genotyped) between all the populations.

#' @param filter.ind.heterozygosity (optional, numerical) Filter the 
#' individual heterozygosity (averaged accross markers). This argument is useful
#' to detect mixed samples. The individuals to remove will standout in the plot.
#' Default: \code{filter.ind.heterozygosity = NULL}.
#' e.g to remove individuals with heterozygosity higher than 10% use:
#' \code{filter.ind.heterozygosity = 0.1}.

#' @param plot.ind.heterozygosity (optional, logical) Plot the distribution 
#' (boxplot and jitter plot) of individual heterozygosity. 
#' Default: \code{plot.ind.heterozygosity = FALSE}.

#' @param filter.reproducibility (optional, numerical) Filter the \code{RepAvg} 
#' column in the data set. Default: \code{filter.reproducibility = NULL}.
#' e.g to keep markers with reproducibility >= 99%,
#' use: \code{filter.reproducibility = 0.99}.

#' @param plot.reproducibility (optional, logical) Plot the distribution 
#' (violin plot) of reproducibility. Default: \code{plot.reproducibility = FALSE}.

#' @param filter.coverage.high (optional, numerical) Filter the upper bound of the
#' \code{AvgCountSnp} column in the data set. Default: \code{filter.coverage.high = NULL}.
#' e.g to keep markers with coverage <= 150 (depth of coverage),
#' use : \code{filter.coverage.high = 150}.

#' @param filter.coverage.low (optional, numerical) Filter the lower bound of the 
#' \code{AvgCountSnp} column in the data set. Default: \code{filter.coverage.low = NULL}.
#' e.g to keep markers with coverage >= 10 (depth of coverage),
#' use : \code{filter.coverage.low = 10}.
#' 
#' @param plot.coverage (optional, logical) Plot the coverage distribution 
#' (violin plot). Default: \code{plot.coverage = FALSE}.

#' @param filter.call.rate (optional, numerical) Filter the \code{CallRate} 
#' column in the data set. Default: \code{filter.call.rate = NULL}. e.g to keep 
#' markers genotyped in more than 95% of the individuals use :
#' \code{filter.call.rate = 0.95}

#' @param plot.call.rate (optional, logical) Plot the distribution 
#' (violin plot) of call rate. Default: \code{plot.call.rate = FALSE}.

#' @param filter.ind.missing.geno (optional, numerical) Filter with the individual's
#' proportion of missing genotype. Similar to call rate.
#' e.g to keep individuals genotyped at >= 0.90 of the markers, use:
#' \code{filter.ind.missing.geno = 0.90}.
#' Default: \code{filter.ind.missing.geno = NULL}.
#' @param plot.ind.missing.geno (optional, logical) Plot the distribution 
#' (violin plot) of an individual's genotype proportion. 
#' Default: \code{plot.ind.missing.geno = FALSE}.

#' @param filter.markers.missing.ind (optional, numerical) Filter markers based 
#' on the proportion of individuals genotyped.
#' e.g to keep markers with >= 0.95 of the genotyped individuals, use:
#' \code{filter.markers.missing.ind = 0.95}.
#' Default: \code{filter.markers.missing.ind = NULL}.
#' @param plot.markers.missing.ind (optional, logical) Plot the distribution 
#' (violin plot) of markers genotype individuals proportion. 
#' Default: \code{plot.markers.missing.ind = FALSE}.

#' @param filter.snp.ld (optional, character) With anonymous markers from
#' reduce representation library like RADseq/GBS/DArT de novo discovery, 
#' you can explore the impact of the number of snp on the read (haplotype/locus).
#' This argument is used to minimize linkage disequilibrium (LD) on short reads,
#' by choosing among these options: 
#' \code{"1snp"} to keep only reads with 1 SNP/reads,
#' \code{"2snp"} to keep only reads with at most 2 SNP/reads,
#' \code{"random"} to select 1 SNP/reads randomly, 
#' \code{"first"} to select the first SNP on all reads, 
#' \code{"last"} to select the last SNP on all reads 
#' (last SNP are usually associated to higher error rate...). 
#' Default: \code{snp.ld = NULL}.
#' Note: for long linkage detection use PLINK linkage disequilibrium based SNP 
#' pruning.

#' @param plot.number.snp.reads (optional, logical) Plot the distribution of SNP
#' per read. 
#' Default: \code{plot.number.snp.reads = FALSE}.




#' @param maf.thresholds (string, double, optional) String with 
#' local/populations and global/overall maf thresholds, respectively.
#' Default: \code{maf.thresholds = NULL}. 
#' e.g. \code{maf.thresholds = c(0.05, 0.1)} for a local maf threshold 
#' of 0.05 and a global threshold of 0.1.

#' @param maf.pop.num.threshold (integer, optional) When maf thresholds are used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. Default: \code{maf.pop.num.threshold = 1}

#' @param maf.operator (character, optional) \code{maf.operator = "AND"} or default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?


#' @param plink (optional, logical). To have the filtered data set output as 
#' plink \code{tped/tfam} file, use \code{plink = TRUE}. 
#' Default: \code{plink = FALSE}

#' @param genind (optional, logical). To have the filtered data set output as 
#' genind object to use inside adegenet, use \code{genind = TRUE}. 
#' Default: \code{genind = FALSE}

#' @param genepop (optional, logical). To have the filtered data set output as 
#' a genepop file, use \code{genepop = TRUE}. 
#' Default: \code{genepop = FALSE}.

#' @param structure (optional, logical). To have the filtered data set output as 
#' a structure file, use \code{structure = TRUE}. 
#' Default: \code{structure = FALSE}.

#' @param filename (optional) The name of the file written to the directory.
#' No file extension at the end. Default: \code{filename = NULL}

#' @param imputation.method Should a map-independent imputation of markers be
#' computed. Available choices are: (1) \code{FALSE} for no imputation.
#' (2) \code{"max"} to use the most frequent category for imputations.
#' (3) \code{"rf"} using Random Forest algorithm. 
#' Default: \code{imputation.method = NULL}.
#' @param impute (character) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
#' @param imputations.group \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by population. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default = \code{"populations"}.
#' @param num.tree The number of trees to grow in Random Forest. Default is 100.
#' @param iteration.rf The number of iterations of missing data algorithm
#' in Random Forest. Default is 10.
#' @param split.number Non-negative integer value used to specify
#' random splitting in Random Forest. Default is 100.
#' @param verbose Logical. Should trace output be enabled on each iteration
#' in Random Forest ? Default is \code{FALSE}.
#' @param parallel.core (optional) The number of cores for OpenMP shared-memory parallel
#' programming of Random Forest imputations. For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.

#' @return Depending on arguments selected, several files are written to the your
#' working directory or \code{folder}
#' The output in your global environment is a list. To view the assignment results
#' \code{$assignment} to view the ggplot2 figure \code{$plot.assignment}. 
#' See example below.

#' @export
#' @rdname dart2df_genind_plink
#' @import dplyr
#' @import parallel
#' @import stringi
#' @import adegenet
#' @importFrom purrr discard
#' @importFrom data.table fread

#' @examples
#' \dontrun{
#' testing
#' }

#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, et al. (2007)
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. 
#' American Journal of Human Genetics. 81: 559–575. doi:10.1086/519795
#' @references Jombart T, Devillard S, Balloux F. (2010)
#' Discriminant analysis of principal components: a new method for the analysis 
#' of genetically structured populations.
#' BMC Genet. 11: 94. doi:10.1186/1471-2156-11-94.
#' @references Jombart T, Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis 
#' of genome-wide SNP data. 
#' Bioinformatics. 27: 3070-3071. doi:10.1093/bioinformatics/btr521
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("ID", "CloneID", "SnpPosition", "CallRate", "AvgCountRef", "AvgCountSnp", 
      "RepAvg", "NOT_USEFUL", "SNP", "CALL_RATE", "AVG_COUNT_REF", 
      "AVG_COUNT_SNP", "REP_AVG", "NEW_ID", "SNP_N", "ALLELE_NAME", "ALLELE_NUMBER",
      "ALLELES_COUNT", "GENOTYPED_PROP", "MISSING_IND_PROP", "AlleleID", "HET_NUMBER",
      "HET_PERCENT", "HET_PROP"
    )
  )
}

dart2df_genind_plink <- function(data,
                                 strata,
                                 pop.levels = NULL,
                                 blacklist.id = NULL,
                                 pop.select = NULL,
                                 filter.monomorphic = TRUE,
                                 common.markers = TRUE,
                                 filter.ind.heterozygosity = NULL,
                                 plot.ind.heterozygosity = FALSE,
                                 filter.reproducibility = NULL,
                                 plot.reproducibility = FALSE,
                                 filter.coverage.high = NULL,
                                 filter.coverage.low = NULL,
                                 plot.coverage = FALSE,
                                 filter.call.rate = NULL,
                                 plot.call.rate = FALSE,
                                 filter.ind.missing.geno = NULL,
                                 plot.ind.missing.geno = FALSE,
                                 filter.markers.missing.ind = NULL,
                                 plot.markers.missing.ind = FALSE,
                                 filter.snp.ld = NULL,
                                 plot.number.snp.reads = FALSE,
                                 maf.thresholds = NULL,
                                 maf.pop.num.threshold = 1,
                                 maf.operator = "OR",
                                 plink = FALSE,
                                 genind = FALSE,
                                 genepop = FALSE,
                                 structure = FALSE,
                                 filename = NULL,
                                 imputation.method = NULL,
                                 impute = "genotype",
                                 imputations.group = "populations",
                                 num.tree = 100,
                                 iteration.rf = 10,
                                 split.number = 100,
                                 verbose = FALSE,
                                 parallel.core = detectCores()-1,
                                 ...) {
  
  cat("#######################################################################\n")
  cat("#################### stackr: dart2df_genind_plink #####################\n")
  cat("#######################################################################\n")
  message("Importing data ...")
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (missing(strata)) stop("strata file missing")
  
  # Filename ------------------------------------------------------------------
  # Create a folder based on filename to save the output files *****************
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename <- stri_join("stackr_default_filename", file.date, sep = "_")
    file.date <- NULL #unused object
  }
  # Strata file ------------------------------------------------------------------
  strata.df <- read_tsv(file = strata, col_names = TRUE)
  
  # Import data ---------------------------------------------------------------
  colnames.keeper <- c(c("AlleleID", "SNP", "SnpPosition", "CallRate", "AvgCountRef", "AvgCountSnp", "RepAvg"), strata.df$INDIVIDUALS)
  
  input <- suppressWarnings(
    data.table::fread(
      input = data,
      sep = "\t",
      stringsAsFactors = FALSE, 
      header = TRUE,
      na.strings = "-",
      strip.white = TRUE,
      select = colnames.keeper,
      # drop = c("AlleleID", "ClusterTempIndex", "ClusterSize", "AlleleSeqDist", "SnpPosition", "CallRate", "OneRatioRef","OneRatioSnp", "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", "AvgPIC", "AvgCountRef", "AvgCountSnp", "RatioAvgCountRefAvgCountSnp", "FreqHetsMinusFreqMinHom", "AlleleCountsCorrelation", "aggregateTagsTotal", "DerivedCorrMinusSeedCorr", "RepRef", "RepSNP", "RepAvg", "PicRepRef", "PicRepSNP", "TotalPicRepRefTest", "TotalPicRepSnpTest", "AlleleSequence", "AlleleSequenceREF", "AlleleSequenceSNP"),
      # drop = c("AlleleID", "AlleleSequence", "AlleleSequenceREF", "AlleleSequenceSNP", "OneRatioRef", "OneRatioSnp", "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", "AvgPIC"),
      showProgress = TRUE,
      verbose = FALSE
    ) %>% 
      as_data_frame() %>%
      rename(LOCUS = AlleleID, POS = SnpPosition, CALL_RATE = CallRate, AVG_COUNT_REF = AvgCountRef, AVG_COUNT_SNP = AvgCountSnp, REP_AVG = RepAvg) %>% 
      arrange(LOCUS, POS)
  )  
  
  # Determine the type of DArT file
  binary <- anyDuplicated(input$LOCUS)
  
  # Screen for duplicate names -------------------------------------------------
  remove.list <- c("LOCUS", "SNP", "POS", "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
  individuals.df <- data_frame(INDIVIDUALS = purrr::discard(.x = colnames(input), .p = colnames(input) %in% remove.list))
  duplicate.individuals <- length(individuals.df$INDIVIDUALS) - n_distinct(individuals.df$INDIVIDUALS)
  if (duplicate.individuals == 0) {
    message("Duplicate individual names in the data: no")
  } else {
    stop(stri_paste("Duplicated individuals names found in the data set.\nNumber of duplicate names = ", duplicate.individuals))
  }
  # removing unused object
  remove.list <- NULL
  individuals.df <- NULL
  duplicate.individuals <- NULL
  
  # Tidying data ---------------------------------------------------------------
  if (binary != 2) {
    message("Tidying the dataset...")
    input <- input %>% 
      tidyr::separate(col = LOCUS, into = c("LOCUS", "NOT_USEFUL"), sep = "\\|", extra = "drop") %>%  
      select(-NOT_USEFUL) %>%
      tidyr::separate(col = SNP, into = c("NOT_USEFUL", "KEEPER"), sep = ":", extra = "drop") %>% 
      select(-NOT_USEFUL) %>%
      tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>% 
      mutate(MARKERS = stri_paste(LOCUS, POS, sep = "_"))
    
    input <- data.table::melt.data.table(
      data = as.data.table(input), 
      id.vars = c("MARKERS", "LOCUS", "POS", "REF", "ALT", "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG"), 
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = "GT"
    ) %>% 
      as_data_frame() %>% 
      mutate(
        REF = stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("01", "02", "03", "04"), vectorize_all = FALSE), # replace nucleotide with numbers
        ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("01", "02", "03", "04"), vectorize_all = FALSE), # replace nucleotide with numbers
        GT = as.character(GT),
        GT = stri_replace_all_fixed(str = GT, pattern = c("0", "1", "2"), replacement = c("REF_REF", "ALT_ALT", "REF_ALT"), vectorize_all = FALSE),
        GT = stri_replace_na(str= GT, replacement = "0_0")
      )
  }
  if (binary == 2) {
    message("Tidying DArT binary data set")
    input <- suppressWarnings(
      input %>% 
        tidyr::separate(col = SNP, into = c("NOT_USEFUL", "KEEPER"), sep = ":", extra = "drop") %>% 
        select(-NOT_USEFUL) %>%
        tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">", extra = "drop") %>%
        tidyr::unite(MARKERS, LOCUS, POS, sep = "_", remove = FALSE )
    )
    
    # necessary to deal with the duplication of lines because of the GT in 2 lines
    grouping.column <- input %>% 
      ungroup() %>% 
      select(MARKERS, LOCUS, POS, REF, ALT, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG) %>% 
      filter(!is.na(REF) | !is.na(ALT)) %>% 
      distinct(MARKERS, LOCUS, POS, REF, ALT, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG, .keep_all = TRUE) %>% 
      mutate(
        REF = stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("01", "02", "03", "04"), vectorize_all = FALSE), # replace nucleotide with numbers
        ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("01", "02", "03", "04"), vectorize_all = FALSE)# replace nucleotide with numbers
      ) 
    
    # Data tidying
    input <- input %>% 
      select(-c(LOCUS, POS, REF, ALT, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG)) %>% 
      mutate(
        ALLELE_NAME = rep("A", n()),
        ALLELE_NUMBER = rep(1:2, each = 1, times = n()/2)
      ) %>% 
      tidyr::unite(ALLELE, c(ALLELE_NAME, ALLELE_NUMBER), sep = "") %>%   
      tidyr::gather(INDIVIDUALS, GENOTYPE, -c(MARKERS, ALLELE)) %>% 
      tidyr::spread(data = ., ALLELE, GENOTYPE) %>% 
      tidyr::unite(GENOTYPE, c(A1, A2)) %>% 
      inner_join(grouping.column, by = c("MARKERS")) %>% 
      mutate(
        GT = ifelse(GENOTYPE == "1_0", "REF_REF",
                    ifelse(GENOTYPE == "0_1", "ALT_ALT", 
                           ifelse(GENOTYPE == "1_1", "REF_ALT", "0_0")))
      ) %>% 
      select(MARKERS, LOCUS, POS, REF, ALT, INDIVIDUALS, GT, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG)
    
    grouping.column <- NULL # remove unused object
  }
  
  # get the number of markers before filters
  locus.before.filters <- n_distinct(input$LOCUS)
  snp.before.filters <- n_distinct(input$MARKERS)
  
  # Strata file ----------------------------------------------------------------
  input <- input %>%
    left_join(strata.df, by = "INDIVIDUALS")
  
  if (ncol(strata.df) == 3) {
    input <- input %>%
      select(-INDIVIDUALS) %>% 
      rename(INDIVIDUALS = NEW_ID)
  }
  
  
  # pop.levels -------------------------------------------------------------------
  if (is.null(pop.levels)) {
    input <- input %>%
      mutate(POP_ID = factor(POP_ID))
  } else {
    input <- input %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered =TRUE))
  }
  
  # Import blacklist id ********************************************************
  if (is.null(blacklist.id)) { # No blacklist of ID
    message("Blacklisted individuals: no")
  } else { # With blacklist of ID
    blacklist.id <- read_tsv(blacklist.id, col_names = TRUE)
    message(stri_paste("Blacklisted individuals: yes (", length(blacklist.id$INDIVIDUALS), " ind.)"))
    message("Filtering with blacklist of individuals")
    input <- suppressWarnings(anti_join(input, blacklist.id, by = "INDIVIDUALS"))
  }
  
  # pop.select -----------------------------------------------------------------
  if (!is.null(pop.select)) {
    message(stri_join(length(pop.select), "population(s) selected", sep = " "))
    input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
  }
  
  # Filter monomorphic markers  ---------------------------------------------------
  if (filter.monomorphic == TRUE) {
    # screen for monomorphic
    blacklist.monomorphic  <- input %>%
      select(MARKERS, GT) %>%
      filter(GT != "0_0") %>%
      distinct(MARKERS, GT) %>% 
      group_by(MARKERS) %>%
      tally %>% 
      filter(n == 1) %>% 
      select(MARKERS) %>% 
      arrange(MARKERS)
    
    # Remove the markers from the dataset
    message(stri_paste("Filter monomorphic: ", n_distinct(blacklist.monomorphic$MARKERS), " markers deleted"))
    if (length(blacklist.monomorphic$MARKERS > 0)) {
      write_tsv(x = blacklist.monomorphic, path = "blacklist.monomorphic.tsv", col_names = TRUE)
      input <- anti_join(input, blacklist.monomorphic, by = "MARKERS")
    }
  }
  
  # Filter Heterozygote/mixed individuals diagnostic ---------------------------
  if (!is.null(filter.ind.heterozygosity)) {
    # filter.ind.heterozygosity <- 0.06 # test
    
    # Create a new df with heterozygote info
    het.ind <- input %>% 
      select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
      group_by(INDIVIDUALS, POP_ID) %>% 
      summarise(
        GENOTYPED = length(GT[GT != "0_0"]),
        HET_NUMBER = length(GT[GT == "REF_ALT"]),
        HET_PROP = HET_NUMBER/GENOTYPED,
        HET_PERCENT = round(HET_PROP*100, 2)
      )
    
    blacklist.individuals  <- het.ind %>%
      filter(HET_PROP > filter.ind.heterozygosity) %>%
      distinct(INDIVIDUALS)
    
    # Remove the individuals from the dataset
    message(stri_paste("Filter individual heterozygosity: ", length(blacklist.individuals$INDIVIDUALS), " individual(s) deleted"))
    if (length(blacklist.individuals$INDIVIDUALS > 0)) {
      write_tsv(x = blacklist.individuals, path = "blacklist.individuals.heterozygosity.tsv", col_names = TRUE)
      input <- anti_join(input, blacklist.individuals, by = "INDIVIDUALS")
    }
    
    
    if (plot.ind.heterozygosity) {
      jitterplot.ind.het <- ggplot(data = het.ind, aes(x = POP_ID, y = HET_PERCENT, colour = POP_ID)) + 
        geom_jitter() + 
        labs(y = "Mean Heterozygosity (percent)")+
        labs(x = "Individuals")+
        labs(colour = "Populations") +
        # theme_bw()+
        theme(
          legend.position = "none",
          # panel.grid.minor.x = element_blank(), 
          # panel.grid.major.y = element_blank(), 
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
          axis.text.x = element_text(size = 10, family = "Helvetica"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
          axis.text.y = element_text(size = 8, family = "Helvetica")
        )
      
      boxplot.ind.het <- ggplot(data = het.ind, aes(x = POP_ID, y = HET_PERCENT, colour = POP_ID)) + 
        geom_boxplot() + 
        labs(y = "Mean Heterozygosity (percent)")+
        labs(x = "Individuals")+
        labs(colour = "Populations") +
        # theme_bw()+
        theme(
          legend.position = "none",
          # panel.grid.minor.x = element_blank(), 
          # panel.grid.major.y = element_blank(), 
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
          axis.text.x = element_text(size = 10, family = "Helvetica"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
          axis.text.y = element_text(size = 8, family = "Helvetica")
        )
      
    } 
    if (!(plot.ind.heterozygosity)) {
      plot.ind.heterozygosity <- "not selected"
    }
  }
  
  if (is.null(filter.ind.heterozygosity) & (plot.ind.heterozygosity)) {
    # Create a new df with heterozygote info
    het.ind <- input %>% 
      select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
      group_by(INDIVIDUALS, POP_ID) %>% 
      summarise(
        GENOTYPED = length(GT[GT != "0_0"]),
        HET_NUMBER = length(GT[GT == "REF_ALT"]),
        HET_PROP = HET_NUMBER/GENOTYPED,
        HET_PERCENT = round(HET_PROP*100, 2)
      )
    
    jitterplot.ind.het <- ggplot(data = het.ind, aes(x = POP_ID, y = HET_PERCENT, colour = POP_ID)) + 
      geom_jitter() + 
      labs(y = "Mean Heterozygosity (percent)")+
      labs(x = "Individuals")+
      labs(colour = "Populations") +
      # theme_bw()+
      theme(
        legend.position = "none",
        # panel.grid.minor.x = element_blank(), 
        # panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = element_text(size = 10, family = "Helvetica"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = element_text(size = 8, family = "Helvetica")
      )
    
    boxplot.ind.het <- ggplot(data = het.ind, aes(x = POP_ID, y = HET_PERCENT, colour = POP_ID)) + 
      geom_boxplot() + 
      labs(y = "Mean Heterozygosity (percent)")+
      labs(x = "Individuals")+
      labs(colour = "Populations") +
      # theme_bw()+
      theme(
        legend.position = "none",
        # panel.grid.minor.x = element_blank(), 
        # panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = element_text(size = 10, family = "Helvetica"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = element_text(size = 8, family = "Helvetica")
      )
    
  }
  
  # Filter common markers between all populations  -------------------------------
  if (common.markers == TRUE) {
    # get the pop number
    pop.number <- n_distinct(input$POP_ID)
    
    # filter
    filter <- input %>% 
      ungroup() %>% 
      filter(GT != "0_0") %>%
      distinct(MARKERS, POP_ID) %>% 
      group_by(MARKERS) %>%
      tally() %>% 
      filter(n == pop.number) %>% 
      arrange(MARKERS) %>%
      distinct(MARKERS)
    
    whitelist.filter <- filter %>% distinct(MARKERS)
    
    blacklist.common.markers <- input %>% 
      distinct(MARKERS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter common markers: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers deleted"))
    if (length(blacklist.common.markers$MARKERS > 0)) {
      write_tsv(x = blacklist.common.markers, path = "blacklist.common.markers.tsv", col_names = TRUE)
      input <- suppressWarnings(input %>% semi_join(filter, by = "MARKERS"))
    }
  }
  
  # Filtering reproducibility  ---------------------------------------------------
  if (!is.null(filter.reproducibility)) {
    filter <- input %>% 
      filter(REP_AVG >= filter.reproducibility)
    
    whitelist.filter <- filter %>% select(MARKERS, LOCUS, POS) %>% distinct(MARKERS, LOCUS, POS)
    
    blacklist.markers.reproducibility <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter reproducibility: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers deleted"))
    if (length(blacklist.markers.reproducibility$MARKERS > 0)) {
      write_tsv(x = blacklist.markers.reproducibility, path = "blacklist.markers.reproducibility.tsv", col_names = TRUE)
    }
    if (plot.reproducibility) {
      data.combined <- bind_rows(
        data.before <- input %>% 
          select(POP_ID, INDIVIDUALS, REP_AVG) %>% 
          mutate(GROUP = rep("before", n())),
        data.after <- filter %>% 
          select(POP_ID, INDIVIDUALS, REP_AVG) %>% 
          mutate(GROUP = rep("after", n()))
      ) %>% 
        mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
      
      plot.repro <- ggplot(data.combined, aes(x = factor(POP_ID), y = REP_AVG, na.rm = TRUE))+
        geom_violin(trim = TRUE)+
        geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
        stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
        labs(x = "Sampling sites")+
        labs(y = "Markers reproducibility average")+
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
          legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )+
        facet_grid(~GROUP)
    } 
    
    if (!(plot.reproducibility)) {
      plot.repro <- "not selected"
    }
    input <- filter
  }
  
  if (is.null(filter.reproducibility) & (plot.reproducibility)) {
    data.combined <- input %>% 
      select(POP_ID, INDIVIDUALS, REP_AVG) %>% 
      mutate(GROUP = rep("before filter", n()))
    
    plot.repro <- ggplot(data.combined, aes(x = factor(POP_ID), y = REP_AVG, na.rm = TRUE))+
      geom_violin(trim = TRUE)+
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
      labs(x = "Sampling sites")+
      labs(y = "Markers reproducibility average before filter")+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
      )+
      facet_grid(~GROUP)
  } 
  
  # Filtering coverage --------------------------------------------------------
  input.before.filter <- input
  
  # high bound
  if (!is.null(filter.coverage.high)) {
    filter <- input %>% 
      filter(AVG_COUNT_SNP <= filter.coverage.high)
    
    whitelist.filter <- filter %>% distinct(MARKERS, LOCUS, POS)
    
    blacklist.markers.coverage.high <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter coverage high: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers deleted"))
    if (length(blacklist.markers.coverage.high$MARKERS > 0)) {
      write_tsv(x = blacklist.markers.coverage.high, path = "blacklist.markers.coverage.high.tsv", col_names = TRUE)
      input <- filter
    }
  }
  # lower bound
  if (!is.null(filter.coverage.low)) {
    filter <- input %>% 
      filter(AVG_COUNT_SNP >= filter.coverage.low)
    
    whitelist.filter <- filter %>% distinct(MARKERS, LOCUS, POS)
    
    blacklist.markers.coverage.low <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter coverage low: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers deleted"))
    if (length(blacklist.markers.coverage.low$MARKERS > 0)) {
      write_tsv(x = blacklist.markers.coverage.low, path = "blacklist.markers.coverage.low.tsv", col_names = TRUE)
      input <- filter
    }
  }
  
  if (!is.null(filter.coverage.high) | !is.null(filter.coverage.low) & (plot.coverage)) {
    data.combined <- bind_rows(
      data.before <- input.before.filter %>% 
        select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>% 
        mutate(GROUP = rep("before", n())),
      data.after <- input %>% 
        select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>% 
        mutate(GROUP = rep("after", n()))
    ) %>% 
      mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
    
    plot.cov <- ggplot(data.combined, aes(x = factor(POP_ID), y = AVG_COUNT_SNP, na.rm = TRUE))+
      geom_violin(trim = TRUE)+
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
      labs(x = "Sampling sites")+
      labs(y = "Markers coverage")+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
      )+
      facet_grid(~GROUP)
  } 
  if (is.null(filter.coverage.high) & is.null(filter.coverage.low) & (plot.coverage)) {
    data.combined <- input.before.filter %>% 
      select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>% 
      mutate(GROUP = rep("before filter", n()))
    
    plot.cov <- ggplot(data.combined, aes(x = factor(POP_ID), y = AVG_COUNT_SNP, na.rm = TRUE))+
      geom_violin(trim = TRUE)+
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
      labs(x = "Sampling sites")+
      labs(y = "Markers coverage before filter")+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
      )+
      facet_grid(~GROUP)
  } 
  if (!(plot.coverage)) {
    plot.cov <- "not selected"
  }  
  
  # Remove unused objects
  input.before.filter <- NULL
  
  # Filtering call rate ---------------------------------------------------------
  if (!is.null(filter.call.rate)) {
    filter <- input %>% 
      filter(CALL_RATE >= filter.call.rate)
    
    whitelist.filter <- filter %>% distinct(MARKERS, LOCUS, POS)
    
    blacklist.call.rate <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter call rate: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers deleted"))
    if (length(blacklist.call.rate$MARKERS > 0)) {
      write_tsv(x = blacklist.call.rate, path = "blacklist.call.rate.tsv", col_names = TRUE)
    }
    if (plot.call.rate) {
      data.combined <- bind_rows(
        data.before <- input %>% 
          select(POP_ID, INDIVIDUALS, CALL_RATE) %>% 
          mutate(GROUP = rep("before", n())),
        data.after <- filter %>% 
          select(POP_ID, INDIVIDUALS, CALL_RATE) %>% 
          mutate(GROUP = rep("after", n()))
      ) %>% 
        mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
      
      plot.cr <- ggplot(data.combined, aes(x = factor(POP_ID), y = CALL_RATE, na.rm = TRUE))+
        geom_violin(trim = TRUE)+
        geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
        stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
        labs(x = "Sampling sites")+
        labs(y = "Markers call rate")+
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
          legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )+
        facet_grid(~GROUP)
      
    }
    if(!(plot.call.rate)) {
      plot.cr <- "not selected"
    }
    input <- filter
  }
  if (is.null(filter.call.rate) & (plot.call.rate)) {
    data.combined <- input %>% 
      select(POP_ID, INDIVIDUALS, CALL_RATE) %>% 
      mutate(GROUP = rep("before filter", n()))
    
    plot.cr <- ggplot(data.combined, aes(x = factor(POP_ID), y = CALL_RATE, na.rm = TRUE))+
      geom_violin(trim = TRUE)+
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
      labs(x = "Sampling sites")+
      labs(y = "Markers call rate before filter")+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
      )+
      facet_grid(~GROUP)
  }
  
  # Filtering genotyped individuals --------------------------------------------
  # filter.ind.missing.geno = NULL,
  # plot.ind.missing.geno = FALSE,
  
  if (!is.null(filter.ind.missing.geno)) {
    # filter.ind.missing.geno <- 0.90 # test
    
    # Create a new df with genotyped prop.
    ind.geno.prop <- input %>% 
      group_by(INDIVIDUALS) %>% 
      summarise(
        GENOTYPED_PROP = length(GT[GT != "0_0"])/length(GT)
      )
    
    # merge with dataset
    input <- input %>% 
      full_join(ind.geno.prop, by = "INDIVIDUALS")
    
    # filter
    filter <- input %>% 
      filter(GENOTYPED_PROP >= filter.ind.missing.geno)
    
    
    whitelist.filter <- filter %>% distinct(INDIVIDUALS)
    
    blacklist.ind.missing.geno <- input %>% 
      distinct(INDIVIDUALS) %>% 
      filter(!INDIVIDUALS %in% whitelist.filter$INDIVIDUALS)
    
    message(stri_paste("Filter individuals with less than ", filter.ind.missing.geno, " missing genotype prop: ", n_distinct(input$INDIVIDUALS) - n_distinct(filter$INDIVIDUALS), " individuals removed"))
    if (length(blacklist.ind.missing.geno$INDIVIDUALS > 0)) {
      write_tsv(x = blacklist.ind.missing.geno, path = "blacklist.ind.missing.geno.tsv", col_names = TRUE)
    }
    
    if (plot.ind.missing.geno) {
      data.combined <- bind_rows(
        data.before <- input %>% 
          select(POP_ID, INDIVIDUALS, GENOTYPED_PROP) %>% 
          mutate(GROUP = rep("before", n())),
        data.after <- filter %>% 
          select(POP_ID, INDIVIDUALS, GENOTYPED_PROP) %>% 
          mutate(GROUP = rep("after", n()))
      ) %>% 
        mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
      
      ind.missing.geno.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = GENOTYPED_PROP, na.rm = TRUE))+
        geom_violin(trim = TRUE)+
        geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
        stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
        labs(x = "Sampling sites")+
        labs(y = "Individual's genotyped proportion")+
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
          legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )+
        facet_grid(~GROUP)
    }
    
    if(!(plot.ind.missing.geno)) {
      ind.missing.geno.plot <- "not selected"
    }
    input <- filter
  }
  
  if (is.null(filter.ind.missing.geno) & (plot.ind.missing.geno)) {
    # Create a new df with genotyped prop.
    ind.geno.prop <- input %>% 
      group_by(INDIVIDUALS) %>% 
      summarise(
        GENOTYPED_PROP = length(GT[GT != "0_0"])/length(GT)
      )
    
    # merge with dataset
    input <- input %>% 
      full_join(ind.geno.prop, by = "INDIVIDUALS")
    
    data.combined <- input %>% 
      select(POP_ID, INDIVIDUALS, GENOTYPED_PROP) %>% 
      mutate(GROUP = rep("before filter", n()))
    
    ind.missing.geno.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = GENOTYPED_PROP, na.rm = TRUE))+
      geom_violin(trim = TRUE)+
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
      labs(x = "Sampling sites")+
      labs(y = "Individual's genotyped proportion")+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
      )+
      facet_grid(~GROUP)
  }
  
  # Filtering filter.markers.missing.ind  -------------------------------------------
  # filter.markers.missing.ind = NULL,
  # plot.markers.missing.ind = FALSE,
  if (!is.null(filter.markers.missing.ind)) {
    # filter.markers.missing.ind <- 0.95 # test
    
    # Create a new df with marker.missing.ind prop.
    marker.missing.ind.prop <- input %>% 
      group_by(MARKERS) %>% 
      summarise(
        MISSING_IND_PROP = length(GT[GT != "0_0"])/length(GT)
      )
    
    # merge with dataset
    input <- input %>% 
      full_join(marker.missing.ind.prop, by = "MARKERS")
    
    # filter
    filter <- input %>% 
      filter(MISSING_IND_PROP >= filter.markers.missing.ind)
    
    whitelist.filter <- filter %>% distinct(MARKERS)
    
    blacklist.marker.missing.ind.prop <- input %>% 
      distinct(MARKERS) %>% 
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter markers with less than ", filter.markers.missing.ind, " missing genotype ind: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers removed"))
    if (length(blacklist.marker.missing.ind.prop$MARKERS > 0)) {
      write_tsv(x = blacklist.marker.missing.ind.prop, path = "blacklist.marker.missing.ind.tsv", col_names = TRUE)
    }
    
    if (plot.markers.missing.ind) {
      data.combined <- bind_rows(
        data.before <- input %>% 
          select(POP_ID, MARKERS, MISSING_IND_PROP) %>% 
          mutate(GROUP = rep("before", n())),
        data.after <- filter %>% 
          select(POP_ID, MARKERS, MISSING_IND_PROP) %>% 
          mutate(GROUP = rep("after", n()))
      ) %>% 
        mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
      
      markers.missing.ind.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = MISSING_IND_PROP, na.rm = TRUE))+
        geom_violin(trim = TRUE)+
        geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
        stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
        labs(x = "Sampling sites")+
        labs(y = "Marker's genotyped proportion")+
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
          legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )+
        facet_grid(~GROUP)
    }
    
    if(!(plot.markers.missing.ind)) {
      markers.missing.ind.plot <- "not selected"
    }
    input <- filter
  }
  
  if (is.null(filter.markers.missing.ind) & (plot.markers.missing.ind)) {
    data.combined <- input %>% 
      select(POP_ID, MARKERS, MISSING_IND_PROP) %>% 
      mutate(GROUP = rep("before filter", n()))
    
    markers.missing.ind.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = MISSING_IND_PROP, na.rm = TRUE))+
      geom_violin(trim = TRUE)+
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
      labs(x = "Sampling sites")+
      labs(y = "Marker's genotyped proportion")+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
      )+
      facet_grid(~GROUP)
  }
  
  
  # snp.ld  --------------------------------------------------------------------
  
  # get the number of SNP per reads/haplotypes/locus
  number.snp.reads <- input %>% 
    group_by (LOCUS) %>%
    summarise(SNP_N = n_distinct(POS))
  if (plot.number.snp.reads) {
    number.snp.reads.plot <- ggplot(number.snp.reads, aes(factor(SNP_N)))+
      geom_bar()+
      labs(x="Number of SNP per haplotypes (reads)")+
      labs(y="Distribution (number)")+
      theme(axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
            axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
            legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
            legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
            strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold"))
  } 
  
  if (!(plot.number.snp.reads)) {
    number.snp.reads.plot <- "not selected"
  }
  
  # filter snp per reads -------------------------------------------------------
  # 1 snp
  if (filter.snp.ld == "1snp") {
    whitelist.filter <- number.snp.reads %>% 
      filter(SNP_N == 1) %>% 
      distinct(LOCUS) %>% 
      arrange(LOCUS)
    
    blacklist.snp.per.reads <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!LOCUS %in% whitelist.filter$LOCUS)
    write_tsv(x = blacklist.snp.per.reads, path = "blacklist.1snp.per.reads.only.tsv", col_names = TRUE)
    
    input <- input %>% 
      semi_join(whitelist.filter, by = c("LOCUS"))
  }
  
  # 2snp
  if(filter.snp.ld == "2snp") {
    whitelist.filter <- number.snp.reads %>% 
      filter(SNP_N <= 2) %>% 
      distinct(LOCUS) %>% 
      arrange(LOCUS)
    
    blacklist.snp.per.reads <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!LOCUS %in% whitelist.filter$LOCUS)
    write_tsv(x = blacklist.snp.per.reads, path = "blacklist.max2snp.per.reads.only.tsv", col_names = TRUE)
    
    input <- input %>% 
      semi_join(whitelist.filter, by = c("LOCUS"))
  }
  
  # Random selection
  if (filter.snp.ld == "random") {
    snp.locus <- input %>% distinct(LOCUS, POS)
    
    whitelist.filter <- snp.locus %>%
      group_by(LOCUS) %>%
      sample_n(size = 1, replace = FALSE) %>% 
      tidyr::unite(MARKERS, c(LOCUS, POS))
    
    blacklist.snp.per.reads <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    write_tsv(x = blacklist.snp.per.reads, path = "blacklist.1snp.per.reads.random.tsv", col_names = TRUE)
    
    input <- input %>% 
      semi_join(whitelist.filter, by = c("MARKERS"))
  }
  
  # Last SNP on the read (where most error usually occurs)
  if (filter.snp.ld == "last") {
    snp.locus <- input %>% distinct(LOCUS, POS)
    
    whitelist.filter <- snp.locus %>%
      group_by(LOCUS) %>%
      summarise(POS = max(POS)) %>% 
      tidyr::unite(MARKERS, c(LOCUS, POS))
    
    blacklist.snp.per.reads <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    write_tsv(x = blacklist.snp.per.reads, path = "blacklist.1snp.per.reads.last.tsv", col_names = TRUE)
    
    input <- input %>% 
      semi_join(whitelist.filter, by = c("MARKERS"))
  }
  
  # First SNP on the read
  if (filter.snp.ld == "first") {
    snp.locus <- input %>% distinct(LOCUS, POS)
    
    whitelist.filter <- snp.locus %>%
      group_by(LOCUS) %>%
      summarise(POS = min(POS)) %>% 
      tidyr::unite(MARKERS, c(LOCUS, POS))
    
    blacklist.snp.per.reads <- input %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    write_tsv(x = blacklist.snp.per.reads, path = "blacklist.1snp.per.reads.first.tsv", col_names = TRUE)
    
    input <- input %>% 
      semi_join(whitelist.filter, by = c("MARKERS"))
  }
  
  # Filter Minor Allele Frequency* ***********************************************
  # maf.thresholds <- c(0.05, 0.1) # test
  if (!is.null(maf.thresholds)) { # with MAF
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
    
    message("Calculating global and local MAF, this may take some time on large data set")
    
    maf.local <- input %>%
      filter(GT != "0_0") %>%
      group_by(MARKERS, POP_ID, REF, ALT) %>%
      summarise(
        N = as.numeric(n()),
        PQ = as.numeric(length(GT[GT == "REF_ALT"])),
        QQ = as.numeric(length(GT[GT == "ALT_ALT"]))
      ) %>%
      mutate(MAF_LOCAL = ((QQ * 2) + PQ) / (2 * N))
    
    maf.global <- maf.local %>%
      group_by(MARKERS) %>%
      summarise_each_(funs(sum), vars = c("N", "PQ", "QQ")) %>%
      mutate(MAF_GLOBAL = ((QQ * 2) + PQ) / (2 * N)) %>%
      select(MARKERS, MAF_GLOBAL)
    
    maf.data <- maf.global %>%
      left_join(maf.local, by = c("MARKERS")) %>%
      select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL) %>% 
      arrange(MARKERS, POP_ID)
    
    
    maf.local <- NULL
    maf.global <- NULL
    
    write_tsv(x = maf.data, path = "maf.data.tsv", col_names = TRUE, append = FALSE)
    message("MAF table (maf.data.tsv) written in your working directory")
    
    # update the vcf with the maf info
    if (maf.operator == "OR") {
      input.maf <- maf.data %>%
        group_by(MARKERS, POP_ID) %>%
        summarise(
          MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
          MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
        ) %>%
        filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
        group_by(MARKERS) %>%
        tally() %>%
        filter(n >= maf.pop.num.threshold) %>%
        select(MARKERS) %>%
        left_join(input, by = "MARKERS") %>%
        arrange(MARKERS, POP_ID)
    } else { # AND operator between local and global maf
      input.maf <- maf.data %>%
        group_by(MARKERS, POP_ID) %>%
        summarise(
          MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
          MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
        ) %>%
        filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
        group_by(MARKERS) %>%
        tally() %>%
        filter(n >= maf.pop.num.threshold) %>%
        select(MARKERS) %>%
        left_join(input, by = "MARKERS") %>%
        arrange(MARKERS, POP_ID)
    }
    
    message(stri_paste("Filter MAF: ", n_distinct(input$MARKERS) - n_distinct(input.maf$MARKERS), " markers deleted"))
    input <- input.maf
    
    # unused object
    input.maf <- NULL 
  } # End of MAF filters
  
  # Writing to working directory the filtered data frame -----------------------
  
  # Whitelist 
  whitelist.markers <- input %>% select(MARKERS, LOCUS, POS) %>% distinct(MARKERS, .keep_all = TRUE)
  write_tsv(x = whitelist.markers, path = "whitelist.markers.tsv", col_names = TRUE)
  message("Writing the whitelist of markers: whitelist.markers.tsv")
  
  # filtered tidy data 
  input.filtered.df <- input %>% 
    select(MARKERS, LOCUS, POS, REF, ALT, GT, POP_ID, INDIVIDUALS, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG) %>% 
    mutate(
      REF = stri_pad_left(str = REF, pad = "0", width = 3),
      ALT = stri_pad_left(str = ALT, pad = "0", width = 3),
      GENOTYPE = ifelse(GT == "REF_REF", stri_join(REF, REF, sep = ""),
                        ifelse(GT == "ALT_ALT",  stri_join(ALT, ALT, sep = ""), 
                               ifelse(GT == "REF_ALT",  stri_join(REF, ALT, sep = ""), "000000")))
    ) %>% 
    select(-GT) %>% 
    arrange(MARKERS, POP_ID, INDIVIDUALS)  
  
  # write to working directory
  filtered.data.name <- stri_paste(filename, "_tidy.tsv", sep = "")
  message(stri_paste("Writing the tidy and filtered data set: ", filtered.data.name, "\nWorking directory: ", getwd()))
  write_tsv(x = input.filtered.df, path = filtered.data.name, col_names = TRUE)
  
  input <- input %>% 
    select(MARKERS, LOCUS, POS, REF, ALT, GT, POP_ID, INDIVIDUALS, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG) %>% 
    arrange(MARKERS, POP_ID, INDIVIDUALS) 
  
  input <- input %>% 
    select(MARKERS, REF, ALT, GT, POP_ID, INDIVIDUALS) %>% 
    arrange(MARKERS, POP_ID, INDIVIDUALS) %>% 
    mutate(
      REF = stri_pad_left(str = REF, pad = "0", width = 3),
      ALT = stri_pad_left(str = ALT, pad = "0", width = 3),
      GENOTYPE = ifelse(GT == "REF_REF", stri_join(REF, REF, sep = ""),
                        ifelse(GT == "ALT_ALT",  stri_join(ALT, ALT, sep = ""), 
                               ifelse(GT == "REF_ALT",  stri_join(REF, ALT, sep = ""), "000000")))
    ) %>% 
    select(MARKERS, -REF, -ALT, POP_ID, INDIVIDUALS, GT = GENOTYPE)
  
  input$POP_ID <- droplevels(input$POP_ID)
  
  # Imputations **************************************************************
  if (!is.null(imputation.method)) {
    message("Preparing the data for imputations")
    
    input.imp <- stackr::stackr_imputations_module(
      data = input, 
      imputation.method = imputation.method, 
      impute = impute, 
      imputations.group = imputations.group, 
      num.tree = num.tree, 
      iteration.rf = iteration.rf, 
      split.number = split.number, 
      verbose = verbose, 
      parallel.core = parallel.core, 
      filename = NULL
    )
    input.imp$POP_ID <- droplevels(input.imp$POP_ID)
    
    # write to working directory
    filtered.data.name <- stri_paste(filename, "_tidy_imputed", ".tsv", sep = "")
    message(stri_paste("Writing the tidy, filtered and imputed data set: ", filtered.data.name, "\nWorking directory: ", getwd()))
    write_tsv(x = input.imp, path = filtered.data.name, col_names = TRUE)
  } # End imputations
  
  if (is.null(imputation.method)) {
    genind.data.imp <- "not selected"
    input.imp <- "not selected"
  }
  
  
  # PLINK ------------------------------------------------------------------------
  write_plink <- function(x, filename) {
    tped <- x %>% 
      arrange(INDIVIDUALS) %>% 
      mutate(
        COL1 = rep("0", n()),
        COL3 = rep("0", n()),
        COL4 = rep("0", n())
      ) %>% 
      select(COL1, MARKERS, COL3, COL4, INDIVIDUALS, GT) %>% 
      mutate(
        A1 = stri_sub(str = GT, from = 1, to = 3),
        A2 = stri_sub(str = GT, from = 4, to = 6)
      ) %>% 
      select(-GT) %>% 
      tidyr::gather(ALLELES, GENOTYPE, -c(COL1, MARKERS, COL3, COL4, INDIVIDUALS)) %>%
      mutate(
        GENOTYPE = as.character(as.numeric(GENOTYPE)),
        GENOTYPE = stri_pad_left(GENOTYPE, width = 2, pad = "0")
      ) %>%  
      arrange(INDIVIDUALS, ALLELES) %>% 
      tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_") %>%
      group_by(COL1, MARKERS, COL3, COL4) %>% 
      tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GENOTYPE) %>% 
      arrange(MARKERS)
    
    tfam <- x %>%
      distinct(POP_ID, INDIVIDUALS) %>% 
      arrange(INDIVIDUALS) %>% 
      mutate(
        COL3 = rep("0",n()),
        COL4 = rep("0",n()),
        COL5 = rep("0",n()),
        COL6 = rep("-9",n())
      )
    write_delim(x = tped, path = stri_paste(filename, ".tped", sep = ""), col_names = FALSE, delim = " ")
    write_delim(x = tfam, path = stri_paste(filename, ".tfam", sep = ""), col_names = FALSE, delim = " ")
  } # end write_plink
  
  if (plink) {
    message("Generating the PLINK tped and tfam files")
    write_plink (x = input, filename = filename)
    
    if (!is.null(imputation.method)) {
      message("Generating the PLINK tped and tfam files: with imputations")
      write_plink (x = input.imp, filename = stri_paste(filename, "_imputed"))
    } # end plink imputed
    
    
  } else {
    message("PLINK output not selected")
  } # end plink
  
  
  # genind object for adegenet -------------------------------------------------
  prepare_genind <- function(x) {
    genind.prep <- x %>% 
      select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
      #faster than: tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
      mutate(
        A1 = stri_sub(str = GT, from = 1, to = 3),
        A2 = stri_sub(str = GT, from = 4, to = 6)
      ) %>% 
      select(-GT) %>% 
      tidyr::gather(data = ., key = ALLELES, 
                    value = GT, 
                    -c(MARKERS, INDIVIDUALS, POP_ID)
      ) %>% # just miliseconds longer than data.table.melt so keeping this one for simplicity
      filter(GT != "000") # remove missing "000"
    
    # this reintroduce the missing, but with NA
    genind.prep <- data.table::dcast.data.table(
      data = as.data.table(genind.prep), 
      formula = POP_ID + INDIVIDUALS + ALLELES ~ MARKERS, 
      value.var = "GT") %>% 
      as_data_frame() %>% 
      plyr::colwise(.fun = factor, exclude = NA)(.) %>% 
      mutate(INDIVIDUALS = as.character(INDIVIDUALS))
    
    # The next part is longer than it used to be with VCF file only, 
    # but it as the advantage of working and simplifying the use for other file type.
    genind.prep <- suppressWarnings(mutate_each(tbl = genind.prep, funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)))
    
    genind.prep <- tidyr::gather(data = genind.prep, 
                                 key = MARKERS, 
                                 value = GT, 
                                 -c(INDIVIDUALS, POP_ID, ALLELES)
    ) %>% # faster than data.table.melt...
      mutate(GT = stri_replace_na(str = GT, replacement = "000")) %>%
      filter(GT != "000") %>%
      select(-ALLELES) %>%
      group_by(POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
      tally %>% # count alleles, longest part of the block
      ungroup()
    
    genind.prep <- genind.prep %>%
      mutate(MARKERS_ALLELES = stri_paste(MARKERS, GT, sep = ":")) %>%  # faster then: tidyr::unite(MARKERS_ALLELES, MARKERS, GT, sep = ":", remove = TRUE)
      select(-GT, -MARKERS) %>% 
      arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES)
    
    genind.prep <- data.table::dcast.data.table(
      data = as.data.table(genind.prep), 
      formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES, 
      value.var = "n") %>% 
      as_data_frame()
    
    genind.prep <- tidyr::gather(data = genind.prep, key = MARKERS_ALLELES, value = COUNT, -c(INDIVIDUALS, POP_ID)) %>% 
      tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = ":", remove = TRUE) %>% 
      mutate(COUNT = as.numeric(stri_replace_na(str = COUNT, replacement = "0"))) %>% 
      group_by(INDIVIDUALS, MARKERS) %>%
      mutate(MAX_COUNT_MARKERS = max(COUNT, na.rm = TRUE)) %>%
      ungroup() %>% 
      mutate(COUNT = ifelse(MAX_COUNT_MARKERS == 0, "erase", COUNT)) %>%
      select(-MAX_COUNT_MARKERS) %>% 
      mutate(COUNT = replace(COUNT, which(COUNT == "erase"), NA)) %>% 
      arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES)
    
    genind.prep <- genind.prep %>%
      mutate(MARKERS_ALLELES = stri_paste(MARKERS, ALLELES, sep = ".")) %>%  # faster then: tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE)
      select(-MARKERS, -ALLELES) %>% 
      mutate(
        POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
        POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
      )
    
    genind.prep <- data.table::dcast.data.table(
      data = as.data.table(genind.prep), 
      formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES, 
      value.var = "COUNT") %>% 
      as_data_frame() %>%
      arrange(POP_ID, INDIVIDUALS)
  } # End prepare genind
  
  if (genind) {
    message("Preparing data for adegenet genind object...")
    genind.prep <- prepare_genind(x = input)
    
    # genind construction: no imputation -----------------------------------------
    # genind arguments common to all data.type
    message("Building the genind object")
    ind <- genind.prep$INDIVIDUALS
    pop <- genind.prep$POP_ID
    genind.df <- genind.prep %>% ungroup() %>% 
      select(-c(INDIVIDUALS, POP_ID))
    rownames(genind.df) <- ind
    loc.names <- colnames(genind.df)
    strata <- genind.prep %>% ungroup() %>% distinct(INDIVIDUALS, POP_ID)
    
    # genind constructor
    prevcall <- match.call()
    no.imputation <- genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = NULL)
    
    if (!is.null(imputation.method)) {
      message("Preparing imputed data for adegenet genind object...")
      genind.prep.imp <- prepare_genind(x = input.imp)
      
      # genind construction: with imputations ----------------------------------
      message("Building the imputed genind object")
      ind <- genind.prep.imp$INDIVIDUALS
      pop <- genind.prep.imp$POP_ID
      genind.df <- genind.prep.imp %>%
        ungroup() %>% 
        select(-c(INDIVIDUALS, POP_ID))
      rownames(genind.df) <- ind
      loc.names <- colnames(genind.df)
      
      strata <- genind.prep.imp %>% 
        ungroup() %>% 
        distinct(INDIVIDUALS, POP_ID)
      
      # genind constructor
      prevcall <- match.call()
      imputed  <- adegenet::genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = NULL)

      ind <- NULL
      pop <- NULL
      genind.df <- NULL
      loc.names <- NULL
      strata <- NULL
      prevcall <- NULL
    } else {
      
    }# end genind imp
    imputed <- "not selected"
  } else {
    no.imputation <- "not selected"
    imputed <- "not selected"
  } # end genind
  
  # genepop --------------------------------------------------------------------
  
  if (genepop) {
    message("Generating the genepop file")
    stackr::write_genepop(data = input, genepop.header = "stackr::dart2df_genind_plink: no imputations", markers.line.format = "line", filename = stri_paste(filename, ".gen", sep = ""))

    if (!is.null(imputation.method)) {
      message("Generating the genepop file: with imputations")
      stackr::write_genepop(data = input.imp, genepop.header = "stackr::dart2df_genind_plink: imputed data", markers.line.format = "line", filename = stri_paste(filename, "_imputed.gen", sep = ""))
    }
  } else {
    message("genepop output not selected")
  }
  
  # STRUCTURE --------------------------------------------------------------------
  
  if (structure) {
    message("Generating the structure file")
    stackr::write_structure(data = input, markers.line = TRUE, filename = stri_paste(filename, ".str", sep = ""))

    if (!is.null(imputation.method)) {
      message("Generating the structure file: with imputations")
      stackr::write_structure(data = input.imp, markers.line = TRUE, filename = stri_paste(filename, "_imputed.str", sep = ""))
    }
  } else {
    message("structure output not selected")
  }
  
  input <- tidyr::separate(data = input, col = MARKERS, into = c("LOCUS", "POS"), sep = "_", remove = FALSE, extra = "drop")
  
  # Results ----------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stri_paste("The number of markers removed by the filters:\nSNP: ", snp.before.filters - n_distinct(input$MARKERS), "\nLOCUS: ", locus.before.filters - n_distinct(input$LOCUS)))
  message("The number of markers before -> after the filters")
  message(stri_paste("SNP: ", snp.before.filters, " -> ", n_distinct(input$MARKERS)))
  message(stri_paste("LOCUS: ", locus.before.filters, " -> ", n_distinct(input$LOCUS)))
  cat("#######################################################################\n")
  
  res <- list()
  res$jitterplot.ind.het <- jitterplot.ind.het
  res$boxplot.ind.het <- boxplot.ind.het
  res$plot.reproducibility <- plot.repro
  res$plot.coverage <- plot.cov
  res$plot.call.rate <- plot.cr
  res$plot.ind.missing.geno <- ind.missing.geno.plot
  res$plot.markers.missing.ind <- markers.missing.ind.plot
  res$plot.number.snp.reads <- number.snp.reads.plot
  res$individual.heterozygosity <- het.ind
  res$maf.data <- maf.data
  res$whitelist.markers <- whitelist.markers
  res$data.filtered <- input.filtered.df
  res$data.filtered.imputed <- input.imp
  res$genind.data <- no.imputation
  res$genind.data.imputed <- imputed
  return(res)
}
