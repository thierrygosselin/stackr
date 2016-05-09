# Import, filter and transform a dart output file to different format

#' @name dart2df_genind_plink
#' @title swiss army knife tool to prepare DArT output file for population genetics analysis.
#' @description Import, filter and transform a DArT output file to different format: 
#' data frame of genotypes, genind object and/or PLINK \code{tped/tfam} format.
#' @param data DArT output file in wide format or binary format tipically 
#' used by CSIRO genomic projects.

#' @param strata A tab delimited file with columns header:
#' \code{INDIVIDUALS} and \code{POP_ID}. 
#' Use the column \code{POP_ID} for any hierarchical grouping. If a third column
#' named \code{NEW_ID} is used, this column will be used to replace the
#' \code{INDIVIDUALS} column in the main data file.

#' @param pop.levels (optional string) A character string with your populations ordered.

#' @param pop.select (optional string) Keep specific populations. 
#' Default = \code{NULL} for no selection and keep all population.
#' e.g. \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population samples.
#' 
#' @param filter.monomorphic (optional, logical) Should the monomorphic markers present in 
#' the dataset be filtered out ? Default: \code{filter.monomorphic = TRUE}.

#' @param common.markers (optional, logical) Default: \code{common.markers = TRUE}, 
#' will only keep markers in common (genotyped) between all the populations.

#' @param filter.reproducibility (optional, numerical) Filter the \code{RepAvg} 
#' column in the data set. Default: \code{filter.reproducibility = NULL}.
#' e.g to keep markers with reproducibility >= 99%,
#' use : \code{filter.reproducibility = 0.99}.

#' @param plot.reproducibility (optional, logical) Plot the distribution 
#' (violin plot) of reproducibility. Default: \code{plot.reproducibility = NULL}.

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
#' (violin plot). Default: \code{plot.coverage = NULL}.

#' @param filter.call.rate (optional, numerical) Filter the \code{CallRate} 
#' column in the data set. Default: \code{filter.call.rate = NULL}. e.g to keep 
#' markers genotyped in more than 95% of the individuals use :
#' \code{filter.call.rate = 0.95}
#' @param plot.call.rate (optional, logical) Plot the distribution 
#' (violin plot) of call rate. Default: \code{plot.call.rate = NULL}.



#' @param filter.snp.ld (optional, character) With anonymous markers from
#' reduce representation library like RADseq/GBS/DArT de novo discovery, 
#' you can explore the impact of the number of snp on the read (haplotype/locus).
#' This argument is used to minimize linkage disequilibrium (LD) on short reads,
#' by choosing among these options: 
#' \code{"1snp"} to keep only reads with 1 SNP/reads,
#' \code{"2snp"} to keep only reads with at most 2 SNP/reads,
#'\code{"random"} to select 1 SNP/reads randomly, 
#'  \code{"first"} to select the first SNP on all reads, 
#' \code{"last"} to select the last SNP on all reads 
#' (last SNP are usually associated to higher error rate...). 
#' Default: \code{snp.ld = NULL}.
#' Note: for long linkage detection use PLINK linkage disequilibrium based SNP 
#' pruning.

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
#' @param genind (optional, logical). To have the filtered data set output as 
#' genind object to use inside adegenet, use \code{genind = TRUE}.

#' @param filename (optional) The name of the file written to the directory.
#' No file extension at the end.

#' @param imputation.method Should a map-independent imputations of markers be
#' computed. Available choices are: (1) \code{FALSE} for no imputation.
#' (2) \code{"max"} to use the most frequent category for imputations.
#' (3) \code{"rf"} using Random Forest algorithm. 
#' Default: \code{imputation.method = FALSE}.
#' @param impute (character) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
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
      "AVG_COUNT_SNP", "REP_AVG", "NEW_ID", "SNP_N", "ALLELE_NAME", "ALLELE_NUMBER"
    )
  )
}

dart2df_genind_plink <- function(data,
                                 strata,
                                 pop.levels,
                                 pop.select,
                                 filter.monomorphic,
                                 common.markers,
                                 filter.reproducibility,
                                 plot.reproducibility,
                                 filter.coverage.high,
                                 filter.coverage.low,
                                 plot.coverage,
                                 filter.call.rate,
                                 plot.call.rate,
                                 filter.snp.ld,
                                 plot.number.snp.reads,
                                 maf.thresholds,
                                 maf.pop.num.threshold,
                                 maf.operator,
                                 plink,
                                 genind,
                                 filename,
                                 imputation.method = FALSE,
                                 impute = "genotypes",
                                 imputations.group = "populations",
                                 num.tree = 100,
                                 iteration.rf = 10,
                                 split.number = 100,
                                 verbose = FALSE,
                                 parallel.core = NULL,
                                 ...) {
  
  cat("#######################################################################\n")
  cat("#################### stackr: dart2df_genind_plink #####################\n")
  cat("#######################################################################\n")
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (missing(strata)) stop("strata file missing")
  if (missing(pop.levels)) pop.levels <- NULL
  if (missing(pop.select)) pop.select <- NULL
  if (missing(filter.monomorphic)) filter.monomorphic <- TRUE
  if (missing(common.markers)) common.markers <- TRUE
  if (missing(filter.reproducibility)) filter.reproducibility <- NULL
  if (missing(plot.reproducibility)) plot.reproducibility <- NULL
  if (missing(filter.coverage.high)) filter.coverage.high <- NULL
  if (missing(filter.coverage.low)) filter.coverage.low <- NULL
  if (missing(plot.coverage)) plot.coverage <- NULL
  if (missing(filter.call.rate)) filter.call.rate <- NULL
  if (missing(plot.call.rate)) plot.call.rate <- NULL
  if (missing(filter.snp.ld)) filter.snp.ld <- NULL
  if (missing(plot.number.snp.reads)) plot.number.snp.reads <- NULL
  if (missing(maf.thresholds)) maf.thresholds <- NULL
  if (missing(maf.pop.num.threshold)) maf.pop.num.threshold <- 1
  if (missing(maf.operator)) maf.operator <- "OR"
  if (missing(plink)) plink <- NULL
  if (missing(genind)) genind <- NULL
  if (missing(filename)) filename <- NULL
  if (missing(imputation.method)) imputation.method <- FALSE
  if (missing(imputations.group)) imputations.group <- "populations"
  if (imputation.method != FALSE & missing(impute)) stop("impute argument is necessary")
  if (imputation.method == FALSE & missing(impute)) impute <- NULL
  if (missing(num.tree)) num.tree <- 100
  if (missing(iteration.rf)) iteration.rf <- 10
  if (missing(split.number)) split.number <- 100
  if (missing(verbose)) verbose <- FALSE
  if (missing(parallel.core) | is.null(parallel.core)) parallel.core <- detectCores()-1
  if (filter.monomorphic == FALSE) filter.monomorphic <- NULL
  if (common.markers == FALSE) common.markers <- NULL
  if (plink == FALSE) plink <- NULL
  if (genind == FALSE) genind <- NULL
  if (plot.reproducibility == FALSE) plot.reproducibility <- NULL
  if (plot.coverage == FALSE) plot.coverage <- NULL
  if (plot.call.rate == FALSE) plot.call.rate <- NULL
  if (plot.number.snp.reads == FALSE) plot.number.snp.reads <- NULL
  
  
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
  
  # Import data ---------------------------------------------------------------
  input <- suppressWarnings(
    data.table::fread(
      input = data,
      sep = "\t",
      stringsAsFactors = FALSE, 
      header = TRUE,
      na.strings = "-",
      strip.white = TRUE,
      drop = c("AlleleID", "AlleleSequence", "AlleleSequenceREF", "AlleleSequenceSNP", "OneRatioRef", "OneRatioSnp", "FreqHomRef", "FreqHomSnp", "FreqHets", "PICRef", "PICSnp", "AvgPIC"),
      showProgress = TRUE,
      verbose = FALSE
    ) %>% 
      as_data_frame() %>%
      rename(LOCUS = CloneID, POS = SnpPosition, CALL_RATE = CallRate, AVG_COUNT_REF = AvgCountRef, AVG_COUNT_SNP = AvgCountSnp, REP_AVG = RepAvg) %>% 
      arrange(LOCUS, POS)
  )  
  
  # Determine the type of DArT file
  binary <- anyDuplicated(input$LOCUS)
  
  # Screen for duplicate -----------------------------------------------------
  remove.list <- c("LOCUS", "SNP", "POS", "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
  individuals.df <- data_frame(INDIVIDUALS = purrr::discard(.x = colnames(input), .p = colnames(input) %in% remove.list))
  duplicate.individuals <- length(individuals.df$INDIVIDUALS) - n_distinct(individuals.df$INDIVIDUALS)
  if (duplicate.individuals == 0) {
    message("Duplicate individual in the data: no")
  } else {
    stop(stri_paste("Duplicated individuals id found in the data set.\nNumber of duplicates= ", duplicate.individuals))
  }
  # removing unused object
  remove.list <- NULL
  individuals.df <- NULL
  duplicate.individuals <- NULL
  
  # Tidying data ---------------------------------------------------------------
  if (binary == 0) {
    message("Tidying the dataset")
    input <- input %>% 
      # working on the column to keep the interesting info
      tidyr::separate(col = LOCUS, into = c("LOCUS", "NOT_USEFUL"), sep = "\\|", extra = "drop") %>%  
      select(-NOT_USEFUL) %>% # discard the 'NOT_USEFUL' column
      tidyr::separate(col = SNP, into = c("NOT_USEFUL", "KEEPER"), sep = ":", extra = "drop") %>% 
      select(-NOT_USEFUL) %>% # discard the 'NOT_USEFUL' column 
      tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>% 
      # merge LOCUS and POS columns into MARKERS
      # gather individuals into 1 columns to work in long format
      tidyr::unite(MARKERS, LOCUS, POS, sep = "_", remove = FALSE ) %>%
      tidyr::gather(INDIVIDUALS, GT, -c(MARKERS, LOCUS, POS, REF, ALT, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG)) %>% 
      # replace nucleotide with numbers
      mutate(
        REF = stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("01", "02", "03", "04"), vectorize_all = FALSE), # replace nucleotide with numbers
        ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("01", "02", "03", "04"), vectorize_all = FALSE), # replace nucleotide with numbers
        GT = ifelse(GT == 0, "REF_REF",
                    ifelse(GT == 1, "ALT_ALT", 
                           ifelse(GT == 2, "REF_ALT", "0_0"))),
        GT = stri_replace_na(str = GT, replacement = "0_0")
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
      distinct(MARKERS, LOCUS, POS, REF, ALT, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG) %>% 
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
  
  
  # Strata file ------------------------------------------------------------------
  strata.df <- read_tsv(file = strata, col_names = TRUE)
  
  input <- input %>%
    mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
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
  
  # pop.select ------------------------------------------------------------------
  if (!is.null(pop.select)) {
    message(stri_join(length(pop.select), "population(s) selected", sep = " "))
    input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
  }
  
  # Filter monomorphic markers  ---------------------------------------------------
  if (filter.monomorphic == TRUE) {
    # screen for monomorphic
    blacklist.monomorphic  <- input %>%
      select(LOCUS, GT) %>%
      filter(GT != "0_0") %>%
      group_by(LOCUS, GT) %>%
      distinct %>% 
      group_by(LOCUS) %>%
      tally %>% 
      filter(n == 1) %>% 
      select(LOCUS) %>% 
      arrange(LOCUS)
    
    # Remove the markers from the dataset
    message(stri_paste("Filter monomorphic: ", n_distinct(blacklist.monomorphic$LOCUS), " markers deleted"))
    if (length(blacklist.monomorphic$LOCUS > 0)) {
      write_tsv(x = blacklist.monomorphic, path = "blacklist.monomorphic.tsv", col_names = TRUE)
      input <- anti_join(input, blacklist.monomorphic, by = "LOCUS")
    }
  }
  # Filter common markers between all populations  -------------------------------
  if (common.markers == TRUE) {
    # get the pop number
    pop.number <- n_distinct(input$POP_ID)
    
    # filter
    filter <- input %>% 
      ungroup() %>% 
      filter(GT != "0_0") %>%
      select(MARKERS, POP_ID) %>% 
      distinct(MARKERS, POP_ID) %>% 
      group_by(MARKERS) %>%
      tally() %>% 
      filter(n == pop.number) %>% 
      arrange(MARKERS) %>%
      select(MARKERS) %>%
      distinct(MARKERS)
    
    whitelist.filter <- filter %>% select(MARKERS) %>% distinct(MARKERS)
    
    blacklist.common.markers <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
      distinct(MARKERS, LOCUS, POS) %>%
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
      select(MARKERS, LOCUS, POS) %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter reproducibility: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers deleted"))
    if (length(blacklist.markers.reproducibility$MARKERS > 0)) {
      write_tsv(x = blacklist.markers.reproducibility, path = "blacklist.markers.reproducibility.tsv", col_names = TRUE)
    }
    if (!is.null(plot.reproducibility)) {
      data.combined <- bind_rows(
        data.before <- input %>% 
          select(POP_ID, INDIVIDUALS, REP_AVG) %>% 
          mutate(GROUP = rep("before", n())),
        data.after <- filter %>% 
          select(POP_ID, INDIVIDUALS, REP_AVG) %>% 
          mutate(GROUP = rep("after", n()))
      ) %>% 
        mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
      
      reproducibility.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = REP_AVG, na.rm = TRUE))+
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
    if (is.null(plot.reproducibility)) {
      reproducibility.plot <- "not selected"
    }
    input <- filter
  }
  
  if (is.null(filter.reproducibility) & !is.null(plot.reproducibility)) {
    data.combined <- input %>% 
      select(POP_ID, INDIVIDUALS, REP_AVG) %>% 
      mutate(GROUP = rep("before filter", n()))
    
    reproducibility.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = REP_AVG, na.rm = TRUE))+
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
    
    whitelist.filter <- filter %>% select(MARKERS, LOCUS, POS) %>% distinct(MARKERS, LOCUS, POS)
    
    blacklist.markers.coverage.high <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
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
    
    whitelist.filter <- filter %>% select(MARKERS, LOCUS, POS) %>% distinct(MARKERS, LOCUS, POS)
    
    blacklist.markers.coverage.low <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter coverage low: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers deleted"))
    if (length(blacklist.markers.coverage.low$MARKERS > 0)) {
      write_tsv(x = blacklist.markers.coverage.low, path = "blacklist.markers.coverage.low.tsv", col_names = TRUE)
      input <- filter
    }
  }
  
  if (!is.null(filter.coverage.high) | !is.null(filter.coverage.low) & !is.null(plot.coverage)) {
    data.combined <- bind_rows(
      data.before <- input.before.filter %>% 
        select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>% 
        mutate(GROUP = rep("before", n())),
      data.after <- input %>% 
        select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>% 
        mutate(GROUP = rep("after", n()))
    ) %>% 
      mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
    
    coverage.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = AVG_COUNT_SNP, na.rm = TRUE))+
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
  if (is.null(filter.coverage.high) & is.null(filter.coverage.low) & !is.null(plot.coverage)) {
    data.combined <- input.before.filter %>% 
      select(POP_ID, INDIVIDUALS, AVG_COUNT_SNP) %>% 
      mutate(GROUP = rep("before filter", n()))
    
    coverage.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = AVG_COUNT_SNP, na.rm = TRUE))+
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
  if (is.null(plot.coverage)) {
    coverage.plot <- "not selected"
  }  
  
  # Remove unused objects
  input.before.filter <- NULL
  
  # Filtering call rate ---------------------------------------------------------
  if (!is.null(filter.call.rate)) {
    filter <- input %>% 
      filter(CALL_RATE >= filter.call.rate)
    
    whitelist.filter <- filter %>% select(MARKERS, LOCUS, POS) %>% distinct(MARKERS, LOCUS, POS)
    
    blacklist.call.rate <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    
    message(stri_paste("Filter call rate: ", n_distinct(input$MARKERS) - n_distinct(filter$MARKERS), " markers deleted"))
    if (length(blacklist.call.rate$MARKERS > 0)) {
      write_tsv(x = blacklist.call.rate, path = "blacklist.call.rate.tsv", col_names = TRUE)
    }
    if (!is.null(plot.call.rate)) {
      data.combined <- bind_rows(
        data.before <- input %>% 
          select(POP_ID, INDIVIDUALS, CALL_RATE) %>% 
          mutate(GROUP = rep("before", n())),
        data.after <- filter %>% 
          select(POP_ID, INDIVIDUALS, CALL_RATE) %>% 
          mutate(GROUP = rep("after", n()))
      ) %>% 
        mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
      
      call.rate.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = CALL_RATE, na.rm = TRUE))+
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
    if(is.null(plot.call.rate)) {
      call.rate.plot <- "not selected"
    }
    input <- filter
  }
  if (is.null(filter.call.rate) & !is.null(plot.call.rate)) {
    data.combined <- input %>% 
      select(POP_ID, INDIVIDUALS, CALL_RATE) %>% 
      mutate(GROUP = rep("before filter", n()))
    
    call.rate.plot <- ggplot(data.combined, aes(x = factor(POP_ID), y = CALL_RATE, na.rm = TRUE))+
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
  
  # snp.ld  --------------------------------------------------------------------
  
  # get the number of SNP per reads/haplotypes/locus
  number.snp.reads <- input %>% 
    group_by (LOCUS) %>%
    summarise(SNP_N = n_distinct(POS))
  if (!is.null(plot.number.snp.reads)) {
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
  
  if (is.null(plot.number.snp.reads)) {
    number.snp.reads.plot <- "not selected"
  }
  
  # filter snp per reads -------------------------------------------------------
  # 1 snp
  if (filter.snp.ld == "1snp") {
    whitelist.filter <- number.snp.reads %>% 
      filter(SNP_N == 1) %>% 
      select(LOCUS) %>% 
      distinct(LOCUS) %>% 
      arrange(LOCUS)
    
    blacklist.snp.per.reads <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
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
      select(LOCUS) %>% 
      distinct(LOCUS) %>% 
      arrange(LOCUS)
    
    blacklist.snp.per.reads <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!LOCUS %in% whitelist.filter$LOCUS)
    write_tsv(x = blacklist.snp.per.reads, path = "blacklist.max2snp.per.reads.only.tsv", col_names = TRUE)
    
    input <- input %>% 
      semi_join(whitelist.filter, by = c("LOCUS"))
  }
  
  # Random selection
  if (filter.snp.ld == "random") {
    snp.locus <- input %>% select(LOCUS, POS) %>% distinct(LOCUS, POS)
    
    whitelist.filter <- snp.locus %>%
      group_by(LOCUS) %>%
      sample_n(size = 1, replace = FALSE) %>% 
      tidyr::unite(MARKERS, c(LOCUS, POS))
    
    blacklist.snp.per.reads <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    write_tsv(x = blacklist.snp.per.reads, path = "blacklist.1snp.per.reads.random.tsv", col_names = TRUE)
    
    input <- input %>% 
      semi_join(whitelist.filter, by = c("MARKERS"))
  }
  
  # Last SNP on the read (where most error usually occurs)
  if (filter.snp.ld == "last") {
    snp.locus <- input %>% select(LOCUS, POS) %>% distinct(LOCUS, POS)
    
    whitelist.filter <- snp.locus %>%
      group_by(LOCUS) %>%
      summarise(POS = max(POS)) %>% 
      tidyr::unite(MARKERS, c(LOCUS, POS))
    
    blacklist.snp.per.reads <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
      distinct(MARKERS, LOCUS, POS) %>%
      filter(!MARKERS %in% whitelist.filter$MARKERS)
    write_tsv(x = blacklist.snp.per.reads, path = "blacklist.1snp.per.reads.last.tsv", col_names = TRUE)
    
    input <- input %>% 
      semi_join(whitelist.filter, by = c("MARKERS"))
  }
  
  # First SNP on the read
  if (filter.snp.ld == "first") {
    snp.locus <- input %>% select(LOCUS, POS) %>% distinct(LOCUS, POS)
    
    whitelist.filter <- snp.locus %>%
      group_by(LOCUS) %>%
      summarise(POS = min(POS)) %>% 
      tidyr::unite(MARKERS, c(LOCUS, POS))
    
    blacklist.snp.per.reads <- input %>% 
      select(MARKERS, LOCUS, POS) %>% 
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
  whitelist.markers <- input %>% select(MARKERS, LOCUS, POS) %>% distinct(MARKERS)
  write_tsv(x = whitelist.markers, path = "whitelist.markers.tsv", col_names = TRUE)
  message("Writing the whitelist of markers: whitelist.markers.tsv")
  
  # filtered tidy data 
  input.filtered.df <- input %>% 
    select(MARKERS, LOCUS, POS, REF, ALT, GT, POP_ID, INDIVIDUALS, CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG) %>% 
    mutate(
      REF = stri_pad_left(str = REF, pad = "0", width = 3),
      ALT = stri_pad_left(str = ALT, pad = "0", width = 3),
      GENOTYPE = ifelse(GT == "REF_REF", stri_join(REF, REF, sep = "/"),
                  ifelse(GT == "ALT_ALT",  stri_join(ALT, ALT, sep = "/"), 
                         ifelse(GT == "REF_ALT",  stri_join(REF, ALT, sep = "/"), "000/000")))
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
  
  # Imputations **************************************************************
  if (imputation.method != "FALSE") {
    message("Preparing the data for imputations")
    
    if (impute == "genotype") {
      input.prep <- input %>%
        select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
        mutate(
          GT = stri_replace_all_fixed(GT, pattern = "0_0", replacement = "NA", vectorize_all = FALSE),
          GT = replace(GT, which(GT == "NA"), NA)
        ) %>%
        group_by(INDIVIDUALS, POP_ID) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>%
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
    }
    
    if (impute == "allele") {
      input.prep <- input %>%
        select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
        mutate(
          GT = stri_replace_all_fixed(GT, pattern = "0_0", replacement = "NA", vectorize_all = FALSE),
          GT = replace(GT, which(GT == "NA"), NA)
        ) %>% 
        tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
        tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
        mutate(GT = replace(GT, which(GT == "NA"), NA)) %>%
        group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
    }
    
    # keep stratification
    strata.df.impute <- input.prep %>% 
      select(INDIVIDUALS, POP_ID) %>% 
      distinct(INDIVIDUALS, POP_ID)
    
    # Imputation with Random Forest
    if (imputation.method == "rf") {
      # Parallel computations options
      options(rf.cores = parallel.core, mc.cores = parallel.core)
      
      # imputations using Random Forest with the package randomForestSRC
      impute_genotype_rf <- function(x) {
        randomForestSRC::impute.rfsrc(data = x,
                                      ntree = num.tree,
                                      nodesize = 1,
                                      nsplit = split.number,
                                      nimpute = iteration.rf,
                                      do.trace = verbose)
      } # End of imputation function
      
      # Random Forest by pop
      if (imputations.group == "populations") {
        message("Imputations computed by populations, take a break...")
        df.split.pop <- split(x = input.prep, f = input.prep$POP_ID) # slip data frame by population
        pop.list <- names(df.split.pop) # list the pop
        imputed.dataset <-list() # create empty list
        
        # Function to go through the populations
        impute_rf_pop <- function(pop.list, ...){
          sep.pop <- df.split.pop[[pop.list]]
          sep.pop <- suppressWarnings(
            plyr::colwise(factor, exclude = NA)(sep.pop)
          )
          # message of progress for imputations by population
          message(paste("Completed imputations for pop ", pop.list, sep = ""))
          # imputed.dataset[[i]] <- impute_markers_rf(sep.pop) # test with foreach
          imputed.dataset <- impute_genotype_rf(sep.pop)
          imputed.dataset <- suppressWarnings(
            plyr::colwise(as.character, exclude = NA)(imputed.dataset)
          )
          return(imputed.dataset)
        } # End impute_rf_pop
        
        input.imp <- list()
        input.imp <- parallel::mclapply(
          X = pop.list, 
          FUN = impute_rf_pop, 
          mc.preschedule = FALSE, 
          mc.silent = FALSE, 
          mc.cores = parallel.core
        )
        
        # Compiling the results
        message("Compiling imputations results")
        input.imp <- suppressWarnings(bind_rows(input.imp))
        
        # Second round of imputations (globally) to remove introduced NA 
        # In case that some pop don't have the markers
        input.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(input.imp)) # Make the columns factor
        input.imp <- impute_genotype_rf(input.imp) # impute globally
        input.imp <- plyr::colwise(as.character, exclude = NA)(input.imp)
        
        # reintroduce the stratification
        input.imp <- suppressWarnings(
          left_join(strata.df.impute, 
                    input.imp %>% 
                      select(-POP_ID)
                    , by = "INDIVIDUALS") %>% 
            arrange(POP_ID, INDIVIDUALS) %>% 
            ungroup()
        )
        
        # dump unused objects
        df.split.pop <- NULL
        pop.list <- NULL
        sep.pop <- NULL
        imputed.dataset <- NULL
        input.prep <- NULL
        
      } # End imputation RF populations
      # Random Forest global
      if (imputations.group == "global") { # Globally (not by pop_id)
        message("Imputations computed globally, take a break...")
        input.prep <- plyr::colwise(factor, exclude = NA)(input.prep)
        input.imp <- impute_genotype_rf(input.prep)
        input.imp <- plyr::colwise(as.character, exclude = NA)(input.imp)
        
        # reintroduce the stratification
        input.imp <- suppressWarnings(
          left_join(strata.df.impute, 
                    input.imp %>% 
                      select(-POP_ID)
                    , by = "INDIVIDUALS") %>% 
            arrange(POP_ID, INDIVIDUALS) %>% 
            ungroup()
        )
        input.prep <- NULL # remove unused object
      } # End imputation RF global
      
      
      # data prep
      if (impute == "genotype") {
        input.imp <- suppressWarnings(
          input.imp %>%
            tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID))
        )
      }
      if (impute == "allele") {
        input.imp <- suppressWarnings(
          input.imp %>%
            tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES))
        )
      }
    } # End imputation RF
    
    # Imputation using the most common genotype
    if (imputation.method == "max") { # End imputation max
      if (imputations.group == "populations") {
        message("Imputations computed by populations")
        
        if (impute == "genotype"){
          input.imp <- suppressWarnings(
            input.prep %>%
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
              ungroup()
          )
        }
        
        if (impute == "allele"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
              group_by(MARKERS, POP_ID) %>%
              mutate(
                GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                GT = replace(GT, which(GT == "NA"), NA)
              ) %>%
              # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
              # will take the global observed values by markers for those cases.
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              ungroup()
          )
        }
        
        input.prep <- NULL # remove unused object
        
      } # End imputation max populations 
      if (imputations.group == "global") {
        # Globally (not by pop_id)
        message("Imputations computed globally")
        if (impute == "genotype"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              ungroup()
          )
        }
        
        if (impute == "allele"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              ungroup()
          )
        }
        
        input.prep <- NULL # remove unused object
      } # End imputation max global 
    } # End imputations max
    
    # prepare the imputed dataset for plink and adegenet
    message("Preparing imputed data set...")
    
    input.keep <- input %>% select(MARKERS, REF, ALT) %>%  distinct(MARKERS, REF, ALT)
    
    if (impute == "genotype") {
      input.imp <- input.imp %>% 
        inner_join(input.keep, by = c("MARKERS"))
    }
    if (impute == "allele") {
      input.imp <- input.imp %>% 
        tidyr::spread(data = ., key = ALLELES, value = GT) %>%
        tidyr::unite(GT, A1, A2, sep = "_", remove = TRUE) %>% 
        inner_join(input.keep, by = c("MARKERS"))
    }
    
    # filtered tidy data imputed -----------------------------------------------
    input.imp.df <- input.imp %>% 
      select(MARKERS, REF, ALT, GT, POP_ID, INDIVIDUALS) %>% 
      arrange(MARKERS, POP_ID, INDIVIDUALS) %>% 
      mutate(
        REF = stri_pad_left(str = REF, pad = "0", width = 3),
        ALT = stri_pad_left(str = ALT, pad = "0", width = 3),
        GENOTYPE = ifelse(GT == "REF_REF", stri_join(REF, REF, sep = "/"),
                    ifelse(GT == "ALT_ALT",  stri_join(ALT, ALT, sep = "/"), 
                           ifelse(GT == "REF_ALT",  stri_join(REF, ALT, sep = "/"), "000/000")))
      ) %>% 
      select(MARKERS, REF, ALT, POP_ID, INDIVIDUALS, GENOTYPE)
    
    # write to working directory
    filtered.data.name <- stri_paste(filename, "_tidy_imputed", ".tsv", sep = "")
    message(stri_paste("Writing the tidy, filtered and imputed data set: ", filtered.data.name, "\nWorking directory: ", getwd()))
    write_tsv(x = input.imp.df, path = filtered.data.name, col_names = TRUE)

    
  } # End imputations
  
  if (imputation.method == FALSE) {
    genind.data.imp <- "not selected"
    input.imp.df <- "not selected"
  }
  
  
  # PLINK ------------------------------------------------------------------------
  if (plink == TRUE) {
    message("Generating the PLINK tped and tfam files")
    tped <- input %>% 
      arrange(INDIVIDUALS) %>% 
      mutate(
        COL1 = rep("0", n()),
        COL3 = rep("0", n()),
        COL4 = rep("0", n())
      ) %>% 
      select(COL1, MARKERS, COL3, COL4, REF, ALT, INDIVIDUALS, GT) %>% 
      mutate(
        GENOTYPE = ifelse(GT == "REF_REF", stri_join(REF, REF, sep = "_"),
                          ifelse(GT == "ALT_ALT",  stri_join(ALT, ALT, sep = "_"), 
                                 ifelse(GT == "REF_ALT",  stri_join(REF, ALT, sep = "_"), "0_0")))
      ) %>%
      select(COL1, MARKERS, COL3, COL4, INDIVIDUALS, GENOTYPE) %>%
      tidyr::separate(col = GENOTYPE, into = c("A1", "A2"), sep = "_") %>% 
      tidyr::gather(ALLELES, GENOTYPE, -c(COL1, MARKERS, COL3, COL4, INDIVIDUALS)) %>% 
      arrange(INDIVIDUALS, ALLELES) %>% 
      tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_") %>%
      group_by(COL1, MARKERS, COL3, COL4) %>% 
      tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GENOTYPE) %>% 
      arrange(MARKERS)
    
    tfam <- input %>%
      select(POP_ID, INDIVIDUALS) %>% 
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
    
    if (imputation.method != "FALSE") {
      message("Generating the imputed PLINK tped and tfam files")
      tped <- input.imp %>% 
        arrange(INDIVIDUALS) %>% 
        mutate(
          COL1 = rep("0", n()),
          COL3 = rep("0", n()),
          COL4 = rep("0", n())
        ) %>% 
        select(COL1, MARKERS, COL3, COL4, REF, ALT, INDIVIDUALS, GT) %>% 
        mutate(
          GENOTYPE = ifelse(GT == "REF_REF", stri_join(REF, REF, sep = "_"),
                            ifelse(GT == "ALT_ALT",  stri_join(ALT, ALT, sep = "_"), 
                                   ifelse(GT == "REF_ALT",  stri_join(REF, ALT, sep = "_"), "0_0")))
        ) %>%
        select(COL1, MARKERS, COL3, COL4, INDIVIDUALS, GENOTYPE) %>%
        tidyr::separate(col = GENOTYPE, into = c("A1", "A2"), sep = "_") %>% 
        tidyr::gather(ALLELES, GENOTYPE, -c(COL1, MARKERS, COL3, COL4, INDIVIDUALS)) %>% 
        arrange(INDIVIDUALS, ALLELES) %>% 
        tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_") %>%
        group_by(COL1, MARKERS, COL3, COL4) %>% 
        tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GENOTYPE) %>% 
        arrange(MARKERS)
      
      tfam <- input.imp %>%
        select(POP_ID, INDIVIDUALS) %>% 
        distinct(POP_ID, INDIVIDUALS) %>% 
        arrange(INDIVIDUALS) %>% 
        mutate(
          COL3 = rep("0",n()),
          COL4 = rep("0",n()),
          COL5 = rep("0",n()),
          COL6 = rep("-9",n())
        )
      
      write_delim(x = tped, path = stri_paste(filename, "_imputed", ".tped", sep = ""), col_names = FALSE, delim = " ")
      write_delim(x = tfam, path = stri_paste(filename, "_imputed", ".tfam", sep = ""), col_names = FALSE, delim = " ") 
    } # end plink imputed
    
    
  } # end plink
  
  if (is.null(plink)) {
    message("PLINK output not selected")
  }
  
  # genind object for adegenet ---------------------------------------------------
  if (genind == TRUE) {
    message("Generating the adegenet genind object")
    genind.prep <- input %>%
      select(MARKERS, INDIVIDUALS, POP_ID, GT) %>% 
      mutate(
        GT = as.character(GT),
        GT = stri_replace_all_fixed(str = GT, pattern = c("REF_REF", "ALT_ALT", "REF_ALT"), replacement = c("2_0", "0_2", "1_1"), vectorize_all = FALSE),
        GT = replace(GT, which(GT == "0_0"), NA)
      ) %>% 
      arrange(MARKERS, POP_ID) %>%
      tidyr::separate(col = GT, into = c("001", "002"), sep = "_", extra = "drop", remove = TRUE) %>%
      tidyr::gather(key = ALLELES, value = COUNT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% # make tidy
      tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>% 
      group_by(POP_ID, INDIVIDUALS) %>%
      tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>%
      ungroup () %>% 
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
        POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
      ) %>% 
      arrange(POP_ID, INDIVIDUALS)
    
    # genind arguments common to all data.type
    ind <- genind.prep$INDIVIDUALS
    pop <- genind.prep$POP_ID
    genind.df <- genind.prep %>% 
      ungroup() %>% 
      select(-c(INDIVIDUALS, POP_ID))
    rownames(genind.df) <- ind
    loc.names <- colnames(genind.df)
    strata <- genind.prep %>% 
      ungroup() %>% 
      select(INDIVIDUALS, POP_ID) %>% 
      distinct(INDIVIDUALS, POP_ID)
    
    # genind constructor
    prevcall <- match.call()
    genind.data <- adegenet::genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = ~POP_ID/INDIVIDUALS)
    
    if (imputation.method != FALSE) {
      message("Generating the adegenet imputed genind object")
      genind.prep.imp <- input.imp %>%
        select(MARKERS, INDIVIDUALS, POP_ID, GT) %>% 
        mutate(
          GT = as.character(GT),
          GT = stri_replace_all_fixed(str = GT, pattern = c("REF_REF", "ALT_ALT", "REF_ALT"), replacement = c("2_0", "0_2", "1_1"), vectorize_all = FALSE),
          GT = replace(GT, which(GT == "0_0"), NA)
        ) %>%
        arrange(MARKERS, POP_ID) %>%
        tidyr::separate(col = GT, into = c("001", "002"), sep = "_", extra = "drop", remove = TRUE) %>%
        tidyr::gather(key = ALLELES, value = COUNT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% # make tidy
        tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>% 
        group_by(POP_ID, INDIVIDUALS) %>%
        tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>%
        ungroup () %>%
        mutate(
          INDIVIDUALS = as.character(INDIVIDUALS),
          POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
          POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
          ) %>% 
        arrange(POP_ID, INDIVIDUALS)
      
      # genind arguments common to all data.type
      ind <- genind.prep.imp$INDIVIDUALS
      pop <- genind.prep.imp$POP_ID
      genind.df.imp <- genind.prep.imp %>% 
        ungroup() %>% 
        select(-c(INDIVIDUALS, POP_ID))
      rownames(genind.df.imp) <- ind
      loc.names <- colnames(genind.df.imp)
      strata <- genind.prep.imp %>% 
        ungroup() %>% 
        select(INDIVIDUALS, POP_ID) %>% 
        distinct(INDIVIDUALS, POP_ID)
      
      # genind constructor
      prevcall <- match.call()
      genind.data.imp <- adegenet::genind(tab = genind.df.imp, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = ~POP_ID/INDIVIDUALS)
    } # end genind imp

  } # end genind
  
  if (is.null(genind)) {
    genind.data <- "not selected"
    genind.data.imp <- "not selected"
  }
  
  # Results ----------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stri_paste("The number of markers removed by the filters:\nSNP: ", snp.before.filters - n_distinct(input$MARKERS), "\nLOCUS: ", locus.before.filters - n_distinct(input$LOCUS)))
  message("The number of markers before -> after the filters")
  message(stri_paste("SNP: ", snp.before.filters, " -> ", n_distinct(input$MARKERS)))
  message(stri_paste("LOCUS: ", locus.before.filters, " -> ", n_distinct(input$LOCUS)))
  cat("#######################################################################\n")
  
  res <- list()
  res$plot.reproducibility <- reproducibility.plot
  res$plot.coverage <- coverage.plot
  res$plot.call.rate <- call.rate.plot
  res$plot.number.snp.reads <- number.snp.reads.plot
  res$maf.data <- maf.data
  res$whitelist.markers <- whitelist.markers
  res$data.filtered <- input.filtered.df
  res$data.filtered.imputed <- input.imp.df
  res$genind.data <- genind.data
  res$genind.data.imputed <- genind.data.imp
  return(res)
}
