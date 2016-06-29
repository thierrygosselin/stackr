#' @name filter_individual
#' @title Individual filter
#' @description Filter individuals genotyped at a marker from a tidy data 
#' set (long format) of any of these file format: 
#' vcf, plink (tped/tfam), stacks haplotype file, genind, 
#' genepop, data frame in wide format. The function uses 
#' \code{\link[stackr]{tidy_genomic_data}} and 
#' \code{\link[stackr]{read_long_tidy_wide}} to load the file. For filtering
#' you can consider the overall number of individuals (no concept of populations here),
#' or the number of genotyped individuals per pop. The threshold can be a fixed 
#' number of individuals, a proportion or a percentage.

#' @param data 6 options: vcf, plink, genind, genepop, 
#' and a data frame in wide or long/tidy format.
#' The function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}} and 
#' \code{\link[stackr]{tidy_genomic_data}}.

#' @param vcf.metadata (optional, logical) For the VCF file, with 
#' default: \code{vcf.metadata = FALSE}, 
#' only the genotype information is kept.
#' \code{vcf.metadata = TRUE}, the metadata contained in the 
#' \code{FORMAT} field will be kept in the tidy data file. 

#' @param interactive.filter (optional, logical) Do you want the filtering session to 
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is 
#' asked to see figures of distribution before making decisions for filtering.

#' @param approach (character). 
#' The approach to filter a marker, is it based on the overall number of 
#' genotyped individuals (no concept of populations here), 
#' \code{approach = "overall"}, or 
#' the number of genotyped individuals per population \code{approach = "pop"}.
#' Default: \code{approach = "pop"}.

#' @param ind.threshold The individual threshold, proportion, percentage or 
#' number e.g. 0.70, 70, 15.
#' Default: \code{ind.threshold = 0.5}.

#' @param percent Is the threshold a percentage? TRUE or FALSE.
#' This argument is necessary to distinguish percentage from integer individual
#' threshold (e.g. 70 percent or 70 individuals).
#' Default: \code{percent = FALSE}.

#' @param prob.pop.threshold (integer, optional) Useful to incorporate problematic 
#' populations dragging down polymorphism discovery, but still wanted for analysis.
#' If individuals with missing genotypes are not managed upstream, 
#' use this threshold to allow variance in the number of populations passing 
#' the \code{ind.threshold} argument. e.g. with \code{prob.pop.threshold = 2},
#' you tolerate a maximum of 2 populations failing to pass the 
#' \code{ind.threshold}. Manage after the problematic populations 
#' downstream with blacklist of individuals and/or missing data imputations. This
#' argument is not necessary with \code{approach = "overall"}.
#' Default: \code{prob.pop.threshold = 0}. See details for more info.

#' @param filename (optional) Name of the filtered tidy data frame file 
#' written to the working directory (ending with \code{.tsv})
#' Default: \code{filename = NULL}.

#' @param strata (optional for data frame and PLINK files, 
#' required for VCF and haplotypes files) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' Default: \code{strata = NULL}.

#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by chromosome and/or locus and/or by snp.
#' The whitelist is in the working directory (e.g. "whitelist.txt").
#' de novo CHROM column with 'un' need to be changed to 1. 
#' In the VCF, the column ID is the LOCUS identification.
#' Default \code{whitelist.markers = NULL} for no whitelist of markers.

#' @param monomorphic.out (optional) Should the monomorphic 
#' markers present in the dataset be filtered out ? 
#' Default: \code{monomorphic.out = TRUE}.

#' @param blacklist.genotype (optional) Useful to erase genotype with below 
#' average quality, e.g. genotype with more than 2 alleles in diploid likely 
#' sequencing errors or genotypes with poor genotype likelihood or coverage. 
#' The blacklist as a minimum of 2 column headers (markers and individuals). 
#' Markers can be 1 column (CHROM or LOCUS or POS), 
#' a combination of 2 (e.g. CHROM and POS or CHROM and LOCUS or LOCUS and POS) or 
#' all 3 (CHROM, LOCUS, POS) The markers columns must be designated: CHROM (character
#' or integer) and/or LOCUS (integer) and/or POS (integer). The id column designated
#' INDIVIDUALS (character) columns header. The blacklist must be in the working 
#' directory (e.g. "blacklist.genotype.txt"). For de novo VCF, CHROM column 
#' with 'un' need to be changed to 1. 
#' Default: \code{blacklist.genotype = NULL} for no blacklist of 
#' genotypes to erase.

#' @param snp.ld (optional) \strong{For VCF file only}. 
#' SNP short distance linkage disequilibrium pruning. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 3 options: \code{"random"} selection, \code{"first"} or
#' \code{"last"} SNP on the same short read/haplotype. For long distance linkage
#' disequilibrium pruning, see details below.
#' Default: \code{snp.ld = NULL}.

#' @param common.markers (optional) Logical. \code{common.markers = TRUE}, 
#' will only keep markers in common (genotyped) between all the populations. 
#' Default: \code{common.markers = FALSE}

#' @param max.marker An optional integer useful to subsample marker number in 
#' large PLINK file. e.g. if the data set 
#' contains 200 000 markers and \code{max.marker = 10000} 10000 markers are
#' subsampled randomly from the 200000 markers. Use \code{whitelist.markers} to
#' keep specific markers.
#' Default: \code{max.marker = NULL}.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

#' @param pop.levels (optional, string) This refers to the levels in a factor. In this 
#' case, the id of the pop.
#' Use this argument to have the pop ordered your way instead of the default 
#' alphabetical or numerical order. e.g. \code{pop.levels = c("QUE", "ONT", "ALB")} 
#' instead of the default \code{pop.levels = c("ALB", "ONT", "QUE")}. 
#' Default: \code{pop.levels = NULL}.

#' @param pop.labels (optional, string) Use this argument to rename/relabel
#' your pop or combine your pop. e.g. To combine \code{"QUE"} and \code{"ONT"} 
#' into a new pop called \code{"NEW"}:
#' (1) First, define the levels for your pop with \code{pop.levels} argument: 
#' \code{pop.levels = c("QUE", "ONT", "ALB")}. 
#' (2) then, use \code{pop.labels} argument: 
#' \code{pop.levels = c("NEW", "NEW", "ALB")}.#' 
#' To rename \code{"QUE"} to \code{"TAS"}:
#' \code{pop.labels = c("TAS", "ONT", "ALB")}.
#' Default: \code{pop.labels = NULL}. If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.

#' @param pop.select (string, optional) Selected list of populations for 
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#'and \code{ONT} population samples (out of 20 pops).
# Default: \code{pop.select = NULL} 


#' @rdname filter_individual
#' @export
#' @import stringi
#' @import dplyr
#' @import readr

#' @details 
#' \strong{Interactive version}
#' 
#' There is 2 steps in the interactive version to visualize and filter
#' the data based on the number of genotyped individuals:
#' 
#' Step 1. Impact of individual threshold on marker discovery
#' 
#' Step 2. Choose the filtering approach and thresholds
#' 
#' 
#' 
#' \strong{prob.pop.threshold}
#' 
#' Used to my functions... then you've seen the \code{pop.num.threshold} here and there.
#' \code{prob.pop.threshold}, is different because the number of populations 
#' genotyped vary across markers, which make difficult the use 
#' of \code{pop.num.threshold}. 
#' e.g. If only 2 populations out of 10 are genotyped for a marker and you use 
#' \code{pop.num.threshold = 0.6} this could lead to inconsistensis. Thus,
#' it's easier to use the number of population you are willing to accept 
#' failling to conform to the threshold. 


#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.ind
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $strata
#' \item $filters.parameters
#' }
#' 
#' With \code{interactive.filter = TRUE}, a list with 2 additionnal objects is created.
#' \enumerate{
#' \item $plot.ind.threshold
#' \item $ind.threshold.helper.table
#' }
#' 
#' The object can be isolated in separate object outside the list by 
#' following the example below.

#' @examples
#' \dontrun{
#' turtle.ind <- filter_individual(
#' data = "turtle.vcf",
#' strata = "turtle.strata.tsv",
#' approach = "pop",
#' ind.thresholds = 50,
#' percent = TRUE,
#' prob.pop.threshold = 2,
#' filename = "tidy.data.turtle.tsv"
#' ) 

#' 
#' #If interactive.filter = TRUE, a list is created and to view the filtered tidy data:
#' tidy.data <- turtle.ind$tidy.filtered.ind
#' 
#' #Inside the same list, to isolate the blacklist.genotypes:
#' bg <- turtle.ind$blacklist.genotypes
#' 
#' # The remaining argument are used in tidy_genomic_data during import and allow
#' # the user to apply some filtering or selection before doing the filtering.
#' }

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("GENOTYPED", "IND_THRESHOLD", "N_IND", "PERCENT", "POP_GENOTYPED", "PROB_POP")
  )
}

filter_individual <- function(
  data,
  vcf.metadata = FALSE,
  interactive.filter = TRUE,
  approach = "pop",
  ind.threshold = 0.5,
  prob.pop.threshold = 0,
  percent = FALSE,
  filename = NULL, 
  blacklist.id = NULL, 
  blacklist.genotype = NULL, 
  whitelist.markers = NULL, 
  monomorphic.out = TRUE, 
  max.marker = NULL,
  snp.ld = NULL, 
  common.markers = FALSE,
  strata = NULL, 
  pop.levels = NULL, 
  pop.labels = NULL, 
  pop.select = NULL
) {
  
  
  cat("#######################################################################\n")
  cat("##################### stackr::filter_individual #######################\n")
  cat("#######################################################################\n")
  
  # manage missing arguments -----------------------------------------------------  
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter == TRUE){
    message("Interactive mode: on")
    message("2 steps to visualize and filter the data based on the number of genotyped individuals:")
    message("Step 1. Impact of individual threshold on marker discovery")
    message("Step 2. Choose the filtering approach and thresholds")
    
    # Folder -------------------------------------------------------------------
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    path.folder <- stri_join(getwd(),"/", "filter_individual_", file.date, sep = "")
    dir.create(file.path(path.folder))
    
    message(stri_join("Folder created: \n", path.folder))
    file.date <- NULL #unused object
  } else {
    path.folder <- getwd()
  }
  
  # Filter parameter file ------------------------------------------------------
  message("Parameters used in this run will be store in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if(length(filters.parameters) == 0) {
    filters.parameters <- data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("Created a file to store filters parameters: filters_parameters.tsv")
  } else {
    message("Using the filters parameters file found in the directory: \nfilters_parameters.tsv")
  }
  
  # File type detection ********************************************************
  if(is.genind(data)){
    data.type <- "genind.file"
    message("File type: genind object")
  } else {
    data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
    if (identical(data.type, "##fileformat=VCF") | stri_detect_fixed(str = data, pattern = ".vcf")) {
      data.type <- "vcf.file"
      message("File type: VCF")
    }
    if (stri_detect_fixed(str = data, pattern = ".tped")) {
      data.type <- "plink.file"
      message("File type: PLINK")
      if (!file.exists(stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        stop("Missing tfam file with the same prefix as your tped")
      }
    } 
    if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stri_detect_fixed(str = data.type, pattern = "MARKERS")) {
      data.type <- "df.file"
      message("File type: data frame of genotypes")
    }
    
    
    if (stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      data.type <- "haplo.file"
      message("File type: haplotypes from stacks")
      # if (is.null(blacklist.genotype)) {
      #   stop("blacklist.genotype file missing. 
      #        Use stackr's missing_genotypes function to create this blacklist")
      # }
    }
    if (stri_detect_fixed(str = data, pattern = ".gen")) {
      data.type <- "genepop.file"
      message("File type: genepop")
    } 
  } # end file type detection
  
  # import data ----------------------------------------------------------------
  if (is.vector(data)){
    message("Using input file in your directory")
    
    input <- stackr::tidy_genomic_data(
      data = data, 
      vcf.metadata = vcf.metadata,
      blacklist.id = blacklist.id, 
      blacklist.genotype = blacklist.genotype, 
      whitelist.markers = whitelist.markers, 
      monomorphic.out = monomorphic.out, 
      max.marker = max.marker,
      snp.ld = snp.ld, 
      common.markers = common.markers,
      strata = strata, 
      pop.levels = pop.levels, 
      pop.labels = pop.labels, 
      pop.select = pop.select,
      filename = NULL
    )
    
    # create a strata.df
    strata.df <- input %>% 
      select(INDIVIDUALS, POP_ID) %>% 
      distinct(INDIVIDUALS, .keep_all = TRUE)
    strata <- strata.df
    pop.levels <- levels(input$POP_ID)
    pop.labels <- pop.levels
  } else {
    message("Using input file from your global environment")
    
    input <- read_long_tidy_wide(data = data)
    strata.df <- input %>% 
      select(INDIVIDUALS, POP_ID) %>% 
      distinct(INDIVIDUALS, .keep_all = TRUE)
    strata <- strata.df
    pop.levels <- levels(input$POP_ID)
    pop.labels <- pop.levels
  }
  
  # prepare filter, table and figure----------------------------------------------
  ind.total <- n_distinct(input$INDIVIDUALS) # total number of individuals
  pop.number <- n_distinct(input$POP_ID) # number of pop
  
  # individuals per pop
  ind.pop <- input %>% 
    distinct(POP_ID, INDIVIDUALS) %>% 
    group_by(POP_ID) %>% 
    tally %>% 
    rename(N_IND = n)
  
  # input genotyped
  input.genotyped <- input %>% 
    filter(GT != "000000")
  
  # overall genotyped individuals
  overall <- input.genotyped %>% 
    select(MARKERS, INDIVIDUALS) %>% 
    group_by(MARKERS) %>% 
    tally %>% 
    rename(GENOTYPED = n) %>% 
    mutate(PERCENT = ceiling(GENOTYPED/ind.total*100))
  
  # number of pop. genotyped per marker
  pop.genotyped.marker <- input.genotyped %>% 
    distinct(MARKERS, POP_ID) %>% 
    group_by(MARKERS) %>% 
    tally %>% 
    rename(POP_GENOTYPED = n)
  
  # genotyped individuals per pop
  pop <- input.genotyped %>% 
    select(MARKERS, INDIVIDUALS, POP_ID) %>% 
    group_by(MARKERS, POP_ID) %>% 
    tally %>% 
    rename(GENOTYPED = n) %>%
    inner_join(ind.pop, by = "POP_ID") %>%
    mutate(PERCENT = ceiling(GENOTYPED/N_IND*100))
  
  input.genotyped <- NULL # unused object
  
  # Step 1. Impact of individual threshold on marker discovery------------------
  if (interactive.filter == TRUE) {
    message("Step 1. Impact of individual threshold on marker discovery")
    message("Show the line graph to inspect the change in the number of markers in relation to the individual percentage thresholds (y/n)): ")
    line.graph <- as.character(readLines(n = 1))
    if (line.graph == "y") {
      threshold.helper.overall <- overall %>% 
        ungroup() %>% 
        summarise(
          `10` = length(PERCENT[PERCENT >= 10]),
          `20` = length(PERCENT[PERCENT >= 20]),
          `30` = length(PERCENT[PERCENT >= 30]),
          `40` = length(PERCENT[PERCENT >= 40]),
          `50` = length(PERCENT[PERCENT >= 50]),
          `60` = length(PERCENT[PERCENT >= 60]),
          `70` = length(PERCENT[PERCENT >= 70]),
          `80` = length(PERCENT[PERCENT >= 80]),
          `90` = length(PERCENT[PERCENT >= 90]),
          `100` = length(PERCENT[PERCENT == 100])
        ) %>% 
        tidyr::gather(data = ., key = IND_THRESHOLD, value = MARKER_NUMBER) %>% 
        mutate(POP_ID = rep("OVERALL", n()))
      
      threshold.helper.pop <- pop %>%
        group_by(POP_ID) %>% 
        summarise(
          `10` = length(PERCENT[PERCENT >= 10]),
          `20` = length(PERCENT[PERCENT >= 20]),
          `30` = length(PERCENT[PERCENT >= 30]),
          `40` = length(PERCENT[PERCENT >= 40]),
          `50` = length(PERCENT[PERCENT >= 50]),
          `60` = length(PERCENT[PERCENT >= 60]),
          `70` = length(PERCENT[PERCENT >= 70]),
          `80` = length(PERCENT[PERCENT >= 80]),
          `90` = length(PERCENT[PERCENT >= 90]),
          `100` = length(PERCENT[PERCENT == 100])
        ) %>% 
        tidyr::gather(data = ., key = IND_THRESHOLD, value = MARKER_NUMBER, -POP_ID)
      
      mean.pop <- threshold.helper.pop %>% 
        group_by(IND_THRESHOLD) %>% 
        summarise(
          MARKER_NUMBER = round(mean(MARKER_NUMBER), 0)
        ) %>% 
        mutate(POP_ID = rep("MEAN_POP", n()))
      
      threshold.helper <- bind_rows(threshold.helper.pop, mean.pop, threshold.helper.overall) %>% 
        mutate(
          IND_THRESHOLD = as.numeric(IND_THRESHOLD),
          POP_ID = factor(POP_ID, levels = c(levels(input$POP_ID), "MEAN_POP", "OVERALL"), ordered = TRUE)
        )
      
      # Set the breaks for the figure
      max.markers <- n_distinct(input$MARKERS)
      # max.markers <- 658
      if (max.markers >= 1000) {
        y.breaks.by <- plyr::round_any(max.markers/10, 100, ceiling)
        y.breaks.max <- plyr::round_any(max.markers, 1000, ceiling)
        y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
      } else {
        y.breaks.by <- plyr::round_any(max.markers/10, 10, ceiling)
        y.breaks.max <- plyr::round_any(max.markers, 100, ceiling)
        y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
      }
      
      plot.ind.threshold <- ggplot(threshold.helper, aes(x = IND_THRESHOLD, y = MARKER_NUMBER)) +
        geom_line() +
        geom_point(size = 2, shape = 21, fill = "white") +
        scale_x_continuous(name = "Individual threshold (percent)", breaks = seq(10, 100, by = 10)) +
        scale_y_continuous(name = "Number of markers", breaks = y.breaks, limits = c(0, y.breaks.max)) +
        theme(
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        ) +
        facet_grid(~POP_ID)
      print(plot.ind.threshold)
      # save
      ggsave(stri_paste(path.folder, "/plot.ind.threshold.pdf"), width = pop.number*4, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stri_paste(path.folder, "/plot.ind.threshold.png"), width = pop.number*4, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the line graph of individual threshold on marker discovery were saved in this directory: \n", path.folder))
      
      # Helper table for individual thresholds
      ind.threshold.helper.table <- threshold.helper %>% 
        group_by(POP_ID) %>% 
        spread(data = ., key = IND_THRESHOLD, MARKER_NUMBER)
      
      message("A table (ind.threshold.helper.table) to help you view the relation between individual thresholds and marker discovery, with your dataset:")
      print(ind.threshold.helper.table)
      message("First column: POP_ID with remaining columns the individual thresholds in percent with the value = the number of markers discovered")
      message("The last 2 rows, MEAN_POP is the mean across your populations and OVERALL is if you had 1 large population")
      write_tsv(x = ind.threshold.helper.table, path = stri_paste(path.folder, "/", "ind.threshold.helper.table.tsv"))
      message(stri_paste("ind.threshold.helper.table was written in this directory: \n", path.folder))
    }
  }
  
  # Step 2. Choose the filtering approach and thresholds------------------------
  # 2 approach: filtering with the overall n. ind. ("overall") or by pop ("pop")
  if (interactive.filter == TRUE) {
    message("Step 2. Choose the filtering approach and thresholds.
The approach to filter a marker, is it based on the overall number of genotyped individuals or 
the number of genotyped individuals per pop ? (overall or pop):")
    approach <- as.character(readLines(n = 1))
    
    message("Enter the individual threshold (number, proportion or percentage).
e.g. enter 10 (for 10 individuals), 0.8 for proportion and 80 for percentage.")
    ind.threshold <- as.numeric(readLines(n = 1))
    
    message("The value you just entered for the individual threshold, is it a percentage? (TRUE/FALSE)")
    percent <- as.character(readLines(n = 1))
    
    if (approach == "pop") {
      message("Tolerance for deviation: look at the plot produced ealier and if you see some populations dragging down
the number of markers for certain percentage thresholds, you have 3 options:\n
1. remove the population (use pop.select argument to keep the desired populations)
2. remove individuals with too many missing genotypes (use blacklist.id argument)
3. use the next threshold to allow variance. Manage the missing values with blacklist of individuals and/or
missing data imputations.\n
Enter the number of problematic population that you allow to deviate from the threshold:")
      prob.pop.threshold <- as.numeric(readLines(n = 1))
    }
  }
  
  # Filtering ------------------------------------------------------------------
  
  # some discrepencies need to be highllighted here. If you have entire pop not genotyped for a markers
  # this will compute them when doing the filtering:
  # summarise(GENOTYPED = length(INDIVIDUALS[GT != "000000"]))
  # if we first remove the individuals not genotyped with :
  # filter(GT != "000000")
  # the pop not genotyped are not accounted for. And this is what we want here.
  # filter_populations take care the ungenotyped pop and common.markers make sure that for certain analysis 
  # you can have common markers or not.
  # so here we focus on when pop got a marker is it at 50% 60% 70%  ... genotyped?
  if (approach == "overall") {
    if (stri_detect_fixed(percent, "F") & ind.threshold > 1) {
      # ind.threshold <- 100
      threshold.id <- "(Using a fixed threshold)"
      prob.pop.threshold <- "NA"
      approach <- "overall individuals (no pop)"
      filter <- overall %>% 
        filter(GENOTYPED >= ind.threshold) %>% 
        distinct(MARKERS)
    } else if (stri_detect_fixed(percent, "F") & stri_detect_fixed(ind.threshold, ".") & ind.threshold < 1) {
      # ind.threshold <- 0.7
      threshold.id <- "(proportion)"
      prob.pop.threshold <- "NA"
      approach <- "overall individuals (no pop)"
      filter <- overall %>% 
        filter(PERCENT >= ind.threshold*100) %>% 
        distinct(MARKERS)
    } else { # percent
      # ind.threshold <- 70
      threshold.id <- "(percent)"
      prob.pop.threshold <- "NA"
      approach <- "overall individuals (no pop)"
      filter <- overall %>% 
        filter(PERCENT >= ind.threshold) %>% 
        distinct(MARKERS)
    }
  } else { # approach by pop
    if (stri_detect_fixed(percent, "F") & ind.threshold > 1) {
      message("Using a fixed threshold")
      threshold.id <- "(ind.)"
      approach <- "individuals by pop"
      # ind.threshold <- 15
      # prob.pop.threshold <- 3
      filter <- pop %>%
        ungroup() %>% 
        filter(GENOTYPED >= ind.threshold) %>% 
        group_by(MARKERS) %>%
        tally() %>% 
        inner_join(pop.genotyped.marker, by = "MARKERS") %>% 
        mutate(PROB_POP = POP_GENOTYPED-n) %>% 
        filter(PROB_POP <= prob.pop.threshold) %>% 
        distinct(MARKERS)
    } else if (stri_detect_fixed(percent, "F") & stri_detect_fixed(ind.threshold, ".") & ind.threshold < 1) {
      message("Using a proportion threshold...")
      threshold.id <- "(proportion)"
      approach <- "individuals by pop"
      # ind.threshold <- 0.6
      # prob.pop.threshold <- 3
      filter <- pop %>% 
        ungroup() %>% 
        filter(PERCENT >= ind.threshold*100) %>% 
        group_by(MARKERS) %>%
        tally() %>% 
        inner_join(pop.genotyped.marker, by = "MARKERS") %>% 
        mutate(PROB_POP = POP_GENOTYPED-n) %>% 
        filter(PROB_POP <= prob.pop.threshold) %>% 
        distinct(MARKERS)
    } else { 
      message("Using a percentage threshold...")
      threshold.id <- "(percent)"
      approach <- "individuals by pop"
      # ind.threshold <- 0.6
      # prob.pop.threshold <- 3
      filter <- pop %>% 
        ungroup() %>% 
        filter(PERCENT >= ind.threshold) %>% 
        group_by(MARKERS) %>%
        tally() %>% 
        inner_join(pop.genotyped.marker, by = "MARKERS") %>% 
        mutate(PROB_POP = POP_GENOTYPED-n) %>% 
        filter(PROB_POP <= prob.pop.threshold) %>% 
        distinct(MARKERS)
    }
  }
  
  # Apply the filter to the tidy data
  filter <- left_join(x = filter, input, by = "MARKERS")
  
  
  # Update filters.parameters SNP ----------------------------------------------
  if (data.type == "vcf.file") {
    # markers.chrom.locus.pos <- input %>% distinct(MARKERS, LOCUS)
    # filter <- left_join(filter, markers.chrom.locus.pos, by = "MARKERS")
    snp.before <-as.integer(n_distinct(input$MARKERS))
    snp.after <-as.integer(n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(n_distinct(input$MARKERS) - n_distinct(filter$MARKERS))
    locus.before <-as.integer(n_distinct(input$LOCUS))
    locus.after <-as.integer(n_distinct(filter$LOCUS))
    locus.blacklist <- as.integer(n_distinct(input$LOCUS) - n_distinct(filter$LOCUS))
  } else if (data.type == "haplo.file") {
    snp.before <- as.character("NA")
    snp.after <- as.character("NA")
    snp.blacklist <- as.character("NA")
    locus.before <-as.integer(n_distinct(input$MARKERS))
    locus.after <-as.integer(n_distinct(filter$MARKERS))
    locus.blacklist <- as.integer(n_distinct(input$MARKERS) - n_distinct(filter$MARKERS))
  } else {
    snp.before <-as.integer(n_distinct(input$MARKERS))
    snp.after <-as.integer(n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(n_distinct(input$MARKERS) - n_distinct(filter$MARKERS))
    locus.before <- as.character("NA")
    locus.after <- as.character("NA")
    locus.blacklist <- as.character("NA")
  }
  
  markers.before <- stri_paste(snp.before, locus.before, sep = "/")
  markers.after <- stri_paste(snp.after, locus.after, sep = "/")
  markers.blacklist <- stri_paste(snp.blacklist, locus.blacklist, sep = "/")
  
  filters.parameters <- data_frame(
    FILTERS = c("Individuals", rep(as.character(""), 2)),
    PARAMETERS = c("approach", "ind.threshold", "prob.pop.threshold"), 
    VALUES = c(approach, paste(">=", ind.threshold, " ", threshold.id), prob.pop.threshold), 
    BEFORE = c("", "", markers.before),
    AFTER = c("", "", markers.after),
    BLACKLIST = c("", "", markers.blacklist),
    UNITS = c("", "genotyped", "SNP/LOCUS"),
    COMMENTS = c("", "", "")
  )
  write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  
  # saving tidy data 
  if (!is.null(filename)) {
    message("Writing the filtered tidy data set in your working directory...")
    # if (!is.null(save.feather)) {
    # feather::write_feather(filter, stri_replace_all_fixed(filename, pattern = ".tsv", replacement = "_feather.tsv", vectorize_all = TRUE))
    # } else {
    write_tsv(filter, filename, append = FALSE, col_names = TRUE)
    # }
  }
  
  # saving whitelist
  message("Writing the whitelist of markers in your working directory\nwhitelist.markers.ind.tsv")
  
  if (data.type == "vcf.file") {
    whitelist.markers <- filter %>% 
      ungroup() %>%
      distinct(CHROM, LOCUS, POS)
  } else {
    whitelist.markers <- filter %>% 
      ungroup() %>%
      select(LOCUS = MARKERS) %>% 
      distinct(LOCUS)
  }
  write_tsv(whitelist.markers, "whitelist.markers.ind.tsv", append = FALSE, col_names = TRUE)
  
  
  # saving blacklist
  message("Writing the blacklist of markers in your working directory\nblacklist.markers.ind.tsv")
  if (data.type == "vcf.file") {
    blacklist.markers <- input %>% 
      ungroup() %>%
      distinct(CHROM, LOCUS, POS) %>% 
      anti_join(whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  } else {
    blacklist.markers <- input %>% 
      ungroup() %>%
      select(LOCUS = MARKERS) %>% 
      distinct(LOCUS) %>% 
      anti_join(whitelist.markers, by = "LOCUS")
  }
  write_tsv(blacklist.markers, "blacklist.markers.ind.tsv", append = FALSE, col_names = TRUE)
  
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stri_paste("approach: ", approach))
  message(stri_paste("ind.threshold: ", ">= ", ind.threshold, " ", threshold.id))
  message(stri_paste("prob.pop.threshold: ", prob.pop.threshold))
  if (data.type != "vcf.file" & data.type != "haplo.file") {
    message(stri_paste("The number of markers removed by the Individual filter: ", n_distinct(input$MARKERS)-n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the Individual filter")
    message(stri_paste("SNP: ", as.integer(n_distinct(input$MARKERS)), " -> ", as.integer(n_distinct(filter$MARKERS))))
  } else if (data.type == "vcf.file") {
    message(stri_paste("The number of markers removed by the Individual filter:\nSNP: ", n_distinct(input$POS) - n_distinct(filter$POS), "\nLOCUS: ", n_distinct(input$LOCUS) - n_distinct(filter$LOCUS)))
    message("The number of markers before -> after the Individual filter")
    message(stri_paste("SNP: ", as.integer(n_distinct(input$POS)), " -> ", as.integer(n_distinct(filter$POS))))
    message(stri_paste("LOCUS: ", as.integer(n_distinct(input$LOCUS)), " -> ", as.integer(n_distinct(filter$LOCUS))))
  } else {# for haplotype file
    message(stri_paste("The number of markers/locus removed by the Individual filter: ", n_distinct(input$MARKERS)-n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the Individual filter")
    message(stri_paste("MARKERS/LOCUS: ", as.integer(n_distinct(input$MARKERS)), " -> ", as.integer(n_distinct(filter$MARKERS))))
  }
  cat("#######################################################################\n")
  res <- list()
  res$tidy.filtered.ind <- filter
  res$whitelist.markers <- whitelist.markers
  res$blacklist.markers <- blacklist.markers
  res$strata <- strata
  res$filters.parameters <- filters.parameters
  if (interactive.filter){
    res$plot.ind.threshold <- plot.ind.threshold
    res$ind.threshold.helper.table <- ind.threshold.helper.table
  }
  return(res)
}
