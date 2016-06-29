#' @name filter_population
#' @title Population filter
#' @description Filter markes based on populations genotyped. Use a tidy data 
#' set (long format) of any of these file format: 
#' vcf, plink (tped/tfam), stacks haplotype file, genind, 
#' genepop, data frame in wide format. The function uses 
#' \code{\link[stackr]{tidy_genomic_data}} and 
#' \code{\link[stackr]{read_long_tidy_wide}} to load the file. For filtering
#' The threshold can be a fixed number of population, a proportion or a percentage.

#' @param data 6 options: vcf, plink, genind, genepop, 
#' and a data frame in wide or long/tidy format.
#' The function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}} and 
#' \code{\link[stackr]{tidy_genomic_data}}.

#' @param interactive.filter (optional, logical) Do you want the filtering session to 
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is 
#' asked to see figures of distribution before making decisions for filtering.

#' @param pop.threshold The population threshold, proportion, percentage or 
#' number e.g. 0.70, 70, 15.
#' Default: \code{pop.threshold = 100}.

#' @param percent Is the threshold a percentage? TRUE or FALSE.
#' This argument is necessary to distinguish percentage from integer population
#' threshold (e.g. 50 percent or 50 populations).
#' Default: \code{percent = TRUE}.

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


#' @rdname filter_population
#' @export
#' @import stringi
#' @import dplyr
#' @import readr

#' @details 
#' \strong{Interactive version}
#' 
#' There is 2 steps in the interactive version to visualize and filter
#' the data based on the population representation:
#' 
#' Step 1. Impact of population threshold on marker discovery
#' 
#' Step 2. Choose the filtering threshold


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
#' \item $plot.pop.threshold
#' \item $pop.threshold.helper.table
#' }
#' 
#' The object can be isolated in separate object outside the list by 
#' following the example below.

#' @examples
#' \dontrun{
#' turtle.pop <- filter_population(
#' data = "turtle.vcf",
#' strata = "turtle.strata.tsv",
#' pop.thresholds = 100,
#' percent = TRUE,
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

filter_population <- function(
  data,  
  interactive.filter = TRUE,
  pop.threshold = 100,
  percent = TRUE,
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
  cat("##################### stackr::filter_population #######################\n")
  cat("#######################################################################\n")
  
  # manage missing arguments -----------------------------------------------------  
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter == TRUE){
    message("Interactive mode: on")
    message("2 steps to visualize and filter the data based on the number of genotyped individuals:")
    message("Step 1. Impact of population threshold on marker discovery")
    message("Step 2. Choose the filtering threshold")
    
    # Folder -------------------------------------------------------------------
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    path.folder <- stri_join(getwd(),"/", "filter_population_", file.date, sep = "")
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
      vcf.metadata = FALSE,
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
  filter.prep <- input %>% 
    filter(GT != "000000") %>% 
    distinct(MARKERS, POP_ID) %>% 
    group_by(MARKERS) %>% 
    tally %>% 
    rename(POP_GENOTYPED = n)
  
  pop.threshold.helper.table <- filter.prep %>% 
    group_by(POP_GENOTYPED) %>% 
    summarise(MARKER_NUMBER = length(MARKERS))
  
  
  
  # Step 1. Impact of population threshold on marker discovery------------------
  if (interactive.filter == TRUE) {
    message("Step 1. Impact of population threshold on marker discovery")
    message("Show the line graph to inspect the change in the number of markers in relation to the population percentage thresholds (y/n)): ")
    line.graph <- as.character(readLines(n = 1))
    if (line.graph == "y") {
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
      
      plot.pop.threshold <- ggplot(pop.threshold.helper.table, aes(x = POP_GENOTYPED, y = MARKER_NUMBER)) +
        geom_line() +
        geom_point(size = 2, shape = 21, fill = "white") +
        scale_x_continuous(name = "Number of populations genotyped") +
        scale_y_continuous(name = "Number of markers", breaks = y.breaks, limits = c(0, y.breaks.max)) +
        theme(
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )
      # plot.pop.threshold
      print(plot.pop.threshold)
      # save
      ggsave(stri_paste(path.folder, "/plot.pop.threshold.pdf"), width = length(pop.threshold.helper.table$POP_GENOTYPED)*2, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stri_paste(path.folder, "/plot.pop.threshold.png"), width = length(pop.threshold.helper.table$POP_GENOTYPED)*2, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the line graph of populations threshold on marker discovery were saved in this directory: \n", path.folder))
      
      
      # Helper table for individual thresholds
      message("A table (pop.threshold.helper.table) to help you view the relation between population threshold and marker discovery, with your dataset:")
      print(pop.threshold.helper.table)
      write_tsv(x = pop.threshold.helper.table, path = stri_paste(path.folder, "/", "pop.threshold.helper.table.tsv"))
      message(stri_paste("pop.threshold.helper.table was written in this directory: \n", path.folder))
    }
  }
  
  # Step 2. Choose the filtering threshold -------------------------------------
  if (interactive.filter == TRUE) {
    message("Step 2. Choose the filtering threshold.")
    message("Enter the population threshold (number, proportion or percentage).
e.g. enter 10 (for 10 populations), 0.8 for proportion and 80 for percentage.")
    pop.threshold <- as.numeric(readLines(n = 1))
    
    message("The value you just enter for the pop threshold, is it a percentage? (TRUE/FALSE)")
    percent <- as.character(readLines(n = 1))
  }
  
  # Filtering ------------------------------------------------------------------
  pop.number <- n_distinct(input$POP_ID) # number of pop
  
  filter.prep <- filter.prep %>% 
    mutate(PERCENT = ceiling(POP_GENOTYPED/pop.number*100))
  
  if (stri_detect_fixed(percent, "F") & pop.threshold > 1) {
    message("Using a fixed threshold")
    threshold.id <- "(pop.)"
    # pop.threshold <- 12
    filter <- filter.prep %>%
      ungroup() %>% 
      filter(POP_GENOTYPED >= pop.threshold) %>%
      distinct(MARKERS)
  } else if (stri_detect_fixed(percent, "F") & stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    message("Using a proportion threshold...")
    threshold.id <- "(proportion)"
    # pop.threshold <- 0.6
    filter <- filter.prep %>% 
      ungroup() %>% 
      filter(PERCENT >= pop.threshold*100) %>%
      distinct(MARKERS)
  } else { 
    message("Using a percentage threshold...")
    threshold.id <- "(percent)"
    # pop.threshold <- 100
    filter <- filter.prep %>% 
      ungroup() %>% 
      filter(PERCENT >= pop.threshold) %>% 
      distinct(MARKERS)
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
    FILTERS = c("Populations"),
    PARAMETERS = c("pop.threshold"), 
    VALUES = c(stri_paste(pop.threshold, " ", threshold.id)), 
    BEFORE = c(markers.before),
    AFTER = c(markers.after),
    BLACKLIST = c(markers.blacklist),
    UNITS = c("SNP/LOCUS"),
    COMMENTS = c("")
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
  message("Writing the whitelist of markers in your working directory\nwhitelist.markers.pop.tsv")
  
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
  write_tsv(whitelist.markers, "whitelist.markers.pop.tsv", append = FALSE, col_names = TRUE)
  
  
  # saving blacklist
  message("Writing the blacklist of markers in your working directory\nblacklist.markers.pop.tsv")
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
  write_tsv(blacklist.markers, "blacklist.markers.pop.tsv", append = FALSE, col_names = TRUE)
  
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stri_paste("pop.threshold: ", ">= ", pop.threshold))
  if (data.type != "vcf.file" & data.type != "haplo.file") {
    message(stri_paste("The number of markers removed by the Population filter: ", n_distinct(input$MARKERS)-n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the Population filter")
    message(stri_paste("SNP: ", as.integer(n_distinct(input$MARKERS)), " -> ", as.integer(n_distinct(filter$MARKERS))))
  } else if (data.type == "vcf.file") {
    message(stri_paste("The number of markers removed by the Population filter:\nSNP: ", n_distinct(input$POS) - n_distinct(filter$POS), "\nLOCUS: ", n_distinct(input$LOCUS) - n_distinct(filter$LOCUS)))
    message("The number of markers before -> after the Population filter")
    message(stri_paste("SNP: ", as.integer(n_distinct(input$POS)), " -> ", as.integer(n_distinct(filter$POS))))
    message(stri_paste("LOCUS: ", as.integer(n_distinct(input$LOCUS)), " -> ", as.integer(n_distinct(filter$LOCUS))))
  } else {# for haplotype file
    message(stri_paste("The number of markers/locus removed by the Population filter: ", n_distinct(input$MARKERS)-n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the Population filter")
    message(stri_paste("MARKERS/LOCUS: ", as.integer(n_distinct(input$MARKERS)), " -> ", as.integer(n_distinct(filter$MARKERS))))
  }
  cat("#######################################################################\n")
  res <- list()
  res$tidy.filtered.pop <- filter
  res$whitelist.markers <- whitelist.markers
  res$blacklist.markers <- blacklist.markers
  res$strata <- strata
  res$filters.parameters <- filters.parameters
  if (interactive.filter){
    res$plot.pop.threshold <- plot.pop.threshold
    res$pop.threshold.helper.table <- pop.threshold.helper.table
  }
  return(res)
}
