# Minor Allele Frequency
#' @title MAF filter
#' @description Minor Allele Frequency filter from a tidy data set (long format)
#' of any of these file format: 
#' vcf, plink (tped/tfam), stacks haplotype file, genind, 
#' genepop, data frame in wide format. The function uses 
#' \code{\link[stackr]{tidy_genomic_data}} and 
#' \code{\link[stackr]{read_long_tidy_wide}} to load the file.

#' @param data 6 options: vcf, plink, genind, genepop, 
#' and a data frame in wide or long/tidy format.
#' The function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}} and 
#' \code{\link[stackr]{tidy_genomic_data}}.

#' @param interactive.filter (optional, logical) Do you want the filtering session to 
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is 
#' asked to see figures of distribution before making decisions for filtering.

#' @param maf.approach (character, optional). 
#' \code{maf.approach = "haplotype"} : looks at the minimum MAF found on the 
#' read/haplotype. Using this option will discard all the markers/snp on 
#' that read based on the thresholds chosen. This method is only available 
#' for VCF and haplotype files, or tidy data frame from those file types.
#' \code{maf.approach = "SNP"} : treats all the SNP on the same 
#' haplotype/read as independent. Doesn't work with haplotype file, 
#' but does work for all other file type.
#' Default is \code{maf.approach = "SNP"}.

#' @param maf.thresholds (string, double) String with 
#' local/populations and global/overall maf thresholds, respectively.
#' e.g. \code{maf.thresholds = c(0.05, 0.1)} for a local maf threshold 
#' of 0.05 and a global threshold of 0.1.

#' @param maf.operator (character, optional) \code{maf.operator = "AND"} or 
#' default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?


#' @param maf.pop.num.threshold (integer, optional) When maf thresholds are used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. 
#' Default: \code{maf.pop.num.threshold = 1}

#' @param filename (optional) Name of the filtered tidy data frame file 
#' written to the working directory (ending with \code{.tsv})
#' Default: \code{filename = NULL}.

# @param save.feather (optional) Use the package
# \href{https://github.com/wesm/feather}{feather} to save the data frame (very fast).
# Default: \code{save.feather = NULL}.

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

#' @param common.markers (optional) Logical. Default: \code{common.markers = TRUE}, 
#' will only keep markers in common (genotyped) between all the populations.

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



#' @rdname filter_maf
#' @export
#' @import stringi
#' @import dplyr
#' @import readr
# @importFrom feather write_feather

#' @details To help choose a threshold for the local and global MAF, look
#' at the function \link{diagnostic_maf}, or use the interactive version.
#' 
#' There is 3 steps in the interactive version:
#' 
#' Step 1. Global MAF: Inspecting the MAF globally
#' 
#' Step 2. Local MAF: Inspecting the MAF at the populations level
#' 
#' Step 3. Filtering markers based on the different MAF arguments

#' @return With \code{interactive.filter = FALSE}, a list in the global environment,
#' with 7 objects:
#' \enumerate{
#' \item $tidy.filtered.maf
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $blacklist.markers
#' \item $maf.data.thresholds
#' \item $strata
#' \item $filters.parameters
#' }
#' 
#' With \code{interactive.filter = TRUE}, a list with 5 additionnal objects is created.
#' \enumerate{
#' \item $violinplot.maf.global <- violinplot.maf.global
#' \item $violinplot.maf.local
#' \item $plot.distribution.maf.locall
#' \item $maf.global.summary
#' \item $maf.helper.table
#' }
#' 
#' The object can be isolated in separate object outside the list by 
#' following the example below.

#' @examples
#' \dontrun{
#' # The minumum
#' turtle.maf <- filter_maf(
#' data = "turtle.vcf",
#' strata = "turtle.strata.tsv",
#' maf.thresholds = c(0.04, 0.02)
#' ) 
#' #This will use the default: interactive version, 
#' maf.approach = "SNP", 
#' maf.operator = "OR",
#' maf.pop.num.threshold = 1
#' 
#' 
#' #If interactive.filter = TRUE, a list is created and to view the filtered tidy data:
#' tidy.data <- turtle.maf$tidy.filtered.maf
#' 
#' #Inside the same list, to isolate the blacklist.genotypes:
#' bg <- turtle.maf$blacklist.genotypes
#' 
#' # The remaining argument are used in tidy_genomic_data during import and allow
#' # the user to apply some filtering or selection before doing the MAF filtering.
#' }


filter_maf <- function(
  data,
  interactive.filter = TRUE,
  maf.approach = "SNP",
  maf.thresholds,
  maf.operator = "OR",
  maf.pop.num.threshold = 1,
  filename = NULL,
  #save.feather = NULL,
  blacklist.id = NULL, 
  blacklist.genotype = NULL, 
  whitelist.markers = NULL, 
  monomorphic.out = TRUE, 
  max.marker = NULL,
  snp.ld = NULL, 
  common.markers = NULL,
  strata = NULL, 
  pop.levels = NULL, 
  pop.labels = NULL, 
  pop.select = NULL
) {

  
  cat("#######################################################################\n")
  cat("######################### stackr::filter_maf ##########################\n")
  cat("#######################################################################\n")
  
  # manage missing arguments -----------------------------------------------------  
  if (missing(data)) stop("Input file missing")
  if (!interactive.filter & missing(maf.thresholds)) stop("The required maf.thresholds argument is missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter == TRUE){
    message("Interactive mode: on")
    message("3 steps to visualize and filter the data based on MAF:")
    message("Step 1. Global MAF: Inspecting the MAF globally")
    message("Step 2. Local MAF: Inspecting the MAF at the populations level")
    message("Step 3. Filtering markers based on the different MAF arguments")
    
    # Folder -------------------------------------------------------------------
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    path.folder <- stri_join(getwd(),"/", "filter_maf_", file.date, sep = "")
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
  
  if(data.type == "haplo.file") {
    message("With stacks haplotype file the maf.approach is automatically set to: haplotype")
    maf.approach <- "SNP" 
    # confusing, but because the haplotpe file doesn't have snp info, only locus info
    # it's treated as markers/snp info and filtered the same way as the approach by SNP.
    # but it's really by haplotype
  }
  
  if (maf.approach == "haplotype") {
    if (data.type != "vcf.file" & data.type != "haplo.file") {
      stop("The haplotype approach during MAF filtering is for VCF and
         stacks haplotypes file, only. Use the snp approach for the other file types")
    }
  }
  
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
      distinct(INDIVIDUALS)
    strata <- strata.df
    pop.levels <- levels(input$POP_ID)
    pop.labels <- pop.levels
  } else {
    message("Using input file from your global environment")
    
    input <- read_long_tidy_wide(data = data)
    strata.df <- input %>% 
      select(INDIVIDUALS, POP_ID) %>% 
      distinct(INDIVIDUALS)
    strata <- strata.df
    pop.levels <- levels(input$POP_ID)
    pop.labels <- pop.levels
  }
  # # File detection to speed up the MAF calculations with VCF file --------------
  # if ("GT_VCF" %in% colnames(input)) {
  #   data.type <- "vcf.file"
  # } else {
  #   data.type <- NULL
  # }
  
  # Detection if data is summarized for MAF or not -----------------------------
  if ("MAF_LOCAL" %in% colnames(input)) {
    summarize.data <- TRUE
  } else {
    message("Summarizing the data by populations and globally")
    summarize.data <- FALSE
  }
  
  # Summarize data for MAF if required -----------------------------------------
  if (!summarize.data){
    message("Calculating global and local MAF on large data set may take some time...")
    
    if (data.type == "vcf.file") {
      maf.local <- input %>%
        filter(GT_VCF != "./.") %>%
        group_by(MARKERS, POP_ID, REF, ALT) %>%
        summarise(
          N = as.numeric(n()),
          PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
          QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
        ) %>%
        mutate(MAF_LOCAL = ((QQ * 2) + PQ) / (2 * N))
      
      maf.global <- maf.local %>%
        group_by(MARKERS) %>%
        summarise_each_(funs(sum), vars = c("N", "PQ", "QQ")) %>%
        mutate(MAF_GLOBAL = ((QQ * 2) + PQ) / (2 * N)) %>%
        select(MARKERS, MAF_GLOBAL)
      
      maf.data <- maf.global %>%
        left_join(maf.local, by = c("MARKERS")) %>%
        select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
      
      maf.local <- NULL
      maf.global <- NULL
    } else { # not vcf file
      
      # We split the alleles here to prep for MAF
      maf.data <- input %>%
        select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        tidyr::separate(data = ., col = GT, into = .(A1, A2), sep = 3, remove = TRUE) %>% 
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
        filter(GT != "000")
      
      maf.data <- maf.data %>%
        group_by(MARKERS, GT, POP_ID) %>%
        tally %>%
        arrange(MARKERS, GT) %>% 
        group_by(MARKERS, GT) %>%
        mutate(sum.pop = sum(n)) %>% 
        group_by(MARKERS) %>%
        mutate(MAF_GLOBAL = min(sum.pop)/sum(n)) %>% 
        group_by(MARKERS, POP_ID) %>%
        mutate(MAF_LOCAL = n/sum(n)) %>% 
        arrange(MARKERS, POP_ID, GT) %>% 
        group_by(MARKERS, POP_ID) %>% 
        filter(n == min(n)) %>% 
        distinct(MARKERS, POP_ID) %>% 
        select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
    }# end maf calculations
  } # end summarize data
  
  # Step 1. Global MAF: Inspecting the MAF globally----------------------------
  if (interactive.filter == TRUE){
    message("Step 1. Global MAF: Inspecting the MAF globally")
    # plot_1: Violin plot global MAF individuals and pop
    message("Show the violin plot for the global MAF (y/n)): ")
    violinplot <- as.character(readLines(n = 1))
    if (violinplot == "y") {
      
      message("Generating violin plot may take some time...")
      # plot
      OVERALL <- NULL
      global.data <- maf.data %>% 
        group_by(MARKERS) %>% 
        select(MARKERS, MAF_GLOBAL) %>% 
        distinct(MARKERS, MAF_GLOBAL) %>%
        mutate(OVERALL = rep("overall", n()))
      
      maf.global.summary <- global.data %>% 
        ungroup %>% 
        summarise(
          MEAN = mean(MAF_GLOBAL, na.rm = TRUE),
          MEDIAN = stats::median(MAF_GLOBAL, na.rm = TRUE),
          RANGE = stri_paste(round(min(MAF_GLOBAL, na.rm = TRUE), 4), " - ", round(max(MAF_GLOBAL, na.rm = TRUE), 4))
        )
      
      violinplot.maf.global <- ggplot(data = global.data, aes(x = OVERALL, y = MAF_GLOBAL, na.rm = TRUE))+
        geom_violin(trim = TRUE)+
        geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
        stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
        # scale_y_continuous(name = "Global MAF", breaks = c(0, 0.01, 0.02, 0.))
        labs(x = "Overall")+
        labs(y = "Global MAF")+
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.text.x = element_blank(), 
          # axis.text.x = element_text(size = 8, family = "Helvetica"), 
          legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )
      print(violinplot.maf.global)
      # save
      ggsave(stri_paste(path.folder, "/maf.global.pdf"), width = 10, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stri_paste(path.folder, "/maf.global.png"), width = 10, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the violin plot for the global MAF were saved in this directory: \n", path.folder))
      
      global.data <- NULL # unused object
      
      write_tsv(x = maf.global.summary, path = "maf.global.summary.tsv")
      message(stri_paste("The global MAF mean: ", round(maf.global.summary$MEAN, 4)))
      message(stri_paste("The global MAF median: ", round(maf.global.summary$MEDIAN, 4)))
      message(stri_paste("The global MAF range: ", maf.global.summary$RANGE))
      message(stri_paste("maf.global.summary.tsv was saved in this directory: \n", path.folder))
    }
  } # end global maf
  
  # Step 2. Local MAF: Inspecting the MAF at the populations level
  if (interactive.filter == TRUE) {
    message("Step 2. Local MAF: Inspecting the MAF at the population level")
    # plot_2: Violin plot local MAF
    message("Show the violin plot for the local MAF (y/n)): ")
    violinplot <- as.character(readLines(n = 1))
    if (violinplot == "y") {
      pop.number <- length(levels(maf.data$POP_ID))
      
      message("Generating violin plot may take some time...")
      # plot
      # maf.local.summary <- maf.data %>% 
      #   ungroup %>% 
      #   summarise(
      #     MEAN = mean(MAF_LOCAL, na.rm = TRUE),
      #     MEDIAN = stats::median(MAF_LOCAL, na.rm = TRUE),
      #     RANGE = stri_paste(round(min(MAF_LOCAL, na.rm = TRUE), 4), " - ", round(max(MAF_GLOBAL, na.rm = TRUE), 4))
      #   )
      
      violinplot.maf.local <- ggplot(data = maf.data, aes(x = POP_ID, y = MAF_LOCAL, na.rm = TRUE))+
        geom_violin(trim = TRUE)+
        geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
        stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
        labs(x = "Populations/Groupings")+
        labs(y = "Local/populations MAF")+
        theme(
          legend.position = "none",
          axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
          axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
          # axis.text.x = element_blank(), 
          axis.text.x = element_text(size = 8, family = "Helvetica"),
          legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
          strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
        )
      print(violinplot.maf.local)
      # save
      ggsave(stri_paste(path.folder, "/maf.local.violinplot.pdf"), width = pop.number, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stri_paste(path.folder, "/maf.local.violinplot.png"), width = pop.number, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the violin plot for the global MAF were saved in this directory: \n", path.folder))
    }
  }
  
  # plot_3: Distribution of local MAF
  if (interactive.filter == TRUE) {
    message("Show site frequency spectrum (y/n): ")
    spectrum <- as.character(readLines(n = 1))
    if(spectrum == "y") {
      pop.number <- length(levels(maf.data$POP_ID))
      plot.distribution.maf.local <- ggplot(data = maf.data, aes(x = MAF_LOCAL, na.rm = FALSE))+
        geom_line(aes(y = ..scaled.., color = POP_ID), stat = "density", adjust = 1)+ # pop colored
        #   scale_colour_manual(name ="Sampling sites", values = colour_palette_sites.pink)+
        scale_x_continuous(breaks = seq(0,1, by =0.1))+
        # labels = c("0", "0.02", "0.05", "0.10", "0.20", "0.50", "1.00"),
        # limits = c(0,1))+
        labs(x = "Minor Allele Frequency (MAF)")+
        labs(y = "Density of MARKERS (scaled)")+
        expand_limits(y = 0)+
        theme(
          axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
          axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"),
          legend.position = "none",
          axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
          # legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
          # legend.text = element_text(size = 12, family = "Helvetica", face = "bold"),
          strip.text.y = element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"), 
          strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
        )+
        facet_grid(~POP_ID)
      print(plot.distribution.maf.local)
      ggsave(stri_paste(path.folder, "/maf.local.spectrum.pdf"), width = pop.number * 5, height = 15, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stri_paste(path.folder, "/maf.local.spectrum.png"), width = pop.number * 5, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the local maf spectrum plot were saved in this directory: \n", path.folder))
    }
  }
  
  # Helper table for global and local MAF
  if (interactive.filter == TRUE) {
    message("Show and write to directory a helper table for MAF (y/n): ")
    helper.table <- as.character(readLines(n = 1))
    if(helper.table == "y") {
      maf.helper.table <- input %>% 
        group_by(POP_ID, INDIVIDUALS) %>% 
        distinct(POP_ID, INDIVIDUALS) %>%
        ungroup %>% 
        group_by(POP_ID) %>% 
        tally
      
      n.ind <- n_distinct(input$INDIVIDUALS)
      TOTAL <- data_frame(POP_ID = "TOTAL/GLOBAL", n = n.ind)
      maf.helper.table <- bind_rows(maf.helper.table, TOTAL) %>% 
        mutate(
          MAF_1_IND = 1/n,
          MAF_2_IND = 2/n,
          MAF_3_IND = 3/n,
          MAF_4_IND = 4/n,
          MAF_5_IND = 5/n,
          MAF_6_IND = 6/n,
          MAF_7_IND = 7/n,
          MAF_8_IND = 8/n,
          MAF_9_IND = 9/n,
          MAF_10_IND = 10/n
        )
      # remove unused objects
      TOTAL <- NULL
      n.ind <- NULL
      
      message("A table (maf.helper.table) to help you view the LOCAL and GLOBAL MAF in terms of individuals, with your dataset:")
      print(maf.helper.table)
      # View(maf.helper.table)
      message("First and second columns: POP_ID and sample size, with total sample size in the last row. 
The remaining columns are the corresponding local and global (last row) MAF of your dataset
from 1 to 10 individuals genotyped as heterozygote(s) for the minor allele")
      write_tsv(x = maf.helper.table, path = stri_paste(path.folder, "/", "maf.helper.table.tsv"))
      message(stri_paste("maf.helper.table.tsv was written in this directory: \n", path.folder))
    }
  }# end helper table
  
  # Thresholds selection -----------------------------------------------------------------
  
  # maf.approach
  if (interactive.filter) {
    message("Step 3. Filtering markers based on the different MAF arguments")
    message("The maf.approach:\n
maf.approach = \"haplotype\" : looks at the minimum MAF found on the read/haplotype. 
This will discard all the markers/snp on that read based on the thresholds chosen.
This method is only available for VCF and haplotype files or tidy data frame from
those file type.\n
maf.approach = \"SNP\" : treats all the SNP on the same haplotype/read as independent.
Doesn't work with haplotype file, but does work for all other file type.")
    message("Enter the maf.approach:")
    maf.approach <- as.character(readLines(n = 1))
    if (!maf.approach %in% c("SNP", "haplotype")) stop("maf.approach: SNP or haplotype")
  }
  
  # maf.thresholds
  if (interactive.filter) {
    message("The maf.thresholds:\n")
    message("Enter the maf local threshold:")
    maf.local.threshold <- as.character(readLines(n = 1))
    message("Enter the maf global threshold:")
    maf.global.threshold <- as.character(readLines(n = 1))
  }
  
  # maf.operator
  if (interactive.filter) {
    message("The maf.operator:\n 
When filtering, to keep the markers, do you want the local \"AND\" the global MAF to pass the thresholds,
or ... you want the local \"OR\" the global MAF to pass the thresholds?")
    message("Enter the maf.operator:")
    maf.operator <- as.character(readLines(n = 1))
    if (!maf.operator %in% c("OR", "AND")) stop("maf.operator: either OR/AND")
  }
  
  # maf.pop.num.threshold
  if (interactive.filter) {
    message("The maf.pop.num.threshold:\n 
How many populations are required to pass the thresholds to keep the locus?\n
Example: if you have 10 populations and choose maf.pop.num.threshold = 3,
3 populations out of 10 are required to have LOCAL and/or GLOBAL MAF higher than the
thresholds entered.")
    message("Enter the maf.pop.num.threshold:")
    maf.pop.num.threshold <- as.numeric(readLines(n = 1))
  }
  
  if (!interactive.filter) {
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
  }
  
  
  # Update the maf.data with pass or not filter based on threshold -------------
  maf.data.thresholds <- maf.data %>% 
    mutate(
      OR = ifelse((MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold), "pass", "pruned"),
      AND = ifelse((MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold), "pass", "pruned")
    )
  
  write_tsv(x = maf.data.thresholds, 
            path = "maf.data.tsv",
            col_names = TRUE, 
            append = FALSE
  )
  message(stri_paste("The MAF summary statistics (maf.data.tsv), written in this directory: \n", path.folder))
  
  # Filtering ------------------------------------------------------------------
  if (maf.approach == "haplotype") {
    filter <- tidyr::separate(data = maf.data, 
                              col = MARKERS, 
                              into = c("CHROM", "LOCUS", "POS"), 
                              sep = "_", 
                              remove = FALSE, 
                              extra = "warn"
    )
    
    if (maf.operator == "OR") {
      filter <- filter %>%
        group_by(LOCUS, POP_ID) %>%
        summarise(
          MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
          MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
        ) %>%
        filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
        group_by(LOCUS) %>%
        tally() %>%
        filter(n >= maf.pop.num.threshold) %>%
        select(LOCUS) %>%
        left_join(input, by = "LOCUS") %>%
        arrange(LOCUS, POP_ID)
    } else { # AND operator between local and global maf
      filter <- filter %>%
        group_by(LOCUS, POP_ID) %>%
        summarise(
          MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
          MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
        ) %>%
        filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
        group_by(LOCUS) %>%
        tally() %>%
        filter(n >= maf.pop.num.threshold) %>%
        select(LOCUS) %>%
        left_join(input, by = "LOCUS") %>%
        arrange(LOCUS, POP_ID)
    }
    # filter <- filter %>% select(-c(CHROM, LOCUS, POS))
  } # end maf haplotype approach
  
  if (maf.approach == "SNP") { # SNP approach
    if (maf.operator == "OR") {
      filter <- maf.data %>%
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
      filter <- maf.data %>%
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
  } # end maf snp approach
  
  # unused object
  maf.data <- NULL    
  
  # Update filters.parameters SNP ----------------------------------------------
  if (data.type == "vcf.file") {
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
    FILTERS = c("Minor Allele Frequency", rep(as.character(""), 4)),
    PARAMETERS = c("maf.approach", "maf.local.threshold", "maf.global.threshold", "maf.operator", "maf.pop.num.threshold"), 
    VALUES = c(maf.approach, paste(">=", maf.local.threshold), paste(">=", maf.global.threshold), maf.operator, paste(">=", maf.pop.num.threshold)), 
    BEFORE = c("", "", "", "", markers.before),
    AFTER = c("", "", "", "", markers.after),
    BLACKLIST = c("", "", "", "", markers.blacklist),
    UNITS = c("", "", "", "", "SNP/LOCUS"),
    COMMENTS = c("", "", "", "", "")
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
  message("Writing the whitelist of markers in your working directory\nwhitelist.markers.maf.tsv")
  
  if (data.type == "vcf.file") {
    whitelist.markers <- filter %>% 
      ungroup() %>%
      select(CHROM, LOCUS, POS) %>% 
      distinct(CHROM, LOCUS, POS)
  } else {
    whitelist.markers <- filter %>% 
      ungroup() %>%
      select(LOCUS = MARKERS) %>% 
      distinct(LOCUS)
  }
  write_tsv(whitelist.markers, "whitelist.markers.maf.tsv", append = FALSE, col_names = TRUE)
  
  
  # saving blacklist
  message("Writing the blacklist of markers in your working directory\nblacklist.markers.maf.tsv")
  if (data.type == "vcf.file") {
    blacklist.markers <- input %>% 
      ungroup() %>%
      select(CHROM, LOCUS, POS) %>% 
      distinct(CHROM, LOCUS, POS) %>% 
      anti_join(whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  } else {
    blacklist.markers <- input %>% 
      ungroup() %>%
      select(LOCUS = MARKERS) %>% 
      distinct(LOCUS) %>% 
      anti_join(whitelist.markers, by = "LOCUS")
  }
  write_tsv(blacklist.markers, "blacklist.markers.maf.tsv", append = FALSE, col_names = TRUE)
  
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stri_paste("maf.approach: ", maf.approach))
  message(stri_paste("maf.thresholds: ", "local = ", maf.local.threshold, ", global = ", maf.global.threshold))
  message(stri_paste("maf.operator: ", maf.operator))
  message(stri_paste("maf.pop.num.threshold: ", maf.pop.num.threshold))
  if (data.type != "vcf.file" & data.type != "haplo.file") {
    message(stri_paste("The number of markers removed by the MAF filter: ", n_distinct(input$MARKERS)-n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the MAF filter")
    message(stri_paste("SNP: ", as.integer(n_distinct(input$MARKERS)), " -> ", as.integer(n_distinct(filter$MARKERS))))
  } else if (data.type == "vcf.file") {
    message(stri_paste("The number of markers removed by the MAF filter:\nSNP: ", n_distinct(input$POS) - n_distinct(filter$POS), "\nLOCUS: ", n_distinct(input$LOCUS) - n_distinct(filter$LOCUS)))
    message("The number of markers before -> after the MAF filter")
    message(stri_paste("SNP: ", as.integer(n_distinct(input$POS)), " -> ", as.integer(n_distinct(filter$POS))))
    message(stri_paste("LOCUS: ", as.integer(n_distinct(input$LOCUS)), " -> ", as.integer(n_distinct(filter$LOCUS))))
  } else {# for haplotype file
    message(stri_paste("The number of markers/locus removed by the MAF filter: ", n_distinct(input$MARKERS)-n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the MAF filter")
    message(stri_paste("MARKERS/LOCUS: ", as.integer(n_distinct(input$MARKERS)), " -> ", as.integer(n_distinct(filter$MARKERS))))
  }
  cat("#######################################################################\n")
  res <- list()
  res$tidy.filtered.maf <- filter
  res$whitelist.markers <- whitelist.markers
  res$blacklist.markers <- blacklist.markers
  res$blacklist.markers <- blacklist.markers
  res$maf.data.thresholds <- maf.data.thresholds
  res$strata <- strata
  res$filters.parameters <- filters.parameters
  if (interactive.filter){
    res$violinplot.maf.global <- violinplot.maf.global
    res$violinplot.maf.local <- violinplot.maf.local
    res$plot.distribution.maf.local <- plot.distribution.maf.local
    res$maf.global.summary <- maf.global.summary
    res$maf.helper.table <- maf.helper.table
  }
  return(res)
}
