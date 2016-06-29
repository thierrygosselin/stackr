#' @name filter_genotype_likelihood
#' @title Genotype likelihood interactive filter
#' @description Filter genotypes and markers based on genotype likelihood using
#' a tidy data set (long format) and any VCF file both with GL information.
#' The function uses 
#' \code{\link[stackr]{tidy_genomic_data}} and 
#' \code{\link[stackr]{read_long_tidy_wide}} to load the file.

#' @param data vcf file or a data frame in long/tidy format.
#' The function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}} and 
#' \code{\link[stackr]{tidy_genomic_data}}.

#' @param interactive.filter (optional, logical) Do you want the filtering session to 
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is 
#' asked to see figures of distribution before making decisions for filtering 
#' with the genotype likelihood.

#' @param strata (optional for data frame and PLINK files, 
#' required for VCF and haplotypes files) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' Default: \code{strata = NULL}.

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



#' @param approach Character. By \code{approach = "SNP"} or 
#' by \code{approach = "haplotype"}. 
#' The function will consider the SNP or haplotype GL statistics to filter the marker. 
#' Default: \code{approach = "haplotype"}.

#' @param gl.ind.threshold (Integer, optional if interactive session) 
#' Threshold number of individual's genotype likelihood.
#' @param gl.mean.threshold (Integer, optional if interactive session) 
#' Threshold number of mean genotype likelihood. 
#' @param gl.min.threshold (Integer, optional if interactive session) 
#' Threshold number of min genotype likelihood. 
#' @param gl.diff.threshold (Integer, optional if interactive session)
#' Threshold number of diff genotype likelihood,
#' the difference between the max and min GL over a loci/read/haplotype.

#' @param pop.threshold (optional if interactive session)
#' A threshold number: proportion, percentage
#' or fixed number e.g. 0.50, 50 or 5.
#' 
#' @param percent (logical, optional if interactive session). 
#' Is the pop.threshold argument a percentage? TRUE/FALSE.
#' This argument is necessary to distinguish percentage from integer 
#' for population threshold, (e.g. 5 percent or 5 populations).

#' @param filename (optional) Name of the filtered data set, 
#' written to the working directory.

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


#' @details There is 4 steps in the interactive version:  
#' 
#' 1. gl_individuals_populations: the user is asked to inspect 
#' the genotype likelihood at the individuals and populations levels"), 
#' 
#' 2. blacklist_genotypes: the option is given to blacklist individual genotypes
#' based on low quality GL, 
#' 
#' 3. gl_markers: the user is asked to inspecting the genotype likelihood 
#' at the marker level and 
#' 
#' 4. blacklist_markers: the option is given to blacklist markers genotypes 
#' based on low quality GL found at the marker level.
#' 
#' Using the haplotype approach: The summary statistics are averaged
#' and nested: SNP -> individuals -> population -> loci. e.g. the mean GL is the average
#' genotype likelihood for all individuals of pop x for loci x.
#' The gl.diff.threshold is the difference between the max and min GL found for 
#' a locus. This argument is given to help your spot big difference in GL for 
#' a loci. Markers with small differences have more stability in the estimate.
#' @return With \code{interactive.filter = FALSE}, the function returns a 
#' filtered tidy vcf data frame inside the global environment. 
#' In the working directory 5 files are written:
#' 1. filters_parameters.tsv is updated or created.
#' The filter parameters, values and returned results 
#' (markers numbers or blacklist genotypes) are inside that file. 
#' 2. blacklist.genotypes.gl.tsv, 3. blacklist.markers.gl.tsv, 
#' 4. whitelist.markers.gl.tsv and 
#' 5. vcf.tidy.paralogs.id.gl.tsv (the same filtered tidy vcf data frame).
#' 
#' With \code{interactive.filter = TRUE}, a list with 20 objects is created. 
#' The information range from summary tables to plots. The objects names are 
#' found by using \code{names(list.name)}. The object can be isolated in separate
#' object outside the list by following the example below.
#' @rdname filter_genotype_likelihood
#' @export
#' @import dplyr
#' @import readr
#' @importFrom data.table fread

#' @examples
#' \dontrun{
#' If I have 10 populations, the 3 examples below will give the same output results:
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(
#' data = "beluga.vcf", 
#' interactive.filter = TRUE,
#' approach = "haplotype", 
#' gl.ind.threshold = 0,
#' gl.mean.threshold = 10, 
#' gl.min.threshold = 0, 
#' gl.diff.threshold = 200, 
#' pop.threshold = 50, 
#' percent = TRUE)
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(
#' data = "beluga.vcf", 
#' interactive.filter = FALSE,
#' approach = "haplotype",
#' gl.ind.threshold = 0,
#' gl.mean.threshold = 10, 
#' gl.min.threshold = 0, 
#' gl.diff.threshold = 200, 
#' pop.threshold = 0.5, 
#' percent = FALSE)
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(
#' data = "beluga.vcf", 
#' approach = "haplotype", 
#' gl.ind.threshold = 0,
#' gl.mean.threshold = 10, 
#' gl.min.threshold = 0, 
#' gl.diff.threshold = 100, 
#' pop.threshold = 5, 
#' percent = FALSE)
#' 
#' If interactive.filter = TRUE, a list is created and to view the filtered tidy vcf:
#' tidy.data <- beluga.vcf.tidy.gl$tidy.vcf.filtered.gl
#' 
#' Inside the same list, to isolate the blacklist.genotypes:
#' bg <- beluga.vcf.tidy.gl$blacklist.genotypes
#' }

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1"){
  utils::globalVariables(
    c('GL', 'GL_MAX', 'GL_MIN', 'MIN_REF', 'MIN_ALT', 'READ_DEPTH_MAX', 
      'GL_MEAN', 'GL_DIFF', 'POP_ID', 'BLACKLIST', '..scaled..', 
      'GENOTYPE_LIKELIHOOD_GROUP', 'VALUE', 'ALLELE_COVERAGE_RATIO', 
      "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")
  )
}

filter_genotype_likelihood <- function (
  data,
  interactive.filter = TRUE,
  strata = NULL, 
  pop.levels = NULL, 
  pop.labels = NULL, 
  pop.select = NULL,
  approach = "haplotype",
  gl.ind.threshold = NULL,
  gl.mean.threshold = NULL,
  gl.min.threshold = NULL,
  gl.diff.threshold = NULL,
  pop.threshold = NULL,
  percent = NULL,
  blacklist.id = NULL, 
  blacklist.genotype = NULL, 
  whitelist.markers = NULL, 
  monomorphic.out = TRUE, 
  max.marker = NULL,
  snp.ld = NULL, 
  common.markers = NULL,
  filename = NULL
) {
  
  cat("#######################################################################\n")
  cat("################# stackr: filter_genotype_likelihood ##################\n")
  cat("#######################################################################\n")
  
  
  # manage missing arguments -----------------------------------------------------  
  if (missing(data)) stop("missing input file")
  
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter == TRUE){
    message("With the interactive mode turn on, we will go through 4 steps to visualize and filter the data based on genotype likelihood")
    message("Step 1. gl_individuals_populations: Inspecting the genotype likelihood at the individuals and populations levels")
    message("Step 2. blacklist_genotypes: blacklisting (erasing) or not individual genotypes based on low quality GL")
    message("Step 3. gl_markers: Inspecting the genotype likelihood at the marker level")
    message("Step 4. blacklist_markers: Blacklisting (erasing) or not markers genotypes based on low quality GL found at the marker level")
    
    # Folder -------------------------------------------------------------------
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    path.folder <- stri_join(getwd(),"/", "filter_gl_", file.date, sep = "")
    dir.create(file.path(path.folder))
    
    message(stri_join("Folder created: ", path.folder))
    file.date <- NULL #unused object
  }
  # Filter parameter file ------------------------------------------------------
  message("Parameters used in this run will be store in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if(length(filters.parameters) == 0) {
    filters.parameters <- data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("Created a file to store filters parameters: filters_parameters.tsv")
  } else {
    message("Using the filters parameters file found in the directory")
  }
  
  # File type detection ********************************************************
  if(is.genind(data)) stop("vcf file with GL information is necessary with this function")
  if (stri_detect_fixed(str = data, pattern = ".tped")) stop("vcf file with GL information is necessary with this function")
  if (stri_detect_fixed(str = data.type, pattern = "Catalog")) stop("vcf file with GL information is necessary with this function")
  if (stri_detect_fixed(str = data, pattern = ".gen")) stop("vcf file with GL information is necessary with this function")
  
  data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
  if (identical(data.type, "##fileformat=VCF") | stri_detect_fixed(str = data, pattern = ".vcf")) {
    data.type <- "vcf.file"
    message("File type: VCF")
  }
  if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stri_detect_fixed(str = data.type, pattern = "MARKERS")) {
    data.type <- "df.file"
    message("File type: data frame of genotypes")
  }
  
  
  # import data ----------------------------------------------------------------
  if (is.vector(data)){
    message("Using input file in your directory")
    
    input <- stackr::tidy_genomic_data(
      data = data, 
      vcf.metadata = TRUE,
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
  
  # pop number
  pop.number <- n_distinct(data$POP_ID) # get the number of population
  
  # function--------------------------------------------------------------------
  plot_violinplot_genotype_likelihood_individuals<- function(data) {
    plot <- ggplot(data, aes(x = factor(POP_ID), y = GL, na.rm = TRUE))+
      geom_violin(trim = TRUE)+
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
      labs(x = "Sampling sites")+
      labs(y = "Genotype likelihood of individuals")+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
      )
    return(plot)
  }
  
  
  ## Step 1: gl_individuals_populations ----------------------------------------
  if (interactive.filter == TRUE){
    message("Step 1. gl_individuals_populations: Inspecting the genotype likelihood at the individuals and populations levels")
    path.folder.step1 <- stri_paste(path.folder, "/01_gl_individuals_populations")
    dir.create(path.folder.step1)
    summary <- summary_genotype_likelihood(tidy.vcf = data, pop.levels = pop.levels, approach = approach, folder = path.folder.step1)
  }
  
  
  # plot_1: Violin plot GL individuals and pop
  if (interactive.filter == TRUE){
    message("Show the violin plot of individuals genotype likelihood (y/n)): ")
    violinplot <- as.character(readLines(n = 1))
    if (violinplot == "y") {
      message("Generating violin plot may take some time...")
      # plot
      genotype.likelihood.violin.plot.individuals <- plot_violinplot_genotype_likelihood_individuals(data = data)
      print(genotype.likelihood.violin.plot.individuals)
      # save
      ggsave(stri_paste(path.folder.step1, "/genotype.likelihood.violin.plot.individuals.pdf"), width = pop.number, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stri_paste(path.folder.step1, "/genotype.likelihood.violin.plot.individuals.png"), width = pop.number, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.individuals) were saved in this directory: ", path.folder.step1))
    }
  }
  ## Step 2: blacklist_genotypes -------------------------------------------------
  # Erasing genotypes below threshold
  if (interactive.filter == TRUE){
    message("Step 2. blacklist_genotypes: blacklisting (erasing) or not individual genotypes based on low quality GL")
    message("Do you want to erase genotypes below a defined threshold (y/n) ?" )
    erase.genotype <- readLines(n = 1)
    if (erase.genotype == "y") {
      message("Enter the genotype likelihood threshold values, below this threshold (GL < threshold), genotypes are erased:")
      gl.ind.threshold <- as.integer(readLines(n = 1))
      path.folder.step2 <- stri_paste(path.folder, "/02_blacklist_genotypes")
      dir.create(path.folder.step2)
    }
  }
  
  if (is.null(gl.ind.threshold)) {
    data.filter.gl.individuals <- data
    blacklist.genotypes <- NULL
  } else {
    
    data.filter.gl.individuals <- data %>% 
      mutate(
        BLACKLIST = ifelse(GL < gl.ind.threshold, "erase", "keep"),
        BLACKLIST = stri_replace_na(str = BLACKLIST, replacement = "keep")
      )
    
    blacklist.genotypes <- data.filter.gl.individuals %>% 
      filter(BLACKLIST == "erase") %>% 
      select(CHROM, LOCUS, POS, POP_ID, INDIVIDUALS)
    
    data.filter.gl.individuals <- suppressWarnings(
      data.filter.gl.individuals %>% 
        mutate(
          GT = ifelse(BLACKLIST == "erase", "./.", GT),
          GL = as.numeric(ifelse(BLACKLIST == "erase", "NA", GL)),
          READ_DEPTH = as.numeric(ifelse(BLACKLIST == "erase", "NA", READ_DEPTH)),
          ALLELE_REF_DEPTH = as.numeric(ifelse(BLACKLIST == "erase", "NA", ALLELE_REF_DEPTH)),
          ALLELE_ALT_DEPTH = as.numeric(ifelse(BLACKLIST == "erase", "NA", ALLELE_ALT_DEPTH)),
          ALLELE_COVERAGE_RATIO = as.numeric(ifelse(BLACKLIST == "erase", "NA", ALLELE_COVERAGE_RATIO))
        ) %>% 
        select(-BLACKLIST)
    )
    
    # interesting stats.
    erased.genotype.number <- length(blacklist.genotypes$INDIVIDUALS)
    total.genotype.number <- length(data$GT[data$GT != "./."])
    percentage <- paste(round(((erased.genotype.number/total.genotype.number)*100), 6), "%", sep = " ")
    cat("######################## ERASING GENOTYPES ############################\n")
    message(stri_paste("Total number of genotypes: ", total.genotype.number))
    message(stri_paste("Blacklisted genotypes bsed on GL: ", erased.genotype.number))
    message(stri_paste("Percentage erased: ", percentage))
    cat("#######################################################################\n")
    
    write_tsv(x = blacklist.genotypes, path = "blacklist.genotypes.gl.tsv", col_names = TRUE)
    message("Writing the blacklist genotypes based on GL information in your working directory\nblacklist.genotypes.gl.tsv")
    
    # Update filters.parameters
    message("Updating the file storing the filters parameters: filters_parameters.tsv")
    filters.parameters <- data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("gl.ind.threshold"), VALUES = as.integer(gl.ind.threshold), BEFORE = as.character(total.genotype.number), AFTER = as.character(total.genotype.number-erased.genotype.number), BLACKLIST = erased.genotype.number, UNITS = "genotypes", COMMENTS = as.character("NA"))
    write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  }
  
  if (interactive.filter == TRUE){
    message("Create summary and show the violin plot of updated individuals genotype likelihood,\ni.e. after erasing lower quality genotypes (y/n)): ")
    summary.blacklist.genotypes <- as.character(readLines(n = 1))
    if (summary.blacklist.genotypes == "y") {
      # updated summary
      summary.blacklist.genotypes <- summary_genotype_likelihood(tidy.vcf = data.filter.gl.individuals, pop.levels = pop.levels, approach = approach, folder = path.folder.step2)
      
      # plot_2: updated plot after erasing genotypes
      message("Generating the updated violin plot may take some time...")
      genotype.likelihood.violin.plot.individuals.updated <- plot_violinplot_genotype_likelihood_individuals(data = data.filter.gl.individuals)
      print(genotype.likelihood.violin.plot.individuals.updated)
      # save
      ggsave(stri_paste(path.folder.step2, "/genotype.likelihood.violin.plot.individuals.updated.pdf"), width = pop.number, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stri_paste(path.folder.step2, "/genotype.likelihood.violin.plot.individuals.updated.png"), width = pop.number, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.individuals.updated) were saved in this directory: ", path.folder.step2))
    }
    # plot_3: before and after figure
    message("Create a combined violin plot of the individuals genotype likelihood with\nbefore and after erasing lower quality genotypes (y/n)): ")
    combined.plot <- as.character(readLines(n = 1))
    if (combined.plot == "y") {
      data.combined <- bind_rows(
        data.select <- data %>% 
          select(POP_ID, INDIVIDUALS, GL) %>% 
          mutate(GROUP = rep("before", n())),
        data.filter.gl.individuals.select <- data.filter.gl.individuals %>% 
          select(POP_ID, INDIVIDUALS, GL) %>% 
          mutate(GROUP = rep("after", n()))
      ) %>% 
        mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
      
      # combined plot
      message("Generating combined violin plot, may take some time...")
      genotype.likelihood.violin.plot.individuals.combined <- plot_violinplot_genotype_likelihood_individuals(data = data.combined) 
      genotype.likelihood.violin.plot.individuals.combined$facet <- facet_grid(~GROUP)
      print(genotype.likelihood.violin.plot.individuals.combined)
      # save
      ggsave(stri_paste(path.folder.step2, "/genotype.likelihood.violin.plot.individuals.combined.pdf"), width = pop.number*2, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stri_paste(path.folder.step2, "/genotype.likelihood.violin.plot.individuals.combined.png"), width = pop.number*2, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.individuals.combined) were saved in this directory: ", path.folder.step2))
    }
  }
  
  ## Step 3: gl_markers -------------------------------------------------------
  # We inspect at the locus level, by pop and generate the figures
  if (interactive.filter == TRUE){
    message("Step 3. gl_markers: Inspecting the genotype likelihood at the marker level")
    path.folder.step3 <- stri_paste(path.folder, "/03_gl_markers")
    dir.create(path.folder.step3)
    message("Show the density distribution plot (y/n)): ")
    density.distribution <- as.character(readLines(n = 1))
    
    # plot_4:Density distribution of genotype likelihood summary of loci
    if (density.distribution == "y") {
      genotype.likelihood.density.distribution.figure <- plot_density_distribution_genotype_likelihood(data = summary.blacklist.genotypes$gl.summary.marker.pop, aes.colour = aes(y = ..scaled.., color = GENOTYPE_LIKELIHOOD_GROUP), adjust.bin = 1) + facet_grid(POP_ID~GENOTYPE_LIKELIHOOD_GROUP, scales = "free")
      print(genotype.likelihood.density.distribution.figure)
      ggsave(stri_paste(path.folder.step3, "/genotype.likelihood.density.distribution.figure.pdf"), width = 15, height = pop.number*1.5, dpi = 600, units = "cm", useDingbats = FALSE)
      ggsave(stri_paste(path.folder.step3, "/genotype.likelihood.density.distribution.figure.png"), width = 15, height = pop.number*1.5, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the plot (genotype.likelihood.density.distribution.figure) were saved in this directory: ", path.folder.step3))
    }
    
    # plot_5: violin plot
    message("Show the violin plot (y/n)): ")
    violinplot <- as.character(readLines(n = 1))
    if (violinplot == "y") {
      genotype.likelihood.violin.plot.figure <- plot_boxplot_genotype_likelihood(data = summary$gl.summary.marker.pop) + facet_wrap(~GENOTYPE_LIKELIHOOD_GROUP, nrow = 1, ncol = 5, scales = "free")
      print(genotype.likelihood.violin.plot.figure)
      ggsave(stri_paste(path.folder.step3, "/genotype.likelihood.violin.plot.figure.pdf"), width = pop.number*2, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
      ggsave(stri_paste(path.folder.step3, "/genotype.likelihood.violin.plot.figure.png"), width = pop.number*2, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.figure) were saved in this directory: ", path.folder.step3))
      message("Look for inconsistencies, patterns and trends between your populations")
    }
  }
  
  if (interactive.filter == FALSE){
    if (approach == "haplotype"){
      message("Approach selected for GL statistics: haplotype")
      data.sum <- data.filter.gl.individuals %>%
        group_by(LOCUS, POP_ID) %>% # at the population level
        summarise(
          GL_MEAN = mean(GL, na.rm = TRUE),
          GL_MEDIAN = stats::median(GL, na.rm = TRUE),
          GL_MIN = min(GL, na.rm = TRUE),
          GL_MAX = max(GL, na.rm = TRUE),
          GL_DIFF = GL_MAX - GL_MIN
        ) %>% group_by(LOCUS, POP_ID)
    } else {
      message("Approach selected for GL statistics: SNP")
      data.sum <- data.filter.gl.individuals %>%
        group_by(LOCUS, POS, POP_ID) %>% # at the population level
        summarise(
          GL_MEAN = mean(GL, na.rm = TRUE),
          GL_MEDIAN = stats::median(GL, na.rm = TRUE),
          GL_MIN = min(GL, na.rm = TRUE),
          GL_MAX = max(GL, na.rm = TRUE),
          GL_DIFF = GL_MAX - GL_MIN
        ) %>% group_by(LOCUS, POS, POP_ID)
    }
  } else { # interactive
    if (approach == "haplotype"){
      data.sum <- summary.blacklist.genotypes$gl.summary.marker.pop %>%
        group_by(LOCUS, POP_ID) %>% 
        tidyr::spread(data = ., key = GENOTYPE_LIKELIHOOD_GROUP, value = VALUE)
    } else {
      data.sum <- summary.blacklist.genotypes$gl.summary.marker.pop %>%
        group_by(LOCUS, POS, POP_ID) %>% 
        tidyr::spread(data = ., key = GENOTYPE_LIKELIHOOD_GROUP, value = VALUE)
    }
  }
  
  ## Step 4: blacklist_markers ---------------------------------------------------
  # At this point interactive or not, we have the same input data
  if (interactive.filter == TRUE){
    message("Step 4. blacklist_markers: Blacklisting (erasing) or not markers genotypes based on low quality GL found at the marker level")
    message("Enter the population threshold (proportion, percentage or fixed) to keep the marker\ni.e. that whatever the locus GL statistics, it will need to pass in x (proportion, percentage or fixed) pop to be kept:")
    pop.threshold <- as.numeric(readLines(n = 1))
    path.folder.step4 <- stri_paste(path.folder, "/04_blacklist_markers")
    dir.create(path.folder.step4)
  }
  
  # percent ?
  if (interactive.filter == TRUE){
    message("Is the pop.threshold entered above a percentage (TRUE/FALSE) ?")
    percent <- as.character(readLines(n = 1))
  }
  
  if(stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message(stri_paste("Using a proportion threshold of , ", pop.threshold))
    threshold.id <- "proportion"
  } else if (stri_detect_fixed(percent, "T")) {
    multiplication.number <- 100/pop.number
    message(stri_paste("Using a percentage threshold of ", pop.threshold))
    threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message(stri_paste("Using a fixed threshold of ", pop.threshold))
    threshold.id <- "fixed"
  }
  
  # gl.mean.threshold
  if (interactive.filter == TRUE){  
    message(stri_paste("Enter the mean genotype likelihood threshold number.\ni.e. markers with gl.mean.threshold >= threshold value in ", pop.threshold, " (", threshold.id, ") ", "pop, will be kept"))
    message("To turn off this filter, enter: off")
    message("gl.mean.threshold:")
    gl.mean.threshold <- readLines(n = 1)
  }
  if (gl.mean.threshold == "off") gl.mean.threshold <- NULL
  
  if (!is.null(gl.mean.threshold)) {
    gl.mean.threshold <- as.integer(gl.mean.threshold)
    data.sum <- data.sum %>%
      filter(GL_MIN >= gl.mean.threshold)
    
    # Update filters.parameters
    filters.parameters <- data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("gl.mean.threshold"), VALUES = as.integer(gl.mean.threshold), BEFORE = as.character("NA"), AFTER = as.character("NA"), BLACKLIST = as.character("NA"), UNITS = "NA", COMMENTS = as.character("NA"))
    write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  }
  
  
  # gl.min.threshold
  if (interactive.filter == TRUE){
    message(stri_paste("Enter the minimum genotype likelihood threshold number.\ni.e. markers with gl.min.threshold >= threshold value in ", pop.threshold, " (", threshold.id, ") ", "pop, will be kept"))
    message("To turn off this filter, enter: off")
    message("gl.min.threshold:")
    gl.min.threshold <- readLines(n = 1)
  }
  
  if (gl.min.threshold == "off") gl.min.threshold <- NULL
  
  if (!is.null(gl.min.threshold)) {
    gl.min.threshold <- as.integer(gl.min.threshold)
    data.sum <- data.sum %>%
      filter(GL_MIN >= gl.min.threshold)
    # Update filters.parameters
    filters.parameters <- data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("gl.min.threshold"), VALUES = as.integer(gl.min.threshold), BEFORE = as.character("NA"), AFTER = as.character("NA"), BLACKLIST = as.character("NA"), UNITS = "NA", COMMENTS = as.character("NA"))
    write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  }
  
  # gl.diff.threshold
  if (interactive.filter == TRUE){
    message(stri_paste("Enter the difference between the max and min genotype likelihood threshold number.\ni.e. markers with gl.diff.threshold <= threshold value in ", pop.threshold, " (", threshold.id, ") ", "pop, will be kept"))
    message("To turn off this filter, enter: off")
    message("gl.diff.threshold:")
    gl.diff.threshold <- readLines(n = 1)
  }
  
  if (gl.diff.threshold == "off") gl.diff.threshold <- NULL
  
  if (!is.null(gl.diff.threshold)) {
    gl.diff.threshold <- as.integer(gl.diff.threshold)
    data.sum <- data.sum %>%
      filter(GL_DIFF <= gl.diff.threshold)
    # Update filters.parameters
    filters.parameters <- data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("gl.diff.threshold"), VALUES = as.integer(gl.diff.threshold), BEFORE = as.character("NA"), AFTER = as.character("NA"), BLACKLIST = as.character("NA"), UNITS = "NA", COMMENTS = as.character("NA"))
    write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  }
  
  # Filter the data set
  if (approach == "haplotype"){
    message("Approach selected: haplotype")
    filter <- data.sum %>%
      group_by(LOCUS) %>%
      tally() %>% # Globally accross loci
      filter((n * multiplication.number) >= as.numeric(pop.threshold)) %>%
      select(LOCUS) %>% 
      left_join(data.filter.gl.individuals, by = "LOCUS") %>%
      arrange(LOCUS, POS, POP_ID)
  } else {
    message("Approach selected: SNP")
    filter <- data.sum %>%
      group_by(LOCUS, POS) %>%
      tally() %>% # Globally accross loci
      filter((n * multiplication.number) >= as.numeric(pop.threshold)) %>%
      select(LOCUS, POS) %>%
      left_join(data.filter.gl.individuals, by = c("LOCUS", "POS")) %>%
      arrange(LOCUS, POS, POP_ID)
  }
  
  # Update filters.parameters SNP
  filters.parameters <- data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("pop.threshold"), VALUES = stri_paste(pop.threshold, " (", threshold.id, ")"), BEFORE = as.integer(n_distinct(data.filter.gl.individuals$POS)), AFTER = as.integer(n_distinct(filter$POS)), BLACKLIST = n_distinct(data.filter.gl.individuals$POS) - n_distinct(filter$POS), UNITS = "SNP", COMMENTS = as.character("NA"))
  write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  
  # Update filters.parameters LOCUS
  filters.parameters <- data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("pop.threshold"), VALUES = stri_paste(pop.threshold, " (", threshold.id, ")"), BEFORE = as.integer(n_distinct(data.filter.gl.individuals$LOCUS)), AFTER = as.integer(n_distinct(filter$LOCUS)), BLACKLIST = n_distinct(data.filter.gl.individuals$LOCUS) - n_distinct(filter$LOCUS), UNITS = "LOCUS", COMMENTS = as.character("NA"))
  write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  
  # saving tidy data 
  if (!is.null(filename)) {
    message("Writing the tidy vcf file filtered with GL information in your working directory...")
    write_tsv(filter, filename, append = FALSE, col_names = TRUE)
  }
  
  # saving whitelist
  whitelist.markers <- filter %>% 
    ungroup() %>%
    distinct(CHROM, LOCUS, POS)
  message("Writing the whitelist of markers based on GL information in your working directory\nwhitelist.markers.gl.tsv")
  write_tsv(whitelist.markers, "whitelist.markers.gl.tsv", append = FALSE, col_names = TRUE)
  
  # saving blacklist
  blacklist.markers <- data %>% 
    ungroup() %>%
    distinct(CHROM, LOCUS, POS) %>% 
    anti_join(whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  message("Writing the blacklist of markers based on GL information in your working directory\nblacklist.markers.gl.tsv")
  write_tsv(blacklist.markers, "blacklist.markers.gl.tsv", append = FALSE, col_names = TRUE)
  
  if (interactive.filter == TRUE){
    message("Summary of individuals genotype likelihood, i.e. after filtering the markers")
    summary.blacklist.markers <- summary_genotype_likelihood(tidy.vcf = filter, pop.levels = pop.levels, approach = approach, folder = path.folder.step4)
    
    
    # before and after figure
    data.combined <- bind_rows(
      summary$gl.summary.marker.pop %>% 
        mutate(GROUP = rep("before", n())),
      summary.blacklist.markers$gl.summary.marker.pop %>% 
        mutate(GROUP = rep("after", n()))
    ) %>% 
      mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
    
    message("Show the density distribution plot before/after filters (y/n)): ")
    density.distribution <- as.character(readLines(n = 1))
    # plot_6: Density distribution of genotype likelihood summary of loci.
    if (density.distribution == "y") {
      genotype.likelihood.density.distribution.figure.before.after.filters <- plot_density_distribution_genotype_likelihood(data = data.combined, aes.colour = aes(y = ..scaled.., color = GENOTYPE_LIKELIHOOD_GROUP), adjust.bin = 1) + facet_grid(POP_ID+GROUP~GENOTYPE_LIKELIHOOD_GROUP, scales = "free")
      print(genotype.likelihood.density.distribution.figure.before.after.filters)
      ggsave(stri_paste(path.folder.step4, "/genotype.likelihood.density.distribution.figure.before.after.filters.pdf"), width = 20, height = pop.number*3, dpi = 600, units = "cm", useDingbats = FALSE)
      ggsave(stri_paste(path.folder.step4, "/genotype.likelihood.density.distribution.figure.before.after.filters.png"), width = 20, height = pop.number*3, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the plot (genotype.likelihood.density.distribution.figure.before.after.filters) were saved in this directory: ", path.folder.step4))
    }
    
    # plot_7: violin plot 
    message("Show the violin plot before/after filters (y/n)): ")
    violin.plot <- as.character(readLines(n = 1))
    if (violin.plot == "y") {
      # combined plot
      message("Generating combined violin plot, may take some time...")
      genotype.likelihood.violin.plot.figure.before.after.filters <- plot_boxplot_genotype_likelihood(data = data.combined) + facet_grid(GROUP~GENOTYPE_LIKELIHOOD_GROUP, scales = "fixed")
      print(genotype.likelihood.violin.plot.figure.before.after.filters)
      ggsave(stri_paste(path.folder.step4, "/genotype.likelihood.violin.plot.figure.before.after.filters.pdf"), width = pop.number*2, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
      ggsave(stri_paste(path.folder.step4, "/genotype.likelihood.violin.plot.figure.before.after.filters.png"), width = pop.number*2, height = 10, dpi = 300, units = "cm")
      message(stri_paste("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.figure.before.after.filters) were saved in this directory: ", path.folder.step4))
    }
  }
  
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stri_paste("The number of markers removed by the GL filter:\nSNP: ", n_distinct(data.filter.gl.individuals$POS) - n_distinct(filter$POS), "\nLOCUS: ", n_distinct(data.filter.gl.individuals$LOCUS) - n_distinct(filter$LOCUS)))
  message("The number of markers before -> after the GL filter")
  message(stri_paste("SNP: ", as.integer(n_distinct(data.filter.gl.individuals$POS)), " -> ", as.integer(n_distinct(filter$POS))))
  message(stri_paste("LOCUS: ", as.integer(n_distinct(data.filter.gl.individuals$LOCUS)), " -> ", as.integer(n_distinct(filter$LOCUS))))
  cat("#######################################################################\n")
  
  if (interactive.filter == TRUE){
    res <- list()
    res$gl.summary.individuals <- summary$gl.individuals
    res$gl.summary.marker.pop <- summary$gl.summary.marker.pop
    res$gl.summary.pop <- summary$gl.summary.pop
    res$gl.summary.individuals.blacklist.genotypes <- summary.blacklist.genotypes$gl.individuals
    res$gl.summary.marker.pop.blacklist.genotypes <- summary.blacklist.genotypes$gl.summary.marker.pop
    res$gl.summary.pop.blacklist.genotypes <- summary.blacklist.genotypes$gl.summary.pop
    res$genotype.likelihood.violin.plot.individuals <- genotype.likelihood.violin.plot.individuals
    res$genotype.likelihood.violin.plot.individuals.updated <- genotype.likelihood.violin.plot.individuals.updated
    res$genotype.likelihood.violin.plot.individuals.combined <- genotype.likelihood.violin.plot.individuals.combined
    res$genotype.likelihood.density.distribution.figure <- genotype.likelihood.density.distribution.figure
    res$genotype.likelihood.violin.plot.figure <- genotype.likelihood.violin.plot.figure
    res$gl.summary.individuals.filter.markers <- summary.blacklist.markers$gl.individuals
    res$gl.summary.marker.pop.filter.markers <- summary.blacklist.markers$gl.summary.marker.pop
    res$gl.summary.pop.filter.markers <- summary.blacklist.markers$gl.summary.pop
    res$genotype.likelihood.density.distribution.figure.before.after.filters <- genotype.likelihood.density.distribution.figure.before.after.filters
    res$genotype.likelihood.violin.plot.figure.before.after.filters <- genotype.likelihood.violin.plot.figure.before.after.filters
    res$blacklist.genotypes <- blacklist.genotypes
    res$blacklist.markers <- blacklist.markers
    res$whitelist.markers <- whitelist.markers
    res$tidy.vcf.filtered.gl <- filter
  } else {
    res <- filter
  }
  
  return(res)
}




