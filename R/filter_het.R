# Observed heterozygosity 
#' @title Heterozygosity filter
#' @description Observed Heterozygosity based filtering. 
#' The filter arguments of \code{filter_het} allows you to test rapidly,
#' with the use of blacklists of markers and individuals and whitelists of markers,
#' if departure from realistic expectations of heterozygosity statistics 
#' are a problem in downstream analysis.

#' \enumerate{
#' \item Highlight outliers individual's heterozygosity for a quick
#' diagnostic of mixed samples. 
#' 
#' \item The observed heterozygosity in the dataset: an assembly artefact, 
#' a genotyping problem, a problem of population groupings 
#' or a reliable signal of biological polymorphism? 
#' Detect assembly artifact or genotyping problem (e.g. under/over-splitting loci)
#' by looking at marker's observed heterozygosity statistics by population 
#' or overall.
#' 
#' \item Use haplotype or snp level statistics. When the haplotype approach
#' is selected, consistencies of heterozygosity statistics along the read are 
#' highlighted.
#' 
#' \item Interactive approach help by visualizing the data before making
#' a decision on thresholds.
#' }

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data

#' @param interactive.filter (optional, logical) Do you want the filtering session to 
#' be interactive. With the default, the user is asked to see figures of 
#' distribution before making decisions for filtering with heterozygosity statistics.
#' Default: \code{interactive.filter == TRUE}.

#' @param ind.heterozygosity.threshold (double, optional)
#' Blacklist individuals based on observed heterozygosity of individuals 
#' (averaged across markers).
#' The value for the threshold is a proportion (0 to 1, where 1 will turn off 
#' the filter).
#' Individuals with mean heterozygosity higher (>) than the threshold
#' will be blacklisted.
#' Default: \code{ind.heterozygosity.threshold = 1}.

#' @param het.approach (Character string). First value, \code{"SNP"} or
#' \code{"haplotype"}. The haplotype approach considers the statistic
#' consistency on the read (locus/haplotype). The major difference: 
#' the haplotype approach results in blacklisting the entire locus/haplotype
#' with all the SNPs on the read.
#' With the SNP approach, SNPs are independently analyzed and blacklisted.
#' The second value, will use the statistics by population \code{"pop"} or 
#' will consider the data overall, as 1 large group \code{"overall"}.
#' Default: \code{het.approach = c("SNP", "overall")}.

#' @param het.threshold Number Biallelic markers usually max 0.5. But departure
#' from this value is common. The higher the proportion threshold, the more relaxed
#' the filter is. With default there is no filtering.
#' Default: \code{het.threshold = 1}.

#' @param het.dif.threshold Number (0 - 1). For \code{het.approach = "haplotype"} only.
#' You can set a threshold for the difference in het along your read. 
#' Set the number your willing to tolerate on the same read/haplotype. 
#' e.g. if you have 2 SNP on a read/haplotype and on as a het of 0.9 and the other
#' 0.1 and you set \code{het.dif.threshold = 0.3}, this markers will be blacklisted.
#' You should strive to have similar statistics along short read like RAD.
#' The higher the proportion threshold, the more relaxed the filter is.
#' With default, there is no filtering and all the range of differences are allowed.
#' Default: \code{het.dif.threshold = 1}.

#' @param outlier.pop.threshold (integer, optional) Useful to incorporate problematic 
#' populations dragging down polymorphism discovery, but still wanted for analysis.
#' Use this threshold to allow variance in the number of populations passing 
#' the thresholds described above.
#' e.g. with \code{outlier.pop.threshold = 2},
#' you tolerate a maximum of 2 populations failing the 
#' \code{het.threshold} and/or \code{het.dif.threshold}.
#' Manage outlier markers, individuals and populations 
#' downstream with blacklists and whitelists produced by the function.
#' Default: \code{outlier.pop.threshold = 1}. See details for more info.

#' @param helper.tables (logical) Output tables that show
#' the number of markers blacklisted or whitelisted based on a series of 
#' automatic thresholds to guide decisions. When \code{interactive.filter == TRUE},
#' helper tables are written to the directory.
#' Default: \code{helper.tables = FALSE}. 

#' @param filename Name of the tidy data set written to the working directory (optional).
#' e.g. "tasmanian.devil.tidy.het.tsv". 
#' Default: \code{filename = NULL}

#' @details 
#' \strong{Interactive version}
#' 
#' There are 4 steps in the interactive version to visualize and filter
#' the data based heterozygosity statistics:
#' 
#' Step 1. Individual's heterozygosity: outliers that might represent mixed samples
#' Step 2. Blacklist outliers based on a proportion threshold of mean heterozygosity
#' Step 3. Observed heterozygosity statistics per populations and overall
#' Step 4: Blacklist markers based on observed heterozygosity
#' 
#' 
#' 
#' \strong{outlier.pop.threshold}
#' 
#' If your a regular stackr user, you've seen the \code{pop.num.threshold}.
#' \code{outlier.pop.threshold}, is different and requires more thinking,
#' because the number of populations genotyped potentially vary across markers,
#' which makes the use of \code{pop.num.threshold} less optimal.
#' e.g. If only 4 populations out of 10 are genotyped, for a marker you want 
#' to keep, using \code{pop.num.threshold = 0.6} will lead to unwanted results
#' and inconsistensis.
#' It's easier to tolerate outliers with this new approach: 
#' \code{outlier.pop.threshold = 2}.

#' @return The function returns inside the global environment a list with
#' 14 objects, the objects names are found by using \code{names(your.list.name)}:
#' 
#' \enumerate{
#' \item filtered tidy data frame: \code{$tidy.filtered.het}
#' \item whitelist of markers:\code{$whitelist.markers}
#' \item the strata:\code{$strata}
#' \item the filters parameters used:\code{$filters.parameters}
#' \item the individual's heterozigosity:\code{$individual.heterozigosity}
#' \item the blacklisted individuals based on the individual's heterozigosity:\code{$blacklist.ind.het}
#' \item a list containing the helper tables:\code{$helper.table.het}
#' \item the boxplot of individual heterozygosity:\code{$individual.heterozygosity.boxplot}
#' \item the manhattan plot of individual heterozygosity:\code{$individual.heterozygosity.manhattan.plot}
#' \item the boxplot of observed heterozygosity averaged across markers and pop:\code{$markers.pop.heterozygosity.boxplot}
#' \item the density plot of observed heterozygosity averaged across markers and pop:\code{$markers.pop.heterozygosity.density.plot}
#' \item the manhattan plot of observed heterozygosity averaged across markers and pop:\code{$markers.pop.heterozygosity.manhattan.plot}
#' \item whitelist of markers:\code{$whitelist.markers}
#' \item whitelist of markers:\code{$whitelist.markers}
#' }
#' 
#' 
#' In the working directory, output is conditional to \code{interactive.filter} argument

#' @rdname filter_het
#' @export

#' @import ggplot2
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame
#' @importFrom tidyr complete gather unite spread

#' @seealso \link{plot_density_distribution_het}


filter_het <- function(
  data,
  interactive.filter = TRUE,
  ind.heterozygosity.threshold = 1,
  het.approach = c("haplotype", "overall"),
  het.threshold = 1,
  het.dif.threshold = 1,
  outlier.pop.threshold = 1,
  helper.tables = FALSE,
  filename = NULL,
  vcf.metadata = FALSE,
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
  cat("######################### stackr::filter_het ##########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  # manage missing arguments -----------------------------------------------------  
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  if (!is.null(pop.select)) {
    pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter) {
    message("Interactive mode: on, 4 steps of data visualization and filtering based on heterozygosity:\n")
    message("Step 1. Individual's heterozygosity: outliers that might represent mixed samples")
    message("Step 2. Blacklist outliers based on a proportion threshold of mean heterozygosity")
    message("Step 3. Observed heterozygosity statistics per populations and overall")
    message("Step 4: Blacklist markers based on observed heterozygosity\n\n")
  }
  # Folder -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
  file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
  file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
  
  path.folder <- stringi::stri_join(getwd(),"/", "filter_het_", file.date, sep = "")
  dir.create(file.path(path.folder))
  
  message(stringi::stri_join("Folder created: ", path.folder))
  file.date <- NULL #unused object
  
  # Filter parameter file ------------------------------------------------------
  message("\nParameters used in this run will be store in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if (length(filters.parameters) == 0) {
    filters.parameters <- tibble::data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("Created a file to store filters parameters: filters_parameters.tsv\n")
  } else {
    message("Using the filters parameters file found in the directory\n")
  }
  # File type detection----------------------------------------------------------
  data.type <- stackr::detect_genomic_format(data)
  
  if (data.type == "haplo.file") {
    message("With stacks haplotype file the approach is automatically set to: haplotype")
    # counter intuitive, but there is no snp or haplotype info in a stacks haplotypes file
    het.approach[1] <- "SNP"
  }
  
  # import data ----------------------------------------------------------------
  message("Importing data ...")
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
    dplyr::select(INDIVIDUALS, POP_ID) %>% 
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)
  
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  pop.number <- dplyr::n_distinct(input$POP_ID) # number of pop
  
  # double check the approach vs the file used
  if (!tibble::has_name(input, "CHROM") & het.approach[1] == "haplotype") {
    stop("The haplotype approach during HET filtering requires LOCUS and SNP 
information")
  }
  
  # highlight heterozygote
  het.summary <- input %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::mutate(
      HET = dplyr::if_else(
        stringi::stri_sub(GT, 1, 3) != stringi::stri_sub(GT, 4, 6), 1, 0
      )
    )
  
  # input <- NULL # no longer needed
  
  
  # Step 1. Highlight individual's heterozygosity  -----------------------------
  # Heterozygosity at the individual level before looking at the markers level per population
  # It's a good way to do outlier diagnostic ... mixed individuals
  
  # Create a new df with heterozygote info
  het.ind <- het.summary %>% 
    dplyr::group_by(INDIVIDUALS, POP_ID) %>% 
    dplyr::summarise(
      GENOTYPED = n(),
      HET_NUMBER = length(HET[HET == 1]),
      HET_PROP = HET_NUMBER / GENOTYPED
    )
  
  # manhattan.max.breaks <- max(het.ind$HET_PROP)
  # manhattan.max.breaks <- 0.26
  
  
  
  individual.heterozygosity.manhattan.plot <- ggplot(data = het.ind, aes(x = POP_ID, y = HET_PROP, colour = POP_ID)) + 
    geom_jitter() + 
    labs(y = "Individual's Mean Heterozygosity (proportion)") +
    labs(x = "Populations") +
    labs(colour = "Populations") +
    # scale_y_continuous(name = , breaks = )
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    )
  
  individual.heterozygosity.boxplot <- ggplot(data = het.ind, aes(x = POP_ID, y = HET_PROP, colour = POP_ID)) + 
    geom_boxplot() + 
    labs(y = "Individual's Mean Heterozygosity (proportion)") +
    labs(x = "Populations") +
    labs(colour = "Populations") +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    )
  
  if (interactive.filter) {
    message("\nStep 1. Individual's heterozygosity: outliers that might represent mixed samples\n")
    print(individual.heterozygosity.manhattan.plot)
    # save
    ggsave(stringi::stri_join(path.folder, "/individual.heterozygosity.manhattan.plot.pdf"), width = pop.number * 2, height = 10, dpi = 600, units = "cm", useDingbats = F)
    ggsave(stringi::stri_join(path.folder, "/individual.heterozygosity.manhattan.plot.png"), width = pop.number * 2, height = 10, dpi = 300, units = "cm")
    message(stringi::stri_join("2 versions (pdf and png) of the plot (individual.heterozygosity.manhattan.plot) were saved in this directory:\n", path.folder))
  }
  
  
  if (interactive.filter) {
    message("\n\nShow the box plot of individual's heterozygosity (y/n): ")
    boxplot <- as.character(readLines(n = 1))
    if (boxplot == "y") {
      print(individual.heterozygosity.boxplot)
      # save
      ggsave(stringi::stri_join(path.folder, "/individual.heterozygosity.boxplot.pdf"), width = pop.number * 2, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder, "/individual.heterozygosity.boxplot.png"), width = pop.number * 2, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (individual.heterozygosity.boxplot) were saved in this directory:\n", path.folder))
    }
  }
  
  ## Step 2: Blacklist outlier individuals -------------------------------------
  # Blacklist individuals based a threshold of mean heterozygosity
  if (interactive.filter) {
    ind.heterozygosity.threshold <- 2#to make sure we go through selection of threshold
    while (isTRUE(ind.heterozygosity.threshold > 1)) {
      message("\n\nStep 2. Blacklist outliers based on a proportion threshold of mean heterozygosity\n")
      message("If you want to blacklist individuals based on the mean heterozygosity,
use this ind.heterozygosity.threshold filter.
The value for the ind.heterozygosity.threshold is a proportion (0 to 1, where 1 will turn off the filter).
Individuals with mean heterozygosity higher (>) than the threshold will be blacklisted.
Enter the value (proportion, e.g. 0.34):")
      ind.heterozygosity.threshold <- as.numeric(readLines(n = 1))
    }
  }
  
  if (is.double(ind.heterozygosity.threshold)) {
    blacklist.ind.het  <- dplyr::ungroup(het.ind) %>%
      dplyr::filter(HET_PROP > ind.heterozygosity.threshold) %>%
      dplyr::distinct(INDIVIDUALS)
  } else {
    ind.heterozygosity.threshold <- NULL
  }
  
  # Remove the individuals from the dataset
  message(stringi::stri_join("Filter individual's heterozygosity: ", length(blacklist.ind.het$INDIVIDUALS), " individual(s) deleted"))
  if (length(blacklist.ind.het$INDIVIDUALS > 0)) {
    readr::write_tsv(
      x = blacklist.ind.het,
      path = stringi::stri_join(path.folder, "/blacklist.individuals.heterozygosity.tsv"),
      col_names = TRUE
    )
    message("Blacklist (blacklist.individuals.heterozygosity) written in the working directory")
    het.summary <- dplyr::anti_join(het.summary, blacklist.ind.het, by = "INDIVIDUALS")
  }
  
  # Step 3. Markers observed heterozygosity statistics per populations and overall----
  
  # decide approach
  # Haplotype or SNP approach ...
  if (interactive.filter) {
    message("\n\nStep 3. Markers observed heterozygosity statistics per populations and overall\n")
    
    het.approach[1] <- "not defined"
    while (isTRUE(het.approach[1] != "haplotype" & het.approach[1] != "SNP")) {# to make sure the answer is ok 
      message("het.approach argument: haplotype or SNP ?\n")
      message("Decide on the best approach to filter the data based on oserved heterozygosity\n")
      message("The haplotype approach considers the statistic consistency on the read (locus/haplotype).
The major difference: the haplotype approach results in blacklisting the entire
locus/haplotype with all the SNPs on the read. With the SNP approach, SNPs are 
independently analyzed and blacklisted.\n")
      message("Enter if you want the haplotype or SNP approach (haplotype/SNP):")
      het.approach[1] <- as.character(readLines(n = 1))
    }
  }
  # overall or by populations...
  # het.approach[2] <- "test" #test
  if (interactive.filter) {
    het.approach[2] <- "not defined"
    while (isTRUE(het.approach[2] != "overall" & het.approach[2] != "pop")) {# to make sure the answer is ok
      message("het.approach argument: overall or by pop ?\n")
      message("Decide on the best approach to filter the data based on oserved heterozygosity\n")
      message("The overall approach doesn't look at the observed heterozygosity by populations.
It's just 1 large population of sample.
If you're not sure about your population structure and/or population sampling,
use the overall approach.\n")
      message("Enter if you want the overall or pop approach (overall/pop):")
      het.approach[2] <- as.character(readLines(n = 1))
    }
  }
  
  if (tibble::has_name(het.summary, "LOCUS") & het.approach[1] == "haplotype") {
    # By pop
    het.summary.pop <- het.summary %>%
      dplyr::group_by(MARKERS, LOCUS, POP_ID) %>%
      dplyr::summarise(HET_O = as.numeric(length(HET[HET == 1]) / n())) %>% 
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::summarise(
        HET_MEAN = mean(HET_O),
        HET_MAX = max(HET_O),
        HET_MIN = min(HET_O),
        HET_DIF = HET_MAX - HET_MIN
      )
    
    # overall
    het.summary.overall <- het.summary %>%
      dplyr::group_by(MARKERS, LOCUS) %>%
      dplyr::summarise(HET_O = as.numeric(length(HET[HET == 1]) / n())) %>% 
      dplyr::group_by(LOCUS) %>%
      dplyr::summarise(
        HET_MEAN = mean(HET_O),
        HET_MAX = max(HET_O),
        HET_MIN = min(HET_O),
        HET_DIF = HET_MAX - HET_MIN
      ) %>% 
      dplyr::mutate(POP_ID = rep("OVERALL", n()))
    
    
    het.summary.tidy <- suppressWarnings(
      dplyr::bind_rows(het.summary.pop, het.summary.overall) %>% 
        tidyr::gather(
          data = ., 
          key = HET_GROUP, 
          value = VALUE, 
          -c(LOCUS, POP_ID)
        ) %>% 
        dplyr::rename(MARKERS = LOCUS) %>% 
        dplyr::mutate(
          HET_GROUP = factor(
            HET_GROUP, 
            levels = c("HET_MEAN", "HET_MIN", "HET_MAX", "HET_DIF"), 
            ordered = TRUE
          )
        )
    )
  } else {
    # By pop
    het.summary.pop <- het.summary %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(HET_MEAN = as.numeric(length(HET[HET == 1]) / n()))
    
    # Overall
    het.summary.overall <- het.summary %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(HET_MEAN = as.numeric(length(HET[HET == 1]) / n())) %>% 
      dplyr::mutate(POP_ID = rep("OVERALL", n()))
    
    het.summary.tidy <- suppressWarnings(
      dplyr::bind_rows(het.summary.pop, het.summary.overall) %>% 
        tidyr::gather(
          data = ., 
          key = HET_GROUP, 
          value = VALUE, 
          -c(MARKERS, POP_ID)
        )
    )
  }
  # Tidy use for figures
  
  # # unused object
  het.summary <- het.summary.pop 
  het.summary.pop <- NULL 
  
  
  if (!is.null(pop.levels)) {
    het.summary.tidy <- het.summary.tidy %>%
      dplyr::mutate(POP_ID = factor(POP_ID, levels = c(pop.labels, "OVERALL"), ordered = TRUE)) %>% 
      dplyr::arrange(POP_ID)
  }
  
  
  
  # Visualization --------------------------------------------------------------
  markers.pop.heterozygosity.density.plot <- ggplot(het.summary.tidy, aes(x = VALUE, na.rm = FALSE)) +
    geom_line(aes(y = ..scaled..), stat = "density", adjust = 0.3) +
    labs(x = "Observed Heterozygosity") +
    labs(y = "Density of SNP (scaled)") +
    expand_limits(y = 0) +
    theme(
      axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"), 
      axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 12, family = "Helvetica", face = "bold"),
      strip.text.y = element_text(angle = 0, size = 12, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
    ) +
    facet_grid(POP_ID ~ HET_GROUP)
  
  markers.pop.heterozygosity.manhattan.plot <- ggplot(data = het.summary.tidy, aes(x = POP_ID, y = VALUE, colour = POP_ID)) + 
    geom_jitter() + 
    labs(y = "Observed Heterozygosity") +
    labs(x = "Populations") +
    labs(colour = "Populations") +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_text(size = 10, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    ) +
    facet_grid(~ HET_GROUP)
  
  
  markers.pop.heterozygosity.boxplot <- ggplot(data = het.summary.tidy, aes(x = POP_ID, y = VALUE, colour = POP_ID)) + 
    geom_boxplot() + 
    labs(y = "Observed Heterozygosity") +
    labs(x = "Populations") +
    labs(colour = "Populations") +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_text(size = 10, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    ) +
    facet_grid(~ HET_GROUP, scales = "free_y")
  
  het.summary.tidy <- NULL #unused object
  
  # helper tables ----------------------------------------------------------------
  
  if (helper.tables) {
    message("Generating helper table(s)...")
    if (tibble::has_name(het.summary, "LOCUS") &  het.approach[1] == "haplotype") {# by Haplotype
      n.markers <- dplyr::n_distinct(het.summary$LOCUS)
      
      # by pop
      het.helper <- dplyr::select(het.summary, LOCUS, POP_ID, HET_MAX, HET_DIF) %>% 
        dplyr::filter(POP_ID != "OVERALL") %>% 
        dplyr::mutate(
          MAX_0.4 = dplyr::if_else(HET_MAX <= 0.4, "whitelist", "blacklist"),
          MAX_0.5 = dplyr::if_else(HET_MAX <= 0.5, "whitelist", "blacklist"),
          MAX_0.6 = dplyr::if_else(HET_MAX <= 0.6, "whitelist", "blacklist"),
          MAX_0.7 = dplyr::if_else(HET_MAX <= 0.7, "whitelist", "blacklist"),
          MAX_0.8 = dplyr::if_else(HET_MAX <= 0.8, "whitelist", "blacklist"),
          MAX_0.9 = dplyr::if_else(HET_MAX <= 0.9, "whitelist", "blacklist")
        ) %>% 
        tidyr::gather(
          data = .,
          key = MAX_THRESHOLD,
          value = MAX_OUTLIERS,
          MAX_0.4:MAX_0.9
        ) %>%
        dplyr::mutate(
          DIF_0.1 = dplyr::if_else(HET_DIF <= 0.1, "whitelist", "blacklist"),
          DIF_0.2 = dplyr::if_else(HET_DIF <= 0.2, "whitelist", "blacklist"),
          DIF_0.3 = dplyr::if_else(HET_DIF <= 0.3, "whitelist", "blacklist"),
          DIF_0.4 = dplyr::if_else(HET_DIF <= 0.4, "whitelist", "blacklist"),
          DIF_0.5 = dplyr::if_else(HET_DIF <= 0.5, "whitelist", "blacklist"),
          DIF_0.6 = dplyr::if_else(HET_DIF <= 0.6, "whitelist", "blacklist"),
          DIF_0.7 = dplyr::if_else(HET_DIF <= 0.7, "whitelist", "blacklist"),
          DIF_0.8 = dplyr::if_else(HET_DIF <= 0.8, "whitelist", "blacklist"),
          DIF_0.9 = dplyr::if_else(HET_DIF <= 0.9, "whitelist", "blacklist")
        ) %>%
        tidyr::gather(
          data = .,
          key = DIF_THRESHOLD,
          value = DIF_OUTLIERS,
          DIF_0.1:DIF_0.9
        ) %>%
        tidyr::unite(
          data = .,
          col = MAX_DIF_THRESHOLD, MAX_THRESHOLD, DIF_THRESHOLD,
          sep = "_",
          remove = FALSE
        ) %>% 
        dplyr::mutate(
          MAX_DIF_OUTLIERS = dplyr::if_else(
            MAX_OUTLIERS == "whitelist" & DIF_OUTLIERS == "whitelist",
            "whitelist",
            "blacklist")
        ) %>% 
        dplyr::group_by(LOCUS, MAX_DIF_THRESHOLD, MAX_THRESHOLD, DIF_THRESHOLD) %>% 
        dplyr::summarise(
          MAX_DIF_OUTLIERS = length(MAX_DIF_OUTLIERS[MAX_DIF_OUTLIERS == "whitelist"]),
          MAX_OUTLIERS = length(MAX_OUTLIERS[MAX_OUTLIERS == "whitelist"]),
          DIF_OUTLIERS = length(DIF_OUTLIERS[DIF_OUTLIERS == "whitelist"])
        ) %>% 
        dplyr::ungroup(.) %>%
        dplyr::mutate(
          MAX_DIF_OUTLIERS = pop.number - MAX_DIF_OUTLIERS,
          MAX_OUTLIERS = pop.number - MAX_OUTLIERS,
          DIF_OUTLIERS = pop.number - DIF_OUTLIERS
        )
      
      max.threshold <- dplyr::distinct(
        het.helper, LOCUS, MAX_THRESHOLD, MAX_OUTLIERS
      ) %>%
        dplyr::group_by(MAX_THRESHOLD, MAX_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          MAX_OUTLIERS,
          nesting(MAX_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(MAX_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>% 
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>% 
        dplyr::group_by(MAX_THRESHOLD, MAX_OUTLIERS) %>% 
        tidyr::spread(data = ., key = MAX_THRESHOLD, value = MARKERS, fill = 0)
      
      dif.threshold <- dplyr::distinct(
        het.helper, LOCUS, DIF_THRESHOLD, DIF_OUTLIERS
      ) %>%
        dplyr::group_by(DIF_THRESHOLD, DIF_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          DIF_OUTLIERS,
          nesting(DIF_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(DIF_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>% 
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>% 
        dplyr::group_by(DIF_THRESHOLD, DIF_OUTLIERS) %>% 
        tidyr::spread(data = ., key = DIF_THRESHOLD, value = MARKERS, fill = 0)
      
      
      max.dif.threshold.combined <- dplyr::distinct(
        het.helper, LOCUS, MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS
      ) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          MAX_DIF_OUTLIERS,
          nesting(MAX_DIF_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>% 
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>% 
        dplyr::group_by(MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS) %>% 
        tidyr::spread(data = ., key = MAX_DIF_THRESHOLD, value = MARKERS, fill = 0)
      
      # overall
      het.helper.overall <- dplyr::select(het.summary.overall, LOCUS, HET_MAX, HET_DIF) %>% 
        dplyr::mutate(
          MAX_0.4 = dplyr::if_else(HET_MAX <= 0.4, "whitelist", "blacklist"),
          MAX_0.5 = dplyr::if_else(HET_MAX <= 0.5, "whitelist", "blacklist"),
          MAX_0.6 = dplyr::if_else(HET_MAX <= 0.6, "whitelist", "blacklist"),
          MAX_0.7 = dplyr::if_else(HET_MAX <= 0.7, "whitelist", "blacklist"),
          MAX_0.8 = dplyr::if_else(HET_MAX <= 0.8, "whitelist", "blacklist"),
          MAX_0.9 = dplyr::if_else(HET_MAX <= 0.9, "whitelist", "blacklist")
        ) %>% 
        tidyr::gather(
          data = .,
          key = MAX_THRESHOLD,
          value = MAX_OUTLIERS,
          MAX_0.4:MAX_0.9
        ) %>%
        dplyr::mutate(
          DIF_0.1 = dplyr::if_else(HET_DIF <= 0.1, "whitelist", "blacklist"),
          DIF_0.2 = dplyr::if_else(HET_DIF <= 0.2, "whitelist", "blacklist"),
          DIF_0.3 = dplyr::if_else(HET_DIF <= 0.3, "whitelist", "blacklist"),
          DIF_0.4 = dplyr::if_else(HET_DIF <= 0.4, "whitelist", "blacklist"),
          DIF_0.5 = dplyr::if_else(HET_DIF <= 0.5, "whitelist", "blacklist"),
          DIF_0.6 = dplyr::if_else(HET_DIF <= 0.6, "whitelist", "blacklist"),
          DIF_0.7 = dplyr::if_else(HET_DIF <= 0.7, "whitelist", "blacklist"),
          DIF_0.8 = dplyr::if_else(HET_DIF <= 0.8, "whitelist", "blacklist"),
          DIF_0.9 = dplyr::if_else(HET_DIF <= 0.9, "whitelist", "blacklist")
        ) %>%
        tidyr::gather(
          data = .,
          key = DIF_THRESHOLD,
          value = DIF_OUTLIERS,
          DIF_0.1:DIF_0.9
        ) %>%
        tidyr::unite(
          data = .,
          col = MAX_DIF_THRESHOLD, MAX_THRESHOLD, DIF_THRESHOLD,
          sep = "_",
          remove = FALSE
        ) %>% 
        dplyr::mutate(
          MAX_DIF_OUTLIERS = dplyr::if_else(
            MAX_OUTLIERS == "whitelist" & DIF_OUTLIERS == "whitelist",
            "whitelist",
            "blacklist")
        ) %>% 
        dplyr::group_by(LOCUS, MAX_DIF_THRESHOLD, MAX_THRESHOLD, DIF_THRESHOLD) %>% 
        dplyr::summarise(
          MAX_DIF_OUTLIERS = length(MAX_DIF_OUTLIERS[MAX_DIF_OUTLIERS == "whitelist"]),
          MAX_OUTLIERS = length(MAX_OUTLIERS[MAX_OUTLIERS == "whitelist"]),
          DIF_OUTLIERS = length(DIF_OUTLIERS[DIF_OUTLIERS == "whitelist"])
        ) %>% 
        dplyr::ungroup(.) %>%
        dplyr::mutate(
          MAX_DIF_OUTLIERS = 1 - MAX_DIF_OUTLIERS,
          MAX_OUTLIERS = 1 - MAX_OUTLIERS,
          DIF_OUTLIERS = 1 - DIF_OUTLIERS
        )
      
      max.threshold.overall <- dplyr::distinct(
        het.helper.overall, LOCUS, MAX_THRESHOLD, MAX_OUTLIERS
      ) %>%
        dplyr::group_by(MAX_THRESHOLD, MAX_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          MAX_OUTLIERS,
          nesting(MAX_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(MAX_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>% 
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>% 
        dplyr::group_by(MAX_THRESHOLD, MAX_OUTLIERS) %>% 
        tidyr::spread(data = ., key = MAX_THRESHOLD, value = MARKERS, fill = 0) %>% 
        dplyr::filter(MAX_OUTLIERS  == 0) %>% 
        dplyr::ungroup(.) %>% 
        dplyr::select(-MAX_OUTLIERS)
      
      dif.threshold.overall <- dplyr::distinct(
        het.helper.overall, LOCUS, DIF_THRESHOLD, DIF_OUTLIERS
      ) %>%
        dplyr::group_by(DIF_THRESHOLD, DIF_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          DIF_OUTLIERS,
          nesting(DIF_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(DIF_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>% 
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>% 
        dplyr::group_by(DIF_THRESHOLD, DIF_OUTLIERS) %>% 
        tidyr::spread(data = ., key = DIF_THRESHOLD, value = MARKERS, fill = 0) %>% 
        dplyr::filter(DIF_OUTLIERS  == 0) %>% 
        dplyr::ungroup(.) %>% 
        dplyr::select(-DIF_OUTLIERS)
      
      
      max.dif.threshold.combined.overall <- dplyr::distinct(
        het.helper.overall, LOCUS, MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS
      ) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(
          data = .,
          MAX_DIF_OUTLIERS,
          nesting(MAX_DIF_THRESHOLD),
          fill = list(n = 0)
        ) %>%
        dplyr::group_by(MAX_DIF_THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>% 
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>% 
        dplyr::group_by(MAX_DIF_THRESHOLD, MAX_DIF_OUTLIERS) %>% 
        tidyr::spread(data = ., key = MAX_DIF_THRESHOLD, value = MARKERS, fill = 0) %>% 
        dplyr::filter(MAX_DIF_OUTLIERS  == 0) %>% 
        dplyr::ungroup(.) %>% 
        dplyr::select(-MAX_DIF_OUTLIERS)
      
      helper.table.het <- list(
        helper.table.het.max.threshold = max.threshold, 
        helper.table.het.max.threshold.overall = max.threshold.overall, 
        helper.table.het.dif.threshold = dif.threshold,
        helper.table.het.dif.threshold.overall = dif.threshold.overall,
        helper.table.het.max.dif.threshold.combined = max.dif.threshold.combined,
        helper.table.het.max.dif.threshold.combined.overall = max.dif.threshold.combined.overall
        
      )
      if (interactive.filter) {
        readr::write_tsv(
          x = max.threshold, 
          path = stringi::stri_join(path.folder, "/helper.table.het.max.threshold.tsv")
        )
        readr::write_tsv(
          x = max.threshold.overall, 
          path = stringi::stri_join(path.folder, "/helper.table.het.max.threshold.overall.tsv")
        )
        
        readr::write_tsv(
          x = dif.threshold, 
          path = stringi::stri_join(path.folder, "/helper.table.het.dif.threshold.tsv")
        )
        readr::write_tsv(
          x = dif.threshold.overall, 
          path = stringi::stri_join(path.folder, "/helper.table.het.dif.threshold.overall.tsv")
        )
        readr::write_tsv(
          x = max.dif.threshold.combined, 
          path = stringi::stri_join(path.folder, "/helper.table.het.max.dif.threshold.tsv")
        )
        readr::write_tsv(
          x = max.dif.threshold.combined.overall, 
          path = stringi::stri_join(path.folder, "/helper.table.het.max.dif.threshold.overall.tsv")
        )
      }
      n.markers <- het.helper <- max.threshold <- max.threshold.overall <- dif.threshold <- dif.threshold.overall <- max.dif.threshold.combined <- max.dif.threshold.combined.overall <- NULL
    } else {# by SNP
      n.markers <- dplyr::n_distinct(het.summary$MARKERS)
      
      # by pop
      helper.table.het.pop <- dplyr::select(het.summary, MARKERS, POP_ID, HET_MEAN) %>% 
        dplyr::filter(POP_ID != "OVERALL") %>% 
        dplyr::mutate(
          `0.4` = dplyr::if_else(HET_MEAN <= 0.4, "whitelist", "blacklist"),
          `0.5` = dplyr::if_else(HET_MEAN <= 0.5, "whitelist", "blacklist"),
          `0.6` = dplyr::if_else(HET_MEAN <= 0.6, "whitelist", "blacklist"),
          `0.7` = dplyr::if_else(HET_MEAN <= 0.7, "whitelist", "blacklist"),
          `0.8` = dplyr::if_else(HET_MEAN <= 0.8, "whitelist", "blacklist"),
          `0.9` = dplyr::if_else(HET_MEAN <= 0.9, "whitelist", "blacklist")
        ) %>% 
        tidyr::gather(data = ., key = THRESHOLD, value = OUTLIERS, `0.4`:`0.9`) %>%
        dplyr::group_by(MARKERS, THRESHOLD, THRESHOLD) %>% 
        dplyr::summarise(OUTLIERS = length(OUTLIERS[OUTLIERS == "whitelist"])) %>% 
        dplyr::ungroup(.) %>%
        dplyr::mutate(OUTLIERS = pop.number - OUTLIERS) %>% 
        dplyr::group_by(THRESHOLD, OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., OUTLIERS, nesting(THRESHOLD), fill = list(n = 0)) %>%
        dplyr::group_by(THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = n.markers - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>% 
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>% 
        dplyr::group_by(THRESHOLD, OUTLIERS) %>% 
        tidyr::spread(data = ., key = THRESHOLD, value = MARKERS, fill = 0)
      
      # by overall
      helper.table.het.overall <- dplyr::select(het.summary.overall, MARKERS, HET_MEAN) %>% 
        dplyr::mutate(
          `0.4` = dplyr::if_else(HET_MEAN <= 0.4, "whitelist", "blacklist"),
          `0.5` = dplyr::if_else(HET_MEAN <= 0.5, "whitelist", "blacklist"),
          `0.6` = dplyr::if_else(HET_MEAN <= 0.6, "whitelist", "blacklist"),
          `0.7` = dplyr::if_else(HET_MEAN <= 0.7, "whitelist", "blacklist"),
          `0.8` = dplyr::if_else(HET_MEAN <= 0.8, "whitelist", "blacklist"),
          `0.9` = dplyr::if_else(HET_MEAN <= 0.9, "whitelist", "blacklist")
        ) %>% 
        tidyr::gather(data = ., key = THRESHOLD, value = OUTLIERS, `0.4`:`0.9`) %>%
        dplyr::group_by(MARKERS, THRESHOLD, THRESHOLD) %>% 
        dplyr::summarise(OUTLIERS = length(OUTLIERS[OUTLIERS == "whitelist"])) %>% 
        dplyr::ungroup(.) %>%
        dplyr::mutate(OUTLIERS = pop.number - OUTLIERS) %>% 
        dplyr::group_by(THRESHOLD, OUTLIERS) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., OUTLIERS, nesting(THRESHOLD), fill = list(n = 0)) %>%
        dplyr::group_by(THRESHOLD) %>%
        dplyr::mutate(
          WHITELIST = cumsum(n),
          BLACKLIST = 1 - WHITELIST,
          MARKERS = stringi::stri_join(WHITELIST, " (", BLACKLIST, ")")
        ) %>% 
        dplyr::select(-c(n, WHITELIST, BLACKLIST)) %>% 
        dplyr::group_by(THRESHOLD, OUTLIERS) %>% 
        tidyr::spread(data = ., key = THRESHOLD, value = MARKERS, fill = 0) %>% 
        dplyr::filter(OUTLIERS  == 0) %>% 
        dplyr::ungroup(.) %>% 
        dplyr::select(-OUTLIERS)
      
      helper.table.het <- list(
        helper.table.het.pop = helper.table.het.pop,
        helper.table.het.overall = helper.table.het.overall
      )
      
      if (interactive.filter) {
        readr::write_tsv(
          x = helper.table.het.pop,
          path = stringi::stri_join(path.folder, "/helper.table.het.threshold.pop.tsv")
        )
        readr::write_tsv(
          x = helper.table.het.overall,
          path = stringi::stri_join(path.folder, "/helper.table.het.threshold.overall.tsv")
        )
      }
      n.markers <- NULL
    }# End by SNP approach
  } else {
    helper.table.het <- "not selected"
  }# End helper table
  
  
  if (interactive.filter) {
    print(markers.pop.heterozygosity.density.plot)
    # save
    ggsave(stringi::stri_join(path.folder, "/markers.pop.heterozygosity.density.plot.pdf"), width = pop.number * 2, height = pop.number * 2, dpi = 600, units = "cm", useDingbats = F)
    ggsave(stringi::stri_join(path.folder, "/markers.pop.heterozygosity.density.plot.png"), width = pop.number * 2, height = pop.number * 2, dpi = 300, units = "cm")
    message(stringi::stri_join("2 versions (pdf and png) of the plot (markers.pop.heterozygosity.density.plot) were saved in this directory:\n", path.folder))
  }
  
  if (interactive.filter) {
    message("Show the manhattan plot of markers observed heterozygosity per populations and overall (y/n): ")
    manhattan.plot <- as.character(readLines(n = 1))
    if (manhattan.plot == "y") {
      message("Rendering the plot may take some time depending on the number of markers and populations...")
      print(markers.pop.heterozygosity.manhattan.plot)
      # save
      ggsave(stringi::stri_join(path.folder, "/markers.pop.heterozygosity.manhattan.plot.pdf"), width = pop.number * 2, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder, "/markers.pop.heterozygosity.manhattan.plot.png"), width = pop.number * 2, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (markers.pop.heterozygosity.manhattan.plot) were saved in this directory:\n", path.folder))
    }
  }
  
  if (interactive.filter) {
    message("Show the boxplot of markers observed heterozygosity per populations and overall (y/n): ")
    boxplot <- as.character(readLines(n = 1))
    if (boxplot == "y") {
      message("Rendering the plot may take some time depending on the number of markers and populations...")
      print(markers.pop.heterozygosity.boxplot)
      # save
      ggsave(stringi::stri_join(path.folder, "/markers.pop.heterozygosity.boxplot.pdf"), width = pop.number * 2, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder, "/markers.pop.heterozygosity.boxplot.png"), width = pop.number * 2, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (markers.pop.heterozygosity.boxplot) were saved in this directory:\n", path.folder))
    }
  }
  
  # Step 4: Blacklist markers based on observed heterozygosity------------------
  if (interactive.filter) {
    message("\n\nStep 4: Blacklist markers based on observed heterozygosity")
  }
  
  if (het.approach[1] == "haplotype") {
    message("Approach selected: haplotype")
    # het.threshold
    if (interactive.filter) {
      het.threshold <- 2
      while (isTRUE(het.threshold > 1)) {
        message("het.threshold argument:\n
With the haplotype approach, the het.threshold argument is independent of the
number of SNP/locus and is for filtering the mean maximum observed heterozygosity
of the haplotype. Anything higher (>) than the threshold will be blacklisted.\n")
        message("Example 1: LOCUS 9865 as 2 SNP on the 100 bp read, snp1 HET_OBS = 0.2 and snp2 HET_OBS = 0.9.
The max observed heterozygosity is 0.9, consequently, if you choose a threshold
of 0.5, all SNPs on this locus will be blacklisted, including snp1.\n")
        message("Example 2: LOCUS 2377 as 1 SNP on the 100 bp read, snp1 HET_OBS = 0.6.
The max observed heterozygosity is 0.6, consequently, if you choose a threshold
of 0.5, the SNP/locus/haplotype are blacklisted.\n")
        
        message("Markers will be discarded from the dataset if they don't pass the 2 next filters
(het.dif.threshold and het.pop.threshold).
Enter the het.threshold value (0 to 1, where 0.1 is very strict and 1 turns off the filter):")
        het.threshold <- as.double(readLines(n = 1))
      }
    }
    # het.dif.threshold
    if (interactive.filter) {
      het.dif.threshold <- 2
      while (isTRUE(het.dif.threshold > 1)) {
        message("het.dif.threshold argument:\n
Markers with 1 SNP/read are not affected by this threshold.
Locus/haplotypes with > 1 SNP are inspected for consistencies of observed heterozygosity along the read.
Locus (read/haplotypes and all SNP on it) are blacklisted if the difference in observed het > threshold.\n")
        message("Example: LOCUS 7654 as 2 SNP, snp1 HET_OBS = 0.1 and snp2 HET_OBS = 0.5.
If you choose a threshold of 0.2 all SNP on this locus will be blacklisted.
Choosing a threshold of 0.6 will not blacklist the SNPs and it's locus.
The locus and it's SNPs will be discarded if they don't pass the remaining filter (het.pop.threshold)\n")
        message("Enter the het.dif.threshold value (0 to 1, where 0.1 is very strict and 1 turn off the filter):")
        het.dif.threshold <- as.double(readLines(n = 1))
      }
    }
    #outlier.pop.threshold
    if (interactive.filter) {
      if (het.approach[2] == "pop") {
        outlier.pop.threshold <- pop.number + 10
        while (isTRUE(outlier.pop.threshold > pop.number)) {
          message("outlier.pop.threshold argument:\n
This filter works by counting the number of populations that didn't pass
the previous 2 filters (het.threshold, het.dif.threshold).
If the number of populations is higher (>) than the threshold, the SNPs and it's locus are discarded.\n")
          message("Useful to incorporate problematic populations dragging down 
polymorphism discovery, but still wanted for analysis.
Use this threshold to allow variance in the number of populations passing 
the previous thresholds. Blacklist and whitelist produced by the filter 
allows to manage outlier markers, individuals and populations.\n")
          message("Example: with a outlier.pop.threshold = 2, you tolerate a
maximum of 2 outlier populations (failing het.threshold and/or het.dif.threshold).
The lower the number, the more severe the filtering, while entering the 
number of populations in the dataset turns off the filter.\n")
          message("Enter the outlier.pop.threshold:")
          
          outlier.pop.threshold <- as.numeric(readLines(n = 1))
        }
      }
    }
    
    # het.threshold <- 0.5 # test
    # het.dif.threshold <- 0.5# test
    # outlier.pop.threshold <- 2# test
    if (het.approach[2] == "pop") {
      filter <- dplyr::ungroup(het.summary) %>%
        dplyr::group_by(LOCUS) %>% 
        dplyr::summarise(OUTLIERS = length(POP_ID[HET_DIF > het.dif.threshold | HET_MAX > het.threshold])) %>%
        dplyr::filter(OUTLIERS <= outlier.pop.threshold) %>%
        dplyr::select(LOCUS) %>%
        dplyr::left_join(input, by = "LOCUS") %>%
        dplyr::arrange(LOCUS, POP_ID)
    } else {
      filter <- dplyr::ungroup(het.summary.overall) %>%
        dplyr::group_by(LOCUS) %>% 
        dplyr::filter(HET_DIF <= het.dif.threshold) %>%
        dplyr::filter(HET_MAX <= het.threshold) %>%
        dplyr::select(LOCUS) %>%
        dplyr::left_join(input, by = "LOCUS") %>%
        dplyr::arrange(LOCUS, POP_ID)
    }
  } else {# by SNP
    message("Approach selected: SNP")
    # het.threshold
    if (interactive.filter) {
      het.threshold <- 2
      while (isTRUE(het.threshold > 1)) {
        message("het.threshold argument:\n
With the SNP approach, the het.threshold argument is independent of the
number of SNP/locus and is for filtering the observed heterozygosity
of each SNP independently of it's haplotype.\n")
        message("Example 1: LOCUS 9865 as 2 SNP on the 100 bp read, snp1 HET_OBS = 0.2 and snp2 HET_OBS = 0.9.
If you choose a threshold of 0.5, only snp2 is blacklisted.\n")
        message("Markers will be discarded from the dataset if they don't pass
the next filters (het.pop.threshold).
Enter the het.threshold value (0 to 1, where 0.1 is very strict and 1 turns off the filter):")
        het.threshold <- as.double(readLines(n = 1))
      }
    }
    #outlier.pop.threshold
    if (interactive.filter) {
      if (het.approach[2] == "pop") {
        outlier.pop.threshold <- pop.number + 10
        while (isTRUE(outlier.pop.threshold > pop.number)) {
          message("outlier.pop.threshold argument:\n
This filter works by counting the number of populations that didn't pass
the previous filter (het.threshold).
If the number of populations is higher (>) than the threshold, the SNPs and it's locus are discarded.\n")
          message("Useful to incorporate problematic populations dragging down 
polymorphism discovery, but still wanted for analysis.
Use this threshold to allow variance in the number of populations passing 
the previous thresholds. Blacklist and whitelist produced by the filter 
allows to manage outlier markers, individuals and populations.\n")
          message("Example: with a outlier.pop.threshold = 2, you tolerate a
maximum of 2 outlier populations (failing het.threshold and/or het.dif.threshold).
The lower the number, the more severe the filtering, while entering the 
number of populations in the dataset turns off the filter.\n")
          message("Enter the outlier.pop.threshold:")
          
          outlier.pop.threshold <- as.numeric(readLines(n = 1))
          
        }
      }
    }
    # het.threshold <- 0.5 # test
    # outlier.pop.threshold <- 2# test
    if (het.approach[2] == "pop") {
      filter <- dplyr::ungroup(het.summary) %>%
        dplyr::group_by(MARKERS) %>% 
        dplyr::summarise(OUTLIERS = length(POP_ID[HET_MEAN > het.threshold])) %>% #using HET_MEAN is the same with snp approach
        dplyr::filter(OUTLIERS <= outlier.pop.threshold) %>%
        dplyr::select(MARKERS) %>%
        dplyr::left_join(input, by = "MARKERS") %>%
        dplyr::arrange(MARKERS, POP_ID)
    } else {
      filter <- dplyr::ungroup(het.summary.overall) %>%
        dplyr::group_by(MARKERS) %>% 
        dplyr::filter(HET_MEAN <= het.threshold) %>%
        dplyr::select(MARKERS) %>%
        dplyr::left_join(input, by = "MARKERS") %>%
        dplyr::arrange(MARKERS, POP_ID)
    }
  }# end snp approach to filtering
  
  # Update filters.parameters SNP ----------------------------------------------
  # Prepare a list of markers and number of markers before filtering
  if (tibble::has_name(het.summary, "LOCUS") & het.approach[1] == "haplotype") {
    snp.before <- dplyr::n_distinct(input$POS)
    locus.before <- dplyr::n_distinct(input$LOCUS)
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(snp.before - snp.after)
    locus.after <- as.integer(dplyr::n_distinct(filter$LOCUS))
    locus.blacklist <- as.integer(locus.before - locus.after)
  } else {
    snp.before <- dplyr::n_distinct(input$MARKERS)
    snp.after <- as.integer(dplyr::n_distinct(filter$MARKERS))
    snp.blacklist <- as.integer(snp.before - snp.after)
    locus.before <- as.character("NA")
    locus.after <- as.character("NA")
    locus.blacklist <- as.character("NA")
  }
  
  ind.before <- dplyr::n_distinct(input$INDIVIDUALS)
  ind.blacklisted <- length(blacklist.ind.het$INDIVIDUALS)
  ind.after <- ind.before - ind.blacklisted
  
  markers.before <- stringi::stri_join(snp.before, locus.before, sep = "/")
  markers.after <- stringi::stri_join(snp.after, locus.after, sep = "/")
  markers.blacklist <- stringi::stri_join(snp.blacklist, locus.blacklist, sep = "/")
  
  if (tibble::has_name(het.summary, "LOCUS")) {
    markers.df <- dplyr::distinct(input, CHROM, LOCUS, POS)
  } else {
    markers.df <- dplyr::distinct(input, MARKERS)
  }
  
  
  if (het.approach[2] == "overall") outlier.pop.threshold <- "using overall"
  
  filters.parameters <- tibble::data_frame(
    FILTERS = c("Observed Heterozygosity", rep(as.character(""), 4)),
    PARAMETERS = c("ind.mean.het", "het.approach", "het.threshold", "het.dif.threshold", "outlier.pop.threshold"), 
    VALUES = c(ind.heterozygosity.threshold, paste(het.approach, collapse = " and "), paste("<=", het.threshold), paste("<=", het.dif.threshold), paste("<=", outlier.pop.threshold)), 
    BEFORE = c(ind.before, "", "", "", markers.before),
    AFTER = c(ind.after, "", "", "", markers.after),
    BLACKLIST = c(ind.blacklisted, "", "", "", markers.blacklist),
    UNITS = c("individual" , "", "", "", "SNP/LOCUS"),
    COMMENTS = c("", "", "", "", "")
  )
  readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  
  # saving filtered tidy data --------------------------------------------------
  # filename <- "test.tidy.tsv"#test
  if (!is.null(filename)) {
    message("Writing the filtered tidy data set in your working directory...")
    readr::write_tsv(filter, paste0(path.folder,"/", filename), append = FALSE, col_names = TRUE)
  }
  # saving whitelist -----------------------------------------------------------
  message("Writing the whitelist of markers in your working directory\nwhitelist.markers.het.tsv")
  
  if (tibble::has_name(het.summary, "LOCUS")) {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(CHROM, LOCUS, POS)
  } else {
    whitelist.markers <- dplyr::ungroup(filter) %>%
      dplyr::distinct(MARKERS)
  }
  readr::write_tsv(whitelist.markers, paste0(path.folder,"/whitelist.markers.het.tsv"), append = FALSE, col_names = TRUE)
  
  
  # saving blacklist -----------------------------------------------------------
  message("Writing the blacklist of markers in your working directory\nblacklist.markers.het.tsv")
  if (tibble::has_name(het.summary, "LOCUS")) {
    blacklist.markers <- dplyr::anti_join(markers.df, whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  } else {
    blacklist.markers <- dplyr::anti_join(markers.df, whitelist.markers, by = "MARKERS")
  }
  readr::write_tsv(blacklist.markers, paste0(path.folder,"/blacklist.markers.het.tsv"), append = FALSE, col_names = TRUE)
  
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stringi::stri_join("ind.heterozygosity.threshold: ", ind.heterozygosity.threshold))
  message(stringi::stri_join("Blacklisted individuals: ", ind.blacklisted))
  message(stringi::stri_join("het.approach: ", paste(het.approach, collapse = " and ")))
  message(stringi::stri_join("het.threshold: ", het.threshold))
  message(stringi::stri_join("het.dif.threshold: ", het.dif.threshold))
  message(stringi::stri_join("outlier.pop.threshold: ", outlier.pop.threshold))
  if (tibble::has_name(het.summary, "LOCUS") & het.approach[1] == "haplotype") {
    message(stringi::stri_join("The number of markers removed by the HET filter:\nSNP: ", snp.before - dplyr::n_distinct(filter$POS), "\nLOCUS: ", locus.before - dplyr::n_distinct(filter$LOCUS)))
    message("The number of markers before -> after the HET filter")
    message(stringi::stri_join("SNP: ", snp.before, " -> ", as.integer(dplyr::n_distinct(filter$POS))))
    message(stringi::stri_join("LOCUS: ", locus.before, " -> ", as.integer(dplyr::n_distinct(filter$LOCUS))))
  } else {# for haplotype file
    message(stringi::stri_join("The number of markers/locus removed by the HET filter: ", snp.before - dplyr::n_distinct(filter$MARKERS)))
    message("The number of markers before -> after the HET filter")
    message(stringi::stri_join("MARKERS/LOCUS: ", snp.before, " -> ", as.integer(dplyr::n_distinct(filter$MARKERS))))
  }
  if (!interactive.filter) {
    timing <- proc.time() - timing
    message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
  }
  cat("############################## completed ##############################\n")
  res <- list()
  res$tidy.filtered.het <- filter
  res$whitelist.markers <- whitelist.markers
  res$blacklist.markers <- blacklist.markers
  res$strata <- strata.df
  res$filters.parameters <- filters.parameters
  res$individual.heterozigosity <- het.ind
  res$blacklist.ind.het <- blacklist.ind.het
  res$helper.table.het <- helper.table.het 
  res$individual.heterozygosity.boxplot <- individual.heterozygosity.boxplot
  res$individual.heterozygosity.manhattan.plot <- individual.heterozygosity.manhattan.plot
  res$markers.pop.heterozygosity.boxplot <- markers.pop.heterozygosity.boxplot
  res$markers.pop.heterozygosity.density.plot <- markers.pop.heterozygosity.density.plot
  res$markers.pop.heterozygosity.manhattan.plot <- markers.pop.heterozygosity.manhattan.plot
  return(res)
}
