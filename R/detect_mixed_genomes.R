# detect mixed genomes
#' @title Detect mixed genomes
#' @description Highlight outliers individual's heterozygosity for a quick
#' diagnostic of mixed samples.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{stackr} \code{\link{tidy_genomic_data}}.



#' @param ind.heterozygosity.threshold (double, optional)
#' Blacklist individuals based on observed heterozygosity of individuals 
#' (averaged across markers).
#' The value for the threshold is a proportion (0 to 1), where 1 is similar to the
#' default, \code{ind.heterozygosity.threshold = NULL}, and turn off the filter
#' (the function will only output the plots and table of heterozygosity).
#' Individuals with mean heterozygosity higher (>) than the threshold
#' will be blacklisted. 

#' @return The function returns inside the global environment a list with
#' 4 objects:
#' 
#' \enumerate{

#' \item the individual's heterozigosity:\code{$individual.heterozigosity}
#' \item the blacklisted individuals based on the individual's heterozigosity:\code{$blacklist.ind.het}
#' \item the boxplot of individual heterozygosity:\code{$individual.heterozygosity.boxplot}
#' \item the manhattan plot of individual heterozygosity:\code{$individual.heterozygosity.manhattan.plot}
#' }
#' 
#' 
#' In the working directory, output is conditional to \code{interactive.filter} argument

#' @rdname detect_mixed_genomes
#' @export

#' @import ggplot2
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame
#' @importFrom tidyr complete gather unite spread nesting
#' @importFrom stats median sd

#' @seealso \link{plot_density_distribution_het}

#' @examples
#' \dontrun{
#' #Step1: highlight outlier individuals, the simplest way to run:
#' outlier.ind.het <- detect_mixed_genomes(data = "wombat_tidy.tsv")
#' 
#' # this example without threshold will not produce a blacklist of individuals.
#' 
#' # to look at the table with individual's heterozygosity
#' outlier.ind.het$individual's heterozygosity
#' 
#' # To view the manhattan plot:
#' outlier.ind.het$individual.heterozygosity.manhattan.plot
#' 
#' # To view the box plot
#' outlier.ind.het$$individual.heterozygosity.boxplot
#' 
#' # To save the boxplot:
#' ggsave("individual.heterozygosity.boxplot.pdf", width = 15, height = 10, dpi = 600, units = "cm", useDingbats = F)
#' # prefer a PNG:
#' ggsave("individual.heterozygosity.boxplot.png", width = 15, height = 10, dpi = 300, units = "cm")
#' 
#' 
#' # Based on the look of the distribution using both jitter and boxplot,
#' # choose a threshold to blacklist the outliers and re-run the function.
#' 
#' This can be done in one step with the interactive mode
#' in \code{\link{filter_het}}.
#' 
#' readr::write_tsv(x = blacklist.ind.het, path = "blacklist.individuals.heterozygosity.tsv", col_names = TRUE)
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_mixed_genomes <- function(
  data,
  ind.heterozygosity.threshold = NULL
) {
  cat("#######################################################################\n")
  cat("#################### stackr::detect_mixed_genomes #####################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  
  # manage missing arguments ---------------------------------------------------
  if (missing(data)) stop("missing data argument")
  
  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- stackr::tidy_wide(data = data, import.metadata = FALSE)
  } else {
    input <- data
  }
  
  # check genotype column naming
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input), 
    pattern = "GENOTYPE", 
    replacement = "GT", 
    vectorize_all = FALSE
  )
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # # strata
  # strata <- input %>% 
  #   ungroup %>% 
  #   distinct(POP_ID, INDIVIDUALS)
  # 
  
  # highlight heterozygote
  het.summary <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, GT) %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::mutate(
      HET = dplyr::if_else(
        stringi::stri_sub(GT, 1, 3) != stringi::stri_sub(GT, 4, 6), 1, 0
      )
    )
  
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
    ) %>%
    dplyr::arrange(POP_ID, HET_PROP) %>% 
    dplyr::ungroup(.)
  
  het.ind.overall <- dplyr::mutate(.data = het.ind, POP_ID = as.character(POP_ID)) %>%
    dplyr::bind_rows(dplyr::mutate(.data = het.ind, POP_ID = rep("OVERALL", n()))) %>%
    dplyr::mutate(POP_ID = factor(POP_ID, levels = c(levels(het.summary$POP_ID), "OVERALL")))
    
  
  # Get stats...
  
  het.ind.stats <- het.ind.overall %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      HET_MEAN = mean(HET_PROP, na.rm = TRUE),
      HET_MEDIAN = stats::median(HET_PROP, na.rm = TRUE),
      HET_SD = stats::sd(HET_PROP, na.rm = TRUE),
      HET_MIN = min(HET_PROP, na.rm = TRUE),
      HET_MAX = max(HET_PROP, na.rm = TRUE)
    ) %>% 
    dplyr::mutate_if(.tbl = ., .predicate = is.numeric, .funs = round, digits = 4) %>% 
    tidyr::unite(data = ., HET_RANGE, HET_MIN, HET_MAX, sep = " - ") %>% 
    dplyr::arrange(POP_ID, HET_MEAN)
  
  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }
  y.breaks.by <- rounder(max(het.ind$HET_PROP)/10, 0.001, ceiling)
  y.breaks.max <- rounder(max(het.ind$HET_PROP), 0.001, ceiling)
  y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)

  
  individual.heterozygosity.manhattan.plot <- ggplot(data = het.ind.overall, aes(x = POP_ID, y = HET_PROP, colour = POP_ID)) + 
    geom_jitter() + 
    labs(y = "Individual's Mean Heterozygosity (proportion)") +
    # labs(x = "Populations") +
    # labs(colour = "Populations") +
    scale_y_continuous(name = waiver(), breaks = y.breaks, limits = c(0, y.breaks.max), expand = c(0.06, 0)) +
    # theme_minimal() +
    theme(
      legend.position = "none",
      # panel.grid.major.y = element_line(linetype = "solid"),
      # panel.grid.minor.y = element_line(linetype = "longdash", size = 1),
      # panel.background = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      # axis.text.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    ) +
    geom_hline(mapping = aes(yintercept = HET_MEAN), het.ind.stats, linetype = "dotted", size = 0.6) + #mean
    # geom_hline(mapping = aes(yintercept = HET_sig_minus), het.ind.stats.pop, linetype = "dashed") + #3 sigma -
    # geom_hline(mapping = aes(yintercept = HET_sig_plus), het.ind.stats.pop, linetype = "dashed") + #3 sigma +
    facet_grid(~ POP_ID, switch = "x", scales = "free")
  # individual.heterozygosity.manhattan.plot
  
  individual.heterozygosity.boxplot <- ggplot(data = het.ind.overall, aes(x = POP_ID, y = HET_PROP, colour = POP_ID)) + 
    geom_boxplot() + 
    labs(y = "Individual's Mean Heterozygosity (proportion)") +
    labs(x = "Populations") +
    labs(colour = "Populations") +
    scale_y_continuous(name = waiver(), breaks = y.breaks, limits = c(0, y.breaks.max), expand = c(0.06, 0)) +
    theme_classic() +
    # theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    )
  # individual.heterozygosity.boxplot
  
  ## Step 2: Blacklist outlier individuals -------------------------------------
  # Blacklist individuals based a threshold of mean heterozygosity
  if (!is.null(ind.heterozygosity.threshold)) {
    # ind.heterozygosity.threshold <- c(0.035, 0.10)
    threshold.min <- ind.heterozygosity.threshold[1]
    threshold.max <- ind.heterozygosity.threshold[2]
    
    blacklist.ind.het  <- dplyr::ungroup(het.ind) %>%
      dplyr::filter(HET_PROP > threshold.max | HET_PROP < threshold.min) %>% 
      dplyr::distinct(INDIVIDUALS)
    
    message(stringi::stri_join("Filter individual's heterozygosity: ", length(blacklist.ind.het$INDIVIDUALS), " individual(s) blacklisted"))

  } else {
    blacklist.ind.het <- "ind.heterozygosity.threshold is necessary to get a blacklist of individuals"
  }
  message(stringi::stri_join("Computation time: ", round((proc.time() - timing)[[3]]), " sec"))
  cat("############################## completed ##############################\n")
  res <- list(
    individual.heterozygosity = het.ind,
    individual.heterozygosity.statistics = het.ind.stats,
    blacklist.ind.het = blacklist.ind.het,
    individual.heterozygosity.boxplot = individual.heterozygosity.boxplot,
    individual.heterozygosity.manhattan.plot = individual.heterozygosity.manhattan.plot
  )
  return(res)
}
