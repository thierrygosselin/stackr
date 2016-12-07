# Detect duplicate genomes

#' @name detect_duplicate_genomes
#' @title Compute pairwise genome similarity or distance between individuals 
#' to highligh potential duplicate individuals
#' @description The function can compute two methods 
#' to highligh potential duplicate individuals.
#' \enumerate{
#' \item distance between individuals and/or
#' \item pairwise genome similarity
#' }

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{stackr} \code{\link{tidy_genomic_data}}.

#' @param distance.method (character) The distance measure used inside \code{stats::dist} 
#' (<= 30000 markers) or \code{amap::Dist} (> 30000 markers). 
#' This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary".
#' Using \code{distance.method = NULL} will not run this method.
#' Default: \code{distance.method = "manhattan"}. This is very fast
#' compared to the genome similarity method. It uses allele counts and the codes
#' are tailored for biallelic and multiallelic markers.

#' @param genome (logical) Computes pairwise genome similarity in parallel.
#' The proportion of the shared genotypes is averaged across shared markers between 
#' each pairwise comparison. This method makes filtering easier because the 
#' threshold is more intuitive with the plots produced, but it's much longer 
#' to run, even in parallel, so better to run overnight. 
#' Default: \code{genome = FALSE}.

#' @param parallel.core (optional) The number of core for parallel computation.
#' Default: \code{parallel.core = parallel::detectCores()-1}.

#' @return A list with potentially 8 objects: 
#' \code{$distance }: results of the distance method
#' \code{$distance.stats}: Summary statistics of the distance method
#' \code{$pairwise.genome.similarity}: results of the genome method
#' \code{$genome.stats}: Summary statistics of the genome method
#' \code{$violin.plot.distance}: violin plot showing the distribution of pairwise distances
#' \code{$manhattan.plot.distance}: same info different visual with manhattan plot
#' \code{$violin.plot.genome}: violin plot showing the distribution of pairwise genome similarities
#' \code{$manhattan.plot.genome}: same info different visual with manhattan plot
#' 
#' Saved in the working directory:
#' individuals.pairwise.dist.tsv, individuals.pairwise.distance.stats.tsv, 
#' individuals.pairwise.genome.similarity.tsv, individuals.pairwise.genome.stats.tsv

#' @details
#' Strategically, run the default first (\code{distance.method},
#' but no \code{genome})
#' 
#' \strong{\code{distance.method} argument is fast, but...}
#' 
#' you don't know if the observed comparison (close or distant)
#' is influenced by missing values/the number of markers in common
#' between the pair compared. This is something that needs to be considered.
#' Be suspicious of a \emph{distant outlier} from the same pop pairwise comparison,
#' and similarly, be suspicious of a \emph{close outlier} from a different pop
#' pairwise comparisons.
#' 
#' If there is no outliers, don't bother with the other method (\code{genome = TRUE}).
#' 
#' 
#' \strong{\code{genome = TRUE} argument is slower, but...}
#' If you see outliers with the first run, take the time to run the function
#' with \code{genome = TRUE}. Because this option is much better at detecting
#' duplicated individuals and it also shows the impact of \strong{missingness}
#' or the number of \strong{shared markers} between comparisons.
#' 
#' Your outlier duo could well be the result of one of the individual having
#' an extremely low number genotypes...


#' @export
#' @rdname detect_duplicate_genomes
#' @importFrom stringi stri_paste stri_replace_all_fixed
#' @importFrom dplyr arrange rename select group_by filter mutate rename_ filter_ bind_cols bind_rows summarise n_distinct intersect
#' @importFrom utils combn
#' @importFrom stats na.omit var median quantile dist
#' @importFrom amap Dist
#' @importFrom reshape2 melt
#' @importFrom lazyeval interp
#' @importFrom readr write_tsv
#' @importFrom parallel detectCores
#' @importFrom purrr flatten map
#' @importFrom tibble as_data_frame has_name remove_rownames column_to_rownames
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light


#' @examples
#' \dontrun{
#' # First run and simplest way (if you have the tidy df):
#' dup <- stackr::detect_duplicate_genomes(data = "wombat_tidy.tsv")
#' 
#' #If you need a tidy df:
#' dup <- stackr::tidy_genomic_data(
#' data = "wombat_tidy.tsv",
#' strata = "wombat.strata.tsv",
#' vcf.metadata = FALSE
#' ) %>%
#' stackr::detect_duplicate_genomes(data = .)
#' 
#' # This will use by defaul:
#' distance.method = "manhattan"
#' genome = FALSE
#' #parallel.core = all my CPUs - 1 
#' 
#' # To view the manhattan plot:
#' dup$manhattan.plot.distance
#' 
#' # to view the data stats
#' dup.data.stats <- dup$distance.stats
#' 
#' # to view the data
#' dup.data <- dup$distance
#' 
#' # Based on the look of the distribution using both manhattan and boxplot, 
#' # I can filter the dataset to highlight potential duplicates: 
#' dup.filtered <- dplyr::filter(.data = dup.data, DISTANCE_RELATIVE < 0.2)
#' 
#' # To run the distance (with euclidean distance instead of the default manhattan, 
#' # and also carry the second analysis (with the genome method):
#' dup <- stackr::tidy_genomic_data(
#' data = "wombat_tidy.tsv",
#' strata = "wombat.strata.tsv",
#' vcf.metadata = FALSE
#' ) %>%
#' stackr::detect_duplicate_genomes(data = ., distance.method = "euclidean", genome = TRUE)
#' 
#' # to view the data of the genome data
#' dup.data <- dup$pairwise.genome.similarity
#' 
#' # Based on the look of the distribution using both manhattan and boxplot, 
#' # I can filter the dataset based on 98% of identical genotype proportion, 
#' # to highlight potential duplicates: 
#' dup.filtered <- dplyr::filter(.data = dup.data, PROP_IDENTICAL > 0.98)
#' 
#' # Get the list of duplicates id
#' dup.list.names <- data.frame(INDIVIDUALS = unique(c(dup.filtered$ID1, dup.filtered$ID2)))
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_duplicate_genomes <- function(
  data,
  distance.method = "manhattan",
  genome = FALSE,
  parallel.core = parallel::detectCores() - 1
) {
  cat("\n\n")
  cat("###############################################################################\n")
  cat("###################### stackr::detect_duplicate_genomes #######################\n")
  cat("###############################################################################\n")
  timing <- proc.time()
  
  # Manage missing arguments ---------------------------------------------------
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
  
  # strata
  strata <- dplyr::ungroup(input) %>% 
    dplyr::distinct(POP_ID, INDIVIDUALS)
  
  # New list to prepare for results
  res <- list()
  
  # Preparing data for comparisons ---------------------------------------------
  message("Preparing data for analysis")
  
  # GT_BIN available
  if (!is.null(distance.method) & tibble::has_name(input, "GT_BIN")) {
    input.prep <- dplyr::ungroup(input) %>%
      dplyr::select(MARKERS, INDIVIDUALS, ALT = GT_BIN) %>%
      dplyr::mutate(REF = 2 - ALT) %>%
      tidyr::gather(data = ., key = ALLELES, value = n, -c(MARKERS,INDIVIDUALS)) %>% 
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
      dplyr::select(-MARKERS, -ALLELES) %>% 
      dplyr::arrange(MARKERS_ALLELES,INDIVIDUALS)
  }
  
  # GT_BIN NOT available
  if (!tibble::has_name(input, "GT_BIN")) {
    # Allele count
    missing.geno <- dplyr::select(.data = input, MARKERS, INDIVIDUALS, GT) %>% 
      dplyr::filter(GT == "000000") %>%
      dplyr::select(-GT) %>% 
      dplyr::mutate(MISSING = rep("blacklist", n()))
    
    input.prep <- dplyr::ungroup(input) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>% 
      dplyr::filter(GT != "000000") %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
        A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
      ) %>% 
      dplyr::select(-GT) %>% 
      tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS)) %>% 
      dplyr::arrange(MARKERS, INDIVIDUALS, GT) %>%
      dplyr::count(x = ., INDIVIDUALS, MARKERS, GT) %>% 
      dplyr::ungroup(.) %>% 
      tidyr::complete(data = ., INDIVIDUALS, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
      dplyr::anti_join(missing.geno, by = c("MARKERS", "INDIVIDUALS")) %>% 
      dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, GT, sep = ".")) %>%
      dplyr::select(-MARKERS, -GT) %>%
      dplyr::arrange(MARKERS_ALLELES, INDIVIDUALS)
    
    missing.geno <- NULL # unused object
  }#end preparing data
  
  # Computing distance ---------------------------------------------------------
  # distance.method <- "euclidean"
  # parallel.core <- 8
  
  if (!is.null(distance.method)) {
    message(stringi::stri_join("Computing ", distance.method, " distances between individuals..."))
    
    res$distance <- distance_individuals(
      x = input.prep,
      strata = strata,
      distance.method = distance.method, 
      parallel.core = parallel.core
    )
    
    readr::write_tsv(
      x = res$distance,
      path = "individuals.pairwise.dist.tsv",
      col_names = TRUE
    )
    
    # Stats
    message("Generating summary statistics")
    res$distance.stats <- res$distance %>% 
      dplyr::summarise(
        MEAN = mean(DISTANCE_RELATIVE, na.rm = TRUE),
        MEDIAN = stats::median(DISTANCE_RELATIVE, na.rm = TRUE),
        SE = round(sqrt(stats::var(DISTANCE_RELATIVE, na.rm = TRUE)/length(stats::na.omit(DISTANCE_RELATIVE))), 2),
        MIN = round(min(DISTANCE_RELATIVE, na.rm = TRUE), 2),
        MAX = round(max(DISTANCE_RELATIVE, na.rm = TRUE), 2),
        QUANTILE25 = stats::quantile(DISTANCE_RELATIVE, 0.25), # quantile25
        QUANTILE75 = stats::quantile(DISTANCE_RELATIVE, 0.75)#, # quantile75
        # OUTLIERS_LOW = QUANTILE25 - (1.5 * (QUANTILE75 - QUANTILE25)), # outliers : below the outlier boxplot
        # OUTLIERS_HIGH = QUANTILE75 + (1.5 * (QUANTILE75 - QUANTILE25)) # outliers : higher the outlier boxplot
      )
    readr::write_tsv(
      x = res$distance.stats, 
      path = "individuals.pairwise.distance.stats.tsv", 
      col_names = TRUE
    )
    
    message("Generating plots")
    # violin plot
    res$violin.plot.distance <- ggplot2::ggplot(
      data = res$distance, 
      ggplot2::aes(x = PAIRWISE, y = DISTANCE_RELATIVE, na.rm = TRUE)
    ) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = "black") +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(y = "Distance (relative)\n <- close     distant ->") +
      ggplot2::labs(x = "Pairwise comparisons") +
      ggplot2::theme(
        # legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        # panel.grid.major.y = element_blank(), 
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
      )
    
    # Manhattan plot
    res$manhattan.plot.distance <- ggplot2::ggplot(
      data = res$distance,
      ggplot2::aes(x = PAIRWISE, y = DISTANCE_RELATIVE, colour = POP_COMP)
    ) + 
      ggplot2::geom_jitter(alpha = 0.3) + 
      ggplot2::labs(y = "Distance (relative)\n <- close     distant ->") +
      ggplot2::labs(x = "Pairwise comparisons") +
      ggplot2::labs(colour = "Population comparisons") +
      ggplot2::scale_colour_manual(values = c("#0571b0", "black")) +
      # ggplot2::scale_y_reverse() +
      ggplot2::theme_light() +
      ggplot2::theme(
        # legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        # panel.grid.major.y = element_blank(), 
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
      )
  } # end distance method
  
  # Compute genome similarity -------------------------------------------------
  if (genome) {
    
    # If GT_BIN available, we need a new input.prep (not the same as dist method)
    if (tibble::has_name(input, "GT_BIN")) {
      input.prep <- dplyr::filter(.data = input, !is.na(GT_BIN))
    }
    
    input <- NULL
    
    # list of id
    id.list <- unique(input.prep$INDIVIDUALS) # id list
    
    # all combination of individual pair
    id.pairwise <- utils::combn(unique(id.list), 2, simplify = FALSE) 
    list.pair <- 1:length(id.pairwise)

    message(stringi::stri_join("Pairwise comparisons: ", length(list.pair)))
    message("Starting scan for duplicate genomes, take a break...")
    
    # list.pair <- 1:5 # test
    
    pairwise.genome.similarity <- .stackr_parallel(
      X = list.pair, 
      FUN = pairwise_genome_similarity, 
      mc.preschedule = FALSE, 
      mc.silent = FALSE, 
      mc.cleanup = TRUE,
      mc.cores = parallel.core,
      id.pairwise = id.pairwise,
      input.prep = input.prep
    )
    pairwise.genome.similarity <- dplyr::bind_rows(pairwise.genome.similarity)
    
    # Include population info with strata
    ID1.pop <- suppressWarnings(
      pairwise.genome.similarity %>% 
        dplyr::select(INDIVIDUALS = ID1) %>% 
        dplyr::inner_join(strata, by = "INDIVIDUALS") %>% 
        dplyr::select(ID1_POP = POP_ID))
    
    ID2.pop <- suppressWarnings(
      pairwise.genome.similarity %>% 
        dplyr::select(INDIVIDUALS = ID2) %>% 
        dplyr::inner_join(strata, by = "INDIVIDUALS") %>% 
        dplyr::select(ID2_POP = POP_ID))
    
    pairwise.genome.similarity <- dplyr::bind_cols(
      pairwise.genome.similarity, ID1.pop, ID2.pop
    ) %>% 
      dplyr::mutate(
        POP_COMP = ifelse(ID1.pop == ID2.pop, "same pop", "different pop"),
        POP_COMP = factor(POP_COMP, levels = c("same pop", "different pop"), ordered = TRUE),
        PAIRWISE = rep("pairwise comparison", n()),
        METHOD = rep("genome similarity", n())
      )
    
    res$pairwise.genome.similarity <- pairwise.genome.similarity
    
    readr::write_tsv(
      x = pairwise.genome.similarity, 
      path = "individuals.pairwise.genome.similarity.tsv", 
      col_names = TRUE
    )
    
    # Stats
    message("Generating summary statistics")
    res$genome.stats <- pairwise.genome.similarity %>% 
      dplyr::summarise(
        MEAN = mean(PROP_IDENTICAL, na.rm = TRUE),
        MEDIAN = stats::median(PROP_IDENTICAL, na.rm = TRUE),
        SE = round(sqrt(stats::var(PROP_IDENTICAL, na.rm = TRUE)/length(stats::na.omit(PROP_IDENTICAL))), 2),
        MIN = round(min(PROP_IDENTICAL, na.rm = TRUE), 2),
        MAX = round(max(PROP_IDENTICAL, na.rm = TRUE), 2),
        QUANTILE25 = stats::quantile(PROP_IDENTICAL, 0.25), # quantile25
        QUANTILE75 = stats::quantile(PROP_IDENTICAL, 0.75)#, # quantile75
        # OUTLIERS_LOW = QUANTILE25 - (1.5 * (QUANTILE75 - QUANTILE25)), # outliers : below the outlier boxplot
        # OUTLIERS_HIGH = QUANTILE75 + (1.5 * (QUANTILE75 - QUANTILE25)) # outliers : higher the outlier boxplot
      )
    readr::write_tsv(
      x = res$genome.stats, 
      path = "individuals.pairwise.genome.stats.tsv", 
      col_names = TRUE
    )
    
    # Visualization ------------------------------------------------------------
    message("Generating the plots")
    
    # violin plot
    res$violin.plot.genome <- ggplot2::ggplot(
      data = pairwise.genome.similarity, 
      ggplot2::aes(x = PAIRWISE, y = PROP_IDENTICAL, na.rm = TRUE)
    ) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = "black") +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      ggplot2::labs(y = "Genome similarity (proportion)") +
      ggplot2::labs(x = "Pairwise comparison") +
      ggplot2::theme(
        # legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        # panel.grid.major.y = element_blank(), 
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
      )
    # Manhattan plot
    res$manhattan.plot.genome <- ggplot2::ggplot(
      data = pairwise.genome.similarity, 
      ggplot2::aes(x = PAIRWISE, y = PROP_IDENTICAL, colour = POP_COMP, size = MARKERS_COMMON)
    ) + 
      ggplot2::geom_jitter(alpha = 0.3) +
      ggplot2::labs(y = "Genome similarity (proportion)") +
      ggplot2::labs(x = "Pairwise comparisons") +
      ggplot2::labs(colour = "Population comparisons") +
      ggplot2::scale_colour_manual(values = c("#0571b0", "black")) +
      ggplot2::scale_size_area(name = "Markers in common", max_size = 5) +
      ggplot2::theme_light() +
      ggplot2::theme(
        # legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        # panel.grid.major.y = element_blank(), 
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
      )
  } # end genome method
  
  # RESULTS --------------------------------------------------------------------
  if (genome) message("A table: pairwise.genome.similarity.tsv was written in the working directory")
  "individuals.pairwise.dist.tsv"
  cat("################################### RESULTS ###################################\n")
  message("Object in the list (if all arguments are selected):\n
$distance                         # Distance method results
$distance.stats                   # Summary statistics of the distance method
$pairwise.genome.similarity       # Genome method results
$genome.stats                     # Summary statistics of the genome method\n\n
Visualization:
$violin.plot.distance
$manhattan.plot.distance
$violin.plot.genome
$manhattan.plot.genome\n
Saved in the working directory:
individuals.pairwise.dist.tsv
individuals.pairwise.distance.stats.tsv
individuals.pairwise.genome.similarity.tsv
individuals.pairwise.genome.stats.tsv
")
  message(stringi::stri_join("Working directory: ", getwd()))
  message(stringi::stri_join("Computation time: ", round((proc.time() - timing)[[3]]), " sec"))
  cat("############################## completed ##############################\n")
  return(res)
} # end function detect_duplicate_genomes




# Internal functions: ---------------------------------------------------------


# distance method --------------------------------------------------------------
#' @title Distance individuals
#' @description distance method
#' @rdname distance_individuals
#' @export
#' @keywords internal
distance_individuals <- function(
  x,
  strata = NULL,
  distance.method = "manhattan",
  parallel.core = parallel::detectCores() - 1
) {
  # Prep data 
  dist.computation <- suppressWarnings(
    dplyr::ungroup(x) %>%
      data.table::as.data.table() %>%
      data.table::dcast.data.table(
        data = .,
        formula = INDIVIDUALS ~ MARKERS_ALLELES, 
        value.var = "n"
      ) %>%
      tibble::as_data_frame(.) %>% 
      tibble::remove_rownames(.) %>% 
      tibble::column_to_rownames(df = ., var = "INDIVIDUALS")
  )
  # rownames(dist.computation) <- dist.computation[["INDIVIDUALS"]]
  # dist.computation[["INDIVIDUALS"]] <- NULL

  # compute distance
  # gain in speed between the 2 is very small on small data set
  if (dplyr::n_distinct(x$MARKERS_ALLELES) > 60000) {
    dist.computation <- suppressWarnings(
      amap::Dist(
        x = dist.computation, 
        method = distance.method, 
        nbproc = parallel.core
      )
    )
    
  } else {
    dist.computation <- stats::dist(
      x = dist.computation, 
      method = distance.method
    )
  }
  
  # melt the dist matrice into a data frame
  dist.computation <- stats::na.omit(
    reshape2::melt(
      data = as.matrix(dist.computation), 
      varnames = c("ID1", "ID2"), 
      value.name = "DISTANCE",
      na.rm = TRUE)[reshape2::melt(upper.tri(as.matrix(dist.computation), diag = FALSE))$value,]
  ) %>% 
    tibble::as_data_frame() %>%
    dplyr::ungroup(.) %>% 
    dplyr::mutate(DISTANCE_RELATIVE = DISTANCE/max(DISTANCE)) %>% 
    dplyr::arrange(DISTANCE)
  
  # Include population info with strata
  ID1.pop <- suppressWarnings(
    dplyr::select(.data = dist.computation, INDIVIDUALS = ID1) %>% 
      dplyr::inner_join(strata, by = "INDIVIDUALS") %>% 
      dplyr::select(ID1_POP = POP_ID))
  
  ID2.pop <- suppressWarnings(
    dplyr::select(.data = dist.computation, INDIVIDUALS = ID2) %>% 
      dplyr::inner_join(strata, by = "INDIVIDUALS") %>% 
      dplyr::select(ID2_POP = POP_ID))
  
  dist.computation <- dplyr::bind_cols(dist.computation, ID1.pop, ID2.pop) %>% 
    dplyr::mutate(
      POP_COMP = ifelse(ID1.pop == ID2.pop, "same pop", "different pop"),
      POP_COMP = factor(POP_COMP, levels = c("same pop", "different pop"), ordered = TRUE),
      PAIRWISE = rep("pairwise", n()),
      METHOD = rep(distance.method, n())
    )
  x <- ID1.pop <- ID2.pop <- NULL
  return(dist.computation)
}#End distance_individuals


# pairwise genome similarity method---------------------------------------------
#' @title Pairwise genome similarity
#' @description genome method
#' @rdname pairwise_genome_similarity
#' @export
#' @keywords internal
pairwise_genome_similarity <- function(list.pair, id.pairwise = NULL, input.prep = NULL, ...) {
  # list.pair <- 2
  id.select <- stringi::stri_join(purrr::flatten(id.pairwise[list.pair]))
  id1 <- id.select[[1]]
  id2 <- id.select[[2]]
  
  if (tibble::has_name(input.prep, "GT_BIN")) {
    # filtered dataset for the 2 ind.
    input.select <- dplyr::filter(
      .data = input.prep,
      !is.na(GT_BIN) & INDIVIDUALS %in% id.select
    ) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT_BIN)
    
    # genotypes & markers
    id1.data <- dplyr::filter(.data = input.select, INDIVIDUALS %in% id1) %>%
      dplyr::select(-INDIVIDUALS)
    
    id2.data <- dplyr::filter(.data = input.select, INDIVIDUALS %in% id2) %>%
      dplyr::select(-INDIVIDUALS)
    
    # markers
    id1.markers <- dplyr::select(id1.data, MARKERS)
    
    id2.markers <- dplyr::select(id2.data, MARKERS)
    
    # output comparison
    genome.comparison <- tibble::data_frame(
      ID1 = id1,
      ID2 = id2,
      MARKERS_COMMON = nrow(dplyr::intersect(x = id1.markers, y = id2.markers)),
      IDENTICAL = nrow(dplyr::intersect(x = id1.data, y = id2.data)),
      DIFFERENT = MARKERS_COMMON - IDENTICAL,
      PROP_IDENTICAL = IDENTICAL / MARKERS_COMMON
    )
    #unused objets:
    input.select <- id1.data <- id2.data <- id1.markers <- id2.markers <- NULL
  } else {
    input.select <- dplyr::filter(
      .data = input.prep, 
      INDIVIDUALS %in% id.select & n != 0
      ) %>% 
      dplyr::mutate(
        MARKERS = stringi::stri_sub(str = MARKERS_ALLELES, from = 1, to = -5)
      )

    # genotypes & markers
    id1.data <- dplyr::filter(.data = input.select, INDIVIDUALS %in% id1) %>%
      dplyr::select(-INDIVIDUALS)
    
    id2.data <- dplyr::filter(.data = input.select, INDIVIDUALS %in% id2) %>%
      dplyr::select(-INDIVIDUALS)

    # markers
    id1.markers <- dplyr::distinct(id1.data, MARKERS)
    id2.markers <- dplyr::distinct(id2.data, MARKERS)
    
    # output comparison
    genome.comparison <- tibble::data_frame(
      ID1 = id1,
      ID2 = id2,
      MARKERS_COMMON = nrow(dplyr::intersect(x = id1.markers, y = id2.markers)),
      IDENTICAL = nrow(dplyr::intersect(x = id1.data, y = id2.data) %>% dplyr::distinct(MARKERS)),
      DIFFERENT = MARKERS_COMMON - IDENTICAL,
      PROP_IDENTICAL = IDENTICAL / MARKERS_COMMON
    )
    #unused objets:
    input.select <- id1.data <- id2.data <- id1.markers <- id2.markers <- NULL
  }
  return(genome.comparison)
} # end duplicate pairwise
