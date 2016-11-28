# Detect duplicate genomes

#' @name detect_duplicate_genomes
#' @title Compute pairwise genome similarity or distance between individuals 
#' to highligh potential duplicate individuals
#' @description The function can compute two methods 
#' to highligh potential duplicate individuals (1. distance between individuals
#' or 2. pairwise genome similarity).

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{stackr} \code{\link{tidy_genomic_data}}.

#' @param distance.method (character) The distance measure used inside \code{stats::dist} 
#' (<= 30000 markers) or \code{amap::Dist} (> 30000 markers). 
#' This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary".
#' Using \code{distance.method = NULL} will not run this method.
#' Default: \code{distance.method = "manhattan"}. This is very fast (max 5 min)
#' compared to the genome similarity method. 

#' @param genome (logical) Computes pairwise genome similarity in parallel.
#' The proportion of the shared genotypes is averaged across shared markers between 
#' each pairwise comparison. This method makes filtering easier because the 
#' threshold is more intuitive with the plots produced, but it's much longer 
#' to run, even in parallel, so better to run overnight. 
#' Default: \code{genome = FALSE}.

#' @param parallel.core (optional) The number of core for parallel computation.
#' Default: \code{parallel.core = detectCores()-1}.

#' @return A list with potentially 8 objects: 
#' \code{$distance }: results of the distance method
#' \code{$distance.stats}: Summary statistics of the distance method
#' \code{$pairwise.genome.similarity}: results of the genome method
#' \code{$genome.stats}: Summary statistics of the genome method
#' \code{$violin.plot.distance}: violin plot showing the distribution of pairwise distances
#' \code{$jitter.plot.distance}: same info different visual with jitter plot
#' \code{$violin.plot.genome}: violin plot showing the distribution of pairwise genome similarities
#' \code{$jitter.plot.genome}: same info different visual with jitter plot
#' 
#' Saved in the working directory:
#' individuals.pairwise.dist.tsv, individuals.pairwise.distance.stats.tsv, 
#' individuals.pairwise.genome.similarity.tsv, individuals.pairwise.genome.stats.tsv

#' @export
#' @rdname detect_duplicate_genomes
#' @importFrom stringi stri_paste stri_replace_all_fixed
#' @importFrom dplyr rename select group_by filter mutate rename_ filter_ bind_cols bind_rows summarise
#' @importFrom utils combn
#' @importFrom stats na.omit
#' @importFrom amap Dist
#' @importFrom data.table fread
#' @importFrom reshape2 melt
#' @importFrom lazyeval interp
#' @importFrom readr write_tsv
#' @importFrom parallel detectCores
#' @importFrom purrr flatten map

#' @examples
#' \dontrun{
#' # to get pairwise distance only, the simplest way to run:
#' dup <- detect_duplicate_genomes(data = "wombat_tidy.tsv")
#' # This will use by defaul \code{distance.method = "manhattan"}, 
#' \code{genome = FALSE}, and all my CPU -1 as default for \code{parallel.core}
#' 
#' # To view the jitter plot:
#' dup$jitter.plot.distance
#' 
#' # to view the data stats
#' dup.data.stats <- dup$distance.stats
#' 
#' # to view the data
#' dup.data <- dup$distance
#' 
#' # Based on the look of the distribution using both jitter and boxplot, 
#' I can filter the dataset to highlight potential duplicates: 
#' dup.filtered <- filter(.data = dup.data, DISTANCE < 3000000)
#' 
#' # To run the distance (with euclidean distance instead of the default manhattan, 
#' # with the genome methods:
#' dup <- detect_duplicate_genomes(
#' data = "wombat_tidy.tsv", 
#' distance.method = "euclidean",
#' genome = TRUE
#' )
#' # to view the data of the genome data
#' dup.data <- dup$pairwise.genome.similarity
#' 
#' # Based on the look of the distribution using both jitter and boxplot, 
#' # I can filter the dataset based on 98% of identical genotype proportion, 
#' # to highlight potential duplicates: 
#' dup.filtered <- filter(.data = dup.data, PROP_IDENTICAL > 0.98)
#' 
#' # Get the list of duplicates id
#' dup.list.names <- data.frame(INDIVIDUALS = unique(c(dup.filtered$ID1, dup.filtered$ID2)))
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_duplicate_genomes <- function(
  data,
  distance.method = "manhattan",
  genome = FALSE,
  parallel.core = detectCores() - 1) {
  
  cat("###############################################################################\n")
  cat("######################## stackr: detect_duplicate_genomes ########################\n")
  cat("###############################################################################\n")
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
  
  # strata
  strata <- input %>% 
    ungroup %>% 
    distinct(POP_ID, INDIVIDUALS)
  
  # Functions -------------------------------------------------------------------
  # distance method
  distance_individuals <- function(x, distance.method, parallel.core) {
    if (!requireNamespace("tibble")) warning("tibble not installed")
    
    # Prep data 
    dist.computation <- suppressWarnings(
      x %>% 
        dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
        dplyr::group_by(INDIVIDUALS) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
        ungroup %>% 
        tibble::remove_rownames(.) %>% 
        tibble::column_to_rownames(df = ., var = "INDIVIDUALS")
    )    
    
    # compute distance
    # gain in speed between the 2 is very small on small data set
    if (n_distinct(x$MARKERS) > 30000) {
      if (requireNamespace("amap")) {
        dist.computation <- suppressWarnings(amap::Dist(
          x = dist.computation, 
          method = distance.method, 
          nbproc = parallel.core
        )
        )
      } else {
        warning("amap not installed, using stats::dist instead")
        dist.computation <- stats::dist(
          x = dist.computation, 
          method = distance.method
        )
      }
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
      as_data_frame() %>% 
      arrange(DISTANCE)
    
    # Include population info with strata
    ID1.pop <- suppressWarnings(
      dist.computation %>% 
        dplyr::select(INDIVIDUALS = ID1) %>% 
        dplyr::inner_join(strata, by = "INDIVIDUALS") %>% 
        dplyr::select(ID1_POP = POP_ID))
    
    ID2.pop <- suppressWarnings(
      dist.computation %>% 
        dplyr::select(INDIVIDUALS = ID2) %>% 
        dplyr::inner_join(strata, by = "INDIVIDUALS") %>% 
        dplyr::select(ID2_POP = POP_ID))
    
    dist.computation <- dplyr::bind_cols(dist.computation, ID1.pop, ID2.pop) %>% 
      dplyr::mutate(
        POP_COMP = ifelse(ID1.pop == ID2.pop, "same pop", "different pop"),
        POP_COMP = factor(POP_COMP, levels = c("same pop", "different pop"), ordered = TRUE),
        PAIRWISE = rep("pairwise", n()),
        METHOD = rep(distance.method, n())
      )
    return(dist.computation)
  }
  
  # pairwise genome similarity method
  pairwise_genome_similarity <- function(list.pair, ...) {
    # list.pair <- 2
    id.select <- stringi::stri_join(purrr::flatten(id.pairwise[list.pair]))
    
    id1 <- id.select[[1]]
    id2 <- id.select[[2]]
    
    identical.gt <- input %>%
      dplyr::filter_(lazyeval::interp(~ INDIVIDUALS == as.name(id1) | INDIVIDUALS == as.name(id2))) %>% 
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>% 
      tidyr::spread(data = ., INDIVIDUALS, GT) %>%
      dplyr::rename_(ID1 = as.name(id1), ID2 = as.name(id2)) %>%
      dplyr::mutate(
        ID1 = stringi::stri_replace_all_fixed(
          str = as.character(ID1), 
          pattern = c("/", ":", "_", "-", "."), 
          replacement = "", 
          vectorize_all = FALSE),
        ID2 = stringi::stri_replace_all_fixed(
          str = as.character(ID2), 
          pattern = c("/", ":", "_", "-", "."), 
          replacement = "", 
          vectorize_all = FALSE)
      ) %>% 
      dplyr::filter(ID1 != "000000" | ID2 != "000000") %>% 
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        IDENTICAL_GT = ifelse(ID1 == ID2, "IDENTICAL", "DIFFERENT")
      ) %>%
      dplyr::group_by(IDENTICAL_GT) %>% 
      tally %>% 
      tidyr::spread(data = ., IDENTICAL_GT, n) %>% 
      dplyr::mutate(
        TOTAL_MARKERS_GENOTYPED = IDENTICAL + DIFFERENT,
        PROP_IDENTICAL = IDENTICAL/TOTAL_MARKERS_GENOTYPED,
        ID1 = id1,
        ID2 = id2
      ) %>% 
      dplyr::select(ID1, ID2, IDENTICAL, DIFFERENT, TOTAL_MARKERS_GENOTYPED, PROP_IDENTICAL)
    return(identical.gt)
  } # end duplicate pairwise
  
  # New list to prepare for results
  res <- list()
  
  # Compute distance between individuals --------------------------------------
  if (!is.null(distance.method)) {
    message(stringi::stri_join("Computing ", distance.method, " distances between individuals..."))
    # distance.method <- "euclidean"
    # distance.method <- "manhattan"
    # parallel.core <- 8
    
    dist.computation <- distance_individuals(
      x = input, 
      distance.method = distance.method, 
      parallel.core = parallel.core
    )
    res$distance <- dist.computation
    readr::write_tsv(
      x = dist.computation,
      path = "individuals.pairwise.dist.tsv",
      col_names = TRUE
    )
    
    
    # Stats
    message("Generating summary statistics")
    distance.stats <- dist.computation %>% 
      dplyr::summarise(
        MEAN = mean(DISTANCE, na.rm = TRUE),
        MEDIAN = stats::median(DISTANCE, na.rm = TRUE),
        SE = round(sqrt(stats::var(DISTANCE, na.rm = TRUE)/length(stats::na.omit(DISTANCE))), 2),
        MIN = round(min(DISTANCE, na.rm = TRUE), 2),
        MAX = round(max(DISTANCE, na.rm = TRUE), 2),
        QUANTILE25 = stats::quantile(DISTANCE, 0.25), # quantile25
        QUANTILE75 = stats::quantile(DISTANCE, 0.75)#, # quantile75
        # OUTLIERS_LOW = QUANTILE25 - (1.5 * (QUANTILE75 - QUANTILE25)), # outliers : below the outlier boxplot
        # OUTLIERS_HIGH = QUANTILE75 + (1.5 * (QUANTILE75 - QUANTILE25)) # outliers : higher the outlier boxplot
      )
    res$distance.stats <- distance.stats
    readr::write_tsv(
      x = distance.stats, 
      path = "individuals.pairwise.distance.stats.tsv", 
      col_names = TRUE
    )
    
    message("Generating the plots")
    # violin plot
    res$violin.plot.distance <- ggplot(
      data = dist.computation, 
      aes(x = PAIRWISE, y = DISTANCE, na.rm = TRUE)
    ) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = "black") +
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      labs(y = "Distance\n <- close     distant ->") +
      labs(x = "Pairwise comparisons") +
      theme(
        # legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = element_text(size = 8, family = "Helvetica")
      )
    
    # Jitter plot
    res$jitter.plot.distance <- ggplot(data = dist.computation, aes(x = PAIRWISE, y = DISTANCE, colour = POP_COMP)) + 
      geom_jitter(alpha = 0.3) + 
      labs(y = "Distance\n <- distant     close ->") +
      labs(x = "Pairwise comparisons") +
      labs(colour = "Population comparisons") +
      scale_colour_manual(values = c("#0571b0", "black")) +
      scale_y_reverse() +
      theme_light() +
      theme(
        # legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = element_text(size = 8, family = "Helvetica")
      )
  } # end distance method
  
  # Compute pairwise search for duplicate --------------------------------------
  if (genome) {
    # list of id
    id.list <- unique(input$INDIVIDUALS) # id list
    
    # all combination of individual pair
    id.pairwise <- utils::combn(unique(id.list), 2, simplify = FALSE) 
    list.pair <- 1:length(id.pairwise)
    # list.pair <- 5
    # parallel.core<-8
    message("Starting scan for duplicate genomes, take a break...")
    pairwise.genome.similarity <- .stackr_parallel(
      X = list.pair, 
      FUN = pairwise_genome_similarity, 
      mc.preschedule = FALSE, 
      mc.silent = FALSE, 
      mc.cleanup = TRUE,
      mc.cores = parallel.core
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
    genome.stats <- pairwise.genome.similarity %>% 
      dplyr::summarise(
        MEAN = mean(PROP_IDENTICAL, na.rm = TRUE),
        MEDIAN = stats::median(PROP_IDENTICAL, na.rm = TRUE),
        SE = round(sqrt(stats::var(PROP_IDENTICAL, na.rm = TRUE)/length(na.omit(PROP_IDENTICAL))), 2),
        MIN = round(min(PROP_IDENTICAL, na.rm = TRUE), 2),
        MAX = round(max(PROP_IDENTICAL, na.rm = TRUE), 2),
        QUANTILE25 = stats::quantile(PROP_IDENTICAL, 0.25), # quantile25
        QUANTILE75 = stats::quantile(PROP_IDENTICAL, 0.75)#, # quantile75
        # OUTLIERS_LOW = QUANTILE25 - (1.5 * (QUANTILE75 - QUANTILE25)), # outliers : below the outlier boxplot
        # OUTLIERS_HIGH = QUANTILE75 + (1.5 * (QUANTILE75 - QUANTILE25)) # outliers : higher the outlier boxplot
      )
    res$genome.stats <- genome.stats
    readr::write_tsv(
      x = genome.stats, 
      path = "individuals.pairwise.genome.stats.tsv", 
      col_names = TRUE
    )
    
    # Visualization ------------------------------------------------------------
    message("Generating the plots")
    
    # violin plot
    res$violin.plot.genome <- ggplot(
      data = pairwise.genome.similarity, 
      aes(x = PAIRWISE, y = PROP_IDENTICAL, na.rm = TRUE)
    ) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = "black") +
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      labs(y = "Genome similarity (proportion)") +
      labs(x = "Pairwise comparison") +
      theme(
        # legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = element_text(size = 8, family = "Helvetica")
      )
    # jitter plot
    res$jitter.plot.genome <- ggplot(
      data = pairwise.genome.similarity, 
      aes(x = PAIRWISE, y = PROP_IDENTICAL, colour = POP_COMP)
    ) + 
      geom_jitter(alpha = 0.3) + 
      labs(y = "Genome similarity (proportion)") +
      labs(x = "Pairwise comparisons") +
      labs(colour = "Population comparisons") +
      scale_colour_manual(values = c("#0571b0", "black")) +
      theme_light() +
      theme(
        # legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = element_text(size = 8, family = "Helvetica")
      )
  } # end genome method
  
  # RESULTS --------------------------------------------------------------------
  if (genome) message("A table: pairwise.genome.similarity.tsv was written in the working directory")
  "individuals.pairwise.dist.tsv"
  cat("################################### RESULTS ###################################\n")
  message("Object in the list (if all arguments are selected):\n
$distance                         # results of the distance method
$distance.stats                   # Summary statistics of the distance method
$pairwise.genome.similarity       # results of the genome method
$genome.stats                     # Summary statistics of the genome method\n
Visualization:
$violin.plot.distance
$jitter.plot.distance
$violin.plot.genome
$jitter.plot.genome\n
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
