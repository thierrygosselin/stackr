# Find duplicate individual

#' @name find_duplicate_id
#' @title Compute pairwise genome similarity to highligh potential duplicate 
#' individuals
#' @description The function computes pairwise genome similarity 
#' to highligh potential duplicate individuals.

#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. \code{\link{tidy_genomic_data}} can
#' transform numerous genomic data format in tidy data frames.

#' @param parallel.core (optional) The number of core for parallel computation 
#' of pairwise Fst. 
#' If not selected \code{detectCores()-1} is used as default.

#' @return A list with 2 objects: 
#' \code{$duplicate.info.summary}: a data frame with the summary of the similarity analysis.
#' \code{$similarity.plot}: boxplot of the similarity analysis.

#' @export
#' @rdname find_duplicate_id
#' @import stringi
#' @import dplyr
#' @import utils
#' @importFrom purrr map
#' @importFrom data.table fread

#' @examples
#' \dontrun{
#' dup <- find_duplicate_id(data = "wombat_tidy.tsv", parallel.core = 8)
#' 
#' # To view the data frame:
#' dup.df <- dup$duplicate.info.summary
#' 
#' # To view the boxplot:
#' bp <- dup$similarity.plot
#' 
#' # Based on the look of the distribution using the boxplot, I choose to filter 
#' # the resuls with a threshold of 90%:
#' dup.filtered <-  %>% filter(PROP_IDENTICAL > 0.90)
#' 
#' # Get the list of duplicates id
#' dup.list.names <- data.frame(INDIVIDUALS = unique(c(dup.filtered$ID1, dup.filtered$ID2)))
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1"){
  utils::globalVariables(
    c("ID1", "ID2", "IDENTICAL_GT", "IDENTICAL", "DIFFERENT", 
      "TOTAL_MARKERS_GENOTYPED", "PROP_IDENTICAL", "PAIRWISE")
  )
}

find_duplicate_id <- function(data, parallel.core = detectCores()-1) {
  cat("#######################################################################\n")
  cat("###################### stackr: find_duplicate_id ######################/n")
  cat("#######################################################################\n")
  
  # manage missing arguments ---------------------------------------------------
  if (missing(data)) stop("missing data argument")
  if (missing(parallel.core) | is.null(parallel.core)) parallel.core <- detectCores()-1
  
  # Import data ----------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data)
  
  
  # Compute pairwise search for duplicate ----------------------------------------
  # list of id
  id.list <- unique(input$INDIVIDUALS) # id list
  
  # all combination of individual pair
  id.pairwise <- combn(unique(id.list), 2, simplify = FALSE) 
  list.pair <- 1:length(id.pairwise)
  
  # Fst for all pairwise populations
  duplicate_pairwise <- function(list.pair, ...) {
    # list.pair <- 2
    id.select <- stri_paste(flatten(id.pairwise[list.pair]))
    
    id1 <- id.select[[1]]
    id2 <- id.select[[2]]
    
    identical.gt <- input %>%
      filter_(lazyeval::interp(~ INDIVIDUALS == as.name(id1) | INDIVIDUALS == as.name(id2))) %>% 
      select(MARKERS, INDIVIDUALS, GT) %>% 
      tidyr::spread(data = ., INDIVIDUALS, GT) %>%
      rename_(ID1 = as.name(id1), ID2 = as.name(id2)) %>%
      mutate(
        ID1 = stri_replace_all_fixed(
          str = as.character(ID1), 
          pattern = c("/", ":", "_", "-", "."), 
          replacement = "", 
          vectorize_all = FALSE),
        ID2 = stri_replace_all_fixed(
          str = as.character(ID2), 
          pattern = c("/", ":", "_", "-", "."), 
          replacement = "", 
          vectorize_all = FALSE)
      ) %>% 
      filter(ID1 != "000000" | ID2 != "000000") %>% 
      group_by(MARKERS) %>%
      mutate(
        IDENTICAL_GT = ifelse(ID1 == ID2, "IDENTICAL", "DIFFERENT")
      ) %>%
      group_by(IDENTICAL_GT) %>% 
      tally %>% 
      tidyr::spread(data = ., IDENTICAL_GT, n) %>% 
      mutate(
        TOTAL_MARKERS_GENOTYPED = IDENTICAL + DIFFERENT,
        PROP_IDENTICAL = IDENTICAL/TOTAL_MARKERS_GENOTYPED,
        ID1 = id1,
        ID2 = id2
      ) %>% 
      select(ID1, ID2, IDENTICAL, DIFFERENT, TOTAL_MARKERS_GENOTYPED, PROP_IDENTICAL)
    return(identical.gt)
  } # end duplicate pairwise
  
  # list.pair <- 5
  # parallel.core<-8
  message("Starting scan for duplicate genome, take a break...")
  duplicate.info <- parallel::mclapply(
    X = list.pair, 
    FUN = duplicate_pairwise, 
    mc.preschedule = FALSE, 
    mc.silent = FALSE, 
    mc.cores = parallel.core
  )
  duplicate.info.summary <- bind_rows(duplicate.info)
  
  res <- list()
  res$duplicate.info.summary
  
  write_tsv(x = duplicate.info.summary, path = "duplicate.info.summary.tsv", col_names = TRUE)
  
  # Visualization --------------------------------------------------------------
  duplicate.info.summary.mod <- duplicate.info.summary %>% 
    mutate(PAIRWISE = rep("pairwise comparison", n()))
  
  similarity.plot <- ggplot(duplicate.info.summary.mod, aes(x = PAIRWISE, y = PROP_IDENTICAL, na.rm = TRUE))+
    # geom_violin(trim = TRUE)+
    geom_boxplot(width = 0.5)+
    stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
    labs(x = "Pairwise comparison")+
    labs(y = "Genome similarity proportion")+
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = element_blank(), 
      # axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
  
  res$similarity.plot <- similarity.plot
  return(res)
} # end function find_duplicate_id
