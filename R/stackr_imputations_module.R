# stackr impute module

#' @name stackr_imputations_module
#' @title Map-independent imputation of missing genotype
#' 
#' @description Map-independent imputation of missing genotype using Random Forest 
#' or the most frequent category. Impute genotypes or alleles. Used internally 
#' in \href{https://github.com/thierrygosselin/assigner}{assigner} and
#' \href{https://github.com/thierrygosselin/stackr}{stackr} and 
#' might be of interest for users.

#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. To import, the function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}}. See details for more info.

#' @param imputation.method (character, optional) 
#' Methods available for map-independent imputations of missing genotype: 
#' (1) \code{"max"} to use the most frequent category for imputations.
#' (2) \code{"rf"} using Random Forest algorithm. 
#' Default: no imputation \code{imputation.method = NULL}.

#' @param impute (character, optional) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
#' Default: \code{"genotype"}.

#' @param imputations.group (character, optional) \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default = \code{"populations"}.

#' @param num.tree (integer, optional) The number of trees to grow in Random Forest. 
#' Default: \code{num.tree = 100}.

#' @param iteration.rf (integer, optional) The number of iterations of missing data algorithm
#' in Random Forest. 
#' Default: \code{iteration.rf = 10}.

#' @param split.number (integer, optional) Non-negative integer value used to specify
#' random splitting in Random Forest. 
#' Default: \code{split.number = 100}.

#' @param verbose (logical, optional) Should trace output be enabled on each iteration
#' in Random Forest ? 
#' Default: \code{verbose = FALSE}.

#' @param parallel.core (optional) The number of core for OpenMP shared-memory parallel
#' programming of Random Forest imputations. For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.

#' @param filename (optional) The file name for the \strong{imputed} 
#' tidy data frame written to the working directory.
#' With missing argument or default: \code{filename = NULL}, the imputed 
#' tidy data is in the global environment only 
#' (i.e. not written in the working directory).

#' @param ... other parameters passed to the function.


#' @return The output in your global environment is the imputed tidy data frame.
#' If \code{filename} is provided, the imputed tidy data frame is also 
#' written to the working directory.


#' @export
#' @rdname stackr_imputations_module
#' @import parallel
#' @import stringi
#' @import dplyr
#' @importFrom purrr map
#' @importFrom purrr flatten
#' @importFrom purrr keep
#' @importFrom data.table fread
#' @importFrom data.table melt.data.table
#' @importFrom data.table as.data.table

#' @examples
#' \dontrun{
#' # The simplest form:
#' wolf.imputed.data <- stackr_imputations_module(data = "wolf.tidy.dataset.tsv")
#' This will impute missing genotype by populations using random forest.
#' The number of tree, iterations.rf, split.number and parallel.core used 
#' will be the default.
#' }

#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# # required to pass the R CMD check and have 'no visible binding for global variable'
# if (getRversion() >= "2.15.1") {
#   utils::globalVariables(
#     c("DP", "AD", "vcf.headers", "GT_VCF", "INDIVIDUALS2", ""
#     )
#   )
# }

stackr_imputations_module <- function(
  data, 
  imputation.method = NULL,
  impute = "genotype",
  imputations.group = "populations",
  num.tree = 100,
  iteration.rf = 10,
  split.number = 100,
  verbose = FALSE,
  parallel.core = detectCores()-1,
  filename = NULL,
  ...) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
  if (!"MARKERS" %in% colnames(input) & "LOCUS" %in% colnames(input)) {
    input <- rename(.data = input, MARKERS = LOCUS)
  }
  
  input <- select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT)
  
  # Imputations ***************************************************************
  message("Preparing the data for imputations")
  
  if (impute == "genotype") {
    input.prep <- input %>%
      mutate(
        GT = stri_replace_all_fixed(GT, pattern = "000000", replacement = "NA", vectorize_all = FALSE),
        GT = replace(GT, which(GT == "NA"), NA)
      ) %>%
      group_by(INDIVIDUALS, POP_ID) %>% 
      tidyr::spread(data = ., key = MARKERS, value = GT) %>%
      ungroup() %>% 
      arrange(POP_ID, INDIVIDUALS)
  }
  if (impute == "allele") {
    input.prep <- input %>%
      tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
      tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%  
      mutate(
        GT = stri_replace_all_fixed(GT, pattern = "000", replacement = "NA", vectorize_all = FALSE),
        GT = replace(GT, which(GT == "NA"), NA)
      ) %>%
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
        mc.cleanup = TRUE,
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

  # prepare the imputed dataset for gsi_sim or adegenet
  message("Preparing imputed data set...")
  if (impute == "allele") {
    input.imp <- input.imp %>% 
      tidyr::spread(data = ., key = ALLELES, value = GT) %>% 
      tidyr::unite(data = ., GT, A1, A2, sep = "")
  }
  
  # results --------------------------------------------------------------------
  # Write to working directory
  if (!is.null(filename)) {
  message(stri_paste("Writing the imputed tidy data to the working directory: \n"), filename)
  write_tsv(x = input.imp, path = filename, col_names = TRUE)
  }
  return(input.imp)
} # End imputations
