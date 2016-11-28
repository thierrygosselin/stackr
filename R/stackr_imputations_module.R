# stackr impute module

#' @name stackr_imputations_module
#' @title Map-independent imputation of missing genotype
#' 
#' @description Map-independent imputation of missing genotype using Random Forest 
#' or the most frequent category. Impute genotypes or alleles. Used internally 
#' in \href{https://github.com/thierrygosselin/assigner}{assigner} and
#' \href{https://github.com/thierrygosselin/stackr}{stackr} and 
#' might be of interest for users. 
#' 
#' The goal of \code{missing_visualization} is to provide a simple solution for a complicated
#' problem: missing observations in genotype datasets. Before running this function
#' to populate the original dataset with synthetic data you should try 
#' \code{\link[stackr]{missing_visualization}} to detect patterns of missingness.
#' Follow the vignette for more info. 

#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. To import, the function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{tidy_wide}}. See details for more info.

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
#' @importFrom parallel detectCores
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate mutate_all summarise left_join inner_join right_join anti_join semi_join full_join summarise_each_ funs bind_rows

#' @importFrom tidyr spread gather separate
#' @importFrom purrr map flatten keep
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom randomForestSRC impute.rfsrc
#' @importFrom stringi stri_replace_all_fixed stri_replace_na stri_join
#' @importFrom tibble has_name


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

stackr_imputations_module <- function(
  data, 
  imputation.method = NULL,
  impute = "genotype",
  imputations.group = "populations",
  num.tree = 100,
  iteration.rf = 10,
  split.number = 100,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  filename = NULL,
  ...
) {
  
  # for timing
  timing <- proc.time()
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  
  # Import data ---------------------------------------------------------------
  input <- stackr::tidy_wide(data = data, import.metadata = TRUE)
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (!tibble::has_name(input, "MARKERS") && tibble::has_name(input, "LOCUS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # scan for the columnn CHROM and keep the info to include it after imputations
  # detect if data is biallelic
  if (tibble::has_name(input, "CHROM")) {
    marker.meta <- dplyr::distinct(.data = input, MARKERS, CHROM, LOCUS, POS)
  } else {
    marker.meta <- NULL
  }
  
  # scan for REF allele column
  if (tibble::has_name(input, "REF")) {
    ref.column <- TRUE
    biallelic <- TRUE
  } else {
    ref.column <- FALSE
    biallelic <- detect_biallelic_markers(data = input)
  }
  
  # select the column for the imputations
  input <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT)
  
  # output the proportion of missing genotypes BEFORE imputations
  still.na <- dplyr::summarise(.data = input, MISSING = round(length(GT[GT == "000000"])/length(GT), 3))
  message(stringi::stri_join("Proportion of missing genotypes before imputations: ", still.na$MISSING))
  
  # Imputations ----------------------------------------------------------------
  message("Preparing the data for imputations")
  
  # prepare the data in working with genotype or alleles
  if (impute == "genotype") {
    input.prep <- input %>%
      dplyr::mutate(GT = replace(GT, which(GT == "000000"), NA)) %>%
      dplyr::group_by(INDIVIDUALS, POP_ID) %>% 
      tidyr::spread(data = ., key = MARKERS, value = GT) %>%
      dplyr::ungroup(.) %>% 
      dplyr::arrange(POP_ID, INDIVIDUALS)
  }
  if (impute == "allele") {
    input.prep <- input %>%
      tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
      tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%  
      dplyr::mutate(GT = replace(GT, which(GT == "000"), NA)) %>%
      dplyr::group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
      tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
      dplyr::ungroup(.) %>% 
      dplyr::arrange(POP_ID, INDIVIDUALS)
  }
  
  # keep stratification
  strata.df.impute <- dplyr::distinct(.data = input.prep, INDIVIDUALS, POP_ID)
  
  # Imputation with Random Forest  ---------------------------------------------
  if (imputation.method == "rf") {
    # Parallel computations options
    options(rf.cores = parallel.core, mc.cores = parallel.core)
    
    # imputations using Random Forest with the package randomForestSRC
    impute_genotype_rf <- function(x) {
      randomForestSRC::impute.rfsrc(
        data = data.frame(x),
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
      imputed.dataset <- list() # create empty list
      
      # Function to go through the populations
      impute_rf_pop <- function(pop.list, ...){
        sep.pop <- df.split.pop[[pop.list]]
        sep.pop <- dplyr::mutate_all(.tbl = sep.pop, .funs = factor, exclude = NA)
        imputed.dataset <- impute_genotype_rf(sep.pop)
        imputed.dataset <- dplyr::mutate_all(.tbl = imputed.dataset, .funs = as.character, exclude = NA)
        # message of progress for imputations by population
        message(paste("Completed imputations for pop ", pop.list, sep = ""))
        return(imputed.dataset)
      } # End impute_rf_pop
      
      input.imp <- list()
      input.imp <- .stackr_parallel(
        X = pop.list,
        FUN = impute_rf_pop,
        mc.preschedule = FALSE,
        mc.silent = FALSE,
        mc.cleanup = TRUE,
        mc.cores = parallel.core
      )
      
      # Compiling the results
      message("Compiling imputations results")
      input.imp <- suppressWarnings(dplyr::bind_rows(input.imp))
      
      # reintroduce the stratification
      input.imp <- suppressWarnings(
        dplyr::left_join(
          strata.df.impute,
          dplyr::select(.data = input.imp, -POP_ID)
          , by = "INDIVIDUALS"
        ) %>% 
          dplyr::arrange(POP_ID, INDIVIDUALS) %>% 
          dplyr::ungroup(.)
      )
      
      # remove unused objects
      df.split.pop <- pop.list <- sep.pop <- imputed.dataset <- input.prep <- NULL
    } # End imputation RF populations
    
    # Random Forest global
    if (imputations.group == "global") { # Globally (not by pop_id)
      message("Imputations computed globally, take a break...")
      input.prep <- dplyr::mutate_all(.tbl = input.prep, .funs = factor, exclude = NA)
      input.imp <- impute_genotype_rf(input.prep)
      input.imp <- dplyr::mutate_all(.tbl = input.imp, .funs = as.character, exclude = NA)
      
      # reintroduce the stratification
      input.imp <- suppressWarnings(
        dplyr::left_join(strata.df.impute,
                         dplyr::select(.data = input.imp, -POP_ID)
                         , by = "INDIVIDUALS") %>% 
          dplyr::arrange(POP_ID, INDIVIDUALS) %>% 
          dplyr::ungroup(.)
      )
      input.prep <- NULL # remove unused object
    } # End imputation RF global
    
    # data prep
    if (impute == "genotype") {
      input.imp <- suppressWarnings(
        tidyr::gather(
          data = input.imp,
          key = MARKERS, value = GT,
          -c(INDIVIDUALS, POP_ID)
        )
      )
    }
    if (impute == "allele") {
      input.imp <- suppressWarnings(
        tidyr::gather(
          data = input.imp,
          key = MARKERS, value = GT,
          -c(INDIVIDUALS, POP_ID, ALLELES)
        )
      )
    }
    
  } # End imputation RF
  
  # Imputation using the most common genotype ----------------------------------
  if (imputation.method == "max") { # End imputation max
    if (imputations.group == "populations") {
      message("Imputations computed by populations")
      
      if (impute == "genotype") {
        input.imp <- suppressWarnings(
          input.prep %>%
            tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
            dplyr::group_by(MARKERS, POP_ID) %>%
            dplyr::mutate(
              GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
              GT = replace(GT, which(GT == "NA"), NA)
            ) %>%
            # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
            # will take the global observed values by markers for those cases.
            # dplyr::group_by(MARKERS) %>%
            # dplyr::mutate(GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
            dplyr::ungroup(.)
        )
      }
      if (impute == "allele") {
        input.imp <- suppressWarnings(
          input.prep %>%
            tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
            dplyr::group_by(MARKERS, POP_ID) %>%
            dplyr::mutate(
              GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
              GT = replace(GT, which(GT == "NA"), NA)
            ) %>%
            # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
            # will take the global observed values by markers for those cases.
            # dplyr::group_by(MARKERS) %>%
            # dplyr::mutate(GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
            dplyr::ungroup(.)
        )
      }
      input.prep <- NULL # remove unused object
    } # End imputation max populations 
    if (imputations.group == "global") {
      # Globally (not by pop_id)
      message("Imputations computed globally")
      if (impute == "genotype") {
        input.imp <- suppressWarnings(
          input.prep %>%
            tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
            dplyr::ungroup(.)
        )
      }
      
      if (impute == "allele") {
        input.imp <- suppressWarnings(
          input.prep %>%
            tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
            dplyr::ungroup(.)
        )
      }
      
      input.prep <- NULL # remove unused object
    } # End imputation max global 
  } # End imputations max
  
  # results --------------------------------------------------------------------
  message("Tidying the imputed data set")
  
  # back to a tidy format for the allele dataset
  if (impute == "allele") {
    input.imp <- input.imp %>% 
      tidyr::spread(data = ., key = ALLELES, value = GT) %>% 
      tidyr::unite(data = ., GT, A1, A2, sep = "")
  }
  
  # Replace NA by 000000 in GT column
  input.imp$GT <- stringi::stri_replace_na(str = input.imp$GT, replacement = "000000")
  
  
  # Compute REF/ALT allele... might have change depending on prop of missing values
  if (ref.column) {
    message("Adjusting REF/ALT alleles to account for imputations...")
    input.temp <- stackr::change_alleles(data = input.imp)
    input.imp <- input.temp$input
    
  } # end computing REF/ALT
  
  # Integrate marker.meta columns
  if (!is.null(marker.meta)) {
    input.imp <- dplyr::left_join(input.imp, marker.meta, by = "MARKERS") #%>% dplyr::select(MARKERS, CHROM, LOCUS, POS, REF, ALT, POP_ID, INDIVIDUALS, GT, GT_VCF, GT_BIN)
  }
  
  # Write to working directory
  if (!is.null(filename)) {
    message(stringi::stri_join("Writing the imputed tidy data to the working directory: \n"), filename)
    readr::write_tsv(x = input.imp, path = filename, col_names = TRUE)
  }
  
  # output the proportion of missing genotypes after imputations
  still.na <- dplyr::summarise(.data = input.imp, MISSING = round(length(GT[GT == "000000"])/length(GT), 3))
  message(stringi::stri_join("Proportion of missing genotypes after imputations: ", still.na$MISSING))
  
  # for timing
  timing <- proc.time() - timing
  message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
  return(input.imp)
} # End imputations
