# stackr impute module

#' @name stackr_imputations_module
#' @title Map-independent imputation of missing genotype
#' 
#' @description Used internally in \href{https://github.com/thierrygosselin/assigner}{assigner} and
#' \href{https://github.com/thierrygosselin/stackr}{stackr} and 
#' might be of interest for users.
#' The goal of this module is to provide a simple solution for
#' a complicated problem: missing observations in RADseq genomic datasets.
#' This function will compute \strong{map-independent imputations} of missing
#' genotypes.
#' 
#' \strong{Key features:}
#' 
#' \itemize{
#' \item Random Forests (rf) or the most observed genotypes (~mean/mode).
#' \item Imputations conducted by populations or globally.
#' \item Haplotype/SNP approach: correlation among SNPs is accounted for during
#' rf imputation, i.e. imputation is automatically conducted by haplotype
#' when marker meta-information is avaialble (chromosome, locus and position,
#' usually from VCF files). The alternative, considers all the markers independent
#' and imputation is conducted by SNPs.
#' \item Genotype likelihood (GL): when the GL info is detected
#' (GL column in FORMAT field of VCF files),
#' genotypes with higher likelihood will have higher probability during bootstrap
#' samples of trees in random Forests. Note that the use of genotype likelihoods
#' in the form of normalized, phred-scaled likelihoods (PL, e.g. from GATK)
#' are not recognized, yet... it's still under development.
#' \item Predictive mean matching: the rf option uses a fast k-nearest neighbor
#' (KNN) searching algorithms (see argument documentation and details below).
#' \item Optimized for speed: the package
#' \href{https://github.com/imbs-hl/ranger}{ranger}
#' (see Wright and Ziegler, 2016) provides a fast C++ version
#' of the original implementation of rf from Breiman (2001), and
#' imputations of genotypes are conducted in parallel across CPUs.
#' A progress bar is now available to see if you have time for a coffee break!
#' }
#'
#'
#' Before running this function to populate the original dataset with synthetic
#' data you should try \code{\link[stackr]{missing_visualization}}
#' to detect patterns of missingness.
#' Follow the \href{https://www.dropbox.com/s/4zf032g6yjatj0a/vignette_missing_data_analysis.nb.html?dl=0}{vignette}
#' for more info.


#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. To import, the function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{tidy_wide}}.
#' \emph{See details of this function for more info}.

#' @param imputation.method (character, optional) 
#' Methods available for map-independent imputations of missing genotype: 
#' (1) \code{imputation.method = "max"} to use the most frequent category for 
#' the imputations.
#' (2) \code{imputation.method = "rf"} using Random Forests algorithm.
#' \code{imputation.method = NULL} will return the original dataset, without
#' imputation.
#' Default: \code{imputation.method = "rf"}. 

#' @param imputations.group (character, optional) \code{"global"} or
#' \code{"populations"}.
#' Should the imputations be computed globally or by populations. Note that
#' imputing genotype globally can create huge bias for example by
#' introducing foreign genotypes in some populations.
#' Default = \code{"populations"}.

#' @param num.tree (integer, optional) The number of trees to grow with the 
#' Random Forests approach.
#' Default: \code{num.tree = 50}.

#' @param pred.mean.matching (integer, optional) Used in conjunction with 
#' random Forests. Number of candidate non-missing
#' value to sample from during the predictive mean matching step.
#' A fast k-nearest neighbor searching algorithms is used with this approach.
#' \code{pred.mean.matching = 3} will use 3 nighbors.
#' Default: \code{pred.mean.matching = 0}, avoids this step. 

#' @param random.seed (integer, optional) For reproducibility, set an integer
#' that will be used to initialize the random generator. With default,
#' a random number is generated.
#' Default: \code{random.seed = NULL}.

#' @param verbose (optional, logical) When \code{verbose = TRUE} 
#' the function is a little more chatty during execution.
#' Default: \code{verbose = TRUE}.

#' @param parallel.core (optional) The number of core used for parallel
#' execution of Random Forests imputations.
#' Markers are imputed in parallel, populations are processed sequentially.
#' Default: \code{parallel::detectCores() - 1}.

#' @param filename (optional) The file name for the \strong{imputed} 
#' tidy data frame written to the working directory.
#' With missing argument or default: \code{filename = NULL}, the imputed 
#' tidy data is in the global environment only 
#' (i.e. not written in the working directory).

#' @return The output in your global environment is the imputed tidy data frame.
#' If \code{filename} is provided, the imputed tidy data frame is also 
#' written to the working directory. The original data is returned for markers
#' with \emph{all} or \emph{no} NA.

#' @details
#' \strong{haplotype/SNP approach:}
#' 
#' The \strong{haplotype approach} is automatically used when markers meta-information
#' is detected (chromosome/CHROM, locus/ID and SNP/POS columns, usually from a VCF file).
#' Missing genotypes from SNPs on the same locus or same RADseq read is undertaken
#' simulteneously to account for the correlation of the linked SNPs.
#' Alternatively, a \strong{snp approach} is used, and the SNP are considered
#' independent. Imputations of genotypes is then conducted for each marker separately.
#' 
#' 
#' \strong{Predictive mean matching:}
#' 
#' Random Forests already behave like a nearest neighbor
#' classifier, with adaptive metric. Now we have the option to conduct
#' predictive mean matching on top of the prediction based missing value
#' imputation.PMM tries to raise the variance in the resulting conditional
#' distributions to a realistic level.
#' The closest k predicted values are identified by a fast
#' k-nearest neighbour approach using \code{\link[FNN]{knnx.index}} function
#' wrapped in the package \code{\link[missRanger]{pmm}}.
#' Returned value correspond to the mean value.

#' @note
#' 
#' \strong{Reference genome or linkage map available ?}
#' 
#' Numerous approaches are available and more appropriate, please search
#' the literature
#' (\href{https://online.papersapp.com/collections/05d6e65a-73c9-49e6-9c75-289a818f76f3/share}{references}).
#' 
#' 
#' \strong{Deprecated arguments:}
#' 
#' \itemize{
#' \item \code{impute} is no longer available.
#' Imputing using \code{impute = "allele"} option was wrong because it
#' was using F1 genotypes for imputations. Now imputation is only conducted at
#' the genotype level.
#' \item \code{iteration.rf} is no longer necessary and iteration is now set to a
#' maximum of 10 000 wich is more than enough. Most RADseq dataset with Random
#' Forests approach will reach consensus before 10 or 15 iterations.
#' \item \code{split.number} is automatically set.
#' }

#' @seealso
#' \href{https://github.com/mayer79/missRanger}{missRanger}
#' 
#' \href{https://github.com/imbs-hl/ranger}{ranger}
#' 
#' \href{https://github.com/stekhoven/missForest}{missForest}
#' 
#' \href{https://github.com/kogalur/randomForestSRC}{randomForestSRC}


#' @export
#' @rdname stackr_imputations_module
#' @importFrom parallel detectCores
#' @importFrom dplyr distinct group_by ungroup rename arrange tally filter filter_ select select_ one_of mutate mutate_all summarise left_join funs bind_rows
#' @importFrom tidyr gather unite
#' @importFrom purrr map flatten keep discard flatten_chr flatten_dbl flatten_lgl
#' @importFrom stringi stri_replace_na
#' @importFrom tibble has_name as_data_frame
#' @importFrom stats predict reformulate
#' @importFrom lazyeval interp
#' @importFrom ranger ranger
#' @importFrom missRanger pmm


#' @examples
#' \dontrun{
#' # The simplest way to run when you have a tidy dataset:
#' 
#' wolf.imputed <- stackr::stackr_imputations_module(data = "wolf.tidy.dataset.tsv")
#' 
#' # This will impute the missing genotypes by population using random Forests.
#' # The remaining arguments will be the defaults.
#' 
#' # When you start with a vcf file you can use magrittr %>% to `pipe` the
#' # result, below an example with more arguments offered by the functions:
#' 
#' wolf.imp <- stackr::tidy_genomic_data(
#'     data = "batch_1.vcf",
#'     strata = "strata.wolf.10pop.tsv",
#'     vcf.metadata = TRUE,
#'     whitelist.markers = "whitelist.loci.txt",
#'     verbose = TRUE) %>%
#' stackr::stackr_imputations_module(
#'     data = wolf.tidy, pred.mean.matching = 3, parallel.core = 32)
#' }

#' @references Wright, M. N. & Ziegler, A. (2016).
#' ranger: A Fast Implementation of Random Forests for High Dimensional Data
#' in C++ and R.
#' Journal of Statistical Software, in press. http://arxiv.org/abs/1508.04409.
#' @references Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

stackr_imputations_module <- function(
  data,
  imputation.method = "rf",
  imputations.group = "populations",
  num.tree = 50,
  pred.mean.matching = 0,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1,
  random.seed = NULL,
  filename = NULL
) {
  
  if (verbose) {
    cat("\n\n")
    cat("###############################################################################\n")
    cat("##################### stackr::stackr_imputations_module #######################\n")
    cat("###############################################################################\n")
  }
  # for timing
  timing <- proc.time()
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  
  # Set seed for sampling reproducibility
  if (is.null(random.seed)) {
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }
  
  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    input <- stackr::tidy_wide(data = data, import.metadata = TRUE)
  } else {
    input <- data
  }
  
  # output the proportion of missing genotypes BEFORE imputations
  na.before <- dplyr::summarise(.data = input, MISSING = round(length(GT[GT == "000000"])/length(GT), 3)) %>%
    purrr::flatten_dbl(.) %>% format(., scientific = FALSE)
  message("Proportion of missing genotypes before imputations: ", na.before)
  
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
  if (tibble::has_name(input, "GL")) {
    input <- dplyr::select(
      .data = input,
      dplyr::one_of("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS", "GT", "GL")
    )
    if (verbose) message("Genotype likelihood (GL) column detected")
  }
  
  # keep stratification
  strata.before <- dplyr::distinct(.data = input, INDIVIDUALS, POP_ID)
  
  # prepare data
  input <- dplyr::mutate(.data = input, GT = replace(GT, which(GT == "000000"), NA))
  
  # Imputation with Random Forests  ---------------------------------------------
  if (imputation.method == "rf") {
    if (verbose) message("Using Random Forests algorith")
    
    # SNP or haplotype approach ?
    if (tibble::has_name(input, "LOCUS") & tibble::has_name(input, "POS")) {
      input <- tidyr::unite(data = input, col = CHROM_LOCUS, CHROM, LOCUS) %>% 
        dplyr::select(-POS) # not necessary if we keep MARKERS
      if (verbose) message("Imputations computed by haplotype/read")
    } else {
      if (verbose) message("Imputations computed by markers (independent)")
    }
    
    # Random Forests by pop
    if (imputations.group == "populations") {
      if (verbose) message("Imputations computed by populations, take a break...")
      
      # Replace missing by NA and split data frame by population
      df.split.pop <- input %>% split(f = .$POP_ID)
      
      input.imp <- list() # to store results
      input.imp <- purrr::map(.x = df.split.pop,
                              .f = .stackr_imputer,
                              imputations.group = imputations.group,
                              num.tree = num.tree,
                              pred.mean.matching = pred.mean.matching,
                              random.seed = random.seed,
                              parallel.core = parallel.core) %>% 
        dplyr::bind_rows(.)
      df.split.pop <- NULL # remove unused objects
    } # End imputation RF populations
    
    # Random Forests global
    if (imputations.group == "global") { # Globally/overall
      if (verbose) message("Imputations computed globally, take a break...")
      input.imp <- list() # to store results
      input.imp <- .stackr_imputer(data = input,
                                   imputations.group = imputations.group,
                                   num.tree = num.tree,
                                   pred.mean.matching = pred.mean.matching,
                                   random.seed = random.seed,
                                   parallel.core = parallel.core)
    } # End imputation RF global
    
    # Reintroduce the stratification
    input.imp <- suppressWarnings(
      dplyr::left_join(strata.before,
                       dplyr::select(input.imp, MARKERS, INDIVIDUALS, GT)
                       , by = "INDIVIDUALS") %>% 
        dplyr::arrange(POP_ID, INDIVIDUALS) %>% 
        dplyr::ungroup(.)
    )
  } # End imputation RF
  
  # Imputations using the most common genotype ----------------------------------
  if (imputation.method == "max") { # End imputation max
    if (verbose) message("Using the most observed genotype per marker for imputations")
    if (imputations.group == "populations") {
      if (verbose) message("Imputations computed by populations")
      
      input.imp <- dplyr::select(input.imp, MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                      GT = replace(GT, which(GT == "NA"), NA)) %>%
        dplyr::ungroup(.)
      
    } # End imputation max populations 
    if (imputations.group == "global") {
      # Globally (not by pop_id)
      if (verbose) message("Imputations computed globally")
      
      input.imp <- dplyr::select(input.imp, MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(GT = stringi::stri_replace_na(str = GT, replacement = max(GT, na.rm = TRUE))) %>%
        dplyr::ungroup(.)
      
    } # End imputation max global
  } # End imputations max
  
  # prep results ---------------------------------------------------------------
  # Replace NA by 000000 in GT column
  input.imp$GT <- stringi::stri_replace_na(str = input.imp$GT, replacement = "000000")
  
  # Compute REF/ALT allele... might have change depending on prop of missing values
  if (ref.column) {
    if (verbose) message("Adjusting REF/ALT alleles to account for imputations...")
    input.imp2 <- stackr::change_alleles(data = input.imp)$input
  } # end computing REF/ALT
  
  # Integrate marker.meta columns
  if (!is.null(marker.meta)) {
    input.imp <- dplyr::left_join(input.imp, marker.meta, by = "MARKERS")
  }
  
  # Imputations reliability ------------------------------------------------------
  # Highlight potential problems with imputations
  
  
  # blacklist.markers.imputations <- dplyr::filter(.data = input, is.na(GT)) %>%
  #   dplyr::distinct(MARKERS)
  # prob.markers <- nrow(blacklist.markers.imputations)
  # if (prob.markers > 0) {
  #   n.markers <- dplyr::n_distinct(input$MARKERS)
  #   prop.prob.markers <- format(
  #     round(x = prob.markers/n.markers, digits = nchar(as.character(n.markers))),
  #     scientific = FALSE)
  #   prop.prob.markers
  # reintroduce the marker meta info (CHROM, LOCUS, POS)
  # if (!is.null(marker.meta)) {
  #   blacklist.markers.imputations <- dplyr::left_join(blacklist.markers.imputations, marker.meta, by = "MARKERS")
  # }
  
  # e.g. when only 2 genotypes present, and different, in the training set
  # blacklist problematic markers
  imputation.reliability <- FALSE
  if (imputation.reliability) {
    samples.per.pop <- dplyr::filter(.data = input, GT != "000000") %>%
      dplyr::distinct(POP_ID, INDIVIDUALS) %>% 
      dplyr::group_by(POP_ID) %>% 
      dplyr::tally(.)
    
    min.30 <- round(min(samples.per.pop$n) * 0.3, 0)
    
    info <- dplyr::filter(.data = input, GT != "000000") %>% 
      dplyr::distinct(MARKERS, POP_ID, GT) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n == 1) %>%
      dplyr::select(MARKERS, POP_ID) %>% 
      dplyr::ungroup(.) %>%
      dplyr::left_join(input, by = c("MARKERS", "POP_ID")) %>% 
      dplyr::filter(GT != "000000") %>% 
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::tally(.) %>% 
      dplyr::filter(n <= min.30) %>%
      dplyr::ungroup(.)
    
    if (nrow(info) > 0) {
      blacklist.markers.imputations <- dplyr::distinct(.data = info, MARKERS)
      if (!is.null(marker.meta)) {
        blacklist.markers.imputations <- dplyr::left_join(blacklist.markers.imputations, marker.meta, by = "MARKERS")
      }
      # readr::write_tsv(x = blacklist.markers.imputations, path = "blacklist.markers.imputations.tsv")
    }
  }#End imputation.reliability
  
  # Write to working directory
  if (!is.null(filename)) {
    if (verbose) message("Writing the imputed tidy data to the working directory: \n", filename)
    readr::write_tsv(x = input.imp, path = filename, col_names = TRUE)
  }
  
  # Missing after imputation:
  na.after <- dplyr::summarise(.data = input.imp, MISSING = round(length(GT[GT == "000000"])/length(GT), 3)) %>%
    purrr::flatten_dbl(.) %>% format(., scientific = FALSE)
  message("Proportion of missing genotypes after imputations: ", na.after)
  
  if (verbose) {
    # output the proportion of missing genotypes after imputations
    timing <- proc.time() - timing
    message("Computation time: ", round(timing[[3]]), " sec")
    cat("############################## completed ##############################\n")
  }
  return(input.imp)
} # End imputations

# Internal Function -----------------------------------------------------------

# stackr_imputer -----------------------------------------------------------
#' @title .stackr_imputer
#' @description imputations using Ranger package and predictive mean matching
#' @rdname .stackr_imputer
#' @export
#' @keywords internal

.stackr_imputer <- function(
  data,
  imputations.group,
  num.tree = 100,
  pred.mean.matching = 0,
  random.seed = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  # data <- df.split.pop[[7]]#test
  
  # message of progress for imputations by population
  if (imputations.group == "populations") {
    message("Imputations of pop: ", unique(data$POP_ID))
  }
  
  if (tibble::has_name(data, "CHROM_LOCUS")) {
    data <- data %>% split(f = .$CHROM_LOCUS)
  } else {
    data <- data %>% split(f = .$MARKERS)
  }
  
  data.imp <- list() # initiate list
  data.imp <- .stackr_parallel(
    X = data,
    FUN = .impute_genotypes,
    num.tree = num.tree,
    pred.mean.matching = pred.mean.matching,
    random.seed = random.seed,
    mc.cores = parallel.core
  ) %>% 
    dplyr::bind_rows(.)
  return(data.imp)
} #End .stackr_imputer

# impute_genotypes -------------------------------------------------------------
#' @title .impute_genotypes
#' @description imputations using Ranger package and predictive mean matching of missRanger
#' @rdname .impute_genotypes
#' @export
#' @keywords internal

.impute_genotypes <- function(
  data,
  num.tree = 100,
  pred.mean.matching = 0,
  random.seed = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  maxiter <-  10000
  # data <- data[["1_11204"]]#locus specific test
  # data <- data[[7]]# random test
  data <- tibble::as_data_frame(data)
  
  # GL info
  if (tibble::has_name(data, "GL")) {
    case.weights <- dplyr::filter(.data = data, !is.na(GL)) %>%
      dplyr::select(GL) %>%
      purrr::flatten_dbl(.)
    data <- dplyr::select(data, -GL)
  } else {
    case.weights <- NULL
  }
  
  # Dont' waist time imputing... some screening first
  all.missing <- nrow(dplyr::filter(.data = data, !is.na(GT)))
  
  if (all.missing == 0 || !anyNA(data$GT)) {
    # All NA or NO NA will return the original data
    imputed.dataset <- data
  } else {
    # to myself:
    # some test to find best method when > 1 SNP/haplotype
    # n.snp <- dplyr::n_distinct(data$MARKERS)
    # if (tibble::has_name(data, "CHROM_LOCUS")) {
    #   if (n.snp > 1) {
    #     data <- dplyr::group_by(.data = data, CHROM_LOCUS, POP_ID, INDIVIDUALS) %>% 
    #       tidyr::spread(data = ., key = MARKERS, value = GT)
    #   }
    # }
    
    # mutate columns to factors
    data <- dplyr::mutate_all(.tbl = data, .funs = factor, exclude = NA)
    
    if (tibble::has_name(data, "CHROM_LOCUS")) {
      data <- dplyr::select(.data = data, -CHROM_LOCUS)
      grouping.col <- "MARKERS"
      col.select <- c("MARKERS", "GT")
    } else {
      grouping.col <- "."
      col.select <- "GT"
    }
    complete <- dplyr::filter(data, !is.na(GT))
    # complete.id <- dplyr::select(.data = complete, -dplyr::one_of(col.select))
    complete.data <- dplyr::select(.data = complete, dplyr::one_of(col.select))
    missing <- dplyr::filter(data, is.na(GT))
    # missing.id <- dplyr::select(.data = missing, -dplyr::one_of(col.select))
    missing.data <- dplyr::select(.data = missing, dplyr::one_of(col.select))
    
    missing <- complete <- NULL
    # test
    # grouping.col <- purrr::keep(.x = colnames(data), .p = colnames(data) %in% c("CHROM_LOCUS", "POP_ID", "INDIVIDUALS", "MARKERS"))
    # grouping.col <- "MARKERS"
    # genotype.col <- purrr::discard(
    #   .x = colnames(data),
    #   .p = colnames(data) %in% c("MARKERS", "CHROM", "LOCUS", "POS", "CHROM_LOCUS", "INDIVIDUALS", "POP_ID", "GL", "PL"))
    # complete.data <- dplyr::filter_(
    #   .data = data,
    #   lazyeval::interp(~!is.na(genotype.col), genotype.col = as.name(genotype.col))
    # )
    # missing.data <- dplyr::filter_(
    #   .data = data,
    #   lazyeval::interp(~is.na(genotype.col), genotype.col = as.name(genotype.col))
    # )
    # missing.vector <- tibble::as_data_frame(is.na(data)) %>% 
    #   dplyr::select(-dplyr::one_of(grouping.col)) %>%
    #   purrr::flatten_lgl(.)
    # formula : stats::reformulate(grouping.col, response = genotype.col)
    
    i <- 1
    pred.error <- 1
    names(pred.error) <- "GT"
    oob.error <- TRUE # initiate
    while (oob.error && i <= maxiter) {
      data.last <- data
      pred.error.last <- pred.error
      genotype.col <- "GT"
      ranger.res <- ranger::ranger(
        formula = stats::reformulate(grouping.col, response = genotype.col),
        data = complete.data,
        num.trees = num.tree, 
        case.weights = case.weights, # used with GL...
        num.threads = parallel.core,
        seed = random.seed)
      
      predicted <- stats::predict(ranger.res$forest, missing.data)$predictions
      
      # predictive mean matching ---------------------------------------------
      if (pred.mean.matching > 0) {
        
        ytrain <- dplyr::select(complete.data, GT) %>% 
          dplyr::mutate(GT = as.character(GT)) %>% 
          purrr::flatten_chr(.)
        
        predicted <- missRanger::pmm(
          xtrain = ranger.res$predictions,
          xtest = predicted,
          ytrain = ytrain,
          k = pred.mean.matching)
      }
      
      data <- dplyr::select(.data = missing, -GT) %>%
        tibble::add_column(GT = predicted) %>% 
        dplyr::bind_rows(complete) %>%
        dplyr::arrange(MARKERS, INDIVIDUALS)
      
      pred.error[["GT"]] <- ranger.res$prediction.error
      
      if (is.nan(pred.error[[genotype.col]])) {
        pred.error[[genotype.col]] <- 0
      }
    }
    # update out-of-bag error
    oob.error <- mean(pred.error) < mean(pred.error.last)
    i <- i + 1 # update iteration
  }
  
  if (i == maxiter && oob.error || i == 2) {
    imputed.dataset <- data
  } else {
    imputed.dataset <- data.last
  }

imputed.dataset <- dplyr::mutate_all(.tbl = imputed.dataset,
                                     .funs = as.character, exclude = NA)
return(imputed.dataset)
}#End .impute_genotypes
