#' @name generate_missing

#' @title Generate missing data

#' @description
#' Generate missing genotypes following a compound Dirichlet-multinomial distribution.

#' @param data 6 options: vcf (to make vcf population ready, see details below),
#' plink, stacks haplotype file, genind, genlight, genepop, 
#' and a data frame in wide format. 
#' \emph{See details} in \code{\link[stackr]{tidy_genomic_data}}.

#' @param filename (optional) The filename prefix for the objet in the global environment
#' or the working directory. Default: \code{filename = NULL}. A default name will be used,
#' customized with the output file(s) selected.

#' @inheritParams genomic_converter
#' @inheritParams stackr_imputations_module

#' @param average.read.depth (integer) Desired average read depth at a marker in an individual.
#' Default: \code{average.read.depth = 10}.
#' @param min.reads (integer) Minimum number of simulated reads to call a SNP.
#' Default: \code{min.reads = 6}.
#' @param alpha.individuals (integer) Shape parameter for gamma distribution over individuals.
#' Default: \code{alpha.individuals = 5}.
#' @param alpha.marker (integer) Shape parameter for gamma distribution over loci.
#' Default: \code{alpha.marker = 5}.
#' @param random.seed (integer) Random seed number for reproducibility. 
#' With default: \code{random.seed = NULL} the number is generated automatically and randomly.

# @param filename (optional) The name of the tidy data frame written to the directory.
# Use the extension ".tsv" at the end. 
# Several info will be appended to the name of the file.

#' @param ... (upcomming) other parameters passed to the function \code{\link[stackr]{tidy_genomic_data}} 
#' for the input file and \code{\link[stackr]{genomic_converter}} for the output file parameter.

#' @return In the global environment, a list with the tidy data set, the random.seed and function.call.
#' In the working directory, the output file with format selected.

#' @examples
#' \dontrun{
#' The simplest form of the function:
#' 
#' datamissing <- generate_missing(
#' data = "plink.tped",
#' output = "genepop"
#' )
#' }


#' @export
#' @rdname generate_missing
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename if_else mutate full_join
#' @importFrom stats rgamma rmultinom
#' @importFrom tibble data_frame as_data_frame
#' @importFrom tidyr unnest
#' @keywords internal

#' @author Eric C. Anderson \email{eric.anderson@@noaa.gov}, Greg L. Owens \email{gregory.owens@@alumni.ubc.ca} and Thierry Gosselin \email{thierrygosselin@@icloud.com}

generate_missing <- function(
  data,
  output,
  filename = NULL,
  average.read.depth = 10,
  min.reads = 6,
  alpha.individuals = 5,
  alpha.marker = 5,
  random.seed = NULL,
  imputation.method = NULL,
  impute = "genotype",
  imputations.group = "populations", 
  num.tree = 100, 
  iteration.rf = 10, 
  split.number = 100, 
  verbose = FALSE, 
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  
  
  # Import data ----------------------------------------------------------------
  message("WARNING: This function is still under testing, use with caution and report bugs\n\n")
  
  # Checking for missing and/or default arguments
  if (missing(data)) stop("Input file is missing")
  
  # input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  input <- stackr::tidy_genomic_data(
    data = data, 
    ...
  )
  
  # store function call
  function.call <- match.call()
  
  # For long tidy format, switch LOCUS to MARKERS column name, if found MARKERS not found
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # Individuals
  number.individuals <- dplyr::n_distinct(input$INDIVIDUALS)
  individuals.list <- dplyr::distinct(input, INDIVIDUALS) %>% purrr::flatten_chr(.)
  message(paste0("Number of individuals used ", number.individuals))
  
  # Number of markers
  number.markers <- dplyr::n_distinct(input$MARKERS)
  markers.list <- dplyr::distinct(input, MARKERS) %>% purrr::flatten_chr(.)
  message(paste0("Number of markers used ", number.markers))
  
  # Number of populations
  # number.populations <- dplyr::n_distinct(input$POP_ID)
  # pop.list <- dplyr::distinct(input, POP_ID) %>% purrr::flatten_chr(.)
  # message(paste0("Number of populations used ", number.populations))
  
  
  # Set seed -------------------------------------------------------------------
  if (is.null(random.seed)) {
    message("Generating random seed number")
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }
  
  # Simulation -----------------------------------------------------------------
  # simulate the number of reads per individual. Note that this is 
  # one way to simulate a compound Dirichlet-multinomial.
  
  # Simulate the Dirichlet as a bunch of gammas scaled by their sum, and then
  # set that as the cell probs in a multinomial.
  
  # the scale is set here only so things are large enough that there
  # is not a likely chance of underflow
  message("Simulating a compound Dirichlet-multinomial for the number of reads per individuals")
  gamma.read <- stats::rgamma(
    n = number.individuals, 
    shape = alpha.individuals, 
    scale = number.individuals * number.markers * average.read.depth / alpha.individuals
  )
  gamma.read.dirichlet <- gamma.read / sum(gamma.read)
  
  total.reads.per.individuals <- stats::rmultinom(
    n = 1,
    size = number.individuals * number.markers * average.read.depth, 
    prob = gamma.read.dirichlet
  ) %>% 
    tibble::as_data_frame(.) %>% 
    dplyr::rename(TOTAL_READ = V1) %>% 
    dplyr::mutate(INDIVIDUALS = individuals.list)
  
  gamma.read.dirichlet <- gamma.read <- NULL # unused arguments
  
  # then simulate the number of reads per locus within each individual. 
  # For this, we can use a CDM, and bind it altogether into a tidy data frame.
  # We want each locus to have its own characteristic rate at which reads come off of it,
  # so we simulate those first, and then use those to apportion reads from each individual
  
  message("Simulating a compound Dirichlet-multinomial for the number of reads per locus/marker")
  
  loc.gammas <- stats::rgamma(
    n = number.markers, 
    shape = alpha.marker, 
    scale = 100 / alpha.marker
  )
  # a dirichlet r.v. is a vector of gammas with common scale (scaled to sum to one)
  loc.dirichlet <- loc.gammas / sum(loc.gammas)
  
  # Function to apply to each ind and markers
  reads_per_markers <- function(total.reads, markers.list) {
    if (total.reads > 0) {
      res <- stats::rmultinom(n = 1, size = total.reads, prob = loc.dirichlet) %>% 
        tibble::as_data_frame(.) %>% 
        dplyr::rename(READ_DEPTH = V1) %>% 
        dplyr::mutate(MARKERS = markers.list)
    } else {
      res <- tibble::data_frame(READ_DEPTH = rep(0, length(markers.list))) %>%
        dplyr::mutate(MARKERS = markers.list)
    }
    return(res)
  }
  
  data.missing <- total.reads.per.individuals %>%
    dplyr::group_by(INDIVIDUALS) %>% 
    dplyr::mutate(READ_PER_MARKERS = purrr::map(.x = TOTAL_READ, .f = reads_per_markers, markers.list)) %>% 
    tidyr::unnest(.)
  
  data.missing <- data.missing %>% 
    dplyr::full_join(input, by = c("MARKERS", "INDIVIDUALS")) %>% 
    dplyr::mutate(GT = dplyr::if_else(READ_DEPTH < min.reads, "000000", GT)) %>% 
    dplyr::ungroup(.)
  
  # output ---------------------------------------------------------------------
  output <- stackr::genomic_converter(
    data = dplyr::select(.data = data.missing, MARKERS, POP_ID, INDIVIDUALS, GT),
    output = output, 
    monomorphic.out = FALSE, 
    common.markers = FALSE, 
    filename = filename, 
    verbose = FALSE,
    imputation.method = imputation.method,
    impute = impute,
    imputations.group = imputations.group, 
    num.tree = num.tree, 
    iteration.rf = iteration.rf, 
    split.number = split.number, 
    parallel.core = parallel.core,
    ...
  )
  
  if (is.null(imputation.method)) {
    output$tidy.data.imp <- "Imputations not selected"
  }
  return(res = list(tidy.data = data.missing, random.seed = random.seed, function.call = function.call, output = output$tidy.data.imp))
} # End of generate_missing function
