#' @name run_rxstacks
#' @title Run STACKS rxstacks module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/rxstacks.php}{rxstacks}
#' module inside R!


#' @param P (Character) Path to the directory containing STACKS files.
#' Default: \code{P = "06_ustacks_cstacks_sstacks"}.
#' Inside the folder \code{06_ustacks_cstacks_sstacks}, you should have:
#' \itemize{
#'   \item \strong{4 files for each samples:} The sample name is the prefix of
#'   the files ending with:
#' \code{.alleles.tsv.gz, .models.tsv.gz, .snps.tsv.gz, .tags.tsv.gz}.
#' Those files are created in the
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' module.
#' }

#' @param o (Character) Output path to write results.
#' Default: \code{o = "/07_rxstacks_cstacks_sstacks_populations"}.

#' @param strata (optional) Path to a strata file to get
#' the rxstacks results by population. With default, the results is presented
#' overall samples.
#' The strata file is a tab delimited file with 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column can be any hierarchical grouping.
#' If you have already run
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data,
#' the strata file is similar to a stacks \emph{population map file},
#' make sure you
#' have the required column names (\code{INDIVIDUALS} and \code{STRATA}).
#' Default: \code{strata = NULL}.



#' @param t (Integer) enable parallel execution with num_threads threads.
#' Default: \code{t = parallel::detectCores() - 1}

#' @param b (Integer) Database/batch ID of the input catalog to consider.
#' Advice: don't modify the default.
#' Default: \code{b = "guess"}.

#' @param lnl_filter (Logical) Filter catalog loci based on the mean log likelihood of
#' the catalog locus in the population.
#' Default: \code{lnl_filter = TRUE}.

#' @param lnl_lim (Double) Minimum log likelihood required to keep a catalog locus.
#' Default: \code{lnl_lim = NULL}.

#' @param lnl_dist (Logical) Print distribution of mean log likelihoods for
#' catalog loci AFTER the corrections. If you need the log likelihoods BEFORE the
#' correction use the output from \code{\link{run_sstacks}}, or if not generated,
#' use the function \code{\link{summary_catalog_log_lik}}.
#' Default: \code{lnl_dist = TRUE}.

#' @param conf_filter (Logical) Filter confounded loci.
#' Default: \code{conf_filter = TRUE}.

#' @param conf_lim (Double) Proportion of loci in population that must be
#' confounded relative to the catalog locus.
#' Default: \code{conf_lim = 0.75}.


#' @param prune_haplo (Logical) Prune out non-biological haplotypes unlikely to
#' occur in the population.
#' Default: \code{prune_haplo = TRUE}.

#' @param max_haplo (Integer) Only consider haplotypes for pruning if they occur
#' in fewer than \code{max_haplo_cnt samples}.
#' Default: \code{max_haplo = NULL}.

#' @param model_type Either 'snp' (default), 'bounded', or 'fixed'.
#' Default: \code{model_type = "snp"}.

#' @param alpha For the SNP or Bounded SNP model,
#' Chi square significance level required to call
#' a heterozygote or homozygote, either 0.1, 0.05.
#' Default: \code{alpha = 0.1}.
#' @param bound_low For the bounded SNP model, lower bound for epsilon,
#' the error rate, between 0 and 1.0.
#' Default: \code{bound_low = 0}.
#' @param bound_high For the bounded SNP model, upper bound for epsilon,
#' the error rate, between 0 and 1.0.
#' Default: \code{bound_high = 1}.

#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @param verbose Extended logging, including coordinates of all changed
#' nucleotides (forces single-threaded execution).
#' Advice: don't modify the default.
#' Default: \code{verbose = FALSE}



#' @rdname run_rxstacks



#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr group_by summarise bind_rows distinct arrange mutate ntile
#' @importFrom readr read_tsv
#' @importFrom stats median
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot geom_histogram labs facet_wrap theme element_text
#' @importFrom tibble add_column


#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/rxstacks.php}{rxstacks}
#' returns a \code{.matches.tsv.gz file for each sample}

#' @details \strong{Computer or server problem during the rxstacks ? Power outage? No problem:}
#' Restart the function as it was, the function will pick up at the sample it
#' was correcting before the problem.

#' @examples
#' \dontrun{
#' # First: get the log-likelihood distribution of catalog loci from
#' # \code{\link{run_sstacks}} output.
#' # If not run during \code{\link{run_sstacks}},
#' use the function \code{\link{summary_catalog_log.lik}}.
#' # Run rxstacks with appropriate filter thresholds:
#' rxstacks.tiger <- stackr::run_rxstacks(
#' lnl_filter = TRUE,
#' lnl_lim = -10,
#' lnl_dist = FALSE,
#' model_type = "bounded",
#' bound_low = 0, bound_high = 0.15,
#' strata = "population_map_tiger.txt")
#' # get all sorts of tables and figures:
#' rxstacks.tiger$rxstacks.correction.ind
#' rxstacks.tiger$rxstacks.correction.overall
#' rxstacks.tiger$rxstacks.correction.overall.proportion
#' rxstacks.tiger$rxstacks.correction.filled.bar.plot
#' rxstacks.tiger$rxstacks.correction.bar.plot
#' # Get the distribution of log-likelihood of catalog loci after correction:
#' rxstacks.tiger$log.likelihood.fig
#' }

#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/rxstacks.php}{rxstacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.


run_rxstacks <- function(
  P = "06_ustacks_cstacks_sstacks",
  o = "07_rxstacks_cstacks_sstacks_populations",
  strata = NULL,
  t = parallel::detectCores() - 1,
  b = "guess",
  lnl_filter = TRUE,
  lnl_lim = NULL,
  lnl_dist = TRUE,
  conf_filter = TRUE,
  conf_lim = 0.75,
  prune_haplo = TRUE,
  max_haplo = NULL,
  model_type = "snp",
  alpha = 0.1,
  bound_low = 0,
  bound_high = 1,
  h = FALSE,
  verbose = FALSE
) {

  cat("#######################################################################\n")
  cat("######################## stackr::run_rxstacks #########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  res <- list() # return results in this list
  parallel.core <- t

  # Check directory ------------------------------------------------------------
  if (!dir.exists("07_rxstacks_cstacks_sstacks_populations")) {
    dir.create("07_rxstacks_cstacks_sstacks_populations")
  }
  if (!dir.exists("09_log_files")) dir.create("09_log_files")


  # rxstacks options ------------------------------------------------------------
  if (b == "guess") {
    b <- ""
  } else {
    b <- stringi::stri_join("-b ", b)
  }
  ustacks.folder <- P
  rxstacks.folder <- o
  P <- stringi::stri_join("-P ", shQuote(P))
  o <- stringi::stri_join("-o ", shQuote(o))
  t <- stringi::stri_join("-t ", t)


  # Filtering options ----------------------------------------------------------

  if (lnl_filter) {
    lnl_filter <- stringi::stri_join("--lnl_filter ")
  } else {
    lnl_filter <- ""
  }

  if (is.null(lnl_lim)) {
    lnl_lim <- ""
  } else {
    lnl_lim <- stringi::stri_join("--lnl_lim ", lnl_lim)
  }

  if (lnl_dist) {
    lnl_dist <- stringi::stri_join("--lnl_dist ")
    print.log.likelihood <- TRUE
  } else {
    lnl_dist <- ""
    print.log.likelihood <- FALSE
  }

  if (conf_filter) {
    conf_filter <- stringi::stri_join("--conf_filter ")
  } else {
    conf_filter <- ""
  }

  if (is.null(conf_lim)) {
    conf_lim <- ""
  } else {
    conf_lim <- stringi::stri_join("--conf_lim ", conf_lim)
  }

  if (prune_haplo) {
    prune_haplo <- stringi::stri_join("--prune_haplo ")
  } else {
    prune_haplo <- ""
  }

  if (is.null(max_haplo)) {
    max_haplo <- ""
  } else {
    max_haplo <- stringi::stri_join("--max_haplo ", max_haplo)
  }

  # Model options --------------------------------------------------------------
  alpha <- stringi::stri_join("--alpha ", alpha)

  if (model_type == "bounded") {
    bound_low <- stringi::stri_join("--bound_low ", bound_low)
    bound_high <- stringi::stri_join("--bound_high ", bound_high)
  } else {
    bound_low <- ""
    bound_high <- ""
  }

  if (model_type == "fixed") {
    bc_err_freq <- stringi::stri_join("--bc_err_freq ", bc_err_freq)
  } else {
    bc_err_freq <- ""
  }

  model_type <- stringi::stri_join("--model_type ", model_type)



  # Logging Options: ------------------------------------------------------------
  # --verbose: extended logging, including coordinates of all changed nucleotides
  # (forces single-threaded execution).
  if (verbose) {
    verbose <- stringi::stri_join("--verbose ")
  } else {
    verbose <- ""
  }


  if (h) {
    h <- stringi::stri_join("-h ")
  } else {
    h <- ""
  }

  # logs file ------------------------------------------------------------------
  file.date.time <- stringi::stri_replace_all_fixed(
    str = Sys.time(),
    pattern = " EDT", replacement = "") %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("-", " ", ":"),
      replacement = c("", "@", ""),
      vectorize_all = FALSE
    ) %>%
    stringi::stri_sub(str = ., from = 1, to = 13)
  # command --------------------------------------------------------------------
  command.arguments <- paste(
    P, o, t, b, lnl_filter, lnl_lim, lnl_dist, conf_filter, conf_lim,
    prune_haplo, max_haplo, model_type, alpha, bound_low, bound_high, verbose, h
  )
  rxstacks.log.file <- stringi::stri_join("09_log_files/rxstacks_", file.date.time,".log")
  message(stringi::stri_join("For progress, look in the log file:\n", rxstacks.log.file))
  system2(command = "rxstacks", args = command.arguments, stderr = rxstacks.log.file)
  res$output <- "look inside output folder"

  # summarize rxstacks output  -------------------------------------------------
  rx.sum <- stackr::summary_rxstacks(
    rxstacks.folder = rxstacks.folder,
    strata = strata,
    parallel.core = parallel.core,
    verbose = FALSE,
    ... = file.date.time)

  res$rxstacks.correction.ind <- rx.sum$rxstacks.correction.ind
  res$rxstacks.correction.overall <- rx.sum$rxstacks.correction.overall
  res$rxstacks.correction.overall.proportion <- rx.sum$rxstacks.correction.overall.proportion
  res$rxstacks.correction.filled.bar.plot <- rx.sum$rxstacks.correction.filled.bar.plot
  res$rxstacks.correction.bar.plot <- rx.sum$rxstacks.correction.bar.plot
  res$rxstacks.haplo <- rx.sum$rxstacks.haplo
  res$rxstacks.snps <- rx.sum$rxstacks.snps


  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}# end run_rxstacks
