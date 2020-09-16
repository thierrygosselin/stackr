#' @name extract_catalog_sql_ids
#' @title Extract sample SQL IDs from STACKS catalog file
#' @description This function reads the output of
#' \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks} tags
#' file to generate a tibble of sample SQL ids that were used to generate the catalog.
#'
#' @param x (character, path). The path to the \code{catalog.tags.tsv.gz} file.
#' Default: \code{x = "catalog.tags.tsv.gz"}.

#' @param parallel.core (integer) Enable parallel execution with the number of threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @rdname extract_catalog_sql_ids
#' @export

#' @return The function returns a tibble with the SQL IDs.

#' @examples
#' \dontrun{
#' ids <- stackr::extract_catalog_sql_id(x = "06_ustacks_2_gstacks/catalog.tags.tsv.gz")
#' }


#' @seealso

#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks}

#' \code{\link{run_cstacks}}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

extract_catalog_sql_ids <- function(
  x = "06_ustacks_2_gstacks/catalog.tags.tsv.gz",
  parallel.core = parallel::detectCores() - 1
  ) {

  # Test
  # x = "06_ustacks_2_gstacks/catalog.tags.tsv.gz"
  # parallel.core=12
  cat("#######################################################################\n")
  cat("################# stackr::extract_catalog_sql_ids #####################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  opt.change <- getOption("width")
  options(width = 70)

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  filename <- stringi::stri_join("catalog_sqlids_", file.date, ".tsv")
  message("Reading catalog tags file...")
  tags.id <- readr::read_tsv(file = x, col_types = "____c____", skip = 1L, col_names = "ID")

  clean_tags <- function(x) {
    clean.x <- stringi::stri_extract_all_regex(str = x, pattern = "[0-9]+_", simplify = TRUE) %>%
      stringi::stri_replace_all_fixed(str = ., pattern = "_", replacement = "", vectorize_all = FALSE) %>%
      unique %>%
      as.integer
    return(clean.x)
  }


  n.tags <- nrow(tags.id)
  if (n.tags < 10) {
    tags.core <- n.tags
  } else {
    tags.core <- 10
  }
  cpu.rounds <- n.tags/tags.core/parallel.core
  tags.id <- tibble::add_column(
    .data = tags.id, SPLIT_VEC = as.integer(floor((parallel.core * cpu.rounds * (1:n.tags - 1) / n.tags) + 1))
  ) %>%
    dplyr::group_by(SPLIT_VEC) %>%
    dplyr::group_split(.tbl = ., .keep = FALSE)
  message("Extracting and summarizing SQL IDs information...")
  sql.ids <- .stackr_parallel(
      X = tags.id,
      FUN = clean_tags,
      mc.cores = parallel.core
      ) %>%
    unlist %>%
    unique %>%
    sort
  tags.id <- cpu.rounds <- tags.core <- n.tags <- NULL


  sql.ids <- tibble::tibble(CATALOG_SQL_IDS = sql.ids)
  message("Number of sample used to build the catalog: ", nrow(sql.ids))

  if (file.exists("08_stacks_results")) {
    filename <- file.path("08_stacks_results", filename)
    readr::write_tsv(x = sql.ids, path = filename)
    message("File written: ", filename)
  }

  timing <- proc.time() - timing
  options(width = opt.change)
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(sql.ids)
} #End extract_catalog_sql_ids


