#' @name extract_catalog_sql_ids
#' @title Extract sample SQL IDs from STACKS catalog file
#' @description This function reads the output of
#' \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks} tags
#' file to generate a tibble of sample SQL ids that were used to generate the catalog.
#'
#' @param x (character, path). The path to the \code{catalog.tags.tsv.gz} file.
#' Default: \code{x = "catalog.tags.tsv.gz"}.

#' @rdname extract_catalog_sql_ids
#' @export

#' @return The function returns a tibble with the SQL IDs.

#' @examples
#' \dontrun{
#' ids <- stackr::extract_catalog_sql_id(x = "catalog.tags.tsv.gz")
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

extract_catalog_sql_ids <- function(x = "catalog.tags.tsv.gz") {

  # Test
  # x = "catalog.tags.tsv.gz"

  cat("#######################################################################\n")
  cat("################# stackr::extract_catalog_sql_ids #####################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  opt.change <- getOption("width")
  options(width = 70)

  if (missing(x)) stop("catalog tags file is required")

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  filename <- stringi::stri_join("catalog_sample_id_", file.date, ".tsv")
  tags.id <- readr::read_tsv(file = x, col_types = "____c____", skip = 1L, col_names = "ID")

  clean_tags <- function(x) {
    clean.x <- stringi::stri_extract_all_regex(str = x, pattern = "[0-9]+_", simplify = TRUE) %>%
      stringi::stri_replace_all_fixed(str = ., pattern = "_", replacement = "", vectorize_all = FALSE) %>%
      unique %>%
      as.integer
    return(clean.x)
  }

  tags.id <- tibble::tibble(
    SAMPLE_CATALOG =
      tags.id %>%
      dplyr::mutate(
        ID = purrr::map(.x = ID, .f = clean_tags)
      ) %>%
      # dplyr::select(SAMPLE_ID) %>%
      unlist(x = .) %>%
      unique %>%
      sort
  )
  message("Number of sample used to build the catalog: ", nrow(tags.id))

  if (file.exists("08_stacks_results")) {
    filename <- file.path("08_stacks_results", filename)
    readr::write_tsv(x = tags.id, path = filename)
    message("File written: ", filename)
  }

  timing <- proc.time() - timing
  options(width = opt.change)
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(tags.id)
} #End extract_catalog_sql_ids


