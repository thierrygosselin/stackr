#' @name run_cstacks
#' @title Run STACKS cstacks module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks}
#' module inside R! The function runs a summary of the log file automatically
#' at the end (\code{\link{summary_cstacks}}). In the event of a power outage,
#' computer or cluster crash, just re-run the function. The function will start
#' over from the last catalog generated.

#' @param P path to the directory containing STACKS files.
#' Default: \code{P = "06_ustacks_2_gstacks"}.
#' Inside the folder \code{06_ustacks_2_gstacks}, you should have:
#' \itemize{
#'   \item \strong{4 files for each samples:} The sample name is the prefix of
#'   the files ending with:
#' \code{.alleles.tsv.gz, .models.tsv.gz, .snps.tsv.gz, .tags.tsv.gz}.
#' Those files are created in the
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' module.
#' }

#' @param o Output path to write catalog.
#' Default: \code{o = "06_ustacks_2_gstacks"}

#' @param M path to a population map file (Required when P is used).
#' Default: \code{M = "02_project_info/population.map.catalog.tsv"}.

#' @param n number of mismatches allowed between sample loci when build the catalog.
#' Default: \code{n = 1}

#' @param parallel.core Enable parallel execution with num_threads threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}

#' @param catalog.path This is for the "Catalog editing" part in cstacks where
#' you can provide the path to an existing catalog.
#' cstacks will add data to this existing catalog.
#' With default: \code{catalog.path = NULL} or with a supplied path, the function
#' The function scan automatically for the presence of a catalog inside the input folder.
#' If none is found, a new catalog is created.
#' If your catalog is not in the input folder, supply a path here.
#' e.g. \code{catalog.path = ~/catalog_folder}, the catalog files are inside the
#' P folder along the samples files and detected automatically.
#' If a catalog is detected in the input folder,
#' the samples in the \code{sample.list} argument
#' will be added in this catalog. The catalog is made of 3 files:
#' \code{catalog.alleles.tsv.gz, catalog.snps.tsv.gz, catalog.tags.tsv.gz}

#' @param max.gaps The number of gaps allowed between stacks before merging.
#' Default: \code{max.gaps = 2}

#' @param min.aln.len The minimum length of aligned sequence in a gapped
#' alignment.
#' Default: \code{min.aln.len = 0.8}

#' @param disable.gapped Disable gapped alignments between stacks.
#' Default: \code{disable.gapped = FALSE} (use gapped alignments).

#' @param k.len Specify k-mer size for matching between between catalog loci
#' (automatically calculated by default).
#' Advice: don't modify.
#' Default: \code{k.len = NULL}

#' @param report.mmatches Report query loci that match more than one catalog locus.
#' Advice: don't modify.
#' Default: \code{report.mmatches = FALSE}

#' @param split.catalog (integer) In how many samples you want to split the
#' catalog population map. This allows to have a backup catalog every
#' \code{split.catalog} samples. Their is obviously a trade-off between the
#' integer use here, the time to initialize an existing catalog and
#' re-starting from zero if everything crash.
#' Default: \code{split.catalog = 20}. Very useful on a personal computer or
#' university computer cluster....


#' @rdname run_cstacks
#' @export
#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}
#' returns a \code{.matches.tsv.gz file for each sample}

#' @details \strong{Computer or server problem during the cstacks ?} Look
#' in the log file to see which individuals remains to be included. Create a
#' new list of individuals to include and use the catalog.path argument to point
#' to the catalog created before the problem.

#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' run_cstacks()
#' # that's it ! Now if you have your own workflow folders, etc. See below.
#' Next example, let say you only want to include 10 individuals/pop and
#' include in the catalog samples with more than 2000000 reads. With the project
#' info file in the global environment:
#' library(tidyverse)
#' individuals.catalog <- project.info.file) %>%
#' filter(RETAINED > 2000000) %>%
#' group_by(POP_ID) %>%
#' sample_n(size = 10, replace = FALSE) %>%
#' ungroup %>%
#' arrange(desc(RETAINED)) %>%
#' distinct(INDIVIDUALS_REP, POP_ID)
#' # Write file to disk
#' readr::write_tsv(x = individuals.catalog,
#' file = "02_project_info/population.map.catalog.tsv")
#' # The next line will give you the list of individuals to include
#' individuals.catalog <- individuals.catalog$INDIVIDUALS_REP
#'
#' # To keep your info file updated with this information:
#' project.info.file <- project.info.file %>%
#' mutate(CATALOG = if_else(INDIVIDUALS_REP %in% individuals.catalog,
#' true = "catalog", false = "not_catalog")
#' )
#' write_tsv(project.info.file, "project.info.catalog.tsv")
#'
#' # Then run the command this way:
#' run_cstacks (
#' P = "06_ustacks_2_gstacks",
#' catalog.path = NULL,
#' n = 1,
#' parallel.core = 32,
#' h = FALSE,
#' max.gaps = 2, min.aln.len = 0.8,
#' k.len = NULL, report.mmatches = FALSE
#' )
#' }

#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

run_cstacks <- function(
    P = "06_ustacks_2_gstacks",
    o = "06_ustacks_2_gstacks",
    M = "02_project_info/population.map.catalog.tsv",
    catalog.path = NULL,
    n = 1,
    parallel.core = parallel::detectCores() - 1,
    max.gaps = 2, min.aln.len = 0.8, disable.gapped = FALSE,
    k.len = NULL, report.mmatches = FALSE,
    split.catalog = 20
) {

  # TEST
  # P = "06_ustacks_2_gstacks"
  # o = "06_ustacks_2_gstacks"
  # M = "02_project_info/population.map.catalog.tsv"
  # n = 3
  # parallel.core = parallel::detectCores() - 1
  # catalog.path = NULL
  # # catalog.path = "/Volumes/THIERRY_MAC/sturgeons_saskatchewan/sturgeon_saskatchewan/catalog_test/catalog_2"
  # max.gaps = 2
  # min.aln.len = 0.8
  # disable.gapped = FALSE
  # k.len = NULL
  # report.mmatches = FALSE
  # h = FALSE
  # split.catalog = 3


  cat("#######################################################################\n")
  cat("######################## stackr::run_cstacks ##########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # Check directory ------------------------------------------------------------
  if (!dir.exists(P)) stop("Missing P directory")
  if (!dir.exists("09_log_files")) dir.create("09_log_files")
  if (!dir.exists("08_stacks_results")) dir.create("08_stacks_results")

  # spliting catalog -----------------------------------------------------------
  if (split.catalog == 0L) split.catalog <- 1L
  pop.map <- readr::read_tsv(file = M, col_names = c("INDIVIDUALS", "STRATA"), col_types = "cc")
  if (nrow(pop.map) < split.catalog) {
    split.catalog <- nrow(pop.map)
  }

  # Crash session --------------------------------------------------------------
  # part added when the computer or current crash...
  if (is.null(catalog.path)) {
    catalog.path <- "06_ustacks_2_gstacks"
    catalog.output <- list.files(path = "08_stacks_results", pattern = "catalog_temp_", full.names = TRUE)
  } else {
    catalog.output <- 0L
  }
  # check the presence of split catalog

  if (length(catalog.output) > 0) {
    catstuff <- split_pop_map(
      pop.map = pop.map,
      split.catalog = split.catalog,
      catalog.path = catalog.path,
      write = FALSE
    )
    catstuff.bk <- catstuff
    catalog.output <- list.files(path = "08_stacks_results", pattern = "catalog_temp_", full.names = TRUE)
    cat.written <- length(purrr::map(.x = catalog.output, .f = list.files, pattern = "catalog.alleles.") %>% purrr::flatten_chr(.))

    if (cat.written < length(catstuff$list.catalog.numbers)) {
      message("Number of temporary catalog: ", cat.written)
      message("Number of catalog split required: ", length(catstuff$list.catalog.numbers))
      check.log <- sort(list.files(path = "09_log_files", pattern = "cstacks", full.names = FALSE))
      if (length(check.log) >= cat.written + 1) {
        n.log <- length(check.log)
        if (identical(check.log[(cat.written + 1)], check.log[(n.log)])) {
          message("Interrupted catalog construction, check log file: ", check.log[(cat.written + 1)])
        } else {
          message("Interrupted catalog construction, check log file: ", check.log[(n.log)])
        }
      }
      message("\nContinuing catalog construction no. ", cat.written + 1)

      catstuff$list.catalog.numbers <- catstuff$list.catalog.numbers[-(1:cat.written)]
      catstuff$catalog.output <- catstuff$catalog.output[-(1:cat.written)]
      catstuff$catalog.path <- catstuff$catalog.path[-(1:cat.written)]
      catstuff$catalog.pop.map <- catstuff$catalog.pop.map[-(1:cat.written)]
    }
  } else {
    catstuff <- split_pop_map(
      pop.map = pop.map,
      split.catalog = split.catalog,
      catalog.path = catalog.path,
      write = TRUE
    )
    catstuff.bk <- catstuff
    message("Number of catalog split: ", length(catstuff$list.catalog.numbers), "\n")
  }

  # run cstacks on split pop map -----------------------------------------------
  purrr::pwalk(
    .l = list(
      catalog.id = catstuff$list.catalog.numbers,
      o = catstuff$catalog.output,
      M = catstuff$catalog.pop.map,
      catalog.path = catstuff$catalog.path
    ),
    .f = run_split_cstacks,
    P = "06_ustacks_2_gstacks",
    n = n,
    parallel.core = parallel.core,
    max.gaps = max.gaps,
    min.aln.len = min.aln.len,
    disable.gapped = disable.gapped,
    k.len = k.len,
    report.mmatches = report.mmatches
  )

  # Copy the last catalog to the output folder ---------------------------------
  message("Moving last catalog generated to folder: ", P)
  file.copy(
    from = list.files(
      path = catstuff.bk$catalog.output[length(catstuff.bk$list.catalog.numbers)],
      full.names = TRUE
    ),
    to = P,
    overwrite = TRUE,
    recursive = TRUE,
    copy.date = TRUE
  )

  # message("\nCatalog files written in: ", P)

  timing <- proc.time() - timing
  message("\nOverall computation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
} # End run_cstacks

# split_pop_map ----------------------------------------------------------------
#' @title split_pop_map
#' @description Split the catalog map file
#' @rdname split_pop_map
#' @export
#' @keywords internal
split_pop_map <- function(pop.map, split.catalog, catalog.path, write = TRUE) {
  n.ind <- nrow(pop.map)

  pop.map.temp <- tibble::add_column(
    .data = pop.map,
    SPLIT_VEC = sort(rep.int(x = 1:ceiling(n.ind/split.catalog), times = split.catalog))[1:n.ind]
  ) %>%
    dplyr::group_by(SPLIT_VEC) %>%
    dplyr::group_split(.tbl = .)

  write_catalog_pop_map <- function(pop.map.temp, write) {
    c.id <- unique(pop.map.temp$SPLIT_VEC)
    c.path <- paste0("08_stacks_results/catalog_temp_", c.id)
    if (!dir.exists(c.path)) dir.create(c.path)
    c.filename <- file.path("02_project_info", paste0("pop_map_catalog_", c.id, ".tsv"))

    if (write) {
      readr::write_tsv(
        x = dplyr::select(pop.map.temp, -SPLIT_VEC),
        file = c.filename,
        col_names = FALSE
      )
    }
    return(res = list(catalog.path = c.path, catalog.pop.map = c.filename))
  } #End write_catalog_pop_map

  catalog.info <- purrr::map(.x = pop.map.temp, .f = write_catalog_pop_map, write = write)

  return(
    list(
      list.catalog.numbers = 1:length(pop.map.temp),
      catalog.output = purrr::map_chr(catalog.info, 1),
      catalog.path = c(catalog.path, purrr::map_chr(catalog.info, 1)[c(1:(length(pop.map.temp) - 1))]),
      catalog.pop.map = purrr::map_chr(catalog.info, 2)
    )
  )
}#End split_pop_map

# run_split_cstacks ----------------------------------------------------------------
#' @title run_split_cstacks
#' @description The function that runs cstacks
#' @rdname run_split_cstacks
#' @export
#' @keywords internal
run_split_cstacks <- function(
    catalog.id = NULL,
    o = "06_ustacks_2_gstacks",
    M = "02_project_info/population.map.catalog.tsv",
    catalog.path = NULL,
    P = "06_ustacks_2_gstacks",
    n = 1,
    parallel.core = parallel::detectCores() - 1,
    max.gaps = 2,
    min.aln.len = 0.8,
    disable.gapped = FALSE,
    k.len = NULL,
    report.mmatches = FALSE
) {
  timing.catalog <- proc.time()

  if (!is.null(catalog.id)) message("Generating catalog: ", catalog.id)

  # Existing catalog -----------------------------------------------------------
  if (is.null(catalog.path)) { # no catalog path, searching in the input path...
    old.catalog <- list.files(path = P, pattern = "catalog")
    # detect.rxstacks.lnls <- list.files(path = P, pattern = "rxstacks_lnls")
    # if (length(detect.rxstacks.lnls) == 1) {
    #   old.catalog <- purrr::discard(.x = old.catalog, .p = old.catalog %in% detect.rxstacks.lnls)
    # }
    if (length(old.catalog) > 0 & length(old.catalog) == 3) {
      message("Found a catalog in the input folder, using files: ")
      message(stringi::stri_join(old.catalog, "\n"))

      catalog.path <- stringi::stri_replace_all_fixed(
        str = old.catalog[1],
        pattern = ".catalog.alleles.tsv.gz",
        replacement = "",
        vectorize_all = FALSE
      )
      catalog.path <- stringi::stri_join(P, "/", catalog.path)
      catalog.path <- stringi::stri_join("--catalog ", shQuote(catalog.path))
    }
    if (length(old.catalog) > 0 & length(old.catalog) < 3) {
      rlang::abort("Incomplete catalog, 3 files are required, see argument documentation")
    }

    if (length(old.catalog) == 0) {
      message("Builing catalog for the first time")
      catalog.path <- ""
    }
  } else {
    old.catalog <- list.files(path = catalog.path, pattern = "catalog")

    # second time because after each pwalk, writting the catalog might take longer...
    if (length(old.catalog) == 0) {
      Sys.sleep(15)
      old.catalog <- list.files(path = catalog.path, pattern = "catalog")
    }


    # old.catalog <- purrr::keep(.x = old.catalog, .p = !old.catalog %in% "catalog.sql.ids.tsv")

    if (length(old.catalog) > 0) {
      if (length(old.catalog) < 3) {
        rlang::abort("Incomplete catalog, 3 files are required, see argument documentation")
      }
      message("Existing catalog: yes")
      message(stringi::stri_join(old.catalog, "\n"))


      # catalog.path <- file.path(
      #   catalog.path,
      #   stringi::stri_replace_all_fixed(
      #     str = old.catalog[1],
      #     pattern = ".alleles.tsv.gz",
      #     replacement = "",
      #     vectorize_all = FALSE
      #   ))
      catalog.path <- file.path(catalog.path, "catalog")
      catalog.path <- stringi::stri_join("--catalog ", shQuote(catalog.path))
    }

    if (length(old.catalog) == 0) {
      message("Builing catalog for the first time")
      catalog.path <- ""
    }
  }


  # cstacks options ------------------------------------------------------------
  # P <- stringi::stri_join("-P ", P)
  # M <- stringi::stri_join("-M ", M)

  # sample for catalog
  sc <- readr::read_tsv(file = M, col_names = c("INDIVIDUALS", "STRATA"), col_types = "cc") %>%
    dplyr::select(INDIVIDUALS) %>%
    purrr::flatten_chr(.)
  sc <- stringi::stri_join("-s ", file.path(P, sc), collapse = " ")

  n <- stringi::stri_join("-n ", n)
  p <- stringi::stri_join("-p ", parallel.core)
  o <- stringi::stri_join("-o ", o)

  # gapped assembly options ---------------------------------------------------
  max.gaps <- stringi::stri_join("--max-gaps ", max.gaps)
  min.aln.len <- stringi::stri_join("--min-aln-len ", min.aln.len)

  if (disable.gapped) {
    disable.gapped <- stringi::stri_join("--disable-gapped ")
  } else {
    disable.gapped <- ""
  }

  # Advanced options -----------------------------------------------------------
  if (is.null(k.len)) {
    k.len <- ""
  } else {
    k.len <- stringi::stri_join("--k-len ", k.len)
  }

  if (report.mmatches) {
    report.mmatches <- stringi::stri_join("--report-mmatches ")
  } else {
    report.mmatches <- ""
  }

  # logs files -----------------------------------------------------------------
  file.date.time <- format(Sys.time(), "%Y%m%d@%H%M")

  cstacks.log.file <- stringi::stri_join("09_log_files/cstacks_", file.date.time,".log")
  message(stringi::stri_join("For progress, look in the log file:\n", cstacks.log.file))


  # command args ---------------------------------------------------------------
  command.arguments <- paste(
    # P, M,
    sc,
    n, p, catalog.path, o,
    max.gaps, min.aln.len, disable.gapped,
    k.len, report.mmatches
  )

  # command
  system2(command = "cstacks", args = command.arguments, stderr = cstacks.log.file)


  # Summary cstacks ------------------------------------------------------------
  sum <- stackr::summary_cstacks(
    cstacks.log = cstacks.log.file,
    verbose = FALSE
  )
  sum <- NULL

  if (!is.null(catalog.id)) {
    timing.catalog <- proc.time() - timing.catalog
    message("\nComputation time to build catalog ", catalog.id, ": ", round(timing.catalog[[3]]), " sec\n")
  }
  timing.catalog <- timing <- NULL

}# end run_split_cstacks
