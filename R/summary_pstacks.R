#' @name summary_pstacks
#' @title Summarize STACKS pstacks files
#' @description This function reads the output of
#' \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/pstacks.php}{pstacks} folder
#' and summarise the information.
#'
#' The information shown can help to choose individuals to include/exclude from
#' the catalog, the next step in
#' \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{https://github.com/thierrygosselin/stackr#stacks-modules-and-radseq-typical-workflow}{pipeline}
#' (see example).
#'

#' @param pstacks.folder (character). The path to the pstacks output files.

#' @param parallel.core (integer, optional) The number of core used for parallel
#' execution.
#' Default: \code{parallel::detectCores() - 1}.

#' @param verbose (logical, optional) Make the function a little more chatty during
#' execution.
#' Default: \code{verbose = FALSE}.

#' @rdname summary_pstacks
#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr group_by tally ungroup summarise summarise_if mutate filter distinct select bind_rows n_distinct arrange
#' @importFrom readr read_tsv
#' @importFrom stats cor sd

#' @return The function returns a summary (data frame) containing:
#' \enumerate{
#' \item INDIVIDUALS: the sample id
#' \item LOCUS_NUMBER: the number of locus
#' \item BLACKLIST_PSTACKS: the number of locus blacklisted by pstacks
#' \item FOR_CATALOG: the number of locus available to generate the catalog
#' \item BLACKLIST_ARTIFACT: the number of artifact genotypes (> 2 alleles, see details)
#' \item FILTERED: the number of locus after artifacts are removed
#' \item HOMOZYGOSITY: the number of homozygous genotypes
#' \item HETEROZYGOSITY: the number of heterozygous genotypes
#' \item MEAN_NUMBER_SNP_LOCUS: the mean number of SNP/locus (excluding the artifact locus)
#' \item MAX_NUMBER_SNP_LOCUS: the max number of SNP/locus observed for this individual (excluding the artifact locus)
#' \item NUMBER_LOCUS_4SNP: the number of locus with 4 or more SNP/locus (excluding the artifact locus)
#' }


#' @examples
#' \dontrun{
#' pstacks.summary <- stackr::summary_pstacks(
#' pstacks.folder = "output_pstacks", verbose = TRUE)
#'
#' # To get the number of snp/locus for a specific individual:
#' snp.info <- pstacks.summary %>%
#' dplyr::filter(INDIVIDUALS == "ID2") %>%
#' dplyr::select(INDIVIDUALS, SNP_LOCUS) %>%
#' tidyr::unnest(.)
#'
#' #Similarly, for the blacklisted locus (artifactual):
#' blacklist <- pstacks.summary %>%
#' dplyr::filter(INDIVIDUALS == "ID2") %>%
#' dplyr::select(INDIVIDUALS, BLACKLIST_ARTIFACT) %>%
#' tidyr::unnest(.)
#'
#' # Decide individuals to include in catalog:
#'
#' # Make the pstacks.summary dataframe population wise by including pop info.
#' # For this we need a strata file ... it's similar to a population map, but with headers.
#' # It's a 2 columns files with INDIVIDUALS and STRATA as column header.
#' strata <- readr::read_tsv("strata.octopus.tsv", col_types = "cc") %>%
#' dplyr::rename(POP_ID = STRATA)
#'
#' # join the dataframes
#' individuals.catalog <- dplyr::left_join(pstacks.summary, strata, by = "INDIVIDUALS") %>%
#' dplyr::arrange(desc(FILTERED))
#'
#' # look at the FILTERED column
#' hist(individuals.catalog$FILTERED)
#' boxplot(individuals.catalog$FILTERED)
#' mean(individuals.catalog$FILTERED)
#'
#' # lets say we want to filter out individuals below 45000 markers
#' individuals.catalog.filtered <- individuals.catalog %>%
#' dplyr::filter(FILTERED > 45000)
#'
#' # number of ind per pop
#' dplyr::select(individuals.catalog.filtered, INDIVIDUALS, POP_ID) %>%
#' dplyr::group_by(POP_ID) %>% dplyr::tally(.)
#'
#' # choose 20% of individuals per pop to represent the catalog
#'
#' individuals.catalog.select <- individuals.catalog.filtered %>%
#' dplyr::group_by(POP_ID) %>%
#' dplyr::sample_frac(tbl = ., size = 0.2, replace = FALSE) %>%
#' dplyr::ungroup(.) %>%
#' dplyr::arrange(desc(FILTERED)) %>%
#' dplyr::distinct(INDIVIDUALS, POP_ID)
#'
#' # Using desc (from large to lower number of markers) ensure that
#' # the catalog is populated with samples with the highest number of loci first.
#' # This speed up the process in cstacks!
#'
#'
#' # Create a file with individuals and pop to use in
#' # cstacks or run_cstacks
#'
#' readr::write_tsv(
#' x = individuals.catalog.select,
#' path = "population.map.octopus.samples.catalog.tsv",
#' col_names = FALSE) # stacks doesn't want column header
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/pstacks.php}{pstacks}

#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks}


#' \code{\link{run_cstacks}}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

summary_pstacks <- function(
  pstacks.folder,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE) {
  if (verbose) cat("#######################################################################\n")
  if (verbose) cat("##################### stackr::summary_pstacks #########################\n")
  if (verbose) cat("#######################################################################\n")
  timing <- proc.time()

  if (missing(pstacks.folder)) stop("pstacks.folder argument required")

  # alleles
  alleles.files <- list.files(
    path = pstacks.folder, pattern = "alleles", full.names = FALSE)

  # models
  models.files <- list.files(
    path = pstacks.folder, pattern = "models", full.names = FALSE)

  n.models <- length(models.files)
  # switch to tags if no models files (same thing, only more info to parse)
  use.tags <- FALSE
  tags.files <- list.files(
    path = pstacks.folder, pattern = "tags", full.names = FALSE)

  if (n.models == 0) {
    models.files <- tags.files
    use.tags <- TRUE
  }

  # snps
  snps.files <- list.files(
    path = pstacks.folder, pattern = "snps", full.names = FALSE)

  # remove catalog alleles
  catalog.file <- list.files(
    path = pstacks.folder, pattern = "catalog", full.names = FALSE)

  if (length(catalog.file) >= 1) {
    alleles.files <- purrr::keep(.x = alleles.files,
                                 .p = !alleles.files %in% catalog.file)
    snps.files <- purrr::keep(.x = snps.files,
                              .p = !snps.files %in% catalog.file)

    models.files <- purrr::keep(.x = models.files,
                                .p = !models.files %in% catalog.file)

    tags.files <- purrr::keep(.x = tags.files,
                              .p = !tags.files %in% catalog.file)

    message("Removing these catalog files from the summary: \n  ",
            stringi::stri_join(catalog.file, collapse = "\n  "))
  }
  n.alleles <- length(alleles.files)
  n.snps <- length(snps.files)
  n.tags <- length(tags.files)
  tags.files <- NULL

  if (n.alleles == n.snps && n.alleles == n.models && n.alleles == n.tags) {
    message("Summarizing ", n.snps, " pstacks (models, snps, tags, alleles) files...")
    sample.name <- stringi::stri_replace_all_fixed(
      str = snps.files,
      pattern = ".snps.tsv.gz",
      replacement = "",
      vectorize_all = FALSE)
  } else {
    message("Unequal numbers of pstacks files: models, snps, tags, alleles")
    message("  alleles: ", n.alleles)
    message("  snps: ", n.snps)
    message("  tags: ", n.tags)
    message("  models: ", n.models)
    message("  These samples will be removed (based on .snps files): ")

    alleles.names <- stringi::stri_replace_all_fixed(
      str = alleles.files,
      pattern = ".alleles.tsv.gz",
      replacement = "",
      vectorize_all = FALSE)

    sample.name <- stringi::stri_replace_all_fixed(
      str = snps.files,
      pattern = ".snps.tsv.gz",
      replacement = "",
      vectorize_all = FALSE)

    missing.samples <- setdiff(alleles.names, sample.name)

    if (length(missing.samples) > 1) {
      missing.samples <- stringi::stri_join(missing.samples, collapse = ", ")
    }
    alleles.names <- NULL
    message("    ", missing.samples)
  }

  alleles.files <- tags.files <- NULL
  opt.change <- getOption("width")
  options(width = 70)

  summarise_pstacks <- function(sample.name, pstacks.folder, use.tags) {
    # sample.name <- stringi::stri_replace_all_fixed(
    #   str = snps.files,
    #   pattern = ".snps.tsv.gz",
    #   replacement = "",
    #   vectorize_all = FALSE)

    # sample.name <- "TOB_70"

    # summary
    if (use.tags) {
      models.summary <- readr::read_tsv(
        file = stringi::stri_join(pstacks.folder, "/", sample.name, ".tags.tsv.gz"),
        col_names = c("SQL_ID", "ID", "LOCUS", "CHROMOSOME", "BASEPAIR", "STRAND", "SEQ_TYPE", "STACK_COMPONENT", "SEQ_ID", "SEQUENCE", "DELEVERAGED_FLAG","BLACKLISTED_FLAG", "LUMBERJACKSTACK_FLAG", "LOG_LIKELIHOOD"),
        col_types = "iiiciccicciiid", na = "-",
        comment = "#")
    } else {
      models.summary <- readr::read_tsv(
        file = stringi::stri_join(pstacks.folder, "/", sample.name, ".models.tsv.gz"),
        col_names = c("SQL_ID", "ID", "LOCUS", "CHROMOSOME", "BASEPAIR", "STRAND", "SEQ_TYPE", "STACK_COMPONENT", "SEQ_ID", "SEQUENCE", "DELEVERAGED_FLAG","BLACKLISTED_FLAG", "LUMBERJACKSTACK_FLAG", "LOG_LIKELIHOOD"),
        col_types = "iiiciccicciiid", na = "-",
        comment = "#")
    }


    summarise.pstacks <- dplyr::filter(models.summary, SEQ_TYPE == "model") %>%
      dplyr::select(LOCUS, SEQUENCE) %>%
      dplyr::mutate(
        SNP_NUMBER = stringi::stri_count_fixed(str = SEQUENCE, pattern = "E"),
        POLYMORPHISM = dplyr::if_else(SNP_NUMBER > 0, "het", "hom"),
        INDIVIDUALS = rep(sample.name, n())) %>%
      dplyr::select(INDIVIDUALS, LOCUS, POLYMORPHISM, SNP_NUMBER) %>%
      dplyr::full_join(
        dplyr::filter(models.summary,SEQ_TYPE == "consensus") %>%
          dplyr::mutate(
            BLACKLIST_PSTACKS = dplyr::if_else(BLACKLISTED_FLAG == 1 | LUMBERJACKSTACK_FLAG == 1, "blacklist", "whitelist"),
            INDIVIDUALS = rep(sample.name, n())) %>%
          dplyr::select(INDIVIDUALS, LOCUS, BLACKLIST_PSTACKS)
        , by = c("INDIVIDUALS", "LOCUS"))

    models.summary <- NULL


    summarise.pstacks <- dplyr::full_join(
      summarise.pstacks,
        # ALLELES
        readr::read_tsv(
          file = stringi::stri_join(pstacks.folder, "/", sample.name, ".alleles.tsv.gz"),
          col_names = c("SQL_ID", "ID", "LOCUS", "HAPLOTYPE", "PERCENT", "COUNT"),
          col_types = "iiicdi", na = "-",
          comment = "#") %>%
          dplyr::group_by(LOCUS) %>%
          dplyr::tally(.) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate(INDIVIDUALS = rep(sample.name, n())) %>%
          dplyr::select(INDIVIDUALS, LOCUS, HAPLOTYPE_NUMBER = n)
      , by = c("INDIVIDUALS", "LOCUS")) %>%
      dplyr::mutate(
        HAPLOTYPE_NUMBER = replace(
          HAPLOTYPE_NUMBER, which(is.na(HAPLOTYPE_NUMBER)), 0),
        BLACKLIST_ARTIFACT = dplyr::if_else(HAPLOTYPE_NUMBER > 2,
                                            "artifact", "not_artifact")) %>%
      dplyr::ungroup(.)

    summary.ind <- dplyr::bind_cols(
      summarise.pstacks %>%
        dplyr::group_by(INDIVIDUALS) %>%
        dplyr::summarise(
          LOCUS_TOTAL = n(),
          BLACKLIST_ARTIFACT = length(LOCUS[BLACKLIST_PSTACKS == "whitelist" & BLACKLIST_ARTIFACT == "artifact"]),
          BLACKLIST_PSTACKS = length(LOCUS[BLACKLIST_PSTACKS == "blacklist"]),
          FOR_CATALOG = LOCUS_TOTAL - BLACKLIST_PSTACKS,
          FILTERED = FOR_CATALOG - BLACKLIST_ARTIFACT) %>%
        dplyr::select(INDIVIDUALS,LOCUS_TOTAL, BLACKLIST_PSTACKS, FOR_CATALOG, BLACKLIST_ARTIFACT, FILTERED),
      summarise.pstacks %>%
        dplyr::filter(BLACKLIST_PSTACKS == "whitelist", BLACKLIST_ARTIFACT == "not_artifact") %>%
        dplyr::summarise(
          HOMOZYGOSITY = length(POLYMORPHISM[POLYMORPHISM == "hom"]),
          HETEROZYGOSITY = length(POLYMORPHISM[POLYMORPHISM == "het"]),
          MEAN_NUMBER_SNP_LOCUS = round(mean(SNP_NUMBER[SNP_NUMBER > 0]), 2),
          MAX_NUMBER_SNP_LOCUS = max(SNP_NUMBER[SNP_NUMBER > 0]),
          NUMBER_LOCUS_4SNP = length(LOCUS[SNP_NUMBER > 3]))
    )
    summarise.pstacks <- sample.name <- NULL
    return(summary.ind)
  }#End summarise_pstacks

  if (length(sample.name) > 1) {
    res <- list()
    res <- .stackr_parallel(
      X = sample.name,
      FUN = summarise_pstacks,
      mc.cores = parallel.core,
      pstacks.folder = pstacks.folder,
      use.tags = use.tags
    ) %>%
      dplyr::bind_rows(.)
  } else {
    res <- summarise_pstacks(
      sample.name = sample.name,
      pstacks.folder = pstacks.folder,
      use.tags = use.tags
    )
  }

  options(width = opt.change)

  mean.summary <- dplyr::summarise_if(.tbl = res, .predicate = is.numeric, .funs = mean)
  n.ind <- dplyr::n_distinct(res$INDIVIDUALS)

  if (n.ind > 1) {
    message("Number of individuals: ", n.ind)
    message("Mean number of:")
    message("  locus +/- \u00B1 SD [range]: ", round(mean.summary$LOCUS_TOTAL), " +/- ", round(stats::sd(res$LOCUS_TOTAL)), " [", stringi::stri_join(range(res$LOCUS_TOTAL), collapse = " - "), "]")
    message("  snp/locus [max], excluding artifact locus: ", round(mean.summary$MEAN_NUMBER_SNP_LOCUS), " [", round(mean.summary$MAX_NUMBER_SNP_LOCUS), "]")
    message("  homozygous locus: ", round(mean.summary$HOMOZYGOSITY))
    message("  heterozygous locus: ", round(mean.summary$HETEROZYGOSITY))
    message("  artifactual locus: ", round(mean.summary$BLACKLIST_ARTIFACT))
    message("  locus with 4 or more SNPs (excluding artifactual locus): ", round(mean.summary$NUMBER_LOCUS_4SNP), "\n\n")

    message("Correlation between number of filtered locus and :")
    message("  heterozygous genotypes: ", round(stats::cor(res$FILTERED, res$HETEROZYGOSITY), 3))
    message("  homozygous genotypes: ", round(stats::cor(res$FILTERED, res$HOMOZYGOSITY), 3))
    message("  artifactual genotypes: ", round(stats::cor(res$FILTERED, res$BLACKLIST_ARTIFACT), 3))
  } else {
    message("Mean number of:")
    message("  locus [range]: ", round(mean.summary$LOCUS_TOTAL), " [", stringi::stri_join(range(res$LOCUS_TOTAL), collapse = " - "), "]")
    message("  snp/locus [max], excluding artifact locus: ", round(mean.summary$MEAN_NUMBER_SNP_LOCUS), " [", round(mean.summary$MAX_NUMBER_SNP_LOCUS), "]")
    message("  homozygous locus: ", round(mean.summary$HOMOZYGOSITY))
    message("  heterozygous locus: ", round(mean.summary$HETEROZYGOSITY))
    message("  artifactual locus: ", round(mean.summary$BLACKLIST_ARTIFACT))
    message("  locus with 4 or more SNPs (excluding artifactual locus): ", round(mean.summary$NUMBER_LOCUS_4SNP), "\n\n")
  }
  timing <- proc.time() - timing
  if (verbose) message("\nComputation time: ", round(timing[[3]]), " sec")
  if (verbose) cat("############################## completed ##############################\n")
  return(res)
}#End summary_pstacks
