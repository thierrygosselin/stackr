#' @name summary_ustacks
#' @title Summarize STACKS ustacks files
#' @description This function reads inside the output of STACKS ustasks folder
#' and summarise the information.

#' @param ustacks.folder (logical). Default: \code{mismatch.testing = FALSE}.

#' @param parallel.core (optional) The number of core used for parallel
#' execution during vcf import.
#' Default: \code{parallel::detectCores() - 1}.

#' @param verbose (optional) Make the function a little more chatty during
#' execution.
#' Default: \code{verbose = FALSE}.

#' @rdname summary_ustacks
#' @export
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
#' @importFrom dplyr group_by tally ungroup summarise summarise_if mutate filter distinct select bind_rows n_distinct arrange
#' @importFrom readr read_tsv
#' @importFrom stats cor

#' @return The function returns a summary (data frame) containing:
#' \enumerate{
#' \item INDIVIDUALS: the sample id
#' \item LOCUS_NUMBER: the number of locus
#' \item BLACKLIST_USTACKS: the number of locus blacklisted by ustacks
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
#' ustacks.summary <- stackr::summary_ustacks(
#' ustacks.folder = "06_ustacks_cstacks_sstacks")
#'
#' # To get the number of snp/locus for a specific individual:
#' snp.info <- ustacks.summary %>%
#' dplyr::filter(INDIVIDUALS == "ID2") %>%
#' dplyr::select(INDIVIDUALS, SNP_LOCUS) %>%
#' tidyr::unnest(.)
#'
#' #Similarly, for the blacklisted locus (artifactual):
#' blacklist <- ustacks.summary %>%
#' dplyr::filter(INDIVIDUALS == "ID2") %>%
#' dplyr::select(INDIVIDUALS, BLACKLIST_ARTIFACT) %>%
#' tidyr::unnest(.)
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

summary_ustacks <- function(
  ustacks.folder,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE) {
  if (verbose) cat("#######################################################################\n")
  if (verbose) cat("##################### stackr::summary_ustacks #########################\n")
  if (verbose) cat("#######################################################################\n")
  timing <- proc.time()

  if (missing(ustacks.folder)) stop("ustacks.folder argument required")

  # alleles
  alleles.files <- list.files(
    path = ustacks.folder, pattern = "alleles", full.names = FALSE)

  # models
  models.files <- list.files(
    path = ustacks.folder, pattern = "models", full.names = FALSE)

  n.models <- length(models.files)
  # switch to tags if no models files (same thing, only more info to parse)

  tags.files <- list.files(
    path = ustacks.folder, pattern = "tags", full.names = FALSE)

  if (n.models == 0) {
    models.files <- tags.files
  }

  # snps
  snps.files <- list.files(
    path = ustacks.folder, pattern = "snps", full.names = FALSE)

  # remove catalog alleles
  catalog.file <- list.files(
    path = ustacks.folder, pattern = "catalog", full.names = FALSE)

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
    message("Summarizing ", length(alleles.files), " ustacks (models, snps, tags, alleles) files...")
  } else {
    message("Unequal numbers of ustacks files: models, snps, tags, alleles")
    message("  alleles: ", n.alleles)
    message("  snps: ", n.snps)
    message("  tags: ", n.tags)
    message("  models: ", n.models)
    message("  Some samples will be removed...selection based on .alleles files")
  }
  snps.files <- tags.files <- NULL
  opt.change <- getOption("width")
  options(width = 70)

  summarise_ustacks <- function(alleles.files, ustacks.folder) {
    sample.name <- stringi::stri_replace_all_fixed(
      str = alleles.files,
      pattern = ".alleles.tsv.gz",
      replacement = "",
      vectorize_all = FALSE)

    # summary
    models.summary <- readr::read_tsv(
      file = stringi::stri_join(ustacks.folder, "/", sample.name, ".models.tsv.gz"),
      col_names = c("SQL_ID", "ID", "LOCUS", "CHROMOSOME", "BASEPAIR", "STRAND", "SEQ_TYPE", "STACK_COMPONENT", "SEQ_ID", "SEQUENCE", "DELEVERAGED_FLAG","BLACKLISTED_FLAG", "LUMBERJACKSTACK_FLAG", "LOG_LIKELIHOOD"),
      col_types = "iiiciccicciiid", na = "-",
      comment = "#")

    summarise.ustacks <- dplyr::filter(models.summary, SEQ_TYPE == "model") %>%
      dplyr::select(LOCUS, SEQUENCE) %>%
      dplyr::mutate(
        SNP_NUMBER = stringi::stri_count_fixed(str = SEQUENCE, pattern = "E"),
        POLYMORPHISM = dplyr::if_else(SNP_NUMBER > 0, "het", "hom"),
        INDIVIDUALS = rep(sample.name, n())) %>%
      dplyr::select(INDIVIDUALS, LOCUS, POLYMORPHISM, SNP_NUMBER) %>%
      dplyr::full_join(
        dplyr::filter(models.summary,SEQ_TYPE == "consensus") %>%
          dplyr::mutate(
            BLACKLIST_USTACKS = dplyr::if_else(BLACKLISTED_FLAG == 1 | LUMBERJACKSTACK_FLAG == 1, "blacklist", "whitelist"),
            INDIVIDUALS = rep(sample.name, n())) %>%
          dplyr::select(INDIVIDUALS, LOCUS, BLACKLIST_USTACKS)
        , by = c("INDIVIDUALS", "LOCUS"))

    models.summary <- NULL


    summarise.ustacks <- dplyr::full_join(
      summarise.ustacks,
        # ALLELES
        readr::read_tsv(
          file = stringi::stri_join(ustacks.folder, "/", sample.name, ".alleles.tsv.gz"),
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
      summarise.ustacks %>%
        dplyr::group_by(INDIVIDUALS) %>%
        dplyr::summarise(
          LOCUS_TOTAL = n(),
          BLACKLIST_ARTIFACT = length(LOCUS[BLACKLIST_USTACKS == "whitelist" & BLACKLIST_ARTIFACT == "artifact"]),
          BLACKLIST_USTACKS = length(LOCUS[BLACKLIST_USTACKS == "blacklist"]),
          FOR_CATALOG = LOCUS_TOTAL - BLACKLIST_USTACKS,
          FILTERED = FOR_CATALOG - BLACKLIST_ARTIFACT) %>%
        dplyr::select(INDIVIDUALS,LOCUS_TOTAL, BLACKLIST_USTACKS, FOR_CATALOG, BLACKLIST_ARTIFACT, FILTERED),
      summarise.ustacks %>%
        dplyr::filter(BLACKLIST_USTACKS == "whitelist", BLACKLIST_ARTIFACT == "not_artifact") %>%
        dplyr::summarise(
          HOMOZYGOSITY = length(POLYMORPHISM[POLYMORPHISM == "hom"]),
          HETEROZYGOSITY = length(POLYMORPHISM[POLYMORPHISM == "het"]),
          MEAN_NUMBER_SNP_LOCUS = round(mean(SNP_NUMBER[SNP_NUMBER > 0]), 2),
          MAX_NUMBER_SNP_LOCUS = max(SNP_NUMBER[SNP_NUMBER > 0]),
          NUMBER_LOCUS_4SNP = length(LOCUS[SNP_NUMBER > 3]))
    ) #%>%
      #dplyr::mutate(LOCUS_INFO = list(summarise.ustacks))
    summarise.ustacks <- sample.name <- NULL
    return(summary.ind)
  }#End allele_summary

  res <- list()
  res <- .stackr_parallel(
    X = alleles.files,
    FUN = summarise_ustacks,
    mc.cores = parallel.core,
    ustacks.folder = ustacks.folder
  ) %>%
    dplyr::bind_rows(.)
  options(width = opt.change)

  mean.summary <- dplyr::summarise_if(.tbl = res, .predicate = is.numeric, .funs = mean)
  n.ind <- dplyr::n_distinct(res$INDIVIDUALS)

  if (n.ind > 1) {
    message("Number of individuals: ", n.ind)
    message("Mean number of:")
    message("  locus +/- \u00B1 SD [range]: ", round(mean.summary$LOCUS_TOTAL), " +/- ", round(sd(res$LOCUS_TOTAL)), " [", stringi::stri_join(range(res$LOCUS_TOTAL), collapse = " - "), "]")
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
}#End summary_ustacks
