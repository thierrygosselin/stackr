## Summary and tables

#' @title Haplotypes file summary
#' @description STACKS batch_x.haplotypes.tsv file summary.
#' The output of the function is a summary table for populations with:
#' \enumerate{
#' \item assembly artifacts and/or sequencing errors: loci with > 2 alleles.
#' Genotypes with more than 2 alleles are erased before estimating the subsequent statistics.
#' \item consensus reads (reads with no variation/polymorphism throughout the dataset).
#' \item monomorphic loci (at the population and overall level)
#' \item polymorphic loci (at the population and overall level)
#' \item haplotypes statistics for the observed and expected homozygosity and
#' heterozygosity (HOM_O, HOM_E, HET_O, HET_E).
#' \item Gene Diversity within populations (Hs) and averaged over all populations.
#' Not to confuse with Expected Heterozygosity(HET_E).
#' Here, Hs statistic includes a correction for sampling bias stemming from
#' sampling a limited number of individuals per population (Nei, 1987).
#' \item Wright’s inbreeding coefficient (Fis)
#' \item Nei's inbreeding coefficient (Gis) that include a correction for sampling bias.
#' \item \strong{IBDG}: a proxy measure of the realized proportion of the genome
#' that is identical by descent
#' \item \strong{FH measure}: an individual-base estimate that is
#' based on the excess in the observed number of homozygous
#' genotypes within an individual relative to the mean number of homozygous
#' genotypes expected under random mating (Keller et al., 2011; Kardos et al., 2015).
#' \item \strong{Pi}: the nucleotide diversity measured here consider the
#' consensus reads in the catalog (no variation between population sequences).
#' }

#' @param data The 'batch_x.haplotypes.tsv' created by STACKS.

#' @param strata A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' If you have already run
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data,
#' the strata file is similar to a stacks `population map file`, make sure you
#' have the required column names (\code{INDIVIDUALS} and \code{STRATA}).

#' @param read.length (number) The length in nucleotide of your reads
#' (e.g. \code{read.length = 100}).


#' @param whitelist.markers (optional) A whitelist of loci with a column name
#' 'LOCUS'. The whitelist is located in the global environment or in the
#' directory (e.g. "whitelist.txt").
#' If a whitelist is used, the files written to the directory will have
#' '.filtered' in the filename.
#' Default: \code{whitelist.markers = NULL}.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the global environment
#' or in the directory (with "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

#' @param keep.consensus (optional) Using whitelist of markers can automatically
#' remove consensus markers from the nucleotide diversity (Pi) calculations.
#' If you know what you are doing, play with this argument
#' to keep consensus markers.
#' Default: \code{keep.consensus = TRUE}.

#' @inheritParams tidy_genomic_data

#' @param parallel.core (optional) The number of core for parallel
#' programming during Pi calculations.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @importFrom stringdist stringdist
#' @importFrom utils combn count.fields
#' @importFrom stats lm na.omit
#' @importFrom stringi stri_replace_all_fixed stri_replace_na stri_join stri_count_fixed
#' @importFrom tibble as_data_frame data_frame add_column add_row
#' @importFrom dplyr select rename n_distinct distinct mutate summarise group_by ungroup arrange left_join full_join semi_join anti_join bind_rows bind_cols if_else
#' @importFrom readr write_tsv read_tsv
#' @importFrom tidyr separate gather
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid stat_smooth ggsave
#' @importFrom purrr flatten_chr map_df

#' @return The function returns a list with:
#' \enumerate{
#' \item $consensus.pop # if consensus reads are found
#' \item $consensus.loci # if consensus reads are found
#' \item $artifacts.pop # if artifacts are found
#' \item $artifacts.loci # if artifacts are found
#' \item $individual.summary: the individual's info for FH and Pi
#' \item $summary: the summary statistics per populations and averaged over all pops.
#' }
#' Also in the list 3 plots (also written in the folder):
#' \enumerate{
#' \item $scatter.plot.Pi.Fh.ind: showing Pi and FH per individuals
#' \item $scatter.plot.Pi.Fh.pop: showing Pi and FH per pop
#' \item $boxplot.pi: showing the boxplot of Pi per pop
#' \item $boxplot.fh: showing the boxplot of FH per pop
#' }
#' use $ to access each #' objects in the list.
#' The function potentially write 4 files in the working directory:
#' blacklist of unique loci with assembly artifacts and/or sequencing errors
#' (per individuals and globally), a blacklist of unique consensus loci
#' and a summary of the haplotype file by population.

#' @examples
#' \dontrun{
#' # The simplest way to run the function:
#' sum <- summary_haplotypes(
#' data = "batch_1.haplotypes.tsv",
#' strata = "strata_brook_charr.tsv",
#' read.length = 90)
#' }


#' @rdname summary_haplotypes
#' @export

#' @references Keller MC, Visscher PM, Goddard ME (2011)
#' Quantification of inbreeding due to distant ancestors and its detection
#'  using dense single nucleotide polymorphism data. Genetics, 189, 237–249.
#' @references Kardos M, Luikart G, Allendorf FW (2015)
#' Measuring individual inbreeding in the age of genomics: marker-based
#' measures are better than pedigrees. Heredity, 115, 63–72.
#' @references Nei M, Li WH (1979)
#' Mathematical model for studying genetic variation in terms of
#' restriction endonucleases.
#' Proceedings of the National Academy of Sciences of
#' the United States of America, 76, 5269–5273.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and
#' Anne-Laure Ferchaud \email{annelaureferchaud@@gmail.com}


summary_haplotypes <- function(
  data,
  strata,
  read.length,
  whitelist.markers = NULL,
  blacklist.id = NULL,
  keep.consensus = TRUE,
  pop.levels = NULL,
  pop.labels = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  cat("#######################################################################\n")
  cat("#################### stackr::summary_haplotypes #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  opt.change <- getOption("width")
  options(width = 70)
  res <- list() # to store results


  if (missing(data)) stop("data argument is missing")
  if (missing(strata)) stop("strata argument is required")
  if (missing(read.length)) stop("read.length argument is required")

  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }

  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")

  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  # Date and time --------------------------------------------------------------
  file.date <- stringi::stri_replace_all_fixed(
    Sys.time(),
    pattern = " EDT", replacement = "") %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("-", " ", ":"), replacement = c("", "@", ""),
      vectorize_all = FALSE) %>%
    stringi::stri_sub(str = ., from = 1, to = 13)

  folder.extension <- stringi::stri_join("summary_haplotypes_", file.date, sep = "")
  path.folder <- stringi::stri_join(getwd(),"/", folder.extension, sep = "")
  dir.create(file.path(path.folder))

  message(stringi::stri_join("Folder created: \n", folder.extension))
  file.date <- NULL #unused object


  # Import haplotype file ------------------------------------------------------
  message("Importing and tidying STACKS haplotype file: ", data)
  number.columns <- max(utils::count.fields(data, sep = "\t"))

  haplotype <- data.table::fread(
    input = data,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = list(character = 1:number.columns),
    verbose = FALSE,
    showProgress = TRUE,
    data.table = FALSE,
    na.strings = "-"
  ) %>%
    tibble::as_data_frame(.) %>%
    dplyr::select(-Cnt)

  if (tibble::has_name(haplotype, "# Catalog ID") || tibble::has_name(haplotype, "Catalog ID")) {
    colnames(haplotype) <- stringi::stri_replace_all_fixed(
      str = colnames(haplotype),
      pattern = c("# Catalog ID", "Catalog ID"), replacement = c("LOCUS", "LOCUS"), vectorize_all = FALSE
    )
  }

  if (tibble::has_name(haplotype, "Seg Dist")) {
    haplotype <- dplyr::select(.data = haplotype, -`Seg Dist`)
  }

  haplotype <- data.table::melt.data.table(
    data = data.table::as.data.table(haplotype),
    id.vars = "LOCUS",
    variable.name = "INDIVIDUALS",
    variable.factor = FALSE,
    value.name = "HAPLOTYPES"
  ) %>%
    tibble::as_data_frame(.)

  number.columns <- NULL

  haplotype$INDIVIDUALS = stringi::stri_replace_all_fixed(
    str = haplotype$INDIVIDUALS,
    pattern = c("_", ":"),
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )

  # Import whitelist of markers-------------------------------------------------
  if (is.null(whitelist.markers)) {
    haplotype <- dplyr::mutate(haplotype, WHITELIST = rep("whitelist", n()))
  } else {
    whitelist.markers <- suppressMessages(readr::read_tsv(whitelist.markers, col_names = TRUE))

    # Check that whitelist has "LOCUS" column
    if (!tibble::has_name(whitelist.markers, "LOCUS")) {
      stop("Whitelist requires a LOCUS column, check argument documentation...")
    }
    whitelist.markers <- dplyr::select(.data = whitelist.markers, LOCUS) %>%
      dplyr::mutate(LOCUS = as.character(LOCUS)) %>%
      dplyr::arrange(as.integer(LOCUS))

    haplotype <- haplotype %>%
      dplyr::mutate(
        WHITELIST = dplyr::if_else(LOCUS %in% whitelist.markers$LOCUS,
                                   "whitelist", "blacklist")
      )
    whitelist.markers <- "whitelist"
  }

  # population levels and strata------------------------------------------------
  message("Making the haplotype file population-wise...")
  if (is.vector(strata)) {
    # message("strata file: yes")
    number.columns.strata <- max(utils::count.fields(strata, sep = "\t"))
    col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
    strata.df <- suppressMessages(readr::read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>%
                                    dplyr::rename(POP_ID = STRATA))
  } else {
    # message("strata object: yes")
    colnames(strata) <- stringi::stri_replace_all_fixed(
      str = colnames(strata),
      pattern = "STRATA",
      replacement = "POP_ID",
      vectorize_all = FALSE
    )
    strata.df <- strata
  }

  # Remove potential whitespace in pop_id
  strata.df$POP_ID <- stringi::stri_replace_all_fixed(
    strata.df$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)

  # Import and filter blacklist id ---------------------------------------------
  if (!is.null(blacklist.id)) {
    blacklist.id <- suppressMessages(readr::read_tsv(blacklist.id, col_names = TRUE))
    message("Filtering with blacklist of individuals: ", nrow(blacklist.id), " individual(s) blacklisted")
    strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")

    blacklist.id$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = blacklist.id$INDIVIDUALS,
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )

    haplotype <- suppressWarnings(dplyr::anti_join(haplotype, blacklist.id, by = "INDIVIDUALS"))
    blacklist.id <- NULL
  }

  # Population levels and strata -----------------------------------------------
  strata.df$INDIVIDUALS = stringi::stri_replace_all_fixed(
    str = strata.df$INDIVIDUALS,
    pattern = c("_", ":"),
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )

  haplotype <- dplyr::left_join(x = haplotype, y = strata.df, by = "INDIVIDUALS")
  strata.df <- NULL

  # using pop.levels and pop.labels info if present
  haplotype <- change_pop_names(data = haplotype, pop.levels = pop.levels, pop.labels = pop.labels)

  pop <- unique(haplotype$POP_ID)

  # Locus with consensus alleles -----------------------------------------------
  message("Scanning for consensus markers...")
  consensus.pop <- haplotype %>%
    dplyr::filter(HAPLOTYPES == "consensus") %>%
    dplyr::distinct(LOCUS, POP_ID) %>%
    dplyr::mutate(CONSENSUS = rep("consensus", times = n())) %>%
    dplyr::arrange(as.numeric(LOCUS), POP_ID)

  if (length(consensus.pop$CONSENSUS) > 0) {
    # Create a list of consensus loci
    blacklist.loci.consensus <- consensus.pop %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(LOCUS) %>%
      dplyr::arrange(as.numeric(LOCUS))
    message("    Generated a file with ", nrow(blacklist.loci.consensus), " consensus loci")
    message("    Details in: blacklist.loci.consensus.txt")
    readr::write_tsv(
      x = blacklist.loci.consensus,
      path = stringi::stri_join(path.folder, "/blacklist.loci.consensus.txt"),
      col_names = TRUE
    )

    haplotype <- dplyr::mutate(
      haplotype,
      CONSENSUS = dplyr::if_else(
        LOCUS %in% blacklist.loci.consensus$LOCUS,
        "consensus", "not_consensus"))

    blacklist.loci.consensus.sum <- blacklist.loci.consensus %>%
      dplyr::ungroup(.) %>%
      dplyr::summarise(CONSENSUS = n()) %>%
      dplyr::select(CONSENSUS)

    consensus <- consensus.pop %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::summarise(CONSENSUS = dplyr::n_distinct(LOCUS))

    res$consensus.pop <- consensus.pop
    res$consensus.loci <- blacklist.loci.consensus
    consensus.pop <- blacklist.loci.consensus <- NULL
  } else {
    res$consensus.pop <- NULL
    res$consensus.loci <- NULL
    blacklist.loci.consensus.sum <- tibble::data_frame(CONSENSUS = as.integer(0))
    consensus <- tibble::data_frame(
      POP_ID = pop,
      CONSENSUS = as.integer(rep(0, length(pop)))
    )
    haplotype <- dplyr::mutate(haplotype, CONSENSUS = rep("not_consensus", n()))
  }

  # samples per pop
  sample.number <- dplyr::distinct(haplotype, INDIVIDUALS, POP_ID) %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::tally(.) %>%
    dplyr::rename(NUM = n) %>%
    tibble::add_row(.data = ., POP_ID = "OVERALL", NUM = sum(.$NUM))

  # Sequencing errors and assembly artifact ------------------------------------
  # Locus with > 2 alleles by individuals
  message("Scanning for assembly artifacts...")
  # haplotype <- haplotype %>%
  #   dplyr::mutate(
  #     ARTIFACTS = stringi::stri_count_fixed(HAPLOTYPES, "/"),
  #     ARTIFACTS = dplyr::if_else(ARTIFACTS > 1, "artifact", "not_artifact"))

  haplotype <- haplotype %>%
    dplyr::mutate(
      POLYMORPHISM = stringi::stri_count_fixed(HAPLOTYPES, "/"),
      POLYMORPHISM = dplyr::if_else(
        POLYMORPHISM > 1, "artifact",
        dplyr::if_else(
          POLYMORPHISM == 1, "het", "hom", missing = "missing")),
      POLYMORPHISM = dplyr::if_else(CONSENSUS == "consensus", "consensus", POLYMORPHISM)
    ) %>%
    dplyr::select(POP_ID, INDIVIDUALS, LOCUS, HAPLOTYPES, WHITELIST, POLYMORPHISM)

  # artifacts.ind <- dplyr::filter(
  #   haplotype,
  #   WHITELIST == "whitelist",
  #   CONSENSUS == "not_consensus")

  artifacts.ind <- dplyr::filter(
    haplotype,
    WHITELIST == "whitelist",
    POLYMORPHISM != "consensus")

  total.genotype.number.haplo <- length(stats::na.omit(artifacts.ind$HAPLOTYPES))

  # artifacts.ind <- artifacts.ind %>%
  #   dplyr::filter(ARTIFACTS == "artifact") %>%
  #   dplyr::arrange(LOCUS, POP_ID, INDIVIDUALS) %>%
  #   dplyr::select(LOCUS, POP_ID, INDIVIDUALS, HAPLOTYPES, ARTIFACTS)
  #
  artifacts.ind <- artifacts.ind %>%
    dplyr::filter(POLYMORPHISM == "artifact") %>%
    dplyr::arrange(LOCUS, POP_ID, INDIVIDUALS) %>%
    dplyr::select(LOCUS, POP_ID, INDIVIDUALS, HAPLOTYPES, ARTIFACTS = POLYMORPHISM)


  if (nrow(artifacts.ind) > 0) {
    erased.genotype.number <- length(artifacts.ind$HAPLOTYPES)
    percent.haplo <- paste(round(((erased.genotype.number/total.genotype.number.haplo)*100), 2), "%", sep = " ")

    # Write the list of locus, individuals with artifacts
    if (is.null(whitelist.markers)) {
      filename.artifacts.ind <- "blacklist.loci.artifacts.ind.txt"
    } else {
      filename.artifacts.ind <- "blacklist.loci.artifacts.ind.filtered.txt"
    }
    readr::write_tsv(
      x = artifacts.ind,
      path = stringi::stri_join(path.folder, "/", filename.artifacts.ind),
      col_names = TRUE)
    message("    Generated a file with ",
            nrow(artifacts.ind),
            " artifact loci by individuals")
    message("    Details in: ", filename.artifacts.ind)

    # res$artifacts.pop <- dplyr::ungroup(artifacts.ind) %>%
    #   dplyr::distinct(LOCUS, POP_ID, ARTIFACTS) %>%
    #   dplyr::arrange(as.numeric(LOCUS), POP_ID)

    res$artifacts.pop <- dplyr::ungroup(artifacts.ind) %>%
      dplyr::distinct(LOCUS, POP_ID, ARTIFACTS) %>%
      dplyr::arrange(as.numeric(LOCUS), POP_ID)

    artifacts <- res$artifacts.pop %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::tally(.) %>%
      dplyr::rename(ARTIFACTS = n)

    res$artifacts.loci <- dplyr::ungroup(res$artifacts.pop) %>%
      dplyr::distinct(LOCUS) %>%
      dplyr::arrange(as.numeric(LOCUS))

    blacklist.loci.artifacts.sum <- res$artifacts.loci %>%
      dplyr::ungroup(.) %>%
      dplyr::summarise(ARTIFACTS = n()) %>%
      dplyr::select(ARTIFACTS)

    # Write the unique list of paralogs blacklisted to a file
    if (is.null(whitelist.markers)) {
      filename.paralogs <- "blacklist.loci.artifacts.txt"
    } else {
      filename.paralogs <- "blacklist.loci.artifacts.filtered.txt"
    }
    readr::write_tsv(
      x = res$artifacts.loci,
      path = stringi::stri_join(path.folder, "/", filename.paralogs),
      col_names = TRUE
    )

    message("    Generated a file with ", nrow(res$artifacts.loci),
            " artifact loci")
    message("    Details in: ", filename.paralogs)

    message("    Out of a total of ", total.genotype.number.haplo, " genotypes")
    message("      ", percent.haplo, " (", erased.genotype.number, ")"," outlier genotypes will be erased")
    haplotype <- suppressWarnings(
      haplotype %>%
        dplyr::mutate(
          HAPLOTYPES = dplyr::if_else(POLYMORPHISM == "artifact", NA_character_, HAPLOTYPES)
        ))
  } else {
    res$artifacts.pop <- res$artifacts.loci <- NULL
    artifacts <- tibble::data_frame(POP_ID = pop, ARTIFACTS = as.integer(rep(0, length(pop))))
    blacklist.loci.artifacts.sum <- tibble::data_frame(ARTIFACTS = as.integer(0))
    # haplotype <- dplyr::mutate(haplotype, ARTIFACTS = rep("not_artifact", n()))
  }
  artifacts.ind <- NULL


  # bind info
  consensus.artifacts <- dplyr::full_join(artifacts, consensus, by = "POP_ID")
  artifacts <- consensus <- NULL

  # Summary dataframe by individual---------------------------------------------
  message("Genome-Wide Identity-By-Descent calculations (FH)...")

  # haplotype.filtered <- haplotype %>%
  #   dplyr::filter(
  #     WHITELIST == "whitelist",
  #     CONSENSUS == "not_consensus",
  #     ARTIFACTS == "not_artifact")

  haplotype.filtered <- haplotype %>%
    dplyr::filter(WHITELIST == "whitelist", POLYMORPHISM != "not_consensus")

  n.individuals <- dplyr::n_distinct(haplotype.filtered$INDIVIDUALS)
  n.markers <- dplyr::n_distinct(haplotype.filtered$LOCUS)

  # summary.ind <- dplyr::mutate(
  #   .data = haplotype.filtered,
  #   IND_LEVEL_POLYMORPHISM = dplyr::if_else(
  #     stringi::stri_count_fixed(HAPLOTYPES, "/") == 1,
  #     "het", "hom", missing = "missing"))

  # summary <- dplyr::mutate(
  #   .data = haplotype.filtered,
  #   IND_LEVEL_POLYMORPHISM = dplyr::if_else(
  #     stringi::stri_count_fixed(HAPLOTYPES, "/") == 1,
  #     "het", "hom", missing = "missing"))

  # Some could start averaging by individuals accross markers and then by populations...
  # GenoDive, Hierfstat and Nei computes by locus...

  if (n.markers < 5000) {
    summary.loc.pop <- haplotype.filtered %>%
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::summarise(
        HOM = length(POLYMORPHISM[POLYMORPHISM == "hom"]),
        HET = length(POLYMORPHISM[POLYMORPHISM == "het"]),
        N_GENOT = HOM + HET,
        HOM_O = HOM / N_GENOT,
        HET_O = HET / N_GENOT
      ) %>%
      dplyr::mutate(
        NN = 2 * N_GENOT,
        N_INV = N_GENOT / (N_GENOT - 1)
      ) %>%
      dplyr::ungroup(.)
  } else {
    loc_pop_stats <- function(x) {
      res <- dplyr::group_by(x, LOCUS, POP_ID) %>%
        dplyr::summarise(
          HOM = length(POLYMORPHISM[POLYMORPHISM == "hom"]),
          HET = length(POLYMORPHISM[POLYMORPHISM == "het"]),
          N_GENOT = HOM + HET,
          HOM_O = HOM / N_GENOT,
          HET_O = HET / N_GENOT
        ) %>%
        dplyr::mutate(
          NN = 2 * N_GENOT,
          N_INV = N_GENOT / (N_GENOT - 1)
        ) %>%
        dplyr::ungroup(.)
      return(res)
    }#End individual_stats

    split.vec <- dplyr::distinct(haplotype.filtered, LOCUS) %>%
      dplyr::mutate(SPLIT_VEC = as.integer(floor((parallel.core * 3 * (1:n.markers - 1) / n.markers) + 1)))

    summary.loc.pop <- haplotype.filtered %>%
      dplyr::left_join(split.vec, by = "LOCUS") %>%
      split(x = ., f = .$SPLIT_VEC) %>%
      .stackr_parallel(
        X = .,
        FUN = loc_pop_stats,
        mc.cores = parallel.core
      ) %>%
      dplyr::bind_rows(.)

    split.vec <- NULL
  }

  # HS Gene Diversity within populations ---------------------------------------
  # not the same as Expected Heterozygosity

  # haplotype.filtered.genotyped <- haplotype.filtered %>% dplyr::filter(!is.na(HAPLOTYPES))

  separate_haplo <- function(x) {
    res <- dplyr::select(x, LOCUS, POP_ID, INDIVIDUALS, HAPLOTYPES)  %>%
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"),
        sep = "/", extra = "drop", remove = TRUE
      ) %>%
      dplyr::mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
      dplyr::select(-INDIVIDUALS) %>%
      tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID)) %>%
      dplyr::select(-ALLELE_GROUP)
    return(res)
  }#End separate_haplo

  n.row <- nrow(haplotype.filtered)
  split.vec <- as.integer(floor((parallel.core * 20 * (1:n.row - 1) / n.row) + 1))
  n.row <- NULL

  # split.vec.locus <- dplyr::distinct(haplotype.filtered, LOCUS) %>%
  #   dplyr::mutate(
  #     SPLIT_VEC = factor(as.integer(floor((parallel.core * 20 * (1:n() - 1) / n()) + 1)))
  #   )

  haplotype.filtered.split <- split(x = haplotype.filtered, f = split.vec) %>%
    .stackr_parallel_mc(
      X = ., FUN = separate_haplo, mc.cores = parallel.core) %>%
    dplyr::bind_rows(.)

  hs <- haplotype.filtered.split %>%
    dplyr::group_by(LOCUS, POP_ID) %>%
    dplyr::count(ALLELES) %>%
    dplyr::summarise(HOM_E = sum((n / sum(n))^2)) %>%
    dplyr::ungroup(.) %>%
    dplyr::full_join(dplyr::select(summary.loc.pop, LOCUS, POP_ID, HET_O, NN, N_INV), by = c("LOCUS", "POP_ID")) %>%
    dplyr::group_by(LOCUS, POP_ID) %>%
    dplyr::mutate(HS = N_INV * (1 - HOM_E - HET_O / NN)) %>%
    dplyr::ungroup(.)

  summary.pop <- summary.loc.pop %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      HOM_O = mean(HOM_O, na.rm = TRUE),
      HET_O = mean(HET_O, na.rm = TRUE)
    ) %>%
    dplyr::full_join(
      dplyr::group_by(.data = hs, POP_ID) %>%
        dplyr::summarise(
          HOM_E = mean(HOM_E, na.rm = TRUE),
          HET_E = (1 - HOM_E),
          HS = mean(HS, na.rm = TRUE))
      , by = "POP_ID")

  summary.locus <- summary.loc.pop %>%
    dplyr::group_by(LOCUS) %>%
    dplyr::summarise(
      HOM_O = mean(HOM_O, na.rm = TRUE),
      HET_O = mean(HET_O, na.rm = TRUE)
    ) %>%
    dplyr::full_join(
      dplyr::group_by(.data = hs, LOCUS) %>%
        dplyr::summarise(
          HOM_E = mean(HOM_E, na.rm = TRUE),
          HET_E = (1 - HOM_E),
          HS = mean(HS, na.rm = TRUE))
      , by = "LOCUS")

  hs <- summary.loc.pop <- NULL

  summary.ind <- haplotype.filtered %>%
    dplyr::group_by(INDIVIDUALS, POP_ID) %>%
    dplyr::summarise(
      HOM = length(POLYMORPHISM[POLYMORPHISM == "hom"]),
      HET = length(POLYMORPHISM[POLYMORPHISM == "het"]),
      N_GENOT = HOM + HET,
      HOM_O = HOM / N_GENOT,
      HET_O = HET / N_GENOT
    ) %>%
    dplyr::ungroup(.)

  haplotype.filtered <- NULL

  # required functions
  # freq_hom <- function(x) {
  #   res <- dplyr::select(x, -SPLIT_VEC) %>%
  #     dplyr::group_by(LOCUS, POP_ID, ALLELES) %>%
  #     dplyr::summarise(
  #       FREQ_ALLELES = length(ALLELES)/mean(DIPLO),
  #       HOM_E = FREQ_ALLELES * FREQ_ALLELES
  #     ) %>%
  #     dplyr::select(-FREQ_ALLELES) %>%
  #     dplyr::ungroup(.)
  #   return(res)
  # }#End freq_hom

  # freq.hom.pop <- haplotype.filtered %>%
  #   dplyr::filter(!is.na(HAPLOTYPES))

  # diplo <- freq.hom.pop %>%
  #   dplyr::distinct(LOCUS, POP_ID, INDIVIDUALS) %>%
  #   dplyr::group_by(LOCUS, POP_ID) %>%
  #   dplyr::tally(.) %>%
  #   dplyr::ungroup(.) %>%
  #   dplyr::mutate(n = 2 * n) %>%
  #   dplyr::rename(DIPLO = n)

  # n.row <- nrow(haplotype.filtered.split)
  # split.vec <- as.integer(floor((parallel.core * 20 * (1:n.row - 1) / n.row) + 1))
  # n.row <- NULL

  # split.vec.locus <- dplyr::distinct(haplotype.filtered.split, LOCUS) %>%
  #   dplyr::mutate(
  #     SPLIT_VEC = factor(as.integer(floor((parallel.core * 20 * (1:n() - 1) / n()) + 1)))
  #   )

  # # freq.hom.pop <- split(x = freq.hom.pop, f = split.vec) %>%
  # #   .stackr_parallel_mc(
  # #     X = ., FUN = separate_haplo, mc.cores = parallel.core) %>%
  # #   dplyr::bind_rows(.) %>%
  # #   dplyr::left_join(diplo, by = c("LOCUS", "POP_ID")) %>%
  # #   dplyr::left_join(split.vec.locus, by = "LOCUS") %>%
  # #   split(x = ., f = .$SPLIT_VEC) %>%
  # #   .stackr_parallel_mc(
  # #     X = ., FUN = freq_hom, mc.cores = parallel.core) %>%
  # #   dplyr::bind_rows(.) %>%
  # #   dplyr::group_by(LOCUS, POP_ID) %>%
  # #   dplyr::summarise(HOM_E = sum(HOM_E)) %>%
  # #   dplyr::group_by(POP_ID) %>%
  # #   dplyr::summarise(HOM_E = mean(HOM_E))
  #
  #
  # freq.hom.pop <- hs %>%
  #   dplyr::group_by(POP_ID) %>%
  #   dplyr::summarise(HOM_E = mean(HOM_E))

  split.vec <- NULL
  # diplo <- split.vec.locus <- split.vec <- NULL

  # IBDg with FH ---------------------------------------------------------------
  # ok so here we need the info by id

  # fh.i <- suppressWarnings(
  #   summary.ind %>%
  #     dplyr::full_join(freq.hom.pop, by = "POP_ID") %>%
  #     dplyr::mutate(FH = ((HOM_O - HOM_E)/(N_GENOT - HOM_E))) %>%
  #     dplyr::select(INDIVIDUALS, POP_ID, FH)
  # )

  res$individual.summary <- suppressWarnings(
    summary.ind %>%
      dplyr::full_join(dplyr::select(summary.pop, -HOM_O, -HET_O), by = "POP_ID") %>%
      dplyr::mutate(FH = ((HOM_O - HOM_E)/(N_GENOT - HOM_E))) %>%
      dplyr::select(INDIVIDUALS, POP_ID, FH)
  )
  summary.ind <- NULL

  summary.pop <- summary.pop %>%
    dplyr::full_join(
      dplyr::group_by(.data = res$individual.summary, POP_ID) %>% dplyr::summarise(FH = mean(FH, na.rm = TRUE))
      , by = "POP_ID") %>%
    dplyr::mutate(
      FIS = dplyr::if_else(
        HET_O == 0, 0, round(((HET_E - HET_O) / HET_E), 6)),
      GIS = dplyr::if_else(
        HET_O == 0, 0, round(1 - HET_O / HS, 6))
    ) %>%
    dplyr::select(POP_ID, HOM_O, HOM_E, HET_O, HET_E, HS, FIS, GIS, FH)

  summary.locus <- summary.locus %>%
    dplyr::mutate(
      FIS = dplyr::if_else(
        HET_O == 0, 0, round(((HET_E - HET_O) / HET_E), 6)),
      GIS = dplyr::if_else(
        HET_O == 0, 0, round(1 - HET_O / HS, 6))
    )

  summary.locus <- NULL

  summary.overall <- suppressWarnings(dplyr::full_join(
    sample.number,
    dplyr::bind_rows(
      summary.pop,
      dplyr::summarise_if(
        .tbl = summary.pop, .predicate = is.numeric, .funs = mean, na.rm = TRUE) %>%
        tibble::add_column(.data = ., POP_ID = "OVERALL", .before = 1)
    ), by = "POP_ID"))

  sample.number <- summary.pop <- NULL
  # fh.tot <- fh.i %>%
  #   dplyr::summarise(
  #     HOM_O = round(mean(HOM_O), 6),
  #     HOM_E = round(mean(HOM_E), 6),
  #     HET_O = round(mean(1 - HOM_O), 6),
  #     HET_E = round(mean(1 - HOM_E), 6),
  #     FIS = ifelse(HET_O == 0, 0, round(((HET_E - HET_O) / HET_E), 6)),
  #     FH = mean(FH)
  #   )

  # fh.tot <- tibble::data_frame(POP_ID = "OVERALL") %>%
  #   dplyr::bind_cols(fh.tot)

  # fh.res <- suppressWarnings(dplyr::bind_rows(fh.pop, fh.tot) %>% dplyr::select(-POP_ID))
  res$summary <- suppressWarnings(
    dplyr::bind_rows(
      haplotype.filtered.split %>%
        dplyr::distinct(LOCUS, POP_ID, ALLELES) %>%
        dplyr::group_by(LOCUS, POP_ID) %>%
        dplyr::summarise(ALLELES_COUNT = length(ALLELES)) %>%
        dplyr::group_by(POP_ID) %>%
        dplyr::summarise(
          MONOMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT == 1]),
          POLYMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT >= 2])
        ) %>%
        dplyr::full_join(consensus.artifacts, by = "POP_ID"),
      haplotype.filtered.split %>%
        dplyr::distinct(LOCUS, ALLELES) %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::summarise(ALLELES_COUNT_OVERALL = length(ALLELES)) %>%
        dplyr::ungroup(.) %>%
        dplyr::summarise(
          MONOMORPHIC = length(ALLELES_COUNT_OVERALL[ALLELES_COUNT_OVERALL == 1]),
          POLYMORPHIC = length(ALLELES_COUNT_OVERALL[ALLELES_COUNT_OVERALL >= 2])
        ) %>%
        dplyr::bind_cols(blacklist.loci.consensus.sum) %>%
        dplyr::bind_cols(blacklist.loci.artifacts.sum) %>%
        tibble::add_column(.data = ., POP_ID = "OVERALL", .before = 1)))

  haplotype.filtered.split <- NULL

  if (keep.consensus) {
    res$summary <- res$summary %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::mutate(
        TOTAL = MONOMORPHIC + POLYMORPHIC + CONSENSUS + ARTIFACTS,
        MONOMORPHIC_PROP = round(MONOMORPHIC / TOTAL, 4),
        POLYMORPHIC_PROP = round(POLYMORPHIC / TOTAL, 4),
        CONSENSUS_PROP = round(CONSENSUS / TOTAL, 4),
        ARTIFACTS_PROP = round(ARTIFACTS / TOTAL, 4)
      )
  } else {
    res$summary <- res$summary %>%
      dplyr::select(-CONSENSUS) %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::mutate(
        TOTAL = MONOMORPHIC + POLYMORPHIC + ARTIFACTS,
        MONOMORPHIC_PROP = round(MONOMORPHIC / TOTAL, 4),
        POLYMORPHIC_PROP = round(POLYMORPHIC / TOTAL, 4),
        ARTIFACTS_PROP = round(ARTIFACTS / TOTAL, 4)
      )
  }
  res$summary <- suppressWarnings(
    dplyr::full_join(res$summary, summary.overall, by = "POP_ID"))

  consensus.artifacts <- blacklist.loci.consensus.sum <- blacklist.loci.artifacts.sum <- summary.overall <- NULL

  # Nei & Li 1979 Nucleotide Diversity -----------------------------------------
  message("Nucleotide diversity (Pi):")
  message("    Read length used: ", read.length)

  # Pi: by individuals----------------------------------------------------------
  message("    Pi calculations: individuals...")
  if (keep.consensus) {
    pi.data <- dplyr::filter(
      haplotype,
      !is.na(HAPLOTYPES),
      WHITELIST == "whitelist" | !POLYMORPHISM %in% "artifact")
  } else {
    pi.data <- dplyr::filter(
      haplotype,
      !is.na(HAPLOTYPES),
      WHITELIST == "whitelist",
      !POLYMORPHISM %in% c("artifact", "consensus"))
  }

  separate_haplo <- function(x) {
    res <- dplyr::select(x, LOCUS, POP_ID, INDIVIDUALS, HAPLOTYPES)  %>%
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"),
        sep = "/", extra = "drop", remove = TRUE
      ) %>%
      dplyr::mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2))
    return(res)
  }#End separate_haplo

  n.row <- nrow(pi.data)
  split.vec <- as.integer(floor((parallel.core * 20 * (1:n.row - 1) / n.row) + 1))
  n.row <- NULL

  pi.data <- split(x = pi.data, f = split.vec) %>%
    .stackr_parallel_mc(
      X = ., FUN = separate_haplo, mc.cores = parallel.core) %>%
    dplyr::bind_rows(.)

  split.vec <- NULL

  res$individual.summary <- dplyr::full_join(
    res$individual.summary,
    pi.data %>%
      dplyr::mutate(
        PI = (stringdist::stringdist(a = ALLELE1, b = ALLELE2, method = "hamming")) / read.length
      ) %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::summarise(PI = mean(PI, na.rm = TRUE))
    , by = "INDIVIDUALS")

  # Pi function ----------------------------------------------------------------

  pi <- function(data, read.length) {
    # data <- dplyr::filter(data, LOCUS == "2")#test
    y <- dplyr::select(data, ALLELES) %>% purrr::flatten_chr(.)

    if (length(unique(y)) <= 1) {
      pi <- tibble::data_frame(PI = as.numeric(0))
    } else {

      #1 Get all pairwise comparison
      allele.pairwise <- utils::combn(unique(y), 2)

      #2 Calculate pairwise nucleotide mismatches
      pairwise.mismatches <- apply(allele.pairwise, 2, function(z) {
        stringdist::stringdist(a = z[1], b = z[2], method = "hamming")
      })

      #3 Calculate allele frequency
      allele.freq <- table(y) / length(y)

      #4 Calculate nucleotide diversity from pairwise mismatches and allele frequency
      pi <- apply(allele.pairwise, 2, function(y) allele.freq[y[1]] * allele.freq[y[2]])
      pi <- tibble::data_frame(PI = sum(pi * pairwise.mismatches) / read.length)
    }
    return(pi)
  }#End pi

  pi_pop <- function(data, read.length, parallel.core) {
    pop <- unique(data$POP_ID)
    message("    Pi calculations for pop: ", pop)
    # data <- df.split.pop[["DD"]]
    # data <- df.split.pop[["SKY"]]

    pi.pop <- split(x = data, f = data$LOCUS)
    pi.pop <- .stackr_parallel(
      X = pi.pop,
      FUN = pi,
      mc.cores = parallel.core,
      read.length = read.length
    ) %>%
      dplyr::bind_rows(.) %>%
      dplyr::summarise(PI_NEI = mean(PI, na.rm = TRUE)) %>%
      tibble::add_column(.data = ., POP_ID = pop, .before = "PI_NEI")

    return(pi.pop)
  }#End pi_pop

  # Pi: by pop------------------------------------------------------------------
  message("    Pi calculations: populations...")
  pi.data <- pi.data %>%
    tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, INDIVIDUALS, POP_ID))

  pi.pop <- pi.data %>%
    split(x = ., f = .$POP_ID) %>%
    purrr::map_df(
      .x = ., .f = pi_pop,
      read.length = read.length, parallel.core = parallel.core
    )

  # Pi: overall  ---------------------------------------------------------------
  message("    Pi calculations: overall")
  pi.overall <- split(x = pi.data, f = pi.data$LOCUS)
  pi.data <- NULL
  pi.overall <- .stackr_parallel(
    X = pi.overall,
    FUN = pi,
    mc.cores = parallel.core,
    read.length = read.length
  ) %>%
    dplyr::bind_rows(.) %>%
    dplyr::summarise(PI_NEI = mean(PI, na.rm = TRUE)) %>%
    tibble::add_column(.data = ., POP_ID = "OVERALL", .before = "PI_NEI")

  # Combine the pop and overall data -------------------------------------------
  pi.overall <- suppressWarnings(
    dplyr::bind_rows(pi.pop, pi.overall))
  pi.pop <- NULL
  res$summary <- suppressWarnings(dplyr::full_join(res$summary, pi.overall, by = "POP_ID"))
  pi.overall <- NULL

  if (is.null(whitelist.markers)) {
    filename.sum <- "haplotype.catalog.loci.summary.pop.tsv"
  } else {
    filename.sum <- "haplotype.catalog.loci.summary.pop.filtered.tsv"
  }
  readr::write_tsv(x = res$summary, path = stringi::stri_join(path.folder, "/",filename.sum))

  # Figures --------------------------------------------------------------------
  res$scatter.plot.Pi.Fh.ind <- dplyr::filter(res$individual.summary, POP_ID != "OVERALL") %>%
    ggplot2::ggplot(data = ., ggplot2::aes(x = FH, y = PI)) +
    ggplot2::geom_point(ggplot2::aes(colour = POP_ID)) +
    ggplot2::stat_smooth(method = stats::lm, level = 0.95, fullrange = FALSE, na.rm = TRUE) +
    ggplot2::labs(x = "Individual IBDg (FH)") +
    ggplot2::labs(y = "Individual nucleotide diversity (Pi)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  ggplot2::ggsave(
    filename = stringi::stri_join(path.folder, "/scatter.plot.Pi.Fh.ind.pdf"),
    plot = res$scatter.plot.Pi.Fh.ind,
    width = 20, height = 15,
    dpi = 600, units = "cm", useDingbats = FALSE)

  res$scatter.plot.Pi.Fh.pop <- dplyr::filter(res$summary, POP_ID != "OVERALL") %>%
    ggplot2::ggplot(data = ., ggplot2::aes(x = FH, y = PI_NEI)) +
    ggplot2::geom_point(ggplot2::aes(colour = POP_ID)) +
    ggplot2::stat_smooth(method = stats::lm, level = 0.95, fullrange = FALSE, na.rm = TRUE) +
    ggplot2::labs(x = "Individual IBDg (FH)") +
    ggplot2::labs(y = "Individual nucleotide diversity (Pi)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  ggplot2::ggsave(
    filename = stringi::stri_join(path.folder, "/scatter.plot.Pi.Fh.pop.pdf"),
    plot = res$scatter.plot.Pi.Fh.pop,
    width = 20, height = 15,
    dpi = 600, units = "cm", useDingbats = FALSE)


  res$boxplot.pi <- dplyr::filter(res$individual.summary, POP_ID != "OVERALL") %>%
    ggplot2::ggplot(data = ., ggplot2::aes(x = factor(POP_ID), y = PI, na.rm = TRUE)) +
    ggplot2::geom_violin(trim = F) +
    ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    ggplot2::labs(x = "Sampling sites") +
    ggplot2::labs(y = "Individual nucleotide diversity (Pi)") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  ggplot2::ggsave(
    filename = stringi::stri_join(path.folder, "/boxplot.pi.pdf"),
    plot = res$boxplot.pi,
    width = 20, height = 15,
    dpi = 600, units = "cm", useDingbats = FALSE)

  res$boxplot.fh <- dplyr::filter(res$individual.summary, POP_ID != "OVERALL") %>%
    ggplot2::ggplot(data = ., ggplot2::aes(x = factor(POP_ID), y = FH, na.rm = T)) +
    ggplot2::geom_violin(trim = F) +
    ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    ggplot2::labs(x = "Sampling sites") +
    ggplot2::labs(y = "Individual IBDg (FH)") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  ggplot2::ggsave(
    filename = stringi::stri_join(path.folder, "/boxplot.fh.pdf"),
    plot = res$boxplot.fh,
    width = 20, height = 15,
    dpi = 600, units = "cm", useDingbats = FALSE)

  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message("Number of populations: ", length(pop))
  message("Number of individuals: ", n.individuals)
  message("Number of loci in the catalog: ", n.markers)
  if (!is.null(res$artifacts.pop)) {
    message("    impacted by assembly artifacts and/or sequencing errors (loci > 2 alleles) = ", dplyr::n_distinct(res$artifacts.pop$LOCUS))
  }
  message("    with consensus alleles = ", dplyr::n_distinct(res$consensus.pop$LOCUS))
  if (!is.null(res$artifacts.pop)) {
    message("Number of artefactual genotypes/total genotypes (percent): ", erased.genotype.number, "/", total.genotype.number.haplo, " (", percent.haplo, ")")
  }
  message("\nFiles written in: ")
  if (!is.null(res$artifacts.pop)) message(filename.artifacts.ind)
  if (!is.null(res$artifacts.pop)) message(filename.paralogs)
  message(filename.sum)
  if (!is.null(res$consensus.loci)) message("blacklist.loci.consensus.txt")
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  options(width = opt.change)
  return(res)
}
