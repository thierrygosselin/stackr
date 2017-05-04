## Summary and tables

#' @title Haplotypes file summary
#' @description STACKS batch_x.haplotypes.tsv file summary.
#' The output of the function is a summary table for populations with:
#' \enumerate{
#' \item assembly artifacts and/or sequencing errors: loci with > 2 alleles.
#' Genotypes with more than 2 alleles are erased before estimating the subsequent statistics.
#' \item consensus loci
#' \item monomorphic loci
#' \item polymorphic loci
#' \item haplotypes statistics for the observed and expected homozygosity and
#' heterozygosity.
#' \item Wright’s inbreeding coefficient (Fis)
#' \item \strong{IBDG}: a proxy measure of the realized proportion of the genome
#' that is identical by descent
#' \item \strong{FH measure}: based on the excess in the observed number of homozygous
#' genotypes within an individual relative to the mean number of homozygous
#' genotypes expected under random mating (Keller et al., 2011; Kardos et al., 2015).
#' \item \strong{Pi}: the nucleotide diversity measured here consider the
#' consensus loci in the catalog (no variation between population sequences).
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

#' @inheritParams tidy_genomic_data

#' @param parallel.core (optional) The number of core for parallel
#' programming during Pi calculations.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @importFrom stringdist stringdist
#' @importFrom utils combn count.fields
#' @importFrom stats lm
#' @importFrom stringi stri_replace_all_fixed stri_replace_na stri_join stri_count_fixed
#' @importFrom tibble as_data_frame data_frame add_column
#' @importFrom dplyr select rename n_distinct distinct mutate summarise group_by ungroup arrange left_join full_join semi_join anti_join bind_rows bind_cols if_else
#' @importFrom readr write_tsv read_tsv
#' @importFrom tidyr separate gather
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid stat_smooth
#' @importFrom purrr flatten_chr map_df

#' @return The function returns a list with:
#' \enumerate{
#' \item $summary
#' \item $artifacts.pop # if artifacts are found
#' \item $artifacts.loci # if artifacts are found
#' \item $consensus.pop # if consensus locus are found
#' \item $consensus.loci # if consensus locus are found
#' \item $fh.pi.individuals: the individual's info for FH and Pi
#' }
#' Also in the list 3 plots:
#' \enumerate{
#' \item $scatter.plot
#' \item $boxplot.pi
#' \item $boxplot.fh
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
  pop.levels = NULL,
  pop.labels = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  cat("#######################################################################\n")
  cat("#################### stackr::summary_haplotypes #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

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

  # Import haplotype file ------------------------------------------------------
  message("Importing STACKS haplotypes file: ", data)
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
  if (!is.null(whitelist.markers)) { # no Whitelist
    whitelist.markers <- suppressMessages(readr::read_tsv(whitelist.markers, col_names = TRUE))
    columns.names.whitelist <- colnames(whitelist.markers)
    if ("CHROM" %in% columns.names.whitelist) {
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
    if ("LOCUS" %in% columns.names.whitelist) {
      whitelist.markers$LOCUS <- as.character(whitelist.markers$LOCUS)
    }
    if ("POS" %in% columns.names.whitelist) {
      whitelist.markers$POS <- as.character(whitelist.markers$POS)
    }
    # haplo.file
    whitelist.markers <- dplyr::select(.data = whitelist.markers, LOCUS) %>%
      dplyr::arrange(as.integer(LOCUS))
    columns.names.whitelist <- colnames(whitelist.markers)
  }

  # Import blacklist id --------------------------------------------------------
  if (!is.null(blacklist.id)) { # No blacklist of ID
    blacklist.id <- suppressMessages(readr::read_tsv(blacklist.id, col_names = TRUE))
  }

  # population levels and strata------------------------------------------------
  message("Making STACKS haplotypes file population-wise...")
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

  # filtering the strata if blacklist id available
  if (!is.null(blacklist.id)) {
    strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
  }
  # Remove potential whitespace in pop_id
  strata.df$POP_ID <- stringi::stri_replace_all_fixed(strata.df$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)

  # # Check with strata and pop.levels/pop.labels
  # if (!is.null(strata) & !is.null(pop.levels)) {
  #   if (length(levels(factor(strata.df$POP_ID))) != length(pop.levels)) {
  #     stop("The number of groups in your strata file must match the number of groups in pop.levels")
  #   }
  # }
  # Filter with whitelist of markers -------------------------------------------
  if (!is.null(whitelist.markers)) {
    message("Filtering with whitelist of loci: ", nrow(whitelist.markers), " loci")
    haplotype <- suppressWarnings(dplyr::semi_join(haplotype, whitelist.markers, by = columns.names.whitelist))
  }

  # Filter with blacklist of individuals ---------------------------------------
  if (!is.null(blacklist.id)) {
    message("Filtering with blacklist of individuals: ", nrow(blacklist.id), " individual(s)")

    blacklist.id$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = blacklist.id$INDIVIDUALS,
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )

    haplotype <- suppressWarnings(dplyr::anti_join(haplotype, blacklist.id, by = "INDIVIDUALS"))
  }

  # Population levels and strata -----------------------------------------------
  strata.df$INDIVIDUALS = stringi::stri_replace_all_fixed(
    str = strata.df$INDIVIDUALS,
    pattern = c("_", ":"),
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )

  haplotype <- dplyr::left_join(x = haplotype, y = strata.df, by = "INDIVIDUALS")

  # using pop.levels and pop.labels info if present
  haplotype <- change_pop_names(data = haplotype, pop.levels = pop.levels, pop.labels = pop.labels)

  # Locus with consensus alleles -----------------------------------------------
  message("Scanning for consensus markers...")

  # consensus.pop <- haplotype %>%
  #   dplyr::mutate(CONSENSUS = stringi::stri_count_fixed(HAPLOTYPES, "consensus")) %>%
  #   dplyr::group_by(LOCUS, POP_ID) %>%
  #   dplyr::summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
  #   dplyr::filter(CONSENSUS_MAX > 0)

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
    message("    Generated a file with ", nrow(blacklist.loci.consensus), " consensus loci: blacklist.loci.consensus.txt")
    readr::write_tsv(
      x = blacklist.loci.consensus,
      path = "blacklist.loci.consensus.txt",
      col_names = TRUE
    )

    blacklist.loci.consensus.sum <- blacklist.loci.consensus %>%
      dplyr::ungroup(.) %>%
      dplyr::summarise(CONSENSUS = n()) %>%
      dplyr::select(CONSENSUS)

    consensus.pop.sum <- consensus.pop %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::summarise(CONSENSUS = dplyr::n_distinct(LOCUS))
  } else {
    blacklist.loci.consensus <- NULL
    blacklist.loci.consensus.sum <- tibble::data_frame(CONSENSUS = as.integer(0))
    consensus.pop.sum <- tibble::data_frame(
      POP_ID = pop.labels,
      CONSENSUS = as.integer(rep(0, length(pop.labels)))
    )
  }

  # Individuals per pop
  ind.pop <- haplotype %>%
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE) %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::tally(.) %>%
    dplyr::rename(NUM = n)

  ind.tot <- haplotype %>%
    dplyr::distinct(INDIVIDUALS) %>%
    dplyr::tally(.) %>%
    dplyr::rename(NUM = n) %>%
    dplyr::mutate(POP_ID = "OVERALL")

  sample.number <- suppressWarnings(dplyr::bind_rows(ind.pop, ind.tot))

  # unused arguments
  ind.pop <- ind.tot <- NULL

  # Sequencing errors and assembly artifact ---------------------------------------------
  # Locus with > 2 alleles by individuals
  message("Scanning for outlier genotypes with > 2 alleles (assembly artifact and sequencing errors)...")

  total.genotype.number.haplo <- length(haplotype$HAPLOTYPES)

  paralogs.ind <- haplotype %>%
    dplyr::filter(stringi::stri_count_fixed(HAPLOTYPES, "/") > 1) %>%
    dplyr::arrange(LOCUS, POP_ID, INDIVIDUALS) %>%
    dplyr::select(LOCUS, POP_ID, INDIVIDUALS, HAPLOTYPES)

  if (nrow(paralogs.ind) > 0) {
    erased.genotype.number <- length(paralogs.ind$HAPLOTYPES)
    percent.haplo <- paste(round(((erased.genotype.number/total.genotype.number.haplo)*100), 2), "%", sep = " ")

    # Write the list of locus, individuals with paralogs
    if (is.null(whitelist.markers)) {
      readr::write_tsv(x = paralogs.ind, path = "blacklist.loci.artifacts.ind.txt", col_names = TRUE)
      filename.paralogs.ind <- "blacklist.loci.artifacts.ind.txt"
    } else {
      readr::write_tsv(x = paralogs.ind, path = "blacklist.loci.artifacts.ind.filtered.txt", col_names = TRUE)
      filename.paralogs.ind <- "blacklist.loci.artifacts.ind.filtered.txt"
    }
    message("    Generated a file with ", nrow(paralogs.ind), " artifact loci by individuals: ", filename.paralogs.ind)

    paralogs.pop <- paralogs.ind %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(LOCUS, POP_ID) %>%
      dplyr::arrange(as.numeric(LOCUS)) %>%
      dplyr::mutate(PARALOGS = rep("paralogs", times = n()))

    paralogs.pop.sum <- paralogs.pop %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::summarise(ARTIFACTS = dplyr::n_distinct(LOCUS))

    blacklist.loci.paralogs <- paralogs.ind %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(LOCUS) %>%
      dplyr::arrange(as.numeric(LOCUS))

    blacklist.loci.paralogs.sum <- blacklist.loci.paralogs %>%
      dplyr::ungroup(.) %>%
      dplyr::summarise(ARTIFACTS = n()) %>%
      dplyr::select(ARTIFACTS)

    # Write the unique list of paralogs blacklisted to a file
    if (is.null(whitelist.markers)) {
      readr::write_tsv(
        x = blacklist.loci.paralogs,
        path = "blacklist.loci.artifacts.txt", col_names = TRUE
      )
      filename.paralogs <- "blacklist.loci.artifacts.txt"
    } else {
      readr::write_tsv(
        x = blacklist.loci.paralogs,
        path = "blacklist.loci.artifacts.filtered.txt", col_names = TRUE
      )
      filename.paralogs <- "blacklist.loci.artifacts.filtered.txt"
    }
    message("    Generated a file with ", nrow(blacklist.loci.paralogs), " artifact loci: ", filename.paralogs)
  } else {
    paralogs.pop <- NULL
    blacklist.loci.paralogs <- NULL
    pop <- unique(strata.df$POP_ID)
    paralogs.pop.sum <- tibble::data_frame(POP_ID = pop, ARTIFACTS = as.integer(rep(0, length(pop))))
    blacklist.loci.paralogs.sum <- tibble::data_frame(ARTIFACTS = as.integer(0))
  }

  # Haplo filtered paralogs ----------------------------------------------------
  if (length(paralogs.ind$HAPLOTYPES > 0)) {
    message("    Erasing genotypes with > 2 alleles")
    message.haplo.erasing.geno <- stringi::stri_join("    Out of a total of ", total.genotype.number.haplo, " genotypes, ", percent.haplo, " (", erased.genotype.number, ")"," will be erased")
    message(message.haplo.erasing.geno)

    erase.paralogs <- paralogs.ind %>%
      dplyr::select(LOCUS, INDIVIDUALS) %>%
      dplyr::mutate(ERASE = rep("erase", n()))

    haplo.filtered.paralogs <- suppressWarnings(
      haplotype %>%
        dplyr::full_join(erase.paralogs, by = c("LOCUS", "INDIVIDUALS")) %>%
        dplyr::mutate(
          ERASE = stringi::stri_replace_na(str = ERASE, replacement = "ok"),
          HAPLOTYPES = ifelse(ERASE == "erase", NA, HAPLOTYPES)
        ) %>%
        dplyr::select(-ERASE)
    )
  } else {
    haplo.filtered.paralogs <- haplotype
  }

  # Haplo filtered for consensus -----------------------------------------------
  if (length(consensus.pop$CONSENSUS) > 0) {
    haplo.filtered.consensus <- haplotype %>%
      dplyr::filter(!LOCUS %in% consensus.pop$LOCUS)

    # Haplo filtered for consensus and paralogs
    haplo.filtered.consensus.paralogs <- haplo.filtered.paralogs %>%
      dplyr::filter(!LOCUS %in% consensus.pop$LOCUS)
  } else {
    haplo.filtered.consensus <- haplotype
    haplo.filtered.consensus.paralogs <- haplo.filtered.paralogs
  }

  # Summary dataframe by individual---------------------------------------------
  message("Genome-Wide Identity-By-Descent calculations (FH)...")
  n.individuals <- dplyr::n_distinct(haplo.filtered.consensus.paralogs$INDIVIDUALS)

  if (n.individuals < 200) {
    summary.ind <- haplo.filtered.consensus.paralogs %>%
      dplyr::mutate(
        ALLELES_COUNT = stringi::stri_count_fixed(HAPLOTYPES, "/"),
        IND_LEVEL_POLYMORPHISM = dplyr::if_else(ALLELES_COUNT == 1, "het", "hom", missing = "missing")
      ) %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::summarise(
        HOM = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "hom"]),
        HET = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "het"]),
        MISSING = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "missing"]),
        N_GENOT = HOM + HET,
        HOM_O = HOM/N_GENOT,
        HET_O = HET/N_GENOT
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::left_join(strata.df, by = "INDIVIDUALS") %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)
  } else {
    individual_stats <- function(x) {
      res <- dplyr::group_by(x, INDIVIDUALS) %>%
        dplyr::summarise(
          HOM = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "hom"]),
          HET = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "het"]),
          MISSING = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "missing"]),
          N_GENOT = HOM + HET,
          HOM_O = HOM/N_GENOT,
          HET_O = HET/N_GENOT
        ) %>%
        dplyr::ungroup(.)
      return(res)
    }#End individual_stats

    split.vec <- dplyr::distinct(strata.df, INDIVIDUALS) %>%
      dplyr::mutate(SPLIT_VEC = as.integer(floor((parallel.core * 3 * (1:n.individuals - 1) / n.individuals) + 1)))

    summary.ind <- haplo.filtered.consensus.paralogs %>%
      dplyr::mutate(
        ALLELES_COUNT = stringi::stri_count_fixed(HAPLOTYPES, "/"),
        IND_LEVEL_POLYMORPHISM = dplyr::if_else(ALLELES_COUNT == 1, "het", "hom", missing = "missing")
      ) %>%
      dplyr::left_join(split.vec, by = "INDIVIDUALS") %>%
      split(x = ., f = .$SPLIT_VEC) %>%
      .stackr_parallel(
        X = .,
        FUN = individual_stats,
        mc.cores = parallel.core
      ) %>%
      dplyr::bind_rows(.) %>%
      dplyr::left_join(strata.df, by = "INDIVIDUALS") %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)

    split.vec <- NULL
  }
    # system.time(freq.alleles.loci.pop <- suppressWarnings(
  #   haplo.filtered.consensus.paralogs %>%
  #     dplyr::filter(!is.na(HAPLOTYPES)) %>%
  #     dplyr::group_by(LOCUS, POP_ID) %>%
  #     dplyr::mutate(DIPLO = length(INDIVIDUALS) * 2) %>%
  #     tidyr::separate(
  #       col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"),
  #       sep = "/", extra = "drop", remove = FALSE
  #     ) %>%
  #     dplyr::mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
  #     dplyr::select(-HAPLOTYPES, -INDIVIDUALS) %>%
  #     tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID, DIPLO)) %>%
  #     dplyr::group_by(LOCUS, POP_ID, ALLELES) %>%
  #     dplyr::summarise(
  #       FREQ_ALLELES = length(ALLELES)/mean(DIPLO),
  #       HOM_E = FREQ_ALLELES * FREQ_ALLELES
  #     ) %>%
  #     dplyr::select(-FREQ_ALLELES) %>% dplyr::arrange(LOCUS, POP_ID) %>%
  #     dplyr::ungroup(.)
  # ))

  # faster in parallel...
  # required functions
  separate_haplo <- function(x) {
    res <- x %>%
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"),
        sep = "/", extra = "drop", remove = FALSE
      ) %>%
      dplyr::mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
      dplyr::select(-HAPLOTYPES, -INDIVIDUALS) %>%
      tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID)) %>%
      dplyr::select(-ALLELE_GROUP)
    return(res)
  }#End separate_haplo
  freq_hom <- function(x) {
    res <- dplyr::select(x, -SPLIT_VEC) %>%
      dplyr::group_by(LOCUS, POP_ID, ALLELES) %>%
      dplyr::summarise(
        FREQ_ALLELES = length(ALLELES)/mean(DIPLO),
        HOM_E = FREQ_ALLELES * FREQ_ALLELES
      ) %>%
      dplyr::select(-FREQ_ALLELES) %>%
      dplyr::ungroup(.)
    return(res)
  }#End freq_hom

  freq.alleles.loci.pop <- haplo.filtered.consensus.paralogs %>%
    dplyr::filter(!is.na(HAPLOTYPES))
  diplo <- freq.alleles.loci.pop %>%
    dplyr::distinct(LOCUS, POP_ID, INDIVIDUALS) %>%
    dplyr::group_by(LOCUS, POP_ID) %>%
    dplyr::tally(.) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(n = 2 * n) %>%
    dplyr::rename(DIPLO = n)
  n.row <- nrow(freq.alleles.loci.pop)
  split.vec <- as.integer(floor((parallel.core * 20 * (1:n.row - 1) / n.row) + 1))
  n.row <- NULL
  split.vec.locus <- dplyr::distinct(freq.alleles.loci.pop, LOCUS) %>%
    dplyr::mutate(
      SPLIT_VEC = factor(as.integer(floor((parallel.core * 20 * (1:n() - 1) / n()) + 1)))
    )
  freq.alleles.loci.pop <- split(x = freq.alleles.loci.pop, f = split.vec) %>%
    .stackr_parallel_mc(
      X = ., FUN = separate_haplo, mc.cores = parallel.core) %>%
    dplyr::bind_rows(.) %>%
    dplyr::left_join(diplo, by = c("LOCUS", "POP_ID")) %>%
    dplyr::left_join(split.vec.locus, by = "LOCUS") %>%
    split(x = ., f = .$SPLIT_VEC) %>%
    .stackr_parallel_mc(
      X = ., FUN = freq_hom, mc.cores = parallel.core) %>%
    dplyr::bind_rows(.)
  diplo <- split.vec.locus <- split.vec <- NULL

  freq.loci.pop <- freq.alleles.loci.pop %>%
    dplyr::group_by(LOCUS, POP_ID) %>%
    dplyr::summarise(HOM_E = sum(HOM_E))

  freq.pop <- freq.loci.pop %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(HOM_E = mean(HOM_E))


  # IBDg with FH ---------------------------------------------------------------
  fh.i <- suppressWarnings(
    summary.ind %>%
      dplyr::full_join(freq.pop, by = "POP_ID") %>%
      dplyr::mutate(FH = ((HOM_O - HOM_E)/(N_GENOT - HOM_E)))
  )
  fh.i.out <- dplyr::select(.data = fh.i, INDIVIDUALS, POP_ID, FH)

  fh.pop <- fh.i %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      HOM_O = round(mean(HOM_O), 6),
      HOM_E = round(mean(HOM_E), 6),
      HET_O = round(mean(1 - HOM_O), 6),
      HET_E = round(mean(1 - HOM_E), 6),
      FIS = ifelse(HET_O == 0, 0, round(((HET_E - HET_O) / HET_E), 6)),
      FH = mean(FH)
    )

  fh.tot <- fh.i %>%
    dplyr::summarise(
      HOM_O = round(mean(HOM_O), 6),
      HOM_E = round(mean(HOM_E), 6),
      HET_O = round(mean(1 - HOM_O), 6),
      HET_E = round(mean(1 - HOM_E), 6),
      FIS = ifelse(HET_O == 0, 0, round(((HET_E - HET_O) / HET_E), 6)),
      FH = mean(FH)
    )

  fh.tot <- tibble::data_frame(POP_ID = "OVERALL") %>%
    dplyr::bind_cols(fh.tot)

  fh.res <- suppressWarnings(dplyr::bind_rows(fh.pop, fh.tot) %>% dplyr::select(-POP_ID))

  # unused arguments
  freq.pop <- summary.ind <- fh.tot <- fh.pop <- NULL

  # Nei & Li 1979 Nucleotide Diversity -----------------------------------------
  message("Nucleotide diversity (Pi):")
  message("    Read length used: ", read.length)

  separate_haplo <- function(x) {
    res <- x %>%
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"),
        sep = "/", extra = "drop", remove = FALSE
      ) %>%
      dplyr::mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2))
  }#End separate_haplo

  pi.data <- dplyr::filter(haplo.filtered.paralogs, !is.na(HAPLOTYPES))
  n.row <- nrow(pi.data)
  split.vec <- as.integer(floor((parallel.core * 20 * (1:n.row - 1) / n.row) + 1))
  n.row <- NULL

  pi.data <- split(x = pi.data, f = split.vec) %>%
    .stackr_parallel_mc(
      X = ., FUN = separate_haplo, mc.cores = parallel.core) %>%
    dplyr::bind_rows(.)

  # Pi: by individuals----------------------------------------------------------
  message("    Pi calculations: individuals...")

  pi.data.i <- pi.data %>%
    dplyr::mutate(
      PI = (stringdist::stringdist(a = ALLELE1, b = ALLELE2, method = "hamming"))/read.length
    ) %>%
    dplyr::group_by(INDIVIDUALS) %>%
    dplyr::summarise(PI = mean(PI))

  pi.i.out <- dplyr::select(.data = pi.data.i, INDIVIDUALS, PI)

  fh.pi.individuals <- dplyr::full_join(fh.i.out, pi.i.out, by = "INDIVIDUALS")

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
      allele.freq <- table(y)/length(y)

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
      mc.cores = 8,
      read.length = read.length
    ) %>%
      dplyr::bind_rows(.) %>%
      dplyr::summarise(PI_NEI = mean(PI)) %>%
      tibble::add_column(.data = ., POP_ID = pop, .before = "PI_NEI")

    return(pi.pop)
  }#End pi_pop

  # Pi: by pop------------------------------------------------------------------
  message("    Pi calculations: populations:")
  pi.data.pop <- pi.data %>%
    tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, INDIVIDUALS, POP_ID))

  pi.pop <- pi.data.pop %>%
    split(x = ., f = .$POP_ID) %>%
    purrr::map_df(
      .x = ., .f = pi_pop,
      read.length = read.length, parallel.core = parallel.core
    )

  # Pi: overall  ---------------------------------------------------------------
  message("    Pi calculations: overall")
  pi.overall <- split(x = pi.data.pop, f = pi.data.pop$LOCUS)
  pi.overall <- .stackr_parallel(
    X = pi.overall,
    FUN = pi,
    mc.cores = parallel.core,
    read.length = read.length
  ) %>%
    dplyr::bind_rows(.) %>%
    dplyr::summarise(PI_NEI = mean(PI)) %>%
    tibble::add_column(.data = ., POP_ID = "OVERALL", .before = "PI_NEI")

  # Combine the pop and overall data -------------------------------------------
  pi.res <- suppressWarnings(dplyr::bind_rows(pi.pop, pi.overall) %>% dplyr::select(-POP_ID))

  pi.pop <- pi.overall <- NULL
  # Summary dataframe by pop ---------------------------------------------------
  message("Working on the summary table...")
  separate_haplo <- function(x) {
    res <- x %>%
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"),
        sep = "/", extra = "drop", remove = TRUE
      ) %>%
      dplyr::mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
      tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID))
  }#End separate_haplo

  summary.prep <- dplyr::filter(haplo.filtered.consensus, !is.na(HAPLOTYPES)) %>%
    dplyr::select(-INDIVIDUALS)
  n.row <- nrow(summary.prep)
  split.vec <- as.integer(floor((parallel.core * 3 * (1:n.row - 1) / n.row) + 1))
  n.row <- NULL

  summary.prep <- split(x = summary.prep, f = split.vec) %>%
    .stackr_parallel_mc(
      X = ., FUN = separate_haplo, mc.cores = parallel.core) %>%
    dplyr::bind_rows(.)
  split.vec <- NULL
  # system.time(summary.prep2 <- suppressWarnings(
  #   haplo.filtered.consensus %>%
  #     dplyr::filter(!is.na(HAPLOTYPES)) %>%
  #     dplyr::select(-INDIVIDUALS) %>%
  #     tidyr::separate(
  #       col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"),
  #       sep = "/", extra = "drop", remove = TRUE
  #     ) %>%
  #     dplyr::mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
  #     tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID))
  # ))

  summary.pop <- suppressWarnings(
    summary.prep %>%
      dplyr::group_by(LOCUS, POP_ID) %>%
      dplyr::distinct(ALLELES, .keep_all = TRUE) %>%
      dplyr::summarise(ALLELES_COUNT = length(ALLELES)) %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::summarise(
        MONOMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT == 1]),
        POLYMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT >= 2])
      ) %>%
      dplyr::full_join(consensus.pop.sum, by = "POP_ID") %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::mutate(TOTAL = MONOMORPHIC + POLYMORPHIC + CONSENSUS) %>%
      dplyr::full_join(paralogs.pop.sum, by = "POP_ID") %>%
      dplyr::mutate(
        MONOMORPHIC_PROP = round(MONOMORPHIC/TOTAL, 4),
        POLYMORPHIC_PROP = round(POLYMORPHIC/TOTAL, 4),
        CONSENSUS_PROP = round(CONSENSUS/TOTAL, 4),
        ARTIFACTS_PROP = round(ARTIFACTS/TOTAL, 4)
      ))

  total <- summary.prep %>%
    dplyr::group_by(LOCUS) %>%
    dplyr::distinct(ALLELES, .keep_all = TRUE) %>%
    dplyr::summarise(ALLELES_COUNT = length(ALLELES)) %>%
    dplyr::summarise(
      MONOMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT == 1]),
      POLYMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT >= 2])
    ) %>%
    dplyr::bind_cols(blacklist.loci.consensus.sum) %>%
    dplyr::mutate(TOTAL = MONOMORPHIC + POLYMORPHIC + CONSENSUS) %>%
    dplyr::bind_cols(blacklist.loci.paralogs.sum) %>%
    dplyr::mutate(
      MONOMORPHIC_PROP = round(MONOMORPHIC/TOTAL, 4),
      POLYMORPHIC_PROP = round(POLYMORPHIC/TOTAL, 4),
      CONSENSUS_PROP = round(CONSENSUS/TOTAL, 4),
      ARTIFACTS_PROP = round(ARTIFACTS/TOTAL, 4)
    )

  total.res <- tibble::data_frame(POP_ID = "OVERALL") %>%
    dplyr::bind_cols(total)

  summary <- suppressWarnings(dplyr::bind_rows(summary.pop, total.res) %>% dplyr::ungroup(.) %>% dplyr::select(-POP_ID))
  summary <- dplyr::bind_cols(sample.number, summary, fh.res, pi.res)

  if (is.null(whitelist.markers)) {
    readr::write_tsv(x = summary, path = "haplotype.catalog.loci.summary.pop.tsv")
    filename.sum <- "haplotype.catalog.loci.summary.pop.tsv"
  } else {
    readr::write_tsv(x = summary, path = "haplotype.catalog.loci.summary.pop.filtered.tsv")
    filename.sum <- "haplotype.catalog.loci.summary.pop.filtered.tsv"
  }


  # Figures --------------------------------------------------------------------
  fh.pi <- pi.data.i %>%
    dplyr::full_join(
      fh.i %>% dplyr::select(INDIVIDUALS, POP_ID, FH)
      , by = "INDIVIDUALS")

  scatter.plot <- ggplot2::ggplot(fh.pi, ggplot2::aes(x = FH, y = PI)) +
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


  boxplot.pi <-   ggplot2::ggplot(fh.pi, ggplot2::aes(x = factor(POP_ID), y = PI, na.rm = T)) +
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

  boxplot.fh <-   ggplot2::ggplot(fh.pi, ggplot2::aes(x = factor(POP_ID), y = FH, na.rm = T)) +
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
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  n.markers <- dplyr::n_distinct(haplotype$LOCUS)
  n.ind <- dplyr::n_distinct(haplotype$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(haplotype$POP_ID)
  message("Number of populations: ", n.pop)
  message("Number of individuals: ", n.ind)
  message("Number of loci: ", n.markers)
  if (!is.null(paralogs.pop)) {
    message("    number of loci in the catalog impacted by assembly artifacts and/or sequencing errors (loci > 2 alleles) = ", dplyr::n_distinct(paralogs.pop$LOCUS))
    message("    number of artefactual genotypes/total genotypes (percent): ", erased.genotype.number, "/", total.genotype.number.haplo, " (", percent.haplo, ")")
  }
  message("    number of loci in the catalog with consensus alleles = ", dplyr::n_distinct(consensus.pop$LOCUS))
  message("Files written in this directory: ", getwd())
  if (!is.null(paralogs.pop)) message(filename.paralogs.ind)
  if (!is.null(paralogs.pop)) message(filename.paralogs)
  message(filename.sum)
  if (!is.null(blacklist.loci.consensus)) message("blacklist.loci.consensus.txt")
  # Results
  results <- list(
    summary = summary,
    artifacts.pop = paralogs.pop,
    artifacts.loci = blacklist.loci.paralogs,
    consensus.pop = consensus.pop,
    consensus.loci = blacklist.loci.consensus,
    fh.pi.individuals = fh.pi.individuals,
    scatter.plot = scatter.plot,
    boxplot.pi = boxplot.pi,
    boxplot.fh = boxplot.fh
  )
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(results)
}
