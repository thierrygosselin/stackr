# Linkage map
# Genotype summary

#' @title Summary of \code{batch_x.genotypes.txt} and
#' \code{batch_x.markers.tsv} files.
#' @description Useful summary of \code{batch_x.genotypes.txt} and
#'  \code{batch_x.markers.tsv} files produced by STACKS genotypes module use
#'  for linkage mapping. Filter segregation distortion and output JoinMap and or
#'  OneMap file.
#' @param genotypes The \code{genotypes = batch_x.genotypes.txt} created by STACKS genotypes
#' module.
#' @param markers The \code{markers = batch_x.markers.tsv} created by STACKS genotypes
#' module.
#' @param filter.monomorphic (optional) Should monomorphic loci be filtered out.
#' Default: \code{filter.monomorphic = TRUE}.
#' @param filter.missing.band (optional) Should loci with missing
#' band be filtered out. Default: \code{filter.missing.band = TRUE}.
#' @param filter.mean.log.likelihood (optional, integer) Apply a mean log likelihood
#' filter to the loci. e.g. filter.mean.log.likelihood = -10.
#' Default: \code{filter.mean.log.likelihood = NULL} (no filter)
#' @param B (optional, integer) The segregation distortion p-value
#' will be computed with a Monte Carlo test with \code{B} replicates,
#' if \code{B} is supplied.
#' For more details, see \code{\link[stats]{chisq.test}}.
#' Default: \code{B = NULL}.
#' @param filter.GOF (optional, integer) Filer value of the goodness-of-fit for segregation
#' distortion. Default: \code{filter.GOF = NULL}.
#' @param filter.GOF.p.value (optional, double) Filer the goodness-of-fit p-value of
#' GOF segregation distortion. Default: \code{filter.GOF.p.value = NULL}.
#' @param ind.genotyped (optional, integer) Filter the number of individual
#' progenies required to keep the marker. Default: \code{ind.genotyped = 1}.

#' @param joinmap (optional) Name of the JoinMap file to write
#' in the working directory. e.g. "joinmap.turtle.family3.loc".
#' Default: \code{joinmap = NULL}.
#' @param onemap (optional) Name of the OneMap file to write
#' in the working directory. e.g. "onemap.turtle.family3.txt".
#' Default: \code{onemap = NULL}.
#' @param filename (optional) The name of the summary file written
#' in the directory. No default.
#' Default: \code{filename = NULL}.

#' @return The function returns a list with the
#' genotypes summary \code{$geno.sum}, joinmap markers \code{$joinmap.sum}
#' and onemap markers \code{$onemap.sum} summary (use $ to access each
#' components). A JoinMap and/or OneMap file can also be exported to
#' the working direcvtory.
#' @details SIGNIFICANCE results pvalue (goodness-of-fit pvalue):
#' <= 0.0001 = ****, <= 0.001 & > 0.0001 = ***,
#' <= 0.01 & > 0.001 = **, <= 0.05 & > 0.01 = *.
#' @rdname summary_genotypes
#' @export
#' @details work in progress...
#' @examples
#' \dontrun{
#' linkage.map.crab <- summary_genotypes(
#' genotypes = "batch_10.markers.tsv",
#' markers = "batch_10.genotypes_300.txt",
#' filter.monomorphic = TRUE,
#' filter.missing.band = TRUE,
#' filter.mean.log.likelihood = -10,
#' B = NULL,
#' ind.genotyped = 300,
#' joinmap = "test.loc",
#' onemap = "test.onemap.txt",
#' filename = "genotypes.summary.tsv"
#' )
#' }
#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

summary_genotypes <- function(
  genotypes,
  markers,
  filter.monomorphic = TRUE,
  filter.missing.band = TRUE,
  filter.mean.log.likelihood = NULL,
  B = NULL,
  filter.GOF = NULL,
  filter.GOF.p.value = NULL,
  ind.genotyped = 1,
  joinmap = NULL,
  onemap = NULL,
  filename = NULL) {

  message("\n\n\nThis function is no longer working with stacks >= v.2.0\n\n\n")
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()
  cat("###############################################################################\n")
  cat("########################### stackr::summary_genotypes #########################\n")
  cat("###############################################################################\n")

  if (missing(genotypes)) stop("batch_x.genotypes.txt file created by STACKS genotypes is required")
  if (missing(markers)) stop("batch_x.markers.tsv file created by STACKS genotypes is required")
  if (filter.monomorphic == FALSE) message("Keeping monomorphic markers for linkage mapping is not recommended")


  message("Importing ", genotypes)
  genotypes.file <- readr::read_tsv(file = genotypes, col_names = c("SQL_ID", "BATCH_ID", "LOCUS", "INDIVIDUALS", "GENOTYPES"), col_types = "iiiic", skip =1, progress = interactive())

  message("Importing ", markers)
  markers.file <- readr::read_tsv(file = markers, col_names = c("SQL_ID", "BATCH_ID", "LOCUS", "TYPE", "TOTAL_GENOTYPES", "MAX", "GENOTYPE_FREQS", "SEGREGATION_DISTORTION", "MEAN_LOG_LIKELIHOOD", "GENOTYPE_MAP", "UNCORRECTED_MARKER"), col_types = c("iiicidcddcc"),skip =1, progress = interactive())

  geno.sum<- genotypes.file %>%
    dplyr::select(LOCUS, INDIVIDUALS, GENOTYPES) %>%
    dplyr::group_by(LOCUS, GENOTYPES) %>%
    dplyr::summarize(COUNT = n()) %>%
    dplyr::group_by(LOCUS) %>%
    tidyr::spread(GENOTYPES, COUNT) %>%
    dplyr::ungroup(.) %>%
    merge(markers.file, by = "LOCUS") %>%
    dplyr::mutate(
      PATTERN = stringi::stri_replace_all_fixed(TYPE, c("--/ab", "aa/ab", "ab/--", "ab/aa","ab/ab", "ab/ac", "ab/cd"),
                                                c("aaxab", "aaxab", "abxaa", "abxaa", "abxab", "abxac", "abxcd"), vectorize_all = FALSE),
      JOINMAP = stringi::stri_replace_all_fixed(TYPE, c("--/ab", "aa/ab", "ab/--", "ab/aa", "ab/ab", "ab/ac", "ab/cd"),
                                                c("nnxnp", "nnxnp", "lmxll", "lmxll", "hkxhk", "efxeg", "abxcd"), vectorize_all = FALSE),
      ONEMAP = stringi::stri_replace_all_fixed(TYPE, c("--/ab", "aa/ab", "ab/--", "ab/aa", "ab/ab", "ab/ac", "ab/cd"),
                                               c("D2.15", "D2.15", "D1.10", "D1.10", "B3.7", "A.2", "A.1"), vectorize_all = FALSE),
      BANDS_EXP = ifelse(PATTERN=="abxac" | PATTERN=="abxcd", 4, ifelse(PATTERN=="abxab", 3, 2))
    ) %>%
    merge(
      genotypes.file %>%
        dplyr::select(LOCUS, INDIVIDUALS, GENOTYPES) %>%
        dplyr::filter (GENOTYPES != "--") %>%
        dplyr::group_by(LOCUS, GENOTYPES) %>%
        dplyr::summarize(COUNT = n()) %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::tally(.) %>%
        dplyr::select(LOCUS, BANDS_OBS=n)
      , by="LOCUS"
    ) %>%
    dplyr::mutate(
      POLYMORPHIC = ifelse(BANDS_OBS > 1, "polymorphic","monomorphic"),
      MISSING_BAND = ifelse(BANDS_OBS < BANDS_EXP, "missing_band", "all_bands")
    ) %>%
    dplyr::ungroup(.)


  if (filter.monomorphic == TRUE) {
    geno.sum <- geno.sum %>% dplyr::filter(POLYMORPHIC == "polymorphic")
  }

  if (filter.missing.band == TRUE) {
    #   filter(PATTERN == "abxcd") %>% # use this for testing
    geno.sum <- dplyr::filter(geno.sum, MISSING_BAND == "all_bands")
  }

  if(!is.null(filter.mean.log.likelihood)) {
    geno.sum <- geno.sum %>% dplyr::filter(MEAN_LOG_LIKELIHOOD >= filter.mean.log.likelihood)
  }


  # Segregation Distortion summary ---------------------------------------------
  if (is.null(B)) {
    sim.p.value <- FALSE
    B <- 1
  } else {# With monte carlo replicates
    sim.p.value <- TRUE
    message("p-values computed with Monte Carlo simulations")
  }

  # abxac and/or abxcd markers
  message("\nGenerating abxac and/or abxcd markers...")
  pattern.A <- geno.sum %>%
    dplyr::filter(PATTERN == "abxac" | PATTERN == "abxcd") %>%
    dplyr::group_by(LOCUS) %>%
    dplyr::do(
      CHISQ=stats::chisq.test(
        c(.$aa,.$ab, .$ac, .$bc),
        p = c(0.25, 0.25, 0.25, 0.25),
        simulate.p.value = sim.p.value, B = B)
    ) %>%
    dplyr::mutate(
      GOF = CHISQ$statistic,
      GOF_PVALUE = CHISQ$p.value
    ) %>%
    dplyr::select(-CHISQ)

  # abxab markers
  message("\nGenerating abxab markers...")

  pattern.B <- geno.sum %>%
    dplyr::filter(PATTERN == "abxab") %>%
    dplyr::group_by(LOCUS) %>%
    dplyr::do(
      CHISQ=stats::chisq.test(
        c(.$aa,.$ab, .$bb),
        p = c(0.25, 0.5, 0.25),
        simulate.p.value = sim.p.value, B = B)
    ) %>%
    dplyr::mutate(
      GOF = CHISQ$statistic,
      GOF_PVALUE = CHISQ$p.value
    ) %>%
    dplyr::select(-CHISQ)

  # aaxab and/or abxaa markers
  message("\nGenerating aaxab and/or abxaa markers")
  pattern.D <- geno.sum %>%
    dplyr::filter(PATTERN=="aaxab" | PATTERN=="abxaa") %>%
    dplyr::group_by(LOCUS) %>%
    dplyr::do(
      CHISQ=stats::chisq.test(
        c(.$aa,.$ab),
        p = c(0.5,0.5),
        simulate.p.value = sim.p.value, B = B)
    ) %>%
    dplyr::mutate(
      GOF = CHISQ$statistic,
      GOF_PVALUE = CHISQ$p.value
    ) %>%
    dplyr::select(-CHISQ)

  chisq.data <- dplyr::bind_rows(pattern.A, pattern.B, pattern.D)

  geno.sum <- dplyr::inner_join(geno.sum, chisq.data, by = "LOCUS") %>%
    dplyr::mutate(
      GOF = round(GOF,2),
      GOF_PVALUE = round(GOF_PVALUE, 4),
      SIGNIFICANCE = ifelse(GOF_PVALUE <= 0.0001, "****",
                            ifelse(GOF_PVALUE <= 0.001 & GOF_PVALUE > 0.0001, "***",
                                   ifelse(GOF_PVALUE <= 0.01 & GOF_PVALUE > 0.001, "**",
                                          ifelse(GOF_PVALUE <= 0.05 & GOF_PVALUE > 0.01, "*", ""))))
    ) %>%
    dplyr::arrange(LOCUS)


  ##filter.GOF = FALSE, filter.GOF.p.value = FALSE, ind.genotyped
  # Filter goodness-of-fit
  if(!is.null(filter.GOF)) {
    message("Segregation distortion: filtering goodness-of-fit")
    geno.sum <- dplyr::filter(geno.sum, GOF < filter.GOF)
  }

  # Filter goodness-of-fit p-value
  if(!is.null(filter.GOF.p.value)) {
    message("Segregation distortion: filtering goodness-of-fit p-value")
    geno.sum <- dplyr::filter(geno.sum, GOF_PVALUE < filter.GOF.p.value)
  }

  # Filter number of individuals genotyped
  message("\nFiltering number of individuals genotyped")
  geno.sum <- dplyr::filter(geno.sum, TOTAL_GENOTYPES > ind.genotyped)

  # pattern summary
  joinmap.sum <- geno.sum %>%
    dplyr::group_by(JOINMAP) %>%
    dplyr::tally(.) %>%
    dplyr::rename(MARKER_NUMBER = n)

  onemap.sum <- geno.sum %>%
    dplyr::group_by(ONEMAP) %>%
    dplyr::tally(.) %>%
    dplyr::rename(MARKER_NUMBER = n)

  # joinmap output--------------------------------------------------------------
  if (!is.null(joinmap)) {
    joinmap.data <- genotypes.file %>%
      dplyr::select(LOCUS, INDIVIDUALS, GENOTYPES) %>%
      dplyr::inner_join(
        geno.sum %>%
          dplyr::select(LOCUS, JOINMAP)
        , by = "LOCUS") %>%
      dplyr::mutate(
        GENOTYPES = ifelse(JOINMAP == "efxeg", (stringi::stri_replace_all_fixed(GENOTYPES, c("aa", "ab", "ac", "bc"), c("ee", "ef", "eg", "fg"), vectorize_all = FALSE)), GENOTYPES),
        GENOTYPES = ifelse(JOINMAP == "hkxhk", (stringi::stri_replace_all_fixed(GENOTYPES, c("aa", "ab", "bb"), c("hh", "hk", "kk"), vectorize_all = FALSE)), GENOTYPES),
        GENOTYPES = ifelse(JOINMAP == "lmxll", (stringi::stri_replace_all_fixed(GENOTYPES, c("aa", "ab"), c("ll", "lm"), vectorize_all = FALSE)), GENOTYPES),
        GENOTYPES = ifelse(JOINMAP == "nnxnp", (stringi::stri_replace_all_fixed(GENOTYPES, c("aa", "ab"), c("nn", "np"), vectorize_all = FALSE)), GENOTYPES),
        JOINMAP = stringi::stri_replace_all_fixed(JOINMAP, c("efxeg", "hkxhk", "lmxll", "nnxnp"), c("<efxeg>", "<hkxhk>", "<lmxll>", "<nnxnp>"), vectorize_all = FALSE)
      ) %>%
      dplyr::group_by(LOCUS, JOINMAP) %>%
      tidyr::spread(INDIVIDUALS, GENOTYPES) %>%
      dplyr::ungroup(.)

    progeny.names <- genotypes.file %>%
      dplyr::distinct(INDIVIDUALS) %>%
      dplyr::arrange(INDIVIDUALS)

    locus.number <- paste("nloc = ", nrow(geno.sum), sep = "")
    ind.number <- paste("nind = ", nrow(progeny.names), sep = "")

    # Now write the output
    joinmap.file <- file(joinmap, "write")
    cat(paste("name = fill this", "popt = fill this with cross type", locus.number , ind.number, sep = "\n"), sep = "\n", file = joinmap.file, append = TRUE)

    utils::write.table(
      joinmap.data, file = joinmap.file, append = TRUE, col.names = FALSE,
      row.names = FALSE, sep = "\t", quote = FALSE)

    suppressWarnings(
      utils::write.table(
        progeny.names, file = joinmap.file, append = TRUE, col.names = "individual names:",
        row.names = FALSE, sep = "\t", quote = FALSE)
    )

    close(joinmap.file)
    message("JoinMap filename written: ", joinmap)
  }

  # onemap output--------------------------------------------------------------
  if (!is.null(onemap)) {
    onemap.data <- genotypes.file %>%
      dplyr::select(LOCUS, INDIVIDUALS, GENOTYPES) %>%
      dplyr::inner_join(
        geno.sum %>%
          dplyr::select(LOCUS, ONEMAP)
        , by="LOCUS") %>%
      dplyr::mutate(
        GENOTYPES = ifelse(ONEMAP == "A.2" & GENOTYPES == "aa", "a",
                           ifelse(ONEMAP == "A.2" & GENOTYPES == "ac", "a",
                                  ifelse(ONEMAP == "A.2" & GENOTYPES == "bc", "ab",
                                         ifelse(ONEMAP == "B3.7" & GENOTYPES == "aa", "a",
                                                ifelse(ONEMAP == "B3.7" & GENOTYPES == "ab", "-",
                                                       ifelse(ONEMAP == "B3.7" & GENOTYPES == "bb", "ab",
                                                              ifelse(ONEMAP == "D1.10" & GENOTYPES == "aa", "a",
                                                                     ifelse(GENOTYPES == "--", "-", GENOTYPES))))))))) %>%
      dplyr::mutate(
        ONEMAP = ifelse(ONEMAP == "A.2", "D1.10", "D1.10")
      ) %>%
      dplyr::group_by(LOCUS, ONEMAP) %>%
      tidyr::spread(INDIVIDUALS, GENOTYPES) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(LOCUS = paste("*", LOCUS, sep="")) %>%
      tidyr::unite(MARKERS, LOCUS, ONEMAP, sep=" ", remove = T) %>%
      dplyr::mutate(MARKERS = paste(.[,1], .[,2], sep = "\t"))

    progeny.names <- genotypes.file %>%
      dplyr::distinct(INDIVIDUALS) %>%
      dplyr::arrange(INDIVIDUALS)


    # Now write the output
    onemap.file <- file(onemap, "write")
    cat(paste(nrow(progeny.names), nrow(geno.sum), sep = "\t"), sep = "\n", file = onemap.file, append = TRUE)
    utils::write.table(onemap.data, file = onemap.file, append = TRUE, col.names = FALSE,
                       row.names = FALSE, sep = ",", quote = FALSE)
    close(onemap.file)
    message("OneMap filename written: ", onemap)
  }

  # Save/Write the file to the working directory--------------------------------
  if (!is.null(filename)) {
    readr::write_tsv(x = geno.sum, file = filename, append = FALSE, col_names = TRUE)
    message("genotypes summary written: ", filename)
  }

  res <- list()
  res$geno.sum <- geno.sum
  res$joinmap.sum <- joinmap.sum
  res$onemap.sum <- onemap.sum

  message("\nComputation time: ", round((proc.time() - timing)[[3]]), " sec")
  cat("################################## completed ##################################\n")
  options(width = opt.change)
  return(res)
}

