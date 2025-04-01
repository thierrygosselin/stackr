#' @name summary_reads
#' @title Summarise the reads for indel and GC content and produce the read depth plot
#' @description Still in dev, but work nicely.
#' Summarise the reads for indel and GC content and produce the read depth plot
#'
#' @param fq.files (character, path) The path to the individual fastq file to check,
#' or the entire fq folder.

#' @param paired.end (logical) Are the files paired-end.
#' Default: \code{paired.end = FALSE}.

#' @param output (path) Write the output in a specific directory.
#' Default: \code{output = "08_stacks_results/02_summary_reads"}.

#' @param read.depth.plot (logical) To get the interesting summaries on
#' the read content, the function is very close to similar to
#' \code{\link{read_depth_plot}}, so you can produce it for each sample with
#' minimal cost on time.
#' Default: \code{read.depth.plot = TRUE}.


#' @inheritParams read_depth_plot

#' @rdname summary_reads
#' @export

#' @return The function returns the read depth groups plot and the read stats overall
#' and by read depth groups.

#' @examples
#' \dontrun{
#' require(vroom)
#' sum <- summary_reads(fq.files = "my_fq_folder")
#' }

summary_reads <- function(
    fq.files,
    paired.end = FALSE,
    output = "08_stacks_results/02_summary_reads",
    read.depth.plot = TRUE,
    min.coverage.fig = 7L,
    parallel.core = parallel::detectCores() - 1
) {
  opt.change <- getOption("width")
  options(width = 70)
  cat("#######################################################################\n")
  cat("####################### stackr::summary_reads #########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  if (!"vroom" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install vroom for this option:\n
                 install.packages("vroom")')
  }

  if (assertthat::is.string(fq.files) && assertthat::is.dir(fq.files)) {
    fq.files <- stackr::list_sample_file(f =  fq.files, full.path = TRUE, recursive = TRUE, paired.end = paired.end)
    message("Analysing ", length(fq.files), " samples...")
  }

  if (!dir.exists(output)) dir.create(output)



  if (length(fq.files) > 1) {
    future::plan(future::multisession, workers = parallel.core)

    p <- progressr::progressor(steps = length(fq.files))
    res <- furrr::future_map(
      .x = fq.files,
      .f = reads_stats,
      read.depth.plot = read.depth.plot,
      min.coverage.fig = min.coverage.fig,
      output = output,
      parallel.core = parallel.core,
      p = p,
      verbose = FALSE
    )


  } else {
    res <- reads_stats(
      fq.files = fq.files,
      read.depth.plot = read.depth.plot,
      min.coverage.fig = min.coverage.fig,
      output = output,
      parallel.core = parallel.core,
      p = NULL,
      verbose = TRUE
    )
  }

  timing <- proc.time() - timing
  options(width = opt.change)
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}#summary_reads



# Internal function required ---------------------------------------------------
#' @title reads_stats
#' @description Main function to summarise the reads
#' @rdname reads_stats
#' @export
#' @keywords internal
reads_stats <- function(
    fq.files,
    read.depth.plot = TRUE,
    min.coverage.fig = 7L,
    output = "08_stacks_results/02_summary_reads",
    parallel.core = parallel::detectCores() - 1,
    p = NULL,
    verbose = TRUE
) {

  if (!is.null(p)) p()
  clean.names <- stackr::clean_fq_filename(basename(fq.files))
  if (verbose) message("Sample name: ", clean.names)

  read.stats <- vroom::vroom(
    file = fq.files,
    col_names = "READS",
    col_types = "c",
    delim = "\t",
    num_threads = parallel.core,
    progress = TRUE
  ) %>%
    dplyr::mutate(
      INFO = seq.int(from = 1L, to = n()),
      SEQ = rep(1:4, n() / 4)
    ) %>%
    dplyr::filter(SEQ == 2L)

  total.sequences <- nrow(read.stats)
  if (verbose) message("Number of reads: ", total.sequences)

  # stats ----------------------------------------------------------------------
  if (verbose) message("Calculating reads stats...")
  read.stats %<>% dplyr::count(READS, name = "DEPTH")

  # detect indel and / or low quality reads...----------------------------------
  # Number of N
  # GC-content/ratio (or guanine-cytosine content), proportion of nitrogenous bases in a DNA
  # DNA with low GC-content is less stable than DNA with high GC-content


  indel <- read.stats %>%
    dplyr::mutate(
      LENGTH = stringi::stri_length(str = READS),
      N = stringi::stri_count_fixed(str = READS, pattern = "N"),
      N_PROP = round(N / LENGTH, 2),
      GC = stringi::stri_count_fixed(str = READS, pattern = "C") +
        stringi::stri_count_fixed(str = READS, pattern = "G"),
      GC_PROP = round(GC / LENGTH, 2)
    )

  indel.stats.by.depth.group <- stats_stackr(data = indel, x = "N", group.by = "DEPTH")
  gc.ratio.by.depth.group <- stats_stackr(data = indel, x = "GC_PROP", group.by = "DEPTH")

  stats.overall <- dplyr::bind_rows(
    stats_stackr(data = indel, x = "N") %>% dplyr::mutate(GROUP = "INDEL", .before = 1L),
    stats_stackr(data = indel, x = "GC_PROP") %>% dplyr::mutate(GROUP = "GC", .before = 1L)
  ) %>%
    tibble::add_column(.data = ., INDIVIDUALS = clean.names, .before = 1L) %>%
    tibble::add_column(.data = ., TOTAL_READS = total.sequences, .before = 2L)

  vroom::vroom_write(x = stats.overall, file = file.path(output, paste0(clean.names, "_stats.overall.tsv")))

  read.stats %<>% dplyr::count(DEPTH, name = "NUMBER_DISTINCT_READS")
  depth.group.levels <- c("low coverage", "target", "high coverage", "distinct reads")
  max.coverage.fig <- suppressWarnings(min(read.stats$DEPTH[read.stats$NUMBER_DISTINCT_READS == 1L]) - 1)

  read.stats %<>%
    dplyr::mutate(
      DEPTH_GROUP = dplyr::case_when(
        NUMBER_DISTINCT_READS == 1L ~ "distinct reads",
        DEPTH < min.coverage.fig ~ "low coverage",
        DEPTH >= min.coverage.fig & DEPTH <= max.coverage.fig ~ "target",
        DEPTH > max.coverage.fig ~ "high coverage"
      ),
      DEPTH_GROUP = factor(x = DEPTH_GROUP, levels = depth.group.levels, ordered = TRUE)
    )
  distinct.sequences <- total.number.distinct.reads <- sum(read.stats$NUMBER_DISTINCT_READS)

  # read_depth_plot ------------------------------------------------------------
  if (read.depth.plot) {
    color.tibble <- tibble::tibble(
      DEPTH_GROUP = c("low coverage", "target", "high coverage", "distinct reads"),
      LABELS = c("low coverage", paste0("target [", min.coverage.fig, " - ", max.coverage.fig, "]"), "high coverage > 1 reads", "high coverage, unique reads"),
      GROUP_COLOR = c("red", "green", "yellow", "orange")
    ) %>%
      dplyr::mutate(
        DEPTH_GROUP = factor(x = DEPTH_GROUP, levels = depth.group.levels, ordered = TRUE)
      )

    hap.read.depth.group.stats <- read.stats %>%
      dplyr::mutate(DISTINCT_READS_DEPTH = DEPTH * NUMBER_DISTINCT_READS) %>%
      dplyr::group_by(DEPTH_GROUP) %>%
      dplyr::summarise(
        NUMBER_READS_PROP = round((sum(DISTINCT_READS_DEPTH) / total.sequences), 4),
        .groups = "drop"
      ) %>%
      dplyr::left_join(color.tibble, by = "DEPTH_GROUP") %>%
      dplyr::mutate(
        LABELS = stringi::stri_join(LABELS, " (", as.character(format(NUMBER_READS_PROP, scientific = FALSE)), ")")
      )
    base_breaks <- function(n = 10){
      function(x) {
        grDevices::axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
      }
    }

    # safely generate the fig
    # with low reads, this generates an error... hence the safe function....
    rdp <- function(read.stats) {
      DEPTH <- NULL
      NUMBER_DISTINCT_READS <- NULL
      rdp <- ggplot2::ggplot(data = read.stats, ggplot2::aes(x = DEPTH, y = NUMBER_DISTINCT_READS)) +
        ggplot2::geom_point(ggplot2::aes(colour = DEPTH_GROUP)) +
        ggplot2::labs(
          title = paste0("Read Depth Groups for sample: ", clean.names),
          subtitle = paste0("Total reads: ", total.sequences),
          x = "Depth of sequencing (log10)",
          y = "Number of distinct reads (log10)"
        ) +
        ggplot2::annotation_logticks() +
        ggplot2::scale_colour_manual(
          name = "Read coverage groups",
          labels = hap.read.depth.group.stats$LABELS,
          values = hap.read.depth.group.stats$GROUP_COLOR
        ) +
        ggplot2::scale_x_log10(breaks = c(1, 5, 10, 25, 50, 75, 100, 250, 500, 1000)) +
        ggplot2::scale_y_log10(breaks = base_breaks()) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(size = 16, face = "bold"),
          legend.title = ggplot2::element_text(size = 16, face = "bold"),
          legend.text = ggplot2::element_text(size = 16, face = "bold"),
          legend.position = c(0.7,0.8)
        )
      rdp
    } #End rdp

    safe_ggplot <- purrr::safely(.f = rdp)

    read.depth.plot <- safe_ggplot(read.stats)

    if (is.null(read.depth.plot$error)) {
      filename.plot <- file.path(output, stringi::stri_join(clean.names, "_hap_read_depth.png"))

      ggplot2::ggsave(
        plot = read.depth.plot$result,
        filename = filename.plot,
        width = 25,
        height = 15,
        dpi = 300,
        units = "cm"
      )
      # read.depth.plot <- read.depth.plot$result

      # results ------------------------------------------------------------------
      read.stats %<>%
        dplyr::left_join(indel.stats.by.depth.group, by = "DEPTH") %>%
        dplyr::left_join(gc.ratio.by.depth.group, by = "DEPTH", suffix = c("_INDEL", "_GC")) %>%
        tibble::add_column(.data = ., INDIVIDUALS = clean.names, .before = 1L) %>%
        tibble::add_column(.data = ., TOTAL_READS = total.sequences, .before = 2L)

      vroom::vroom_write(x = read.stats, file = file.path(output, paste0(clean.names, "_stats.by.depth.groups.tsv")))
      return(clean.names)
    } else {
      problem <- stringi::stri_join(
        "Problem generating read depth plot, probably not enough sequences: ",
        total.sequences, " total sequences"
      )
      if (verbose) message(problem)
      vroom::vroom_write_lines(
        x = problem,
        file = file.path(output, stringi::stri_join(clean.names, "_hap_read_depth_problem.txt")))
      read.depth.plot <- NULL
      return(invisible(NULL))
    }
  }
}# End reads_stats
