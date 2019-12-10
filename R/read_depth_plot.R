#' @name read_depth_plot
#' @title Generate a figure with the read depth groups
#' @description This function reads the fastq file of an individual and generate
#' a figure of read coverage groups.
#'
#' @param fq.file (character, path). The path to the individual fastq file to check.
#' Default: \code{fq.file = "my-sample.fq.gz"}.

#' @param min.coverage.fig (character, path). Minimum coverage used to draw the
#' color on the figure.
#' Default: \code{min.coverage.fig = 7L}.


#' @param parallel.core (integer) Enable parallel execution with the number of threads.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @details
#'
#' 4 read coverage groups are shown:
#' \enumerate{
#' \item distinct reads with low coverage (in red): these reads are likely
#' sequencing errors or uninformative polymorphisms (shared only by a few samples).
#' \item disting reads for a target coverage (in green):
#'   \itemize{
#'   \item Usually represent around 80% of the reads in the FQ file.
#'   \item It’s a safe coverage range to start exploring your data (open for discussion).
#'   \item Lower threshold (default = 7): you can’t escape it, it’s your tolerance
#' to call heterozygote a true heterozygote. You want a minimum coverage for
#' both the reference and the alternative allele. Yes, you can use population
#' information to lower this threshold or use some fancy bayesian algorithm.
#'   \item Higher threshold: is a lot more open for discussion, here it’s the lower
#' limit of another group (the orange, see below for description). Minus 1 bp.
#' }
#' \item distinct reads with high coverage > 1 read depth (in yellow): those are
#' legitimate alleles with high coverage.
#' \item distinct and unique reads with high coverage (in orange): those
#' alleles when they are assembled in locus are usually paralogs, transposable elements, etc.
#' }

#' @rdname read_depth_plot
#' @export

#' @return The function returns the read depth groups plot.

#' @examples
#' \dontrun{
#' require(vroom)
#' check.reads.depth.groups <- read_depth_plot(fq.file = "my-sample.fq.gz")
#' }

#' @references Ilut, D., Nydam, M., Hare, M. (2014).
#' Defining Loci in Restriction-Based Reduced Representation Genomic Data from Non
#' model Species:
#' Sources of Bias and Diagnostics for Optimal Clustering BioMed Research
#' International  2014.
#' https://dx.doi.org/10.1155/2014/675158


read_depth_plot <- function(
  fq.file,
  min.coverage.fig = 7L,
  parallel.core = parallel::detectCores() - 1
) {
  if (!"vroom" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install vroom for this option:\n
                 install.packages("vroom")')
  }
  sample.clean <- stringi::stri_replace_all_fixed(
    str = fq.file,
    pattern = c(".fq.gz", ".fq", ".fasta", ".fastq", ".gzfasta", ".gzfastq", ".fastq.gz", ".FASTQ.gz", ".FASTQ.GZ"),
    replacement = c("", "", "", "", "", "", "", "", ""),
    vectorize_all = FALSE
  )
  message("Sample: ", sample.clean)


  read.stats <- vroom::vroom(
    file = fq.file,
    col_names = "READS",
    col_types = "c",
    num_threads = parallel.core,
    progress = TRUE
  ) %>%
    dplyr::mutate(
      INFO = seq.int(from = 1L, to = n()),
      SEQ = rep(1:4, n() / 4)
    ) %>%
    dplyr::filter(SEQ == 2L)

  total.sequences <- nrow(read.stats)
  message("Number of reads: ", total.sequences)

  read.stats %<>%
    dplyr::group_by(READS) %>%
    dplyr::tally(name = "DEPTH") %>%
    dplyr::group_by(DEPTH) %>%
    dplyr::tally(name = "NUMBER_DISTINCT_READS") %>%
    dplyr::ungroup(.)

  # sum(read.stats$NUMBER_DISTINCT_READS)

  depth.group.levels <- c("low coverage", "target", "high coverage", "distinct reads")

  max.coverage <- read.stats %>%
    dplyr::filter(NUMBER_DISTINCT_READS == 1L)

  MIN_COVERAGE <- min.coverage.fig
  MAX_COVERAGE <- min(max.coverage$DEPTH) - 1

  hap.read.depth <- read.stats %>%
    dplyr::mutate(
      DEPTH_GROUP = dplyr::case_when(
        NUMBER_DISTINCT_READS == 1L ~ "distinct reads",
        DEPTH < MIN_COVERAGE ~ "low coverage",
        DEPTH >= MIN_COVERAGE & DEPTH <= MAX_COVERAGE ~ "target",
        DEPTH > MAX_COVERAGE ~ "high coverage"
      ),
      DEPTH_GROUP = factor(x = DEPTH_GROUP, levels = depth.group.levels, ordered = TRUE)
    )
  distinct.sequences <- total.number.distinct.reads <- sum(hap.read.depth$NUMBER_DISTINCT_READS)

  color.tibble <- tibble::tibble(
    DEPTH_GROUP = c("low coverage", "target", "high coverage", "distinct reads"),
    LABELS = c("low coverage", paste0("target [", MIN_COVERAGE, " - ", MAX_COVERAGE, "]"), "high coverage > 1 reads", "high coverage, unique reads"),
    GROUP_COLOR = c("red", "green", "yellow", "orange")
  ) %>%
    dplyr::mutate(
      DEPTH_GROUP = factor(x = DEPTH_GROUP, levels = depth.group.levels, ordered = TRUE)
    )

  hap.read.depth.group.stats <- hap.read.depth %>%
    dplyr::mutate(DISTINCT_READS_DEPTH = DEPTH * NUMBER_DISTINCT_READS) %>%
    dplyr::group_by(DEPTH_GROUP) %>%
    dplyr::summarise(
      # NUMBER_READS_PROP = round((sum(NUMBER_DISTINCT_READS) / total.number.distinct.reads), 4),
      NUMBER_READS_PROP = round((sum(DISTINCT_READS_DEPTH) / total.sequences), 4)
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::left_join(color.tibble, by = "DEPTH_GROUP") %>%
    dplyr::mutate(
      LABELS = stringi::stri_join(LABELS, " (", as.character(format(NUMBER_READS_PROP, scientific = FALSE)), ")")
    )
  base_breaks <- function(n = 10){
    function(x) {
      grDevices::axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
  }

  read.depth.plot <- ggplot2::ggplot(data = hap.read.depth, ggplot2::aes(x = DEPTH, y = NUMBER_DISTINCT_READS)) +
    ggplot2::geom_point(ggplot2::aes(colour = hap.read.depth$DEPTH_GROUP)) +
    ggplot2::labs(x = "Depth of sequencing (log10)", y = "Number of distinct reads (log10)") +
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
  filename.plot <- stringi::stri_join(sample.clean, "_hap_read_depth.png")
  ggplot2::ggsave(
    plot = read.depth.plot,
    filename = filename.plot,
    width = 25,
    height = 15,
    dpi = 300,
    units = "cm"
  )
  return(read.depth.plot)
}# End read_depth_plot
