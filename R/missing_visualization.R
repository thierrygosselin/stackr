#' @name missing_visualization
#' @title Visualize missing genotypes in genomic data set.
#' @description Use this function to visualize pattern of missing data.
#'  \itemize{
#'    \item \strong{Imput file:} various file format are supported
#'    (see \code{data} argument below).
#'    \item \strong{Filters:} genotypes, markers, individuals and populations can be
#'   filtered and/or selected in several ways using blacklist,
#'   whitelist and other arguments.
#'    \item \strong{IBM-PCoA:} conduct identity-by-missingness analyses using
#'    Principal Coordinates Analysis, PCoA (also called Multidimensional Scaling, MDS).
#'    \item \strong{FH measure vs missingness}: missingness at the individual level
#'    is contrasted against FH, a new measure of IBDg
#'    (Keller et al., 2011; Kardos et al., 2015; Hedrick & Garcia-Dorado, 2016)
#'    FH is based on the excess in the observed number of homozygous
#'    genotypes within an individual relative to the mean number of homozygous
#'    genotypes expected under random mating.
#'    IBDg is a proxy measure of the realized proportion of the genome
#'    that is identical by descent.
#'    Within this function, we're using a modified version of the measure
#'    described in (Keller et al., 2011; Kardos et al., 2015).
#'    The new measure is population-wise and tailored for RADseq data
#'    (see \code{\link[stackr]{ibdg_fh}} for details).
#'    \item \strong{Figures and Tables:} figures and summary tables of
#'    missing information at the marker, individual and population level are
#'    generated.
#'    \item \strong{Blacklist:} create blacklist of individuals based on
#'    desired thresholds of missing genotypes.
#'    \item \strong{Tidy data:} if the filename argument is used, the
#'    function also output the data in a tidy format.
#' }

#' @inheritParams tidy_genomic_data

#' @param strata.select (optional, character) Use this argument to select the column
#' from the strata file to generate the PCoA-IBM plot. More than 1 column you
#' want to visualize, use a string of character
#' e.g. \code{strata.select = c("POP_ID", "LANES", "SEQUENCER", "WATERSHED")} to test
#' 4 grouping columns inside the \code{strata} file.
#' Default: \code{strata.select = "POP_ID"}

#' @param distance.method (character) The distance measure to be used.
#' This must be one of "euclidean", "maximum", "manhattan", "canberra",
#' "binary" or "minkowski". The function uses \code{\link[stats]{dist}}.
#' Default: \code{distance.method = "euclidean"}

#' @param ind.missing.geno.threshold (string) Percentage of missing genotype
#' allowed per individuals.
#' Default:\code{ind.missing.geno.threshold = c(10,20,30,40,50,60,70)}

#' @param filename (optional) Name of the tidy data set,
#' written to the working directory.

#' @return A list is created with several objects: the tidy data,
#' the principal coordinates
#' with eigenvalues of the PCoA, the identity-by-missingness plot, several
#' summary tables and plots of missing information
#' per individuals, populations and markers. Blacklisted id are also included.
#' Whitelists of markers with different level of missingness are also created.
#' A heatmap showing the missing values in black and genotypes in grey provide a
#' general overview of the missing data.

#' @examples
#' \dontrun{
#' Using a  VCF file, the simplest for of the function:
#' ibm.koala <- missing_visualization(
#' data = "batch_1.vcf",
#' strata = "population.map.strata.tsv"
#' )
#' # To see what's inside the list
#' names(ibm.koala)
#' # To view the heatmap:
#' ibm.koala$heatmap
#' # To view the IBM analysis plot:
#' ibm.koala$ibm_plot
#' }

#' @references Legendre, P. and Legendre, L. (1998) Numerical Ecology,
#' 2nd English edition. Amsterdam: Elsevier Science BV.
#' @references Keller MC, Visscher PM, Goddard ME (2011)
#' Quantification of inbreeding due to distant ancestors and its detection
#'  using dense single nucleotide polymorphism data. Genetics, 189, 237–249.
#' @references Kardos M, Luikart G, Allendorf FW (2015)
#' Measuring individual inbreeding in the age of genomics: marker-based
#' measures are better than pedigrees. Heredity, 115, 63–72.
#' @references Hedrick PW, Garcia-Dorado A. (2016)
#' Understanding Inbreeding Depression, Purging, and Genetic Rescue.
#' Trends in Ecology and Evolution. 2016;31: 940-952.
#' doi:10.1016/j.tree.2016.09.005

#' @export
#' @rdname missing_visualization
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid geom_histogram aes_string scale_fill_manual theme_bw stat_smooth geom_boxplot
#' @importFrom dplyr distinct rename arrange mutate select summarise group_by ungroup filter inner_join left_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_replace_all_regex
#' @importFrom utils count.fields
#' @importFrom readr read_tsv write_tsv
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom ape pcoa
#' @importFrom stats dist
#' @importFrom tibble data_frame


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

missing_visualization <- function(
  data,
  vcf.metadata = FALSE,
  strata = NULL,
  strata.select = "POP_ID",
  distance.method = "euclidean",
  ind.missing.geno.threshold = c(10,20,30,40,50,60,70),
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  max.marker = NULL,
  snp.ld = NULL,
  common.markers = FALSE,
  filename = NULL
) {
message("stackr: missing_visualization is now deprecated, please use grur::missing_visualization")
}
