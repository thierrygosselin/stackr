#' @name missing_visualization
#' @title Visualize missing genotypes in genomic data set.
#' @description Use this function to visualize pattern of missing data.
#'  \itemize{
#'    \item \strong{Imput file:} various file format are supported 
#'    (see \code{data} argument below)
#'    \item \strong{Filters:} genotypes, markers, individuals and populations can be 
#'   filtered and/or selected in several ways using blacklist,
#'   whitelist and other arguments
#'    \item \strong{IBM-PCoA:} Conduct identity-by-missingness analyses using 
#'    Principal Coordinates Analysis, PCoA (also called Multidimensional Scaling, MDS) 
#'    \item \strong{Figures and Tables:} create figures and summary tables of 
#'    missing information at the marker, individual and population level. 
#'    \item \strong{Blacklist:} create blacklist of individuals based on 
#'    desired thresholds of missing genotypes
#'    \item \strong{Tidy data:} if the filename argument is used, the 
#'    function also output the data in a tidy format.
#' }

#' @param data 7 options: vcf, plink, stacks haplotype file, genind, genepop, 
#' and a data frame in wide or long/tidy format. 
#' The function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}} and 
#' \code{\link[stackr]{tidy_genomic_data}}.
#' \emph{See details}.

#' @param vcf.metadata (optional, logical) For the VCF file, with \code{vcf.metadata = TRUE}, 
#' the metadata contained in the \code{FORMAT} field will be kept in the tidy 
#' data file. With default: \code{vcf.metadata = FALSE}, only the genotype information is kept.



#' @param strata (optional for data frame and PLINK files, 
#' required for VCF and haplotypes files) A tab delimited file at least 2 columns 
#' with header:
#' \code{INDIVIDUALS} and \code{STRATA}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} and any other columns can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' Default: \code{strata = NULL}. This file is particularly useful to test pattern
#' of missingness stemming from lanes, library or any other grouping.

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

#' @param pop.levels (optional, string) This refers to the levels in a factor. In this 
#' case, the id of the pop.
#' Use this argument to have the pop ordered your way instead of the default 
#' alphabetical or numerical order. e.g. \code{pop.levels = c("QUE", "ONT", "ALB")} 
#' instead of the default \code{pop.levels = c("ALB", "ONT", "QUE")}. 
#' Default: \code{pop.levels = NULL}.

#' @param pop.labels (optional, string) Use this argument to rename/relabel
#' your pop or combine your pop. e.g. To combine \code{"QUE"} and \code{"ONT"} 
#' into a new pop called \code{"NEW"}:
#' (1) First, define the levels for your pop with \code{pop.levels} argument: 
#' \code{pop.levels = c("QUE", "ONT", "ALB")}. 
#' (2) then, use \code{pop.labels} argument: 
#' \code{pop.levels = c("NEW", "NEW", "ALB")}.#' 
#' To rename \code{"QUE"} to \code{"TAS"}:
#' \code{pop.labels = c("TAS", "ONT", "ALB")}.
#' Default: \code{pop.labels = NULL}. If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.

#' @param pop.select (string, optional) Selected list of populations for 
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#'and \code{ONT} population samples (out of 20 pops).
# Default: \code{pop.select = NULL} 

#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by chromosome and/or locus and/or by snp.
#' The whitelist is in the working directory (e.g. "whitelist.txt").
#' de novo CHROM column with 'un' need to be changed to 1. 
#' In the VCF, the column ID is the LOCUS identification.
#' Default \code{whitelist.markers = NULL} for no whitelist of markers.

#' @param monomorphic.out (optional) Should the monomorphic 
#' markers present in the dataset be filtered out ? 
#' Default: \code{monomorphic.out = TRUE}.

#' @param blacklist.genotype (optional) Useful to erase genotype with below 
#' average quality, e.g. genotype with more than 2 alleles in diploid likely 
#' sequencing errors or genotypes with poor genotype likelihood or coverage. 
#' The blacklist as a minimum of 2 column headers (markers and individuals). 
#' Markers can be 1 column (CHROM or LOCUS or POS), 
#' a combination of 2 (e.g. CHROM and POS or CHROM and LOCUS or LOCUS and POS) or 
#' all 3 (CHROM, LOCUS, POS) The markers columns must be designated: CHROM (character
#' or integer) and/or LOCUS (integer) and/or POS (integer). The id column designated
#' INDIVIDUALS (character) columns header. The blacklist must be in the working 
#' directory (e.g. "blacklist.genotype.txt"). For de novo VCF, CHROM column 
#' with 'un' need to be changed to 1. 
#' Default: \code{blacklist.genotype = NULL} for no blacklist of 
#' genotypes to erase.

#' @param snp.ld (optional) \strong{For VCF file only}. 
#' SNP short distance linkage disequilibrium pruning. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 3 options: \code{"random"} selection, \code{"first"} or
#' \code{"last"} SNP on the same short read/haplotype. For long distance linkage
#' disequilibrium pruning, see details below.
#' Default: \code{snp.ld = NULL}.

#' @param common.markers (optional) Logical. Default: \code{common.markers = TRUE}, 
#' will only keep markers in common (genotyped) between all the populations.

#' @param max.marker An optional integer useful to subsample marker number in 
#' large PLINK file. e.g. if the data set 
#' contains 200 000 markers and \code{max.marker = 10000} 10000 markers are
#' subsampled randomly from the 200000 markers. Use \code{whitelist.markers} to
#' keep specific markers.
#' Default: \code{max.marker = NULL}.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

#' @param filename (optional) Name of the tidy data set, 
#' written to the working directory.

#' @return A list is created with several objects: the tidy data, 
#' the principal coordinates 
#' with eigenvalues of the PCoA, the identity-by-missingness plot, several 
#' summary tables and plots of missing information
#' per individuals, populations and markers. Blacklisted id are also included. A
#' heatmap showing the missing values in black and genotypes in grey provide a
#' general overview of the missing data.

#' @details under construction

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

#' @export
#' @rdname missing_visualization
#' @import dplyr
#' @import readr
#' @import stringi
#' @importFrom data.table fread
#' @importFrom data.table melt.data.table
#' @importFrom data.table as.data.table
#' @importFrom ape pcoa
#' @importFrom stats dist
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("V1", "V2", "Missingness", "MISSING_GENOTYPE", "INDIVIDUALS_NUMBER", "PERC")
  )
}

missing_visualization <- function(data,
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
  
  
  
  cat("#######################################################################\n")
  cat("################### stackr: missing_visualization #####################\n")
  cat("#######################################################################\n")
  
  # manage missing arguments -----------------------------------------------------  
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # import data ----------------------------------------------------------------
  if (is.vector(data)){
    message("Using input file in your directory")
    
    input <- stackr::tidy_genomic_data(
      data = data, 
      vcf.metadata = vcf.metadata,
      blacklist.id = blacklist.id, 
      blacklist.genotype = blacklist.genotype, 
      whitelist.markers = whitelist.markers, 
      monomorphic.out = monomorphic.out, 
      max.marker = max.marker,
      snp.ld = snp.ld, 
      common.markers = common.markers,
      strata = strata, 
      pop.select = pop.select,
      filename = filename
    )
  } else {
    message("Using input file from your global environment")
    input <- read_long_tidy_wide(data = data)
  }
  
  if (!"MARKERS" %in% colnames(input) & "LOCUS" %in% colnames(input)) {
    input <- rename(.data = input, MARKERS = LOCUS)
  }
  
  # strata.df --------------------------------------------------------
  message("Including the strata file")
  
  strata.df <- input %>% 
    distinct(INDIVIDUALS, POP_ID) %>% 
    arrange(POP_ID, INDIVIDUALS)
  
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      number.columns.strata <- max(count.fields(strata, sep = "\t"))
      col.types <- stri_paste(rep("c", number.columns.strata), collapse = "")
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>% 
        rename(POP_ID = STRATA)
    } else {
      colnames(strata) <- stri_replace_all_fixed(str = colnames(strata), 
                                                 pattern = "STRATA", 
                                                 replacement = "POP_ID", 
                                                 vectorize_all = FALSE
      )
    }
  }
  strata.df <- strata.df %>% 
    mutate(
      INDIVIDUALS = stri_replace_all_fixed(
        str = INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      ),
      POP_ID = stri_replace_all_fixed(POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)
    )
  input <- suppressWarnings(left_join(x= input, y = strata.df, by = c("INDIVIDUALS", "POP_ID")))
  
  if(is.null(pop.levels)) { # no pop.levels
    input <- mutate(.data = input, POP_ID = factor(POP_ID))
  } else { # with pop.levels
    message("reordering your populations based on pop.levels...")
    input <- mutate(
      .data = input,
      POP_ID = factor(
        stri_replace_all_regex(
          POP_ID, stri_paste("^", pop.levels, "$", sep = ""), pop.labels, vectorize_all = FALSE),
        levels = unique(pop.labels), ordered = TRUE
      )
    )
  }
  
  # create an empty list to store results
  res <- list()
  
  
  # Identity-by-missingness (IBM) analysis -------------------------------------
  # MultiDimensional Scaling analysis (MDS) - Principal Coordinates Analysis (PCoA)
  message("Principal Coordinate Analysis (PCoA)...")
  
  input.prep <- input %>% 
    mutate(GT = ifelse(GT == "000000", "0", "1"), GT = as.numeric(GT))
  
  input.pcoa <- input.prep %>% 
    select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
    group_by(POP_ID, INDIVIDUALS) %>% 
    tidyr::spread(data = ., key = MARKERS, value = GT)
  
  # we need rownames for this
  rownames(input.pcoa) <- input.pcoa$INDIVIDUALS
  input.pcoa <- input.pcoa %>% ungroup() %>% select(-POP_ID, -INDIVIDUALS)
  
  # euclidean distances between the rows
  # distance.method <- "euclidean"
  d <- stats::dist(x = input.pcoa, method = distance.method)
  
  # alternative tested
  # d <- vegan::vegdist(x = input.pcoa, method = distance.method) # longer than stats::dist
  
  # for metric PCoA/MDS
  ibm <- ape::pcoa(D = d) # With Legendre's ape
  
  # alternative tested giving the same results:
  # ibm <- stats::cmdscale(d, eig = TRUE, k = 2) 
  
  # for non-metric PCoA/MDS
  # ibm <- MASS::isoMDS(d, k=2) # k is the number of dim
  
  # alternative: sammon
  # ibm <- MASS::sammon(d, k =2)
  
  # all gives the same results...
  
  # prep the data for figure:
  # for MASS::isoMDS and stats::cmdscale
  # pcoa.df <- data_frame(INDIVIDUALS = rownames(ibm$points), V1 = ibm$points[,1], V2 = ibm$points[,2]) %>% 
  #   inner_join(strata.df, by = "INDIVIDUALS")
  
  # with vegan and ape
  pcoa.df <- data_frame(INDIVIDUALS = rownames(ibm$vectors), V1 = ibm$vectors[,1], V2 = ibm$vectors[,2]) %>% 
    inner_join(strata.df, by = "INDIVIDUALS")
  
  if (!is.null(pop.levels)) {
    pcoa.df <- pcoa.df %>% 
      mutate(
        POP_ID = factor(
          stri_replace_all_regex(
            POP_ID, 
            stri_paste("^", pop.levels, "$", sep = ""), 
            pop.labels, 
            vectorize_all = FALSE),
          levels = unique(pop.labels), ordered = TRUE
        )
      ) %>% 
      arrange(POP_ID)
  }
  
  message("Generating Identity by missingness plot")
  # strata.select <- c("POP_ID", "WATERSHED", "BARRIER", "LANES")
  for (i in strata.select) {
    # i <- "WATERSHED"
    # IBM Figures
    
    ibm.plot <- ggplot(pcoa.df, aes(x = V1, y = V2), environment = environment())+
      geom_point(aes_string(colour = i))+
      labs(title = stri_paste("Principal Coordinates Analysis (PCoA)\n Identity by Missing (IBM) with strata = ", i))+
      labs(x = "PC1")+
      labs(y = "PC2")+
      theme(axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
            axis.title.y = element_text(size = 12, family = "Helvetica", face = "bold"),
            legend.title = element_text(size = 12, family = "Helvetica", face = "bold"), 
            legend.text = element_text(size = 12, family = "Helvetica", face = "bold"), 
            strip.text.x = element_text(size = 12, family = "Helvetica", face = "bold")
      )
    ibm_plot_name <- stri_paste("ibm.plot.strata.", i)
    res[[ibm_plot_name]] <- ibm.plot
  }
  
  # Heatmap --------------------------------------------------------------------
  heatmap.data <- input.prep %>% 
    mutate(
      GT = as.character(GT),
      Missingness = stri_replace_all_regex(GT, pattern = c("^0$", "^1$"), replacement = c("missing", "genotyped"), vectorize_all = FALSE)
    )
  
  heatmap <- ggplot(heatmap.data,(aes(y = MARKERS, x = as.character(INDIVIDUALS)))) +
    geom_tile(aes(fill = Missingness)) +
    scale_fill_manual(values = c("grey", "black"))+
    labs(y = "Markers")+
    labs(x = "Individuals")+
    theme_bw()+
    theme(
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y = element_blank(), 
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_blank(), 
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_blank()
    )
  # heatmap
  
  # Missing summary ------------------------------------------------------------
  message("Generating missing information summary tables and plot")
  
  # Individuals-----------------------------------------------------------------
  message("Missingness per individuals")
  missing.genotypes.ind <- input.prep %>% 
    select(MARKERS, INDIVIDUALS, POP_ID, GT) %>% 
    group_by(INDIVIDUALS, POP_ID) %>% 
    summarise(
      MISSING_GENOTYPE = length(GT[GT == 0]),
      MARKER_NUMBER = length(MARKERS),
      PERC = round((MISSING_GENOTYPE/MARKER_NUMBER)*100, 2)
    ) %>% 
    ungroup() %>% 
    arrange(POP_ID, INDIVIDUALS)
  
  
  # Figure Individuals
  missing.genotypes.ind.plot <- ggplot(data = missing.genotypes.ind, aes(x = INDIVIDUALS, y = PERC, colour = POP_ID)) + 
    geom_point() + 
    labs(y = "Missing genotypes (percent)")+
    labs(x = "Individuals")+
    labs(colour = "Populations") +
    # theme_bw()+
    theme(
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y = element_blank(), 
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_blank(), 
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    )
  # missing.genotypes.ind.plot
  
  # Blacklist individuals-------------------------------------------------------
  # ind.missing.geno.threshold <- c(10, 20,30,50,70)
  for (i in ind.missing.geno.threshold) {
    # i <- 30
    blacklist.id.missing.geno <- missing.genotypes.ind %>% 
      filter(PERC >= i) %>%
      ungroup() %>% 
      select(INDIVIDUALS)
    if (length(blacklist.id.missing.geno$INDIVIDUALS) > 0) {
      blacklist_name <- stri_paste("blacklist.id.missing.", i)
      res[[blacklist_name]] <- blacklist.id.missing.geno
      write_tsv(blacklist.id.missing.geno, stri_paste(blacklist_name, ".tsv"))
    }
  }
  
  # Populations-----------------------------------------------------------------
  message("Missingness per populations")
  missing.genotypes.pop <- missing.genotypes.ind %>% 
    select(INDIVIDUALS, POP_ID, PERC) %>% 
    group_by(POP_ID) %>% 
    summarise(PERC = round(mean(PERC, na.rm = TRUE), 2))
  
  # Figure Populations
  missing.genotypes.pop.plot <- ggplot(data = missing.genotypes.ind, aes(x = POP_ID, y = PERC)) + 
    geom_boxplot()+
    labs(y = "Missing genotypes (percent)")+
    labs(x = "Populations")+
    # theme_bw()+
    theme(
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y = element_blank(), 
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_text(size = 8, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    )
  missing.genotypes.pop.plot
  
  # Markers---------------------------------------------------------------------
  message("Missingness per markers")
  
  missing.genotypes.markers <- input.prep %>% 
    select(MARKERS, INDIVIDUALS, GT) %>% 
    group_by(MARKERS) %>% 
    summarise(
      MISSING_GENOTYPE = length(GT[GT == 0]),
      INDIVIDUALS_NUMBER = length(INDIVIDUALS),
      PERC = round((MISSING_GENOTYPE/INDIVIDUALS_NUMBER)*100, 2)
    ) %>% 
    ungroup() %>% 
    arrange(MARKERS)
  
  # Figure markers
  missing.genotypes.markers.plot <- ggplot(data = missing.genotypes.markers, aes(x = PERC)) + 
    geom_histogram()+
    labs(x = "Missing genotypes (percent)")+
    labs(y = "Markers (number")+
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
  # missing.genotypes.markers.plot
  
  # Results --------------------------------------------------------------------
  res$tidy.data <- input
  res$tidy.data.binary <- input.prep
  res$vectors <- ibm$vectors
  res$heatmap <- heatmap
  res$missing.genotypes.ind <- missing.genotypes.ind
  res$missing.genotypes.ind.plot <- missing.genotypes.ind.plot
  res$missing.genotypes.pop <- missing.genotypes.pop
  res$missing.genotypes.pop.plot <- missing.genotypes.pop.plot
  res$missing.genotypes.markers <- missing.genotypes.markers
  res$missing.genotypes.markers.plot <- missing.genotypes.markers.plot
  return(res)
}
