# Make genomic input file tidy

#' @name tidy_genomic_data
#' @title Transform common genomic dataset format in a tidy data frame
#' @description Transform genomic data set produced by massive parallel sequencing pipeline (e.g.GBS/RADseq, 
#' SNP chip, etc) into a tidy format. The use of blacklist and whitelist along 
#' several filtering options are available to prune the dataset. 
#' Several arguments are available to make your data population-wise and easily 
#' rename the pop id.
#' Used internally in \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data 6 options: vcf (to make vcf population ready, see details below),
#' plink, stacks haplotype file, genind, genlight, genepop, 
#' and a data frame in wide format. \emph{See details}.

#' @param vcf.metadata (optional, logical) For the VCF file, with \code{vcf.metadata = FALSE}, 
#' only the genotype information is kept.
#' With default: \code{vcf.metadata = TRUE}, the metadata contained in the 
#' \code{FORMAT} field will be kept in the tidy data file. 

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
#' The blacklist has a minimum of 2 column headers (markers and individuals). 
#' Markers can be 1 column (CHROM or LOCUS or POS), 
#' a combination of 2 (e.g. CHROM and POS or CHROM and LOCUS or LOCUS and POS) or 
#' all 3 (CHROM, LOCUS, POS). The markers columns must be designated: CHROM (character
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

#' @param maf.thresholds (string, double, optional) String with 
#' local/populations and global/overall maf thresholds, respectively.
#' e.g. \code{maf.thresholds = c(0.05, 0.1)} for a local maf threshold 
#' of 0.05 and a global threshold of 0.1. Available for VCF, PLINK and data frame 
#' files.
#' Default: \code{maf.thresholds = NULL}. 

#' @param maf.approach (character, optional). 
#' \code{maf.approach = "haplotype"} : looks at the minimum MAF found on the 
#' read/haplotype. Using this option will discard all the markers/snp on 
#' that read based on the thresholds chosen. This method is only available 
#' for VCF and haplotype files, or tidy data frame from those file types.
#' \code{maf.approach = "SNP"} : treats all the SNP on the same 
#' haplotype/read as independent. Doesn't work with haplotype file, 
#' but does work for all other file type.
#' Default is \code{maf.approach = "SNP"}.

#' @param maf.operator (character, optional) \code{maf.operator = "AND"} or 
#' default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?

#' @param maf.pop.num.threshold (integer, optional) When maf thresholds are used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. Default: \code{maf.pop.num.threshold = 1}

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

#' @param strata (optional for data frame and PLINK files, 
#' required for VCF and haplotypes files) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' If you have already run 
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data, 
#' the strata file is similar to a stacks `population map file`, make sure you 
#' have the required column names  (\code{INDIVIDUALS} and \code{STRATA}).
#' Default: \code{strata = NULL}.

#' @param pop.select (string, optional) Selected list of populations for 
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#' and \code{ONT} population samples (out of 20 pops).
#' Default: \code{pop.select = NULL} 


#' @param filename (optional) The file name for the tidy data frame
#' written to the working directory.
#' Default: \code{filename = NULL}, the tidy data is 
#' in the global environment only (i.e. not written in the working directory).

#' @details 
#' \strong{Long distance SNP linkage disequilibrium pruning}
#' #' If you have markers position on a genome or a linkage map,
#' you can go further in removing linked markers by using \pkg{SNPRelate} or
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK}, \emph{linkage 
#' disequilibrium based SNP pruning} option.
#' 
#' \strong{Input files:}
#' \enumerate{
#' \item VCF file (e.g. \code{data = "batch_1.vcf"}). 
#' To make the VCF population ready, you need the \code{strata} argument.
#' 
#' \item haplotype file created in STACKS (e.g. \code{data = "batch_1.haplotypes.tsv"}).
#' To make the haplotype file population ready, you need the \code{strata} argument.
#' 
#' \item Data frame
#' To discriminate the long from the wide format, 
#' the function \pkg{stackr} \code{\link[stackr]{read_long_tidy_wide}} searches 
#' for \code{MARKERS or LOCUS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns: 
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals), 
#' the remaining columns are the markers in separate columns storing genotypes.
#' 
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info. 
#' (e.g. from a VCF see \pkg{stackr} \code{\link[stackr]{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID}, 
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#' 
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 of 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 of 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#' 
#' \emph{How to get a tidy data frame ?}
#' \pkg{stackr} \code{\link[stackr]{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame.
#' 
#' \item PLINK file in 
#' \code{tped/tfam} format (e.g. \code{data =  "data.assignment.tped"}). 
#' The first 2 columns of the \code{tfam} file will be used for the 
#' \code{strata} argument below, unless a new one is provided. 
#' Columns 1, 3 and 4 of the \code{tped} are discarded. The remaining columns 
#' correspond to the genotype in the format \code{01/04} 
#' where \code{A = 01, C = 02, G = 03 and T = 04}. For \code{A/T} format, use 
#' PLINK or bash to convert.
#' Use \href{http://vcftools.sourceforge.net/}{VCFTOOLS} with \code{--plink-tped} 
#' to convert very large VCF file. For \code{.ped} file conversion to 
#' \code{.tped} use \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK} 
#' with \code{--recode transpose},
#' 
#' \item \code{\link[adegenet]{genind}} object from \code{\link[adegenet]{adegenet}}.
#' 
#' \item genepop data file (e.g. \code{data = "kiwi_data.gen"}). Here, the function can only use
#' alleles encoded with 3 digits.
#' }

#' @return The output in your global environment is a tidy data frame. 
#' If \code{filename} is provided, the tidy data frame is also 
#' written to the working directory.

#' @export
#' @rdname tidy_genomic_data
#' @import parallel
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join summarise_each_ funs
#' @importFrom adegenet genind2df
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_replace_all_regex stri_sub stri_pad_left stri_count_fixed stri_replace_na 
#' @importFrom stats var median quantile
#' @importFrom purrr map flatten keep discard
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom plyr colwise
#' @importFrom tidyr spread gather unite separate
#' @importFrom utils count.fields
#' @importFrom readr write_tsv read_tsv


#' @examples
#' \dontrun{
#' tidy.vcf <- tidy_genomic_data(
#' data = "batch_1.vcf",
#' whitelist.markers = "whitelist.vcf.txt",
#' snp.ld = NULL,
#' common.markers = TRUE,
#' blacklist.id = "blacklist.id.treefrog.tsv",
#' strata = "strata.treefrog.tsv",
#' pop.levels = c("PAN", "COS")
#' )
#' }


#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, et al. 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. 
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795
#' @references Jombart T, Devillard S, Balloux F. 
#' Discriminant analysis of principal components: a new method for the analysis 
#' of genetically structured populations. 
#' BMC Genet. 2010:11: 94. doi:10.1186/1471-2156-11-94
#' @references Jombart T, Ahmed I. adegenet 1.3-1: new tools for the analysis 
#' of genome-wide SNP data. 
#' Bioinformatics. 2011:27: 3070–3071. doi:10.1093/bioinformatics/btr521
#' @references Raymond M. & Rousset F, (1995). 
#' GENEPOP (version 1.2): population genetics software for exact tests 
#' and ecumenicism. 
#' J. Heredity, 86:248-249

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_genomic_data <- function(
  data, 
  vcf.metadata = TRUE,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  blacklist.genotype = NULL,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  max.marker = NULL,
  blacklist.id = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  pop.select = NULL,
  filename = NULL
) {
  
  
  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) stop("Input file missing")
  
  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  
  if (!is.null(pop.select)) {
    pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  # File type detection----------------------------------------------------------
  data.type <- detect_genomic_format(data)
  
  if (data.type == "haplo.file") {
    message("With stacks haplotype file the maf.approach is automatically set to: haplotype")
    maf.approach <- "SNP"
    # confusing, but because the haplotpe file doesn't have snp info, only locus info
    # it's treated as markers/snp info and filtered the same way as the approach by SNP.
    # but it's really by haplotype
  }
  
  if (maf.approach == "haplotype") {
    if (data.type != "vcf.file" | data.type != "haplo.file") {
      stop("The haplotype approach during MAF filtering is for VCF and
           stacks haplotypes file, only. Use the snp approach for the other file types")
    }
  }
  
  
  # Strata argument required for VCF and haplotypes files-----------------------
  if (data.type == "haplo.file" | data.type == "vcf.file") {
    if (is.null(strata)) stop("strata argument is required")
  }
  
  # Import whitelist of markers-------------------------------------------------
  if (is.null(whitelist.markers)) { # no Whitelist
    # message("Whitelist of markers: no")
  } else {# with Whitelist of markers
    # message("Whitelist of markers: yes")
    suppressMessages(whitelist.markers <- readr::read_tsv(whitelist.markers, col_names = TRUE))
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
    if (data.type == "haplo.file") {
      whitelist.markers <- dplyr::select(.data = whitelist.markers, LOCUS)
      columns.names.whitelist <- colnames(whitelist.markers)
    }
  }
  
  # Import blacklist id --------------------------------------------------------
  if (is.null(blacklist.id)) { # No blacklist of ID
    # message("Blacklisted individuals: no")
  } else {# With blacklist of ID
    # message("Blacklisted individuals: yes")
    suppressMessages(blacklist.id <- readr::read_tsv(blacklist.id, col_names = TRUE))
  }
  
  # population levels and strata------------------------------------------------
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      # message("strata file: yes")
      number.columns.strata <- max(utils::count.fields(strata, sep = "\t"))
      col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
      suppressMessages(strata.df <- readr::read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>% 
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
  }
  
  # Check with strata and pop.levels/pop.labels
  if (!is.null(strata) & !is.null(pop.levels)) {
    if (length(levels(factor(strata.df$POP_ID))) != length(pop.levels)) {
      stop("The number of groups in your strata file must match the number of groups in pop.levels")
    }
  }
  
  
  # Import VCF-------------------------------------------------------------------
  if (data.type == "vcf.file") { # VCF
    message("Importing the VCF...")
    
    input <- data.table::fread(
      input = data,
      sep = "\t",
      stringsAsFactors = FALSE, 
      header = TRUE,
      # Automatically filter with blacklist of id
      # drop = c("QUAL", "FILTER", "INFO", blacklist.id$INDIVIDUALS),
      drop = c("QUAL", "FILTER", "INFO"),
      skip = "CHROM",
      showProgress = TRUE,
      verbose = FALSE
    ) %>% 
      tibble::as_data_frame() %>% 
      dplyr::rename(LOCUS = ID, CHROM = `#CHROM`) %>%
      dplyr::mutate(
        CHROM = stringi::stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1"),
        POS = as.character(POS),
        LOCUS = as.character(LOCUS)
      )
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Detect the format fields
    format.headers <- unlist(stringi::stri_split_fixed(str = input$FORMAT[1], pattern = ":"))
    input <- input %>% dplyr::select(-FORMAT) # no longer needed
    
    # Tidying the VCF to make it easy to work on the data for conversion
    message("Making the VCF population wise")
    input <- data.table::melt.data.table(
      data = data.table::as.data.table(input), 
      id.vars = c("CHROM", "LOCUS", "POS", "REF", "ALT"), 
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = "FORMAT_ID"
    ) %>% 
      tibble::as_data_frame() %>% 
      dplyr::mutate(
        INDIVIDUALS = stringi::stri_replace_all_fixed(
          str = INDIVIDUALS, 
          pattern = c("_", ":"), 
          replacement = c("-", "-"), 
          vectorize_all = FALSE)
      )
    
    # filter blacklisted individuals
    if (!is.null(blacklist.id)) {
      blacklist.id$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = blacklist.id$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"), 
        vectorize_all = FALSE
      )
      
      input <- dplyr::filter(.data = input, !INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
    }
    
    # population levels and strata
    strata.df$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = strata.df$INDIVIDUALS, 
      pattern = c("_", ":"), 
      replacement = c("-", "-"), 
      vectorize_all = FALSE
    )
    
    input <- dplyr::left_join(x = input, y = strata.df, by = "INDIVIDUALS")
    
    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
      input$POP_ID <- droplevels(input$POP_ID)
    }
    
    # Tidy VCF
    message("Tidying the vcf...")
    input <- tidyr::separate(
      data = input, FORMAT_ID, format.headers, sep = ":", extra = "drop"
    )
    
    if (vcf.metadata) {
      # GL cleaning
      if (length(unlist(stringi::stri_extract_all_fixed(str = format.headers, pattern = "GL", omit_no_match = TRUE))) > 0) {
        message("Fixing GL column...")
        input <- input %>% 
          dplyr::mutate( 
            GL = suppressWarnings(as.numeric(stringi::stri_replace_all_fixed(GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all = FALSE)))
          )
      } # end cleaning GL column
      
      # Cleaning DP and changing name to READ_DEPTH
      if (length(unlist(stringi::stri_extract_all_fixed(str = format.headers, pattern = "DP", omit_no_match = TRUE))) > 0) {
        message("Fixing DP (READ_DEPTH) column...")
        input <- input %>% 
          dplyr::rename(READ_DEPTH = DP) %>% 
          dplyr::mutate(
            READ_DEPTH = suppressWarnings(as.numeric(stringi::stri_replace_all_regex(READ_DEPTH, "^0$", "NA", vectorize_all = FALSE))),
            READ_DEPTH = ifelse(GT == "./.", NA, READ_DEPTH)
          )
      } # end cleaning READ_DEPTH (DP) column
      
      # Cleaning AD (ALLELES_DEPTH)
      if (length(unlist(stringi::stri_extract_all_fixed(str = format.headers, pattern = "AD", omit_no_match = TRUE))) > 0) {
        message("Splitting AD columns into allele coverage info...")
        input <- input %>% 
          tidyr::separate(AD, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"), sep = ",", extra = "drop") %>% 
          dplyr::mutate(
            ALLELE_REF_DEPTH = suppressWarnings(as.numeric(stringi::stri_replace_all_regex(ALLELE_REF_DEPTH, "^0$", "NA", vectorize_all = TRUE))),
            ALLELE_ALT_DEPTH = suppressWarnings(as.numeric(stringi::stri_replace_all_regex(ALLELE_ALT_DEPTH, "^0$", "NA", vectorize_all = TRUE)))
            # dplyr::mutate coverage ratio for allelic imbalance
            # ALLELE_COVERAGE_RATIO = suppressWarnings(
            # as.numeric(ifelse(GT == "./." | GT == "0/0" | GT == "1/1", "NA",
            # ((ALLELE_ALT_DEPTH - ALLELE_REF_DEPTH)/(ALLELE_ALT_DEPTH + ALLELE_REF_DEPTH)))))
          )
      } # end cleaning AD column
    } else {
      input <- input %>% 
        dplyr::select(CHROM, LOCUS, POS, REF, ALT, POP_ID, INDIVIDUALS, GT)
    }# end metadata section
    
    # recoding genotype and creating a new column combining CHROM, LOCUS and POS 
    input <- input %>%
      tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE) %>%
      dplyr::mutate(
        REF = stringi::stri_replace_all_fixed(
          str = REF, 
          pattern = c("A", "C", "G", "T"), 
          replacement = c("001", "002", "003", "004"), 
          vectorize_all = FALSE
        ), # replace nucleotide with numbers
        ALT = stringi::stri_replace_all_fixed(
          str = ALT, pattern = c("A", "C", "G", "T"), 
          replacement = c("001", "002", "003", "004"), 
          vectorize_all = FALSE
        ),# replace nucleotide with numbers
        GT = stringi::stri_replace_all_fixed(str = GT, pattern = "|", replacement = "/", vectorized_all = FALSE),
        GT_VCF = GT,
        GT = ifelse(GT == "0/0", stringi::stri_join(REF, REF, sep = ""),
                    ifelse(GT == "1/1",  stringi::stri_join(ALT, ALT, sep = ""),
                           ifelse(GT == "0/1", stringi::stri_join(REF, ALT, sep = ""),
                                  ifelse(GT == "1/0", stringi::stri_join(ALT, REF, sep = ""), "000000")
                           )
                    )
        )
      ) %>% 
      tibble::as_data_frame()
    
    # Experimental (for genlight object)
    input$GT_BIN <- stringi::stri_replace_all_fixed(str = input$GT_VCF, pattern = c("0/0", "1/1", "0/1", "1/0", "./."), replacement = c("0", "2", "1", "1", NA), vectorize_all = FALSE)
    
    # re-computing the REF/ALT allele
    if (!is.null(pop.select) || !is.null(blacklist.id)) {
      message("Adjusting REF/ALT alleles to account for filters...")
      input <- ref_alt_alleles(data = input, monomorphic.out = monomorphic.out)
    } # end re-computing the REF/ALT allele
    
    # Re ordering columns
    if (vcf.metadata) {
      # Re order columns
      common.colnames <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS", "GT_VCF", "GT_BIN", "REF", "ALT")
      vcf.headers <- colnames(input)
      metadata.colnames <- purrr::discard(.x = colnames(input), .p = vcf.headers %in% common.colnames)
      input <- input[c(common.colnames, metadata.colnames)]
      
    } else {
      input <- input %>% dplyr::select(MARKERS, CHROM, LOCUS, POS, POP_ID, INDIVIDUALS, GT_VCF, GT_BIN, REF, ALT, GT)
    }
    # for other part of the script
    biallelic <- TRUE
  } # End import VCF
  
  # Import PLINK ---------------------------------------------------------------
  if (data.type == "plink.file") { # PLINK
    message("Importing the PLINK files...")
    tfam <- data.table::fread(
      input = stringi::stri_replace_all_fixed(
        str = data, 
        pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE
      ),
      sep = " ", 
      header = FALSE, 
      stringsAsFactors = FALSE,
      verbose = FALSE,
      select = c(1,2),
      colClasses = list(character = c(1,2)),
      col.names = c("POP_ID", "INDIVIDUALS"),
      showProgress = TRUE, 
      data.table = FALSE) %>% 
      tibble::as_data_frame() %>%
      dplyr::mutate(
        # remove unwanted sep in individual name and replace with "-"
        INDIVIDUALS = stringi::stri_replace_all_fixed(
          str = INDIVIDUALS, 
          pattern = c("_", ":"), 
          replacement = c("-", "-"), 
          vectorize_all = FALSE),
        # remove potential whitespace in tfam pop id column
        POP_ID = stringi::stri_replace_all_fixed(
          POP_ID,
          pattern = " ", 
          replacement = "_", 
          vectorize_all = FALSE)
      )
    
    # if no strata tfam = strata.df
    if (is.null(strata)) {
      strata.df <- tfam
      
      # Check with strata and pop.levels/pop.labels
      if (!is.null(pop.levels)) {
        if (length(levels(factor(strata.df$POP_ID))) != length(pop.levels)) {
          stop("The number of groups in your tfam file must match the number of groups in pop.levels")
        }
      }
    } else {
      # remove unwanted sep in individual name and replace with "-"
      strata.df$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = strata.df$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
    }
    
    tped.header.prep <- tfam %>% 
      dplyr::select(INDIVIDUALS) %>%
      dplyr::mutate(
        NUMBER = seq(1, n()),
        ALLELE1 = rep("A1", n()), ALLELE2 = rep("A2", n())
      ) %>%
      tidyr::gather(ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, NUMBER)) %>%
      dplyr::arrange(NUMBER) %>% 
      dplyr::select(-ALLELES_GROUP) %>% 
      tidyr::unite(INDIVIDUALS_ALLELES, c(INDIVIDUALS, ALLELES), sep = "_", remove = FALSE) %>% 
      dplyr::arrange(NUMBER) %>% 
      dplyr::mutate(NUMBER = seq(from = (1 + 4), to = n() + 4)) %>% 
      dplyr::select(-ALLELES)
    
    tped.header.names <- c("LOCUS", tped.header.prep$INDIVIDUALS_ALLELES)
    tped.header.integer <- c(2, tped.header.prep$NUMBER)
    
    if (!is.null(blacklist.id)) { # using the blacklist of individuals
      # remove unwanted sep in individual name and replace with "-"
      blacklist.id$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = blacklist.id$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"), 
        vectorize_all = FALSE
      )
      
      whitelist.id <- tped.header.prep %>% 
        dplyr::anti_join(blacklist.id, by = "INDIVIDUALS") %>% 
        dplyr::arrange(NUMBER)
      tped.header.names <- c("LOCUS", whitelist.id$INDIVIDUALS_ALLELES)
      tped.header.integer <- c(2, whitelist.id$NUMBER)
      
      strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
    }
    
    # import PLINK
    input <- data.table::fread( 
      input = data, 
      sep = " ", 
      header = FALSE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      select = tped.header.integer,
      col.names = tped.header.names,
      showProgress = TRUE,
      data.table = FALSE) %>% 
      tibble::as_data_frame() %>% 
      dplyr::mutate(LOCUS = as.character(LOCUS))
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(
        dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist)
      )
    }
    
    # To reduce the size of the dataset we subsample the markers with max.marker
    if (!is.null(max.marker)) {
      message("Using the max.marker to reduce the size of the dataset")
      input <- dplyr::sample_n(tbl = input, size = max(as.numeric(max.marker)), replace = FALSE)
      
      max.marker.subsample.select <- input %>% 
        dplyr::distinct(LOCUS, .keep_all = TRUE) %>% 
        dplyr::arrange(LOCUS)
      
      readr::write_tsv(# save results
        x = max.marker.subsample.select, 
        path = "max.marker.subsample.select.tsv")
    }
    
    # Make tidy
    message("Tidying the PLINK file ...")
    # Filling GT and new separating INDIVIDUALS from ALLELES
    # combining alleles
    input <- data.table::melt.data.table(
      data = data.table::as.data.table(input), 
      id.vars = "LOCUS", 
      variable.name = "INDIVIDUALS_ALLELES", 
      value.name = as.character("GT"), 
      variable.factor = FALSE, 
      value.factor = FALSE
    ) %>% 
      tibble::as_data_frame() %>% 
      dplyr::mutate(GT = stringi::stri_pad_left(str = GT, width = 3, pad = "0")) %>% 
      tidyr::separate(data = ., col = INDIVIDUALS_ALLELES, into = c("INDIVIDUALS", "ALLELES"), sep = "_")
    
    input <- data.table::dcast.data.table(
      data = data.table::as.data.table(input), 
      formula = LOCUS + INDIVIDUALS ~ ALLELES, 
      value.var = "GT"
    ) %>% 
      tibble::as_data_frame() %>% 
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>% 
      dplyr::select(LOCUS, INDIVIDUALS, GT)
    
    # population levels and strata
    message("Integrating the tfam/strata file...")
    
    input <- dplyr::left_join(x = input, y = strata.df, by = "INDIVIDUALS")
    
    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
      input$POP_ID <- droplevels(input$POP_ID)
    }
    
    # removing untyped markers across all-pop
    remove.missing.gt <- input %>%
      dplyr::select(LOCUS, GT) %>%
      dplyr::filter(GT != "000000")
    
    untyped.markers <- dplyr::n_distinct(input$LOCUS) - dplyr::n_distinct(remove.missing.gt$LOCUS)
    if (untyped.markers > 0) {
      message(paste0("Number of marker with 100 % missing genotypes: ", untyped.markers))
      input <- suppressWarnings(
        dplyr::semi_join(input, 
                         remove.missing.gt %>% 
                           dplyr::distinct(LOCUS, .keep_all = TRUE), 
                         by = "LOCUS")
      )
    }
    
    # Unused objects
    tped.header.prep <- tped.header.integer <- tped.header.names <- remove.missing.gt <- NULL
    
    # detect if plink file was biallelic
    biallelic <- detect_biallelic_markers(input)
    
    # give vcf style genotypes
    if (biallelic) {
      input <- ref_alt_alleles(data = input, monomorphic.out = monomorphic.out)
    }
    
  } # End import PLINK
  
  # Import genepop--------------------------------------------------------------
  if (data.type == "genepop.file") {
    message("Tidying the genepop file ...")
    data <- stackr::tidy_genepop(data = data, tidy = TRUE)
    data.type <- "df.file"
  }
  
  # Import DF-------------------------------------------------------------------
  if (data.type == "df.file") { # DATA FRAME OF GENOTYPES
    input <- read_long_tidy_wide(data = data, import.metadata = vcf.metadata)
    
    # For long tidy format, switch LOCUS to MARKERS column name, if found MARKERS not found
    if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
      input <- dplyr::rename(.data = input, MARKERS = LOCUS)
    }
    
    # Change individuals names containing special character
    input$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = input$INDIVIDUALS, 
      pattern = c("_", ":"), 
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      blacklist.id$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = blacklist.id$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    # population levels and strata
    if (!is.null(strata)) {
      strata.df$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = strata.df$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
      
      input <- input %>%
        dplyr::select(-POP_ID) %>% 
        dplyr::left_join(strata.df, by = "INDIVIDUALS")
    }
    
    # Change potential problematic POP_ID space
    input$POP_ID = stringi::stri_replace_all_fixed(
      input$POP_ID, 
      pattern = " ", 
      replacement = "_", 
      vectorize_all = FALSE
    )
    
    # Check with strata and pop.levels/pop.labels
    if (!is.null(pop.levels)) {
      if (length(levels(factor(input$POP_ID))) != length(pop.levels)) {
        stop("The number of groups in your POP_ID column file must match the number of groups in pop.levels")
      }
    }
    
    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    }
    
    # detect if plink file was biallelic
    biallelic <- detect_biallelic_markers(input)
    
    # give vcf style genotypes
    if (biallelic) {
      input <- ref_alt_alleles(data = input, monomorphic.out = monomorphic.out)
    }
  } # End import data frame of genotypes
  
  # Import haplo---------------------------------------------------------------
  if (data.type == "haplo.file") { # Haplotype file
    message("Importing the stacks haplotype file")
    number.columns <- max(utils::count.fields(data, sep = "\t"))
    
    input <- data.table::fread(
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
      tibble::as_data_frame() %>% 
      dplyr::select(-Cnt)
    
    if (tibble::has_name(input, "# Catalog ID") || tibble::has_name(input, "Catalog ID")) {
      colnames(input) <- stringi::stri_replace_all_fixed(
        str = colnames(input), 
        pattern = c("# Catalog ID", "Catalog ID"), replacement = c("LOCUS", "LOCUS"), vectorize_all = FALSE
      )
    }
    
    if (tibble::has_name(input, "Seg Dist")) {
      input <- dplyr::select(.data = input, -`Seg Dist`)
    }
    
    input <- data.table::melt.data.table(
      data = data.table::as.data.table(input), 
      id.vars = "LOCUS", 
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = "GT"
    ) %>% 
      tibble::as_data_frame()
    
    number.columns <- NULL
    
    # remove consensus markers
    message("Scanning for consensus markers")
    consensus.markers <- input %>%
      dplyr::filter(GT == "consensus") %>% 
      dplyr::distinct(LOCUS, .keep_all = TRUE)
    
    if (length(consensus.markers$LOCUS) > 0) {
      input <- suppressWarnings(dplyr::anti_join(input, consensus.markers, by = "LOCUS"))
    }
    message(stringi::stri_join("Consensus markers removed: ", dplyr::n_distinct(consensus.markers$LOCUS)))
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      
      blacklist.id$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = blacklist.id$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
      
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    # population levels and strata
    strata.df$INDIVIDUALS = stringi::stri_replace_all_fixed(
      str = strata.df$INDIVIDUALS, 
      pattern = c("_", ":"), 
      replacement = c("-", "-"), 
      vectorize_all = FALSE
    )
    
    input <- dplyr::left_join(x = input, y = strata.df, by = "INDIVIDUALS")
    
    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    }
    
    # removing errors and potential paralogs (GT with > 2 alleles)
    message("Scanning for genotypes with more than 2 alleles")
    input <- input %>%
      dplyr::mutate(POLYMORPHISM = stringi::stri_count_fixed(GT, "/"))
    
    blacklist.paralogs <- input %>%
      dplyr::filter(POLYMORPHISM > 1) %>% 
      dplyr::select(LOCUS, INDIVIDUALS)
    
    message(stringi::stri_join("Number of genotypes with more than 2 alleles: ", length(blacklist.paralogs$LOCUS), sep = ""))
    
    if (length(blacklist.paralogs$LOCUS) > 0) {
      input <- input %>% 
        dplyr::mutate(GT = ifelse(POLYMORPHISM > 1, NA, GT)) %>% 
        dplyr::select(-POLYMORPHISM)
      
      readr::write_tsv(blacklist.paralogs, "blacklist.genotypes.paralogs.tsv")
    }
    
    tidy.haplo <- input %>% dplyr::select(LOCUS, INDIVIDUALS, GT_HAPLO = GT)
    
    # recode genotypes
    message("Recoding haplotypes in 6 digits format")
    # input.bk <- input
    input <- suppressWarnings(
      input %>% 
        tidyr::separate(
          col = GT, into = c("A1", "A2"), sep = "/", extra = "drop", remove = TRUE
        ) %>%
        dplyr::mutate(A2 = ifelse(is.na(A2), A1, A2)) %>% 
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(LOCUS, POP_ID, INDIVIDUALS)) %>% 
        tidyr::spread(data = ., key = LOCUS, value = GT)
    )
    
    input.id.col <- dplyr::select(.data = input, POP_ID, INDIVIDUALS, ALLELES)
    
    input.variable <- input %>% dplyr::select(-c(POP_ID, INDIVIDUALS, ALLELES)) %>% 
      plyr::colwise(factor, exclude = NA)(.) %>% 
      plyr::colwise(as.numeric)(.)
    
    input <- dplyr::bind_cols(input.id.col, input.variable)
    
    input <- data.table::melt.data.table(
      data = data.table::as.data.table(input), 
      id.vars = c("INDIVIDUALS", "POP_ID", "ALLELES"), 
      variable.name = "LOCUS",
      variable.factor = FALSE,
      value.name = "GT"
    ) %>% 
      tibble::as_data_frame() %>% 
      dplyr::mutate(
        GT = as.character(GT),
        GT = stringi::stri_pad_left(str = GT, width = 3, pad = "0"),
        GT = stringi::stri_replace_na(str = GT, replacement = "000")
      ) %>% 
      tidyr::spread(data = ., key = ALLELES, value = GT) %>% 
      tidyr::unite(data = ., GT, A1, A2, sep = "") 
    
    input <- dplyr::inner_join(input, tidy.haplo, by = c("LOCUS", "INDIVIDUALS")) %>%
      dplyr::select(LOCUS, POP_ID, INDIVIDUALS, GT, GT_HAPLO) %>%
      dplyr::arrange(LOCUS, POP_ID, INDIVIDUALS, GT, GT_HAPLO)
    
    # change the filename and strata.df here
  } # End import haplotypes file
  
  # Import GENIND--------------------------------------------------------------
  # load("/Users/thierry/Documents/skipjack/fst.test.RData")
  if (data.type == "genind.file") { # DATA FRAME OF GENOTYPES
    # data = skipjack.genind
    input <- adegenet::genind2df(data) %>% 
      tibble::rownames_to_column("INDIVIDUALS") %>% 
      dplyr::rename(POP_ID = pop)
    
    message("Tidying the genind object ...")
    # scan for the number of character coding the allele
    allele.sep <- input %>% dplyr::select(-INDIVIDUALS, -POP_ID)
    allele.sep <- unique(nchar(allele.sep[!is.na(allele.sep)]))
    
    if (length(allele.sep) > 1) {
      stop("The number of character/integer string coding the allele is not identical accross markers")
    }
    
    input <- input %>% 
      tidyr::gather(key = LOCUS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% 
      tidyr::separate(
        data = ., col = GT, into = c("A1", "A2"), 
        sep = allele.sep/2, remove = TRUE, extra = "drop"
      ) %>% 
      dplyr::mutate(
        A1 = stringi::stri_pad_left(str = A1, pad = "0", width = 3),
        A2 = stringi::stri_pad_left(str = A2, pad = "0", width = 3)
      ) %>% 
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>% 
      dplyr::mutate(GT = replace(GT, which(GT == "NANA"), "000000"))
    
    # remove unwanted sep in id and pop.id names
    input <- input %>% 
      dplyr::mutate(
        INDIVIDUALS = stringi::stri_replace_all_fixed(
          str = INDIVIDUALS, 
          pattern = c("_", ":"), 
          replacement = c("-", "-"), 
          vectorize_all = FALSE),
        POP_ID = stringi::stri_replace_all_fixed(
          POP_ID,
          pattern = " ", 
          replacement = "_", 
          vectorize_all = FALSE)
      )
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      
      blacklist.id$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = blacklist.id$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
      
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    # population levels and strata
    if (!is.null(strata)) {
      input <- input %>%
        dplyr::select(-POP_ID) %>% 
        dplyr::mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        dplyr::left_join(strata.df, by = "INDIVIDUALS")
    }
    
    # Change potential problematic POP_ID space
    input$POP_ID = stringi::stri_replace_all_fixed(
      input$POP_ID, 
      pattern = " ", 
      replacement = "_", 
      vectorize_all = FALSE
    )
    
    # Check with strata and pop.levels/pop.labels
    if (!is.null(pop.levels)) {
      if (length(levels(factor(input$POP_ID))) != length(pop.levels)) {
        stop("The number of groups in your POP_ID column must match the number of groups in pop.levels")
      }
    }
    
    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    }
    
    # detect if plink file was biallelic
    biallelic <- detect_biallelic_markers(input)
    
    # give vcf style genotypes
    if (biallelic) {
      input <- ref_alt_alleles(data = input, monomorphic.out = monomorphic.out)
    }
    
    # Now the genind and genepop are like ordinary data frames
    data.type <- "df.file" # for subsequent steps
    
  } # End tidy genind
  
  # Arrange the id and create a strata after pop select ------------------------
  input$INDIVIDUALS <- stringi::stri_replace_all_fixed(
    str = input$INDIVIDUALS, 
    pattern = c("_", ":"), 
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )
  
  strata.df <- input %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(POP_ID, INDIVIDUALS)
  
  # Blacklist genotypes --------------------------------------------------------
  if (is.null(blacklist.genotype)) { # no Whitelist
    message("Erasing genotype: no")
  } else {
    message("Erasing genotype: yes")
    suppressMessages(
      blacklist.genotype <- readr::read_tsv(blacklist.genotype, col_names = TRUE) %>% 
        dplyr::mutate(
          INDIVIDUALS = stringi::stri_replace_all_fixed(
            str = INDIVIDUALS, 
            pattern = c("_", ":"), 
            replacement = c("-", "-"),
            vectorize_all = FALSE
          )
        )
      )
      blacklist.genotype <- suppressWarnings(
        plyr::colwise(as.character, exclude = NA)(blacklist.genotype)
      )
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
      
      if ("CHROM" %in% columns.names.blacklist.genotype) {
        columns.names.blacklist.genotype$CHROM <- as.character(columns.names.blacklist.genotype$CHROM)
      }
    
    if (data.type == "haplo.file") {
      blacklist.genotype <- dplyr::select(.data = blacklist.genotype, INDIVIDUALS, LOCUS)
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # control check to keep only individuals in the strata.df
    blacklist.genotype <- suppressWarnings(
      blacklist.genotype  %>% 
        dplyr::filter(INDIVIDUALS %in% strata.df$INDIVIDUALS)
    )
    # control check to keep only whitelisted markers from the blacklist of genotypes
    if (!is.null(whitelist.markers)) {
      blacklist.genotype <- blacklist.genotype
      message("Control check to keep only whitelisted markers present in the blacklist of genotypes to erase.")
      # updating the whitelist of markers to have all columns that id markers
      if (data.type == "vcf.file") {
        whitelist.markers.ind <- input %>% dplyr::distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      } else {
        whitelist.markers.ind <- input %>% dplyr::distinct(LOCUS, INDIVIDUALS)
      }
      
      # updating the blacklist.genotype
      blacklist.genotype <- suppressWarnings(dplyr::semi_join(whitelist.markers.ind, blacklist.genotype, by = columns.names.blacklist.genotype))
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # Create a column names
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    
    # Add one column that will allow to include the blacklist in the dataset 
    # by x column(s) of markers
    blacklist.genotype <- dplyr::mutate(.data = blacklist.genotype, ERASE = rep("erase", n()))
    
    input <- suppressWarnings(
      input %>%
        dplyr::full_join(blacklist.genotype, by = columns.names.blacklist.genotype) %>%
        dplyr::mutate(ERASE = stringi::stri_replace_na(str = ERASE, replacement = "ok"))
    )
    
    input <- input %>% 
      dplyr::mutate(GT = ifelse(ERASE == "erase", "000000", GT)) %>% 
      dplyr::select(-ERASE)
  } # End erase genotypes
  
  # dump unused object
  blacklist.id <- whitelist.markers <- whitelist.markers.ind <- blacklist.genotype <- NULL
  
  # SNP LD  --------------------------------------------------------------------
  if (!is.null(snp.ld)) {
    if (!tibble::has_name(input, "POS")) {
      stop("snp.ld is only available for VCF file, use stackr package for 
             haplotype file and create a whitelist, for other file type, use 
             PLINK linkage disequilibrium based SNP pruning option")
    }
    message("Minimizing LD...")
    snp.locus <- input %>% dplyr::distinct(LOCUS, POS)
    
    # Random selection
    if (snp.ld == "random") {
      snp.select <- snp.locus %>%
        dplyr::group_by(LOCUS) %>%
        sample_n(size = 1, replace = FALSE)
      message(stringi::stri_join("Number of original SNP = ", dplyr::n_distinct(snp.locus$POS), "\n", "Number of SNP randomly selected to keep 1 SNP per read/haplotype = ", dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS)))
    }
    
    # Fist SNP on the read
    if (snp.ld == "first") {
      snp.select <- snp.locus %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::summarise(POS = min(POS))
      message(stringi::stri_join("Number of original SNP = ", dplyr::n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS)))
    }
    
    # Last SNP on the read
    if (snp.ld == "last") {
      snp.select <- snp.locus %>%
        dplyr::group_by(LOCUS) %>%
        dplyr::summarise(POS = max(POS))
      message(stringi::stri_join("Number of original SNP = ", dplyr::n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS)))
    }
    
    # filtering the VCF to minimize LD
    input <- input %>% dplyr::semi_join(snp.select, by = c("LOCUS", "POS"))
    message("Filtering the tidy VCF to minimize LD by keeping only 1 SNP per short read/haplotype")
  } # End of snp.ld control
  
  # Unique markers id ----------------------------------------------------------
  # we want to keep LOCUS in the vcf, but not in the other type of input file
  # if (data.type != "vcf.file") {
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input), 
      pattern = "LOCUS", 
      replacement = "MARKERS", 
      vectorize_all = FALSE
    )
  }
  
  # Removing special characters in markers id ----------------------------------
  input <- input %>% 
    dplyr::mutate(
      MARKERS = stringi::stri_replace_all_fixed(
        str = as.character(MARKERS), 
        pattern = c("/", ":", "-", "."), 
        replacement = "_", 
        vectorize_all = FALSE)
    )
  
  # Markers in common between all populations (optional) -----------------------
  if (common.markers) { # keep only markers present in all pop
    input <- keep_common_markers(input)
  } # End common markers
  
  # Removing monomorphic markers------------------------------------------------
  if (monomorphic.out & !biallelic) {
    message("Removing monomorphic markers: yes")
    mono.out <- discard_monomorphic_markers(input)
    mono.markers <- mono.out$blacklist.momorphic.markers
    if (dplyr::n_distinct(mono.markers$MARKERS) > 0) {
      if (data.type == "haplo.file") mono.markers <- dplyr::rename(.data = mono.markers, LOCUS = MARKERS)
      input <- mono.out$input
      readr::write_tsv(mono.markers, "blacklist.momorphic.markers.tsv")
    }
  } # End monomorphic out
  
  # Minor Allele Frequency filter ----------------------------------------------
  # maf.thresholds <- c(0.05, 0.1) # test
  if (!is.null(maf.thresholds)) { # with MAF
    maf.info <- stackr_maf_module(
      data = input, 
      maf.thresholds = maf.thresholds, 
      maf.pop.num.threshold = maf.pop.num.threshold, 
      maf.approach = maf.approach, 
      maf.operator = maf.operator
    )
    
    input <- maf.info$input
    # maf.data <- maf.info$maf.data
    maf.info <- NULL
  } # End of MAF filters
  
  # Write to working directory
  if (!is.null(filename)) {
    message(stringi::stri_join("Writing the tidy data to the working directory: \n"), filename)
    readr::write_tsv(x = input, path = filename, col_names = TRUE)
  }
  
  res <- input
  return(res)
} # tidy genomic data
