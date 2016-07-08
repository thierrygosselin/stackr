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
#' plink, stacks haplotype file, genind, genepop, 
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
#' @import stringi
#' @import adegenet
#' @import dplyr
#' @importFrom stats var median quantile
#' @importFrom purrr map
#' @importFrom purrr flatten
#' @importFrom purrr keep
#' @importFrom purrr discard
#' @importFrom data.table fread
#' @importFrom data.table melt.data.table
#' @importFrom data.table as.data.table


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

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("DP", "AD", "vcf.headers", "GT_VCF", "INDIVIDUALS2", "ALLELE_REF_DEPTH",
      "ALLELE_ALT_DEPTH", "GT_BIN", "GT_HAPLO")
  )
}

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
  filename = NULL,
  ...) {
  
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  
  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  if (!is.null(pop.labels)) {
    pop.labels <- stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  if (!is.null(pop.select)) {
    pop.select <- stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  # File type detection ********************************************************
  if(is.genind(data)){
    data.type <- "genind.file"
    message("File type: genind object")
  } else {
    data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
    if (identical(data.type, "##fileformat=VCF") | stri_detect_fixed(str = data, pattern = ".vcf")) {
      data.type <- "vcf.file"
      message("File type: VCF")
    }
    if (stri_detect_fixed(str = data, pattern = ".tped")) {
      data.type <- "plink.file"
      message("File type: PLINK")
      if (!file.exists(stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        stop("Missing tfam file with the same prefix as your tped")
      }
    } 
    if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stri_detect_fixed(str = data.type, pattern = "MARKERS")| stri_detect_fixed(str = data.type, pattern = "LOCUS")) {
      data.type <- "df.file"
      message("File type: data frame of genotypes")
    }
    
    
    if (stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      data.type <- "haplo.file"
      message("File type: haplotypes from stacks")
      # if (is.null(blacklist.genotype)) {
      #   stop("blacklist.genotype file missing. 
      #        Use stackr's missing_genotypes function to create this blacklist")
      # }
    }
    if (stri_detect_fixed(str = data, pattern = ".gen")) {
      data.type <- "genepop.file"
      message("File type: genepop")
    } 
  } # end file type detection
  
  if(data.type == "haplo.file") {
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
  
  # Strata argument required for VCF and haplotypes files **********************
  if (data.type == "haplo.file" | data.type == "vcf.file") {
    if (is.null(strata)) stop("strata argument is required")
  }
  
  # Import whitelist of markers ************************************************
  if (is.null(whitelist.markers)) { # no Whitelist
    message("Whitelist of markers: no")
  } else { # with Whitelist of markers
    message("Whitelist of markers: yes")
    whitelist.markers <- read_tsv(whitelist.markers, col_names = TRUE)
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
      whitelist.markers <- select(.data = whitelist.markers, LOCUS)
      columns.names.whitelist <- colnames(whitelist.markers)
    }
    # # for df.file, plink.file, genepop.file and genind objct
    # if (data.type %in% c("genind.file", "plink.file", "df.file", "genepop.file") {
    #   whitelist.markers <- whitelist.markers %>% select(MARKERS = LOCUS)
    #   columns.names.whitelist <- colnames(whitelist.markers)
    # }
  }
  
  # Import blacklist id ********************************************************
  if (is.null(blacklist.id)) { # No blacklist of ID
    message("Blacklisted individuals: no")
  } else { # With blacklist of ID
    message("Blacklisted individuals: yes")
    blacklist.id <- read_tsv(blacklist.id, col_names = TRUE)
  }
  
  # population levels and strata ***********************************************
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      message("strata file: yes")
      number.columns.strata <- max(count.fields(strata, sep = "\t"))
      col.types <- stri_paste(rep("c", number.columns.strata), collapse = "")
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>% 
        rename(POP_ID = STRATA)
    } else {
      message("strata object: yes")
      colnames(strata) <- stri_replace_all_fixed(str = colnames(strata), 
                                                 pattern = "STRATA", 
                                                 replacement = "POP_ID", 
                                                 vectorize_all = FALSE
      )
      strata.df <- strata
    }
    
    # filtering the strata if blacklist id available
    if (!is.null(blacklist.id)) {
      strata.df <- anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
    }
    strata.df$POP_ID <- stri_replace_all_fixed(strata.df$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)
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
      as_data_frame() %>% 
      rename(LOCUS = ID, CHROM = `#CHROM`) %>%
      mutate(
        CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1"),
        POS = as.character(POS),
        LOCUS = as.character(LOCUS)
      )
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Detect the format fields
    format.headers <- unlist(stri_split_fixed(str = input$FORMAT[1], pattern = ":"))
    input <- input %>% select(-FORMAT) # no longer needed
    
    # Tidying the VCF to make it easy to work on the data for conversion
    message("Making the VCF population wise")
    # input <- input %>%
    # tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, LOCUS, POS, REF, ALT)) # Gather individuals in 1 colummn
    
    # filter blacklisted individuals
    
    input <- data.table::melt.data.table(
      data = as.data.table(input), 
      id.vars = c("CHROM", "LOCUS", "POS", "REF", "ALT"), 
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = "FORMAT_ID"
    ) %>% 
      as_data_frame() %>% 
      mutate(
        INDIVIDUALS = stri_replace_all_fixed(
          str = INDIVIDUALS, 
          pattern = c("_", ":"), 
          replacement = c("-", "-"), 
          vectorize_all = FALSE)
      ) %>% 
      filter(!INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
    
    # population levels and strata  --------------------------------------------
    # if (!is.null(blacklist.id)) {
    #   strata.df <- anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
    # }
    input <- left_join(x= input, y = strata.df, by = "INDIVIDUALS")
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # Tidy VCF
    message("Tidying the vcf...")
    input <- tidyr::separate(
      data = input, FORMAT_ID, format.headers, sep = ":", extra = "drop"
    )
    
    if (vcf.metadata) {
      # GL cleaning
      if (length(unlist(stri_extract_all_fixed(str = format.headers, pattern = "GL", omit_no_match = TRUE))) > 0) {
        message("Fixing GL column...")
        input <- input %>% 
          mutate( 
            GL = suppressWarnings(as.numeric(stri_replace_all_fixed(GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all = FALSE)))
          )
      } # end cleaning GL column
      
      # Cleaning DP and changing name to READ_DEPTH
      if (length(unlist(stri_extract_all_fixed(str = format.headers, pattern = "DP", omit_no_match = TRUE))) > 0) {
        message("Fixing DP (READ_DEPTH) column...")
        input <- input %>% 
          rename(READ_DEPTH = DP) %>% 
          mutate(
            READ_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(READ_DEPTH, "^0$", "NA", vectorize_all=FALSE))),
            READ_DEPTH = ifelse(GT == "./.", NA, READ_DEPTH)
          )
      } # end cleaning READ_DEPTH (DP) column
      
      # Cleaning AD (ALLELES_DEPTH)
      if (length(unlist(stri_extract_all_fixed(str = format.headers, pattern = "AD", omit_no_match = TRUE))) > 0) {
        message("Splitting AD columns into allele coverage info...")
        input <- input %>% 
          tidyr::separate(AD, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"), sep = ",", extra = "drop") %>% 
          mutate(
            ALLELE_REF_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(ALLELE_REF_DEPTH, "^0$", "NA", vectorize_all = TRUE))),
            ALLELE_ALT_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(ALLELE_ALT_DEPTH, "^0$", "NA", vectorize_all = TRUE)))
            # Mutate coverage ratio for allelic imbalance
            # ALLELE_COVERAGE_RATIO = suppressWarnings(
            # as.numeric(ifelse(GT == "./." | GT == "0/0" | GT == "1/1", "NA",
            # ((ALLELE_ALT_DEPTH - ALLELE_REF_DEPTH)/(ALLELE_ALT_DEPTH + ALLELE_REF_DEPTH)))))
          )
      } # end cleaning AD column
    } else {
      input <- input %>% 
        select(CHROM, LOCUS, POS, REF, ALT, POP_ID, INDIVIDUALS, GT)
    }# end metadata section
    
    # recoding genotype and creating a new column combining CHROM, LOCUS and POS 
    input <- input %>%
      tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE) %>%
      mutate(
        REF= stri_replace_all_fixed(
          str = REF, 
          pattern = c("A", "C", "G", "T"), 
          replacement = c("001", "002", "003", "004"), 
          vectorize_all = FALSE
        ), # replace nucleotide with numbers
        ALT = stri_replace_all_fixed(
          str = ALT, pattern = c("A", "C", "G", "T"), 
          replacement = c("001", "002", "003", "004"), 
          vectorize_all = FALSE
        ),# replace nucleotide with numbers
        GT_VCF = GT,
        GT = ifelse(GT == "0/0", stri_join(REF, REF, sep = ""),
                    ifelse(GT == "1/1",  stri_join(ALT, ALT, sep = ""),
                           ifelse(GT == "0/1", stri_join(REF, ALT, sep = ""),
                                  ifelse(GT == "1/0", stri_join(ALT, REF, sep = ""), "000000")
                           )
                    )
        )
      ) %>% 
      as_data_frame()
    
    # Experimental (for genlight object)
    input$GT_BIN <- stri_replace_all_fixed(str = input$GT_VCF, pattern = c("0/0", "1/1", "0/1", "1/0", "./."), replacement = c("0", "2", "1", "1", NA), vectorize_all = FALSE)
    
    # re-computing the REF/ALT allele-----------------------------------------------
    if (!is.null(pop.select) || !is.null(blacklist.id)) {
      message("Adjusting REF/ALT alleles to account for filters...")
      
      input.select <- input %>%
        select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
        #faster than: tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
        mutate(
          A1 = stri_sub(str = GT, from = 1, to = 3),
          A2 = stri_sub(str = GT, from = 4, to = 6)
        ) %>%
        select(-GT)
      
      new.ref.alt.alleles <- input.select %>% 
        tidyr::gather(
          data = ., key = ALLELES, 
          value = GT, 
          -c(MARKERS, INDIVIDUALS, POP_ID)
        ) %>% # just miliseconds longer than data.table.melt so keeping this one for simplicity
        filter(GT != "000") %>%  # remove missing "000"
        group_by(MARKERS, GT) %>%
        tally %>% 
        ungroup() %>% 
        mutate(ALLELE_NUMBER = rep(c("A1", "A2"), each = 1, times = n()/2)) %>% 
        group_by(MARKERS) %>%
        mutate(
          PROBLEM = ifelse(n[ALLELE_NUMBER == "A1"] == n[ALLELE_NUMBER == "A2"], "equal_number", "ok"),
          GROUP = ifelse(n == max(n), "REF", "ALT")
        ) %>% 
        group_by(MARKERS, GT) %>% 
        mutate(
          ALLELE = ifelse(PROBLEM == "equal_number" & ALLELE_NUMBER == "A1", "REF", 
                          ifelse(PROBLEM == "equal_number" & ALLELE_NUMBER == "A2", "ALT", 
                                 GROUP)
          )
        ) %>% 
        select(MARKERS, GT, ALLELE) %>%
        group_by(MARKERS) %>% 
        tidyr::spread(data = ., ALLELE, GT) %>% 
        select(MARKERS, REF, ALT)
      
      old.ref.alt.alleles <- distinct(.data = input, MARKERS, REF, ALT)
      
      ref.alt.alleles <- full_join(
        old.ref.alt.alleles,
        new.ref.alt.alleles %>% 
          rename(REF_NEW = REF, ALT_NEW = ALT)
        , by = "MARKERS") %>% 
        mutate(
          REF_ALT_CHANGE = if_else(REF == REF_NEW, "identical", "different")
          ) %>% 
        filter(REF_ALT_CHANGE == "different")
      
      message(stri_paste("Number of markers with REF/ALT change = ", length(ref.alt.alleles$MARKERS)))
      
      ref.alt.alleles.change <- full_join(input.select, new.ref.alt.alleles, by = "MARKERS") %>% 
        mutate(
          A1 = replace(A1, which(A1 == "000"), NA),
          A2 = replace(A2, which(A2 == "000"), NA),
          GT_VCF_A1 = if_else(A1 == REF, "0", "1", missing = "."),
          GT_VCF_A2 = if_else(A2 == REF, "0", "1", missing = "."),
          GT_VCF = stri_paste(GT_VCF_A1, GT_VCF_A2, sep = "/"),
          GT_BIN = stri_replace_all_fixed(
            str = GT_VCF, 
            pattern = c("0/0", "1/1", "0/1", "1/0", "./."), 
            replacement = c("0", "2", "1", "1", NA), 
            vectorize_all = FALSE
          )
        ) %>% 
        select(MARKERS, INDIVIDUALS, REF, ALT, GT_VCF, GT_BIN)
      
      input <- left_join(
        input %>% 
          select(-c(REF, ALT, GT_VCF, GT_BIN)), 
        ref.alt.alleles.change, 
        by = c("MARKERS", "INDIVIDUALS")
        )
      
      # remove unused object
      input.select <- NULL
      new.ref.alt.alleles <- NULL 
      old.ref.alt.alleles <- NULL 
      ref.alt.alleles <- NULL 
      ref.alt.alleles.change <- NULL 
    } # end re-computing the REF/ALT allele
    
    # Re ordering columns-------------------------------------------------------
    if (vcf.metadata) {
      # Re order columns
      common.colnames <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS", "GT_VCF", "GT_BIN", "REF", "ALT")
      vcf.headers <- colnames(input)
      metadata.colnames <- purrr::discard(.x = colnames(input), .p = vcf.headers %in% common.colnames)
      input <- input[c(common.colnames, metadata.colnames)]
      
    } else {
      input <- input %>% select(MARKERS, CHROM, LOCUS, POS, POP_ID, INDIVIDUALS, GT_VCF, GT_BIN, REF, ALT, GT)
    }
  } # End import VCF
  
  # Import PLINK ---------------------------------------------------------------
  if (data.type == "plink.file") { # PLINK
    message("Importing the PLINK files...")
    tfam <- data.table::fread(
      input = stri_replace_all_fixed(
        str = data, 
        pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE
      ),
      sep = " ", 
      header = FALSE, 
      stringsAsFactors = FALSE,
      verbose = FALSE,
      select = c(1,2),
      colClasses=list(character = c(1,2)),
      col.names = c("POP_ID", "INDIVIDUALS"),
      showProgress = TRUE, 
      data.table = FALSE) %>% 
      as_data_frame() %>%
      mutate(# remove "_" in individual name and replace with "-"
        INDIVIDUALS = stri_replace_all_fixed( str = INDIVIDUALS, 
                                              pattern = c("_", ":"), 
                                              replacement = c("-", "-"), 
                                              vectorize_all = FALSE),
        POP_ID = stri_replace_all_fixed(POP_ID, 
                                        pattern = " ", 
                                        replacement = "_", 
                                        vectorize_all = FALSE)
      )
    
    # if no strata tfam = strata.df
    if (is.null(strata)) {
      strata.df <- tfam
    } else {
      strata.df <- mutate(.data = strata.df, 
                          INDIVIDUALS = stri_replace_all_fixed(
                            str = INDIVIDUALS, 
                            pattern = c("_", ":"), 
                            replacement = c("-", "-"),
                            vectorize_all = FALSE
                          )
      )
    }
    
    tped.header.prep <- tfam %>% 
      select(INDIVIDUALS) %>%
      mutate(
        NUMBER = seq(1,n()),
        ALLELE1 = rep("A1", n()), ALLELE2 = rep("A2", n())
      ) %>%
      tidyr::gather(ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, NUMBER)) %>%
      arrange(NUMBER) %>% 
      select(-ALLELES_GROUP) %>% 
      tidyr::unite(INDIVIDUALS_ALLELES, c(INDIVIDUALS, ALLELES), sep = "_", remove = FALSE) %>% 
      arrange(NUMBER) %>% 
      mutate(NUMBER = seq(from = (1 + 4), to = n() + 4)) %>% 
      select(-ALLELES)
    
    tped.header.names <- c("LOCUS", tped.header.prep$INDIVIDUALS_ALLELES)
    tped.header.integer <- c(2, tped.header.prep$NUMBER)
    
    if (!is.null(blacklist.id)) { # using the blacklist of individuals
      blacklist.id <- mutate(
        .data = blacklist.id, 
        INDIVIDUALS = stri_replace_all_fixed(# remove "_" in individual name and replace with "-"
          str = INDIVIDUALS, 
          pattern = c("_", ":"), replacement = c("-", "-"), vectorize_all = FALSE
        )
      )
      
      whitelist.id <- tped.header.prep %>% 
        anti_join(blacklist.id, by = "INDIVIDUALS") %>% 
        arrange(NUMBER)
      tped.header.names <- c("LOCUS", whitelist.id$INDIVIDUALS_ALLELES)
      tped.header.integer <- c(2, whitelist.id$NUMBER)
      
      strata.df <- anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
    }
    
    input <- data.table::fread( # import PLINK
      input = data, 
      sep = " ", 
      header = FALSE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      select = tped.header.integer,
      col.names = tped.header.names,
      showProgress = TRUE,
      data.table = FALSE) %>% 
      as_data_frame() %>% 
      mutate(LOCUS = as.character(LOCUS))
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(
        semi_join(input, whitelist.markers, by = columns.names.whitelist)
      )
    }
    
    # To reduce the size of the dataset we subsample the markers with max.marker
    if (!is.null(max.marker)) {
      message("Using the max.marker to reduce the size of the dataset")
      input <- sample_n(tbl = input, size = max(as.numeric(max.marker)), replace = FALSE)
      
      max.marker.subsample.select <- input %>% 
        distinct(LOCUS, .keep_all = TRUE) %>% 
        arrange(LOCUS)
      
      write_tsv(# save results
        x = max.marker.subsample.select, 
        path = "max.marker.subsample.select.tsv")
    }
    
    # Make tidy
    message("Tidying the PLINK file ...")
    # Filling GT and new separating INDIVIDUALS from ALLELES
    # combining alleles
    input <- data.table::melt.data.table(
      data = as.data.table(input), 
      id.vars = "LOCUS", 
      variable.name = "INDIVIDUALS_ALLELES", 
      value.name = as.character("GT"), 
      variable.factor = FALSE, 
      value.factor = FALSE
    ) %>% 
      as_data_frame() %>% 
      mutate(GT = stri_pad_left(str = GT, width = 3, pad = "0")) %>% 
      tidyr::separate(data = ., col = INDIVIDUALS_ALLELES, into = c("INDIVIDUALS", "ALLELES"), sep = "_")
    
    input <- data.table::dcast.data.table(
      data = as.data.table(input), 
      formula = LOCUS + INDIVIDUALS ~ ALLELES, 
      value.var = "GT"
    ) %>% 
      as_data_frame() %>% 
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>% 
      select(LOCUS, INDIVIDUALS, GT)
    
    # population levels and strata  ----------------------------------------------
    message("Integrating the tfam/strata file...")
    
    input <- left_join(x= input, y = strata.df, by = "INDIVIDUALS")
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # removing untyped markers across all-pop
    remove.missing.gt <- input %>%
      select(LOCUS, GT) %>%
      filter(GT != "000000")
    
    untyped.markers <- n_distinct(input$LOCUS) - n_distinct(remove.missing.gt$LOCUS)
    if (untyped.markers > 0) {
      message(paste0("Number of marker with 100 % missing genotypes: ", untyped.markers))
      input <- suppressWarnings(
        semi_join(input, 
                  remove.missing.gt %>% 
                    distinct(LOCUS, .keep_all = TRUE), 
                  by = "LOCUS")
      )
    }
    
    # Unused objects
    tped.header.prep <- NULL
    tped.header.integer <- NULL
    tped.header.names <- NULL
    remove.missing.gt <- NULL
    
  } # End import PLINK
  
  # Import DF-------------------------------------------------------------------
  if (data.type == "df.file") { # DATA FRAME OF GENOTYPES
    input <- read_long_tidy_wide(data = data)
    
    if ("MARKERS" %in% colnames(input)) {
      input <- rename(.data = input, LOCUS = MARKERS)
    }
    
    
    # Change individuals names containing special character
    input$INDIVIDUALS <- stri_replace_all_fixed(
      str = input$INDIVIDUALS, 
      pattern = c("_", ":"), 
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      blacklist.id$INDIVIDUALS <- stri_replace_all_fixed(
        str = blacklist.id$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
      input <- suppressWarnings(anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    # population levels and strata  --------------------------------------------
    if (!is.null(strata)) {
      strata.df$INDIVIDUALS <- stri_replace_all_fixed(
        str = strata.df$INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
      
      input <- input %>%
        select(-POP_ID) %>% 
        left_join(strata.df, by = "INDIVIDUALS")
    }
    
    # Change potential problematic POP_ID space
    input$POP_ID = stri_replace_all_fixed(input$POP_ID, 
                                          pattern = " ", 
                                          replacement = "_", 
                                          vectorize_all = FALSE)
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
  } # End import data frame of genotypes
  
  # Import haplo---------------------------------------------------------------
  if (data.type == "haplo.file") { # Haplotype file
    message("Importing the stacks haplotype file")
    number.columns <- max(count.fields(data, sep = "\t"))
    
    input <- data.table::fread(
      input = data, 
      sep = "\t", 
      header = TRUE, 
      stringsAsFactors = FALSE,
      colClasses=list(character=1:number.columns),
      verbose = FALSE,
      showProgress = TRUE,
      data.table = FALSE, 
      na.strings = "-"
    ) %>% 
      as_data_frame() %>% 
      select(-Cnt) %>% 
      rename(LOCUS = `Catalog ID`)
    
    input <- data.table::melt.data.table(
      data = as.data.table(input), 
      id.vars = "LOCUS", 
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = "GT"
    ) %>% 
      as_data_frame()
    
    number.columns <- NULL
    
    # remove consensus markers -------------------------------------------------
    message("Scanning for consensus markers")
    consensus.markers <- input %>%
      filter(GT == "consensus") %>% 
      distinct(LOCUS, .keep_all = TRUE)
    
    if (length(consensus.markers$LOCUS) > 0) {
      input <- suppressWarnings(anti_join(input, consensus.markers, by = "LOCUS"))
    }
    message(stri_paste("Consensus markers removed: ", n_distinct(consensus.markers$LOCUS)))
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      input <- suppressWarnings(anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    # population levels and strata  --------------------------------------------
    input <- left_join(x= input, y = strata.df, by = "INDIVIDUALS")
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # removing errors and potential paralogs (GT with > 2 alleles) -------------
    message("Scanning for genotypes with more than 2 alleles")
    input <- input %>%
      mutate(POLYMORPHISM = stri_count_fixed(GT, "/"))
    
    blacklist.paralogs <- input %>%
      filter(POLYMORPHISM > 1) %>% 
      select(LOCUS, INDIVIDUALS)
    
    message(stri_paste("Number of genotypes with more than 2 alleles: ", length(blacklist.paralogs$LOCUS), sep = ""))
    
    if (length(blacklist.paralogs$LOCUS) > 0){
      input <- input %>% 
        mutate(GT = ifelse(POLYMORPHISM >1, NA, GT)) %>% 
        select(-POLYMORPHISM)
      
      write_tsv(blacklist.paralogs, "blacklist.genotypes.paralogs.tsv")
    }
    
    tidy.haplo <- input %>% select(LOCUS, INDIVIDUALS, GT_HAPLO = GT)
    
    # recode genotypes
    message("Recoding haplotypes in 6 digits format")
    # input.bk <- input
    input <- suppressWarnings(
      input %>% 
        tidyr::separate(
          col = GT, into = c("A1", "A2"), sep = "/", extra = "drop", remove = TRUE
        ) %>%
        mutate(A2 = ifelse(is.na(A2), A1, A2)) %>% 
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(LOCUS, POP_ID, INDIVIDUALS))%>% 
        tidyr::spread(data = ., key = LOCUS, value = GT)
    )
    
    input.id.col <- select(.data = input, POP_ID, INDIVIDUALS, ALLELES)
    
    input.variable <- input %>% select(-c(POP_ID, INDIVIDUALS, ALLELES)) %>% 
      colwise(factor, exclude = NA)(.) %>% 
      colwise(as.numeric)(.)
    
    input <- bind_cols(input.id.col, input.variable)
    
    input <- data.table::melt.data.table(
      data = as.data.table(input), 
      id.vars = c("INDIVIDUALS", "POP_ID", "ALLELES"), 
      variable.name = "LOCUS",
      variable.factor = FALSE,
      value.name = "GT"
    ) %>% 
      as_data_frame() %>% 
      mutate(
        GT = as.character(GT),
        GT = stri_pad_left(str = GT, width = 3, pad = "0"),
        GT = stri_replace_na(str = GT, replacement = "000")
      ) %>% 
      tidyr::spread(data = ., key = ALLELES, value = GT) %>% 
      tidyr::unite(data = ., GT, A1, A2, sep = "") 
    
    # input <- arrange(.data = input, LOCUS, POP_ID, INDIVIDUALS) %>% filter(GT == "000000")
    # tidy.haplo <- arrange(.data = tidy.haplo, LOCUS, INDIVIDUALS)
    
    input <- inner_join(input, tidy.haplo, by = c("LOCUS", "INDIVIDUALS")) %>%
      select(LOCUS, POP_ID, INDIVIDUALS, GT, GT_HAPLO) %>%
      arrange(LOCUS, POP_ID, INDIVIDUALS, GT, GT_HAPLO)
    
    # change the filename and strata.df here
  } # End import haplotypes file
  
  # Import genepop **************************************************************
  if (data.type == "genepop.file") {
    message("Tidying the genepop file ...")
    
    # data <- "/Users/thierry/Documents/skipjack/skipjack.filtered_imputed.gen"
    data <- adegenet::read.genepop(data, ncode = 3, quiet = TRUE)
    data.type <- "genind.file"
    genind.type <- "genepop"
  }
  
  # Import GENIND **************************************************************
  # load("/Users/thierry/Documents/skipjack/fst.test.RData")
  if (data.type == "genind.file") { # DATA FRAME OF GENOTYPES
    # data = skipjack.genind
    input <- adegenet::genind2df(data) %>% 
      add_rownames("INDIVIDUALS") %>% 
      rename(POP_ID = pop)
    
    if (genind.type == "genepop") {
      input <- tidyr::gather(
        data = input, key = LOCUS, value = GT, -c(INDIVIDUALS, POP_ID)
      )
    } else {
      message("Tidying the genind object ...")
      input <- input %>% 
        tidyr::gather(key = LOCUS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% 
        tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 1, remove = TRUE, extra = "drop") %>% 
        mutate(
          A1 = stri_pad_left(str= A1, pad = "0", width = 3),
          A2 = stri_pad_left(str= A2, pad = "0", width = 3)
        ) %>% 
        tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>% 
        mutate(GT = replace(GT, which(GT == "NANA"), "000000"))
    }
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      input <- suppressWarnings(anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    # population levels and strata  --------------------------------------------
    if (!is.null(strata)) {
      input <- input %>%
        select(-POP_ID) %>% 
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        left_join(strata.df, by = "INDIVIDUALS")
    }
    
    # Change potential problematic POP_ID space
    input$POP_ID = stri_replace_all_fixed(input$POP_ID, 
                                          pattern = " ", 
                                          replacement = "_", 
                                          vectorize_all = FALSE)
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # Now the genind and genepop are like ordinary data frames
    data.type <- "df.file" # for subsequent steps
    
  } # End tidy genind
  
  # Arrange the id and create a strata after pop select ************************
  input <- input %>% 
    mutate(
      INDIVIDUALS = stri_replace_all_fixed(
        str = INDIVIDUALS, 
        pattern = c("_", ":"), 
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
    )
  
  strata.df <- input %>%
    ungroup() %>%
    distinct(POP_ID, INDIVIDUALS)
  
  # Blacklist genotypes ********************************************************
  if (is.null(blacklist.genotype)) { # no Whitelist
    message("Erasing genotype: no")
  } else {
    message("Erasing genotype: yes")
    blacklist.genotype <- read_tsv(blacklist.genotype, col_names = TRUE) %>% 
      mutate(
        INDIVIDUALS = stri_replace_all_fixed(
          str = INDIVIDUALS, 
          pattern = c("_", ":"), 
          replacement = c("-", "-"),
          vectorize_all = FALSE
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
      blacklist.genotype <- select(.data = blacklist.genotype, INDIVIDUALS, LOCUS)
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # control check to keep only individuals in the strata.df
    blacklist.genotype <- suppressWarnings(
      blacklist.genotype  %>% 
        filter(INDIVIDUALS %in% strata.df$INDIVIDUALS)
    )
    # control check to keep only whitelisted markers from the blacklist of genotypes
    if (!is.null(whitelist.markers)) {
      blacklist.genotype <- blacklist.genotype
      message("Control check to keep only whitelisted markers present in the blacklist of genotypes to erase.")
      # updating the whitelist of markers to have all columns that id markers
      if (data.type == "vcf.file"){
        whitelist.markers.ind <- input %>% distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      } else {
        whitelist.markers.ind <- input %>% distinct(LOCUS, INDIVIDUALS)
      }
      
      # updating the blacklist.genotype
      blacklist.genotype <- suppressWarnings(semi_join(whitelist.markers.ind, blacklist.genotype, by = columns.names.blacklist.genotype))
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # Create a column names
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    
    # Add one column that will allow to include the blacklist in the dataset 
    # by x column(s) of markers
    blacklist.genotype <- mutate(.data = blacklist.genotype, ERASE = rep("erase", n()))
    
    input <- suppressWarnings(
      input %>%
        full_join(blacklist.genotype, by = columns.names.blacklist.genotype) %>%
        mutate(ERASE = stri_replace_na(str = ERASE, replacement = "ok"))
    )
    
    input <- input %>% 
      mutate(GT = ifelse(ERASE == "erase", "000000", GT)) %>% 
      select(-ERASE)
  } # End erase genotypes
  
  
  # population levels **********************************************************
  if(is.null(pop.levels)) { # no pop.levels
    input <- mutate(.data = input, POP_ID = factor(POP_ID))
  } else { # with pop.levels
    input <- mutate(
      .data = input,
      POP_ID = factor(
        stri_replace_all_regex(
          POP_ID, stri_paste("^", pop.levels, "$", sep = ""), pop.labels, vectorize_all = FALSE),
        levels = unique(pop.labels), ordered = TRUE
      )
    )
  }
  
  # dump unused object
  blacklist.id <- NULL
  whitelist.markers <- NULL
  whitelist.markers.ind <- NULL
  blacklist.genotype <- NULL
  
  # LD control... keep only 1 SNP per haplotypes/reads (optional) ************
  if (!is.null(snp.ld)) {
    if (data.type != "vcf.file") {
      stop("snp.ld is only available for VCF file, use stackr package for 
             haplotype file and create a whitelist, for other file type, use 
             PLINK linkage disequilibrium based SNP pruning option")
    }
    message("Minimizing LD...")
    snp.locus <- input %>% distinct(LOCUS, POS)
    # Random selection
    if (snp.ld == "random") {
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>%
        sample_n(size = 1, replace = FALSE)
      message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP randomly selected to keep 1 SNP per read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # Fist SNP on the read
    if (snp.ld == "first") {
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>%
        summarise(POS = min(POS))
      message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # Last SNP on the read
    if (snp.ld == "last") {
      snp.select <- snp.locus %>%
        group_by(LOCUS) %>%
        summarise(POS = max(POS))
      message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
    }
    
    # filtering the VCF to minimize LD
    input <- input %>% semi_join(snp.select, by = c("LOCUS", "POS"))
    message("Filtering the tidy VCF to minimize LD by keeping only 1 SNP per short read/haplotype")
  } # End of snp.ld control
  
  # Unique markers id *********************************************************
  # we want to keep LOCUS in the vcf, but not in the other type of input file
  if (data.type != "vcf.file") {
    colnames(input) <- stri_replace_all_fixed(str = colnames(input), 
                                              pattern = "LOCUS", 
                                              replacement = "MARKERS", 
                                              vectorize_all = FALSE)
  }
  
  # Removing special characters in markers id ******************************
  input <- input %>% 
    mutate(
      MARKERS = stri_replace_all_fixed(
        str = as.character(MARKERS), 
        pattern = c("/", ":", "-", "."), 
        replacement = "_", 
        vectorize_all = FALSE)
    )
  
  # Markers in common between all populations (optional) *********************
  if (common.markers) { # keep only markers present in all pop
    message("Using markers common in all populations:")
    pop.number <- n_distinct(input$POP_ID)
    
    pop.filter <- input %>% filter(GT != "000000")
    
    pop.filter <- pop.filter %>% 
      group_by(MARKERS) %>%
      filter(n_distinct(POP_ID) == pop.number) %>%
      arrange(MARKERS) %>%
      distinct(MARKERS)
    
    markers.input <- n_distinct(input$MARKERS)
    markers.in.common <- n_distinct(pop.filter$MARKERS)
    blacklist.markers.common <- markers.input - markers.in.common
    
    message(stri_join("Number of original markers = ", markers.input, 
                      "\n", "Number of markers present in all the populations = ", 
                      markers.in.common, "\n", 
                      "Number of markers removed = ", 
                      blacklist.markers.common)
    )
    
    if (blacklist.markers.common > 0) {
      input <- suppressWarnings(input %>% semi_join(pop.filter, by = "MARKERS"))
    }
    pop.filter <- NULL # ununsed object
    markers.input <- NULL
    markers.in.common <- NULL
    blacklist.markers.common <- NULL
  } # End common markers
  
  # Removing monomorphic markers------------------------------------------------
  if (monomorphic.out) {
    message("Removing monomorphic markers: yes")
    message("Scanning for monomorphic markers...")
    
    mono.markers <- input %>%
      select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
      mutate(
        A1 = stri_sub(GT, 1, 3),
        A2 = stri_sub(GT, 4,6)
      ) %>% 
      select(-GT)
    
    mono.markers <- data.table::melt.data.table(
      data = as.data.table(mono.markers), 
      id.vars = c("MARKERS", "INDIVIDUALS", "POP_ID"), 
      variable.name = "ALLELES",
      variable.factor = FALSE,
      value.name = "GT"
    ) %>% 
      as_data_frame() %>% 
      filter(GT != "000") %>%
      distinct(MARKERS, GT) %>% 
      select(MARKERS) %>% 
      group_by(MARKERS) %>% 
      tally %>%
      filter(n == 1) %>%
      ungroup() %>% 
      select(MARKERS)
    
    # Remove the markers from the dataset
    message(paste0("Number of monomorphic markers removed: ", n_distinct(mono.markers$MARKERS)))
    
    if (length(mono.markers$MARKERS)>0) {
      input <- anti_join(input, mono.markers, by = "MARKERS")
      if(data.type == "haplo.file") mono.markers <- rename(.data = mono.markers, LOCUS = MARKERS)
      write_tsv(mono.markers, "blacklist.momorphic.markers.tsv")
    }
  }
  
  # Minor Allele Frequency filter ********************************************
  # maf.thresholds <- c(0.05, 0.1) # test
  if (!is.null(maf.thresholds)) { # with MAF
    maf.local.threshold <- maf.thresholds[1]
    maf.global.threshold <- maf.thresholds[2]
    message("MAF filter: yes")
    
    if (data.type == "vcf.file") {
      maf.local <- input %>%
        filter(GT_VCF != "./.") %>%
        group_by(MARKERS, POP_ID, REF, ALT) %>%
        summarise(
          N = as.numeric(n()),
          PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
          QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
        ) %>%
        mutate(MAF_LOCAL = ((QQ * 2) + PQ) / (2 * N))
      
      maf.global <- maf.local %>%
        group_by(MARKERS) %>%
        summarise_each_(funs(sum), vars = c("N", "PQ", "QQ")) %>%
        mutate(MAF_GLOBAL = ((QQ * 2) + PQ) / (2 * N)) %>%
        select(MARKERS, MAF_GLOBAL)
      
      maf.data <- maf.global %>%
        left_join(maf.local, by = c("MARKERS")) %>%
        select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
      
      maf.local <- NULL
      maf.global <- NULL
    } # end maf calculations with vcf
    
    if (data.type == "plink.file" | data.type == "df.file") {
      message("Calculating global and local MAF, this may take some time on large data set")
      
      # We split the alleles here to prep for MAF
      maf.data <- input %>%
        select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        mutate(
          A1 = stri_sub(GT, 1, 3),
          A2 = stri_sub(GT, 4,6)
        ) %>% 
        select(-GT)
      
      maf.data <- data.table::melt.data.table(
        data = as.data.table(maf.data), 
        id.vars = c("MARKERS", "INDIVIDUALS", "POP_ID"), 
        variable.name = "ALLELES",
        variable.factor = FALSE,
        value.name = "GT"
      ) %>% 
        as_data_frame() %>% 
        filter(GT != "000")
      
      maf.data <- maf.data %>%
        group_by(MARKERS, GT, POP_ID) %>%
        tally %>%
        arrange(MARKERS, GT) %>% 
        group_by(MARKERS, GT) %>%
        mutate(sum.pop = sum(n)) %>% 
        group_by(MARKERS) %>%
        mutate(
          MAF_GLOBAL = min(sum.pop)/sum(n),
          ALT = ifelse(min(sum.pop), "alt", "ref")
        ) %>%
        group_by(MARKERS, POP_ID) %>%
        mutate(MAF_LOCAL = n/sum(n)) %>% 
        arrange(MARKERS, POP_ID, GT) %>% 
        group_by(MARKERS, POP_ID) %>% 
        filter(n == min(n)) %>% 
        distinct(MARKERS, POP_ID, .keep_all = TRUE) %>% 
        select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
    }# end maf calculations with PLINK or data frame of genotypes
    
    if (data.type == "haplo.file") {
      stop("MAF filtering is only available for VCF file, use stackr
             package and update your whitelist")
    }
    
    write_tsv(x = maf.data, 
              path = "maf.data.tsv",
              col_names = TRUE, 
              append = FALSE
    )
    message("The MAF table was written in your folder")
    
    # # update the vcf with the maf info
    # input <- full_join(input, maf.data, by = c("MARKERS", "POP_ID"))
    if (maf.approach == "haplotype") {
      
      vcf.maf <- tidyr::separate(data = maf.data, 
                                 col = MARKERS, 
                                 into = c("CHROM", "LOCUS", "POS"), 
                                 sep = "__", 
                                 remove = FALSE, 
                                 extra = "warn"
      )
      
      if (maf.operator == "OR") {
        vcf.maf <- vcf.maf %>%
          group_by(LOCUS, POP_ID) %>%
          summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
          group_by(LOCUS) %>%
          tally() %>%
          filter(n >= maf.pop.num.threshold) %>%
          select(LOCUS) %>%
          left_join(input, by = "LOCUS") %>%
          arrange(LOCUS, POP_ID)
      } else { # AND operator between local and global maf
        vcf.maf <- vcf.maf %>%
          group_by(LOCUS, POP_ID) %>%
          summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
          group_by(LOCUS) %>%
          tally() %>%
          filter(n >= maf.pop.num.threshold) %>%
          select(LOCUS) %>%
          left_join(input, by = "LOCUS") %>%
          arrange(LOCUS, POP_ID)
      }
      vcf.maf <- vcf.maf %>% select(-c(CHROM, LOCUS, POS))
    } # end maf haplotype approach
    
    if (maf.approach == "SNP") { # SNP approach
      if (maf.operator == "OR") {
        vcf.maf <- maf.data %>%
          group_by(MARKERS, POP_ID) %>%
          summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
          group_by(MARKERS) %>%
          tally() %>%
          filter(n >= maf.pop.num.threshold) %>%
          select(MARKERS) %>%
          left_join(input, by = "MARKERS") %>%
          arrange(MARKERS, POP_ID)
      } else { # AND operator between local and global maf
        vcf.maf <- maf.data %>%
          group_by(MARKERS, POP_ID) %>%
          summarise(
            MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
            MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
          ) %>%
          filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
          group_by(MARKERS) %>%
          tally() %>%
          filter(n >= maf.pop.num.threshold) %>%
          select(MARKERS) %>%
          left_join(input, by = "MARKERS") %>%
          arrange(MARKERS, POP_ID)
      }
    } # end maf snp approach
    
    
    message(stri_join("The number of MARKERS removed by the MAF filters = ", 
                      n_distinct(input$MARKERS)-n_distinct(vcf.maf$MARKERS), "\n", 
                      "The number of MARKERS before -> after the MAF filters: ", 
                      n_distinct(input$MARKERS)," -> ", n_distinct(vcf.maf$MARKERS), 
                      " MARKERS"))
    
    input <- vcf.maf
    
    # unused object
    vcf.maf <- NULL 
    maf.data <- NULL
  } # End of MAF filters
  
  # Write to working directory
  if (!is.null(filename)) {
    message(stri_paste("Writing the tidy data to the working directory: \n"), filename)
    write_tsv(x = input, path = filename, col_names = TRUE)
  }
  
  res <- input
  return(res)
} # tidy genomic data
