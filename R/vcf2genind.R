# Write a adegenet genind object from STACKS VCF file

#' @name vcf2genind
#' @title Create a \code{adegenet} \code{\link[adegenet]{genind}} object from a \code{STACKS} vcf file
#' @description This function can first filter the vcf file 
#' with a whitelist of loci
#' and a blacklist of individuals (optional). Then it will convert the file
#' to a \code{adegenet} \code{\link[adegenet]{genind}} object.
#' Map-independent imputation using Random Forest or the most frequent category
#' is also available as an option.


#' @param data (biallelic data) 6 options: vcf, plink, genind, genepop, 
#' and a data frame in wide or long/tidy format. \emph{See details}.
#' The function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}} and 
#' \code{\link[stackr]{tidy_genomic_data}}.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

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
#' files. Use stackr for haplotypes files and use the whitelist of markers.
#' Default: \code{maf.thresholds = NULL}. 

#' @param maf.pop.num.threshold (integer, optional) When maf thresholds are used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. Default: \code{maf.pop.num.threshold = 1}

#' @param maf.approach (character, optional). By \code{maf.approach = "SNP"} or 
#' by \code{maf.approach = "haplotype"}.
#' The function will consider the SNP or ID/LOCUS/haplotype/read MAF statistics 
#' to filter the markers.
#' Default is \code{maf.approach = "SNP"}. The \code{haplotype} approach is 
#' restricted to VCF file.

#' @param maf.operator (character, optional) \code{maf.operator = "AND"} or 
#' default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?

#' @param max.marker An optional integer useful to subsample marker number in 
#' large PLINK file. e.g. if the data set 
#' contains 200 000 markers and \code{max.marker = 10000} 10000 markers are
#' subsampled randomly from the 200000 markers. Use \code{whitelist.markers} to
#' keep specific markers.
#' Default: \code{max.marker = NULL}.

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
#' Default: \code{pop.labels = NULL}. When pop.levels is not null and pop.labels
#' is not specified. pop.labels = pop.levels.
#' If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.

#' @param strata (optional for data frame and PLINK files, 
#' required for VCF and haplotypes files) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' Default: \code{strata = NULL}.

#' @param pop.select (string, optional) Selected list of populations for 
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#'and \code{ONT} population samples (out of 20 pops).
# Default: \code{pop.select = NULL} 

#' @param hierarchy (optional) A formula that explicitely defines hierarchical levels 
#' in your strata. See \code{\link[adegenet]{genind}} for details.

#' @param imputation.method (character, optional) 
#' Methods available for map-independent imputations of missing genotype: 
#' (1) \code{"max"} to use the most frequent category for imputations.
#' (2) \code{"rf"} using Random Forest algorithm. 
#' Default: no imputation \code{imputation.method = NULL}.

#' @param impute (character, optional) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
#' Default: \code{"genotype"}.

#' @param imputations.group (character, optional) \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default = \code{"populations"}.

#' @param num.tree (integer, optional) The number of trees to grow in Random Forest. 
#' Default: \code{num.tree = 100}.

#' @param iteration.rf (integer, optional) The number of iterations of missing data algorithm
#' in Random Forest. 
#' Default: \code{iteration.rf = 10}.

#' @param split.number (integer, optional) Non-negative integer value used to specify
#' random splitting in Random Forest. 
#' Default: \code{split.number = 100}.

#' @param verbose (logical, optional) Should trace output be enabled on each iteration
#' in Random Forest ? 
#' Default: \code{verbose = FALSE}.

#' @param parallel.core (optional) The number of core for OpenMP shared-memory parallel
#' programming of Random Forest imputations. For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.


#' @details 
#' \strong{Input files:}
#' \enumerate{
#' \item VCF file (e.g. \code{data = "batch_1.vcf"}). 
#' To make the VCF population ready, you need the \code{strata} argument.
#' 
#' \item haplotype file created in STACKS (e.g. \code{data = "batch_1.haplotypes.tsv"}).
#' To make the haplotype file population ready, you need the \code{strata} argument.
#' 
#' \item Data frame
#' Tab delimitted.
#' \strong{2 genotypes formats are available, both use 3 character per allele:}
#' 6 characters no allele separator: e.g. \code{001002 of 111333} (for heterozygote individual).
#' 6 characters WITH an allele separator: e.g. \code{001/002 of 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#' Missing alleles are coded \code{000}.
#' To discriminate the long from the wide format, 
#' the function \pkg{stackr} \code{\link[stackr]{read_long_tidy_wide}} searches 
#' for columns number, > 20 for wide 
#' (i.e. don't use less than 10 markers in wide format, the function was not designed for that).
#' 
#' Data Frame wide format:
#' The wide format cannot store metadata info.
#' The wide format contains starts with these 2 id columns: 
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals), 
#' the remaining columns are the markers in separate columns storing genotypes.
#' This format requires column numbers to be larger than 20.

#' Data frame long/tidy format:
#' This format requires column numbers to be within the range: 4 min - 20 max.
#' The long format is considered to be a tidy data frame and can store metadata info. 
#' (e.g. from a VCF see \pkg{stackr} \code{\link[stackr]{tidy_genomic_data}}). The 4 columns
#' required in the long format are: \code{INDIVIDUALS}, \code{POP_ID}, 
#' \code{MARKERS} and \code{GENOTYPE or GT}.
#' 
#' Note that the \code{POP_ID} column can be any hierarchical grouping. 
#' See the argument \code{strata} for other means of controlling grouping used 
#' in the assignment.
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
#' \item genepop data file (e.g. \code{data = kiwi_data.gen}). Here, the function can only use
#' alleles encoded with 3 digits.
#' }
#' 
#' \strong{Imputations details:}
#' The imputations using Random Forest requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set 
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals 
#' will require 15 min.

#' @return When no imputation is selected an object of the 
#' class \code{\link[adegenet]{genind}} is returned.
#' When imputation is selected a list with 2 objects is returned
#' and accessed with \code{$no.imputation} or \code{$imputed}.

#' @export
#' @rdname vcf2genind
# @importFrom adegenet df2genind
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread
#' @import parallel
#' @import dplyr
#' @import stringi

#' @examples
#' \dontrun{
#' snowcrab <- vcf2genind(
#' data = "batch_1.vcf",
#' strata = "snowcrab.strata.tsv",
#' whitelist.markers = "whitelist.txt",
#' common.markers = TRUE,
#' blacklist.id = "blacklist.id.lobster.tsv",
#' pop.levels = c("PAN", "COS")
#' imputation.method = "rf",
#' impute = "genotype",
#' imputations.group <- "populations", 
#' num.tree <- 100,
#' iteration.rf <- 10,
#' split.number <- 100,
#' verbose <- FALSE,
#' parallel.core = 12
#' )
#' 
#' A list with 2 genind objects, with and without imputation: 
#' no.imputation <- snowcrab$no.imputation
#' imputed <- snowcrab$imputed
#' }

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1: 
#' new tools for the analysis of genome-wide SNP data. 
#' Bioinformatics, 27, 3070-3071.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.
#' @seealso \code{adegenet} is available on CRAN \url{http://cran.r-project.org/web/packages/adegenet/} and github \url{https://github.com/thibautjombart/}
#' \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# to get rid of notes in build check
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", 
                           "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", 
                           "SAMPLES", "ALLELES", 'A1', 'A2', 'COUNT', 
                           "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", 
                           "POLYMORPHISM", "POLYMORPHISM_MAX", "other", 
                           "strata", "hierarchy", "GROUP", ".", 'MARKERS', 
                           'MARKERS_ALLELES', 'STRATA'
  )
  )
}


vcf2genind <- function(
  data,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  max.marker = NULL,
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  imputation.method = NULL,
  impute = "genotype",
  imputations.group = "populations",
  num.tree = 100,
  iteration.rf = 10,
  split.number = 100,
  verbose = FALSE,
  parallel.core = detectCores()-1
) {
  
  cat("#######################################################################\n")
  cat("######################### stackr::vcf2genind ##########################\n")
  cat("#######################################################################\n")
  # Checking for missing and/or default arguments ##############################
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  if (is.null(imputation.method)) {
    message("vcf2genind: no imputation...")
  } else {
    message("vcf2genind: with imputations...")
  }
  
  # File type detection#########################################################
  if(adegenet::is.genind(data)){
    data.type <- "genind.file"
    message("File type: genind object")
  } else {
    data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
    if (identical(data.type, "##fileformat=VCF") | stri_detect_fixed(str = data, pattern = ".vcf")) {
      data.type <- "vcf.file"
      # message("File type: VCF")
    }
    if (stri_detect_fixed(str = data, pattern = ".tped")) {
      data.type <- "plink.file"
      # message("File type: PLINK")
      if (!file.exists(stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        stop("Missing tfam file with the same prefix as your tped")
      }
    } 
    if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS")) {
      data.type <- "df.file"
      # message("File type: data frame of genotypes")
    }
    if (stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      # data.type <- "haplo.file"
      message("File type: haplotypes from stacks")
      if (is.null(blacklist.genotype)) {
        stop("blacklist.genotype file missing. 
             Use stackr's missing_genotypes function to create this blacklist")
      }
    }
    if (stri_detect_fixed(str = data, pattern = ".gen")) {
      # data.type <- "genepop.file"
      message("File type: genepop")
      message("Multilocus genepop file won't work, only for biallelic markers")
    } 
    
  } # end file type detection
  
  # Strata argument required for VCF and haplotypes files#######################
  if (data.type == "vcf.file" & is.null(strata)) stop("strata argument is required")
  if (data.type == "haplo.file") stop("This function is for biallelic dataset only")
  
  # Import #####################################################################
  input <- stackr::tidy_genomic_data(
    data = data, 
    vcf.metadata = FALSE,
    blacklist.id = blacklist.id, 
    blacklist.genotype = blacklist.genotype, 
    whitelist.markers = whitelist.markers, 
    monomorphic.out = monomorphic.out, 
    max.marker = max.marker,
    snp.ld = snp.ld, 
    common.markers = common.markers, 
    maf.thresholds = maf.thresholds, 
    maf.pop.num.threshold = maf.pop.num.threshold, 
    maf.approach = maf.approach, 
    maf.operator = maf.operator,
    strata = strata, 
    pop.levels = pop.levels, 
    pop.labels = pop.labels, 
    pop.select = pop.select,
    filename = NULL
  )
  
  # create a strata.df
  strata.df <- input %>% 
    select(INDIVIDUALS, POP_ID) %>% 
    distinct(INDIVIDUALS)
  strata <- strata.df
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  # Preparing adegenet genind object ---------------------------------------------
  message("Preparing adegenet genind object")
  genind.prep <- suppressWarnings(
    input %>%
      tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
      tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID))  %>%
      filter(GT != "000") %>% 
      tidyr::spread(data = ., key = MARKERS, value = GT) %>% # this reintroduce the missing, but with NA
      ungroup() %>% 
      plyr::colwise(.fun = factor, exclude = NA)(.)
  )
  
  # The next part is longer than it used to be with VCF file only, 
  # but it as the advantage of working and simplifying the use for other file type.
  genind.prep <- suppressWarnings(
    genind.prep %>%
      ungroup() %>% 
      mutate_each(funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
      ungroup()
  )
    
    genind.prep <- tidyr::gather(data = genind.prep, key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
      mutate(GT = stri_replace_na(str = GT, replacement = "000")) %>%
      filter(GT != "000") %>%
      select(-ALLELES) %>%
      group_by(POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
      tally %>%
      ungroup()
      
    genind.prep <- genind.prep %>%
      tidyr::unite(MARKERS_ALLELES, MARKERS, GT, sep = ":", remove = TRUE) %>%
      arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
      group_by(POP_ID, INDIVIDUALS) %>% 
      tidyr::spread(data = ., key = MARKERS_ALLELES, value = n) %>%
      ungroup()
    
    genind.prep <- tidyr::gather(data = genind.prep, key = MARKERS_ALLELES, value = COUNT, -c(INDIVIDUALS, POP_ID)) %>% 
      tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = ":", remove = TRUE) %>% 
      mutate(COUNT = as.numeric(stri_replace_na(str = COUNT, replacement = "0"))) %>% 
      group_by(INDIVIDUALS, MARKERS) %>%
      mutate(MAX_COUNT_MARKERS = max(COUNT, na.rm = TRUE)) %>%
      ungroup() %>% 
      mutate(COUNT = ifelse(MAX_COUNT_MARKERS == 0, "erase", COUNT)) %>%
      select(-MAX_COUNT_MARKERS) %>% 
      mutate(COUNT = replace(COUNT, which(COUNT == "erase"), NA)) %>% 
      arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES)
    
    genind.prep <- genind.prep %>%  
      tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
      tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>%
      mutate(
        INDIVIDUALS = as.character(INDIVIDUALS),
        POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
        POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
      ) %>% 
      arrange(POP_ID, INDIVIDUALS)

  
  # genind construction: no imputation -----------------------------------------
  # genind arguments common to all data.type
  message("Building the genind object")
  ind <- genind.prep$INDIVIDUALS
  pop <- genind.prep$POP_ID
  genind.df <- genind.prep %>% ungroup() %>% 
    select(-c(INDIVIDUALS, POP_ID))
  rownames(genind.df) <- ind
  loc.names <- colnames(genind.df)
  strata <- genind.prep %>% ungroup() %>% select(INDIVIDUALS, POP_ID) %>% distinct(INDIVIDUALS, POP_ID)
  
  # genind constructor
  prevcall <- match.call()
  no.imputation <- genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = NULL)
  
  # sum <- summary(genind.object) # test
  # sum$NA.perc # test
  
  ind <- NULL
  pop <- NULL
  genind.df <- NULL
  loc.names <- NULL
  strata <- NULL
  prevcall <- NULL
  genind.prep <- NULL
  
  if (is.null(imputation.method)) {
    message("A large 'genind' object (no imputation) was created in your Environment")
  } else if (imputation.method == "max"){
    message("Calculating map-independent imputations using the most frequent allele.")
  } else {
    message("Calculating map-independent imputations using random forest")
  }
  
  # Imputations-----------------------------------------------------------------
  if (!is.null(imputation.method)) {
    
    input.imp <- stackr::stackr_imputations_module(
      data = input, 
      imputation.method = imputation.method, 
      impute = impute, 
      imputations.group = imputations.group, 
      num.tree = num.tree, 
      iteration.rf = iteration.rf, 
      split.number = split.number, 
      verbose = verbose, 
      parallel.core = parallel.core, 
      filename = "dataset.imputed.tsv"
    )
    
    message("Preparing imputed data for adegenet genind object")
    genind.prep.imp <- input.imp %>%
      tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3) %>%  # separate the genotypes into alleles
      tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID))
    
    genind.prep.imp <- suppressWarnings(
      genind.prep.imp %>%
        tidyr::spread(data = ., key = MARKERS, value = GT) %>% # this reintroduce the missing, but with NA
        ungroup() %>% 
        plyr::colwise(.fun = factor, exclude = NA)(.)
    )
    
    genind.prep.imp <- suppressWarnings(
      genind.prep.imp %>%
        ungroup() %>% 
        mutate_each(funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
        ungroup() %>% 
        tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
        select(-ALLELES) %>%
        group_by(POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
        tally %>%
        ungroup() %>%
        tidyr::unite(MARKERS_ALLELES, MARKERS, GT, sep = ":", remove = TRUE) %>%
        arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
        group_by(POP_ID, INDIVIDUALS) %>% 
        tidyr::spread(data = ., key = MARKERS_ALLELES, value = n) %>%
        ungroup() %>%
        tidyr::gather(data = ., key = MARKERS_ALLELES, value = COUNT, -c(INDIVIDUALS, POP_ID)) %>% 
        tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = ":", remove = TRUE) %>% 
        mutate(COUNT = as.numeric(stri_replace_na(str = COUNT, replacement = "0"))) %>% 
        ungroup() %>%
        arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES) %>% 
        tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
        tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>% 
        ungroup () %>%
        mutate(
          INDIVIDUALS = as.character(INDIVIDUALS),
          POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
          POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
        ) %>% 
        arrange(POP_ID, INDIVIDUALS)
    )
    
    # genind construction: no imputation ---------------------------------------
    # 1) the genind without imputations is put in a res list
    res <- list()
    res$no.imputation <- no.imputation
    no.imputation <- NULL # drop unused object
    
    # 2) the genind with imputations
    message("Building the genind object")
    ind <- genind.prep.imp$INDIVIDUALS
    pop <- genind.prep.imp$POP_ID
    genind.df <- genind.prep.imp %>%
      ungroup() %>% 
      select(-c(INDIVIDUALS, POP_ID))
    rownames(genind.df) <- ind
    loc.names <- colnames(genind.df)
    
    strata <- genind.prep.imp %>% 
      ungroup() %>% 
      select(INDIVIDUALS, POP_ID) %>% 
      distinct(INDIVIDUALS, POP_ID)
    
    # genind constructor
    prevcall <- match.call()
    res$imputed  <- adegenet::genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = hierarchy)
    message("A large 'genind' object was created in your Environment (with and without imputations)")
    # sum <- summary(res$imputed) # test
    # sum$NA.perc # test
    
    ind <- NULL
    pop <- NULL
    genind.df <- NULL
    loc.names <- NULL
    strata <- NULL
    prevcall <- NULL
    # genind.prep <- NULL
    # genind.prep.imp <- NULL
    
  } # End imputations

  # outout results -------------------------------------------------------------
  return(res)
}
