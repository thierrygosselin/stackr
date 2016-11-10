# Import and prune data file

#' @name data_pruning
#' @title Data set pruning
#' @description This function enable to prune the data set with different 
#' argument (see below) in order to be prep for subsequent analysis.
#' Various input files are offered. Individuals, populations and
#' markers can be filtered and/or selected in several ways using blacklist,
#' whitelist and other arguments.

#' @param data 6 options: vcf (to make vcf population ready, see details below),
#' plink, stacks haplotype file, genind, genepop, 
#' and a data frame in wide format. \emph{See details}.

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
#'and \code{ONT} population samples (out of 20 pops).
# Default: \code{pop.select = NULL} 


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
#' for "MARKERS" in column names (TRUE = long format).
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
#' \code{MARKERS} and \code{GENOTYPE or GT}. The rest are considered metata info.
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




#' @export
#' @rdname data_pruning
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread
#' @importFrom utils write.table packageVersion

#' @examples
#' \dontrun{
#' data_pruning(data = "plink.tped",
#' pop.levels = c("pop_1", "pop_2", "pop_3", "pop_4", "pop_5", "pop_6", "pop_7"),
#' monomorphic.out <- TRUE,
#' maf.thresholds <- c(0.05, 0.1),
#' maf.pop.num.threshold <- 1,
#' maf.approach <- "SNP",
#' maf.operator <- "OR",
#' common.markers <- TRUE,
#' pop.select <- c("pop_1", "pop_2", "pop_3")
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
#' American Journal of Human Genetics. 2007; 81: 559–575. doi:10.1086/519795
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

data_pruning <- function(data,
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
  
  
  # Create a filename to save the output files ********************************
  # Get date and time to have unique filenaming
  if (is.null(filename)) {
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    if (data.type == "vcf.file") {
      filename <- stri_replace_all_fixed(data, pattern = ".vcf", replacement = (stri_paste("_pruning_", file.date, ".vcf")), vectorize_all = FALSE)
    }
    
    if (data.type == "haplo.file") {
      filename <- stri_replace_all_fixed(data, pattern = ".tsv", replacement = (stri_paste("_pruning_", file.date, ".tsv")), vectorize_all = FALSE)
    }
    
    if (data.type == "df.file") {
      filename <- stri_paste(data, "_pruning_", file.date, ".tsv")
    }
    
    if (data.type == "plink.file") {
      filename.tped <- stri_replace_all_fixed(data, pattern = ".tped", replacement = (stri_paste("_pruning_", file.date, ".tped")), vectorize_all = FALSE)
      filename.tmap <- stri_replace_all_fixed(data, pattern = ".tped", replacement = (stri_paste("_pruning_", file.date, ".tmap")), vectorize_all = FALSE)
    }
    
  }
  
  # Import file ----------------------------------------------------------------
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
    strata = strata, 
    pop.levels = pop.levels, 
    pop.labels = pop.labels, 
    pop.select = pop.select,
    filename = NULL
  )
  
  # create a strata.df
  strata.df <- NULL
  strata.df <- input %>% 
    select(INDIVIDUALS, POP_ID) %>% 
    distinct(INDIVIDUALS, .keep_all = TRUE)
  
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  # back to original format ****************************************************
  if (data.type == "vcf.file") {
    if ("POP_ID" %in% colnames(input)){
      input <- select(.data = input, -POP_ID)
    }
    
    
    input <- input %>% 
      select(-GT) %>% 
      rename(GT = GT_VCF)
    
    info.field <- suppressWarnings(
      input %>% 
        group_by(MARKERS) %>%
        filter(GT != "./.") %>% 
        tally %>% 
        mutate(INFO = stri_paste("NS=", n, sep = "")) %>% 
        select(-n)
    )
    
    output <- suppressWarnings(
      left_join(input, info.field, by = "MARKERS") %>% 
        group_by(MARKERS, INFO) %>% 
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT) %>%
        ungroup() %>% 
        # tidyr::separate(MARKERS, c("CHROM", "ID", "POS"), sep = "_", extra = "warn") %>%
        mutate(
          LOCUS = as.numeric(LOCUS),
          POS = as.numeric(POS),
          QUAL = rep(".", n()),
          FILTER = rep("PASS", n()),
          FORMAT = rep("GT", n())
        ) %>% 
        arrange(CHROM, LOCUS, POS) %>% 
        ungroup() %>%
        select(-MARKERS) %>% 
        select('#CHROM' = CHROM, POS, ID = LOCUS, REF, ALT, QUAL, FILTER, INFO, FORMAT, everything())
    )
    
    # File format
    file.format <- "##fileformat=VCFv4.3"
    file.format <- as.data.frame(file.format)
    utils::write.table(x = file.format, file = filename, sep = " ", append = FALSE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # File date
    file.date <- stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")
    file.date <- stri_paste("##fileDate=", file.date, sep = "")
    utils::write.table(x = file.date, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Source
    file.source <- as.data.frame(stri_paste("##source=stackr v.", utils::packageVersion("stackr"), sep = ""))
    utils::write.table(x = file.source, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Info field 1
    info1 <- '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">'
    info1 <- as.data.frame(info1)
    utils::write.table(x = info1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Format field 1
    format1 <- '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    format1 <- as.data.frame(format1)
    utils::write.table(x = format1, file = filename, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # Write the prunned vcf to the file
    suppressWarnings(
      write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE)
    )
  }
  if (data.type == "haplo.file") {
    cnt.field <- suppressWarnings(
      input %>% 
        group_by(MARKERS) %>%
        filter(GT != "-") %>% 
        tally %>%
        rename(Cnt = n) %>% 
        arrange(as.integer(MARKERS))
    )
    
    output <- suppressWarnings(
      left_join(input, cnt.field, by = "MARKERS") %>% 
        group_by(MARKERS) %>% 
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT) %>%
        ungroup() %>% 
        arrange(as.integer(MARKERS)) %>% 
        rename(`Catalog ID` = MARKERS)
    )
    
    # Write the prunned vcf to the file
    suppressWarnings(
      write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE)
    )
  } # end haplotypes file
  
  if (data.type == "plink.file") {
    # to create a PLINK tped and tfam
    tped <- input %>% 
      mutate(
        COL1 = rep("0", n()),
        COL3 = rep("0", n()),
        COL4 = rep("0", n())
      ) %>% 
      arrange(POP_ID, INDIVIDUALS, INDIVIDUALS_ALLELES) %>% 
      select(COL1, MARKERS, COL3, COL4, INDIVIDUALS_ALLELES, GT) %>%
      group_by(COL1, MARKERS, COL3, COL4) %>% 
      tidyr::spread(data = ., key = INDIVIDUALS_ALLELES, value = GT) %>% 
      arrange(MARKERS)
    
    write_delim(x = tped, path = filename.tped, col_names = FALSE, delim = " ")
    
    # Create a tfam file
    tfam <- strata.df %>% 
      select(POP_ID, INDIVIDUALS) %>% 
      arrange(POP_ID, INDIVIDUALS) %>% 
      mutate(
        COL3 = rep("0",n()),
        COL4 = rep("0",n()),
        COL5 = rep("0",n()),
        COL6 = rep("-9",n())
      )
    
    # the tfam must have the same name as the tped 
    write_delim(x = tfam, path = filename.tmap, col_names = FALSE, delim = " ")
  } # end plink
  if (data.type == "df.file") {
    output <- input %>% 
      select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
      arrange(MARKERS, POP_ID, INDIVIDUALS) %>% 
      group_by(POP_ID, INDIVIDUALS) %>% 
      tidyr::spread(data = ., key = MARKERS, value = GT)
    
    # Write the prunned vcf to the file
    suppressWarnings(
      write_tsv(x = output, path = filename, append = TRUE, col_names = TRUE)
    )
  }# end df
  
} # End data_pruning



