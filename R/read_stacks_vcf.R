#' @name read_stacks_vcf
#' @title Read a VCF file produced by STACKS and transform in tidy format
#' @description Import a VCF file created by STACKS and mofify to a tidy format.
#' @param vcf.file The VCF file created by STACKS.

#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by CHROM and/or locus and/or by snp.
#' The whitelist is in the working directory (e.g. "whitelist.txt").
#' de novo CHROM column with 'un' need to be changed to 1. 
#' Default \code{NULL} for no whitelist of markers.

#' @param blacklist.markers (optional) A blacklist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by CHROM and/or locus and/or by snp.
#' The blacklist is in the working directory (e.g. "blacklist.txt").
#' de novo CHROM column with 'un' need to be changed to 1. 
#' Default \code{NULL} for no blacklist of markers.

#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An character string with your populations ordered.
#' @param pop.labels An optional character string with new populations names.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").
#' @param pop.select (string) Conduct the assignment analysis on a
#' selected list of populations. Default = \code{NULL} for no selection and keep
#' all population.
#' e.g. \code{pop.select = "QUE"} to select QUE population samples.
#' \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population samples.

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
#' with 'un' need to be changed to 1. Default \code{NULL} for no blacklist of 
#' genotypes to erase.

#' @param filename (optional) The name of the file written in the directory.
#' @rdname read_stacks_vcf
#' @export
#' @import dplyr
#' @import readr
#' @details Your individuals are identified in this form : 
#' CHI-QUE-ADU-2014-020, (SPE-POP-MAT-YEA-ID). Then, \code{pop.id.start} = 5
#' and \code{pop.id.end} = 7. With 5 populations: ALK, JAP, NOR, QUE, LAB, 
#' the \code{pop.levels} option could look like this: 
#' c("JAP", "ALK", "QUE", "LAB", "NOR"). Useful for hierarchical clustering.
#' Allele depth was introduced in STACKS v.1.28, so this field won't 
#' be available to prior version.
#' @examples
#' \dontrun{
#' vcf.tidy <- read_stacks_vcf(
#' vcf.file = "batch_1.vcf", 
#' pop.id.start = 5, 
#' pop.id.end = 7, 
#' pop.levels = c("JAP", "ALK", "QUE", "LAB", "NOR"), 
#' whitelist.markers = "whitelist.loci.txt", 
#' blacklist.markers = "blacklist.paralogs.txt", 
#' filename = "vcf.tidy.paralogs.tsv")
#' }
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.
#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1"){
  utils::globalVariables(
    c('QUAL', 'FILTER', 'INFO', 'ID', '#CHROM', 'FORMAT', 'FORMAT_ID', 'REF', 
      'ALT', 'ERASE', 'GT', 'ALLELE_DEPTH', 'GL', 'READ_DEPTH', 
      'ALLELE_REF_DEPTH', 'ALLELE_ALT_DEPTH')
  )
}
read_stacks_vcf <- function(vcf.file, 
                            whitelist.markers, 
                            blacklist.markers,
                            blacklist.genotype,
                            pop.id.start, 
                            pop.id.end, 
                            pop.levels, 
                            pop.labels, 
                            blacklist.id,
                            pop.select,
                            filename) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(vcf.file)) stop("VCF file required")
  if (missing(whitelist.markers)) whitelist.markers <- NULL # no Whitelist
  if (missing(blacklist.markers)) blacklist.markers <- NULL # no Blacklist
  if (missing(blacklist.genotype)) blacklist.genotype <- NULL # no genotype to erase
  if (missing(blacklist.id)) blacklist.id <- NULL # No blacklist of ID
  if (missing(pop.levels)) stop("pop.levels required")
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels
  if (missing(pop.id.start)) stop("pop.id.start required")
  if (missing(pop.id.end)) stop("pop.id.end required")
  if (missing(pop.select)) pop.select <- NULL
  if (missing(filename)) filename <- NULL
  
  # Import/read VCF ************************************************************
  message("Importing the VCF...")
  vcf <- read_delim(
    vcf.file,
    delim = "\t",
    comment = "##",
    progress = interactive()
  ) %>%
    select(-c(QUAL, FILTER, INFO)) %>%
    rename(LOCUS = ID, CHROM = `#CHROM`) %>%
    mutate(
      CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
    )
  
  # Detect STACKS version ******************************************************
  if(stri_detect_fixed(vcf$FORMAT[1], "AD")) {
    stacks.version <- "new"
  } else{
    stacks.version <- "old"
  }
  vcf <- vcf %>% select(-FORMAT)
  
  # Whitelist of markers *******************************************************
  if (is.null(whitelist.markers)) { # no Whitelist
    message("No whitelist of markers to apply to the VCF")
    vcf <- vcf
  } else { # with Whitelist of markers
    message("Filtering the VCF with the whitelist of markers from your directory")
    whitelist.markers <- read_tsv(whitelist.markers, col_names = TRUE)
    columns.names.whitelist <- colnames(whitelist.markers)
    if ("CHROM" %in% columns.names.whitelist){
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
    vcf <- suppressWarnings(
      vcf %>%
        semi_join(whitelist.markers, by = columns.names.whitelist)
    )
  }
  
  # Blacklist of markers *******************************************************
  if (is.null(blacklist.markers)) { # no Blacklist
    message("No blacklist of markers to apply to the VCF")
    vcf <- vcf
  } else { # with Blacklist of markers
    message("Filtering the VCF with the blacklist of markers from your directory")
    blacklist.markers <- read_tsv(blacklist.markers, col_names = TRUE)
    columns.names.blacklist <- colnames(blacklist.markers)
    if ("CHROM" %in% columns.names.blacklist){
      blacklist.markers$CHROM <- as.character(blacklist.markers$CHROM)
    }
    vcf <- suppressWarnings(
      vcf %>%
        anti_join(blacklist.markers, by = columns.names.blacklist)
    )
  }
  
  # Tidying the VCF to make it easy to work on the data for conversion *********
  message("Making the VCF population wise")
  vcf <- vcf %>% 
    tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, LOCUS, POS, REF, ALT)) %>% # Gather individuals in 1 colummn
    mutate( # Make population ready
      POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
      POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered =TRUE),
      INDIVIDUALS =  as.character(INDIVIDUALS)
    )
  
  # Blacklist id ***************************************************************
  if (is.null(blacklist.id)) { # No blacklist of ID
    message("No individual blacklisted")
    vcf <- vcf
  } else { # With blacklist of ID
    message("Using the blacklisted id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = TRUE)
    vcf <- suppressWarnings(
      vcf %>%
        anti_join(blacklist.id, by = "INDIVIDUALS") %>%
        mutate(POP_ID = droplevels(POP_ID))
    )
  }

  # Pop select *****************************************************************
  if (is.null(pop.select)){
    vcf <- vcf
  } else {
    message(stri_join(length(pop.select), "population(s) selected", sep = " "))
    vcf <- suppressWarnings(
      vcf %>%
        filter(POP_ID %in% pop.select)
    )
  }
  
  # Blacklist genotypes ********************************************************
  if (is.null(blacklist.genotype)) { # no Whitelist
    message("No genotype to erase")
    vcf <- vcf
  } else {
    message("Erasing genotype with the blacklist")
    blacklist.genotype <- read_tsv(blacklist.genotype, col_names = TRUE)
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    if ("CHROM" %in% columns.names.blacklist.genotype){
      columns.names.blacklist.genotype$CHROM <- as.character(columns.names.blacklist.genotype$CHROM)
    }
    
    # control check to keep only whitelisted markers from the blacklist of genotypes
    if (!is.null(whitelist.markers)){
      message("Control check to keep only whitelisted markers 
              present in the blacklist of genotypes to erase.")
      # updating the whitelist of markers to have all columns that id markers
      whitelist.markers.ind <- vcf %>% select(CHROM, LOCUS, POS, INDIVIDUALS) %>% distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      # updating the blacklist.genotype
      blacklist.genotype <- suppressWarnings(semi_join(whitelist.markers.ind, blacklist.genotype, by = columns.names.blacklist.genotype))
    } else {
      blacklist.genotype <- blacklist.genotype
    }
    
    # control check to remove blacklisted individuals from the blacklist of genotypes
    if (!is.null(blacklist.id)){
      message("Control check to remove blacklisted individuals 
              present in the blacklist of genotypes to erase.")
      blacklist.genotype <- suppressWarnings(anti_join(blacklist.genotype, blacklist.id, by = "INDIVIDUALS"))
    } else {
      blacklist.genotype <- blacklist.genotype
    }
    
    # Add one column that will allow to include the blacklist in the dataset 
    # by x column(s) of markers
    blacklist.genotype <- mutate(.data = blacklist.genotype, ERASE = rep("erase", n()))
    
    vcf <- suppressWarnings(
      vcf %>%
        full_join(blacklist.genotype, by = c("CHROM", "LOCUS", "POS", "INDIVIDUALS")) %>%
        mutate(
          ERASE = stri_replace_na(str = ERASE, replacement = "ok"),
          GT = ifelse(ERASE == "erase", "./.", GT)
        ) %>% 
        select(-ERASE)
    )
  } # end erase genotypes
  
  # Separate FORMAT and COVERAGE columns ---------------------------------------
  message("Tidying the VCF...")
  
  if(stacks.version == "new"){
    vcf <- suppressWarnings(
      vcf %>%
        tidyr::separate(FORMAT, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"),
                        sep = ":", extra = "warn") %>% 
        tidyr::separate(ALLELE_DEPTH, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
                        sep = ",", extra = "warn")
    )
  } else {
    # stacks version prior to v.1.29 had no Allele Depth field...
    message("Hum....")
    message("you are using an older version of STACKS...")
    message("It's not too late to use the last STACKS version, see STACKS change log for problems associated with older vcf files, for more details see: http://catchenlab.life.illinois.edu/stacks/")
    message("Continuing to work on tidying your VCF, for this time ;)")
    
    vcf <- vcf %>%
      tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "GL"),
                      sep = ":", extra = "warn")
  }
  
  # Work with Mutate on CHROM and GL -------------------------------------------
  message("Fixing columns...")
  vcf <- vcf %>%
    mutate(
      GL = suppressWarnings(as.numeric(stri_replace_all_fixed(GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all=F)))
    ) %>%
    # Mutate read depth
    mutate(
      READ_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(READ_DEPTH, "^0$", "NA", vectorize_all=F)))
      )
  
  # mutate the alleles REF/ALT depth
  if(stacks.version == "new"){
    vcf <- vcf %>%
      mutate(
        ALLELE_REF_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(ALLELE_REF_DEPTH, "^0$", "NA", vectorize_all = TRUE))),
        ALLELE_ALT_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(ALLELE_ALT_DEPTH, "^0$", "NA", vectorize_all = TRUE))),
        # Mutate coverage ratio for allelic imbalance
        ALLELE_COVERAGE_RATIO = suppressWarnings(
          as.numeric(ifelse(GT == "./." | GT == "0/0" | GT == "1/1", "NA",
                            ((ALLELE_ALT_DEPTH - ALLELE_REF_DEPTH)/(ALLELE_ALT_DEPTH + ALLELE_REF_DEPTH)))))
      )
  } else {# stacks version prior to v.1.29 had no Allele Depth field...
    vcf <- vcf 
  }

  # Reorder the columns --------------------------------------------------------
  message("Reordering columns ...")
  if(stacks.version == "new"){
    vcf <- vcf[c("CHROM", "LOCUS", "POS", "REF", "ALT", "POP_ID", "INDIVIDUALS", "GT", "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "ALLELE_COVERAGE_RATIO", "GL")]
  } else {# stacks version prior to v.1.29 had no Allele Depth field...
    vcf <- vcf[c("CHROM", "LOCUS", "POS", "REF", "ALT", "POP_ID", "INDIVIDUALS", "GT", "READ_DEPTH", "GL")]
  }  
  
  # Save/Write the file to the working directory--------------------------------
  if (is.null(filename)) {
    saving <- "Saving was not selected..."
  } else {
    message("Saving the file in your working directory, may take some time...")
    write_tsv(vcf, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  }
  
  # Message at the end ---------------------------------------------------------- 
  invisible(cat(sprintf(
    "%s\n
Working directory:
%s",
    saving, getwd()
  )))
  vcf
}
