#' @title Read a VCF file produced by STACKS and transform in tidy format
#' @description Import a VCF file created by STACKS and mofify to a tidy format.
#' @param vcf.file The VCF file created by STACKS.
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An character string with your populations ordered.
#' @param pop.labels An optional character string with new populations names.
#' @param whitelist An optional filter of loci can be applied to the vcf, using
#' a file in the working directory (e.g. "myfile.txt") or an object
#' in the global environment.
#' @param blacklist (optional) Blacklist file with loci to filter out
#' of the vcf, using a file in the working directory (e.g. "myfile.txt")
#' or an object in the global environment.
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the directory
#'  (e.g. "blacklist.id.30.txt").
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
#' whitelist = "whitelist.loci.txt", 
#' blacklist = "blacklist.loci.paralogs.txt", 
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


read_stacks_vcf <- function(vcf.file, pop.id.start, pop.id.end, pop.levels, pop.labels, whitelist, blacklist, blacklist.id,
                            filename) {
#   setwd("/Users/thierry/Dropbox/partage_thierry_laura/gsi_sim")
#   # setwd("/Users/thierry/Dropbox/esturgeon_dropbox/stacks_populations_2015/01_stacks_populations/populations_8pop")
#   vcf.file <- "batch_1.vcf"
#   pop.id.start <- 1
#   pop.id.end <- 3
#   pop.levels <- c("TRI", "BON", "ANT", "GAS", "CAR", "MAL", "EDN", "CAP", "BRA", "SID", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO")
#   pop.labels <- c("TRI", "BON", "ANT", "GAS", "CAR", "MAL", "MAG", "MAG", "BRA", "DIN", "CAN", "LOB", "BRO", "OFF", "SEA", "BOO", "MAR", "BUZ", "RHO")
#   blacklist <- "blacklist.loci.paralogs.txt"
#   blacklist.id <- "blacklist.id.lobster.tsv"
#   filename <- "vcf.tidy.id.tsv"
#   
  
  X1 <- NULL
  POP_ID <- NULL
  QUAL <- NULL
  FILTER <- NULL
  FORMAT <- NULL
  ID <- NULL
  `#CHROM` <- NULL
  INFO <- NULL
  N <- NULL
  AF <- NULL
  INDIVIDUALS <- NULL
  REF <- NULL
  ALT <- NULL
  READ_DEPTH <- NULL
  REF_FREQ <- NULL
  ALT_FREQ <- NULL
  ALLELE_DEPTH <- NULL
  GT <- NULL
  GL <- NULL
  #   ALLELE_P <- NULL
  #   ALLELE_Q <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL
  
  # Import/read VCF ------------------------------------------------------------- 
  message("Importing the VCF...")
  
  vcf <- read_delim(
    vcf.file, 
    delim = "\t", 
    comment = "##",
    progress = interactive()
  ) %>% 
    select(-c(QUAL, FILTER)) %>%
    rename(LOCUS = ID, CHROM = `#CHROM`)
  
  # detect stacks version
  if(stri_detect_fixed(vcf$FORMAT[1], "AD")) {
    stacks.version <- "new"
  } else{
    stacks.version <- "old"
  }
  
  vcf <- vcf %>% select(-FORMAT)
  
    # Whitelist/Blacklist ---------------------------------------------------------- 
  if (missing(whitelist) == "TRUE") {
    vcf <- vcf
    message("No whitelist to apply to the VCF")
  } else if (is.vector(whitelist) == "TRUE") {
    vcf <- vcf %>% 
      semi_join(
        read_tsv(whitelist, col_names = T),
        by = "LOCUS"
      )
    message("Filtering the VCF with the whitelist from your directory")
  } else {
    vcf <- vcf %>%
      semi_join(whitelist, by = "LOCUS")
    message("Filtering the VCF with the whitelist from your global environment")
  }
  
  if (missing(blacklist) == "TRUE") {
    vcf <- vcf
    message("No blacklist to apply to the VCF")
  } else if (is.vector(blacklist) == "TRUE") {
    message("Filtering the VCF with the blacklist from your directory")
    vcf <- vcf  %>% 
      anti_join(
        read_tsv(blacklist, col_names = T),
        by = "LOCUS"
      )
  } else {
    message("Filtering the VCF with the blacklist from your global environment")
    vcf <- vcf %>% 
      anti_join(blacklist, by = "LOCUS")
  }
  
  # Make VCF tidy-----------------------------------------------------------------
    vcf <- vcf %>%
    tidyr::separate(INFO, c("N", "AF"), sep = ";", extra = "warn") %>%
    mutate(
      N = as.numeric(stri_replace_all_fixed(N, "NS=", "", vectorize_all=F)),
      AF = stri_replace_all_fixed(AF, "AF=", "", vectorize_all=F)
    ) %>%
    tidyr::separate(AF, c("REF_FREQ", "ALT_FREQ"), sep = ",", extra = "warn") %>%
    mutate(
      REF_FREQ = as.numeric(REF_FREQ),
      ALT_FREQ = as.numeric(ALT_FREQ)
    )
  # Gather individuals in 1 colummn --------------------------------------------
  vcf <- tidyr::gather(vcf, INDIVIDUALS, FORMAT, -c(CHROM:ALT_FREQ))
  message("Gathering individuals in 1 column")
  
  # Blacklist id----------------------------------------------------------------
  
  if (missing(blacklist.id) == "TRUE") {
    message("No individual blacklisted")
    blacklist.id <- NULL
    vcf <- vcf
  } else if (is.vector(blacklist.id) == "TRUE") {
    message("Using the blacklisted id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)
    vcf <- suppressWarnings(
      vcf %>% 
        anti_join(blacklist.id, by = "INDIVIDUALS")
    )
  } else {
    message("Using the blacklisted id from your global environment")
    blacklist.id <- blacklist.id
    vcf <- suppressWarnings(
      vcf %>% 
        anti_join(blacklist.id, by = "INDIVIDUALS")
    )
  }
  
  # Separate FORMAT and COVERAGE columns ---------------------------------------
  message("Tidying the VCF...")

  if(stacks.version == "new"){
    vcf <- vcf %>%
      tidyr::separate(FORMAT, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"),
                      sep = ":", extra = "warn")
  } else {
    # stacks version prior to v.1.29 had no Allele Depth field...
    message("Hum....")
    message("you are using an older version of STACKS...")
    message("It's not too late to use the last STACKS version, see STACKS change log for problems associated with older vcf files, for more details see: http://catchenlab.life.illinois.edu/stacks/")
    message("Continuing to work on tidying your VCF, for this time ;)")
    
    vcf <- vcf %>%
      tidyr::separate(FORMAT, c("GT", "READ_DEPTH", "GL"),
                      sep = ":", extra = "warn")
  }
  
  # Work with Mutate on CHROM and GL -------------------------------------------
  message("Fixing columns...")
  
  vcf <- vcf %>%
    mutate(
      CHROM = suppressWarnings(as.numeric(stri_replace_all_fixed(CHROM, "un", "1", vectorize_all=F))),
      GL = suppressWarnings(as.numeric(stri_replace_all_fixed(GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all=F)))
    ) %>%
    # Mutate read depth
    mutate(READ_DEPTH = suppressWarnings(as.numeric(stri_replace_all_regex(READ_DEPTH, "^0$", "NA", vectorize_all=F))))
  
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
  } else {
    # stacks version prior to v.1.29 had no Allele Depth field...
    vcf <- vcf 
  }  
  
  # Populations levels ---------------------------------------------------------
  message("Adding a population column, so that your tidy VCF is population wise...")
  
  if(missing(pop.labels)){
    pop.labels <- pop.levels
  } else {
    pop.labels <- pop.labels
  }
  
  vcf <- suppressWarnings(
    vcf %>% 
      mutate(
        POP_ID = factor(str_sub(INDIVIDUALS, pop.id.start, pop.id.end), 
                        levels = pop.levels, labels = pop.labels, ordered = T)
      ) %>%
      arrange(LOCUS, POS, POP_ID, INDIVIDUALS)
  )
  
  # Reorder the columns --------------------------------------------------------
  message("Reordering columns ...")
  if(stacks.version == "new"){
  # vcf <- vcf[c("CHROM", "LOCUS", "POS", "N", "REF", "ALT", "REF_FREQ", "ALT_FREQ", "POP_ID", "INDIVIDUALS", "GT", "ALLELE_P", "ALLELE_Q", "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "ALLELE_COVERAGE_RATIO", "GL")] # testing
  vcf <- vcf[c("CHROM", "LOCUS", "POS", "N", "REF", "ALT", "REF_FREQ", "ALT_FREQ", "POP_ID", "INDIVIDUALS", "GT", "READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH", "ALLELE_COVERAGE_RATIO", "GL")]
  } else {
    # stacks version prior to v.1.29 had no Allele Depth field...
    vcf <- vcf[c("CHROM", "LOCUS", "POS", "N", "REF", "ALT", "REF_FREQ", "ALT_FREQ", "POP_ID", "INDIVIDUALS", "GT", "READ_DEPTH", "GL")]
  }  
  
  
  # Save/Write the file to the working directory--------------------------------
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory, may take some time...")
    write_tsv(vcf, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected..."
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
