# tidy_stacks_haplotype_vcf
#' @title Make a tidy format of the batch_x.haplotypes.vcf file
#' @description Import and transform in tidy format the batch_x.haplotypes.vcf
#' file produced by STACKS.
#' @param data The stacks haplotype VCF file. e.g. \code{"batch_1.haplotypes.vcf"}
#' 
#' @param vcf.metadata (optional, logical) For the VCF file, 
#' with \code{vcf.metadata = FALSE}, 
#' only the genotype information is kept.
#' With default: \code{vcf.metadata = TRUE}, the metadata contained in the 
#' \code{FORMAT} field will be kept in the tidy data file.

#' @param strata A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. 
#' The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' If you have already run 
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data, 
#' the strata file is similar to a stacks `population map file`, make sure you 
#' have the required column names  (\code{INDIVIDUALS} and \code{STRATA}).

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

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

#' @param filename (optional) The file name for the tidy data frame
#' written to the working directory.
#' Default: \code{filename = NULL}, the tidy data is 
#' in the global environment only (i.e. not written in the working directory).

#' @examples
#' \dontrun{
#' The simplest form of the function
#' tidy.haplo.vcf <- tidy_stacks_haplotypes_vcf(
#' data = "batch_1.haplotypes.vcf",
#' strata = "strata.sturgeon.1pop.tsv"
#' )
#' Look carefully at the other arguments to customize your tidy haplotype vcf.
#' }

#' @rdname tidy_stacks_haplotypes_vcf
#' @export 
#' @import dplyr
#' @import stringi
#' @import tidyr
#' @importFrom data.table fread
#' @importFrom data.table melt.data.table
#' @importFrom data.table as.data.table


#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_stacks_haplotypes_vcf <- function(
  data, 
  vcf.metadata = TRUE, 
  strata, 
  pop.levels = NULL, 
  pop.labels = pop.levels, 
  pop.select = NULL, 
  whitelist.markers = NULL, 
  blacklist.id = NULL, 
  filename = NULL
  ) {
  

  # required to pass the R CMD check and have 'no visible binding for global variable'
  # if (getRversion() >= "2.15.1") {
  #   utils::globalVariables(
  #     c("DP", "AD", "vcf.headers", "GT_VCF", "INDIVIDUALS2", ""
  #     )
  #   )
  # }
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  if (is.null(strata)) stop("strata argument is required")
  
  
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
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
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
  }
  
  # Import VCF ****************************************************************
    message("Importing the Haplotypes VCF...")

    input <- data.table::fread(
      input = data,
      sep = "\t",
      stringsAsFactors = FALSE, 
      header = TRUE,
      # Automatically filter with blacklist of id
      drop = c("QUAL", "FILTER", "INFO", blacklist.id$INDIVIDUALS),
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
    message("Making the Haplotype VCF population wise")
    # input <- input %>%
    # tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, LOCUS, POS, REF, ALT)) # Gather individuals in 1 colummn
    
    input <- data.table::melt.data.table(
      data = as.data.table(input), 
      id.vars = c("CHROM", "LOCUS", "POS", "REF", "ALT"), 
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = "FORMAT_ID"
    ) %>% 
      as_data_frame()
    
    # population levels and strata  --------------------------------------------
    input <- left_join(x= input, y = strata.df, by = "INDIVIDUALS")
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # Tidy VCF
    message("Tidying the haplotype vcf...")
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
          )
      } # end cleaning AD column
    } else {
      input <- input %>% 
        select(CHROM, LOCUS, POS, REF, ALT, POP_ID, INDIVIDUALS, GT)
    }# end metadata section
    
    # recoding genotype and creating a new column combining CHROM, LOCUS and POS 
    input <- input %>%
      tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_", remove = FALSE) %>% 
      mutate(
        GT_VCF = GT,
        GT = stri_replace_all_fixed(str = GT, pattern = "./.", replacement = "-1/-1", vectorize_all = FALSE)
        ) %>% 
      tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = "/", remove = TRUE, extra = "drop") %>% 
      mutate(
        A1 = as.numeric(A1),
        A2 = as.numeric(A2),
        A1 = A1 + 1,
        A2 = A2 + 1,
        A1 = stri_pad_left(str= A1, pad = "0", width = 3),
        A2 = stri_pad_left(str= A2, pad = "0", width = 3)
        ) %>% 
      tidyr::unite(data = ., GT, A1, A2, sep = "")
      
    # Re ordering columns
    if (vcf.metadata) {
      # Re order columns
      common.colnames <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID", "INDIVIDUALS", "GT_VCF", "REF", "ALT")
      vcf.headers <- colnames(input)
      metadata.colnames <- purrr::discard(.x = colnames(input), .p = vcf.headers %in% common.colnames)
      input <- input[c(common.colnames, metadata.colnames)]
    } else {
      input <- input %>% select(MARKERS, CHROM, LOCUS, POS, POP_ID, INDIVIDUALS, GT_VCF, REF, ALT, GT)
    }

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

    # Write to working directory
    if (!is.null(filename)) {
      message(stri_paste("Writing the tidy data to the working directory: \n"), filename)
      write_tsv(x = input, path = filename, col_names = TRUE)
    }
    return(input)
}
