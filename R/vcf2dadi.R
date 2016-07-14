#' @name vcf2dadi
#' @title Create a \code{dadi} SNP input file from a any vcf file.
#' @description This function will create a \code{dadi} SNP input file using a
#' VCF file (Danecek et al. 2011). Missing data can bias demographic inference, 
#' \code{vcf2dadi} was created to address this problem, providing a customizable 
#' imputation framework specifically designed to work with GBS/RAD data.
#' If your VCF is not filtered, you can supply the function a whitelist of loci and a 
#' blacklist of individuals.
#' @inheritParams genomic_converter 
#' @inheritParams tidy_genomic_data
#' @inheritParams stackr_imputations_module

#' @param fasta.outgroup (optional) The fasta output file from STACKS. This file is 
#' required to use an outgroup. Default: \code{fasta.outgroup = NULL}.

#' @param fasta.ingroup (optional) The fasta output file from STACKS. This file is 
#' required to use with an outgroup. Leave empty if no outgroup. 
#' Default: \code{fasta.ingroup = NULL}.

#' @param sumstats.outgroup (optional) The sumstats output file from STACKS 
#' when running STACKS for the outgroup fasta file. This file is 
#' required to use an outgroup. Default: \code{sumstats.outgroup = NULL}.

#' @param sumstats.ingroup (optional) The sumstats output file from STACKS
#' when running STACKS for the ingroup fasta file.This file is 
#' required to use with an outgroup. Leave empty if no outgroup. 
#' Default: \code{sumstats.ingroup = NULL}.

#' @param dadi.input.filename (optional) Name of the \code{dadi} SNP input file 
#' written to the working directory. e.g. \code{dadi.file.tsv}. 
#' Default use date and time to make the file. If used, the file extension
#' need to finish with \code{.tsv or .txt}.


#' @details 
#' \strong{Input data:}
#'  
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
#' (e.g. from a VCF see \pkg{stackr} \code{\link{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID}, 
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#' 
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 or 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 or 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#' 
#' 
#' \strong{Imputations details:}
#' The imputations using Random Forest requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set 
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals 
#' will require 15 min.

#' @export
#' @rdname vcf2dadi
#' @import parallel
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread

#' @return The function returns a list with 1 or 2 objects (w/o imputations): 
#' `$dadi.no.imputation`
#' `$dadi.imputed`
#' The data frame are accessed form the list with `$`.

#' @examples
#' #See vignette `vignette_vcf2dadi` for more info.
#' \dontrun{
#' dadi.no.imputation <- vcf2dadi(
#' data = "batch_1.vcf", 
#' whitelist.markers = "whitelist.loci.txt",
#' strata = "strata.file.tsv",
#' pop.levels = c("PAN", "COS"),
#' common.markers = TRUE, 
#' fasta.ingroup = "batch_1.ingroup.fa", 
#' fasta.outgroup = "batch_1.outgroup.fa", 
#' sumstats.ingroup = "batch_1.sumstats.ingroup.tsv", 
#' sumstats.outgroup = "batch_1.sumstats.outgroup.tsv"
#' )
#' 
#' With Imputations:
#' dadi.files <- vcf2dadi(
#' data = "batch_1.vcf", 
#' whitelist.markers = "whitelist.loci.txt",
#' strata = "strata.file.tsv",
#' pop.levels = c("PAN", "COS"),
#' common.markers = TRUE, 
#' fasta.ingroup = "batch_1.ingroup.fa", 
#' fasta.outgroup = "batch_1.outgroup.fa", 
#' sumstats.ingroup = "batch_1.sumstats.ingroup.tsv", 
#' sumstats.outgroup = "batch_1.sumstats.outgroup.tsv",
#' imputation.method = "max", 
#' impute = "allele", 
#' imputations.group = "populations", 
#' num.tree = 100, 
#' iteration.rf = 10, 
#' split.number = 100, 
#' verbose = FALSE, 
#' parallel.core = 8
#' )
#' # to get the imputed data frame:
#' dadi.imputed.df <- dadi.files$dadi.imputed
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
#' @references Gutenkunst RN, Hernandez RD, Williamson SH, Bustamante CD (2009)
#' Inferring the Joint Demographic History of Multiple Populations 
#' from Multidimensional SNP Frequency Data (G McVean, Ed,). 
#' PLoS genetics, 5, e1000695.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and
#' Anne-Laure Ferchaud \email{annelaureferchaud@@gmail.com} 

vcf2dadi <- function(
  data,
  fasta.ingroup = NULL,
  fasta.outgroup = NULL,
  sumstats.ingroup = NULL,
  sumstats.outgroup = NULL,
  dadi.input.filename = NULL,
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
){
  
  # dadi unicode character: \u2202
  cat("#######################################################################\n")
  cat("########################## stackr::vcf2dadi ###########################\n")
  cat("#######################################################################\n")
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  if (is.null(imputation.method)) {
    message("vcf2dadi: no imputation...")
  } else {
    message("vcf2dadi: with imputations...")
  }
  
  # File type detection ********************************************************
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
  
  # Strata argument required for VCF and haplotypes files **********************
  if (data.type == "vcf.file" & is.null(strata)) stop("strata argument is required")
  if (data.type == "haplo.file") stop("This function is for biallelic dataset only")
  
  # Import input ***************************************************************
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
    distinct(INDIVIDUALS, POP_ID, .keep_all = TRUE)
  strata <- strata.df
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  # Compute count and Minor Allele Frequency -----------------------------------
  # We split the alleles here to prep for MAF
  # need to compute REF/ALT allele for non VCF file
  if (!tibble::has_name(input, "GT_VCF")) {
    ref.alt.alleles.change <- ref_alt_alleles(data = input)
    input <- left_join(input, ref.alt.alleles.change, by = c("MARKERS", "INDIVIDUALS"))
  }
  
  # MAF = the ALT allele in the VCF
  message("Computing the Allele Frequency Spectrum")
  
  input.count <- input %>% 
    group_by(MARKERS, POP_ID, REF, ALT) %>%
    summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT_VCF[GT_VCF == "0/0"])),
      PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
      QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
    ) %>%
    mutate(MAF = ((QQ*2) + PQ)/(2*N)) %>% 
    ungroup() %>%
    mutate(
      REF = stri_replace_all_fixed(str = REF, pattern = c("001", "002", "003", "004"), replacement = c("A", "C", "G", "T"), vectorize_all = FALSE),
      ALT = stri_replace_all_fixed(str = ALT, pattern = c("001", "002", "003", "004"), replacement = c("A", "C", "G", "T"), vectorize_all = FALSE)
    )
  
  # Function to make dadi input  data format ***********************************
  message("Preparing \u2202a\u2202i input SNP data format")
  write_dadi <- function(input, write.imputation, ...){
    # input <- input.count.imp # testing
    # input <- input.count # testing
    # input.bk <- input # testing
    
    if (is.null(fasta.ingroup)) {
      dadi.input <- suppressWarnings(
        input %>%
          group_by(MARKERS, POP_ID, REF, ALT) %>%
          mutate(
            A1 = (2 * PP) + PQ,
            A2 = (2 * QQ) + PQ
          ) %>%
          ungroup() %>% 
          select(POP_ID, Allele1 = REF, A1, Allele2 = ALT, A2, MARKERS) %>%
          tidyr::gather(ALLELE_GROUP, COUNT, -c(POP_ID, MARKERS, Allele1, Allele2)) %>% 
          tidyr::unite(POP, POP_ID, ALLELE_GROUP, sep = "_") %>% 
          group_by(MARKERS, Allele1, Allele2) %>% 
          tidyr::spread(data = ., key = POP, value = COUNT) %>%
          mutate(
            IN_GROUP = rep("---", n()), #version 2
            OUT_GROUP = rep("---", n())
          ) %>% 
          select(IN_GROUP, OUT_GROUP, Allele1, contains("A1"), Allele2, contains("A2"), MARKERS) %>% 
          arrange(MARKERS)
      )
    } # no fasta, no outgroup
    if (!is.null(fasta.ingroup)) {# With outgroup and ingroup fasta file
      message("using outgroup info")
      # Get the list of ref. allele in the vcf of the ingroup
      ref.allele.vcf.ingroup <- input %>% 
        ungroup() %>% 
        distinct(MARKERS, REF, .keep_all = TRUE) %>%
        arrange(MARKERS) %>% 
        tidyr::separate(MARKERS, c("CHROM", "LOCUS", "POS"), sep = "__") %>%
        distinct(CHROM, LOCUS, POS, REF, .keep_all = TRUE) %>%
        mutate(
          CHROM = as.character(stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>%
        arrange(CHROM, LOCUS, POS)
      # View(ref.allele.vcf.ingroup)
      
      # keep the list of markers
      markers <- ref.allele.vcf.ingroup %>% ungroup %>% select(-REF)
      
      # import out- and in- group fasta files ********************************
      fasta.outgroup <- data.table::fread(
        input = fasta.outgroup,
        sep = "\t",
        stringsAsFactors = FALSE,
        showProgress = FALSE,
        verbose = FALSE, 
        header = FALSE,
        col.names = "LOCUS"
      ) %>%
        mutate(
          ID.FILTER = rep(c("LOCUS", "SEQ"), times = n()/2),
          ANCESTRAL = rep("outgroup", times = n())
        )
      
      fasta.ingroup <- data.table::fread(
        input = fasta.ingroup,
        sep = "\t",
        stringsAsFactors = FALSE,
        showProgress = FALSE,
        verbose = FALSE, 
        header = FALSE,
        col.names = "LOCUS"
      ) %>%
        mutate(
          ID.FILTER = rep(c("LOCUS", "SEQ"), times = n()/2),
          ANCESTRAL = rep("ingroup", times = n())
        )
      
      fasta.data <- bind_rows(fasta.outgroup, fasta.ingroup)
      
      fasta.ingroup <- NULL
      fasta.outgroup <- NULL
      
      message("Preparing fasta sequences")
      fasta.prep <- fasta.data %>% 
        filter(ID.FILTER == "LOCUS") %>% 
        select(-ID.FILTER) %>% 
        bind_cols(fasta.data %>% 
                    filter(ID.FILTER == "SEQ") %>% 
                    select(-ID.FILTER, -ANCESTRAL, SEQUENCES = LOCUS)
        ) %>% 
        tidyr::separate(LOCUS, c("LOCUS", "ALLELE"), sep = "_Sample_", extra = "warn" ) %>%
        tidyr::separate(ALLELE, c("GARBAGE", "ALLELE"), sep = "_Allele_", extra = "warn" ) %>%
        tidyr::separate(ALLELE, c("ALLELE", "INDIVIDUALS"), sep = " ", extra = "warn" ) %>% 
        mutate(LOCUS = as.integer(stri_replace_all_fixed(LOCUS, pattern = ">CLocus_", replacement = "", vectorize_all = FALSE))) %>%
        select(-GARBAGE, -INDIVIDUALS) %>% 
        distinct(LOCUS, ALLELE, ANCESTRAL, SEQUENCES, .keep_all = TRUE) 
      
      fasta.data <- NULL
      
      # import out- and in- group sumstats files*****************************
      
      # Note: only one sumstats might be necessary to have the info needed, need checking
      sumstats.outgroup <- data.table::fread(
        input = sumstats.outgroup,
        sep = "\t",
        stringsAsFactors = FALSE, 
        skip = "Batch ID",
        select = c("Chr", "Locus ID", "BP", "Col"),
        showProgress = TRUE,
        verbose = FALSE
      ) %>% 
        as_data_frame() %>%
        select(CHROM = Chr, LOCUS = `Locus ID`, POS = BP, SNP_READ_POS = Col) %>%
        mutate(
          CHROM = as.character(stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>% 
        distinct(CHROM, LOCUS, SNP_READ_POS, .keep_all = TRUE) %>% 
        semi_join(markers, by = c("CHROM", "LOCUS", "POS")) %>% 
        arrange(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        mutate(ANCESTRAL = rep("outgroup", times = n())) %>%
        # the SNP position on the read is +1 for the fasta file
        mutate(SNP_READ_POS = SNP_READ_POS + 1)
      
      sumstats.ingroup <- data.table::fread(
        input = sumstats.ingroup,
        sep = "\t",
        stringsAsFactors = FALSE, 
        skip = "Batch ID",
        select = c("Chr", "Locus ID", "BP", "Col"),
        showProgress = TRUE,
        verbose = FALSE
      ) %>% 
        as_data_frame() %>%
        select(CHROM = Chr, LOCUS = `Locus ID`, POS = BP, SNP_READ_POS = Col) %>%
        mutate(
          CHROM = as.character(stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>% 
        distinct(CHROM, LOCUS, SNP_READ_POS, .keep_all = TRUE) %>% 
        semi_join(markers, by = c("CHROM", "LOCUS", "POS")) %>% 
        arrange(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        mutate(ANCESTRAL = rep("ingroup", times = n())) %>%
        # the SNP position on the read is +1 for the fasta file
        mutate(SNP_READ_POS = SNP_READ_POS + 1)
      
      markers <- NULL # remove unused object
      
      # When REF and ALT allele are equal in number, this can be problematic between ANCESTRAL group of alleles.
      # For those loci, we will set them based on both group (out- and in-group), so that there is no differences
      
      # THINKING: when ingroup REF and ALT equal, set to REF in VCF, and EQUAL in outgroup
      # when outgroup REF and ALT equal, set the REF based on the ingroup VCF
      
      max.length.read <- stri_length(str = fasta.prep[1,4])
      
      message("combining information from the fasta and sumstats file")
      ref.allele.vcf.ingroup.fasta <- ref.allele.vcf.ingroup %>% 
        # The sumstats is used to get the position of the markers along the sequence read
        left_join(
          sumstats.ingroup %>% 
            select(-ANCESTRAL)
          , by = c("CHROM", "LOCUS", "POS")
        ) %>% 
        mutate(SNP_READ_POS = stri_replace_na(SNP_READ_POS, replacement = "NA")) %>%
        filter(SNP_READ_POS != "NA") %>%
        mutate(SNP_READ_POS = as.integer(SNP_READ_POS)) %>% 
        arrange(CHROM, LOCUS, POS) %>%
        ungroup() %>% 
        # get the fasta sequence for those LOCUS...
        left_join(
          fasta.prep %>% 
            filter(ANCESTRAL == "ingroup") %>% 
            select (-ALLELE, -ANCESTRAL)
          , by = "LOCUS"
        ) %>%
        mutate(SNP_READ_POS = stri_replace_na(SNP_READ_POS, replacement = "NA")) %>%
        filter(SNP_READ_POS != "NA") %>%
        # get the flanking bases, left and right, of the SNP
        mutate(
          SNP_READ_POS = as.integer(SNP_READ_POS),
          FASTA_REF = stri_sub(SEQUENCES, from = SNP_READ_POS, to = SNP_READ_POS),
          IN_GROUP = stri_sub(SEQUENCES, from = SNP_READ_POS-1, to = SNP_READ_POS+1)
        ) %>%
        # remove lines with no match between all the alleles in the fasta file and the REF in the VCF
        filter(FASTA_REF == REF) %>%
        group_by(CHROM, LOCUS, POS, REF) %>%
        distinct(CHROM, LOCUS, POS, REF, .keep_all = TRUE) %>% 
        mutate(
          IN_GROUP = ifelse(SNP_READ_POS == max.length.read, stri_pad(IN_GROUP, width = 3, side = "right", pad = "-"), IN_GROUP),
          IN_GROUP = ifelse(SNP_READ_POS == 1, stri_pad(IN_GROUP, width = 3, side = "left", pad = "-"), IN_GROUP)
        )
      
      
      ingroup <- ref.allele.vcf.ingroup.fasta %>% 
        ungroup() %>% 
        select(CHROM, LOCUS, POS, IN_GROUP) %>% 
        # mutate(
        #   POS = stri_pad_left(str = POS, width = 8, pad = "0"),
        #   LOCUS = stri_pad_left(str = LOCUS, width = 8, pad = "0")
        # ) %>%
        arrange(CHROM, LOCUS, POS) %>%
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__")
      
      markers.id <- ref.allele.vcf.ingroup.fasta %>% 
        ungroup %>% 
        select(CHROM, LOCUS, POS, SNP_READ_POS)
      
      markers.id.ingroup.nuc <- ref.allele.vcf.ingroup.fasta %>% 
        ungroup %>% 
        select(CHROM, LOCUS, POS, IN_GROUP)
      
      # remove unused objects
      ref.allele.vcf.ingroup <- NULL
      ref.allele.vcf.ingroup.fasta <- NULL
      sumstats.ingroup <- NULL
      sumstats.outgroup <- NULL
      
      # Same thing but with outgroup
      ref.allele.outgroup.fasta <- markers.id %>% 
        left_join(
          fasta.prep %>% 
            filter(ANCESTRAL == "outgroup") %>% 
            select (-ALLELE, -ANCESTRAL)
          , by = "LOCUS"
        ) %>%
        mutate(SEQUENCES = stri_replace_na(SEQUENCES, replacement = "NA")) %>%
        filter(SEQUENCES != "NA") %>%
        # get the flanking bases, left and right, of the SNP
        mutate(OUT_GROUP = stri_sub(SEQUENCES, from = SNP_READ_POS-1, to = SNP_READ_POS+1)) %>%
        select(CHROM, LOCUS, POS, OUT_GROUP, SNP_READ_POS) %>%
        group_by(CHROM, LOCUS, POS, OUT_GROUP, SNP_READ_POS) %>%
        tally %>% 
        group_by(CHROM, LOCUS, POS) %>%
        filter(n == max(n))
      
      fasta.prep <- NULL # remove unused object
      
      # The problem is that some markers in the outgroup will have number of REF and ALT alleles equals...
      # Making the ancestral allele call ambiguous (50/50 chance of differences...)
      ambiguous.ancestral.allele <- ref.allele.outgroup.fasta %>%
        ungroup() %>%
        select(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        group_by(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        tally %>% 
        filter(n > 1) %>% 
        select(CHROM, LOCUS, POS, SNP_READ_POS) %>% 
        left_join(markers.id.ingroup.nuc, by = c("CHROM", "LOCUS", "POS")) %>% 
        rename(OUT_GROUP = IN_GROUP) %>% 
        ungroup() %>% 
        arrange(CHROM, LOCUS, POS)
      
      markers.id.ingroup.nuc <- NULL # remove unused object
      
      outgroup <- ref.allele.outgroup.fasta %>% 
        select(-n) %>% 
        anti_join(ambiguous.ancestral.allele, by = c("CHROM", "LOCUS", "POS")) %>%
        ungroup() %>% 
        arrange(CHROM, LOCUS, POS) %>% 
        bind_rows(ambiguous.ancestral.allele) %>%
        ungroup() %>% 
        arrange(CHROM, LOCUS, POS) %>%
        mutate(
          OUT_GROUP = ifelse(SNP_READ_POS == max.length.read, stri_pad(OUT_GROUP, width = 3, side = "right", pad = "-"), OUT_GROUP),
          OUT_GROUP = ifelse(SNP_READ_POS == 1, stri_pad(OUT_GROUP, width = 3, side = "left", pad = "-"), OUT_GROUP)
        ) %>% 
        select(CHROM, LOCUS, POS, OUT_GROUP) %>% 
        # mutate(
        #   POS = stri_pad_left(str = POS, width = 8, pad = "0"),
        #   LOCUS = stri_pad_left(str = LOCUS, width = 8, pad = "0")
        # ) %>%
        arrange(CHROM, LOCUS, POS) %>%
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__")
      
      ambiguous.ancestral.allele <- NULL # remove unused object
      ref.allele.outgroup.fasta <- NULL # remove unused object
      
      # common markers between ingroup and outgroup
      markers.ingroup <- ingroup %>% select(MARKERS)
      markers.outgroup <- outgroup %>% select(MARKERS)
      
      markers.ingroup.outgroup.common <- bind_rows(markers.ingroup, markers.outgroup) %>% 
        group_by(MARKERS) %>%
        tally %>% 
        filter(n == 2) %>%
        arrange(MARKERS) %>%
        distinct(MARKERS)
      
      markers.ingroup <- NULL
      markers.outgroup <- NULL
      
      message(stri_join("Number of markers common between in- and out- group = ", n_distinct(markers.ingroup.outgroup.common$MARKERS)))
      dadi.input <- suppressWarnings(
        input %>%
          group_by(MARKERS, POP_ID, REF, ALT) %>%
          mutate(
            A1 = (2 * PP) + PQ,
            A2 = (2 * QQ) + PQ
          ) %>%
          ungroup() %>%
          select(POP_ID, Allele1 = REF, A1, Allele2 = ALT, A2, MARKERS) %>%
          tidyr::gather(ALLELE_GROUP, COUNT, -c(POP_ID, MARKERS, Allele1, Allele2)) %>%
          tidyr::unite(POP, POP_ID, ALLELE_GROUP, sep = "_") %>%
          group_by(MARKERS, Allele1, Allele2) %>%
          tidyr::spread(data = ., key = POP, value = COUNT) %>% 
          select(Allele1, contains("A1"), Allele2, contains("A2"), MARKERS) %>%
          arrange(MARKERS) %>% 
          ungroup %>% 
          semi_join(markers.ingroup.outgroup.common, by = "MARKERS") %>% 
          inner_join(ingroup, by = "MARKERS") %>% 
          inner_join(outgroup, by = "MARKERS") %>% 
          select(IN_GROUP, OUT_GROUP, Allele1, contains("A1"), Allele2, contains("A2"), MARKERS) %>% 
          arrange(MARKERS)
      )
      
      # remove unused objects
      markers.id <- NULL
      markers.ingroup.outgroup.common <- NULL
      ingroup <- NULL
      outgroup <- NULL
    }
    
    # We need to modify the header to remove "_A1" and "_A2" that were necessary for spreading the info accross columns
    header.dadi <- colnames(dadi.input)
    colnames(dadi.input) <- stri_replace_all_fixed(str = header.dadi, pattern = c("_A1", "_A2"), replacement = "", vectorize_all = FALSE)
    
    if(is.null(dadi.input.filename)){
      # when filename is not provided will save the 'dadi.input' with date and time
      file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
      file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "_", ""), vectorize_all = FALSE)
      if (write.imputation == FALSE) {
        dadi.input.filename <- stri_paste("dadi_input_", file.date, ".tsv", sep = "")
      }
      
      if (write.imputation == TRUE) {
        dadi.input.filename.imp <- stri_paste("dadi_input_imputed_", file.date, ".tsv", sep = "")
      }
    } else {
      if (write.imputation == FALSE) {
        dadi.input.filename <- dadi.input.filename
      }
      
      if (write.imputation == TRUE) {
        dadi.input.filename.imp <- stri_replace_all_fixed(
          dadi.input.filename,
          pattern = c(".txt", ".tsv"),
          replacement = c("_imputed.txt", "_imputed.tsv"),
          vectorize_all = FALSE
        )
      }
    }
    file.version <- suppressWarnings(stri_paste("#\u2202a\u2202i SNP input file generated with stackr v.", packageVersion("stackr"), sep = ""))
    file.date <- suppressWarnings(stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = ""))
    file.header.line <- suppressWarnings(as.data.frame(stri_paste(file.version, file.date, sep = " ")))
    
    if (write.imputation == FALSE) {
      write_tsv(x = file.header.line, path = dadi.input.filename, append = FALSE, col_names = FALSE) # write the header line
      write_tsv(x = dadi.input, path = dadi.input.filename, append = TRUE, col_names = TRUE) # write the data frame
      message(stri_paste("\u2202a\u2202i input file name is: ", dadi.input.filename, "\n", "Saved here: ", getwd()))
    }
    
    if (write.imputation == TRUE) {
      write_tsv(x = file.header.line, path = dadi.input.filename.imp, append = FALSE, col_names = FALSE) # write the header line
      write_tsv(x = dadi.input, path = dadi.input.filename.imp, append = TRUE, col_names = TRUE) # write the data frame
      message(stri_paste("\u2202a\u2202i input file name is: ", dadi.input.filename.imp, "\n", "Saved here: ", getwd()))
    }
    
    return(dadi.input)
  } # End function write_dadi
  
  # Create list to store results
  res <- list()
  
  # without imputations (automatic)
  res$dadi.no.imputation <- write_dadi(input = input.count, write.imputation = FALSE)
  
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
      filename = NULL
    )
    
    # transform the imputed dataset  -------------------------------------------
    message("Computing the Allele Frequency Spectrum for the imputed data")
    
    input.count.imp <- input.imp %>% 
      group_by(MARKERS, POP_ID, REF, ALT) %>%
      summarise(
        N = as.numeric(n()),
        PP = as.numeric(length(GT_VCF[GT_VCF == "0/0"])),
        PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
        QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
      ) %>%
      mutate(MAF = ((QQ*2) + PQ)/(2*N))
    
    input.imp <- NULL # remove unused object
    
    # output in dadi
    res$dadi.imputed <- write_dadi(input = input.count.imp, write.imputation = TRUE)
  } # end dadi with imputations
  cat("############################## completed ##############################\n")
  
  return(res)
  } # End vcf2dadi
