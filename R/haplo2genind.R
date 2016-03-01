# Write a adegenet genind object from STACKS haplotypes file

# to get rid of notes in build check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Catalog ID", "Catalog.ID", "Catalog.ID = LOCUS", "Catalog.ID = `Catalog ID`", "Cnt", "HAPLOTYPES", "SAMPLES", "ALLELE", "ALLELE1", "ALLELE2", "GENOTYPE", "NUCLEOTIDES", "INDIVIDUALS", "POP_ID", "POLYMORPHISM", "POLYMORPHISM_MAX", "other", "strata", "hierarchy", "GROUP", ".", "CONSENSUS", "CONSENSUS_MAX", "MAX_COUNT_MARKERS"))


#' @name haplo2genind
#' @title Convert between batch_x.haplotypes.tsv and \code{adegenet} 
#' \code{\link[adegenet]{genind}} object
#' @description This function can first filter the haplotypes file 
#' with a whitelist of loci
#' and a blacklist of individuals (optional). Then it will convert the file
#' to a \code{adegenet} \code{\link[adegenet]{genind}} object.
#' Map-independent imputation using Random Forest or the most frequent category
#' is also available as an option.

#' @param data Haplotype file created in STACKS
#' e.g. \code{data = "batch_1.haplotypes.tsv"}, 

#' @param whitelist.markers (optional) A whitelist containing a LOCUS (integer) 
#' column header. The whitelist is in the working directory (e.g. "whitelist.txt").
#' Default: \code{whitelist.markers = NULL} for no whitelist of markers. 
#' In the stacks haplotype file, the \code{LOCUS = Catalog ID}.

#' @param blacklist.genotype (optional) Useful to erase genotype with below 
#' average quality, e.g. genotype with more than 2 alleles in diploid likely 
#' sequencing errors or genotypes with poor genotype likelihood or coverage. 
#' The blacklist file header should be \code{LOCUS}. If you need to work at the
#' SNP level, use \code{vcf2genind}.

#' @param erase (character) If you don't use the \code{blacklist.genotype} 
#' argument you still need to decide what to do with loci or genotypes with more
#' than 2 alleles. Loci with more than 2 alleles (paralogs and/or sequencing 
#' errors) will be removed with \code{erase = "loci"}. 
#' Default: \code{erase= "genotype"} will erase genotypes with more than 
#' 2 alleles by individual (paralogs and/or sequencing errors), keeping the loci
#' for other individuals.

#' @param common.markers (optional) Logical. Default = \code{FALSE}.
#' With \code{TRUE}, will keep markers genotyped in all the populations.


#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").

#' @param pop.levels (required) A character string with your populations ordered.
#' @param pop.labels (optional) A character string for your populations labels.
#' If you need to rename sampling sites in \code{pop.levels} or combined sites/pop
#' into a different names, here is the place.
#' @param pop.id.start The start of your population id
#' in the name of your individual sample. Your individuals are identified 
#' in this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
#' then, \code{pop.id.start} = 5. If you didn't name your individuals
#' with the pop id in it, use the \code{strata} argument. 
#' @param pop.id.end The end of your population id
#' in the name of your individual sample. Your individuals are identified 
#' in this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
#' then, \code{pop.id.end} = 7. If you didn't name your individuals
#' with the pop id in it, use the \code{strata} argument.
#' @param strata (optional) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. Default: \code{strata = NULL}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping.
#' @param hierarchy (optional) A formula that explicitely defines hierarchical levels 
#' in your strata. See \code{\link[adegenet]{genind}} for details.

#' @param pop.select (string) Conduct the assignment analysis on a
#' selected list of populations. Default = \code{NULL} for no selection and keep
#' all population.
#' e.g. \code{pop.select = "QUE"} to select QUE population samples.
#' \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population samples.

#' @param imputation.method Should a map-independent imputations of markers be
#' computed. Available choices are: (1) \code{FALSE} for no imputation.
#' (2) \code{"max"} to use the most frequent category for imputations.
#' (3) \code{"rf"} using Random Forest algorithm. 
#' Default: \code{imputation.method = FALSE}.
#' @param impute (character) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
#' @param imputations.group \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default = \code{"populations"}.
#' @param num.tree The number of trees to grow in Random Forest. Default is 100.
#' @param iteration.rf The number of iterations of missing data algorithm
#' in Random Forest. Default is 10.
#' @param split.number Non-negative integer value used to specify
#' random splitting in Random Forest. Default is 100.
#' @param verbose Logical. Should trace output be enabled on each iteration
#' in Random Forest ? Default is \code{FALSE}.
#' @param parallel.core (optional) The number of core for OpenMP shared-memory parallel
#' programming of Random Forest imputations. For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.
#' @details You need to have either the \code{pop.id.start} and \code{pop.id.end}
#' or the \code{strata} argument, to identify your populations.
#' The imputations using Random Forest requires more time to compute
#' and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals
#' will require 15 min.

#' @return When no imputation is selected an object of the 
#' class \code{\link[adegenet]{genind}} is returned.
#' When imputation is selected a list with 2 objects is returned
#' and accessed with \code{$no.imputation} or \code{$imputed}.
#' @export
#' @rdname haplo2genind
#' @import dplyr
#' @import parallel
#' @import stringi
#' @import adegenet
#' @importFrom data.table fread

#' @examples
#' \dontrun{
#' snowcrab <- haplo2genind(
#' data = "batch_1.haplotypes.tsv",
#' whitelist.markers = "whitelist.txt",
#' common.markers = TRUE,
#' blacklist.id = "blacklist.id.lobster.tsv",
#' pop.levels = c("PAN", "COS")
#' pop.id.start = 5, pop.id.end = 7,
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
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @seealso \code{adegenet} is available on CRAN \url{http://cran.r-project.org/web/packages/adegenet/} and github \url{https://github.com/thibautjombart/}
#' \code{randomForestSRC} is available on CRAN \url{http://cran.r-project.org/web/packages/randomForestSRC/} and github \url{https://github.com/ehrlinger/randomForestSRC}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

haplo2genind <- function(data, 
                         whitelist.markers = NULL,
                         blacklist.genotype = NULL,
                         erase = "genotype",
                         common.markers = NULL,
                         blacklist.id = NULL,
                         pop.levels,
                         pop.labels,
                         pop.id.start, 
                         pop.id.end,
                         pop.select = NULL,
                         strata = NULL,
                         hierarchy = NULL,
                         imputation.method = FALSE,
                         impute = "genotypes",
                         imputations.group = "populations",
                         num.tree = 100,
                         iteration.rf = 10,
                         split.number = 100,
                         verbose = FALSE,
                         parallel.core = NULL
) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Haplotype file is missing")
  if (missing(whitelist.markers)) whitelist.markers <- NULL # no Whitelist
  if (missing(blacklist.genotype)) blacklist.genotype <- NULL # no genotype to erase
  if (missing(erase)) erase <- "genotype"
  if (missing(common.markers)) common.markers <- FALSE
  if (missing(blacklist.id)) blacklist.id <- NULL # No blacklist of ID
  if (missing(pop.levels)) stop("pop.levels required")
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels
  if (missing(pop.id.start)) pop.id.start <- NULL
  if (missing(pop.id.end)) pop.id.end <- NULL
  if (missing(pop.select)) pop.select <- NULL
  if (missing(strata)) strata <- NULL
  if (missing(hierarchy)) hierarchy <- NULL
  if (missing(imputation.method)) imputation.method <- FALSE
  if (missing(imputations.group)) imputations.group <- "populations"
  if (imputation.method != FALSE & missing(impute)) stop("impute argument is necessary")
  if (imputation.method == FALSE & missing(impute)) impute <- NULL
  if (missing(num.tree)) num.tree <- 100
  if (missing(iteration.rf)) iteration.rf <- 10
  if (missing(split.number)) split.number <- 100
  if (missing(verbose)) verbose <- FALSE
  if (missing(parallel.core) | is.null(parallel.core)) parallel.core <- detectCores()-1
  
    if (imputation.method == "FALSE") {
      message("haplo2genind: no imputation...")
    } else {
      message("haplo2genind: with imputations...")
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
    }
    
    whitelist.markers <- select(.data = whitelist.markers, LOCUS)
    columns.names.whitelist <- colnames(whitelist.markers)
    
    
    # Import blacklist id ********************************************************
    if (is.null(blacklist.id)) { # No blacklist of ID
      message("Blacklisted individuals: no")
    } else { # With blacklist of ID
      message("Blacklisted individuals: yes")
      blacklist.id <- read_tsv(blacklist.id, col_names = TRUE)
    }
    
    # Import data ****************************************************************
    message("Importing the stacks haplotype file")
    input <- data.table::fread(
      input = data, 
      sep = "\t", 
      header = TRUE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      showProgress = TRUE,
      data.table = FALSE
    ) %>% 
      as_data_frame() %>% 
      select(-Cnt) %>% 
      rename(LOCUS = `Catalog ID`) %>%
      tidyr::gather(INDIVIDUALS, GT, -LOCUS)
    
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
    
    if (is.null(strata)){
      input <- input %>%
        mutate( # Make population ready
          POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
          POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered = TRUE),
          INDIVIDUALS =  as.character(INDIVIDUALS)
        )
    } else { # Make population ready with the strata provided
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA)
      
      input <- input %>%
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        left_join(strata.df, by = "INDIVIDUALS") %>% 
        mutate(POP_ID = factor(POP_ID, levels = unique(pop.labels), ordered =TRUE))
    }
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # Blacklist genotypes ********************************************************
    if (is.null(blacklist.genotype)) { # no Whitelist
      message("Erasing genotype: no")
    } else {
      message("Erasing genotype: yes")
      blacklist.genotype <- read_tsv(blacklist.genotype, col_names = TRUE)
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
      if ("CHROM" %in% columns.names.blacklist.genotype) {
        columns.names.blacklist.genotype$CHROM <- as.character(columns.names.blacklist.genotype$CHROM)
      }
      
      blacklist.genotype <- select(.data = blacklist.genotype, INDIVIDUALS, LOCUS)
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
      
      # control check to keep only individuals in pop.select
      if (!is.null(pop.select)) {
        message("Control check to keep only individuals present in pop.select")
        # updating the blacklist.genotype
        if (is.null(strata)){
          blacklist.genotype <- suppressWarnings(
            blacklist.genotype  %>% 
              mutate( # Make population ready
                POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
                POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = F), levels = unique(pop.labels), ordered =T),
                INDIVIDUALS =  as.character(INDIVIDUALS) 
              ) %>% 
              filter(POP_ID %in% pop.select) %>% 
              select(-POP_ID)
          )
        } else {
          blacklist.genotype <- suppressWarnings(
            blacklist.genotype %>%
              mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
              left_join(strata.df, by = "INDIVIDUALS") %>% 
              filter(POP_ID %in% pop.select) %>% 
              select(-POP_ID)
          )
        }
      }
      
      # control check to keep only whitelisted markers from the blacklist of genotypes
      if (!is.null(whitelist.markers)) {
        blacklist.genotype <- blacklist.genotype
        message("Control check to keep only whitelisted markers present in the blacklist of genotypes to erase.")
        # updating the whitelist of markers to have all columns that id markers
        whitelist.markers.ind <- input %>% select(LOCUS, INDIVIDUALS) %>% distinct(LOCUS, INDIVIDUALS)
        
        # updating the blacklist.genotype
        blacklist.genotype <- suppressWarnings(semi_join(whitelist.markers.ind, blacklist.genotype, by = columns.names.blacklist.genotype))
        columns.names.blacklist.genotype <- colnames(blacklist.genotype)
      }
      
      # control check to remove blacklisted individuals from the blacklist of genotypes
      if (!is.null(blacklist.id)) {
        message("Control check to remove blacklisted individuals present in the blacklist of genotypes to erase.")
        blacklist.genotype <- suppressWarnings(anti_join(blacklist.genotype, blacklist.id, by = "INDIVIDUALS"))
        columns.names.blacklist.genotype <- colnames(blacklist.genotype)
      }
      
      # Add one column that will allow to include the blacklist in the dataset 
      # by x column(s) of markers
      blacklist.genotype <- mutate(.data = blacklist.genotype, ERASE = rep("erase", n()))
      
      input <- suppressWarnings(
        input %>%
          full_join(blacklist.genotype, by = columns.names.blacklist.genotype) %>%
          mutate(ERASE = stri_replace_na(str = ERASE, replacement = "ok"))
      )
      
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "-", GT)) %>% 
        select(-ERASE)
      
    } # End erase genotypes
    
    # Paralogs *******************************************************************
    message("Looking for paralogs and/or errors...")
    if (erase == "loci") {
      message("Loci with more than 2 alleles will be removed")
      
      paralogs <- input %>%
        mutate(POLYMORPHISM = stri_count_fixed(GT, "/")) %>%
        group_by(LOCUS) %>%
        summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
        filter(POLYMORPHISM_MAX > 1) %>%
        group_by(LOCUS) %>%
        select(LOCUS) %>%
        distinct(LOCUS)
      
      input <- suppressWarnings(
        input %>%
          filter(subset =! LOCUS %in% paralogs$LOCUS)
      )
      message(stri_join("Found and/or removed", n_distinct(paralogs$LOCUS), "paralogs and/or errors", sep = " "))
    } 
    if (erase == "genotype") {
      message("Erasing genotypes with more than 2 alleles")
      
      # get the number of genotypes...
      haplo.number <- input %>%
        filter(GT != "-") %>%
        select(GT)
      
      input <- input %>% 
        mutate(POLYMORPHISM = stri_count_fixed(GT, "/"))
      
      erased.genotype.number <- length(input$INDIVIDUALS[input$POLYMORPHISM > 1])
      total.genotype.number.haplo <- length(haplo.number$GT)
      percent.haplo <- paste(round(((erased.genotype.number/total.genotype.number.haplo)*100), 4), "%", sep = " ")
      message(stri_paste("Out of a total of ", total.genotype.number.haplo, " genotypes, ", percent.haplo, " (", erased.genotype.number, ")"," will be erased"))
      
      # Erasing genotype with the blacklist
      message("Erasing... Erasing...")
      input <- suppressWarnings(
        input %>%
          mutate(GT = ifelse(POLYMORPHISM > 1, "-", GT)) %>% 
          select(-POLYMORPHISM)
      )
    }
    
    # consensus ******************************************************************
    consensus <- input %>%
      group_by(LOCUS) %>% 
      mutate(CONSENSUS = stri_count_fixed(GT, "consensus")) %>%
      summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
      filter(CONSENSUS_MAX > 0) %>%
      select(LOCUS)
    
    input <- input %>%
      filter(!LOCUS %in% consensus$LOCUS)
    
    # dump unused object
    consensus <- NULL
    blacklist.id <- NULL
    whitelist.markers <- NULL
    whitelist.markers.ind <- NULL
    blacklist.genotype <- NULL 
    haplo.number <- NULL
    
    # Unique markers id **********************************************************
    # Recycled codes from other functions, kept for internal consistencies
    input <- input %>% rename(MARKERS = LOCUS)
    
    # Markers in common between all populations (optional) *********************
    if (common.markers == TRUE) { # keep only markers present in all pop
      message("Using markers common in all populations:")
      pop.number <- n_distinct(input$POP_ID)
      
      pop.filter <- input %>% filter(GT != "-")
      
      pop.filter <- pop.filter %>% 
        group_by(MARKERS) %>%
        filter(n_distinct(POP_ID) == pop.number) %>%
        arrange(MARKERS) %>%
        select(MARKERS) %>%
        distinct(MARKERS)
      
      message(stri_join("Number of original markers = ", n_distinct(input$MARKERS), 
                        "\n", "Number of markers present in all the populations = ", 
                        n_distinct(pop.filter$MARKERS), "\n", 
                        "Number of markers removed = ", 
                        n_distinct(input$MARKERS) - n_distinct(pop.filter$MARKERS))
      )
      input <- suppressWarnings(input %>% semi_join(pop.filter, by = "MARKERS"))
      pop.filter <- NULL # ununsed object
    } # End common markers
    
    # Change the genotype coding  **********************************************
    message("Recoding genotypes")
    input <- suppressWarnings(
      input %>%
        mutate(GT = stri_replace_all_fixed(GT, "-", "000/000", vectorize_all=F)) %>% 
        arrange(MARKERS) %>% 
        tidyr::separate(
          col = GT, into = c("A1", "A2"), 
          sep = "/", extra = "drop", remove = TRUE
        ) %>%
        mutate(
          A2 = stri_replace_na(str = A2, replacement = "000"),
          A2 = ifelse(A2 == "000", A1, A2)
        ) %>%
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
        ungroup()
    )
    
    input <- suppressWarnings(
      input %>%
        filter(GT != "000") %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>% # this reintroduce the missing, but with NA
        ungroup() %>% 
        plyr::colwise(.fun = factor, exclude = NA)(.)
    )
    
    input <- suppressWarnings(
      input %>%
        ungroup() %>% 
        mutate_each(funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
        ungroup() %>%
        tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
        mutate(GT = stri_pad_left(str = GT, width = 3, pad = "0"))
    )
    
    genind.prep <- suppressWarnings(
      input %>% 
        mutate(GT = stri_replace_na(str = GT, replacement = "000")) %>%
        filter(GT != "000") %>%
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
        group_by(INDIVIDUALS, MARKERS) %>%
        mutate(MAX_COUNT_MARKERS = max(COUNT, na.rm = TRUE)) %>%
        ungroup() %>% 
        mutate(COUNT = ifelse(MAX_COUNT_MARKERS == 0, "erase", COUNT)) %>%
        select(-MAX_COUNT_MARKERS) %>% 
        mutate(COUNT = replace(COUNT, which(COUNT == "erase"), NA)) %>% 
        arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES) %>% 
        tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
        tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>% 
        arrange(POP_ID, INDIVIDUALS)
    )
    
    # genind constructor **********************************************************
    ind <- as.character(genind.prep$INDIVIDUALS)
    pop <- genind.prep$POP_ID
    genind.df <- genind.prep %>% ungroup() %>% 
      select(-c(INDIVIDUALS, POP_ID))
    rownames(genind.df) <- ind
    loc.names <- colnames(genind.df)
    strata <- genind.prep %>% ungroup() %>% select(INDIVIDUALS, POP_ID) %>% distinct(INDIVIDUALS, POP_ID)
    prevcall <- match.call()
    res <- genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = hierarchy)
    message("A large 'genind' object was created in your Environment")
    
    # sum <- summary(res) # test
    # sum$NA.perc # test
    
    ind <- NULL
    pop <- NULL
    genind.df <- NULL
    genind.prep <- NULL
    loc.names <- NULL
  
  # Imputations **************************************************************
  if (imputation.method != "FALSE") {
    message("Preparing the data for imputations")
    if (impute == "genotype") {
      input.prep <- input %>%
        tidyr::spread(data = ., key = ALLELES, value = GT) %>%
        tidyr::unite(data = ., GT, A1, A2, sep = "", remove = TRUE) %>% 
        mutate(GT = replace(GT, which(GT == "NANA"), NA)) %>%
        arrange(POP_ID, INDIVIDUALS, MARKERS) %>% 
        group_by(INDIVIDUALS, POP_ID) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>%
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
    }
    if (impute == "allele") {
      input.prep <- input %>%
        group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
    }
    
    # Imputation with Random Forest
    if (imputation.method == "rf") {
      # Parallel computations options
      options(rf.cores = parallel.core, mc.cores = parallel.core)
      
      # imputations using Random Forest with the package randomForestSRC
      impute_genotype_rf <- function(x) {
        randomForestSRC::impute.rfsrc(data = x,
                                      ntree = num.tree,
                                      nodesize = 1,
                                      nsplit = split.number,
                                      nimpute = iteration.rf,
                                      do.trace = verbose)
      } # End of imputation function
      
      # Random Forest by pop
      if (imputations.group == "populations") {
        message("Imputations computed by populations, take a break...")
        df.split.pop <- split(x = input.prep, f = input.prep$POP_ID) # slip data frame by population
        pop.list <- names(df.split.pop) # list the pop
        imputed.dataset <-list() # create empty list
        
        # Function to go through the populations
        impute_rf_pop <- function(pop.list, ...){
          sep.pop <- df.split.pop[[pop.list]]
          sep.pop <- suppressWarnings(
            plyr::colwise(factor, exclude = NA)(sep.pop)
          )
          # message of progress for imputations by population
          message(paste("Completed imputations for pop ", pop.list, sep = ""))
          # imputed.dataset[[i]] <- impute_markers_rf(sep.pop) # test with foreach
          imputed.dataset <- impute_genotype_rf(sep.pop)
          return(imputed.dataset)
        } # End impute_rf_pop
        
        input.imp <- list()
        input.imp <- parallel::mclapply(
          X = pop.list, 
          FUN = impute_rf_pop, 
          mc.preschedule = FALSE, 
          mc.silent = FALSE, 
          mc.cores = parallel.core
        )
        
        # Compiling the results
        message("Compiling imputations results")
        input.imp <- suppressWarnings(bind_rows(input.imp))
        
        # Second round of imputations (globally) to remove introduced NA 
        # In case that some pop don't have the markers
        input.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(input.imp)) # Make the columns factor
        input.imp <- impute_genotype_rf(input.imp) # impute globally
        
        # dump unused objects
        df.split.pop <- NULL
        pop.list <- NULL
        sep.pop <- NULL
        imputed.dataset <- NULL
        input.prep <- NULL
        
      } # End imputation RF populations
      # Random Forest global
      if (imputations.group == "global") { # Globally (not by pop_id)
        message("Imputations computed globally, take a break...")
        input.prep <- plyr::colwise(factor, exclude = NA)(input.prep)
        input.imp <- impute_genotype_rf(input.prep)
        
        input.prep <- NULL # remove unused object
      } # End imputation RF global
      if (impute == "genotype") {
        input.imp <- suppressWarnings(
          input.imp %>% 
            tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% 
            tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
            tidyr::gather(data = ., key = ALLELES, value = GT, -c(INDIVIDUALS, POP_ID, MARKERS)) #%>% 
            # select(-ALLELES)
        )
      }
      if (impute == "allele") {
        input.imp2 <- suppressWarnings(
          input.imp %>% 
            tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) 
          # %>% 
            # select(-ALLELES)
        )
      }
      
    } # End imputation RF
    # Imputation using the most common genotype
    if (imputation.method == "max") { # End imputation max
      if (imputations.group == "populations") {
        message("Imputations computed by populations")
        
        if (impute == "genotype"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
              group_by(MARKERS, POP_ID) %>%
              mutate(
                GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                GT = replace(GT, which(GT == "NA"), NA)
              ) %>%
              # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
              # will take the global observed values by markers for those cases.
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              # group_by(INDIVIDUALS, POP_ID) %>% 
              # tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
        }
        
        if (impute == "allele"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
              group_by(MARKERS, POP_ID) %>%
              mutate(
                GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                GT = replace(GT, which(GT == "NA"), NA)
              ) %>%
              # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
              # will take the global observed values by markers for those cases.
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              # group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
              # tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
        }
        
        input.prep <- NULL # remove unused object
        
      } # End imputation max populations 
      if (imputations.group == "global") {
        # Globally (not by pop_id)
        message("Imputations computed globally")
        if (impute == "genotype"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              # group_by(INDIVIDUALS, POP_ID) %>% 
              # tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
        }
        
        if (impute == "allele"){
          input.imp <- suppressWarnings(
            input.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              # group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
              # tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
        }
        
        input.prep <- NULL # remove unused object
      } # End imputation max global 
    } # End imputations max
    
    # prepare the imputed dataset for adegenet
    message("Preparing imputed data set for adegenet...")
    genind.prep.imp <- suppressWarnings(
      input.imp %>%
        # tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
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
        arrange(POP_ID, INDIVIDUALS)
    )
    input.imp <- NULL
    # results WITH imputations ------------------------------------------------
    # 1) the genind without imputations is modified and put in a new list
    no.imputation <- res
    res <- list()
    res$no.imputation <- no.imputation
    
    # 2) the genind with imputations
    ind <- as.character(genind.prep.imp$INDIVIDUALS)
    pop <- genind.prep.imp$POP_ID
    genind.df <- genind.prep.imp %>% ungroup() %>% 
      select(-c(INDIVIDUALS, POP_ID))
    rownames(genind.df) <- ind
    loc.names <- colnames(genind.df)
    strata <- genind.prep.imp %>% ungroup() %>% select(INDIVIDUALS, POP_ID) %>% distinct(INDIVIDUALS, POP_ID)
    
    # genind constructor
    prevcall <- match.call()
    res$imputed <- genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = NULL)
    
    # sum <- summary(res$imputed) # test
    # sum$NA.perc # test
    
    ind <- NULL
    pop <- NULL
    genind.df <- NULL
    genind.prep <- NULL
    genind.prep.imp <- NULL
    input.imp <- NULL # remove unused object
    message("A large 'genind' object was created in your Environment (with and without imputations)")
  } # End imputations
  # outout results -------------------------------------------------------------
  return(res)
} # End haplo2genind function
