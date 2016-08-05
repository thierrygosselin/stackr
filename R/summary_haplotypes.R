## Summary and tables

#' @title Haplotypes file summary
#' @description STACKS batch_x.haplotypes.tsv file summary.
#' The output of the function is a summary table for populations with:
#' \enumerate{
#' \item putative sequencing errors: loci with > 2 alleles and/or potential 
#' paralogs. Genotypes with more than 2 alleles are 
#' erased before estimating the subsequent statistics.
#' \item consensus loci
#' \item monomorphic loci
#' \item polymorphic loci
#' \item haplotypes statistics for the observed and expected homozygosity and
#' heterozygosity.
#' \item Wright’s inbreeding coefficient (Fis)
#' \item \strong{IBDG}: a proxy measure of the realized proportion of the genome
#' that is identical by descent
#' \item \strong{FH measure}: based on the excess in the observed number of homozygous
#' genotypes within an individual relative to the mean number of homozygous 
#' genotypes expected under random mating (Keller et al., 2011; Kardos et al., 2015).
#' \item \strong{Pi}: the nucleotide diversity measured here consider the 
#' consensus loci in the catalog (no variation between population sequences).
#' }

#' @param data The 'batch_x.haplotypes.tsv' created by STACKS.

#' @param strata A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. 
#' The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' If you have already run 
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data, 
#' the strata file is similar to a stacks `population map file`, make sure you 
#' have the required column names (\code{INDIVIDUALS} and \code{STRATA}).
#' Default: \code{strata = NULL}.

#' @param whitelist.markers (optional) A whitelist of loci with a column name
#' 'LOCUS'. The whitelist is located in the global environment or in the 
#' directory (e.g. "whitelist.txt"). 
#' If a whitelist is used, the files written to the directory will have 
#' '.filtered' in the filename. 
#' Default: \code{whitelist.markers = NULL}.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the global environment 
#' or in the directory (with "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

#' @param pop.levels An optional character string with your populations ordered.
#' Default: \code{pop.levels = NULL}.
#' @param pop.labels An optional character string with new populations names.
#' Default: \code{pop.labels = NULL}.

#' @param read.length (number) The length in nucleotide of your reads 
#' (e.g. \code{read.length = 100}).

#' @param parallel.core (optional) The number of core for parallel
#' programming during Pi calculations. 
#' Default: \code{parallel.core = detectCores() - 1}.

#' @importFrom stringdist stringdist
#' @importFrom utils combn
#' @importFrom stats lm
#' @importFrom parallel mclapply
#' @return The function returns a list with:
#' \enumerate{
#' \item $summary
#' \item $paralogs.pop
#' \item $paralogs.loci
#' \item $consensus.pop
#' \item $consensus.loci
#' \item $fh.pi.individuals: the individual's info for FH and Pi
#' }
#' Also in the list 3 plots:
#' \enumerate{
#' \item $scatter.plot
#' \item $boxplot.pi
#' \item $boxplot.fh
#' }
#' use $ to access each #' objects in the list.
#' The function potentially write 3 files in the working directory:
#' blacklist of unique putative paralogs and unique consensus loci 
#' and a summary of the haplotypes file by population.

#' @examples
#' \dontrun{
#' sum <- summary_haplotypes(
#' data = "batch_1.haplotypes.tsv", 
#' strata = "strata_brook_charr.tsv", 
#' pop.levels = c("SKY", "LIM", "TWE"), 
#' read.length = 90)
#' }


#' @rdname summary_haplotypes
#' @export 

#' @references Keller MC, Visscher PM, Goddard ME (2011)
#' Quantification of inbreeding due to distant ancestors and its detection
#'  using dense single nucleotide polymorphism data. Genetics, 189, 237–249.
#' @references Kardos M, Luikart G, Allendorf FW (2015)
#' Measuring individual inbreeding in the age of genomics: marker-based
#' measures are better than pedigrees. Heredity, 115, 63–72.
#' @references Nei M, Li WH (1979)
#' Mathematical model for studying genetic variation in terms of 
#' restriction endonucleases. 
#' Proceedings of the National Academy of Sciences of 
#' the United States of America, 76, 5269–5273.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and 
#' Anne-Laure Ferchaud \email{annelaureferchaud@@gmail.com}


summary_haplotypes <- function(
  data,
  strata,
  whitelist.markers = NULL, 
  blacklist.id = NULL, 
  pop.levels = NULL,
  pop.labels = NULL,
  read.length,
  parallel.core = detectCores() - 1
) {
  
  cat("#######################################################################\n")
  cat("#################### stackr: summary_haplotypes #######################\n")
  cat("#######################################################################\n")
  
  
  if (missing(data)) stop("data argument is missing")
  if (missing(strata)) stop("strata argument is required")
  if (missing(read.length)) stop("read.length argument is required")
  
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  # Import haplotype file ------------------------------------------------------
  message("Importing STACKS haplotypes file")
  number.columns <- max(count.fields(data, sep = "\t"))
  
  haplotype <- data.table::fread(
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
    as_data_frame() %>% 
    select(-Cnt) %>% 
    rename(LOCUS = `Catalog ID`)
  
  haplotype <- data.table::melt.data.table(
    data = data.table::as.data.table(haplotype), 
    id.vars = "LOCUS", 
    variable.name = "INDIVIDUALS",
    variable.factor = FALSE,
    value.name = "HAPLOTYPES"
  ) %>% 
    as_data_frame()
  
  number.columns <- NULL
  
  haplotype$INDIVIDUALS = stri_replace_all_fixed(
    str = haplotype$INDIVIDUALS, 
    pattern = c("_", ":"), 
    replacement = c("-", "-"), 
    vectorize_all = FALSE
  )
  
  # Import whitelist of markers-------------------------------------------------
  if (!is.null(whitelist.markers)) { # no Whitelist
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
    whitelist.markers <- select(.data = whitelist.markers, LOCUS)
    columns.names.whitelist <- colnames(whitelist.markers)
  }
  
  # Import blacklist id --------------------------------------------------------
  if (!is.null(blacklist.id)) { # No blacklist of ID
    blacklist.id <- read_tsv(blacklist.id, col_names = TRUE)
  }
  
  # population levels and strata------------------------------------------------
  message("Making STACKS haplotypes file population-wise...")
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      # message("strata file: yes")
      number.columns.strata <- max(count.fields(strata, sep = "\t"))
      col.types <- stri_paste(rep("c", number.columns.strata), collapse = "")
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>% 
        rename(POP_ID = STRATA)
    } else {
      # message("strata object: yes")
      colnames(strata) <- stri_replace_all_fixed(
        str = colnames(strata), 
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
    # Remove potential whitespace in pop_id
    strata.df$POP_ID <- stri_replace_all_fixed(strata.df$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  # Check with strata and pop.levels/pop.labels
  if (!is.null(strata) & !is.null(pop.levels)) {
    if (length(levels(factor(strata.df$POP_ID))) != length(pop.levels)) {
      stop("The number of groups in your strata file must match the number of groups in pop.levels")
    }
  }
  # Filter with whitelist of markers -------------------------------------------
  if (!is.null(whitelist.markers)) {
    message("Filtering with whitelist of markers")
    haplotype <- suppressWarnings(semi_join(haplotype, whitelist.markers, by = columns.names.whitelist))
  }
  
  # Filter with blacklist of individuals ---------------------------------------
  if (!is.null(blacklist.id)) {
    message("Filtering with blacklist of individuals")
    
    blacklist.id$INDIVIDUALS <- stri_replace_all_fixed(
      str = blacklist.id$INDIVIDUALS, 
      pattern = c("_", ":"), 
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )
    
    haplotype <- suppressWarnings(anti_join(haplotype, blacklist.id, by = "INDIVIDUALS"))
  }
  
  # Population levels and strata -----------------------------------------------
  strata.df$INDIVIDUALS = stri_replace_all_fixed(
    str = strata.df$INDIVIDUALS, 
    pattern = c("_", ":"), 
    replacement = c("-", "-"), 
    vectorize_all = FALSE
  )
  
  if (is.null(pop.levels)) { # no pop.levels
    strata.df <- mutate(.data = strata.df, POP_ID = factor(POP_ID))
  } else {# with pop.levels
    strata.df <- mutate(
      .data = strata.df,
      POP_ID = factor(
        stri_replace_all_regex(
          POP_ID, stri_paste("^", pop.levels, "$", sep = ""), pop.labels, vectorize_all = FALSE),
        levels = unique(pop.labels), ordered = TRUE
      )
    )
  }
  
  haplotype <- left_join(x = haplotype, y = strata.df, by = "INDIVIDUALS")
  
  # Locus with concensus alleles -----------------------------------------------
  message("Scanning for consensus markers...")
  
  # consensus.pop <- haplotype %>%
  #   mutate(CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus")) %>%
  #   group_by(LOCUS, POP_ID) %>%
  #   summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
  #   filter(CONSENSUS_MAX > 0)
  
  consensus.pop <- haplotype %>%
    filter(HAPLOTYPES == "consensus") %>%
    distinct(LOCUS, POP_ID) %>% 
    mutate(CONSENSUS = rep("consensus", times = n())) %>% 
    arrange(as.numeric(LOCUS), POP_ID)
  
  # Create a list of consensus loci
  if (length(consensus.pop$CONSENSUS) > 0) {
    blacklist.loci.consensus <- consensus.pop %>%
      ungroup %>% 
      distinct(LOCUS) %>%
      arrange(as.numeric(LOCUS))
    
    write_tsv(
      x = blacklist.loci.consensus, 
      path = "blacklist.loci.consensus.txt", 
      col_names = TRUE
    )
  }
  
  # Individuals per pop
  ind.pop <- haplotype %>%
    distinct(INDIVIDUALS, .keep_all = TRUE) %>% 
    group_by(POP_ID) %>%
    tally() %>% 
    rename(NUM = n)
  
  ind.tot <- haplotype %>%
    distinct(INDIVIDUALS) %>% 
    tally() %>% 
    rename(NUM = n) %>% 
    mutate(POP_ID = "OVERALL")
  
  sample.number <- suppressWarnings(bind_rows(ind.pop, ind.tot))
  
  # unused arguments
  ind.pop <- ind.tot <- NULL
  
  # Sequencing errors and Paralogs ---------------------------------------------
  # Locus with > 2 alleles by individuals 
  message("Scanning for genotypes with > 2 alleles (paralogs and sequencing errors)...")
  paralogs.ind <- haplotype %>% 
    filter(stri_count_fixed(HAPLOTYPES, "/") > 1) %>% 
    arrange(LOCUS, POP_ID, INDIVIDUALS) %>% 
    select(LOCUS, POP_ID, INDIVIDUALS, HAPLOTYPES)
  
  if (length(paralogs.ind$HAPLOTYPES > 0)) {
    # Write the list of locus, individuals with paralogs
    if (is.null(whitelist.markers)) {
      write_tsv(x = paralogs.ind, path = "blacklist.loci.paralogs.ind.txt", col_names = TRUE)
      filename.paralogs.ind <- "blacklist.loci.paralogs.ind.txt"
    } else {
      write_tsv(x = paralogs.ind, path = "blacklist.loci.paralogs.ind.filtered.txt", col_names = TRUE)
      filename.paralogs.ind <- "blacklist.loci.paralogs.ind.filtered.txt"
    }
    
    paralogs.pop <- paralogs.ind %>%
      ungroup() %>%
      distinct(LOCUS, POP_ID) %>%
      arrange(as.numeric(LOCUS)) %>% 
      mutate(PARALOGS = rep("paralogs", times = n()))
    
    blacklist.loci.paralogs <- paralogs.ind %>%
      ungroup() %>%
      distinct(LOCUS) %>%
      arrange(as.numeric(LOCUS))
    
    # Write the unique list of paralogs blacklisted to a file
    if (is.null(whitelist.markers)) {
      write_tsv(
        x = blacklist.loci.paralogs,
        path = "blacklist.loci.paralogs.txt", col_names = TRUE
      )
      filename.paralogs <- "blacklist.loci.paralogs.txt"
    } else {
      write_tsv(
        x = blacklist.loci.paralogs, 
        path = "blacklist.loci.paralogs.filtered.txt", col_names = TRUE
      )
      filename.paralogs <- "blacklist.loci.paralogs.filtered.txt"
    }
  }
  
  # Haplo filtered paralogs ----------------------------------------------------
  # haplo.filtered.paralogs <- haplotype %>%
  #   filter(!LOCUS %in% blacklist.loci.paralogs$LOCUS)
  #or
  if (length(paralogs.ind$HAPLOTYPES > 0)) {
    message("Erasing genotypes with > 2 alleles")
    erase.paralogs <- paralogs.ind %>% 
      select(LOCUS, INDIVIDUALS) %>% 
      mutate(ERASE = rep("erase", n()))
    
    haplo.filtered.paralogs <- suppressWarnings(
      haplotype %>%
        full_join(erase.paralogs, by = c("LOCUS", "INDIVIDUALS")) %>%
        mutate(
          ERASE = stri_replace_na(str = ERASE, replacement = "ok"),
          HAPLOTYPES = ifelse(ERASE == "erase", NA, HAPLOTYPES)
        ) %>% 
        select(-ERASE)
    )
  } else {
    haplo.filtered.paralogs <- haplotype
  }
  
  # Haplo filtered for consensus -----------------------------------------------
  if (length(consensus.pop$CONSENSUS) > 0) {
    haplo.filtered.consensus <- haplotype %>%
      filter(!LOCUS %in% consensus.pop$LOCUS)
    
    # Haplo filtered for consensus and paralogs
    haplo.filtered.consensus.paralogs <- haplo.filtered.paralogs %>%
      filter(!LOCUS %in% consensus.pop$LOCUS)
  } else {
    haplo.filtered.consensus <- haplotype
    haplo.filtered.consensus.paralogs <- haplo.filtered.paralogs
  }
  
  # Summary dataframe by individual---------------------------------------------
  message("Genome-Wide Identity-By-Descent calculations (FH)...")
  
  summary.ind <- haplo.filtered.consensus.paralogs %>%
    mutate(ALLELES_COUNT = stri_count_fixed(HAPLOTYPES, "/")) %>% 
    mutate(
      IND_LEVEL_POLYMORPHISM = if_else(ALLELES_COUNT == 1, "het", "hom", missing = "missing")
    ) %>% 
    group_by(INDIVIDUALS) %>%
    summarise(
      HOM = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "hom"]),
      HET = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "het"]),
      MISSING = length(IND_LEVEL_POLYMORPHISM[IND_LEVEL_POLYMORPHISM == "missing"]),
      N_GENOT = HOM + HET,
      HOM_O = HOM/N_GENOT,
      HET_O = HET/N_GENOT
    ) %>% 
    left_join(strata.df, by = "INDIVIDUALS") %>% 
    arrange(POP_ID, INDIVIDUALS)
  
  freq.alleles.loci.pop <- suppressWarnings(
    haplo.filtered.consensus.paralogs %>% 
      filter(!is.na(HAPLOTYPES)) %>% 
      group_by(LOCUS, POP_ID) %>%
      mutate(DIPLO = length(INDIVIDUALS) * 2) %>% 
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
        sep = "/", extra = "drop", remove = F
      ) %>%
      mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
      select(-HAPLOTYPES, -INDIVIDUALS) %>% 
      tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID, DIPLO)) %>%
      group_by(LOCUS, POP_ID, ALLELES) %>% 
      summarise(
        FREQ_ALLELES = length(ALLELES)/mean(DIPLO),
        HOM_E = FREQ_ALLELES * FREQ_ALLELES
      ) %>% 
      select(-FREQ_ALLELES)
  )
  
  freq.loci.pop <- freq.alleles.loci.pop %>% 
    group_by(LOCUS, POP_ID) %>%
    summarise(HOM_E = sum(HOM_E)) 
  
  freq.pop <- freq.loci.pop %>% 
    group_by(POP_ID) %>%
    summarise(HOM_E = mean(HOM_E))
  
  
  # IBDg with FH ---------------------------------------------------------------
  fh.i <- summary.ind %>% 
    full_join(freq.pop, by = "POP_ID") %>% 
    mutate(FH = ((HOM_O - HOM_E)/(N_GENOT - HOM_E)))
  
  fh.i.out <- fh.i %>% select(INDIVIDUALS, POP_ID, FH)
  
  fh.pop <- fh.i %>% 
    group_by(POP_ID) %>% 
    summarise(
      HOM_O = round(mean(HOM_O), 6),
      HOM_E = round(mean(HOM_E), 6),
      HET_O = round(mean(1 - HOM_O), 6),
      HET_E = round(mean(1 - HOM_E), 6),
      FIS = ifelse(HET_O == 0, 0, round(((HET_E - HET_O) / HET_E), 6)),
      FH = mean(FH)
    )
  
  fh.tot <- fh.i %>% 
    summarise(
      HOM_O = round(mean(HOM_O), 6),
      HOM_E = round(mean(HOM_E), 6),
      HET_O = round(mean(1 - HOM_O), 6),
      HET_E = round(mean(1 - HOM_E), 6),
      FIS = ifelse(HET_O == 0, 0, round(((HET_E - HET_O) / HET_E), 6)),
      FH = mean(FH)
    )
  
  fh.tot <- data_frame(POP_ID = "OVERALL") %>%
    bind_cols(fh.tot)
  
  fh.res <- suppressWarnings(bind_rows(fh.pop, fh.tot) %>% select(-POP_ID))
  
  # unused arguments  
  freq.pop <- summary.ind <- fh.tot <- fh.pop <- NULL
  
  # Nei & Li 1979 Nucleotide Diversity -----------------------------------------
  message("Nucleotide diversity (Pi) calculations")
  
  pi.data <- suppressWarnings(
    haplo.filtered.paralogs %>%
      filter(!is.na(HAPLOTYPES)) %>% 
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
        sep = "/", extra = "drop", remove = TRUE
      ) %>%
      mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2))
  )
  # Pi: by individuals----------------------------------------------------------
  message("Pi calculations by individuals...")
  
  pi.data.i <- pi.data %>%
    mutate(
      PI = (stringdist::stringdist(a = ALLELE1, b = ALLELE2, method = "hamming"))/read.length
    ) %>% 
    group_by(INDIVIDUALS) %>% 
    summarise(PI = mean(PI))
  
  pi.i.out <- pi.data.i %>% select(INDIVIDUALS, PI)
  
  fh.pi.individuals <- full_join(fh.i.out, pi.i.out, by = "INDIVIDUALS") 
  
  # Pi: by pop------------------------------------------------------------------
  message("Pi calculations by populations, take a break...")
  pi.data.pop <- pi.data %>% 
    tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, INDIVIDUALS, POP_ID))
  
  df.split.pop <- split(x = pi.data.pop, f = pi.data.pop$POP_ID) # slip data frame by population
  pop.list <- names(df.split.pop) # list the pop
  pi.res <- list() # create empty list
  
  # Pi function ----------------------------------------------------------------
  pi <- function(y, read.length) {
    
    if (length(unique(y)) <= 1) {
      PI <- 0
      PI <- data.frame(PI)
    } else {
      
      #1 Get all pairwise comparison
      allele_pairwise <- utils::combn(unique(y), 2)
      # allele_pairwise
      
      #2 Calculate pairwise nucleotide mismatches
      pairwise_mismatches <- apply(allele_pairwise, 2, function(z) {
        stringdist::stringdist(a = z[1], b = z[2], method = "hamming")
      })
      # pairwise_mismatches
      
      #3 Calculate allele frequency
      allele_freq <- table(y)/length(y)
      # allele_freq
      
      #4 Calculate nucleotide diversity from pairwise mismatches and allele frequency
      pi.prep <- apply(allele_pairwise, 2, function(y) allele_freq[y[1]] * allele_freq[y[2]])
      # pi.prep
      # read.length <- 80
      PI <- sum(pi.prep * pairwise_mismatches) / read.length
      PI <- data.frame(PI)
    }
    return(PI)
  }
  pi_pop <- function(pop.list, read.length) {
    # message of progress for pi calculation by population
    # pop.list <- "SKY" # test
    # pop.pi.calculations <- paste("Pi calculations for pop ", pop.list, sep = "")
    # message(pop.pi.calculations)
    
    pop.data <- df.split.pop[[pop.list]]
    pi.pop.data <- pop.data %>% 
      group_by(LOCUS, POP_ID) %>% 
      do(., pi(y = .$ALLELES, read.length = read.length))
    pi.res[[pop.list]] <- pi.pop.data
  }
  pi.res <- parallel::mclapply(
    X = pop.list, 
    FUN = pi_pop, 
    mc.preschedule = FALSE, 
    mc.silent = FALSE, 
    mc.cleanup = TRUE,
    mc.cores = parallel.core,
    read.length = read.length
  )
  # pi.res <- as.data.frame(bind_rows(pi.res))
  pi.res <- suppressWarnings(
    bind_rows(pi.res) %>%
      group_by(POP_ID) %>% 
      summarise(PI_NEI = mean(PI))
  )
  
  df.split.pop <- NULL
  pop.list <- NULL
  
  # Pi: overall  ---------------------------------------------------------------
  message("Calculating Pi overall")
  pi.overall <- pi.data.pop %>% 
    group_by(LOCUS) %>% 
    do(., pi(y = .$ALLELES, read.length = read.length)) %>% 
    ungroup() %>%
    summarise(PI_NEI = mean(PI))
  
  pi.tot <- data_frame(POP_ID = "OVERALL") %>% bind_cols(pi.overall)
  
  # Combine the pop and overall data
  pi.res <- suppressWarnings(bind_rows(pi.res, pi.tot) %>% select(-POP_ID))
  
  # Summary dataframe by pop ---------------------------------------------------
  message("Working on the summary table")
  
  summary.prep <- suppressWarnings(
    haplo.filtered.consensus %>% 
      filter(!is.na(HAPLOTYPES)) %>%
      select(-INDIVIDUALS) %>%
      tidyr::separate(
        col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
        sep = "/", extra = "drop", remove = TRUE
      ) %>%
      mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
      tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID))
  )
  
  summary.pop <- summary.prep %>%
    group_by(LOCUS, POP_ID) %>%
    distinct(ALLELES, .keep_all = TRUE) %>% 
    summarise(ALLELES_COUNT = length(ALLELES)) %>%
    group_by(POP_ID) %>% 
    summarise(
      MONOMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT == 1]),
      POLYMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT >= 2])
    ) %>%
    full_join(
      consensus.pop %>%
        group_by(POP_ID) %>%
        summarise(CONSENSUS = n_distinct(LOCUS)),
      by = "POP_ID"
    ) %>%
    group_by(POP_ID) %>%
    mutate(TOTAL = MONOMORPHIC + POLYMORPHIC + CONSENSUS) %>%
    full_join(
      paralogs.pop %>%
        group_by(POP_ID) %>%
        summarise(PARALOGS = n_distinct(LOCUS)),
      by = "POP_ID"
    ) %>%
    mutate(
      MONOMORPHIC_PROP = round(MONOMORPHIC/TOTAL, 4),
      POLYMORPHIC_PROP = round(POLYMORPHIC/TOTAL, 4),
      CONSENSUS_PROP = round(CONSENSUS/TOTAL, 4),
      PARALOG_PROP = round(PARALOGS/TOTAL, 4)
    )
  
  
  total <- summary.prep %>%
    group_by(LOCUS) %>%
    distinct(ALLELES, .keep_all = TRUE) %>% 
    summarise(ALLELES_COUNT = length(ALLELES)) %>%
    summarise(
      MONOMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT == 1]),
      POLYMORPHIC = length(ALLELES_COUNT[ALLELES_COUNT >= 2])
    ) %>%
    bind_cols(
      blacklist.loci.consensus %>%
        ungroup %>%
        summarise(CONSENSUS = n()) %>% 
        select(CONSENSUS)
    ) %>%
    mutate(TOTAL = MONOMORPHIC + POLYMORPHIC + CONSENSUS) %>%
    bind_cols(
      blacklist.loci.paralogs %>%
        ungroup %>%
        summarise(PARALOGS = n()) %>% 
        select(PARALOGS)
    ) %>%
    mutate(
      MONOMORPHIC_PROP = round(MONOMORPHIC/TOTAL, 4),
      POLYMORPHIC_PROP = round(POLYMORPHIC/TOTAL, 4),
      CONSENSUS_PROP = round(CONSENSUS/TOTAL, 4),
      PARALOG_PROP = round(PARALOGS/TOTAL, 4)
    )
  
  total.res <- data_frame(POP_ID = "OVERALL") %>%
    bind_cols(total)
  
  summary <- suppressWarnings(bind_rows(summary.pop, total.res) %>% ungroup() %>% select(-POP_ID))
  summary <- bind_cols(sample.number, summary, fh.res, pi.res)
  
  if (is.null(whitelist.markers)) {
    write.table(summary, "haplotype.catalog.loci.summary.pop.tsv", 
                sep = "\t", row.names = F, col.names = T, quote = F)
    filename.sum <- "haplotype.catalog.loci.summary.pop.tsv"
  } else {
    write.table(summary, "haplotype.catalog.loci.summary.pop.filtered.tsv", 
                sep = "\t", row.names = F, col.names = T, quote = F)
    filename.sum <- "haplotype.catalog.loci.summary.pop.filtered.tsv"
  } 
  
  
  # Figures --------------------------------------------------------------------
  fh.pi <- pi.data.i %>% 
    full_join(
      fh.i %>% select(INDIVIDUALS, POP_ID, FH)
      , by = "INDIVIDUALS")
  
  scatter.plot <- ggplot(fh.pi, aes(x = FH, y = PI)) + 
    geom_point(aes(colour = POP_ID)) +
    stat_smooth(method = stats::lm, level = 0.95, fullrange = FALSE, na.rm = TRUE) +
    labs(x = "Individual IBDg (FH)") + 
    labs(y = "Individual nucleotide diversity (Pi)") +
    theme(
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"), 
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.y = element_text(angle = 0, size = 10, family = "Helvetica", face = "bold"), 
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
  
  
  boxplot.pi <-   ggplot(fh.pi, aes(x = factor(POP_ID), y = PI, na.rm = T)) +
    geom_violin(trim = F) +
    geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    labs(x = "Sampling sites") +
    labs(y = "Individual nucleotide diversity (Pi)") +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
  
  
  
  boxplot.fh <-   ggplot(fh.pi, aes(x = factor(POP_ID), y = FH, na.rm = T)) +
    geom_violin(trim = F) +
    geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    labs(x = "Sampling sites") +
    labs(y = "Individual IBDg (FH)") +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stri_paste("The number of loci in the catalog = ", n_distinct(haplotype$LOCUS)))
  message(stri_paste("The number of putative paralogs loci in the catalog (> 2 alleles) = ", n_distinct(paralogs.pop$LOCUS)))
  message(stri_paste("The number of loci in the catalog with consensus alleles = ", n_distinct(consensus.pop$LOCUS)))
  message(stri_paste("3 files were written in this directory: ", getwd()))
  message(filename.paralogs)
  message(filename.sum)
  message("blacklist.loci.consensus.txt")
  cat("#######################################################################\n")
  # Results
  results <- list()
  results$summary <- summary
  results$paralogs.pop <- paralogs.pop
  results$paralogs.loci <- blacklist.loci.paralogs
  results$consensus.pop <- consensus.pop
  results$consensus.loci <- blacklist.loci.consensus
  results$fh.pi.individuals <- fh.pi.individuals
  results$scatter.plot <- scatter.plot
  results$boxplot.pi <- boxplot.pi
  results$boxplot.fh <- boxplot.fh
  cat("############################## completed ##############################\n")
  return(results)
}
