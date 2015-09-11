## Summary and tables

#' @title Haplotypes file summary
#' @description STACKS batch_x.haplotypes.tsv file summary.
#' The output of the function is a summary table for populations with putative
#' paralogs, consensus, monomorphic and polymorphic loci. 
#' The haplotypes statistics for the observed and expected homozygosity and
#' heterozygosity. Wright’s inbreeding coefficient (Fis), and a proxy measure of 
#' the realized proportion of the genome that is identical by descent (IBDG),
#' the FH measure based on the excess in the observed number of homozygous
#' genotypes within an individual relative to the mean number of homozygous 
#' genotypes expected under random mating 
#' (Keller et al., 2011; Kardos et al., 2015). 
#' The nucleotide diversity (Pi) is also given. Pi measured here consider the 
#' consensus loci in the catalog (no variation between population sequences).
#' @param haplotypes.file The 'batch_x.haplotypes.tsv' created by STACKS.
#' @param whitelist.loci (optional) A whitelist of loci with a column name
#' 'LOCUS'. The whitelist is located in the global environment or in the directory (e.g. "whitelist.txt"). 
#' If a whitelist is used, the files written to the directory will have 
#' '.filtered' in the filename. 
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the global environment 
#' or in the directory (with "blacklist.txt").
#' @param pop.id.start The start of your population id 
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id 
#' in the name of your individual sample.
#' @param pop.levels An optional character string with your populations ordered.
#' @param read.length The length in nucleotide of your reads (e.g. 80 or 100).
#' @import stringdist
#' @return The function returns a list with the summary, the paralogs and
#' consensus loci by populations and unique loci and 3 plots (use $ to access each 
#' objects). 
#' Write 3 files in the working directory:
#' blacklist of unique putative paralogs and unique consensus loci 
#' and a summary of the haplotypes file by population.
#' @details If the object for the function is 'haplotype.file.summary' then:
#' 
#' haplo.summary <- haplotype.file.summary$summary
#'
#' paralogs.pop <- haplotype.file.summary$paralogs.pop
#' 
#' paralogs.loci <- haplotype.file.summary$paralogs.loci
#' 
#' consensus.pop <- haplotype.file.summary$consensus.pop
#' 
#' consensus.loci <- haplotype.file.summary$consensus.loci
#' 
#' scatter.plot <- haplotype.file.summary$scatter.plot
#' 
#' boxplot.pi <- haplotype.file.summary$boxplot.pi
#' 
#' boxplot.fh <- haplotype.file.summary$boxplot.fh
#' @rdname summary_haplotypes_v2
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


summary_haplotypes <- function(haplotypes.file, 
                               whitelist.loci = NULL, 
                               blacklist.id = NULL, 
                               pop.id.start, 
                               pop.id.end, 
                               pop.levels, 
                               read.length) {
  
  POP_ID <- NULL
  POLYMORPHISM <- NULL
  POLYMORPHISM_MAX <- NULL
  PARALOGS <- NULL
  CONSENSUS <- NULL
  CONSENSUS_MAX <- NULL
  ALLELES_COUNT <- NULL
  ALLELES_COUNT_SUM <- NULL
  POP_LEVEL_POLYMORPHISM <- NULL
  MONOMORPHIC <- NULL
  POLYMORPHIC <- NULL
  TOTAL <- NULL
  IND_LEVEL_POLYMORPHISM <- NULL
  N_GENOT <- NULL
  ALLELE_GROUP <- NULL
  ALLELES <- NULL
  HOM_O <- NULL
  HOM_E <- NULL
  HET_O <- NULL
  HET_E <- NULL
  FH <- NULL
  PI <- NULL
  HOM <- NULL
  HET <- NULL
  DIPLO <- NULL
  FREQ_ALLELES <- NULL
  
  # Import haplotype file ------------------------------------------------------
  haplotype <- read_tsv(haplotypes.file, col_names = T) %>%
    rename(LOCUS =`Catalog ID`) %>%
    tidyr::gather(INDIVIDUALS, HAPLOTYPES, -c(LOCUS, Cnt)) %>%
    mutate(
      POP_ID = str_sub(INDIVIDUALS, pop.id.start, pop.id.end),
      POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
    )
  
  # Whitelist loci -------------------------------------------------------------
  if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "TRUE") {
    message("Whitelist of loci: from the directory")
    whitelist <- read_tsv(whitelist.loci, col_names = T)
  } else if (missing(whitelist.loci) == "FALSE" & is.vector(whitelist.loci) == "FALSE") {
    message("Whitelist of loci: from your global environment")
    whitelist <- whitelist.loci
  } else {
    message("Whitelist of loci: no")
    whitelist <- NULL
  }
  
  # Blacklisted individuals and consensus loci ---------------------------------
  if (missing(blacklist.id) == "FALSE" & is.vector(blacklist.id) == "TRUE") {
    message("Blacklisted id: file from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)
  } else if (missing(blacklist.id) == "FALSE" & is.vector(blacklist.id) == "FALSE") {
    message("Blacklisted id: object from your global environment")
    blacklist.id <- blacklist.id
  } else {
    message("Blacklisted id: no")
    blacklist.id <- NULL
  }
  
  # Locus with concensus alleles -----------------------------------------------
  message("Looking for consensus...")
  
  if (is.null(blacklist.id) == TRUE) {
    # Without blacklisted individuals
    consensus.pop <- haplotype %>%
      mutate(CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus")) %>%
      group_by(LOCUS, POP_ID) %>%
      summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
      filter(CONSENSUS_MAX > 0) %>%
      mutate(CONSENSUS = rep("consensus", times = n())) %>%
      select(LOCUS, POP_ID, CONSENSUS)
    
  } else {
    # With a blacklist of individuals
    consensus.pop <- haplotype %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      mutate(CONSENSUS = stri_count_fixed(HAPLOTYPES, "consensus")) %>%
      group_by(LOCUS, POP_ID) %>%
      summarise(CONSENSUS_MAX = max(CONSENSUS)) %>%
      filter(CONSENSUS_MAX > 0) %>%
      mutate(CONSENSUS = rep("consensus", times = n())) %>%
      select(LOCUS, POP_ID, CONSENSUS)
  }  
# Create a list of consensus loci
  blacklist.loci.consensus <- consensus.pop %>%
    group_by(LOCUS) %>%
    select (LOCUS) %>%
    distinct(LOCUS) %>%
    arrange(LOCUS)
  
  write.table(blacklist.loci.consensus, 
              "blacklist.loci.consensus.txt",
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  consensus.pop  
  
  if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == TRUE) {
    
    # Combination 1: No whitelist and No blacklist -----------------------------
    haplotype <- haplotype    
    
  } else if (is.null(whitelist.loci) == FALSE & is.null(blacklist.id) == TRUE) {
    
    # Combination 2: Using whitelist, but No blacklist -------------------------
    haplotype <- haplotype %>%
      semi_join(whitelist, by = "LOCUS") %>% 
      arrange(LOCUS)
    
  } else if (is.null(whitelist.loci) == TRUE & is.null(blacklist.id) == FALSE) {
    
    # Combination 3: Using a blacklist of id, but No whitelist -----------------
    haplotype <- haplotype %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(LOCUS)
    
  } else {
    # Combination 4: Using a whitelist and blacklist---------------------------
    haplotype <- haplotype %>%
      semi_join(whitelist, by = "LOCUS") %>%
      mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>%
      anti_join(blacklist.id, by = "INDIVIDUALS") %>%
      arrange(LOCUS)
  }
  
  # Individuals per pop
  ind.pop <- haplotype %>%
    distinct(INDIVIDUALS) %>% 
    group_by(POP_ID) %>%
    tally() %>% 
    rename(NUM = n)
   
  ind.tot <-haplotype %>%
    distinct(INDIVIDUALS) %>% 
    tally() %>% 
    rename(NUM = n)
  
  ind.tot <- data_frame(POP_ID="OVERALL") %>%
    bind_cols(ind.tot)
  
  sample.number <- bind_rows(ind.pop, ind.tot)
  
  # Paralogs... Locus with > 2 alleles by individuals --------------------------
  # Create a blacklist of catalog loci with paralogs
  message("Looking for paralogs...")
  paralogs.pop <- haplotype %>%
    mutate(POLYMORPHISM = stri_count_fixed(HAPLOTYPES, "/")) %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(POLYMORPHISM_MAX = max(POLYMORPHISM)) %>%
    filter(POLYMORPHISM_MAX > 1) %>%
    mutate(PARALOGS = rep("paralogs", times = n())) %>%
    select(LOCUS, POP_ID, PARALOGS)
  
  blacklist.loci.paralogs <- paralogs.pop %>%
    group_by(LOCUS) %>%
    select (LOCUS) %>%
    distinct(LOCUS) %>%
    arrange(LOCUS)
  
  # Write the paralogs blacklisted to a file
   if (missing(whitelist.loci) == "FALSE") {
    write.table(blacklist.loci.paralogs,
                "blacklist.loci.paralogs.filtered.txt",
                sep = "\t", row.names = F, col.names = T, quote = F)
    filename.paralogs <- "blacklist.loci.paralogs.filtered.txt"
  } else {
    write.table(blacklist.loci.paralogs,
                "blacklist.loci.paralogs.txt",
                sep = "\t", row.names = F, col.names = T, quote = F)
    filename.paralogs <- "blacklist.loci.paralogs.txt"
  }

  paralogs.pop
  
  # dump unused object
  whitelist <- NULL
  blacklist.id <- NULL
  
  # Haplo filtered paralogs
  haplo.filtered.paralogs <- haplotype %>%
    filter(!LOCUS %in% blacklist.loci.paralogs$LOCUS)
  
  # Haplo filtered for consensus -----------------------------------------------
  haplo.filtered.consensus <- haplotype %>%
    filter(!LOCUS %in% consensus.pop$LOCUS)
  
  # Haplo filtered for consensus and paralogs-----------------------------------
  haplo.filtered.consensus.paralogs <- haplotype %>%
    filter(!LOCUS %in% consensus.pop$LOCUS)%>%
    filter(!LOCUS %in% blacklist.loci.paralogs$LOCUS)
  
  # Summary dataframe by individual---------------------------------------------
  message("Genome-Wide Identity-By-Descent calculations (FH)...")
  
  
  summary.ind <- haplo.filtered.consensus.paralogs %>%
    mutate(ALLELES_COUNT = stri_count_fixed(HAPLOTYPES, "/")) %>% 
    mutate(
      IND_LEVEL_POLYMORPHISM = ifelse(HAPLOTYPES == "-", "missing",
                                      ifelse(ALLELES_COUNT == 0 & HAPLOTYPES != "-", "hom", "het"))
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
    mutate(POP_ID = factor(str_sub(INDIVIDUALS, pop.id.start, pop.id.end), 
                           levels = pop.levels, ordered = T)) %>% 
    arrange(POP_ID, INDIVIDUALS)
  
  freq.alleles.loci.pop <- haplo.filtered.consensus.paralogs %>% 
    filter(HAPLOTYPES != "-") %>% 
    group_by(LOCUS, POP_ID) %>%
    mutate(DIPLO= length(INDIVIDUALS) *2) %>% 
    tidyr::separate(
      col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
      sep = "/", extra = "drop", remove = F
    ) %>%
    mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
    select(-Cnt, -HAPLOTYPES, -INDIVIDUALS) %>% 
    tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID, DIPLO)) %>%
    group_by(LOCUS, POP_ID, ALLELES) %>% 
    summarise(
      FREQ_ALLELES = length(ALLELES)/mean(DIPLO),
      HOM_E = FREQ_ALLELES * FREQ_ALLELES
    ) %>% 
    select(-FREQ_ALLELES)
  
  freq.loci.pop<- freq.alleles.loci.pop %>% 
    group_by(LOCUS, POP_ID) %>%
    summarise(HOM_E = sum(HOM_E)) 
  
  freq.pop <- freq.loci.pop %>% 
    group_by(POP_ID) %>%
    summarise(HOM_E = mean(HOM_E))
  
  
  # IBDg with FH ---------------------------------------------------------------
  
  fh.i <- summary.ind %>% 
    full_join(freq.pop, by = "POP_ID") %>% 
    mutate(FH = ((HOM_O - HOM_E)/(N_GENOT - HOM_E)))
  
  fh.pop <- fh.i %>% 
    group_by(POP_ID) %>% 
    summarise(
      HOM_O = round(mean(HOM_O), 6),
      HOM_E = round(mean(HOM_E), 6),
      HET_O = round(mean(1-HOM_O), 6),
      HET_E = round(mean(1-HOM_E), 6),
      FIS = ifelse(HET_O == 0, 0, round (((HET_E - HET_O) / HET_E), 6)),
      FH = mean(FH)
    )
  
  fh.tot <- fh.i %>% 
    summarise(
      HOM_O = round(mean(HOM_O), 6),
      HOM_E = round(mean(HOM_E), 6),
      HET_O = round(mean(1-HOM_O), 6),
      HET_E = round(mean(1-HOM_E), 6),
      FIS = ifelse(HET_O == 0, 0, round (((HET_E - HET_O) / HET_E), 6)),
      FH = mean(FH)
    )
  
  fh.tot <- data_frame(POP_ID="OVERALL") %>%
    bind_cols(fh.tot)
  fh.res <- bind_rows(fh.pop, fh.tot) %>% select(-POP_ID)
  
  # fh.i <- NULL
  freq.pop <- NULL
  summary.ind <- NULL
  fh.tot <- NULL
  fh.pop <- NULL
  
  
  # Nei & Li 1979 Nucleotide Diversity -----------------------------------------
  message("Nucleotide diversity (Pi) calculations")
  
  pi.data <- haplo.filtered.paralogs %>%
    select(-Cnt) %>% 
    filter(HAPLOTYPES != "-") %>% 
    tidyr::separate(
      col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
      sep = "/", extra = "drop", remove = T
    ) %>%
    mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2))
  
  # Pi: by individuals----------------------------------------------------------
  message("Pi calculations by individuals...")
  
  pi.data.i <- pi.data %>%
    mutate(
      PI = (stringdist::stringdist(a = ALLELE1, b = ALLELE2, method = "hamming"))/read.length
    ) %>% 
    group_by(INDIVIDUALS) %>% 
    summarise(PI = mean(PI))
  
  
  # Pi: by pop------------------------------------------------------------------
  message("Pi calculations by populations, take a break...")
  
  pi.data.pop <- pi.data %>% 
    tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, INDIVIDUALS, POP_ID))
  
  df.split.pop <- split(x = pi.data.pop, f = pi.data.pop$POP_ID) # slip data frame by population
  pop.list <- names(df.split.pop) # list the pop
  pi.res <-list() # create empty list
  
  for (i in pop.list) {
    # message of progress for pi calculation by population
    pop.pi.calculations <- paste("Pi calculations for pop ", i, sep = "")
    message(pop.pi.calculations)
    
    pi <- function(y, read.length) {
      
      if(length(unique(y)) <= 1){
        PI <- 0
        PI <- data.frame(PI)
      } else{
        
        #1 Get all pairwise comparison
        allele_pairwise <- combn(unique(y), 2)
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
        PI <- sum(pi.prep*pairwise_mismatches)/read.length
        PI <- data.frame(PI)
      }
      return(PI)
    }
    pop.data <- df.split.pop[[i]]
    pi.pop.data <- pop.data %>% 
      group_by(LOCUS, POP_ID) %>% 
      do(., pi(y = .$ALLELES, read.length = read.length))
    pi.res[[i]] <- pi.pop.data
  }
  pi.res <- as.data.frame(bind_rows(pi.res))
  pi.res <- pi.res %>% group_by(POP_ID) %>% summarise(PI_NEI = mean(PI))
  
  df.split.pop <- NULL
  pop.list <- NULL
  
  
  # Pi: overall  ---------------------------------------------------------------
  message("Calculating Pi overall")
  pi.overall <- pi.data.pop %>% 
    group_by(LOCUS) %>% 
    do(., pi(y = .$ALLELES, read.length = read.length)) %>% 
    ungroup() %>%
    summarise(PI_NEI = mean(PI))
  
  pi.tot <- data_frame(POP_ID="OVERALL") %>% bind_cols(pi.overall)
  
  # Combine the pop and overall data
  pi.res <- bind_rows(pi.res, pi.tot) %>% select(-POP_ID)
  
  # Summary dataframe by pop ---------------------------------------------------
  message("Working on the summary table")
  
  summary.prep <- haplo.filtered.consensus %>% 
    filter(HAPLOTYPES != "-") %>%
    select(-Cnt, -INDIVIDUALS) %>%
    tidyr::separate(
      col = HAPLOTYPES, into = c("ALLELE1", "ALLELE2"), 
      sep = "/", extra = "drop", remove = T
    ) %>%
    mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2)) %>%
    tidyr::gather(ALLELE_GROUP, ALLELES, -c(LOCUS, POP_ID))
  
  summary.pop <- summary.prep %>%
    group_by(LOCUS, POP_ID) %>%
    distinct(ALLELES) %>% 
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
    distinct(ALLELES) %>% 
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
  
  total.res <- data_frame(POP_ID="OVERALL") %>%
    bind_cols(total)
  
  summary <- bind_rows(summary.pop, total.res) %>% select(-POP_ID)
  summary <- bind_cols(sample.number, summary, fh.res, pi.res)
  
  if (missing(whitelist.loci) == "FALSE") {
    write.table(summary, "haplotype.catalog.loci.summary.pop.filtered.tsv", 
                sep = "\t", row.names = F, col.names = T, quote = F)
    filename.sum <- "haplotype.catalog.loci.summary.pop.filtered.tsv"
  } else {
    write.table(summary, "haplotype.catalog.loci.summary.pop.tsv", 
                sep = "\t", row.names = F, col.names = T, quote = F)
    filename.sum <- "haplotype.catalog.loci.summary.pop.tsv"
  }
  
  
  # Figures --------------------------------------------------------------------
  fh.pi <- pi.data.i %>% 
    full_join(
      fh.i %>% select(INDIVIDUALS, POP_ID, FH)
      , by = "INDIVIDUALS")
  
  scatter.plot <- ggplot(fh.pi, aes(x = FH, y = PI)) + 
    geom_point(aes(colour = POP_ID)) +
    stat_smooth(method = lm, level = 0.95, fullrange = F, na.rm = T)+
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
  
  
  boxplot.pi <-   ggplot(fh.pi, aes(x = factor(POP_ID), y = PI, na.rm = T))+
    geom_violin(trim = F)+
    geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
    stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
    labs(x = "Sampling sites")+
    labs(y = "Individual nucleotide diversity (Pi)")+
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
  
  
  
  boxplot.fh <-   ggplot(fh.pi, aes(x = factor(POP_ID), y = FH, na.rm = T))+
    geom_violin(trim = F)+
    geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA)+
    stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white")+
    labs(x = "Sampling sites")+
    labs(y = "Individual IBDg (FH)")+
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
    )
  
  invisible(cat(sprintf(
    "The number of loci in the catalog = %s LOCI
The number of putative paralogs loci in the catalog (> 2 alleles) = %s LOCI
The number of loci in the catalog with consensus alleles = %s LOCI
3 files were written in this directory: %s
1. %s
2. blacklist.loci.consensus.txt
3. %s", 
    n_distinct(haplotype$LOCUS),
    n_distinct(paralogs.pop$LOCUS),
    n_distinct(consensus.pop$LOCUS),
    getwd(),
    filename.paralogs,
    filename.sum
  )))
  
  # Results
  results <- list()
  results$summary <- summary
  results$paralogs.pop <- paralogs.pop
  results$paralogs.loci <- blacklist.loci.paralogs
  results$consensus.pop <- consensus.pop
  results$consensus.loci <- blacklist.loci.consensus
  results$scatter.plot <- scatter.plot
  results$boxplot.pi <- boxplot.pi
  results$boxplot.fh <- boxplot.fh
  
  return(results)
  
}








#' @title Import and summarise the batch_x.hapstats.tsv file
#' @description Import and summarise the batch_x.hapstats.tsv file.
#' Necessary preparation for density distribution and box plot figures.
#' @param data The 'batch_x.hapstats.tsv' created by STACKS.
#' @param pop.num The number of populations analysed.
#' @param pop.col.types \code{"integer"} or \code{"character"} used in STACKS populations module?
#' @param pop.integer.equi When Integer was used for your population id,
#' give the character equivalence
#' @param pop.levels A character string with your populations in order.
#' @rdname summary_hapstats
#' @export 

summary_hapstats <- function(data, pop.num, pop.col.types, pop.integer.equi, pop.levels) {
  
  POP_ID <- NULL
  
  skip.lines <- pop.num + 1
  
  if(pop.col.types == "integer"){
    col.types = "iiciiiiddddddc"
  } 
  if(pop.col.types == "character") {
    col.types = "iiciciiddddddc"
  } else {
    col.types = NULL
  }
  hapstats <- read_tsv(data,
                       na = "NA",
                       skip = skip.lines,
                       progress = interactive(),
                       col_names = c("BATCH_ID", "LOCUS", "CHR", "BP", "POP_ID", "N", "HAPLOTYPE_CNT", "GENE_DIVERSITY", "SMOOTHED_GENE_DIVERSITY", "SMOOTHED_GENE_DIVERSITY_PVALUE", "HAPLOTYPE_DIVERSITY", "SMOOTHED_HAPLOTYPE_DIVERSITY", "SMOOTHED_HAPLOTYPE_DIVERSITY_PVALUE", "HAPLOTYPES"),
                       col_types = col.types) %>%
    mutate (
      POP_ID = stri_replace_all_fixed(POP_ID, seq(from = 1, to = pop.num, by = 1), pop.integer.equi, vectorize_all=F),
      POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)
    ) %>%
    arrange(LOCUS, POP_ID)
}




## VCF
#' @title Summary statistics of a tidy VCF by population and markers
#' @description Summarise and prepare the tidy VCF. 
#' Summary, by population and markers (SNP), of frequency of the REF 
#' and the ALT alleles, the observed and the expected heterozygosity 
#' and the inbreeding coefficient. The Global MAF of Loci, 
#' with STACKS GBS/RAD loci = read or de novo haplotypes, 
#' is included and repeated over SNP.
#' @param filename (optional) Name of the file written to the working directory.
#' @param data The tidy VCF file created with \link{read_stacks_vcf}.
#' @rdname summary_stats_vcf_tidy
#' @export

summary_stats_vcf_tidy <- function(data, filename) {
  
  
  GT <- NULL
  GL <- NULL
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  N <- NULL
  HET_O <- NULL
  HOM_O <- NULL
  HET_E <- NULL
  HOM_E <- NULL
  FREQ_ALT <- NULL
  FREQ_REF <- NULL
  GLOBAL_MAF <- NULL
  PP <- NULL
  PQ <- NULL
  QQ <- NULL
  
  
  
  
  vcf.summary <- data %>%
    filter(GT != "./.") %>%
    group_by(LOCUS, POS, POP_ID) %>%
    summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT[GT == "0/0"])),
      PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
      QQ = as.numeric(length(GT[GT == "1/1"]))
    ) %>%
    mutate(
      FREQ_REF = ((PP*2) + PQ)/(2*N),
      FREQ_ALT = ((QQ*2) + PQ)/(2*N),
      HET_O = PQ/N,
      HET_E = 2 * FREQ_REF * FREQ_ALT,
      FIS = ifelse(HET_O == 0, 0, round (((HET_E - HET_O) / HET_E), 6))
    )
  
  global.maf <- vcf.summary %>%
    group_by(LOCUS, POS) %>%
    summarise_each_(funs(sum), vars = c("N", "PP", "PQ", "QQ")) %>%
    mutate(GLOBAL_MAF = (PQ + (2 * QQ)) / (2*N)) %>%
    select(LOCUS, POS, GLOBAL_MAF)
  
  vcf.prep <- global.maf %>%
    left_join(vcf.summary, by = c("LOCUS", "POS"))
  
  vcf.prep <- vcf.prep[c("LOCUS", "POS", "POP_ID", "N", "PP", "PQ", "QQ", "FREQ_REF", "FREQ_ALT", "GLOBAL_MAF", "HET_O", "HET_E", "FIS")]
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(vcf.prep, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected"
  }
  
  return(vcf.prep)
}

#' @title Summary statistics of a tidy VCF by population
#' @description Summarise the tidy VCF. 
#' The populations summary on :  frequency of the REF 
#' and the ALT alleles, the observed and the expected heterozygosity 
#' and the inbreeding coefficient. The Global MAF of Loci, 
#' with STACKS GBS/RAD loci = read or de novo haplotypes, 
#' is included and repeated over SNP.
#' @param filename (optional) Name of the file written to the working directory.
#' @param data The tidy VCF file created with read_stacks_vcf.
#' @rdname summary_stats_pop
#' @export

summary_stats_pop <- function(data, filename) {
  
  
  POP_ID <- NULL
  N <- NULL
  HET_O <- NULL
  HET_E <- NULL
  FREQ_REF <- NULL
  FIS <- NULL
  SNP <- NULL
  LOCUS <- NULL
  
  
  vcf.summary <- data %>%
    group_by(POP_ID) %>%
    summarise(
      SNP = length(unique(POS)),
      LOCUS = length(unique(LOCUS)),
      N = max(N, na.rm = TRUE),
      FREQ_REF = mean(FREQ_REF, na.rm = TRUE),
      HET_O = mean(HET_O, na.rm = TRUE),
      HET_E = mean(HET_E, na.rm = TRUE),
      FIS = mean(FIS, na.rm = TRUE)
    ) %>%
    select(POP_ID, N, SNP, LOCUS, FREQ_REF, HET_O, HET_E, FIS)  
  
  
  if (missing(filename) == "FALSE") {
    message("Saving the file in your working directory...")
    write_tsv(vcf.summary, filename, append = FALSE, col_names = TRUE)
    saving <- paste("Saving was selected, the filename:", filename, sep = " ")
  } else {
    saving <- "Saving was not selected"
  }
  
  return(vcf.summary)
}



## Coverage
#' @title Coverage summary
#' @description This function create a table summary of the important
#' coverage statistics from the tidy vcf created with read_stacks_vcf.
#' @param tidy.vcf.file The tidy VCF file created with read_stacks_vcf.
#' @param pop.levels Character string defining your ordered populations.
#' @param filename Name of the file saved to the working directory.
#' @details The tables contains summary statistics (mean, median, min, max)
#' of read, ref and alt allele coverage. To access 
#' the two tables, use $. The table that summarise by populations was created
#' using average nested: loci -> individuals -> populations.
#' The long format is used for creating figures.
#' @return A list with 2 tables: the long format of loci and populations
#' coverage statistics and the short format by populations.
#' The short-format is more user-friendly and
#' is written to the working directory.
#' @rdname summary_coverage
#' @export

summary_coverage <- function (tidy.vcf.file, pop.levels, filename) {
  
  POP_ID <- NULL
  READ_DEPTH <- NULL
  ALLELE_REF_DEPTH <- NULL
  ALLELE_ALT_DEPTH <- NULL
  INDIVIDUALS <- NULL
  
  
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, col_names = T, col_types = "iiiiccddcdccddddc")
    message("Using the file in your directory")
    
  } else {
    data = tidy.vcf.file
    message("Using the file from your global environment")
    
  }
  
  coverage.sum.loci <- data %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      READ_MEAN = mean(READ_DEPTH, na.rm = T),
      READ_MEDIAN = median(READ_DEPTH, na.rm = T),
      READ_MIN = min(READ_DEPTH, na.rm = T),
      READ_MAX = max(READ_DEPTH, na.rm = T),
      REF_MEAN = mean(ALLELE_REF_DEPTH, na.rm = T),
      REF_MEDIAN = median(ALLELE_REF_DEPTH, na.rm = T),
      REF_MIN = min(ALLELE_REF_DEPTH, na.rm = T),
      REF_MAX = max(ALLELE_REF_DEPTH, na.rm = T),
      ALT_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = T),
      ALT_MEDIAN = median(ALLELE_ALT_DEPTH, na.rm = T),
      ALT_MIN = min(ALLELE_ALT_DEPTH, na.rm = T),
      ALT_MAX = max(ALLELE_ALT_DEPTH, na.rm = T)
    ) %>%
    melt(
      id.vars = c("LOCUS", "POP_ID"),
      #     measure.vars = c(), # if left blank will use all the non id.vars
      variable.name = "COVERAGE_GROUP", 
      value.name = "VALUE"
    )
  if (missing(pop.levels) == "TRUE") {
    coverage <- coverage.sum.loci
  } else {
    coverage <- coverage.sum.loci %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
  }
  
  # by pop
  coverage.sum.pop <- data %>%
    group_by(POP_ID, INDIVIDUALS) %>%
    summarise(
      READ_DEPTH_MEAN = mean(READ_DEPTH, na.rm = T),
      READ_DEPTH_MEDIAN = median(READ_DEPTH, na.rm = T),
      READ_DEPTH_MIN = min(READ_DEPTH, na.rm = T),
      READ_DEPTH_MAX = max(READ_DEPTH, na.rm = T),
      ALLELE_REF_DEPTH_MEAN = mean(ALLELE_REF_DEPTH, na.rm = T),
      ALLELE_REF_DEPTH_MEDIAN = median(ALLELE_REF_DEPTH, na.rm = T),
      ALLELE_REF_DEPTH_MIN = min(ALLELE_REF_DEPTH, na.rm = T),
      ALLELE_REF_DEPTH_MAX = max(ALLELE_REF_DEPTH, na.rm = T),
      ALLELE_ALT_DEPTH_MEAN = mean(ALLELE_ALT_DEPTH, na.rm = T),
      ALLELE_ALT_DEPTH_MEDIAN = median(ALLELE_ALT_DEPTH, na.rm = T),
      ALLELE_ALT_DEPTH_MIN = min(ALLELE_ALT_DEPTH, na.rm = T),
      ALLELE_ALT_DEPTH_MAX = max(ALLELE_ALT_DEPTH, na.rm = T)
    ) %>%
    group_by(POP_ID) %>%
    summarise_each_(funs(mean), vars = c("READ_DEPTH_MEAN", "READ_DEPTH_MEDIAN", "READ_DEPTH_MIN", "READ_DEPTH_MAX", "ALLELE_REF_DEPTH_MEAN", "ALLELE_REF_DEPTH_MEDIAN", "ALLELE_REF_DEPTH_MIN", "ALLELE_REF_DEPTH_MAX", "ALLELE_ALT_DEPTH_MEAN", "ALLELE_ALT_DEPTH_MEDIAN", "ALLELE_ALT_DEPTH_MIN", "ALLELE_ALT_DEPTH_MAX")) %>%
    melt(
      id.vars = c("POP_ID"),
      variable.name = "GENOTYPE_LIKELIHOOD_GROUP", 
      value.name = "VALUE"
    )
  
  coverage.summary.total <- coverage.sum.pop %>%
    summarise_each(funs(mean))
  
  coverage.summary.total[1,1] <- "TOTAL"
  
  if (missing(pop.levels) == "TRUE") {
    coverage.summary.pop.total <- coverage.sum.pop %>%
      rbind(coverage.summary.total)
  } else {
    coverage.summary.pop.total <- coverage.sum.pop %>%
      rbind(coverage.summary.total) %>% 
      mutate(POP_ID = factor(POP_ID, levels = c(pop.levels, "TOTAL"), ordered = T)) %>%
      arrange(POP_ID)
  }
  coverage.summary.pop.total
  
  write.table(coverage.summary.pop.total, filename, sep = "\t", row.names = F, col.names = T, quote = F)
  
  invisible(cat(sprintf(
    "Filename: %s
    Written in this directory: %s", 
    filename,
    getwd()
  )))
  
  # results
  results <- list()
  results$coverage.summary.long <- coverage.sum.loci
  results$coverage.summary.pop <- coverage.summary.pop.total
  
  return(results)
  
}




#' @title Table of low coverage genotypes
#' @description This function create a table summary of the genotypes
#' below a user-define threshold.
#' coverage statistics by populations.
#' @param tidy.vcf.file The tidy VCF file created with read_stacks_vcf.
#' @param pop.levels Character string defining your ordered populations.
#' @param read.depth.threshold The read depth threshold to evaluate.
#' @param filename.low.coverage Filename of the low coverage table written
#' in the working directory.
#' @param filename.low.coverage.imbalance Filename of ...
#' @return a list of 2 tables (accessed with $). The values in the tables
#' represent percentage of samples.
#' @details work in progress....
#' Table 1: low coverage summary $low.coverage.summary (homo- and
#' hetero- zygotes genotypes).
#' Table 2: summary of coverage imbalance between alleles in the heterozygotes.
#' 0/0 : homozygote REF allele.
#' 1/1 : homozygote ALT allele.
#' 0/1 : heterozygote with coverage REF > ALT allele.
#' 1/0 : heterozygote with coverage REF < ALT allele.
#' @rdname table_low_coverage_summary
#' @export


table_low_coverage_summary <- function(tidy.vcf.file,
                                       pop.levels, 
                                       read.depth.threshold,
                                       filename.low.coverage,
                                       filename.low.coverage.imbalance) {
  
  INDIVIDUALS <- NULL
  POP_ID <- NULL
  GT <- NULL
  READ_DEPTH <- NULL
  ALLELE_COVERAGE_RATIO <- NULL
  SAMPLES_NUMBER <- NULL
  TOTAL_NUMBER <- NULL
  IMBALANCE_NUMBER <- NULL
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, 
                     col_names = T, 
                     col_types = "diidccddccccdddddc") %>%
      mutate(INDIVIDUALS = factor(INDIVIDUALS))
    message("Using the file in your directory")
    
  } else {
    data <- tidy.vcf.file
    message("Using the file from your global environment")
    
  }
  
  if (missing(pop.levels) == "TRUE") {
    data <- tidy.vcf.file
    
  } else {
    data <- tidy.vcf.file %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T))
  }
  
  low.coverage.summary <- data %>%
    filter(GT != "./.") %>%
    select(GT, POP_ID) %>%
    group_by(GT, POP_ID) %>%
    summarise(
      TOTAL_NUMBER = n()
    ) %>%
    full_join(
      data %>%
        filter(READ_DEPTH < read.depth.threshold & GT != "./.") %>%
        group_by(GT, POP_ID) %>%
        summarise(
          SAMPLES_NUMBER = n()
        ),
      by = c("GT", "POP_ID")
    ) %>%
    full_join(
      data %>%
        filter(READ_DEPTH < read.depth.threshold & GT != "./.") %>%
        filter(ALLELE_COVERAGE_RATIO != "NA" & ALLELE_COVERAGE_RATIO != 0 ) %>%
        group_by(GT, POP_ID) %>%
        summarise(
          IMBALANCE_NUMBER = n()
        ),
      by = c("GT", "POP_ID")
    ) %>%
    mutate(
      LOW_COVERAGE_PERCENT = round(SAMPLES_NUMBER / TOTAL_NUMBER * 100, 2),
      IMBALANCE_PERCENT = round(IMBALANCE_NUMBER / TOTAL_NUMBER * 100, 2)
    )
  
  low.coverage.summary.table <- low.coverage.summary %>%
    dcast(POP_ID ~ GT, value.var = "LOW_COVERAGE_PERCENT")
  
  write.table(low.coverage.summary.table, filename.low.coverage, sep = "\t", row.names = F, col.names = T, quote = F)
  
  low.coverage.imbalance.summary.table <- low.coverage.summary %>%
    filter(GT != "0/0" & GT != "1/1") %>%
    dcast(POP_ID ~ GT, value.var = "IMBALANCE_PERCENT")
  
  write.table(low.coverage.imbalance.summary.table, filename.low.coverage.imbalance, sep = "\t", row.names = F, col.names = T, quote = F)
  
  invisible(
    cat(
      sprintf(
        "2 files:
%s
%s\n
Written in the directory:
%s",
        filename.low.coverage, filename.low.coverage.imbalance, getwd()
      )))
  
  res <- list()
  res$low.coverage.summary <- low.coverage.summary.table
  res$heterozygote.imbalance <- low.coverage.imbalance.summary.table
  
  return(res)
}




## Genotype likelihood ###

#' @title Genotype likelihood summary
#' @description This function create 2 tables summary of the important
#' genotype likelihood statistics from the tidy vcf created
#' with \link{read_stacks_vcf}.
#' @param tidy.vcf.file The tidy VCF file created with read_stacks_vcf.
#' @param pop.levels Character string defining your ordered populations.
#' @param filename The name of the file written in the directory.
#' @return A list with 2 tables: the long format of loci ($gl.summary.long)
#' and populations genotype likelihood statistics
#' and the short format by populations ($gl.summary.pop).
#' The short-format is more user-friendly and
#' is written to the working directory.
#' @details The table contains summary statistics: mean, median, min, max and 
#' diff (max-min), of genotype likelihood by locus and populations. To access 
#' the two tables, use $. The table that summarise by populations was created
#' using average nested: loci -> individuals -> populations.
#' @rdname summary_genotype_likelihood
#' @export

summary_genotype_likelihood <- function(tidy.vcf.file, pop.levels, filename){
  
  POP_ID <- NULL
  GL <- NULL
  GL_MAX <- NULL
  GL_MIN <- NULL
  INDIVIDUALS <- NULL
  
  
  
  if (is.vector(tidy.vcf.file) == "TRUE") {
    data <- read_tsv(tidy.vcf.file, col_names = T, col_types = "iiiiccddcdccddddc")
    message("Using the file in your directory")
  } else {
    data <- tidy.vcf.file
    message("Using the file from your global environment")
    
  }
  
  GL.loci.pop <- data %>%
    group_by(LOCUS, POP_ID) %>%
    summarise(
      GL_MEAN = mean(GL, na.rm = T),
      GL_MEDIAN = median(GL, na.rm = T),
      GL_MIN = min(GL, na.rm = T),
      GL_MAX = max(GL, na.rm = T),
      GL_DIFF = GL_MAX - GL_MIN
    ) %>%
    melt(
      id.vars = c("LOCUS", "POP_ID"),
      variable.name = "GENOTYPE_LIKELIHOOD_GROUP", 
      value.name = "VALUE"
    )
  
  if (missing(pop.levels) == "TRUE") {
    GL.loci.pop.summary <- GL.loci.pop
  } else {
    GL.loci.pop.summary <- GL.loci.pop %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)) %>%
      arrange(POP_ID)
  }
  
  
  GL.pop <- data %>%
    group_by(POP_ID, INDIVIDUALS) %>%
    summarise(
      GL_MEAN = mean(GL, na.rm = T),
      GL_MEDIAN = median(GL, na.rm = T),
      GL_MIN = min(GL, na.rm = T),
      GL_MAX = max(GL, na.rm = T),
      GL_DIFF = GL_MAX - GL_MIN
    ) %>%
    group_by(POP_ID) %>%
    summarise_each_(funs(mean), vars = c("GL_MEAN", "GL_MEDIAN", "GL_MIN", "GL_MAX", "GL_DIFF")) %>%
    melt(
      id.vars = c("POP_ID"),
      variable.name = "GENOTYPE_LIKELIHOOD_GROUP", 
      value.name = "VALUE"
    )
  
  if (missing(pop.levels) == "TRUE") {
    GL.pop.summary <- GL.pop
  } else {
    GL.pop.summary <- GL.pop %>%
      mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = T)) %>%
      arrange(POP_ID)
  }
  
  GL.pop.summary.table <- GL.pop.summary %>%
    dcast(POP_ID ~ GENOTYPE_LIKELIHOOD_GROUP, value.var = "VALUE")
  
  
  message("Saving the summary table by pop in your working directory")
  write.table(GL.pop.summary.table, 
              filename,
              sep = "\t",
              row.names = F,
              col.names = T,
              quote = F
  )
  
  
  invisible(cat(sprintf(
    "Filename:
%s
Written in the directory:
%s",
    filename, getwd()
  )))
  
  # results
  results <- list()
  results$gl.summary.long <- GL.loci.pop.summary
  results$gl.summary.pop <- GL.pop.summary.table
  
  return(results)
  
}




#' @title Import and summarise the batch_x.phistats.tsv file
#' @description Import and summarise the batch_x.phistats.tsv file.
#' Necessary preparation for density distribution and box plot figures.
#' @param data The 'batch_x.phistats.tsv' created by STACKS.
#' @param skip.lines The number of line without the header 
#' to start reading the data.
#' @rdname summary_phistats
#' @export

summary_phistats <- function(data, skip.lines) {
  
  BATCH_ID <- NULL
  LOCUS <- NULL
  CHR <- NULL
  BP <- NULL
  POP_COUNT <- NULL
  PHI_ST <- NULL
  SMOOTHED_PHI_ST <- NULL
  SMOOTHED_PHI_ST_P_VALUE <- NULL
  PHI_CT <- NULL
  SMOOTHED_PHI_CT <- NULL
  SMOOTHED_PHI_CT_P_VALUE <- NULL
  PHI_SC <- NULL
  SMOOTHED_PHI_SC <- NULL
  SMOOTHED_PHI_SC_P_VALUE <- NULL
  FST_PRIME <- NULL
  SMOOTHED_FST_PRIME <- NULL
  SMOOTHED_FST_PRIME_P_VALUE <- NULL
  D_EST <- NULL
  SMOOTHED_D_EST <- NULL
  SMOOTHED_D_EST_P_VALUE <- NULL
  
  
  
  phistats <- read_tsv(data,
                       na = "NA",
                       skip = skip.lines,
                       progress = interactive(),
                       col_types = "iiciiddddddddddddddd",
                       col_names = c("BATCH_ID", "LOCUS", "CHR", "BP", "POP_COUNT", "PHI_ST", "SMOOTHED_PHI_ST", "SMOOTHED_PHI_ST_P_VALUE", "PHI_CT", "SMOOTHED_PHI_CT", "SMOOTHED_PHI_CT_P_VALUE", "PHI_SC", "SMOOTHED_PHI_SC", "SMOOTHED_PHI_SC_P_VALUE", "FST_PRIME", "SMOOTHED_FST_PRIME", "SMOOTHED_FST_PRIME_P_VALUE", "D_EST", "SMOOTHED_D_EST", "SMOOTHED_D_EST_P_VALUE")
  ) %>%
    select(-c(BATCH_ID, CHR, SMOOTHED_PHI_ST, SMOOTHED_PHI_ST_P_VALUE, SMOOTHED_PHI_CT, SMOOTHED_PHI_CT_P_VALUE, SMOOTHED_PHI_SC, SMOOTHED_PHI_SC_P_VALUE, SMOOTHED_FST_PRIME, SMOOTHED_FST_PRIME_P_VALUE, SMOOTHED_D_EST, SMOOTHED_D_EST_P_VALUE)) %>%
    melt(
      id.vars = c("LOCUS","BP","POP_COUNT"),
      variable.name = c("PHI_ST", "FST_PRIME", "D_EST"),
      value.name = "VALUE"
    )
}

