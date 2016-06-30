#' @name run_sstacks
#' @title Run STACKS sstacks module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks} 
#' module inside R! 
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}
#' requires the 3 files from the catalog, created in 
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks} and 
#' for each sample, 4 files ending with 
#' \code{.alleles.tsv.gz, .models.tsv.gz, .snps.tsv.gz, .tags.tsv.gz}. 
#' Those files are created in the 
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' module.

#' @param input.path Path to input file. 
#' Default: \code{input.path = "06_ustacks_cstacks_sstacks"}

#' @param p Enable parallel execution with num_threads threads. 
#' Default: \code{ p = 4}

#' @param b MySQL ID of this batch. 
#' Default: \code{b = 1}.

#' @param catalog.prefix This is for the \code{c} option in
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}. 
#' \code{c: TSV file from which to load the catalog loci.}
#' Here, you give the prefix of the catalog file and the function take cares of
#' the rest. 
#' Default: \code{catalog.prefix = "batch_1"}. I'll

#' @param sample.list This is for the \code{s} option in
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}. 
#' \code{s: Filename prefix from which to load sample loci}.
#' Here, you have 2 choices: 1. you leave empty and let the function use the 
#' default:\code{sample.list = NULL} which will scan for the files in the 
#' \code{input.path} folder given above. 2. you supply a character string of the
#' samples. This could come from the \code{INDIVIDUALS_REP} column of the 
#' project info file, e.g. \code{sample.list = project.info$INDIVIDUALS_REP}.

#' @param o output path to write results.
#' Default: \code{o = "06_ustacks_cstacks_sstacks"}

#' @param g Base matching on genomic location, not sequence identity.
#' Default: \code{g = FALSE}

#' @param x Don't verify haplotype of matching locus.
#' Default: \code{x = FALSE}

#' @param v Print program version.
#' Default: \code{v = FALSE}

#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @param gapped Gapped assembly options: do you want to preform 
#' gapped alignments between stacks.
#' Default: \code{gapped = TRUE}


#' @rdname run_sstacks
#' @export
#' @import stringi
#' @import dplyr
#' @import readr

#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}
#' returns a \code{.matches.tsv.gz file for each sample}

#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' run_sstacks
#' that's it !
#' # Now if you have your own workflow folders, etc. Enter them like this:
#' run_sstacks (input.path = "/my/input/path", p = 32, b = 2, catalog.prefix = "batch_2", 
#' sample.list = c("ind1", "ind2", ...), o = "/my/output/path", g = FALSE,
#' x = FALSE, gapped = FALSE)
#' }

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("INDIVIDUALS_REP")
  )
}

#' @seealso 
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks} 

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011) 
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. 
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) 
#' Stacks: an analysis tool set for population genomics. 
#' Molecular Ecology, 22, 3124-3140.

# sstacks ----------------------------------------------------------------------
run_sstacks <- function(
  input.path = "06_ustacks_cstacks_sstacks",
  p = 4, 
  b = 1,
  catalog.prefix = "batch_1", 
  sample.list = NULL, 
  o = "06_ustacks_cstacks_sstacks", 
  g = FALSE, 
  x = FALSE,
  v = FALSE, 
  h = FALSE,
  gapped = TRUE
) {  
  
  cat("#######################################################################\n")
  cat("######################## stackr::run_sstacks ##########################\n")
  cat("#######################################################################\n")
  
  # Check directory ------------------------------------------------------------
  if(!dir.exists(input.path)) dir.create(input.path)
  if(!dir.exists("09_log_files")) dir.create("09_log_files")
  
  # sstacks options ------------------------------------------------------------
  
  # p: enable parallel execution with num_threads threads.
  p <- stri_paste("-p ", p)
  
  # b: MySQL ID of this batch.
  b <- stri_paste("-b ", b)
  
  # c: TSV file from which to load the catalog loci.
  # setwd("/Volumes/Drobo/esturgeon_drobo/sturgeon_manitoba/07-stacks_rx")
  # input.path <- "/Volumes/Drobo/esturgeon_drobo/sturgeon_manitoba/07-stacks_rx"
  # catalog.prefix <- "batch_1"
  
  list.catalog.files <- list.files(path = input.path, pattern = catalog.prefix)
  if (length(list.catalog.files) > 0) {
    message("Reading the catalog files...")
  } else {
    stop("Catalog files are required, check the input.path or catalog prefix
        arguments")
  }
  
  c <- stri_paste("-c ", input.path, "/", catalog.prefix)
  
  # o: output path to write results.
  o <- stri_paste("-o ", shQuote(o))
  
  # g: base matching on genomic location, not sequence identity.
  if(g) {
    g <- stri_paste("-g ")
  } else {
    g <- ""
  }
  
  # x: don't verify haplotype of matching locus.
  if(x) {
    x <- stri_paste("-x ")
  } else {
    x <- ""
  }  
  
  # v: print program version.
  if(v) {
    v <- stri_paste("-v ")
  } else {
    v <- ""
  }
  
  # h: display this help messsage.
  if (h) {
    h <- stri_paste("-h ")
  } else {
    h <- ""
  }
  
  # Gapped assembly options ---------------------------------------------------
  # --gapped: preform gapped alignments between stacks.
  if(gapped) {
    gapped <- stri_paste("--gapped ")
  } else {
    gapped <- ""
  } 
  
  # s: filename prefix from which to load sample loci---------------------------
  samples.in.folder <- data_frame(INDIVIDUALS_REP = list.files(input.path)) %>% 
    filter(!grepl("batch", INDIVIDUALS_REP))
  
  samples.matched <- samples.in.folder %>% 
    filter(grepl("matches", INDIVIDUALS_REP)) %>% 
    mutate(INDIVIDUALS_REP = stri_replace_all_fixed(
      str = INDIVIDUALS_REP, 
      pattern = ".matches.tsv.gz", replacement = "", 
      vectorized_all = FALSE)
    ) %>% 
    distinct(INDIVIDUALS_REP)
  message(stri_paste("Samples in folder already matched to the catalog: ", 
                     length(samples.matched$INDIVIDUALS_REP)))
  
  samples.to.match <- samples.in.folder %>% 
    filter(grepl("alleles", INDIVIDUALS_REP)) %>% 
    mutate(INDIVIDUALS_REP = stri_replace_all_fixed(
      str = INDIVIDUALS_REP, 
      pattern = ".alleles.tsv.gz", replacement = "", 
      vectorized_all = FALSE)
    ) %>% 
    distinct(INDIVIDUALS_REP) %>% 
    filter(!INDIVIDUALS_REP %in% samples.matched$INDIVIDUALS_REP)
  
  message(stri_paste("Samples in folder not matched to the catalog: ", length(samples.to.match$INDIVIDUALS_REP)))
  
  
  if (is.null(sample.list)) {
    s <- stri_paste("-s ", shQuote(stri_paste(input.path, "/", samples.to.match$INDIVIDUALS_REP)))
    message(stri_paste("Matching ", length(samples.to.match$INDIVIDUALS_REP), " sample(s) to the catalog..."))
  } else {
    sample.list <- purrr::keep(.x = sample.list, .p = sample.list %in% samples.to.match$INDIVIDUALS_REP)
    sample.list <- stri_paste(input.path, "/", sample.list)
    s <- stri_paste("-s ", shQuote(sample.list))
    message(stri_paste("Matching ", length(sample.list), " sample(s) to the catalog..."))
  }
  
  # logs files -----------------------------------------------------------------
  number.log.files <- length(list.files(path = "09_log_files", pattern = "sstacks"))
  if (number.log.files != 0) {
    log.number <- number.log.files + 1
    log.file <- stri_paste("09_log_files/sstacks_", log.number,".log")
  } else {
    log.file <- "09_log_files/sstacks.log"
  }
  
  # command args ---------------------------------------------------------------
  command.arguments <- c(
    p, b, c, s, o, g, x, v, h, gapped
  )
  
  # command
  system.time(
    system2(
      command = "sstacks", 
      args = command.arguments, 
      stdout = log.file, 
      stderr = log.file
    )
  )
  
  # # transfer back to s3
  # if (transfer.s3) {
  #   cstacks.files.to.s3 <- list.files(path = sample.list.path, pattern = individual, full.names = FALSE)
  #   purrr::walk(.x = cstacks.files.to.s3, .f = copy_s3, from.folder = from.folder, destination.folder = destination.folder)
  # }
  cat("#######################################################################\n")
}# end run_sstacks
