#' @name run_cstacks
#' @title Run STACKS cstacks module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/cstacks.php}{cstacks} 
#' module inside R!
#' Inside the folder \code{06_ustacks_cstacks_sstacks}, you should have:
#' \itemize{
#'   \item \strong{4 files for each samples:} The sample name is the prefix of 
#'   the files ending with:
#' \code{.alleles.tsv.gz, .models.tsv.gz, .snps.tsv.gz, .tags.tsv.gz}.
#' Those files are created in the 
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/ustacks.php}{ustacks}
#' module.
#' }

#' @param input.path Path to input file. 
#' Default: \code{input.path = "06_ustacks_cstacks_sstacks"}

#' @param sample.list This is for the \code{s} option in
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{cstacks}. 
#' \code{s: Filename prefix from which to load loci in the catalog}.
#' Here, you have 2 choices: 1. you leave empty and let the function use the 
#' default:\code{sample.list = NULL} which will scan for the files in the 
#' \code{input.path} folder given above. 2. you supply a character string of the
#' samples. This could come from the \code{INDIVIDUALS_REP} column of the 
#' project info file, e.g. \code{sample.list = project.info$INDIVIDUALS_REP}.

#' @param b MySQL ID of this batch. 
#' Default: \code{b = 1}.

#' @param o output path to write results.
#' Default: \code{o = "06_ustacks_cstacks_sstacks"}

#' @param g Base catalog matching on genomic location, not sequence identity.
#' Default: \code{g = FALSE}
#' 
#' @param m Include tags in the catalog that match to more than one entry.
#' Default: \code{m = FALSE}

#' @param n Number of mismatches allowed between sample tags when generating 
#' the catalog.
#' Default: \code{n = 1}

#' @param p Enable parallel execution with num_threads threads. 
#' Default: \code{p = 4}

#' @param h Display this help messsage.
#' Default: \code{h = FALSE}

#' @param catalog.path This is for the "Catalog editing" part in cstacks where 
#' you can provide the path to an existing catalog.
#' cstacks will add data to this existing catalog.
#' With default: \code{catalog.path = NULL}, the catalog files is inside the 
#' input.path filder along the samples files and detected automatically. 
#' In the catalog files are detected, the samples in the \code{sample.list} argument
#' will be included in this catalog.

#' @param gapped Gapped assembly options: do you want to preform 
#' gapped alignments between stacks.
#' Default: \code{gapped = TRUE}

#' @param max_gaps The number of gaps allowed between stacks before merging.
#' Default: \code{max_gaps = 2}

#' @param min_aln_len The minimum length of aligned sequence in a gapped 
#' alignment.
#' Default: \code{min_aln_len = 0.8}

#' @param k_len Specify k-mer size for matching between between catalog loci 
#' (automatically calculated by default).
#' Default: \code{k_len = NULL}

#' @param report_mmatches Report query loci that match more than one catalog locus.
#' Default: \code{report_mmatches = FALSE}

#' @rdname run_cstacks
#' @export
#' @import stringi
#' @import dplyr
#' @import readr

#' @return \href{http://catchenlab.life.illinois.edu/stacks/comp/sstacks.php}{sstacks}
#' returns a \code{.matches.tsv.gz file for each sample}

#' @details \strong{Computer or server problem during the cstacks ?} Look 
#' in the log file to see which individuals remains to be included. Create a 
#' new list of individuals to include and use the catalog.path argument to point
#' to the catalog created before the problem.

#' @examples
#' \dontrun{
#' # The simplest form of the function:
#' run_cstacks()
#' # that's it ! Now if you have your own workflow folders, etc. See below.
#' Next example, let say you only want to include 10 individuals/pop and 
#' include in the catalog samples with more than 2000000 reads. With the project
#' info file in the global environment:
#' individuals.catalog <- project.info.file) %>% 
#' filter(RETAINED > 2000000) %>% 
#' group_by(POP_ID) %>% 
#' sample_n(size = 10, replace = FALSE) %>% 
#' ungroup %>% 
#' arrange(desc(RETAINED)) %>% 
#' distinct(INDIVIDUALS_REP)
#' # The next line will give you the list of individuals to include
#' individuals.catalog <- individuals.catalog$INDIVIDUALS_REP
#' 
#' # To keep your info file updated with this information:
#' project.info.file <- project.info.file %>% 
#' mutate(CATALOG = if_else(INDIVIDUALS_REP %in% individuals.catalog, 
#' true = "catalog", false = "not_catalog")
#' )
#' write_tsv(project.info.file, "project.info.catalog.tsv")
#' 
#' # Then run the command this way:
#' run_cstacks (
#' sample.list = individuals.catalog, # or 
#' input.path = "06_ustacks_cstacks_sstacks",
#' catalog.path = NULL, 
#' b = 1, 
#' o = "06_ustacks_cstacks_sstacks", 
#' g = FALSE, 
#' m = FALSE, 
#' n = 1,
#' p = 32, 
#' h = FALSE,
#' gapped = TRUE, max_gaps = 2, min_aln_len = 0.8,
#' k_len = NULL, report_mmatches = FALSE
#' )
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

run_cstacks <- function(
  input.path = "06_ustacks_cstacks_sstacks",
  o = "06_ustacks_cstacks_sstacks", 
  b = 1, 
  catalog.path = NULL, 
  sample.list, 
  g = FALSE, 
  m = FALSE, 
  n = 1,
  p = 32, 
  h = FALSE,
  gapped = TRUE, max_gaps = 2, min_aln_len = 0.8,
  k_len = NULL, report_mmatches = FALSE
  # , transfer.s3 = FALSE, 
  # from.folder = NULL, destination.folder = NULL,
) {  
  
  # Check directory ------------------------------------------------------------
  if(!dir.exists(input.path)) dir.create(input.path)
  if(!dir.exists("09_log_files")) dir.create("09_log_files")
  
  # Catalog editing ------------------------------------------------------------
  old.catalog <- list.files(path = input.path, pattern = "batch_")
  if (length(old.catalog) > 0 & length(old.catalog) == 3) {
    message("Found a catalog in the input folder, using files: ")
    message(stri_paste(old.catalog, ", "))
  }
  if (is.null(catalog.path)) {
    catalog.path <- ""
  } else {
    catalog.path <- stri_paste("--catalog ", shQuote(catalog.path))
  }
  
  # cstacks options ------------------------------------------------------------
  b <- stri_paste("-b ", b)
  
  o <- stri_paste("-o ", shQuote(o))
  
  if(g) {
    g <- stri_paste("-g ")
  } else {
    g <- ""
  }
  
  if(m) {
    m <- stri_paste("-m ")
  } else {
    m <- ""
  }  
  
  n <- stri_paste("-n ", n)
  p <- stri_paste("-p ", p)
  
  if (h) {
    h <- stri_paste("-h ")
  } else {
    h <- ""
  }
  
  
  # gapped assembly options ---------------------------------------------------
  if(gapped) {
    gapped <- stri_paste("--gapped ")
  } else {
    gapped <- ""
  } 
  
  max_gaps <- stri_paste("--max_gaps ", max_gaps)
  min_aln_len <- stri_paste("--min_aln_len ", min_aln_len)
  
  # Advanced options ----------------------------------------------------------
  
  if (is.null(k_len)) {
    k_len <- ""
  } else {
    k_len <- stri_paste("--k_len ", k_len)
  }
  
  if(report_mmatches) {
    report_mmatches <- stri_paste("--report_mmatches ")
  } else {
    report_mmatches <- ""
  }
  
  # Samples to include in the catalog ------------------------------------------
  # s: filename prefix from which to load loci into the catalog.
  sample.list <- stri_paste(input.path, "/", sample.list)
  s <- stri_paste("-s ", shQuote(sample.list))
  
  # command args ---------------------------------------------------------------
  command.arguments <- paste(
    sample.list,
    input.path,
    catalog.path, b, o, g, m, n, p, h, gapped, max_gaps, min_aln_len, k_len,
    report_mmatches, s
  )
  
  # command
  system.time(system2(command = "cstacks", args = command.arguments, stdout = "09_log_files/cstacks.log", stderr = "09_log_files/cstacks.log"))
  
  # # transfer back to s3
  # if (transfer.s3) {
  #   cstacks.files.to.s3 <- list.files(path = sample.list.path, pattern = individual, full.names = FALSE)
  #   purrr::walk(.x = cstacks.files.to.s3, .f = copy_s3, from.folder = from.folder, destination.folder = destination.folder)
  # }
}# end run_cstacks
