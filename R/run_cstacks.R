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
#' With default: \code{catalog.path = NULL} or with a supplied path, the function
#' The function scan automatically for the presence of a catalog inside the input folder.
#' If none is found, a new catalog is created. 
#' If your catalog is not in the input folder, supply a path here. 
#' e.g. \code{catalog.path = ~/catalog_folder}
#' 
#' 
#' , the catalog files are inside the 
#' input.path folder along the samples files and detected automatically. 
#' If a catalog is detected in the input folder, 
#' the samples in the \code{sample.list} argument
#' will be added in this catalog. The catalog is made of 3 files: 
#' \code{batch_1.catalog.alleles.tsv.gz, 
#' batch_1.catalog.snps.tsv.gz, 
#' batch_1.catalog.tags.tsv.gz}

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
#' @importFrom stringi stri_join stri_replace_all_fixed

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
  if (!dir.exists(input.path)) dir.create(input.path)
  if (!dir.exists("09_log_files")) dir.create("09_log_files")
  
  # Catalog editing ------------------------------------------------------------
  
  if (is.null(catalog.path)) { # no catalog path, searching in the input path...
    old.catalog <- list.files(path = input.path, pattern = "batch_")
    if (length(old.catalog) > 0 & length(old.catalog) == 3) {
      message("Found a catalog in the input folder, using files: ")
      message(stringi::stri_join(old.catalog, "\n"))
      
      catalog.path <- stringi::stri_replace_all_fixed(
        str = old.catalog[1], 
        pattern = ".catalog.alleles.tsv.gz", 
        replacement = "", 
        vectorize_all = FALSE
      )
      catalog.path <- stringi::stri_join(input.path, "/", catalog.path)
      catalog.path <- stringi::stri_join("--catalog ", shQuote(catalog.path))
    } 
    if (length(old.catalog) > 0 & length(old.catalog) < 3) {
      stop("Incomplete catalog, 3 files are required, see argument documentation")
    }
    
    if (length(old.catalog) == 0) {
      message("Builing catalog for the first time")
      catalog.path <- ""
    }
  } else {
    old.catalog <- list.files(path = catalog.path, pattern = "batch_")
    if (length(old.catalog) > 0 & length(old.catalog) == 3) {
      message("Found the catalog in the catalog path using files: ")
      message(stringi::stri_join(old.catalog, "\n"))
      
      catalog.path <- stringi::stri_replace_all_fixed(
        str = old.catalog[1], 
        pattern = ".catalog.alleles.tsv.gz", 
        replacement = "", 
        vectorize_all = FALSE
      )
      catalog.path <- stringi::stri_join(input.path, "/", catalog.path)
      catalog.path <- stringi::stri_join("--catalog ", shQuote(catalog.path))
    } 
    if (length(old.catalog) > 0 & length(old.catalog) < 3) {
      stop("Incomplete catalog, 3 files are required, see argument documentation")
    }
    
    if (length(old.catalog) == 0) {
      message("Builing catalog for the first time")
      catalog.path <- ""
    }
  }
    
    
    # cstacks options ------------------------------------------------------------
    b <- stringi::stri_join("-b ", b)
    
    o <- stringi::stri_join("-o ", shQuote(o))
    
    if (g) {
      g <- stringi::stri_join("-g ")
    } else {
      g <- ""
    }
    
    if (m) {
      m <- stringi::stri_join("-m ")
    } else {
      m <- ""
    }  
    
    n <- stringi::stri_join("-n ", n)
    p <- stringi::stri_join("-p ", p)
    
    if (h) {
      h <- stringi::stri_join("-h ")
    } else {
      h <- ""
    }
    
    
    # gapped assembly options ---------------------------------------------------
    if (gapped) {
      gapped <- stringi::stri_join("--gapped ")
    } else {
      gapped <- ""
    } 
    
    max_gaps <- stringi::stri_join("--max_gaps ", max_gaps)
    min_aln_len <- stringi::stri_join("--min_aln_len ", min_aln_len)
    
    # Advanced options ----------------------------------------------------------
    
    if (is.null(k_len)) {
      k_len <- ""
    } else {
      k_len <- stringi::stri_join("--k_len ", k_len)
    }
    
    if (report_mmatches) {
      report_mmatches <- stringi::stri_join("--report_mmatches ")
    } else {
      report_mmatches <- ""
    }
    
    # Samples to include in the catalog ------------------------------------------
    # s: filename prefix from which to load loci into the catalog.
    sample.list <- stringi::stri_join(input.path, "/", sample.list)
    s <- stringi::stri_join("-s ", shQuote(sample.list))
    
    # logs files -----------------------------------------------------------------
    file.date.time <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date.time <- stringi::stri_replace_all_fixed(
      file.date.time, 
      pattern = c("-", " ", ":"), 
      replacement = c("", "@", ""), 
      vectorize_all = FALSE
      )
    file.date.time <- stri_sub(file.date.time, from = 1, to = 13)
    
    log.file <- stringi::stri_join("09_log_files/cstacks_", file.date.time,".log")
    message(stringi::stri_join("For progress, look in the log file: ", log.file))
    
    
    # command args ---------------------------------------------------------------
    command.arguments <- paste(
      sample.list,
      input.path,
      catalog.path, b, o, g, m, n, p, h, gapped, max_gaps, min_aln_len, k_len,
      report_mmatches, s
    )
    
    # command
    system2(command = "cstacks", args = command.arguments, stdout = log.file, stderr = log.file)
    
    # # transfer back to s3
    # if (transfer.s3) {
    #   cstacks.files.to.s3 <- list.files(path = sample.list.path, pattern = individual, full.names = FALSE)
    #   purrr::walk(.x = cstacks.files.to.s3, .f = copy_s3, from.folder = from.folder, destination.folder = destination.folder)
    # }
  }# end run_cstacks
  
