#' @name run_process_radtags
#' @title Run STACKS process_radtags module
#' @description Run \href{http://catchenlab.life.illinois.edu/stacks/}{STACKS}
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php}{process_radtags}
#' module inside R. This function was tested on single-end sequences only.
#' If you want to contribute with other types of RAD data, let me know.

#' @param project.info (character, path) Path to the project info file.
#' The file is tab separated and as 3 columns named: \code{LANES}, \code{BRACODES}
#' and \code{INDIVIDUALS}. See details for more info.

#' @param path.seq.lanes (character, path) Path to sequencing lanes if
#' processing single-end sequences.
#' Same as \code{f} in STACKS command line.
#' Default: \code{path.seq.lanes = "03_sequencing_lanes"}.

#' @param paired.files (logical) files contained within the directory are paired.
#' Default: \code{paired.files = FALSE}.

#' @param interleaved (logical) specify that the paired-end reads are interleaved in
#' single files.
#' Default: \code{interleaved = FALSE}.

#' @param input.type (character) Input file type, either \code{"fastq"}, \code{"gzfastq"},
#' \code{"bam"}, or \code{bustard}. Default: \code{input.type = "guess"}.
#' @param output.path (character, path) Path to output the processed files.
#' Default: \code{output.path = "04_process_radtags"}.
#' @param output.type (character) Output file type: \code{"fastq"}, \code{"gzfastq"},
#' \code{"fasta"}, or \code{gzfasta}.
#' Default: \code{output.type = "guess"} (match input type).
#' @param pe.1 (character, path) First input file in a set of
#' paired-end sequences.
#' In stacks it's the argument \code{1}.
#' Default: \code{pe.1 = "FORWARD"}. Corresponding to the \code{FORWARD} column in
#' the project info file.
#' @param pe.2 (character, path) Second input file in a set of
#' paired-end sequences.
#' In stacks it's the argument \code{2}.
#' Default: \code{pe.2 = "REVERSE"}. Corresponding to the \code{REVERSE} column in
#' the project info file.

#' @param clean.data (logical) Clean data, remove any read with an uncalled base.
#' Default: \code{clean.data = TRUE}.
#' @param discard.reads (logical) Discard reads with low quality scores.
#' Default: \code{discard.reads = TRUE}.
#' @param rescue.barcodes (logical) Rescue barcodes and RAD-Tags.
#' Default: \code{rescue.barcodes = TRUE}.
#' @param truncate (integer) Truncate final read length to this value.
#' Default: \code{truncate = 90}.
#' @param capture.discard.reads (logical) Capture discarded reads to a file.
#' Default: \code{capture.discard.reads = FALSE}.


#' @param barcode.inline.null (logical) Barcode is inline with sequence,
#' occurs only on single-end read.
#' Default: \code{barcode.inline.null = TRUE}.
#' @param barcode.index.null Barcode is provided in FASTQ header
#' (Illumina i5 or i7 read).
#' Default: \code{barcode.index.null = FALSE}.
#' @param barcode.null.index Barcode is provded in FASTQ header
#' (Illumina i7 read if both i5 and i7 read are provided).
#' Default: \code{barcode.null.index = FALSE}.
#' @param barcode.inline.inline Barcode is inline with sequence,
#' occurs on single and paired-end read.
#' Default: \code{barcode.inline.inline = FALSE}.
#' @param barcode.index.index Barcode is provded in FASTQ header
#' (Illumina i5 and i7 reads).
#' Default: \code{barcode.index.index = FALSE}.
#' @param barcode.inline.index Barcode is inline with sequence on single-end
#' read and occurs in FASTQ header (from either i5 or i7 read).
#' Default: \code{barcode.inline.index = FALSE}.
#' @param barcode.index.inline Barcode occurs in FASTQ header
#' (Illumina i5 or i7 read) and is inline with single-end sequence
#' (for single-end data) on paired-end read (for paired-end data).
#' Default: \code{barcode.index.inline = FALSE}.

#' @param enzyme (character) Provide the restriction enzyme used
#' (cut site occurs on single-end read).
#' If double-digest use: \code{enzyme = NULL}.
#' Currently supported enzymes include:
#' \code{"aciI"}, \code{"ageI"}, \code{"aluI"}, \code{"apaLI"}, \code{"apeKI"},
#' \code{"apoI"}, \code{"aseI"}, \code{"bamHI"}, \code{"bbvCI"}, \code{"bfaI"},
#' \code{"bfuCI"}, \code{"bgIII"},
#' \code{"bsaHI"}, \code{"bspDI"}, \code{"bstYI"}, \code{"claI"}, \code{"csp6I"},
#' \code{"ddeI"}, \code{"dpnII"}, \code{"eaeI"}, \code{"ecoRI"}, \code{"ecoRV"},
#' \code{"ecoT22I"}, \code{"haeIII"}, \code{"hinP1I"}, \code{"hindIII"}, \code{"hpaII"},
#' \code{"kpnI"}, \code{"mluCI"}, \code{"mseI"}, \code{"mspI"}, \code{"ncoI"},
#' \code{"ndeI"}, \code{"nheI"}, \code{"nlaIII"}, \code{"nspI"}, \code{"notI"},
#' \code{"nsiI"},
#' \code{"pstI"}, \code{"rsaI"}, \code{"sacI"}, \code{"sau3AI"}, \code{"sbfI"},
#' \code{"sexAI"}, \code{"sgrAI"}, \code{"speI"}, \code{"sphI"}, \code{"taqI"},
#' \code{"xbaI"} or \code{"xhoI"}.
#' @param renz.1 (character) When a double digest was used,
#' provide the first restriction enzyme used.
#' @param renz.2 (character) When a double digest was used,
#' provide the second restriction enzyme used (cut site occurs on the paired-end read).

#' @param bestrad (logical) library was generated using BestRAD,
#' check for restriction enzyme on either read and potentially tranpose reads.
#' Default: \code{bestrad = FALSE}.

#' @param adapter.1 (character) Provide adaptor sequence that may occur on the
#' single-end read for filtering.
#' @param adapter.2 (character) Provide adaptor sequence that may occur on the
#' paired-read for filtering.
#' @param adapter.mm (integer) Number of mismatches allowed in the adapter sequence.
#' @param retain.header (logical) Retain unmodified FASTQ headers in the output.
#' Default: \code{retain.header = TRUE}.
#' @param merge (logical) If no barcodes are specified,
#' merge all input files into a single output file.
#' Default: \code{merge = FALSE}.
#' @param filter.illumina (logical) Discard reads that have been marked by
#' Illumina's chastity/purity filter as failing. With Ion Torrent sequencing, you
#' have to use \code{filter.illumina = FALSE}.
#' Default: \code{filter.illumina = TRUE}.
#' @param disable.rad.check (logical) Disable checking if the RAD site is intact.
#' Default: \code{disable.rad.check = FALSE}.

#' @param encoding (character) Specify how quality scores are encoded,
#' \code{"phred33"} for Illumina 1.8+/Sanger or \code{"phred64"} for Illumina 1.3-1.5.
#' Default: \code{encoding = "phred33"}.

#' @param window.size (double) Set the size of the sliding window as a fraction of the
#' read length, between 0 and 1.
#' Default: \code{window.size = 0.15}.

#' @param score.limit (integer) Set the score limit.
#' If the average score within the sliding window drops below this value,
#' the read is discarded.
#' Default: \code{score.limit = 10}.

#' @param len.limit (integer) Specify a minimum sequence length
#' (useful if your data has already been trimmed).
#' Default: \code{len.limit = NULL}.
#' @param barcode.dist.1 (integer) The number of allowed mismatches when
#' rescuing single-end barcodes.
#' Default: \code{barcode.dist.1 = 1}.
#' @param barcode.dist.2 (integer) The number of allowed mismatches when
#' rescuing paired-end barcodes.
#' Default: \code{barcode.dist.2 = NULL} (will default to barcode.dist.1).

#' @param threads (optional) The number of threads to run (stacks, max 24).
#' Default: \code{threads = min(parallel::detectCores(), 2)}.

#' @param parallel.core (optional) The number of core for parallel
#' processing. Much faster to have higher number here than with \code{threads}.
#' Default: \code{floor(parallel::detectCores() / 2)}.


#' @rdname run_process_radtags
#' @export
#' @return For stacks specific output see \href{http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php}{process_radtags}.
#'
#'

#' @details
#' \strong{Step 1. Individual naming}
#'
#' Please take a moment to think about the way you want to name your individuals...
#' Here is my recipe:
#' \enumerate{
#' \item Include metadata: SPECIES, POPULATIONS or SAMPLING SITES, MATURITY, YEAR of sampling, numerical ID with 3 to 4 digits.
#' \item 3 letters in ALL CAPS: capital letters are easier to read and will reduce confusion for other
#' people in deciphering your handwritting or the font you've used. The 3 letters will keep it short.
#' \item only 1 type of separator: the dash (-): why? it's the easiest separator to avoid confusion.
#' Please avoid the underscore (_) as it will sometimes be impossible to tell
#' if it's a whitespace or an underscore in printd lists or even in some codes.
#'
#' \strong{example:}
#'
#' Following this convention: SPECIES-POP-MATURITY-YEAR-ID
#'
#' My ids: STU-QUE-ADU-2017-0001, STU-ONT-JUV-2016-0002
#'
#' \strong{Species: }Sturgeon
#'
#' \strong{Sampling site: }Québec and Ontario
#'
#' \strong{Year of sampling: } 2016 and 2017
#'
#' \strong{MATURITY: } adult and juvenile
#'
#' \strong{ID: } 2 samples 0001 and 0002
#' }
#'
#' \strong{Step 2. project info file:}
#'
#' Create a tab separated file (e.g. in MS Excel) with 3 columns named:
#' \code{LANES}, \code{BRACODES} and \code{INDIVIDUALS}.
#' For each individuals you've just created,
#' give the barcodes and lanes name.
#'
#' \strong{REPLICATES ?}
#'
#' You have replicates? Awesome. stackr makes it easy to keep track of replicates.
#' Use the same name for individual replicates. They will have different barcodes,
#' and can potentially be on different lanes. No problem. stackr will combine
#' fastq file at the end, keeping original replicates intact. However, stackr
#' will be appending integers (e.g. STU-QUE-ADU-2017-0001-1, STU-QUE-ADU-2017-0001-2)
#' at the end of the names you've chosen). Combined replicates
#' will have -R at the end (e.g STU-QUE-ADU-2017-0001-R for the combination of the 2 replicates.)



#' @examples
#' \dontrun{
#' # library(stackr)
#' # If you haven't already build the folders to store all the files:
#' # stackr::build_stackr_workflow_dir()
#' #
#' # run a double digest process_radtags within R:
#' process.radtags.tuna <- stackr::run_process_radtags(
#' project.info = "02_project_info/project.info.tuna.tsv",
#' path.seq.lanes = "03_sequencing_lanes",
#' renz.1 = "pstI",
#' renz.2 = "mspI",
#' adapter.1 = "CGAGATCGGAAGAGCGGG",
#' adapter.mm = 2
#' )
#' # remaining arguments are defaults, so carefully look at them in the doc.
#' }


#' @seealso
#' \href{http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php}{process_radtags}.

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.

# run_process_radtags ----------------------------------------------------------
run_process_radtags <- function(
    project.info,
    path.seq.lanes = "03_sequencing_lanes",
    paired.files = FALSE,
    interleaved = FALSE,
    input.type = "guess",
    output.path = "04_process_radtags",
    output.type = "guess",
    pe.1 = NULL,
    pe.2 = NULL,
    clean.data = TRUE,
    discard.reads = TRUE,
    rescue.barcodes = TRUE,
    truncate = 90,
    capture.discard.reads = FALSE,
    barcode.inline.null = TRUE,
    barcode.index.null = FALSE,
    barcode.null.index = FALSE,
    barcode.inline.inline = FALSE,
    barcode.index.index = FALSE,
    barcode.inline.index = FALSE,
    barcode.index.inline = FALSE,
    enzyme = NULL,
    renz.1 = NULL,
    renz.2 = NULL,
    bestrad = FALSE,
    adapter.1 = NULL,
    adapter.2 = NULL,
    adapter.mm = NULL,
    retain.header = TRUE,
    merge = FALSE,
    filter.illumina = TRUE,
    disable.rad.check = FALSE,
    encoding = "phred33",
    window.size = 0.15,
    score.limit = 10,
    len.limit = NULL,
    barcode.dist.1 = 1,
    barcode.dist.2 = NULL,
    threads = min(parallel::detectCores(), 2),
    parallel.core = floor(parallel::detectCores() / 2)
) {
  cat("#######################################################################\n")
  cat("#################### stackr::run_process_radtags ######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()


  # Check directory ------------------------------------------------------------
  if (!dir.exists(output.path)) dir.create(output.path)
  if (!dir.exists("09_log_files")) dir.create("09_log_files")
  if (!dir.exists("09_log_files/01_process_radtags_log")) dir.create("09_log_files/01_process_radtags_log")
  if (!dir.exists("02_project_info")) dir.create("02_project_info")
  if (!dir.exists("08_stacks_results")) dir.create("08_stacks_results")
  if (!dir.exists("08_stacks_results/01_process_radtags_results")) dir.create("08_stacks_results/01_process_radtags_results")
  if (!dir.exists("02_project_info/01_lanes_barcodes_completed")) dir.create("02_project_info/01_lanes_barcodes_completed")


  # Date and time --------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Import project info file ---------------------------------------------------
  message("Importing project info")
  if (paired.files) {
    project.info.file <- readr::read_tsv(file = project.info, col_types = "cccc")
  } else {
    project.info.file <- readr::read_tsv(file = project.info, col_types = "ccc")
  }

  # Manage duplicate ID
  message("Scanning for duplicate IDs...")
  duplicate.id <- project.info.file %>%
    dplyr::group_by(INDIVIDUALS) %>%
    dplyr::tally(.) %>%
    dplyr::filter(n > 1 & !is.na(INDIVIDUALS))

  if (nrow(duplicate.id) > 0) {
    duplicate.id <- duplicate.id %>%
      dplyr::distinct(INDIVIDUALS) %>%
      dplyr::inner_join(project.info.file, by = "INDIVIDUALS") %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::mutate(REPLICATES = seq(1, n()))
  } else {
    duplicate.id <- project.info.file %>%
      dplyr::mutate(REPLICATES = rep(NA, n()))
  }

  # Manage horrible sequencing lanes names

  # paired-end
  if (paired.files) {
    # here we use LANE but it's really the group of paired-end files
    short.name.lanes <- project.info.file %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(FORWARD, REVERSE) %>%
      dplyr::mutate(
        LANES_SHORT = stringi::stri_join(rep("LANE"), stringi::stri_pad_left(
          str = seq(1, nrow(.)), width = 2, pad = "0"), sep = "_")
      )
    project.info.file <- project.info.file %>%
      dplyr::left_join(duplicate.id, by = c("INDIVIDUALS", "BARCODES", "FORWARD", "REVERSE")) %>%
      dplyr::left_join(short.name.lanes, by = c("FORWARD", "REVERSE")) %>%
      dplyr::arrange(LANES_SHORT, INDIVIDUALS) %>%
      dplyr::mutate(
        INDIVIDUALS_REP = ifelse(!is.na(REPLICATES), stringi::stri_paste(INDIVIDUALS, REPLICATES, sep = "-"), INDIVIDUALS),
        REPLICATES = stringi::stri_replace_na(str = REPLICATES, replacement = 0)
      ) %>%
      dplyr::select(LANES_SHORT, FORWARD, REVERSE, BARCODES, INDIVIDUALS, INDIVIDUALS_REP, REPLICATES) %>%
      dplyr::filter(!is.na(INDIVIDUALS)) %>%
      dplyr::filter(!is.na(BARCODES)) %>%
      dplyr::arrange(LANES_SHORT, INDIVIDUALS)

  } else {
    short.name.lanes <- project.info.file %>%
      dplyr::ungroup(.) %>%
      dplyr::distinct(LANES) %>%
      dplyr::arrange(LANES)

    pad.length <- stringi::stri_length(nrow(short.name.lanes))

    short.name.lanes %<>%
      dplyr::mutate(
        LANES_SHORT = stringi::stri_join(rep("LANE"), stringi::stri_pad_left(
          str = seq(1, nrow(.)), width = pad.length, pad = "0"), sep = "_")
      )

    project.info.file <- project.info.file %>%
      dplyr::left_join(duplicate.id, by = c("INDIVIDUALS", "BARCODES", "LANES")) %>%
      dplyr::left_join(short.name.lanes, by = c("LANES")) %>%
      dplyr::arrange(LANES_SHORT, INDIVIDUALS) %>%
      dplyr::mutate(
        INDIVIDUALS_REP = ifelse(!is.na(REPLICATES), stringi::stri_paste(INDIVIDUALS, REPLICATES, sep = "-"), INDIVIDUALS),
        REPLICATES = stringi::stri_replace_na(str = REPLICATES, replacement = 0)
      ) %>%
      dplyr::select(LANES, LANES_SHORT, BARCODES, INDIVIDUALS, INDIVIDUALS_REP, REPLICATES) %>%
      dplyr::filter(!is.na(INDIVIDUALS)) %>%
      dplyr::filter(!is.na(BARCODES)) %>%
      dplyr::arrange(LANES_SHORT, INDIVIDUALS)
  }


  readr::write_tsv(
    x = project.info.file,
    file = stringi::stri_join("02_project_info/project.info.", file.date, ".tsv")
    )

  sample.per.lanes <- project.info.file %>%
    dplyr::group_by(LANES_SHORT) %>%
    dplyr::tally(.) %>%
    dplyr::arrange(LANES_SHORT)
  # write_tsv(sample.per.lanes, "sample.per.lanes.tsv")

  # get the list of sequencing lane present in the folder
  # lane.list <- list.files(path = path.seq.lanes, full.names = TRUE)
  # check lanes in project info file and directory
  lanes.folder <- list.files(path = path.seq.lanes, full.names = FALSE)
  message("Number of LANES in the folder: ", length(lanes.folder))
  if (paired.files) {
    lanes.todo <- unique(project.info.file$LANES_SHORT)
    if (FALSE %in% (unique(sort(unique(lanes.folder)) %in% sort(unique(c(project.info.file$FORWARD, project.info.file$REVERSE)))))) {
      message("ATTENTION: Lane names don't match between project info file and lanes in folder...")
      message("process_radtags will run the lanes specified in the project info file")
    }
    lane.list <- purrr::keep(.x = lanes.folder, .p = lanes.folder %in% lanes.todo)
  } else {
    lanes.todo <- unique(project.info.file$LANES)
    message("Number of LANES TODO from the project.info.file: ", length(lanes.todo))
    # no.problem <- unique(lanes.todo %in% lanes.folder)
    no.problem <- intersect(lanes.todo, lanes.folder)
    message("Number of LANES in common between the project.info.file and the actual files in folder: ", length(no.problem))

    if (length(lanes.todo) != length(no.problem)) {
      message("\n\nATTENTION: Lane names don't match between project info file and lanes in folder...")
      message("process_radtags will run the lanes specified in the project info file found in the folder")
      problematic.lanes <- tibble::tibble("LANES_PROBLEM" = setdiff(lanes.todo, lanes.folder))
      readr::write_tsv(x = problematic.lanes,file =  stringi::stri_join("02_project_info/problematic.lanes.", file.date, ".tsv"))
    }

    if (length(lanes.todo) != length(lanes.folder)) {
      message("\n\nATTENTION: Lane names don't match between project info file and lanes in folder...")
      message("process_radtags will run the lanes specified in the project info file found in the folder")
      problematic.lanes <- tibble::tibble("LANES_PROBLEM" = setdiff(lanes.folder, lanes.todo))
      readr::write_tsv(x = problematic.lanes, file = stringi::stri_join("02_project_info/problematic.lanes.", file.date, ".tsv"))
    }


    if (!length(no.problem) > 1) {
      rlang::abort("\n\nProblem: Lanes in project info file don't match the ones in the folder....")
    }
    lane.list <- stringi::stri_join(path.seq.lanes, "/", no.problem)
    project.info.file %<>% dplyr::filter(LANES %in% no.problem)
    readr::write_tsv(
      x = project.info.file,
      file = stringi::stri_join("02_project_info/project.info.todo.", file.date, ".tsv")
      )
  }

  process_radtags_lane <- function(
    lane.list,
    file.date = file.date,
    project.info.file = project.info.file,
    path.seq.lanes = "03_sequencing_lanes",
    paired.files = FALSE,
    interleaved = FALSE,
    input.type = "guess",
    output.path = "04_process_radtags",
    output.type = "guess",
    pe.1 = NULL,
    pe.2 = NULL,
    clean.data = TRUE,
    discard.reads = TRUE,
    rescue.barcodes = TRUE,
    truncate = 90,
    capture.discard.reads = FALSE,
    barcode.inline.null = TRUE,
    barcode.index.null = FALSE,
    barcode.null.index = FALSE,
    barcode.inline.inline = FALSE,
    barcode.index.index = FALSE,
    barcode.inline.index = FALSE,
    barcode.index.inline = FALSE,
    enzyme = NULL,
    renz.1 = NULL,
    renz.2 = NULL,
    bestrad = FALSE,
    adapter.1 = NULL,
    adapter.2 = NULL,
    adapter.mm = NULL,
    retain.header = TRUE,
    merge = FALSE,
    filter.illumina = TRUE,
    disable.rad.check = FALSE,
    encoding = "phred33",
    window.size = 0.15,
    score.limit = 10,
    len.limit = NULL,
    barcode.dist.1 = 1,
    barcode.dist.2 = NULL,
    threads = min(parallel::detectCores(), 2)
  ) {
    # lane.list <- lane.list[1]

    # generate a barcode file for the lane -------------------------------------
    if (paired.files) {
      info <- dplyr::ungroup(project.info.file) %>%
        dplyr::filter(LANES_SHORT %in% lane.list)
    } else {
      lane.todo <- basename(lane.list)
      info <- dplyr::ungroup(project.info.file) %>%
        dplyr::filter(LANES %in% lane.todo)
    }
    lane.short <- unique(info$LANES_SHORT)
    barcode.file <- info %>% dplyr::select(BARCODES, INDIVIDUALS_REP)
    barcode.filename <- stringi::stri_join("02_project_info/barcodes_id", "_", lane.short, ".txt")
    lane.completed <- stringi::stri_join("02_project_info/01_lanes_barcodes_completed/barcodes_id", "_", lane.short, ".txt")
    readr::write_tsv(x = barcode.file, file = barcode.filename, col_names = FALSE)

    # process_radtags_options --------------------------------------------------
    if (!paired.files) {
      f <- stringi::stri_join("-f ", shQuote(lane.list))
      P <- ""
      pe.1 <- ""
      pe.2 <- ""
      paired.analysis <- FALSE
    } else {
      f <- ""
      P <- "-P "
      pe.1 <- stringi::stri_join("-1 ", shQuote(file.path(path.seq.lanes, unique(info$FORWARD))))
      pe.2 <- stringi::stri_join("-2 ", shQuote(file.path(path.seq.lanes, unique(info$REVERSE))))
      paired.analysis <- TRUE
    }

    i <- ""
    y <- ""
    I <- ""

    if (input.type != "guess") i <- stringi::stri_join("-i ", shQuote(input.type))
    if (output.type != "guess") y <- stringi::stri_join("-y ", shQuote(output.type))
    if (interleaved) I <- "-I "    # for paired-end data



    temp.dir <- stringi::stri_join(output.path, lane.short, sep = "/")
    dir.create(temp.dir)
    o <- stringi::stri_join("-o ", shQuote(temp.dir))
    b <- stringi::stri_join("-b ", shQuote(barcode.filename))
    c <- ""
    q <- ""
    r <- ""
    D <- ""

    if (clean.data) c <- stringi::stri_join("-c ")
    if (discard.reads) q <- stringi::stri_join("-q ")
    if (rescue.barcodes) r <- stringi::stri_join("-r ")
    if (capture.discard.reads) D <- stringi::stri_join("-D ")
    t <- stringi::stri_join("-t ", truncate)


    # BARCODES OPTIONS -----------------------------------------------------------
    if (barcode.inline.null) {
      barcode.inline.null <- "--inline-null"
    } else {
      barcode.inline.null <- ""
    }
    if (barcode.index.null) {
      barcode.index.null <- "--index-null"
    } else {
      barcode.index.null <- ""
    }
    if (barcode.null.index) {
      barcode.null.index <- "--null-index"
    } else {
      barcode.null.index <- ""
    }

    if (barcode.inline.inline) {
      barcode.inline.inline <- "--inline-inline"
    } else {
      barcode.inline.inline <- ""
    }
    if (barcode.index.index) {
      barcode.index.index <- "--index-index"
    } else {
      barcode.index.index <- ""
    }
    if (barcode.inline.index) {
      barcode.inline.index <- "--inline-index"
    } else {
      barcode.inline.index <- ""
    }
    if (barcode.index.inline) {
      barcode.index.inline <- "--index-inline"
    } else {
      barcode.index.inline <- ""
    }

    # Restriction enzyme options--------------------------------------------------
    if (is.null(enzyme)) {
      enzyme <- ""
    } else {
      enzyme <- stringi::stri_join("-e ", shQuote(enzyme))
    }
    if (is.null(renz.1)) {
      renz.1 <- ""
    } else {
      renz.1 <- stringi::stri_join("--renz-1 ", shQuote(renz.1))
    }
    if (is.null(renz.2)) {
      renz.2 <- ""
    } else {
      renz.2 <- stringi::stri_join("--renz-2 ", shQuote(renz.2))
    }

    # Protocol-specific options:
    # --bestrad: library was generated using BestRAD, check for restriction enzyme on either read and potentially tranpose reads.
    if (bestrad) {
      bestrad <- "--bestrad "
    } else {
      bestrad <- ""
    }

    #  Adapter options: ------------------------------------------------------------
    if (is.null(adapter.1)) {
      adapter.1 <- ""
    } else {
      adapter.1 <- stringi::stri_join("--adapter-1 ", shQuote(adapter.1))
    }
    if (is.null(adapter.2)) {
      adapter.2 <- ""
    } else {
      adapter.2 <- stringi::stri_join("--adapter-2 ", shQuote(adapter.2))
    }
    if (is.null(adapter.mm)) {
      adapter.mm <- ""
    } else {
      adapter.mm <- stringi::stri_join("--adapter-mm ", adapter.mm)
    }
    # Output options--------------------------------------------------------------
    if (retain.header) {
      retain.header <- stringi::stri_join("--retain-header")
    } else {
      retain.header <- ""
    }
    if (merge) {
      merge <- stringi::stri_join("--merge")
    } else {
      merge <- ""
    }

    # Advanced options -----------------------------------------------------------
    if (filter.illumina) {
      filter.illumina <- stringi::stri_join("--filter-illumina")
    } else {
      filter.illumina <- ""
    }
    if (disable.rad.check) {
      disable.rad.check <- disable.rad.check("--disable-rad-check")
    } else {
      disable.rad.check <- ""
    }
    if (is.null(!len.limit)) {
      len.limit <- stringi::stri_join("--len-limit ", len.limit)
    } else {
      len.limit <- ""
    }

    if (rescue.barcodes) {
      if (paired.analysis) {
        if (is.null(barcode.dist.2)) barcode.dist.2 <- barcode.dist.1
        barcode.dist.2 <- stringi::stri_join("--barcode-dist-2 ", barcode.dist.2)
      } else {
        barcode.dist.2 <- ""
      }
      barcode.dist.1 <- stringi::stri_join("--barcode-dist-1 ", barcode.dist.1)
    }

    # modified and new arguments as of stacks v.2.62 ---------------------------
    # number of threads to run (max 24)
    # available.threads <- parallel::detectCores()
    threads <- stringi::stri_join("--threads ", threads)
    # threads <- " "

    # encoding
    encoding <- stringi::stri_join("--encoding ", shQuote(encoding))

    # set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).
    window.size <- stringi::stri_join("--window-size ", window.size)

    # set the phred score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).
    score.limit <- stringi::stri_join("--score-limit ", score.limit)

    # Command ------------------------------------------------------------------
    command.arguments <- paste(
      P,
      I,
      i,
      b,
      o,
      f,
      pe.1,
      pe.2,
      c,
      q,
      r,
      t,
      D,
      y,
      barcode.inline.null,
      barcode.index.null,
      barcode.null.index,
      barcode.inline.inline,
      barcode.index.index,
      barcode.inline.index,
      barcode.index.inline,
      enzyme,
      renz.1, renz.2,
      bestrad,
      adapter.1,
      adapter.2,
      adapter.mm,
      retain.header,
      merge,
      filter.illumina,
      disable.rad.check,
      encoding,
      window.size,
      score.limit,
      len.limit,
      barcode.dist.1,
      barcode.dist.2,
      threads
    )

    # log file -----------------------------------------------------------------
    process.radtags.log.file <- stringi::stri_join("09_log_files/01_process_radtags_logs/process_radtags_", lane.short, "_", file.date, ".log")

    # running the command ------------------------------------------------------
    system2(command = "process_radtags", args = command.arguments, stderr = process.radtags.log.file)

    # results --------------------------------------------------------------------
    # Moving completed lanes
    file.rename(from = barcode.filename, to = lane.completed)


    # Importing log file created by stacks to summarise and rename
    log.file <- list.files(path = temp.dir, pattern = "process_radtags", full.names = TRUE)

    # lane statistics
    lanes.stats <- readr::read_tsv(
      file = log.file,
      col_names = c("DESCRIPTION", "READS"),
      col_types = "ci",
      trim_ws = TRUE,
      skip = 12,
      n_max = 1
    ) %>%
      tibble::add_column("PCT" = "100%") %>%
      dplyr::bind_rows(
        readr::read_tsv(
          file = log.file,
          col_names = c("DESCRIPTION", "READS", "PCT"),
          col_types = "ci",
          trim_ws = TRUE,
          skip = 13,
          n_max = 6
        )
      ) %>%
      dplyr::mutate(
        DESCRIPTION = stringi::stri_trans_toupper(str = DESCRIPTION),
        DESCRIPTION = stringi::stri_replace_all_charclass(DESCRIPTION, "\\p{WHITE_SPACE}", "_"),
        LANES_SHORT = rep(lane.short, dplyr::n())
      ) %>%
      dplyr::group_by(LANES_SHORT) %>%
      tidyr::pivot_wider(data = ., names_from = DESCRIPTION, values_from = c(READS, PCT))

    # save lane stats
    readr::write_tsv(
      x = lanes.stats,
      file = stringi::stri_join(
        "08_stacks_results/01_process_radtags_results/",
        stringi::stri_join("process_radtags_", lane.short, "_stats"), ".tsv"))

    # barcode stats
    temp.file <- suppressWarnings(suppressMessages(readr::read_table(file = log.file, col_names = FALSE)))

    skip.number <- which(stringi::stri_detect_fixed(
      str = temp.file$X1,
      pattern = "BEGIN"))[3] + 4

    barcodes.stats <- suppressMessages(
      suppressWarnings(
        readr::read_delim(
          file = log.file,
          delim = "\t",
          skip = skip.number,
          n_max = nrow(barcode.file),
          col_names = TRUE,
          progress = interactive()
        ) %>%
          dplyr::mutate(LANES_SHORT = rep(lane.short, dplyr::n()))
      )
    )
    # rename and move log file in log folder
    new.log.file.name <- stringi::stri_join(
      "09_log_files/01_process_radtags_logs/process_radtags_summary_log_", lane.short, ".log")
    file.rename(from = log.file, to = new.log.file.name)

    # save barcodes stats
    barcode.res.list.name <- stringi::stri_join(
      "process_radtags_", lane.short, "_barcodes_stats")
    readr::write_tsv(
      x = barcodes.stats,
      file = stringi::stri_join("08_stacks_results/01_process_radtags_results/", barcode.res.list.name, ".tsv"))

    return(barcodes.stats)
  } # end process_radtags_lane


  # run in parallel ------------------------------------------------------------
  message("\nRunning several process_radtags in parallel")
  message("Monitoring progress in real time available in log files")
  names(lane.list) <- lane.list

  # process.radtags.results <- .stackr_parallel_mc(
  #   X = lane.list,
  #   FUN = process_radtags_lane,
  #   mc.preschedule = FALSE,
  #   mc.silent = FALSE,
  #   mc.cleanup = TRUE,
  #   mc.cores = parallel.core,
  #   file.date = file.date,
  #   project.info.file = project.info.file,
  #   path.seq.lanes = path.seq.lanes,
  #   paired.files = paired.files,
  #   interleaved = interleaved,
  #   input.type = input.type,
  #   output.path = output.path,
  #   output.type = output.type,
  #   pe.1 = pe.1,
  #   pe.2 = pe.2,
  #   clean.data = clean.data,
  #   discard.reads = discard.reads,
  #   rescue.barcodes = rescue.barcodes,
  #   truncate = truncate,
  #   capture.discard.reads = capture.discard.reads,
  #   barcode.inline.null = barcode.inline.null,
  #   barcode.index.null = barcode.index.null,
  #   barcode.null.index = barcode.null.index,
  #   barcode.inline.inline = barcode.inline.inline,
  #   barcode.index.index = barcode.index.index,
  #   barcode.inline.index = barcode.inline.index,
  #   barcode.index.inline = barcode.index.inline,
  #   enzyme = enzyme,
  #   renz.1 = renz.1,
  #   renz.2 = renz.2,
  #   bestrad = bestrad,
  #   adapter.1 = adapter.1,
  #   adapter.2 = adapter.2,
  #   adapter.mm = adapter.mm,
  #   retain.header = retain.header,
  #   merge = merge,
  #   filter.illumina = filter.illumina,
  #   disable.rad.check = disable.rad.check,
  #   encoding = encoding,
  #   window.size = window.size,
  #   score.limit = score.limit,
  #   len.limit = len.limit,
  #   barcode.dist.1 = barcode.dist.1,
  #   barcode.dist.2 = barcode.dist.2,
  #   threads = threads
  # )
  process.radtags.results <- stackr_future(
    .x = lane.list,
    .f = process_radtags_lane,
    flat.future = "dfr",
    parallel.core = parallel.core,
    file.date = file.date,
    project.info.file = project.info.file,
    path.seq.lanes = path.seq.lanes,
    paired.files = paired.files,
    interleaved = interleaved,
    input.type = input.type,
    output.path = output.path,
    output.type = output.type,
    pe.1 = pe.1,
    pe.2 = pe.2,
    clean.data = clean.data,
    discard.reads = discard.reads,
    rescue.barcodes = rescue.barcodes,
    truncate = truncate,
    capture.discard.reads = capture.discard.reads,
    barcode.inline.null = barcode.inline.null,
    barcode.index.null = barcode.index.null,
    barcode.null.index = barcode.null.index,
    barcode.inline.inline = barcode.inline.inline,
    barcode.index.index = barcode.index.index,
    barcode.inline.index = barcode.inline.index,
    barcode.index.inline = barcode.index.inline,
    enzyme = enzyme,
    renz.1 = renz.1,
    renz.2 = renz.2,
    bestrad = bestrad,
    adapter.1 = adapter.1,
    adapter.2 = adapter.2,
    adapter.mm = adapter.mm,
    retain.header = retain.header,
    merge = merge,
    filter.illumina = filter.illumina,
    disable.rad.check = disable.rad.check,
    encoding = encoding,
    window.size = window.size,
    score.limit = score.limit,
    len.limit = len.limit,
    barcode.dist.1 = barcode.dist.1,
    barcode.dist.2 = barcode.dist.2,
    threads = threads
  )
  process.radtags.results %<>%
    dplyr::bind_rows(.) %>%
    dplyr::rename_with(.data = ., .fn = toupper) %>%
    dplyr::rename_with(
      .data = .,
      .fn = stringi::stri_replace_all_charclass,
      pattern = "\\p{WHITE_SPACE}",
      replacement = "_",
      vectorize_all = FALSE
    ) %>%
    dplyr::mutate(FILENAME = as.character(FILENAME)) %>%
    dplyr::select(
      INDIVIDUALS_REP = FILENAME,
      BARCODES = BARCODE,
      TOTAL,
      NO_RADTAG = RAD_CUTSITE_NOT_FOUND, # replaced NORADTAG v2.60
      LOW_QUALITY,
      RETAINED_READS,
      tidyselect::everything()
    )


  # transfer the fq.gz file from the separate folder into a new folder

  # get the names of the folder
  folder.list <- list.files(path = output.path, full.names = TRUE)
  # folder.list

  destination.folder <- output.path

  transfer_folder_individuals_fq <- function(folder.list) {
    fq.names <- list.files(path = folder.list)

    transfer_individual_fq <- function(fq.names) {
      fq.names.full.path <- stringi::stri_join(folder.list, fq.names, sep = "/")

      # move to destination folder
      new.fq.path <- stringi::stri_join(destination.folder, fq.names, sep = "/")
      file.rename(from = fq.names.full.path, to = new.fq.path)
    }
    purrr::walk(.x = fq.names, .f = transfer_individual_fq)
  }
  purrr::walk(.x = folder.list, .f = transfer_folder_individuals_fq)

  #Remove foldes of lanes
  remove_lane_folder <- function(folder.list) file.remove(folder.list)
  purrr::walk(.x = folder.list, .f = remove_lane_folder)

  # combine replicates in a new fq file-----------------------------------------
  # message("Scanning for replicates...")
  if (paired.files) {
    project.info.file <- project.info.file %>%
      dplyr::mutate(
        FQ_FILES_F = stringi::stri_join(INDIVIDUALS_REP, ".1.fq.gz"),
        FQ_FILES_R = stringi::stri_join(INDIVIDUALS_REP, ".2.fq.gz")
      ) %>%
      dplyr::inner_join(process.radtags.results, by = "INDIVIDUALS_REP")
  } else {
    project.info.file <- project.info.file %>%
      dplyr::mutate(FQ_FILES = stringi::stri_join(INDIVIDUALS_REP, ".fq.gz")) %>%
      dplyr::inner_join(process.radtags.results, by = c("INDIVIDUALS_REP", "BARCODES", "LANES_SHORT"))
  }
  readr::write_tsv(
    x = project.info.file,
    file = stringi::stri_join("02_project_info/project.info.", file.date, ".tsv")
    )

  # replicates
  replicates <- project.info.file %>%
    dplyr::filter(REPLICATES >= 1) %>%
    dplyr::arrange(INDIVIDUALS_REP) %>%
    dplyr::distinct(INDIVIDUALS_REP, .keep_all = TRUE)

  replicate.presence <- nrow(replicates)

  if (replicate.presence > 0) {
    if (paired.files) {
      replicates.sum <- replicates %>%
        dplyr::group_by(INDIVIDUALS) %>%
        dplyr::summarise(
          TOTAL = sum(TOTAL),
          NO_RADTAG = sum(NO_RADTAG),
          LOW_QUALITY = sum(LOW_QUALITY),
          RETAINED_READS = sum(RETAINED_READS),
          .groups = "keep"
        ) %>%
        dplyr::mutate(
          INDIVIDUALS_REP = stringi::stri_join(INDIVIDUALS, "-R"),
          LANES_SHORT = rep("multiple", n()),
          BARCODES = rep("multiple", n()),
          FORWARD = rep("multiple", n()),
          REVERSE = rep("multiple", n()),
          FQ_FILES_F = stringi::stri_join(INDIVIDUALS_REP, ".1.fq.gz"),
          FQ_FILES_R = stringi::stri_join(INDIVIDUALS_REP, ".2.fq.gz"),
          REPLICATES = rep("R", n())
        )
    } else {
      replicates.sum <- replicates %>%
        dplyr::group_by(INDIVIDUALS) %>%
        dplyr::summarise(
          TOTAL = sum(TOTAL),
          NO_RADTAG = sum(NO_RADTAG),
          LOW_QUALITY = sum(LOW_QUALITY),
          RETAINED_READS = sum(RETAINED_READS),
          .groups = "keep"
        ) %>%
        dplyr::mutate(
          INDIVIDUALS_REP = stringi::stri_join(INDIVIDUALS, "-R"),
          LANES = rep("multiple", n()),
          LANES_SHORT = rep("multiple", n()),
          BARCODES = rep("multiple", n()),
          FQ_FILES = stringi::stri_join(INDIVIDUALS_REP, ".fq.gz"),
          REPLICATES = rep("R", n())
        )
    }


    replicate.individual <- unique(replicates$INDIVIDUALS)

    message("Replicate individuals: ", length(replicate.individual))

    # combine fq files of replicates ---------------------------------------------
    message("Combining replicates...")
    combine_replicates_fq <- function(replicate.individual, replicates) {
      paired <- FALSE # default

      # isolate the replicates
      rep.info <- dplyr::filter(replicates, INDIVIDUALS == replicate.individual)

      # isolate the first replicates
      rep1.fq <- dplyr::filter(rep.info, REPLICATES == 1)

      if (tibble::has_name(rep.info, "FQ_FILES_F")) paired <- TRUE

      # Create a new empty fq files
      if (paired) {
        # forward
        forward.name <- stringi::stri_join(rep1.fq$INDIVIDUALS, "-R.1.fq.gz")
        new.fq.f.path <- file.path(o, forward.name)
        file.create(new.fq.f.path)
        # reverse
        reverse.name <- stringi::stri_join(rep1.fq$INDIVIDUALS, "-R.2.fq.gz")
        new.fq.r.path <- file.path(o, reverse.name)
        file.create(new.fq.r.path)
      } else {
        rep1.fq.path <- stringi::stri_join(o, "/", rep1.fq$FQ_FILES)
        new.fq.path <- stringi::stri_join(o, "/", rep1.fq$INDIVIDUALS, "-R.fq.gz")
        file.create(new.fq.path)
      }



      # get the fq to combined in a vector
      combine_all_rep <- function(fq.files.to.combined, new.fq.path, output.path) {
        fq <- file.path(output.path, fq.files.to.combined)
        file.append(file1 = new.fq.path, file2 = fq)
      }


      if (paired) {
        # forward
        fq.files.forward.to.combined <- rep.info$FQ_FILES_F

        purrr::walk(.x = fq.files.forward.to.combined,
                    .f = combine_all_rep,
                    new.fq.path = new.fq.f.path,
                    output.path = output.path)

        # Reverse
        fq.files.reverse.to.combined <- rep.info$FQ_FILES_R
        purrr::walk(.x = fq.files.reverse.to.combined,
                    .f = combine_all_rep,
                    new.fq.path = new.fq.r.path,
                    output.path = output.path)

      } else {
        fq.files.to.combined <- rep.info$FQ_FILES

        purrr::walk(.x = fq.files.to.combined,
                    .f = combine_all_rep,
                    new.fq.path = new.fq.path,
                    output.path = output.path)
      }

    }

    purrr::walk(.x = replicate.individual,
                .f = combine_replicates_fq, replicates = replicates)

    # combine replicate sum with project.info
    # project.info.file$REPLICATES <- as.character(project.info.file$REPLICATES)

    project.info.file <- dplyr::bind_rows(project.info.file, replicates.sum) %>%
      dplyr::arrange(INDIVIDUALS, INDIVIDUALS_REP)

    message("In project info file:
    replicates have integer (e.g. id-1, id-2, etc.) appended to the original name
    combined replicates have \"-R\" appended to the original name")
  } else {
    message("Replicate individuals: 0")
  }

  # Add SQL_ID column for ustacks ----------------------------------------------
  project.info.file.sqlid <- project.info.file %>%
    dplyr::arrange(INDIVIDUALS_REP) %>%
    dplyr::mutate(
      SQL_ID = seq(1, n())
    )
  message("Updating info file with SQL ID")


  new.info.file.name <- stringi::stri_join("02_project_info/project.info.", file.date, ".tsv")
  readr::write_tsv(x = project.info.file.sqlid, file = new.info.file.name)
  message("New project info, see: ", new.info.file.name)
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  res <- list(project.info = project.info.file.sqlid,
              sample.per.lanes = sample.per.lanes)
  return(res)
} #End run_process_radtags
