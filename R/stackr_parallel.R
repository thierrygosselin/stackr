# mclappply_win ----------------------------------------------------------------
#' @name mclapply_win
#' @title hack to switch function for parallel computation based on OS
#' @description \code{parallel::mclapply} doesn't work on Windows,
#' because forking is not supported.
#' This function defines a socket version of mclapply for windows computer
#' An implementation that switch automatically the parallel process when detecting
#' the os.
#' The code below was inspired from 
#' \pkg{parallel} \code{\link{mclapply}},
#' \href{https://github.com/nathanvan}{Nathan VanHoudnos},
#' \href{https://github.com/kvnkuang/pbmcapply}{Kevin Kuang},
#' \href{https://github.com/psolymos/pbapply}{Peter Solymos} and 
#' \href{https://github.com/EricArcher/}{Eric Archer}.


# @inheritParams parallel::mclapply
# Doesnt work and throws an error for bad markup so have to do it manually until
# parallel fix this bug
#' @param X see \pkg{parallel} \code{\link{mclapply}}
#' @param FUN see \pkg{parallel} \code{\link{mclapply}}
#' @param ... see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.preschedule see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.set.seed see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.silent see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cores see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cleanup see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.allow.recursive see \pkg{parallel} \code{\link{mclapply}}

# @return For mclapply, a list of the same length as X and named by X.
#' @importFrom utils sessionInfo
#' @importFrom parallel detectCores makeCluster clusterExport mclapply parLapply stopCluster
#' @rdname stackr_parallel
#' @export
#' @keywords internal

mclapply_win <- function(
  X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
  mc.silent = FALSE, mc.cores = NULL, mc.cleanup = TRUE, mc.allow.recursive = TRUE
) {
  
  # Create a cluster
  if (is.null(mc.cores)) {
    mc.cores <- parallel::detectCores() - 1
  }
  cl <- parallel::makeCluster(mc.cores)
  
  # We need to find the names of the loaded packages and export them to cluster
  tryCatch(
    {
      loaded.packages <- c(
        utils::sessionInfo()$basePkgs, #Base packages
        names(utils::sessionInfo()$otherPkgs) #Additional packages
      )
      
      #Export the packages to the clusters
      parallel::clusterExport(cl, 'loaded.packages', envir = environment())
      
      # Load the libraries on all the clusters
      parallel::parLapply(
        cl, 1:length(cl), function(xx){
          lapply(loaded.packages, function(yy) {
            require(yy , character.only = TRUE)})
        }
      )
      
      # We want the enclosing environment, not the calling environment
      cluster_export <- function(cl, FUN) {
        env <- environment(FUN)
        while (!identical(env, globalenv())) {
          env <- parent.env(env)
          parallel::clusterExport(cl, ls(all.names = TRUE, envir = env), envir = env)
        }
        parallel::clusterExport(cl, ls(all.names = TRUE, envir = env), envir = env)
      } # End cluster_export
      
      cluster_export(cl, FUN)
      
      # Run the lapply in parallel, with a special case for the ... arguments
      if (length(list(...)) == 0) {
        return(parallel::parLapply(cl = cl, X = X, fun = FUN))
      } else {
        return(parallel::parLapply(cl = cl, X = X, fun = FUN, ...))
      }
    }, finally = {
      parallel::stopCluster(cl) #Stop the cluster
    }
  )#End tryCatch
}#End mclapply_win


# mclapply with progress bar ---------------------------------------------------
#' @title mclapply_progress_bar
#' @description progress bar
#' @param X see \pkg{parallel} \code{\link{mclapply}}
#' @param FUN see \pkg{parallel} \code{\link{mclapply}}
#' @param ... see \pkg{parallel} \code{\link{mclapply}}
# @param mc.preschedule see \pkg{parallel} \code{\link{mclapply}}
# @param mc.set.seed see \pkg{parallel} \code{\link{mclapply}}
# @param mc.silent see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cores see \pkg{parallel} \code{\link{mclapply}}
# @param mc.cleanup see \pkg{parallel} \code{\link{mclapply}}
# @param mc.allow.recursive see \pkg{parallel} \code{\link{mclapply}}
# @param style see \pkg{utils} \code{\link{txtProgressBar}}
#' @rdname mclapply_progress_bar
#' @export
#' @keywords internal
#' @importFrom future plan futureCall value multiprocess
#' @importFrom parallel mclapply
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom pryr named_dots dots
PORT = 6011
mclapply_progress_bar <- function(
  X, FUN, ..., mc.cores = parallel::detectCores() - 1, ignore.interactive = FALSE
) {
  # mc.preschedule = TRUE, mc.set.seed = TRUE,
  # mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE,
  # mc.allow.recursive = TRUE  
  
  # new plan
  first.plan <- future::plan(strategy = "list")
  on.exit(future::plan(strategy = first.plan))
  future::plan(future::multiprocess)
  
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }
  
  # Account for interactive mode
  if (!interactive() & !ignore.interactive) {
    parallel::mclapply(X, FUN, ..., mc.cores = mc.cores)
  }
  
  monitoring <- future::futureCall(function(X, FUN, ..., mc.cores) {
    socket.server <- socketConnection(
      open = "wb", port = PORT, blocking = TRUE, server = TRUE)
    tryCatch(
      result <- parallel::mclapply(X, function(...) {
        res <- FUN(...)
        writeBin(1, socket.server)
        return(res)
      }, ..., mc.cores = mc.cores),
      finally = {close(socket.server)}
    )
    return(result)
  },
  args = list(X, FUN, ..., mc.cores = mc.cores),
  envir = environment(),
  # globals = TRUE
  globals = list(PORT = PORT)
  )
  
  length <- length(X)
  pb <- utils::txtProgressBar(min = 0, max = length, style = 3)
  utils::setTxtProgressBar(pb, 0)
  progress <- 0
  
  # Creating socket client and updating pb
  socket.active <- FALSE
  while (!socket.active) {
    Sys.sleep(0.5)
    try(socket.connection <- socketConnection(
      open = "rb", port = PORT, blocking = TRUE, server = FALSE)
      , silent = TRUE)
    if (exists("socket.connection")) {
      socket.active <- TRUE
      while (progress < length) {
        readBin(socket.connection, "double")
        progress <- progress + 1
        utils::setTxtProgressBar(pb, progress)
      }
      close(socket.connection)
    }
  }
  
  return(future::value(monitoring))
}

# .stackr_parallel--------------------------------------------------------------
# Overwrite the serial version of mclapply on Windows only
# @name .stackr_parallel
# @title Enable parallel execution on Windows
# @description Internal hack to enable parallel execution of \pkg{assigner}
#' functions on Windows.
# @inheritParams parallel::mclapply
#' @return For mclapply, a list of the same length as X and named by X.
# @importFrom parallel detectCores makeCluster clusterExport mclapply parLapply stopCluster
# @importFrom pbmcapply pbmclapply
#' @rdname stackr_parallel
#' @keywords internal
#' @export
.stackr_parallel <- switch(
  Sys.info()[['sysname']],
  Windows = {mclapply_win},
  Linux   = {mclapply_progress_bar},
  # Linux   = {parallel::mclapply},
  # Linux   = {pbmcapply::pbmclapply},
  Darwin  = {mclapply_progress_bar}
  # Darwin  = {parallel::mclapply}
  # Darwin  = {pbmcapply::pbmclapply}
)

