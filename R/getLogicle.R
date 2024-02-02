#' @title Estimation of linearization width for Logicle transformation across a set 
#' of FCS files
#' 
#' @description \code{estimateLogicles} reads FCS and compensates an FCS file
#' and subsequently estimates the linearization width parameter for a logicle 
#' transformation. Note that reading all FCS files can take some time for large 
#' sets of files.
#' 
#' @usage estimateLogicles(files)
#' 
#' @param files List of FCS file paths
#' 
#' @return A dataframe with estimated linearization width for every channel listed
#' in $SPILL slot
#' 
#' @importFrom flowCore read.FCS compensate estimateLogicle
#' @importFrom PeacoQC RemoveMargins
#' @export
estimateLogicles <- function(files){
  parameters <- list()
  for (file in files){
    ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)
    ff <- PeacoQC::RemoveMargins(ff, channels=c("FSC-A", "SSC-A", "FSC-H", "SSC-H"),
                           output="frame")
    ff <- flowCore::compensate(ff, ff@description$SPILL)
    
    # Estimate the linearization width using flowCore's estimateLogicle
    # Sometimes, this fails because of a low value of m
    flag <- FALSE
    tryCatch({
      estimate <- flowCore::estimateLogicle(ff, channels = colnames(ff@description$SPILL))
    }, error = function(e) {
      # Handle the error
      cat("Could not estimate logicle for", file, ":", conditionMessage(e), "\n")
      flag <- TRUE
    })
    
    widths <- list()
    for (channel in colnames(ff@description$SPILL)){
      logicleParameters <- as.list(environment(estimate@transforms[[channel]]@f))
      if (!flag){
        w <- logicleParameters[["w"]][[1]]
      } else {
        # Mark cases where estimateLogicle failed with NA
        w <- NA
      }
      widths[[channel]] <- w
    }
    parameters[[file]] <- widths
  }
  parameters <- do.call(rbind, lapply(names(parameters), function(name) {
    data.frame(file = name, do.call(cbind, parameters[[name]]), check.names = FALSE)
  }))
  rownames(parameters) <- parameters$file
  parameters$file <- NULL
  return(parameters)
}


#' @title Get an optimized logicle transformList object based on a set of FCS files
#' 
#' @description Uses \code{estimateLogicles} to determine the optimal linearization
#' width for every channel in a set of FCS files and builds a transformList object
#' using the median width for every channel
#' 
#' @usage getLogicle(files)
#' 
#' @param files List of FCS file paths
#' @param verbose Option for printing identified parameters
#' 
#' @return A flowCore transformList object
#' 
#' @importFrom flowCore logicleTransform transformList
#' 
#' @seealso \code{\link{estimateLogicles}}
#' 
#' @export
getLogicle <- function(files, verbose = TRUE){
  # Estimate linearization widths
  wParams <- estimateLogicles(files)
  transforms <- list()
  for (channel in colnames(wParams)){
    w <- median(wParams[, channel])
    if (verbose){
      message("Estimated median optimal linearization width for ", channel, " at ",
              round(w, 3))
    }
    transforms[[channel]] <- flowCore::logicleTransform(w = w, 
                                                        transformationId = channel)
    # For some reason, running the the line below is required to construct
    # a proper transformList??? I cannot figure out why...
    # If not called, the same width is used for all channels
    summary(transforms[[channel]])
  }
  tfList <- flowCore::transformList(colnames(wParams), transforms)
  return(tfList)
}
