#' @title Match individually exported gated populations to a full FCS file
#' 
#' @description \code{matchFCSexports} uses a common variable between a full
#' FCS file and individually exported FCS files to re-construct a gating hierarchy.
#' 
#' @param file The original full FCS file
#' @param exports List of FCS files containing exported populations
#' @param labels How to label every export in the output dataframe. By default,
#' uses the full path of every file
#' @param matchVar The variable in each FCS file to use for matching
#' 
#' @return A dataframe indicating for every event whether its present in each
#' exported FCS file
#' 
#' @importFrom flowCore read.FCS
#' @export
matchFCSexports <- function(file, exports, labels = NULL, matchVar = "event_ID"){
  ff <- flowCore::read.FCS(file, truncate_max_range=FALSE)
  if (!matchVar %in% colnames(ff@exprs)){
    stop(matchVar, " not in ", file)
  }
  output <- data.frame(matchVar = ff@exprs[, matchVar])
  colnames(output) <- matchVar
  
  if (is.null(labels)){
    labels <- exports
  }
  
  for (i in seq_along(exports)) {
    export <- flowCore::read.FCS(exports[[i]], truncate_max_range = FALSE)
    if (!matchVar %in% colnames(export@exprs)){
      stop(matchVar, " not in ", export)
    }
    
    output[, labels[i]] <- ifelse(output[, matchVar] %in% export@exprs[, matchVar], 1, 0)

    # Check if the exports contain any events not present in the full file
    newCells <- ifelse(export@exprs[, matchVar] %in% output[, matchVar], 1, 0)
    if (sum(newCells == 0) > 1){
      stop("Found ", sum(newCells == 0), " events in ", exports[[i]], " not in ", file)
    }
    
    if (sum(output[, labels[[i]]]) != nrow(export@exprs)){
      stop("Mismatch merged events and export length for ", exports[[i]])
    }
  }
  return(output)
}
