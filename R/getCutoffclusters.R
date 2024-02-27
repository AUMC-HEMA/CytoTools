#' @title Perform metaclustering based on manually identified MFI cutoffs
#' 
#' @description \code{getCutoffclusters} assesses the MFIs of a given marker across
#' a range of clusters. Using a named list of cutoffs, clusters are assigned
#' as either marker+ or marker-. This labeling is converted into a metaclustering.
#' 
#' @usage getCutoffclusters(MFIData, cutoffs)
#' 
#' @param MFIData Matrix of MFIs of every marker (column) on every cluster (row)
#' @param cutoffs Named list of parameters (corresponding to MFIData columns) and cutoffs
#' 
#' @return A named list containing the patterns, cutoffs and metaclustering
#' 
#' @export
getCutoffclusters <- function(MFIData, cutoffs){
  peakLabels <- data.frame(matrix(ncol = length(cutoffs), nrow = nrow(MFIData)))
  colnames(peakLabels) <- names(cutoffs)
  for (parameter in names(cutoffs)){
    labels <- ifelse(MFIData[, parameter] > cutoffs[[parameter]], "+", "-")
    peakLabels[, parameter] <- paste0(parameter, labels)
  }
  # Create a representation of all marker strings (e.g., "CD3high_CD4high_CD8low")
  fullStrings <- apply(peakLabels, 1, function(row) paste(row, collapse = "/"))
  peakLabels$string <- fullStrings
  # Combine the unique strings into metaclusters
  metaclustering <- as.integer(factor(fullStrings, levels = unique(fullStrings)))
  # Format final output
  peakLabels$cluster <- rownames(peakLabels)
  peakLabels$metacluster <- metaclustering
  output <- list("data" = peakLabels,
                 "cutoffs" = cutoffs,
                 "metaclustering" = metaclustering)
  return(output)
}
