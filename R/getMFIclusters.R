#' @title Label peaks in a numerical vector of MFIs across clusters
#' 
#' @description \code{labelPeaks} assesses the MFIs of a given marker across
#' a range of clusters. The number of peaks (max 3) are detected and labeled
#' as either "unimodal" (1 peak), "low/high" (2 peaks) or "low/medium/high".
#' 
#' @usage labelPeaks(values)
#' 
#' @param values Vector of MFIs for a marker across clusters
#' 
#' @return A named list with cutoffs and labels for all clusters
#' 
#' @importFrom pracma findpeaks
#' @export
labelPeaks <- function(values){
  # Find peaks
  dens <- density(values)
  # Identify peak patterns (maximum 3 per channel)
  peaks <- pracma::findpeaks(dens$y, npeaks = 3)
  if (nrow(peaks) == 1){
    cutoffs <- c()
    labels <- rep("unimodal", length(values))
  } else if (nrow(peaks) == 2){
    cutoffs <- c(dens$x[peaks[1, 4]])
    labels <- ifelse(values > cutoffs[1], "high", "low")
  } else if (nrow(peaks) == 3){
    cutoffs <- c(dens$x[peaks[1, 4]], dens$x[peaks[2, 4]])
    labels <- cut(values, c(-Inf, cutoffs, Inf), labels = c("low", "medium", "high"),
                  include.lowest = TRUE)
  }
  output <- list("cutoffs" = cutoffs,
                 "labels" = labels)
}


#' @title Perform metaclustering based on MFI patterns across clusters
#' 
#' @description \code{getMFIclusters} assesses the MFIs of a given marker across
#' a range of clusters. The number of peaks (max 3) are detected and labeled
#' as either "unimodal" (1 peak), "low/high" (2 peaks) or "low/medium/high". Metaclusters
#' are formed by aggregating patterns across all markers for all clusters.
#' 
#' @usage getMFIclusters(MFIData)
#' 
#' @param MFIData Matrix of MFIs of every marker (column) on every cluster (row)
#' 
#' @return A named list containing the patterns, cutoffs and metaclustering
#' 
#' @export
getMFIclusters <- function(MFIData){
  peakLabels <- data.frame(matrix(ncol = ncol(MFIData), nrow = nrow(MFIData)))
  colnames(peakLabels) <- colnames(MFIData)
  peakCutoffs <- list()
  for (parameter in colnames(MFIData)){
    peakData <- labelPeaks(MFIData[, parameter])
    peakLabels[, parameter] <- peakData$labels
    peakCutoffs[[parameter]] <- peakData$cutoffs
  }
  # Create a representation containing the marker and the peak (e.g., "CD4high")
  markerStrings <- peakLabels
  for (i in seq_along(colnames(MFIData))){
    markerStrings[, i] <- paste0(colnames(MFIData)[i], peakLabels[, i])
  }
  # Create a representation of all marker strings (e.g., "CD3high_CD4high_CD8low")
  fullStrings <- apply(markerStrings, 1, function(row) paste(row, collapse = "_"))
  # Combine the unique strings into metaclusters
  metaclustering <- as.integer(factor(fullStrings, levels = unique(fullStrings)))

  # Format final output
  peakLabels$cluster <- rownames(peakLabels)
  peakLabels$metacluster <- metaclustering
  output <- list("data" = peakLabels,
                 "cutoffs" = peakCutoffs,
                 "metaclustering" = metaclustering)
  return(output)
}
