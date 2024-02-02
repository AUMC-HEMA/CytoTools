#' @title Automated doublet gating
#' 
#' @description \code{gateDoublets} identifies doublets using three methods
#' 
#' @usage gateDoublets(ff)
#' 
#' @param ff A flowframe
#' @param channel1 Name of FSC-A channel (default = "FSC-A")
#' @param channel2 Name of FSC-H channel (default = "FSC-H")
#' @param predictionLevel Value of prediction level used (flowStats)
#' @param nMAD Number of median absolute deviations (PeacoQC)
#' @param nBins Number of bins (Splines)
#' @param nSD Number of standard deviations (Splines)
#' 
#' @return A matrix with a binary (0, 1) column for every method indicating
#' whether a cell is a doublet (1) or singlet (0)
#' 
#' @importFrom flowStats singletGate
#' @importFrom flowCore filter
#' @export
gateDoublets <- function(ff, channel1 = "FSC-A", channel2 = "FSC-H", 
                         predictionLevel = 0.99, nMAD = 4, nBins = 10, nSD = 2){
  # Gate doublets with flowStats singletGate
  boundaries <- flowStats::singletGate(ff, 
                                       area = channel1,
                                       height = channel2,
                                       prediction_level = predictionLevel)
  flowStats_doublet <- !flowCore::filter(x=ff, filter=boundaries)@subSet
  flowStats_doublet <- as.matrix(as.integer(flowStats_doublet))
  colnames(flowStats_doublet) <- 'flowStats_doublet'
  
  # Gate doublets with PeacoQC
  # Code was modified to save memory & compute (PeacoQC always returns flowframe) 
  # https://github.com/saeyslab/PeacoQC/blob/master/R/PeacoQC.R
  # Calculate the ratios
  ratio <- flowCore::exprs(ff)[, channel1] /
    (1+ flowCore::exprs(ff)[, channel2])
  # Define the region that is accepted
  r <- stats::median(ratio)
  r_m <- stats::mad(ratio)
  # Make selection
  selection <- ratio < r + nMAD * r_m
  PeacoQC_doublet <- !selection
  PeacoQC_doublet <- as.matrix(as.integer(PeacoQC_doublet))
  colnames(PeacoQC_doublet) <- 'PeacoQC_doublet'
  
  # Gate doublets using the spline-based classification of Duetz et al. (2021)
  fsc_a <- flowCore::exprs(ff)[, channel1]
  fsc_h <- flowCore::exprs(ff)[, channel2]
  bins <- cut(fsc_a, nBins)
  ratios <- fsc_h / fsc_a
  slope_per_bin <- tapply(ratios, bins, mean)
  expected_values <- fsc_a * slope_per_bin[bins]
  deviations <- abs(fsc_h - expected_values)
  x <- tapply(fsc_a, bins, mean)
  e <- tapply(expected_values, bins, mean)
  d <- tapply(deviations, bins, function(x){mean(x) + nSD * sd(x)})
  y <- e - d
  spl <- splinefun(x, y)
  selection <- fsc_h > spl(fsc_a)
  spline_doublet <- !selection
  spline_doublet <- as.matrix(as.integer(spline_doublet))
  colnames(spline_doublet) <- 'Spline_doublet'

  # Combine output
  output <- cbind(PeacoQC_doublet, flowStats_doublet, spline_doublet)
  return(output)
}
