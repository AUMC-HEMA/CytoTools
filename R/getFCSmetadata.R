#' @title Fast metadata extraction from a list of FCS files
#' 
#' @description \code{getFCSmetadata} uses flowCore to read the headers from FCS
#' files and returns a large dataframe with metadata for all files.
#' 
#' Note that this function can take long for very large sets of files. In that case,
#' run the function for subsets and concatenate the output!
#' 
#' @param files List of FCS files
#' @param commonParams List of common FCS parameters to save
#' @param ignore List of parameters from commonParams to ignore
#' 
#' @importFrom flowCore read.FCSheader
#' @export
getFCSmetadata <- function(files, commonParams = c("$FIL", "GUID", "EXPORT TIME",
                                                   "CYTNUM"),
                           ignore = NULL){
  headers <- flowCore::read.FCSheader(files)
  # Iterate over the headers and convert them to dataframes
  parameterList <- list()
  for (i in 1:length(headers)) {
    header <- headers[[i]]
    # First, identify all the markers that were used in this FCS
    # If markers were not measured, they are often excluded as $P_S parameter
    # but still have a $P_N and $P_V parameter
    pMarkers <- header[grepl("^\\$P[1-9][0-9]?[S]$", names(header))]
    pNames <- header[grepl("^\\$P[1-9][0-9]?[N]$", names(header))]
    pVoltages <- header[grepl("^\\$P[1-9][0-9]?[V]$", names(header))]
    # Keep only the values for markers
    pNames <- pNames[as.character(lapply(names(pMarkers), function(x) gsub("S", "N", x)))]
    pVoltages <- pVoltages[as.character(lapply(names(pMarkers), function(x) gsub("S", "V", x)))]
    # Add all the variables to a dataframe
    parameters <- list()
    for (parameter in commonParams){
      if (!is.null(ignore)){
        if (parameter %in% ignore){
          next
        }
      }
      parameters[parameter] <- header[[parameter]]
    }
    parameters['channel_string'] <- paste(pNames, collapse = "_")
    parameters['marker_string'] <- paste(pMarkers, collapse = "_")
    parameters['voltage_string'] <- paste(pVoltages, collapse = "_")
    parameterList[[files[i]]] <- parameters
  }
  # Convert list of lists to a dataframe
  output <- do.call(rbind, lapply(names(parameterList), function(name) {
    data.frame(Name = name, do.call(cbind, parameterList[[name]]), check.names = FALSE)
  }))  
  return(output)
}
