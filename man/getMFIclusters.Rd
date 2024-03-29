% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMFIclusters.R
\name{getMFIclusters}
\alias{getMFIclusters}
\title{Perform metaclustering based on MFI patterns across clusters}
\usage{
getMFIclusters(MFIData)
}
\arguments{
\item{MFIData}{Matrix of MFIs of every marker (column) on every cluster (row)}
}
\value{
A named list containing the patterns, cutoffs and metaclustering
}
\description{
\code{getMFIclusters} assesses the MFIs of a given marker across
a range of clusters. The number of peaks (max 3) are detected and labeled
as either "unimodal" (1 peak), "low/high" (2 peaks) or "low/medium/high". Metaclusters
are formed by aggregating patterns across all markers for all clusters.
}
