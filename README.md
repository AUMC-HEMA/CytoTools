# CytoTools
A tiny R package containing useful scripts for pre-processing cytometry data. For internal AUMC use. 



## Installation
You can install the package using devtools.

```r
# Install devtools if not available
install.packages("devtools")

devtools::install_github("AUMC-HEMA/CytoTools")
```



## Functionalities

**Quick metadata extraction (getFCSmetadata)**

*getFCSmetadata* is a wrapper function around flowCore's "read.FCSheader" that allows for a super quick extraction of all metadata (channel and marker names, cytometers, voltages, etc.) from a large set of FCS files. 



**Matching exported FCS populations (matchFCSexports)**

Because Infinicyt only allows for exporting gated populations (but does not support exporting gates), we have created a workflow to assign each cell an identifier which can be used to retro-actively assign cells to gates based on exported FCS files in the original FCS file ("*matchFCSexports*"). This allows for evaluating all gated populations together in one file. 



**Automated doublet gating (gateDoublets)**

A number of automated doublet gates have been incorporated in cytometry software. The *gateDoublets* function gates doublets using three different methods and returns a matrix with columns indicating whether a method considers that cell a doublet. 



**Optimized logicle transformations (getLogicle)**

Suboptimal compensation introduces negative events which are problematic for default log-transformation. To solve this the logicle transformation uses a linear transformation for a range of the data. The width of this range can be estimated automatically using flowCore. 

The *getLogicle* estimates the median width across a set of FCS files and returns a transformList object which can be used to transform a flowframe. 



**Perform metaclustering based on cluster MFIs**

The default consensus metaclustering of packages like FlowSOM is often suboptimal, and choices for metacluster counts are hard to determine. When using mostly lineage or cell type markers for clustering (e.g., CD3, CD4, CD8), it is relatively straightforward to separate populations into negative/dim/positive/bright subsets. 

The *getMFIclusters* function can use the MFIs of such markers across N clusters. For each marker, a density estimation is performed, which is used to detect peaks of MFIs across clusters (e.g., negative / positive). This results in a string representation for every cluster (e.g., "CD3low"). These are combined across markers for every cluster (e.g., "CD3high/CD4low/CD8high") and all clusters are grouped together if they have a similar string, resulting in a metaclustering. 

The benefit of this approach is that metaclusters have more unimodal marker distributions, while also being interpretable based on their label. 

```{r}
# Perform FlowSOM
fSOM <- FlowSOM(ff, xdim = 20, ydim = 20) # Large grid for better density patterns
MFIData <- GetClusterMFIs(fSOM, colsUsed = TRUE)

# Get new metaclustering
mc <- getMFIclusters(MFIData)$metaclustering

# Replace existing metaclusters in FlowSOM
fSOM <- UpdateMetaclusters(fSOM, clusterAssignment = mc)
```
