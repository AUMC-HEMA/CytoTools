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

**Matching exported FCS populations (matchFCSexports)**

Because Infinicyt only allows for exporting gated populations (but does not support exporting gates), we have created a workflow to assign each cell an identifier which can be used to retro-actively assign cells to gates based on exported FCS files in the original FCS file ("*matchFCSexports*"). This allows for evaluating all gated populations together in one file. 



**Automated doublet gating (gateDoublets)**

A number of automated doublet gates have been incorporated in cytometry software. The *gateDoublets* function gates doublets using three different methods and returns a matrix with columns indicating whether a method considers that cell a doublet. 



**Optimized logicle transformations (getLogicle)**

Suboptimal compensation introduces negative events which are problematic for default log-transformation. To solve this the logicle transformation uses a linear transformation for a range of the data. The width of this range can be estimated automatically using flowCore. 

The *getLogicle* estimates the median width across a set of FCS files and returns a transformList object which can be used to transform a flowframe. 
