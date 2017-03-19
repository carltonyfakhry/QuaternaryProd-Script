
# QuaternaryProd-Script
A script which computes the Quaternary Dot Product Scoring Statistic (or simply the Quaternary Statistic) for signed and unsigned causal graphs using the *QuaternaryProd* R package over the STRINGdb network.

## Installation of QuaternaryProd
Before running this script, make sure you have *QuaternaryProd* installed. To install *QuaternaryProd*, first make sure you have the latest version of *Rstudio*, *R* and the *devtools* package. Second, you can install the package by running the following in R-studio:
```{R}
library(devtools)
install_github("carltonyfakhry/QuaternaryProd", build_vignettes = TRUE, local = FALSE)
```
## Usage
For an introduction to the Quaternary Dot Product Scoring Statistic, please see 
the Vignette for *QuaternaryProd* using the following:
```{R}
browseVignettes("QuaternaryProd")
```
## Running the script
First, you need to install the *fdrtool* R package. Second, you need to download the files in this repository and change your directory to where the files have been downloaded. Finally, to execute the script over the *e2f3_sig* dataset and the STRINGdb network, run the following in your terminal:

*Rscript script.R ./string.ents ./string.rels ./e2f3_sig.txt ./e2f3_sig_output.txt*

The results will be outputed to *e2f3_sig_output.txt*. The same script can be executed with the other datasets *myc_sig* and *ras_sig*.

