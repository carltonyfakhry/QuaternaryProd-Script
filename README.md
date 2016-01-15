
# QuaternaryProd-Script
A script which computes the Quaternary Product Scoring Statistic for signed and unsigned causal graphs using the *QuaternaryProd* R package.

## Installation of QuaternaryProd
Before running this script, make sure you have *QuaternaryProd* installed. To install *QuaternaryProd*, first make sure you have the latest version of *Rstudio*, *R* and the *devtools* package. Second, you can install the package by running the following in R-studio:
```{R}
library(devtools)
install_github("carltonyfakhry/QuaternaryProd", build_vignettes = TRUE, local = FALSE)
```
## Usage
For an introduction to the Quaternary Product Scoring Statistic and for an example on how to compute it over the publicly available network *Stringdb*, please see 
the *Vignette* for this package using the following:
```{R}
browseVignettes("QuaternaryProd")
```
## Running the script
To run the script over the e2f3_sig dataset run the following in your terminal:

Rscript cre4_run.R ./string.ents ./string.rels ./e2f3_sig.txt ./e2f3_sig_output.txt
