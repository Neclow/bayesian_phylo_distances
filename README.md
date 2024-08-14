# Bayesian distance-based phylogenetics for the genomics era

Code for "Bayesian distance-based phylogenetics for the genomics era" (under review)

## Installation

### Dependencies

R version: 4.4.1
```r
require(devtools)
install.packages(
    c(
        "ape",
        "dplyr",
        "doParallel",
        "GGally",
        "magrittr",
        "phangorn",
        "Rcpp",
        "TreeSim",
    )
)

install_github("sebastianduchene/NELSI")
```
