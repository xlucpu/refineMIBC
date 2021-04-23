# refineMIBC

<!-- badges: start -->
<!-- badges: end -->

This package implements a 120-gene template, that assigns subtype labels according to the multi-omics consensus ensemble of muscle-invasive bladder cancer (MIBC) using nearest template prediction. The consensus ensemble identifies 4 integrative consensus subtypes: basal-inflamed, basal-noninflamed, luminal-excluded, and luminal-desert. This package also deploys a 5-immune-gene classifier to refine each basal-like MIBC as either basal-inflamed or basal-noninflamed by random forest classifier if basal-like classification has been already identified by other approaches (e.g., CMS, PAM, oneNN, Lund, etc.).

## Citation

For now, you can cite the following bioRxiv preprint

## Installation

You may install this package with:

``` r
if (!require("devtools")) 
    install.packages("devtools")
devtools::install_github("xlucpu/refineMIBC")
```

## Example

``` r
# load R package and internal data set
library(refineMIBC)
load(system.file("extdata", "demo.RData", package = "refineMIBC", mustWork = TRUE)) # load example data

# extract expression and TCGA subtype of MIBC
expr <- demo$MIBC.expr
subt <- demo$MIBC.subt

# refine MIBC into four subtypes
iCS  <- refineMIBC(expr       = expr, # expression data of MIBC
                   scaleFlag  = TRUE, # scale the data
                   centerFlag = TRUE) # center the data

head(iCS)
#                      prediction d.basal-inflamed d.basal-noninflamed d.luminal-desert d.luminal-excluded    p.value         FDR
# TCGA-HQ-A5NE-01A basal-inflamed        0.5413788           0.6479540        0.6007494          0.6973374 0.02797203 0.032200358
# TCGA-FD-A5BT-01A basal-inflamed        0.5420301           0.8216466        0.8498781          0.6313805 0.00100000 0.001455882
# TCGA-XF-AAMW-01A basal-inflamed        0.5475961           0.6043551        0.7397574          0.8120306 0.00100000 0.001455882
# TCGA-K4-A6FZ-01A basal-inflamed        0.5055447           0.6362818        0.5627713          0.8083026 0.00100000 0.001455882
# TCGA-BT-A42F-01A basal-inflamed        0.5338917           0.7166934        0.8419884          0.8426489 0.00100000 0.001455882
# TCGA-UY-A8OB-01A basal-inflamed        0.4858243           0.7400020        0.8302346          0.8539241 0.00100000 0.001455882

# refine basal-like MIBC into either basal-inflamed or basal-noninflamed
isbasal <- subt$TCGA == "Basal_squamous"
iBasal  <- refineBasal(expr    = expr,    # expression data of MIBC
                       isBasal = isbasal) # basal-like indicator
head(iBasal)
#              samID   prob       basal
# 1 TCGA-HQ-A5NE-01A 0.6500    inflamed
# 2 TCGA-FD-A5BT-01A 1.0000    inflamed
# 3 TCGA-XF-AAMW-01A 0.4545 noninflamed
# 4 TCGA-K4-A6FZ-01A 0.1965 noninflamed
# 5 TCGA-BT-A42F-01A 0.6350    inflamed
# 6 TCGA-UY-A8OB-01A 0.5295    inflamed
```
