
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BCCLong

<!-- badges: start -->
<!-- badges: end -->

The goal of BCCLong is to compute a Bayesian Consensus Clustering (BCC)
model for mixed-type longitudinal data

## Description

Statistical methods for clustering a single longitudinal trajectory have
been well-developed and widely used in many different medical research
areas. However, it is very common these days to encounter situations
where several longitudinal markers or responses are collected
simultaneously in a study and there is a growing interest to examine how
multiple longitudinal characteristics could collectively contribute to
disaggregating disease heterogeneity. Therefore, the `BCClong` package
has been created. `BCClong` is an R package for performing Bayesian
Consensus Clustering (BCC) model for clustering continuous, discrete and
categorical longitudinal data, which are commonly seen in many clinical
studies \[Lu et al., 2021\]<https://doi.org/10.1002/sim.9225>.

## Installation

You can install the development version of BCC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ZhiwenT/BCClong", build_vignettes = TRUE)
library("BCClong")
```

## Overview

To list all the functions available in the package:

``` r
ls("package:BCClong")
```

Currently, there are 5 function in this package which are
***BCC.multi***, ***BayesT***, ***model.selection.criteria***,
***traceplot***, ***trajplot***.

***BCC.multi*** function performs clustering on mixed-type (continuous,
discrete and categorical) longitudinal markers using Bayesian consensus
clustering method with MCMC sampling and provide a summary statistics
for the computed model. This function will take in a data set and
multiple parameters and output a BCC model with summary statistics.

***BayesT*** function assess the model goodness of fit by calculate the
discrepancy measure T(, ) with following steps (a) Generate T.obs based
on the MCMC samples (b) Generate T.rep based on the posterior
distribution of the parameters (c) Compare T.obs and T.rep, and
calculate the P values.

***model.selection.criteria*** function calculates DIC and WAIC for the
fitted model ***traceplot*** function visualize the MCMC chain for model
parameters ***trajplot*** function plot the longitudinal trajectory of
features by local and global clustering

more information can be found by using the code below after installation

``` r
?BCClong::BCC.multi
?BCClong::BayesT
?BCClong::model.selection.criteria
?BCClong::traceplot
?BCClong::trajplot
```

The package tree structure is provide below

``` r
- BCClong
  |- BCClong.Rproj
  |- DESCRIPTION
  |- NAMESPACE
  |- LICENSE
  |- README
  |- NEWS
  |- data
    |- conRes.rda
    |- epil.rda
    |- epil1.rda
    |- epil2.rda
    |- epil3.rda
    |- example.rda
    |- example1.rda
    |- PBCseqfit.rda
  |- inst
    |-CITATION
  |- man
    |- BayesT.Rd
    |- BCC.multi.Rd
    |- model.selection.criteria.Rd
    |- traceplot.Rd
    |- trajplot.Rd
    |- print.BCC.Rd
    |- plot.BCC.Rd
    |- summary.BCC.Rd
  |- R
    |- bcclong.R
    |- classMethods.R
    |- DiscrepancyMeasure.R
    |- modelSelection.R
    |- RcppExports.R
    |- Traceplot.R
    |- Trajplot.R
  |- src
    |- c_which.h
    |- c_which.cpp
    |- BCC.cpp
    |- Likelihood.cpp
    |- RcppExports.cpp
    |- Makevars
    |- Makevars.win
  |- vignettes
    |- ContinuousData.Rmd
    |- ContinuousData.html
    |- MixedTypeData.Rmd
    |- MixedTypeData.html
```

## Tutorials

For tutorials and plot interpretation, refer to the vignette:

``` r
browseVignettes("BCClong")
```

Three options include a HTMl version, source R markdown file and R code
file. There are two tutorials in this package, one is for dataset with
continuous data only, and the second one is for dataset with mixed type
of data. Tutorial can also be found from the link below. Make sure to
open the html file in browser, the github website only shows the source
code.

For multiple continuous longitudinal markers only:

<https://htmlpreview.github.io/?https://github.com/ZhiwenT/BCClong/blob/main/vignettes/ContinuousData.html>

For multiple mixed type longitudinal markers:

<https://htmlpreview.github.io/?https://github.com/ZhiwenT/BCClong/blob/main/vignettes/MixedTypeData.html>

## Citation for Package

``` r
citation("BCClong")
```

Tan, Z., Shen, C., Lu, Z. (2022) BCClong: an R package for performing
Bayesian Consensus Clustering model for clustering continuous, discrete
and categorical longitudinal data. URL
<https://github.com/ZhiwenT/BCClong>

## References

- \[Lu, Z., & Lou, W. (2021). Bayesian consensus clustering for
  Multivariate Longitudinal Data. *Statistics in Medicine*, 41(1),
  108–127.\]<https://doi.org/10.1002/sim.9225>

- \[Tan, Z., Shen, C., Subbarao, P., Lou, W. and Lu, Z., 2022. A Joint
  Modeling Approach for Clustering Mixed-Type Multivariate Longitudinal
  Data: Application to the CHILD Cohort Study. *arXiv
  preprint*.\]<https://doi.org/10.48550/arXiv.2210.08385>

## Maintainer

- Zhiwen Tan (<z.tan@queensu.ca>).

## Contributions

`BCClong` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/ZhiwenT/BCClong/issues).
