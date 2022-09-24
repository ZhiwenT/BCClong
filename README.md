
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BCCLong

<!-- badges: start -->
<!-- badges: end -->

The goal of BCCLong is to compute a Bayesian Consensus Clustering (BCC)
model for mixed-type longitudinal data

## Description

`BCClong` is an R package for performing Bayesian Consensus Clustering
(BCC) model for clustering continuous, discrete and categorical
longitudinal data, which are commonly seen in many clinical studies [Lu
et al., 2021](https://doi.org/10.1002/sim.9225). Statistical methods for
clustering a single longitudinal trajectory have been well-developed and
widely used in many different medical research areas. However, it is
very common these days to encounter situations where several
longitudinal markers or responses are collected simultaneously in a
study and there is a growing interest to examine how multiple
longitudinal characteristics could collectively contribute to
disaggregating disease heterogeneity.

## Installation

You can install the development version of BCC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ZhiwenT/BCClong", build_vignettes = TRUE)
library("MPLNClust")
```

## Overview

To list all the functions available in the package:

``` r
ls("package:BCClong")
```

## References

-   [Lu, Z., & Lou, W. (2021). Bayesian consensus clustering for
    Multivariate Longitudinal Data. *Statistics in Medicine*, 41(1),
    108–127.](https://doi.org/10.1002/sim.9225)

## Maintainer

-   Zhiwen Tan (<z.tan@queensu.ca>).

## Contributions

`BCClong` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/ZhiwenT/BCClong/issues).
