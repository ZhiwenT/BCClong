
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BCCLong

<!-- badges: start -->
<!-- badges: end -->

The goal of BCCLong is to compute a Bayesian Consensus Clustering (BCC)
model for mixed-type longitudinal data

## Description

`BCClong` is an R package for performing Bayesian Consensus Clustering
(BCC) model for clustering continuous, discrete and categorical
longitudinal data, which are commonly seen in many clinical studies \[Lu
et al., 2021\]. Statistical methods for clustering a single longitudinal
trajectory have been well-developed and widely used in many different
medical research areas. However, it is very common these days to
encounter situations where several longitudinal markers or responses are
collected simultaneously in a study and there is a growing interest to
examine how multiple longitudinal characteristics could collectively
contribute to disaggregating disease heterogeneity.

## Installation

You can install the development version of BCC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ZhiwenT/BCClong")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BCClong)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

## References

-   [Lu, Z., & Lou, W. (2021). Bayesian consensus clustering for
    Multivariate Longitudinal Data. *Statistics in Medicine*, 41(1),
    108–127.](https://doi.org/10.1002/sim.9225)
