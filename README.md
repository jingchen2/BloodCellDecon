
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BloodCellDecon

<!-- badges: start -->
<!-- badges: end -->

The goal of BloodCellDecon is to …

## Installation

You can install the development version of BloodCellDecon from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jingchen2/BloodCellDecon")
```

## Example

This is a basic example which shows you how to deconvolute test samples
in test.beta matrix:

``` r
library(BloodCellDecon)

test.res=estimateCellComposition(test.beta = test.beta, ref.beta.mat = ref.projection.EPIC$ref.beta.mat,projection = ref.projection.EPIC$projection, n.PC = 20,extended = F)
#> [1] "1800 probes found in test data."

p=plot_celltype(test.res*100,test.pd*100,celltype = 'Bmem')
#> `geom_smooth()` using formula = 'y ~ x'
```

<img src="man/figures/README-example-1.png" width="100%" />