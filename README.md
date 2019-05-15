
# SMARTp

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/SMARTp)](https://cran.r-project.org/package=SMARTp) -->

<!-- [![Travis build status](https://travis-ci.org/SMARTp/SMARTp.svg?branch=master)](https://travis-ci.org/SMARTp/SMARTp) -->

<!-- [![Codecov test coverage](https://codecov.io/gh/SMARTp/SMARTp/branch/master/graph/badge.svg)](https://codecov.io/gh/SMARTp/SMARTp?branch=master) -->

<!-- badges: end -->

## Overview

A SMART design for non-surgical treatments of chronic periodontitis with
spatially-referenced and non-randomly missing skewed outcomes.

## Installation

``` r
# Install from CRAN (when available)
install.packages("SMARTp")
# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("bandyopd/SMARTp")
```

## Usage

`library(SMARTp)` will load the following functions:

  - **CAR\_cov\_teeth**, for generating the CAR structure.
  - **MC\_var\_yibar\_mis**, for estimating mean and variance of the
    average change in CAL for each patient by Monte Carlo method.
  - **SampleSize\_SMARTp**, for sample size calculation under a
    clustered SMART design for chronic periodontitis.

See `?SampleSize_SMARTp` for a complete example of how to use this
package.

## Contact

Dipankar Bandyopadhyay, PhD, Professor: <dbandyop@vcu.edu>
