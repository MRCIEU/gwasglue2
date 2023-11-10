
# gwasglue2

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/MRCIEU/gwasglue2/branch/main/graph/badge.svg)](https://app.codecov.io/gh/MRCIEU/gwasglue2?branch=main)
[![R-CMD-check](https://github.com/MRCIEU/gwasglue2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MRCIEU/gwasglue2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This package is part of the [OpenGWAS project](https://gwas.mrcieu.ac.uk). It aims to reduce friction between connecting GWAS summary data sources to a range of analytical tools. See the `Strategy` vignette for more information. 

It aims to replace the original [gwasglue](https://github.com/mrcieu/gwasglue) package but it is still in early development.

## Installation

You can install the development version of **gwasglue2** from [GitHub](https://github.com/) with:

```r
if(!require("remotes"))
{
    install.packages("remotes")
}
remotes::install_github("MRCIEU/gwasglue2")
```
