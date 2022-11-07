<!-- badges: start -->
<!-- [![codecov](https://codecov.io/gh/stan-dev/projpred/branch/master/graph/badge.svg)](https://app.codecov.io/gh/stan-dev/projpred) -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/projpred?color=blue)](https://CRAN.R-project.org/package=projpred)
<!-- badges: end -->

# projpred [<img src="man/figures/logo.svg" align="right" height="139" alt="Stan Logo"/>](https://mc-stan.org)

The R package **projpred** performs the projection predictive variable selection
for various regression models. Usually, the reference model will be an
[**rstanarm**](https://mc-stan.org/rstanarm/) or
[**brms**](https://paul-buerkner.github.io/brms/) fit, but custom reference
models can also be used. Details on supported model types are given in section
["Supported types of reference
models"](https://mc-stan.org/projpred/articles/projpred.html#refmodtypes) of the
main vignette^[The main vignette can be accessed offline by typing
`vignette(topic = "projpred", package = "projpred")` or---more
conveniently---`browseVignettes("projpred")` within R.].

For details on how to cite **projpred**, see the [projpred citation
info](https://cran.r-project.org/web/packages/projpred/citation.html) on
CRAN.^[The citation information can be accessed offline by typing
`print(citation("projpred"), bibtex = TRUE)` within R.]. Further references
(including earlier work that **projpred** is based on) are given in section
["Introduction"](https://mc-stan.org/projpred/articles/projpred.html#introduction)
of the main vignette.

The [vignettes](https://mc-stan.org/projpred/articles/)^[The overview of all
vignettes can be accessed offline by typing `browseVignettes("projpred")` within
R.] (currently, the main vignette is the only one) illustrate how to use the
**projpred** functions in conjunction. Details on the **projpred** functions as
well as some shorter examples may be found in the
[documentation](https://mc-stan.org/projpred/reference/index.html)^[The
documentation can be accessed offline using `?` or `help()` within R.].

## Installation

There are two ways for installing **projpred**: from
[CRAN](https://CRAN.R-project.org/package=projpred) or from
[GitHub](https://github.com/stan-dev/projpred). The GitHub version might be more
recent than the CRAN version, but the CRAN version might be more stable.

### From CRAN

```r
install.packages("projpred")
```

### From GitHub

This requires the [**devtools**](https://devtools.r-lib.org/) package, so if
necessary, the following code will also install **devtools** (from
[CRAN](https://CRAN.R-project.org/package=devtools)):
```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("stan-dev/projpred", build_vignettes = TRUE)
```
To save time, you may omit `build_vignettes = TRUE`.
