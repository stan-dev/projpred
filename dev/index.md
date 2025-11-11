# projpred

The R package **projpred** performs the projection predictive variable
selection for various regression models. Usually, the reference model is
an [**rstanarm**](https://mc-stan.org/rstanarm/) or
[**brms**](https://paulbuerkner.com/brms/) fit, but custom reference
models can also be used. Details on supported model types are given in
section [“Supported types of
models”](https://mc-stan.org/projpred/articles/projpred.html#modtypes)
of the main vignette[¹](#fn1).

For details on how to cite **projpred**, see the [projpred citation
info](https://CRAN.R-project.org/package=projpred/citation.html) on
CRAN[²](#fn2). Further references (including earlier work that
**projpred** is based on) are given in section
[“Introduction”](https://mc-stan.org/projpred/articles/projpred.html#intro)
of the main vignette.

The [vignettes](https://mc-stan.org/projpred/articles/)[³](#fn3)
illustrate how to use the **projpred** functions in conjunction. Details
on the **projpred** functions as well as some shorter examples may be
found in the
[documentation](https://mc-stan.org/projpred/reference/index.html)[⁴](#fn4).

## Installation

There are two ways for installing **projpred**: from
[CRAN](https://CRAN.R-project.org/package=projpred) or from
[GitHub](https://github.com/stan-dev/projpred). The GitHub version might
be more recent than the CRAN version, but the CRAN version might be more
stable.

### From CRAN

``` r
install.packages("projpred")
```

### From GitHub

This requires the [**devtools**](https://devtools.r-lib.org/) package,
so if necessary, the following code will also install **devtools** (from
[CRAN](https://CRAN.R-project.org/package=devtools)):

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("stan-dev/projpred", build_vignettes = TRUE)
```

To save time, you may omit `build_vignettes = TRUE`.

## Contributing to projpred

We welcome contributions! The **projpred** package is under active
development. If you find bugs or have ideas for new features (for us or
yourself to implement) please [open an
issue](https://github.com/stan-dev/projpred/issues) on GitHub. See
[CONTRIBUTING.md](https://github.com/stan-dev/projpred/blob/master/.github/CONTRIBUTING.md)
for more details.

------------------------------------------------------------------------

1.  The main vignette can be accessed offline by typing
    [`vignette(topic = "projpred", package = "projpred")`](https://mc-stan.org/projpred/dev/articles/projpred.md)
    or—more conveniently—`browseVignettes("projpred")` within R.

2.  The citation information can be accessed offline by typing
    `print(citation("projpred"), bibtex = TRUE)` within R.

3.  The overview of all vignettes can be accessed offline by typing
    `browseVignettes("projpred")` within R.

4.  The documentation can be accessed offline using `?` or
    [`help()`](https://rdrr.io/r/utils/help.html) within R.
