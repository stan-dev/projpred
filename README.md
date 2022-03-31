<!-- badges: start -->
<!-- [![codecov](https://codecov.io/gh/stan-dev/projpred/branch/master/graph/badge.svg)](https://app.codecov.io/gh/stan-dev/projpred) -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/projpred?color=blue)](https://CRAN.R-project.org/package=projpred)
<!-- badges: end -->

# projpred [<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" align="right" width="120" alt="Stan Logo"/>](https://mc-stan.org)

The **projpred** R package performs the projection predictive variable selection
for generalized linear models as well as generalized linear and additive
multilevel models. The package is compatible with the
[**rstanarm**](https://mc-stan.org/rstanarm/) and
[**brms**](https://paul-buerkner.github.io/brms/) packages, but custom reference
models can also be used.

The projection predictive variable selection is based on the ideas of Goutis and
Robert (1998) and Dupuis and Robert (2003). The methods implemented in
**projpred** are described in detail in Piironen et al. (2020) and Catalina et
al. (2020). They are evaluated in comparison to many other methods in Piironen
and Vehtari (2017). Type `citation("projpred")` in R (or see the `CITATION`
file) for details on how to cite **projpred**.

Currently, the supported response distributions (objects of class `family` in R)
are `gaussian()`, `binomial()` (via the **brms** package, `brms::bernoulli()` is
also supported), and `poisson()`.

See the [vignette](https://mc-stan.org/projpred/articles/projpred.html) for an
example application. Details on **projpred**'s functions as well as some shorter
examples may be found in the documentation.

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

## References

Catalina, A., Bürkner, P.-C., and Vehtari, A. (2020). Projection predictive
inference for generalized linear and additive multilevel models.
*arXiv:2010.06994*. URL: <https://arxiv.org/abs/2010.06994>.

Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative models
via an entropic explanatory power. *Journal of Statistical Planning and
Inference*, **111**(1-2):77–94. DOI:
[10.1016/S0378-3758(02)00286-0](https://doi.org/10.1016/S0378-3758(02)00286-0).

Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear models:
A Bayesian approach via Kullback–Leibler projections. *Biometrika*,
**85**(1):29–37.

Piironen, J. and Vehtari, A. (2017). Comparison of Bayesian predictive methods
for model selection. *Statistics and Computing*, **27**(3):711-735. DOI:
[10.1007/s11222-016-9649-y](https://doi.org/10.1007/s11222-016-9649-y).

Piironen, J., Paasiniemi, M., and Vehtari, A. (2020). Projective inference in
high-dimensional problems: Prediction and feature selection. *Electronic Journal
of Statistics*, **14**(1):2155-2197. DOI:
[10.1214/20-EJS1711](https://doi.org/10.1214/20-EJS1711).
