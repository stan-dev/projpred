[![Travis-CI Build Status](https://travis-ci.org/stan-dev/projpred.svg?branch=master)](https://travis-ci.org/stan-dev/projpred)
[![codecov](https://codecov.io/gh/stan-dev/projpred/branch/master/graph/badge.svg)](https://codecov.io/gh/stan-dev/projpred)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/projpred?color=blue)](http://cran.r-project.org/web/packages/projpred)

<br>

<div style="text-align:left">
<span><a href="http://mc-stan.org">
<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" width=100 alt="Stan Logo"/> </a><h2><strong>projpred</strong></h2><h4>Projection predictive variable selection</h4></span>
</div>

<br>
The **projpred** package performs projection predictive variable selection for
generalized linear models. Currently **projpred** is most easily compatible with
[**rstanarm**](http://mc-stan.org/rstanarm) but other reference models can also
be used.

The method is described and evaluated in comparison to many other methods in

* Juho Piironen and Aki Vehtari (2017). Comparison of Bayesian predictive methods for model selection. *Statistics and Computing*, 27(3):711-735. doi:10.1007/s11222-016-9649-y. ([online](https://link.springer.com/article/10.1007/s11222-016-9649-y)).

Currently, the supported models (family objects in R) include Gaussian, Binomial
and Poisson families. See the **projpred**
[vignette](http://mc-stan.org/projpred/articles) for examples.



## Getting Started

To get started see the __projpred__ [vignettes](http://mc-stan.org/projpred/articles).


## Installation

Install the latest release from **CRAN**

```r
install.packages("projpred")
```

Install the latest development version from **GitHub**

```r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("stan-dev/projpred", build_vignettes = TRUE)
```
