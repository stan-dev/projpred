[<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" width=100 alt="Stan Logo"/>](http://mc-stan.org)

# projpred

[![Build Status](https://travis-ci.org/stan-dev/projpred.svg?branch=master)](https://travis-ci.org/stan-dev/projpred)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/projpred?color=blue)](http://cran.r-project.org/web/packages/projpred)

An R package to perform projection predictive variable selection for generalized linear models. Compatible with [rstanarm][] but also other reference models can be used. 

The method is described and evaluated in comparison to many other methods in Piironen and Vehtari (2017). 

Currently, the supported models (family objects in R) include Gaussian, Binomial and Poisson families. See the [quickstart-vignette][] for examples.


### Resources

* [mc-stan.org/projpred](http://mc-stan.org/projpred) (online documentation, vignettes)
* [Ask a question](http://discourse.mc-stan.org) (Stan Forums on Discourse)
* [Open an issue](https://github.com/stan-dev/projpred/issues) (GitHub issues for bug reports, feature requests)


### Installation

* Install the latest release from CRAN:

```r
install.packages("rstantools")
```

* Install latest development version from GitHub (requires [devtools](https://github.com/r-lib/devtools) package):

```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
devtools::install_github('stan-dev/projpred', build_vignettes = TRUE)
```
    
### Example

```R
rm(list=ls())
library(projpred)
library(rstanarm)
options(mc.cores = parallel::detectCores())
set.seed(1)

# Gaussian and Binomial examples from the glmnet-package
data('df_gaussian', package = 'projpred')
#data('df_binom', package = 'projpred')

# fit the full model with a sparsifying prior
fit <- stan_glm(y ~ x, family = gaussian(), data = df_gaussian,
                prior = hs(df = 1, global_scale=0.01), iter = 500, seed = 1)
#fit <- stan_glm(y ~ x, family = binomial(), data = df_binom
#                prior = hs(df = 1, global_scale=0.01), iter = 500, seed = 1)


# perform the variable selection
fit_v <- varsel(fit)

# print the results
varsel_stats(fit_v)

# project the parameters for model sizes nv = 3,5 variables 
projs <- project(fit_v, nv = c(3, 5))

# predict using only the 5 most relevant variables
pred <- proj_linpred(fit_v, xnew=df_gaussian$x, nv=5, integrated=T)

# perform cross-validation for the variable selection
fit_cv <- cv_varsel(fit, cv_method='LOO')

# plot the validation results 
varsel_plot(fit_cv)
```


### References

Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative models via an entropic explanatory power. *Journal of Statistical Planning and Inference*, 111(1-2):77–94.

Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear models: a Bayesian approach via Kullback–Leibler projections. *Biometrika*, 85(1):29–37.

Juho Piironen and Aki Vehtari (2017). Comparison of Bayesian predictive methods for model selection. *Statistics and Computing*, 27(3):711-735. doi:10.1007/s11222-016-9649-y. ([online][piironenvehtari]).


  [rstanarm]: https://github.com/stan-dev/rstanarm
  [piironenvehtari]: https://link.springer.com/article/10.1007/s11222-016-9649-y
  [quickstart-vignette]: https://htmlpreview.github.io/?https://github.com/stan-dev/projpred/blob/master/vignettes/quickstart.html

