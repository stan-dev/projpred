[<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" width=100 alt="Stan Logo"/>](https://mc-stan.org)

# projpred

[![Build Status](https://travis-ci.org/stan-dev/projpred.svg?branch=master)](https://travis-ci.org/stan-dev/projpred)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/projpred?color=blue)](https://cran.r-project.org/package=projpred)

An R package to perform projection predictive variable selection for generalized linear models. Compatible with [rstanarm][] and [brms][] but other reference models can also be used. 

The method is described in detail in Piironen et al. (2020) and Catalina et al. (2020) and evaluated in comparison to many other methods in Piironen and Vehtari (2017). 

Currently, the supported models (family objects in R) include Gaussian, Binomial and Poisson families. See the [quickstart-vignette][] for examples on GLMs and [quickstart-glmms-vignette][] for examples on GLMMs.


### Resources

* [mc-stan.org/projpred](https://mc-stan.org/projpred) (online documentation, vignettes)
* [Ask a question](https://discourse.mc-stan.org) (Stan Forums on Discourse)
* [Open an issue](https://github.com/stan-dev/projpred/issues) (GitHub issues for bug reports, feature requests)


### Installation

* Install the latest release from CRAN:

```r
install.packages('projpred')
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

split_structure <- break_up_matrix_term(y ~ x, data = df_gaussian)
df_gaussian <- split_structure$data
formula <- split_structure$formula

# fit the full model with a sparsifying prior
fit <- stan_glm(formula, family = gaussian(), data = df_gaussian,
                prior = hs(df = 1, global_scale=0.01), iter = 500, seed = 1)
#fit <- stan_glm(formula, family = binomial(), data = df_binom
#                prior = hs(df = 1, global_scale=0.01), iter = 500, seed = 1)


# perform the variable selection
vs <- varsel(fit)

# print the results
summary(vs)

# project the parameters for model sizes nterms = 3,5 variables 
projs <- project(vs, nterms = c(3, 5))

# predict using only the 5 most relevant variables
pred <- proj_linpred(vs, newdata=df_gaussian, nterms=5, integrated=TRUE)

# perform cross-validation for the variable selection
cvs <- cv_varsel(fit, cv_method='LOO')

# plot the validation results 
plot(cvs)
```


### References

Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative models via an entropic explanatory power. *Journal of Statistical Planning and Inference*, 111(1-2):77–94.

Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear models: a Bayesian approach via Kullback–Leibler projections. *Biometrika*, 85(1):29–37.

Piironen, Juho and Vehtari, Aki (2017). Comparison of Bayesian predictive methods for model selection. *Statistics and Computing*, 27(3):711-735. doi:10.1007/s11222-016-9649-y. ([online][piironenvehtari]).

Piironen, Juho, Paasiniemi, Markus and Vehtari, Aki (2020). Projective inference in high-dimensional problems: prediction and feature selection. *Electronic Journal of Statistics*, 14(1): 2155-2197 ([Online][projpred]).

Catalina, Alejandro, Bürkner, Paul-Christian, Vehtari, Aki (2020). Projection Predictive Inference for Generalized Linear and Additive Multilevel Models. ([arXiv:2010.06994][new-projpred]).


  [rstanarm]: https://github.com/stan-dev/rstanarm
  [brms]: https://github.com/paul-buerkner/brms
  [piironenvehtari]: https://link.springer.com/article/10.1007/s11222-016-9649-y
  [projpred]: https://projecteuclid.org/euclid.ejs/1589335310
  [new-projpred]: https://arxiv.org/abs/2010.06994
  [quickstart-vignette]: https://htmlpreview.github.io/?https://github.com/stan-dev/projpred/blob/master/vignettes/quickstart.html
  [quickstart-glmms-vignette]: https://htmlpreview.github.io/?https://github.com/stan-dev/projpred/blob/master/vignettes/quickstart_glmm.html

