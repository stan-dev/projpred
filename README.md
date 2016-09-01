# glmproj

An R package to perform projection predictive variable selection for generalized linear models fitted with [rstanarm][]. 

The package uses forward search as a search heuristic, that is, starting from the empty submodel model, adds variables one at a time, each time choosing the variable that decreases the KL-divergence to the full model the most. 

Currently, supported models include gaussian with identity link function, binomial with probit and logit link functions and poisson with log link function.

Installation
------------

    devtools::install_github('paasim/glmproj')
    
Usage
-----

The package provides the following functions:
* varsel
* cv\_varsel
* 

Example
-------

    rm(list=ls())
    library(glmproj)
    library(rstanarm)
    options(mc.cores = parallel::detectCores())
    set.seed(1)

    # Gaussian and Binomial examples from the glmnet-package
    data('QuickStartExample', package = 'glmnet')
    #data('BinomialExample', package = 'glmnet')
    df1 <- list(x = x, y = y)

    # fit the full model with a sparsifying prior
    fit <- stan\_glm(y ~ x, gaussian(), df1, prior = hs(df = 1))
    #fit <- stan_glm(y ~ x, binomial(), df1, prior = hs(df = 1))

    # perform the variable selection
    vars <- varsel(fit, verbose = T)
    # print the results
    vars

    # project the parameters for a model of size 5 and 8
    projection <- project(vars, fit, size = c(5,8))
    projection

    # perform cross-validation for the variable selection
    # this takes some time to complete, especially for the non-gaussian case.
    fits <- cv_fit(fit)
    cv_vars <- cv_varsel(fit, fits, verbose = T)

    # plot the results
    plot(cv_vars)

References
------------
Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative models via an entropic explanatory power. Journal of Statistical Planning and Inference, 111(1-2):77–94.

Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear models: a Bayesian approach via Kullback–Leibler projections. Biometrika, 85(1):29–37.

Juho Piironen and Aki Vehtari (2016). Comparison of Bayesian predictive methods for model selection. Statistics and Computing, ([online][piironenvehtari]).



  [rstanarm]: https://github.com/stan-dev/rstanarm
  [piironenvehtari]: https://link.springer.com/article/10.1007/s11222-016-9649-y

