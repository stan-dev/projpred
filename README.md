# glmproj

An R package to perform projection predictive variable selection for generalized linear models fitted with [rstanarm][]. 

The package uses forward search starting from the empty submodel model, adds variables one at a time, each time choosing the variable that decreases the KL-divergence from the projection to the full model the most. 

Currently, the supported models (family objects in R) include Gaussian, Binomial and Poisson families.

Installation
------------

    devtools::install_github('paasim/glmproj')
    
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
    df1 <- data.frame(x = I(x), y = y)

    # fit the full model with a sparsifying prior
    fit <- stan_glm(y ~ x, gaussian(), df1, prior = hs(df = 3), iter = 500)
    #fit <- stan_glm(y ~ x, binomial(), df1, prior = hs(df = 4), iter = 500)


    # perform the variable selection
    # note that this may take some time for other GLMs than 
    # Gaussian with identity link
    fit_v <- varsel(fit, verbose = T)
    # print the results
    summary_varsel(fit_v)

    # project the parameters for a model of size 5 and 8
    fit_p <- predict_proj(fit_v, nv = c(5,8))
    coef(fit_p)

    # perform cross-validation for the variable selection
    # this takes some time to complete, especially for the non-Gaussian case.
    k_fold <- glmproj::kfold(fit, save_fits = T)
    fit_v <- cv_varsel(fit, k_fold, verbose = T)

    # plot the results
    plot_varsel(fit_v)

References
------------
Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative models via an entropic explanatory power. Journal of Statistical Planning and Inference, 111(1-2):77–94.

Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear models: a Bayesian approach via Kullback–Leibler projections. Biometrika, 85(1):29–37.

Juho Piironen and Aki Vehtari (2016). Comparison of Bayesian predictive methods for model selection. Statistics and Computing, ([online][piironenvehtari]).


  [rstanarm]: https://github.com/stan-dev/rstanarm
  [piironenvehtari]: https://link.springer.com/article/10.1007/s11222-016-9649-y

