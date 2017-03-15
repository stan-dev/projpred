# glmproj

An R package to perform projection predictive variable selection for generalized linear models fitted with [rstanarm][]. 

The method is described and evaluated in comparison to many other methods in Piironen and Vehtari (2016). 

Currently, the supported models (family objects in R) include Gaussian, Binomial and Poisson families. See the examples in the vignettes directory.

Installation
------------

	if (!require(devtools)) {
  		install.packages("devtools")
  		library(devtools)
	}
	devtools::install_github('paasim/glmproj', build_vignettes = TRUE)

    
Example
-------

    rm(list=ls())
    library(glmproj)
    library(rstanarm)
    options(mc.cores = parallel::detectCores())
    set.seed(1)

    # Gaussian and Binomial examples from the glmnet-package
    data('df_gaussian', package = 'glmproj')
    #data('df_binom', package = 'glmproj')

    # fit the full model with a sparsifying prior
    fit <- stan_glm(y ~ x, family = gaussian(), data = df_gaussian,
                    prior = hs(df = 1, global_scale=0.03), iter = 500, seed = 1)
    #fit <- stan_glm(y ~ x, family = binomial(), data = df_binom
    #                prior = hs(df = 1, global_scale=0.03), iter = 500, seed = 1)


    # perform the variable selection
    fit_v <- varsel(fit)
    
    # print the results
    varsel_statistics(fit_v)

    # project the parameters for model sizes nv = 3,5 variables 
    fit_p <- project(fit_v, nv = c(3, 5))
    proj_coef(fit_p)
    
    # perform cross-validation for the variable selection
    fit_cv <- cv_varsel(fit, cv_method='LOO')

    # plot the results
    varsel_plot(fit_cv)

References
------------
Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative models via an entropic explanatory power. Journal of Statistical Planning and Inference, 111(1-2):77–94.

Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear models: a Bayesian approach via Kullback–Leibler projections. Biometrika, 85(1):29–37.

Juho Piironen and Aki Vehtari (2016). Comparison of Bayesian predictive methods for model selection. Statistics and Computing, ([online][piironenvehtari]).


  [rstanarm]: https://github.com/stan-dev/rstanarm
  [piironenvehtari]: https://link.springer.com/article/10.1007/s11222-016-9649-y

