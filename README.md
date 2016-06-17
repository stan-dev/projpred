# glmproj

An R package to perform projection predictive variable selection for generalized linear models fitted with [rstanarm][]. 

The package uses forward search as a search heuristic, that is, starting from a submodel model with only intercept, and adds variables one at a time, each time choosing the variable that decreases the KL-divergence to the full model the most. 

Currently, supported models include gaussian with identity link function, binomial with probit and logit link functions and poisson with log link function.

Installation
------------
    devtools::install_github('paasim/glmproj')
    
    
Usage
-----

The package provides two functions, `glm_proj` and `predict`. The first function can be used to perform variable selection and to project the information in the posterior of the full model onto the submodels. The latter can be used to predict with the submodels. For additional information about the functions (eg. function arguments), use `?glm_proj` or `?predict.glmproj` in R.


Example
-------

    library(glmproj)
    library(rstanarm)
    options(mc.cores = parallel::detectCores())
    
    data('crimedata', package = 'glmproj')
    n <- nrow(crimedata)
	n_train <- 1000
	n_test <- n - n_train
    d <- ncol(crimedata)
    y_train <- crimedata[1:n_train, 1]
	x_train <- crimedata[1:n_train, 2:d]
	x_test <- cbind(rep(1, n_test), crimedata[(n_train+1):n, 2:d])

    family <- gaussian()
    prior <- hs(df = 3)
    fit <- stan_glm(y_train ~ x_train, family = family, prior = prior)

    gproj <- glm_proj(fit)
    
    # predict using 20 variables
    y_pred <- predict(gproj, x = x_test, d = 20)


References
------------
Dupuis, J. A. and Robert, C. P. (2003). Variable selection in qualitative models via an entropic explanatory power. Journal of Statistical Planning and Inference, 111(1-2):77–94.

Goutis, C. and Robert, C. P. (1998). Model choice in generalised linear models: a Bayesian approach via Kullback–Leibler projections. Biometrika, 85(1):29–37.

Juho Piironen and Aki Vehtari (2016). Comparison of Bayesian predictive methods for model selection. Statistics and Computing, ([online][piironenvehtari]).



  [rstanarm]: https://github.com/stan-dev/rstanarm
  [piironenvehtari]: https://link.springer.com/article/10.1007/s11222-016-9649-y

