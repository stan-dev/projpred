# glmproj

An R package to perform projection predictive variable selection for generalized linear models fitted with [rstanarm][].


Installation
------------
    devtools::install_github('paasim/glmproj')


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

    gproj <- glm_proj(fit, d = 20)
    y_pred <- predict(gproj, x_test)
    

[rstanarm]: https://github.com/stan-dev/rstanarm
[glmnet]: https://cran.r-project.org/web/packages/glmnet/index.html

References
------------
Piironen, J. and Vehtari, A. (2015). Comparison of Bayesian predictive methods for model selection. arXiv:1503.08650 ([online][])

  [online]: http://arxiv.org/abs/1503.08650
