#
# The functions in this file are used to compute the elastic net coefficient paths
# for a GLM. The main function is glm_elnet, other functions are auxiliaries.
# The L1-regularized projection path is computed by replacing the actual data y
# by the fit of the full model when calling glm_elnet. glm_elnet uses function glm_elnet_c
# from elnetfun.cpp.
#

pseudo_data <- function(f, y, family, offset=rep(0,length(f)), weights=rep(1.0,length(f)) ) {
    #
    # Returns locations z and weights w (inverse-variances) of the Gaussian pseudo-observations
    # based on the quadratic approximation to the loss function (negative log likelihood) at 
    # when the given fit f = eta = x*beta + beta0. Returns also the deviance at f.
    #
    mu <- family$linkinv(f)
    dmu_df <- family$mu.eta(f)
    z <- (f - offset) + (y - mu)/dmu_df
    w <- (weights * dmu_df^2)/family$variance(mu)
    dev <- sum( family$dev.resids(y, mu, weights) )
    return(list(z=z, w=w, dev=dev))
}


lambda_grid <- function(x, y, family, alpha=1.0, eps=1e-2, nlam=100) {
	#
    # Standard lambda sequence as described in Friedman et al. (2009), section 2.5.
    # The grid will have nlam values, evenly spaced in the log-space between lambda_max
    # and lambda_min. lambda_max is the smallest value for which all the regression
    # coefficients will be zero.
    #
	n <- dim(x)[1]
	obs <- pseudo_data(rep(0,n), y, family)
	
	if (alpha == 0)
	    # initialize ridge as if alpha = 0.01
	    alpha <- 0.01
	
	lambda_max <- max(abs( t(x) %*% (obs$z*obs$w) )) / alpha
	lambda_min <- eps*lambda_max
	loglambda <- seq(log(lambda_min), log(lambda_max), len=nlam)
	return(rev(exp(loglambda)))
}


glm_elnet <- function(x, y, family=gaussian(), nlambda=100, lambda_min_ratio=1e-3,
                      lambda=NULL, alpha=1.0, thresh=1e-6, 
					  qa_updates_max=ifelse(family$family=='gaussian', 1, 100), 
					  pmax=dim(as.matrix(x))[2], pmax_strict=FALSE,
					  weights=NULL, offset=NULL, intercept=TRUE) {
	#
	# Fits GLM with elastic net penalty on the regression coefficients.
	# Computes the whole regularization path.
	# Does not handle any dispersion parameters.
	#
	np <- dim(x)
	if (is.null(np) || (np[2] <= 1)) 
		stop("x should be a matrix with 2 or more columns")
	
	if (is.null(lambda))
		lambda <- lambda_grid(x,y,family,alpha,nlam=nlambda,eps=lambda_min_ratio)
	
	x <- as.matrix(x)
	if (is.null(weights))
		weights <- 1.0
	if (is.null(offset))
		offset <- 0.0
	pseudo_obs <- function(f) {return(pseudo_data(f,y,family,offset=offset,weights=weights))}
	out <- glm_elnet_c(x,pseudo_obs,lambda,alpha,intercept,thresh,qa_updates_max,pmax,pmax_strict)
	return(list( beta=out[[1]], beta0=as.vector(out[[2]]), npasses=out[[3]], 
				 updates_qa=as.vector(out[[4]]), updates_as=as.vector(out[[5]]) ))
}


glm_ridge <- function(x, y, family=gaussian(), lambda=0, thresh=1e-6, 
                      qa_updates_max=ifelse(family$family=='gaussian', 1, 100),
                      weights=NULL, offset=NULL, intercept=TRUE) {
    #
    # Fits GLM with ridge penalty on the regression coefficients.
    # Does not handle any dispersion parameters.
    #
	if (length(x) == 0 && !intercept)
		# null model with no predictors and no intercept
		out <- list( beta=matrix(integer(length=0)), beta0=0, qa_updates=0 )
	else {
		# normal case
		x <- as.matrix(x)
		if (is.null(weights))
			weights <- 1.0
		if (is.null(offset))
			offset <- 0.0
		pseudo_obs <- function(f) {return(pseudo_data(f,y,family,offset=offset,weights=weights))}
		out <- glm_ridge_c(x, pseudo_obs, lambda, intercept, thresh, qa_updates_max)
	}
    return(list( beta=out[[1]], beta0=out[[2]], qa_updates=out[[3]] ))
}








