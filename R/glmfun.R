#
# The functions in this file are used to compute the elastic net coefficient paths
# for a GLM. The main function is glm_elnet, other functions are auxiliaries.
# The L1-regularized projection path is computed by replacing the actual data y
# by the fit of the full model when calling glm_elnet. Uses functions in glmfun.cpp.
#

pseudo_data <- function(f, y, family, offset=rep(0,length(f)), weights=rep(1.0,length(f)), obsvar=0, wprev=NULL) {
  #
  # Returns locations z and weights w (inverse-variances) of the Gaussian pseudo-observations
  # based on the linear approximation to the link function at f = eta = x*beta + beta0,
  # as explained in McGullagh and Nelder (1989). Returns also the deviance at f.
  #
  f <- f + offset
  mu <- family$linkinv(f)
  dmu_df <- family$mu.eta(f)
  z <- (f - offset) + (y - mu)/dmu_df
  if (family$family == 'Student_t') {
  	# Student-t does not belong to the exponential family and thus it has its own
  	# way of computing the observation weights
  	if (is.null(wprev)) {
  		# initialization of the em-iteration; loop recursively until stable initial weights are found
  		wprev <- weights
  		while(T) {
  			wtemp <- projpred:::pseudo_data(f,y,family, offset=offset, wprev=wprev)$w
  			if (max(abs(wtemp-wprev)) < 1e-6)
  				break
  			wprev <- wtemp
  		}
  	}
  	# given the weights from the previous em-iteration, update s2 based on the previous weights and mu,
  	# and then compute new weights w
  	nu <- family$nu
  	s2 <- sum(wprev/sum(weights)*(obsvar+(z-mu)^2))
  	w <- (nu+1)/(nu + 1/s2*(obsvar+(z-mu)^2))
  } else
  	w <- (weights * dmu_df^2)/family$variance(mu)
  dev <- sum( family$dev.resids(y, mu, weights) )
  return(list(z=z, w=w, dev=dev))
}


lambda_grid <- function(x, y, family, offset, weights, obsvar=0, alpha=1.0, 
												eps=1e-2, nlam=100, ret.init.weights=F) {
	#
  # Standard lambda sequence as described in Friedman et al. (2009), section 2.5.
  # The grid will have nlam values, evenly spaced in the log-space between lambda_max
  # and lambda_min. lambda_max is the smallest value for which all the regression
  # coefficients will be zero (assuming alpha > 0, alpha = 0 will be initialized 
	# as if alpha = 0.01).
  #
	n <- dim(x)[1]
	obs <- pseudo_data(rep(0,n), y, family, offset, weights, obsvar=obsvar)

	if (alpha == 0)
	    # initialize ridge as if alpha = 0.01
	    alpha <- 0.01

	lambda_max <- max(abs( t(x) %*% (obs$z*obs$w) )) / alpha
	lambda_min <- eps*lambda_max
	loglambda <- seq(log(lambda_min), log(lambda_max), len=nlam)
	if (ret.init.weights)
		return( list(lambda = rev(exp(loglambda)), w0=obs$w) )
	else
		return(rev(exp(loglambda)))
}


glm_elnet <- function(x, y, family=gaussian(), nlambda=100, lambda_min_ratio=1e-3,
                      lambda=NULL, alpha=1.0, thresh=1e-6,
                      qa_updates_max=ifelse(family$family=='gaussian' &&
                                              family$link=='identity', 1, 100),
                      pmax=dim(as.matrix(x))[2], pmax_strict=FALSE,
                      weights=NULL, offset=NULL, obsvar=0, intercept=TRUE, normalize=TRUE) {
	#
	# Fits GLM with elastic net penalty on the regression coefficients.
	# Computes the whole regularization path.
	# Does not handle any dispersion parameters.
	#
	np <- dim(x)
	if (is.null(np) || (np[2] <= 1))
		stop("x should be a matrix with 2 or more columns")

	# ensure x is in matrix form and fill in missing weights and offsets
	x <- as.matrix(x)
	if (is.null(weights))
		weights <- rep(1.0, nrow(x))
	if (is.null(offset))
		offset <- rep(0.0, nrow(x))

	if (normalize) {
		# normalize the predictor matrix
		mx <- colMeans(x)
		sx <- apply(x,2,'sd')
		x <- scale(x, center=T, scale=T)
	}

	# default lambda-sequence
	if (is.null(lambda)) {
		temp <- lambda_grid(x, y, family, offset, weights, alpha, obsvar=obsvar, nlam=nlambda,
												eps=lambda_min_ratio, ret.init.weights = T)
		lambda <- temp$lambda
		w0 <- temp$w0
	} else
		w0 <- weights
		

	# call the c++-function that serves as the workhorse
	pseudo_obs <- function(f,wprev) {return(pseudo_data(f,y,family,offset=offset,weights=weights,obsvar=obsvar,wprev=wprev))}
	out <- glm_elnet_c(x,pseudo_obs,lambda,alpha,intercept,thresh,qa_updates_max,pmax,pmax_strict,w0)
	beta <- out[[1]]
	beta0 <- as.vector(out[[2]])

	if (normalize) {
		# return the intecept and the coefficients on the original scale
		beta <- sweep(beta, 1, sx, '/')
		beta0 <- beta0 - colSums(sweep(beta, 1, mx, '*'))
	}

	return(list( beta=beta, beta0=beta0, npasses=out[[3]],
				 updates_qa=as.vector(out[[4]]), updates_as=as.vector(out[[5]]) ))
}


glm_ridge <- function(x, y, family=gaussian(), lambda=0, thresh=1e-6,
                      qa_updates_max=ifelse(family$family=='gaussian' &&
                                              family$link=='identity', 1, 100),
                      weights=NULL, offset=NULL, obsvar=0, intercept=TRUE) {
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
			weights <- rep(1.0, nrow(x))
		if (is.null(offset))
			offset <- rep(0.0, nrow(x))
		pseudo_obs <- function(f,wprev) {return(pseudo_data(f,y,family,offset=offset,weights=weights,obsvar=obsvar,wprev=wprev))}
		out <- glm_ridge_c(x, pseudo_obs, lambda, intercept, thresh, qa_updates_max, weights)
	}
  return(list( beta=out[[1]], beta0=out[[2]], qa_updates=out[[3]] ))
}


glm_forward <- function(x, y, family=gaussian(), lambda=0, thresh=1e-6,
                        qa_updates_max=ifelse(family$family=='gaussian' &&
                                                family$link=='identity', 1, 100),
                        weights=NULL, offset=NULL, obsvar=0, intercept=TRUE,
                        pmax=dim(as.matrix(x))[2]) {
  #
  # Runs forward stepwise regression. Does not handle any dispersion parameters.
  #
  if (length(x) == 0 && !intercept)
    # null model with no predictors and no intercept
    return( list( beta=matrix(integer(length=0)), beta0=0, varorder=integer(length=0) ) )
  else {
    # normal case
    x <- as.matrix(x)
    if (is.null(weights))
      weights <- rep(1, nrow(x))
    if (is.null(offset))
      offset <- rep(0.0, nrow(x))
    pseudo_obs <- function(f,wprev) pseudo_data(f,y,family,offset=offset,weights=weights,obsvar=obsvar,wprev=wprev)
    out <- glm_forward_c(x, pseudo_obs, lambda, intercept, thresh, qa_updates_max, pmax, weights)
  }
  return(list( beta=out[[1]], beta0=as.vector(out[[2]]), varorder=as.vector(out[[3]])+1 ))
}





