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
  # as explained in McGullagh and Nelder (1989). Returns also the deviance and its pointwise
  # derivative w.r.t f at the current f (notice though, that this 'deviance' does not contain additional
	# constants, so even when the model fits perfectly to the data, the deviance is not zero).
  #
  mu <- family$linkinv(f+offset)
  dmu_df <- family$mu.eta(f+offset)
  z <- f + (y - mu)/dmu_df
  
  if (family$family == 'Student_t') {
    # Student-t does not belong to the exponential family and thus it has its own
    # way of computing the observation weights
    if (is.null(wprev)) {
      # initialization of the em-iteration; loop recursively until stable initial weights are found
      wprev <- weights
      while(T) {
        wtemp <- pseudo_data(f,y,family, offset=offset, weights=weights, wprev=wprev, obsvar=obsvar)$w
        if (max(abs(wtemp-wprev)) < 1e-6)
          break
        wprev <- wtemp
      }
    }
    # given the weights from the previous em-iteration, update s2 based on the previous weights and mu,
    # and then compute new weights w
    nu <- family$nu
    s2 <- sum(wprev*(obsvar+(y-mu)^2)) / sum(weights) 
    w <- weights*(nu+1)/(nu + 1/s2*(obsvar+(y-mu)^2))
    dev <- sum(family$loss_fun(mu, y, weights, sqrt(s2))) # sum( -2*family$ll_fun(mu, sqrt(s2), y, weights) )
    grad <- weights*2*(mu-y)/(nu*s2) * (nu+1)/(1+(y-mu)^2/(nu*s2)) * dmu_df
    # print(paste0('s2 = ', s2))
    # print(paste0('w = ', w))
    # print(paste0('r^2 = ', (y-mu)^2))
    # print('---------------')
    
    
  } else if (family$family %in% c('gaussian','poisson','binomial')) {
    # exponential family distributions
    w <- (weights * dmu_df^2)/family$variance(mu)
    dev <- sum(family$loss_fun(mu, y, weights)) #sum( -2*family$ll_fun(mu, 1, y, weights) )
    grad <- -2*w*(z-f)
    
  } else {
    stop(sprintf('Don\'t know how to compute quadratic approximation and gradients for family \'%s\'.',
                 family$family))
  }
  
  return(list(z=z, w=w, dev=dev, grad=grad))
}


lambda_grid <- function(x, y, family, offset, weights, intercept, obsvar=0, alpha=1.0, 
												eps=1e-2, nlam=100, ret.all=F) {
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
	
	# find the initial intercept (that is, assuming all coefficients = 0)
	if (intercept) {
	  beta0_old <- 0
	  while(T) {
	    beta0 <- sum(obs$w*obs$z) / sum(obs$w) # intercept
	    if (abs(beta0-beta0_old) < 1e-6)
	      break
	    obs <- pseudo_data(beta0*rep(1,n), y, family, offset, weights, obsvar=obsvar)
	    beta0_old <- beta0
	  } 
	} else
	  beta0 <- 0
	
	resid <- obs$z - beta0 # residual after taking into account the intercept
	lambda_max <- max(abs( t(x) %*% (resid*obs$w) )) / alpha
	lambda_min <- eps*lambda_max
	loglambda <- seq(log(lambda_min), log(lambda_max), len=nlam)
	if (ret.all)
		return( list(lambda = rev(exp(loglambda)), beta0=beta0, w0=obs$w) )
	else
		return(rev(exp(loglambda)))
}



glm_elnet <- function(x, y, family=gaussian(), nlambda=100, lambda_min_ratio=1e-3,
                      lambda=NULL, alpha=1.0, thresh=1e-6,
                      qa_updates_max=ifelse(family$family=='gaussian' &&
                                              family$link=='identity', 1, 100),
                      pmax=dim(as.matrix(x))[2]+1, pmax_strict=FALSE,
                      weights=NULL, offset=NULL, obsvar=0, intercept=TRUE, normalize=TRUE) {
	#
	# Fits GLM with elastic net penalty on the regression coefficients.
	# Computes the whole regularization path.
	# Does not handle any dispersion parameters.
	#
	# np <- dim(x)
	# if (is.null(np) || (np[2] <= 1))
		# stop("x should be a matrix with 2 or more columns")

	# ensure x is in matrix form and fill in missing weights and offsets
	x <- as.matrix(x)
	if (is.null(weights))
		weights <- rep(1.0, nrow(x))
	if (is.null(offset))
		offset <- rep(0.0, nrow(x))

	if (normalize) {
		# normalize the predictor matrix. notice that the variables are centered only if
	  # intercept is used.
	  if (intercept)
		  mx <- colMeans(x)
	  else
	    mx <- rep(0,ncol(x))
		sx <- apply(x,2,'sd')
		x <- scale(x, center=intercept, scale=T)
	}

	# default lambda-sequence, including optimal start point
	if (is.null(lambda)) {
		temp <- lambda_grid(x, y, family, offset, weights, intercept, alpha, 
		                    obsvar=obsvar, nlam=nlambda, eps=lambda_min_ratio, ret.all = T)
		lambda <- temp$lambda
		w0 <- temp$w0
		beta0 <- temp$beta0
	} else {
		beta0 <- 0
		w0 <- weights
	}
		
	# call the c++-function that serves as the workhorse
	pseudo_obs <- function(f,wprev) pseudo_data(f,y,family,offset=offset,weights=weights,obsvar=obsvar,wprev=wprev)
	out <- glm_elnet_c(x,pseudo_obs,lambda,alpha,intercept,thresh,qa_updates_max,pmax,pmax_strict,beta0,w0)
	beta <- out[[1]]
	beta0 <- as.vector(out[[2]])

	if (normalize) {
		# return the intecept and the coefficients on the original scale
		# beta <- sweep(beta, 1, sx, '/')
	  beta <- beta/sx
		# beta0 <- beta0 - colSums(sweep(beta, 1, mx, '*'))
	  beta0 <- beta0 - colSums(mx*beta)
	}

	return(list( beta=beta, beta0=beta0, npasses=out[[3]],
				 updates_qa=as.vector(out[[4]]), updates_as=as.vector(out[[5]]) ))
}


glm_ridge <- function(x, y, family=gaussian(), lambda=0, thresh=1e-9, qa_updates_max=NULL,
                      weights=NULL, offset=NULL, obsvar=0, intercept=TRUE, ls_iter_max=30) {
  #
  # Fits GLM with ridge penalty on the regression coefficients.
  # Does not handle any dispersion parameters.
  #
  if (family$family == 'gaussian' && family$link == 'identity') {
    qa_updates_max <- 1
    ls_iter_max <- 1
  }
  else if (is.null(qa_updates_max))
    qa_updates_max <- 100
  
	if (length(x) == 0 && !intercept)
		# null model with no predictors and no intercept
		out <- list( beta=matrix(integer(length=0)), beta0=0, w=rep(1,length(y)), qa_updates=0 )
	else {
		# normal case
		x <- as.matrix(x)
		if (is.null(weights))
			weights <- rep(1.0, nrow(x))
		if (is.null(offset))
			offset <- rep(0.0, nrow(x))
		
		w0 <- weights 
		pseudo_obs <- function(f,wprev) pseudo_data(f,y,family,offset=offset,weights=weights,obsvar=obsvar,wprev=wprev)
		out <- glm_ridge_c(x, pseudo_obs, lambda, intercept, thresh, qa_updates_max, w0,ls_iter_max)
	}
  return(list( beta=out[[1]], beta0=out[[2]], w=out[[3]], qa_updates=out[[4]] ))
}


glm_forward <- function(x, y, family=gaussian(), lambda=0, thresh=1e-9, qa_updates_max=NULL,
                        weights=NULL, offset=NULL, obsvar=0, intercept=TRUE,
                        pmax=dim(as.matrix(x))[2]) {
  #
  # Runs forward stepwise regression. Does not handle any dispersion parameters.
  #
  if (family$family == 'gaussian' && family$link == 'identity')
    qa_updates_max <- 1
  else if (is.null(qa_updates_max))
    qa_updates_max <- 100
  
  if (length(x) == 0 && !intercept)
    # null model with no predictors and no intercept
    return( list( beta=matrix(integer(length=0)), beta0=0, varorder=integer(length=0) ) )
  else {
    # normal case
    x <- as.matrix(x)
    if (is.null(weights))
      weights <- rep(1.0, nrow(x))
    if (is.null(offset))
      offset <- rep(0.0, nrow(x))
    
    w0 <- weights
    pseudo_obs <- function(f,wprev) pseudo_data(f,y,family,offset=offset,weights=weights,obsvar=obsvar,wprev=wprev)
    out <- glm_forward_c(x, pseudo_obs, lambda, intercept, thresh, qa_updates_max, pmax, w0)
  }
  return(list( beta=out[[1]], beta0=as.vector(out[[2]]), varorder=as.vector(out[[3]])+1 ))
}





