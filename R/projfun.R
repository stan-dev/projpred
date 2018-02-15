# Function handles for the projection
#

project_gaussian <- function(vind, p_full, d_train, family_kl, intercept, regul = 1e-12) {

    x <- d_train$x
    mu <- p_full$mu
    dis <- p_full$dis
    
    if ("weights" %in% names(d_train))
        wobs <- d_train$weights
    else
        wobs <- rep(1.0, NROW(mu))
    if ("weights" %in% names(p_full))
        wsample <- p_full$weights
    else
        wsample <- rep(1.0, NCOL(mu))

    # ensure the weights are normalized
    wobs <- wobs/sum(wobs)
    wsample <- wsample/sum(wsample)

    if (intercept) {
        # add vector of ones to x and transform the variable indices
        x <- cbind(1, x)
        vind <- c(1, vind + 1)
    } else if (length(vind) == 0) {
        # no intercept used and vind is empty, so projecting to the completely
        # null model with eta=0 always
    		pobs <- pseudo_data(0, mu, family_kl, offset=d_train$offset, weights=wobs)
        beta_sub <- matrix(integer(length=0), ncol=NCOL(mu))
        dis_sub <- family_kl$dis_fun(list(mu=pobs$z, var=p_full$var), list(mu=0), pobs$w)
        kl <- weighted.mean(colSums(wobs*pobs$z^2), wsample)
        submodel <- list(kl = kl, weights = wsample, dis = dis_sub, vind = vind,
                      intercept = intercept)
        return(c(submodel, .split_coef(beta_sub, intercept)))
    }

    xp <- x[, vind, drop = F]
    Dp <- dim(xp)[2]
    regulmat <- diag(regul*rep(1.0, Dp), Dp, Dp)

    # Solve the projection equations (with l2-regularization)
    pobs <- pseudo_data(0, mu, family_kl, offset=d_train$offset, weights=wobs) # this will remove the offset
    wsqrt <- sqrt(pobs$w)
    beta_sub <- solve( crossprod(wsqrt*xp)+regulmat, crossprod(wsqrt*xp, wsqrt*pobs$z) )
    musub <- xp%*%beta_sub
    dis_sub <- family_kl$dis_fun(list(mu=pobs$z, var=p_full$var), list(mu=musub), pobs$w)
    kl <- weighted.mean(colSums(wobs*((pobs$z-musub)^2)), wsample) # not the actual kl-divergence, but a reasonable surrogate..
    submodel <- list(kl = kl, weights = wsample, dis = dis_sub)

    # split b to alpha and beta, add it to submodel and return the result
    submodel <- c(submodel, .split_coef(beta_sub, intercept))
    if(length(vind) == 1 && intercept) {
      submodel$vind <- integer(length=0)
    } else {
      submodel$vind <- vind[(1+intercept*1):length(vind)] - intercept*1
    }
    submodel$intercept <- intercept
    return(submodel)
}




project_nongaussian <- function(vind, p_full, d_train, family_kl, intercept,
									regul=1e-9, coef_init=NULL) {
	
	# find the projected regression coefficients for each sample
	xsub <- d_train$x[, vind, drop = F]
	d <- NCOL(xsub)
	n <- NROW(p_full$mu)
	S <- NCOL(p_full$mu)
	
  # loop through each draw and compute the projection for it individually
  beta <- matrix(0, nrow=d, ncol=S)
  alpha <- rep(0, S)
  w <- matrix(nrow=n, ncol=S)
  # beta_init <- beta[,1]
  # beta0_init <- alpha[1]
  for (s in 1:S) {
    out <- glm_ridge(x = xsub, y = p_full$mu[, s, drop = F],
                     family=family_kl, lambda=regul, weights=d_train$weights,
                     offset=d_train$offset, obsvar=p_full$var[,s], intercept=intercept) 
    beta[,s] <- out$beta
    alpha[s] <- out$beta0
    w[,s] <- out$w
    # beta_init <- beta[,s]
    # beta0_init <- alpha[s]
  }
	
	# compute the dispersion parameters and kl-divergences, and combine the results
	submodel <- list()
	mu <- family_kl$mu_fun(xsub, alpha, beta, d_train$offset)
	submodel$dis <- family_kl$dis_fun(p_full, list(mu=mu,w=w), d_train$weights)
	submodel$kl <- weighted.mean(family_kl$kl(p_full, d_train, list(mu=mu,dis=submodel$dis)), p_full$weights)
	submodel$weights <- p_full$weights
	submodel$alpha <- alpha
	submodel$beta <- beta
	submodel$vind <- vind
	submodel$intercept <- intercept
	return(submodel)
}



# function handle for the projection over samples. Gaussian case
# uses analytical solution to do the projection over samples.
.get_proj_handle <- function(family_kl, regul=1e-9) {

    # Use analytical solution for gaussian because it is faster
    if(family_kl$family == 'gaussian' && family_kl$link == 'identity') {
        return(
            function(vind, p_full, d_train, intercept) {
                project_gaussian(vind, p_full, d_train, family_kl, intercept, regul=regul)
        })
    } else {
      # return handle to project_nongaussian with family_kl set accordingly
        return(
            function(vind, p_full, d_train, intercept) {
                project_nongaussian(vind, p_full, d_train, family_kl, intercept, regul=regul)
        })
    }
}


.get_submodels <- function(vind, nv, family_kl, p_full, d_train, intercept, regul) {
    #
    # Project onto given model sizes nv. Returns a list of submodels.
    #
    projfun <- .get_proj_handle(family_kl, regul)

    submodels <- lapply(nv,
        function(nv) {
            if (nv == 0)
                vind <- integer(length=0) # empty
            else
                vind <- vind[1:nv]
            return(projfun(vind, p_full, d_train, intercept))
        })
    return(submodels)
}

