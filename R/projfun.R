# Function handles for the projection
#

project_gaussian <- function(ind, p_full, d_train, intercept, regul = 1e-12) {

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

    if(intercept) {
        # add vector of ones to x and transform the variable indices
        x <- cbind(1, x)
        ind <- c(1, ind + 1)
    } else if (length(ind) == 0) {
        # no intercept used and ind is empty, so projecting to the completely
        # null model with eta=0 always
        beta_sub <- matrix(integer(length=0), ncol=NCOL(mu))
        dis_sub <- sqrt( colSums(wobs*mu^2) + dis^2 )
        kl <- weighted.mean(log(dis_sub) - log(dis), wsample)
        p_sub <- list(kl = kl, weights = wsample, dis = dis_sub, ind = ind)
        return(c(p_sub, .split_coef(beta_sub, intercept)))
    }

    xp <- x[, ind, drop = F]
    Dp <- dim(xp)[2]
    regulmat <- diag(regul*rep(1.0, Dp), Dp, Dp)

    # Solve the projection equations (with l2-regularization)
    w <- sqrt(wobs)
    beta_sub <- solve( crossprod(w*xp)+regulmat, crossprod(w*xp, w*mu) )
    dis_sub <- sqrt( colSums(wobs*(mu - xp%*%beta_sub)^2) + dis^2 )
    kl <- weighted.mean(log(dis_sub) - log(dis), wsample)
    p_sub <- list(kl = kl, weights = wsample, dis = dis_sub)

    # split b to alpha and beta, add it to p_sub and return the result
    p_sub <- c(p_sub, .split_coef(beta_sub, intercept))
    if(length(ind) == 1 && intercept) {
      p_sub$ind <- integer(length=0)
    } else {
      p_sub$ind <- ind[(1+intercept*1):length(ind)] - intercept*1
    }
    p_sub$intercept <- intercept
    return(p_sub)
}


project_nongaussian <- function(ind, p_full, d_train, family_kl, intercept,
									regul=1e-9, coef_init=NULL) {
	
	# find the projected regression coefficients for each sample
	xsub <- d_train$x[, ind, drop = F]
	d <- NCOL(xsub)
	S <- NCOL(p_full$mu)
	beta <- matrix(0, nrow=d, ncol=S)
	alpha <- rep(0, S)
	for (s in 1:S) {
		out <- glm_ridge(x = xsub, y = p_full$mu[, s, drop = F],
						 family=family_kl, lambda=regul, weights=d_train$weights,
						 offset=d_train$offset, intercept=intercept, thresh=1e-6) 
		
		beta[,s] <- out$beta
		alpha[s] <- out$beta0
	}
	
	# compute the dispersion parameters and kl-divergences, and combine the results
	p_sub <- list()
	mu <- family_kl$mu_fun(xsub, alpha, beta, d_train$offset)
	p_sub$dis <- family_kl$dis_fun(p_full, d_train, list(mu=mu))
	p_sub$kl <- weighted.mean(family_kl$kl(p_full, d_train, list(mu=mu)), p_full$weights)
	p_sub$weights <- p_full$weights
	p_sub$alpha <- alpha
	p_sub$beta <- beta
	p_sub$ind <- ind
	p_sub$intercept <- intercept
	return(p_sub)
}


# function handle for the projection over samples. Gaussian case
# uses analytical solution to do the projection over samples.
.get_proj_handle <- function(family_kl, regul=1e-9) {

    # Use analytical solution for gaussian as it is a lot faster
    if(family_kl$family == 'gaussian' && family_kl$link == 'identity') {
        #return(project_gaussian)
        return(
            function(ind, p_full, d_train, intercept) {
                project_gaussian(ind, p_full, d_train, intercept, regul=regul)
        })
    } else {
      # return handle to project_nongaussian with family_kl set accordingly
        return(
            function(ind, p_full, d_train, intercept) {
                project_nongaussian(ind, p_full, d_train, family_kl, intercept, regul=regul)
        })
    }
}


.get_submodels <- function(chosen, nv, family_kl, p_full, d_train, intercept, regul) {
    #
    # Project onto given model sizes nv. Returns a list of submodels.
    #
    projfun <- .get_proj_handle(family_kl, regul)

    p_sub <- lapply(nv,
        function(nv) {
            if (nv == 0)
                ind <- integer(length=0) # empty
            else
                ind <- chosen[1:nv]
            return(projfun(ind, p_full, d_train, intercept))
        })
    return(p_sub)
}

