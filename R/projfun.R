#' Function handles for the projection
#'

project_gaussian <- function(ind, p_full, d_train, intercept, regul = 1e-12) {

    x <- d_train$x
    mu <- p_full$mu
    dis <- p_full$dis
    if ("wobs" %in% names(d_train))
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
        p_sub <- list(kl = kl, weights = wsample, dis = dis_sub)
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
    c(p_sub, .split_coef(beta_sub, intercept))
}


project_nongaussian <- function(chosen, p_full, d_train, family_kl, intercept) {

    # perform the projection over samples
    res <- sapply(1:ncol(p_full$mu), function(s_ind) {
        IRLS(list(mu = p_full$mu[, s_ind, drop = F], dis = p_full$dis[s_ind]),
             list(x = d_train$x[, chosen, drop = F], weights = d_train$weights,
                  offset = d_train$offset), family_kl, intercept)
    })

    # weight the results by sample/cluster weights
    list(kl = weighted.mean(unlist(res['kl',]), p_full$weights),
         weights = p_full$weights,
         dis = unlist(res['dis',]),
         alpha = unlist(res['alpha',]),
         beta = do.call(cbind, res['beta',]))
}



project_nongaussian_new <- function(ind, p_full, d_train, family_kl, intercept=TRUE,
									regul=1e-12, coef_init=NULL) {
	
	# find the projected regression coefficients for each sample
	res <- sapply(1:ncol(p_full$mu),
				  function(s) {
				  	glm_ridge(x = d_train$x[, ind, drop = F], y = p_full$mu[, s, drop = F],
				  			  family=family_kl, lambda=regul, weights=d_train$weights,
				  			  offset=d_train$offset, intercept=intercept) 
				  })
	
	
}


# function handle for the projection over samples. Gaussian case
# uses analytical solution to do the projection over samples.
.get_proj_handle <- function(family_kl) {

    # Use analytical solution for gaussian as it is a lot faster
    if(family_kl$family == 'gaussian' && family_kl$link == 'identity') {
        return(project_gaussian)
    } else {
      # return handle to project_nongaussian with family_kl set accordingly
      return(
        function(chosen, p_full, d_train, intercept) {
          project_nongaussian(chosen, p_full, d_train, family_kl, intercept)
        })
    }
}


.get_submodels <- function(chosen, nv, family_kl, p_full, d_train, intercept) {
    #
    # Project onto given model sizes nv. Returns a list of submodels.
    #
    projfun <- .get_proj_handle(family_kl)

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



