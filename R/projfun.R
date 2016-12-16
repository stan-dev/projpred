



project_gaussian <- function(ind, p_full, d_train, intercept=TRUE, regul=1e-12, coef_init=NULL) {

    x <- d_train$x
    mu <- p_full$mu
    dis <- p_full$dis
    if ("wobs" %in% names(d_train))
        wobs <- d_train$weights
    else
        wobs <- rep(1.0, dim(as.matrix(mu))[1])
    if ("weights" %in% names(p_full))
        wsample <- p_full$weights
    else
        wsample <- rep(1.0, dim(as.matrix(mu))[2])

    if(intercept) {
        # add vector of ones to x and transform the variable indices
        x <- cbind(1, x)
        ind <- c(1, ind + 1)
    }

    xp <- x[, ind, drop = F]
    Dp <- dim(xp)[2]
    regulmat <- diag(regul*rep(1.0, Dp), Dp, Dp)

    # normalize the weights
    wobs <- wobs/sum(wobs)
    wsample <- wsample/sum(wsample)

    # Solve the projection equations (with l2-regularization)
    w <- sqrt(wobs)
    beta_sub <- solve( crossprod(w*xp)+regulmat, crossprod(w*xp, w*mu) )
    dis_sub <- sqrt( colSums(wobs*(mu - xp%*%beta_sub)^2) + dis^2 )
    kl <- weighted.mean(log(dis_sub) - log(dis), wsample)
    p_sub <- list(kl = kl, weights = wsample, dis = dis_sub)

    # split b to alpha and beta, add it to p_sub and return the result
    c(p_sub, .split_coef(beta_sub, intercept))
}




project_nongaussian <- function(chosen, p_full, d_train, family_kl, intercept=TRUE,
                                regul=1e-12, coef_init=NULL) {

    if (is.null(coef_init)) {
        # initialize at the origin
        coef_init <- list( alpha=0, beta=matrix(rep(0,length(chosen)), ncol=1) )
    } else {
        # simply pick up the relevant indices from beta
        coef_init <- within(coef_init, beta <- coef_init$beta[chosen])
    }
    
    # perform the projection over samples
    res <- sapply(1:ncol(p_full$mu), function(s_ind) {
        IRLS(list(mu = p_full$mu[, s_ind, drop = F], dis = p_full$dis[s_ind]),
             list(x = d_train$x[, chosen, drop = F], weights = d_train$weights,
                  offset = d_train$offset), family_kl, intercept, regul, coef_init)
    })

    # weight the results by sample/cluster weights
    list(kl = weighted.mean(unlist(res['kl',]), p_full$weights),
         weights = p_full$weights,
         dis = unlist(res['dis',]),
         alpha = unlist(res['alpha',]),
         beta = do.call(cbind, res['beta',]))
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
        function(chosen, p_full, d_train, intercept, regul=1e-12, coef_init=NULL) {
          project_nongaussian(chosen, p_full, d_train, family_kl, intercept, regul, coef_init)
        })
    }
}

















