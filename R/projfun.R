



# THIS FUNCTION IS CURRENTLY UNUSED
project_gaussian_new <- function(ind, f, sigma2, x, wobs=rep(1.0,dim(as.matrix(x)[1])),
                            wsample=rep(1.0,dim(as.matrix(x)[2])),
                            type='one-to-one', intercept=TRUE, regul=1e-12) 
{
    
    if(intercept) {
        x <- cbind(1, x)
        ind <- c(1, ind + 1)
    }
    
    xp <- x[, ind, drop = F]
    Dp <- dim(xp)[2]
    regulvec <- regul*rep(1.0, Dp) # c((1-intercept)*regul, rep(regul, length(ind) - 1))
    regulmat <- diag(regulvec, length(regulvec), length(regulvec))
    
    # normalize the weights
    wobs <- wobs/sum(wobs)
    wsample <- wsample/sum(wsample)
    
    # Solution for the gaussian case (with l2-regularization)
    w <- sqrt(wobs)
    b <- solve( crossprod(w*xp)+regulmat, crossprod(w*xp, w*f) )
    dis <- sqrt( colSums(wobs*(f - xp%*%b)^2) + sigma2 ) ## check this!!
    kl <- weighted.mean(log(dis) - log(sqrt(sigma2)), wsample)
    p_sub <- list(kl = kl, dis = dis)
    
    # split b to alpha and beta, add it to p_sub and return the result
    c(p_sub, .split_coef(b, intercept))
}




project_gaussian <- function(chosen, p_full, d_train, intercept, regul, coef_init) {
    
    if(intercept) {
        d_train$x <- cbind(1, d_train$x)
        chosen <- c(1, chosen + 1)
    }
    
    regulvec <- c((1-intercept)*regul, rep(regul, length(chosen) - 1))
    regulmat <- diag(regulvec, length(regulvec), length(regulvec))
    
    w <- sqrt(d_train$weights)
    # Solution for the gaussian case (with l2-regularization)
    b <- solve(crossprod(w*d_train$x[, chosen, drop = F]) + regulmat,
               crossprod(w*d_train$x[, chosen, drop = F], w*p_full$mu))
    dis <- sqrt(colMeans(d_train$weights*(
        p_full$mu - d_train$x[, chosen, drop = F]%*%b)^2) + p_full$dis^2)
    p_sub <- list(kl = weighted.mean(log(dis) - log(p_full$dis) +
                                         colSums(b^2*regulvec), p_full$cluster_w),
                  dis = dis)
    # split b to alpha and beta, add it to p_sub and return the result
    c(p_sub, .split_coef(b, intercept))
}


project_nongaussian <- function(chosen, p_full, d_train, intercept, regul, coef_init) {
    
    # perform the projection over samples
    res <- sapply(1:ncol(p_full$mu), function(s_ind) {
        IRLS(list(mu = p_full$mu[, s_ind, drop = F], dis = p_full$dis[s_ind]),
             list(x = d_train$x[, chosen, drop = F], weights = d_train$weights,
                  offset = d_train$offset), family_kl, intercept, regul,
             within(coef_init, beta <- coef_init$beta[chosen]))
    })
    
    # weight the results by sample/cluster weights
    list(kl = weighted.mean(unlist(res['kl',]), p_full$cluster_w),
         alpha = unlist(res['alpha',]),
         beta = do.call(cbind, res['beta',]),
         dis = unlist(res['dis',]))
}


# function handle for the projection over samples. Gaussian case
# uses analytical solution to do the projection over samples.
.get_proj_handle <- function(family_kl) {
    
    # Use analytical solution for gaussian as it is a lot faster
    if(family_kl$family == 'gaussian' && family_kl$link == 'identity') {
        return(project_gaussian)
    } else {
        return(project_nongaussian)
    }
}


