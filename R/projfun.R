

project_gaussian <- function(ind, f, sigma2, x, wobs=rep(1.0,dim(as.matrix(x)[1]),
                            wsample=rep(1.0,dim(as.matrix(x)[2]),
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