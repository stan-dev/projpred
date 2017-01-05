#' Iteratively reweighted least squared
#'
#' Does the projection from the full model to the submodel.
#'
#' Written by hand because glm.fit in R is too slow.
#' - Initialization of b currenly quite ad hoc
#' - rewrite in C?
#'

IRLS <- function(p_full, d_train, family_kl, intercept, regul = 1e-12,
                 eps = 1e-12, max_it = 300) {

  # check if no intercept and no explanatory variables
  if(NCOL(d_train$x) == 0 && !intercept) {
    mu <- family_kl$linkinv(rep(0, NROW(d_train$x)))
    dis <- family_kl$dis_fun(p_full, d_train, list(mu = mu))
    kl <- family_kl$kl(p_full, d_train, list(mu = mu, dis = dis))
    return(c(list(kl = kl, dis = dis),  .split_coef(matrix(0,0,1), intercept)))
  }

  b <- matrix(0, NCOL(d_train$x), 1)
  if(intercept) {
    b <- c(0, b)
    d_train$x <- cbind(1, d_train$x)
  }
  # if initializing b to 0 doesn't produce valid eta, set it to 1
  if(!family_kl$valideta(drop(d_train$x%*%b))) b <- b + 1

  eta <- drop(d_train$x%*%b)
  mu <- family_kl$linkinv(eta + d_train$offset)
  dev <- sum(family_kl$dev.resids(p_full$mu, mu, d_train$weights))
  stepsize <- 1
  if((!is.finite(dev) | !family_kl$valideta(eta) | !family_kl$validmu(mu)) &
     stepsize > eps) {
    b <- b / 2
    stepsize <- stepsize / 2
    eta <- drop(d_train$x%*%b)
    mu <- family_kl$linkinv(eta + d_train$offset)
    dev <- sum(family_kl$dev.resids(p_full$mu, mu, d_train$weights))
  }
  if(stepsize <= eps) stop('Can\'t initialize the projection.')

  regulvec <- c((1-intercept)*regul, rep(regul, max(NCOL(d_train$x) - 1, 0)))
  regulmat <- diag(regulvec, length(regulvec), length(regulvec))

  it <- 1
  while(it < max_it) {
    mu.eta <- family_kl$mu.eta(eta)
    w <- sqrt(d_train$weights*mu.eta^2/family_kl$variance(mu))
    z <- eta - d_train$offset + (p_full$mu-mu)/mu.eta
    b_old <- b
    b <- solve(crossprod(d_train$x*w) + regulmat, crossprod(d_train$x*w, w*z))
    eta <- drop(d_train$x%*%b)
    mu <- family_kl$linkinv(eta + d_train$offset)
    devold <- dev
    dev <- sum(family_kl$dev.resids(p_full$mu, mu, d_train$weights))

    # line search, deviance must be decreasing.
    stepsize <- 1
    while((dev > devold | !family_kl$valideta(eta) | !family_kl$validmu(mu)) &
          stepsize > eps) {
      b <- (b + b_old)/2
      stepsize <- stepsize/2
      eta <- drop(d_train$x%*%b)
      mu <- family_kl$linkinv(eta + d_train$offset)
      dev <- sum(family_kl$dev.resids(p_full$mu, mu, d_train$weights))
    }
    if(stepsize < eps) {
      b <- b_old
      eta <- drop(d_train$x%*%b)
      mu <- family_kl$linkinv(eta + d_train$offset)
      dev <- sum(family_kl$dev.resids(p_full$mu, mu, d_train$weights))
      break
    }

    if(abs(dev - devold)/(0.1 + abs(dev)) < eps) break
    it <- it + 1
  }
  if(it >= max_it) warning("Maximum number of iterations reached.")
  if(!family_kl$valideta(eta) | !family_kl$validmu(mu))
    warning("Numerical problems in the projection.")

  p_sub <- list(mu = mu)
  p_sub$dis <- family_kl$dis_fun(p_full, d_train, p_sub)
  res <- list(kl = family_kl$kl(p_full, d_train, p_sub) + sum(b^2*regulvec))
  res$dis <- p_sub$dis
  # split b to alpha and beta, add it to p_sub and return the result
  c(res, .split_coef(b, intercept))
}

