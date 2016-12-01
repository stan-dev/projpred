#' Iteratively reweighted least squared
#'
#' Does the projection from the full model to the submodel.
#'
#' Written by hand because glm.fit in R is too slow.
#' - eps and max_it cannot be changed currently.
#' - assumes d_train$x has intercept as the first column.
#'

IRLS <- function(p_full, d_train, family_kl, intercept, regul, coef_init,
                 eps = 1e-12, max_it = 300) {
  b <- coef_init$beta
  if(intercept) {
    b <- c(coef_init$alpha, b)
    d_train$x <- cbind(1, d_train$x)
  }

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

  regulvec <- c((1-intercept)*regul, rep(regul, NCOL(d_train$x) - 1))
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
  p_sub$dis <- family_kl$dis(p_full, d_train, p_sub)
  res <- list(kl = family_kl$kl(p_full, d_train, p_sub) + sum(b^2*regulvec))
  res$dis <- p_sub$dis
  # split b to alpha and beta, add it to p_sub and return the result
  c(res, .split_coef(b, intercept))
}

