#' Iteratively reweighted least squared
#'
#' Does the projection from the full model to the submodel.
#'
#' Written by hand because glm.fit in R is too slow.
#'

IRLS <- function(p_full, d_train, b0, args) {
  fam <- args$family_kl
  b <- b0
  eta <- drop(d_train$x%*%b)
  mu <- fam$linkinv(eta + d_train$offset)
  dev <- sum(fam$dev.resids(p_full$mu, mu, d_train$weights))
  stepsize <- 1
  if((!is.finite(dev) | !fam$valideta(eta) | !fam$validmu(mu)) &
     stepsize > args$epsilon) {
    b <- b / 2
    stepsize <- stepsize / 2
    eta <- drop(d_train$x%*%b)
    mu <- fam$linkinv(eta + d_train$offset)
    dev <- sum(fam$dev.resids(p_full$mu, mu, d_train$weights))
  }
  if(stepsize <= args$epsilon) stop('Can\'t initialize the projection.')

  regulvec <- c((1-args$intercept)*args$regul, rep(args$regul, NCOL(d_train$x) - 1))
  regulmat <- diag(regulvec, length(regulvec), length(regulvec))

  it <- 1
  while(it < args$max_it) {
    mu.eta <- fam$mu.eta(eta)
    w <- sqrt(d_train$weights*mu.eta^2/fam$variance(mu))
    z <- eta - d_train$offset + (p_full$mu-mu)/mu.eta
    b_old <- b
    b <- solve(crossprod(d_train$x*w) + regulmat, crossprod(d_train$x*w, w*z))
    eta <- drop(d_train$x%*%b)
    mu <- fam$linkinv(eta + d_train$offset)
    devold <- dev
    dev <- sum(fam$dev.resids(p_full$mu, mu, d_train$weights))

    # line search, deviance must be decreasing.
    stepsize <- 1
    while((dev > devold | !fam$valideta(eta) | !fam$validmu(mu)) &
          stepsize > args$epsilon) {
      b <- (b + b_old)/2
      stepsize <- stepsize/2
      eta <- drop(d_train$x%*%b)
      mu <- fam$linkinv(eta + d_train$offset)
      dev <- sum(fam$dev.resids(p_full$mu, mu, d_train$weights))
    }
    if(stepsize < args$epsilon) {
      b <- b_old
      eta <- drop(d_train$x%*%b)
      mu <- fam$linkinv(eta + d_train$offset)
      dev <- sum(fam$dev.resids(p_full$mu, mu, d_train$weights))
      break
    }

    if(abs(dev - devold)/(0.1 + abs(dev)) < args$epsilon) break
    it <- it + 1
  }
  if(it >= args$max_it) warning("Maximum number of iterations reached.")
  if(!fam$valideta(eta) | !fam$validmu(mu)) warning("Numerical problems in the projection.")

  p_sub <- list(mu = mu)
  p_sub$dis <- fam$dis(p_full, d_train, p_sub)
  res <- list(kl = fam$kl(p_full, d_train, p_sub) + sum(b^2*regulvec), b = b)
  if(fam$family %in% c('gaussian','Gamma')) res$dis <- p_sub$dis
  res
}

