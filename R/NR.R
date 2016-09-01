#' Newton-Rhapson
#'
#' Minimizes the KL divergence between the full model and the submodel
#' using Newton-Rhapson algorithm.
#'
#' Naming convention described in fsel.R.


NR <- function(p_full, d_train, b0, family_kl, max_it = 30, eps = 1e-10) {

  p_sub <- list(b = b0)

  it <- 1

  p_sub$eta <- d_train$x%*%p_sub$b
  p_sub$mu <- family_kl$linkinv(p_sub$eta + d_train$offset)
  p_sub$dis <- family_kl$dis(p_full, d_train, p_sub)
  kl <- family_kl$kl(p_full, d_train, p_sub)

  alpha <- 0.1
  beta <- 0.5

  while(it < max_it) {

    dkl <- family_kl$dkl(p_full, d_train, p_sub, family_kl$dme(p_sub))
    db <- solve(dkl$h, dkl$g)
    decr <- crossprod(dkl$g, db)

    # check convergence
    if(decr*0.5 < eps) break

    step <- 1

    p_sub$eta <- d_train$x%*%(p_sub$b - step*db)
    p_sub$mu <- family_kl$linkinv(p_sub$eta + d_train$offset)
    p_sub$dis <- family_kl$dis(p_full, d_train, p_sub)
    kl_new <- family_kl$kl(p_full, d_train, p_sub)

    # find stepsize
    while(kl_new > kl - (alpha*step)*decr) {
      step <- step*beta
      p_sub$eta <- d_train$x%*%(p_sub$b - step*db)
      p_sub$mu <- family_kl$linkinv(p_sub$eta + d_train$offset)
      p_sub$dis <- family_kl$dis(p_full, d_train, p_sub)
      kl_new <- family_kl$kl(p_full, d_train, p_sub)
    }

    p_sub$b <- p_sub$b - step*db
    kl <- kl_new

    # check convergence
    if(step < 1e-8) break

    it <- it + 1
  }

  if(it >= max_it)
    warning(paste0("Maximum number of iterations reached."))

  res <- list(kl = kl, b = p_sub$b)
  if(family_kl$family %in% c('gaussian','Gamma')) res$dis <- p_sub$dis
  res
}
