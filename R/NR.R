#' Newton-Rhapson
#'
#' Minimizes the KL divergence between the full model and the submodel
#' using Newton-Rhapson algorithm.
#'
#' @param \code{mu_p} Fitted values of the full model.
#' @param \code{x} A model matrix of the selected variables.
#' @param \code{b_p} Sampled estimates of the coefficient of the full model.
#' @param \code{w} Observation weights.
#' @param \code{dis_p} dispersion parameter of the full model.
#' @param \code{funs} Model-specific helper functions.
#' @param \code{max_it} Maximum number of iterations for the algorithm. Defaults to 30.
#' @param \code{eps} Tolerance, when derivative to any direction is at most eps, algorithm stops. Defaults to 1e-10.

NR <- function(p, d, b0, family_kl, max_it = 30, eps = 1e-10) {

  q <- list(b = b0)

  it <- 1

  q$eta <- d$x%*%q$b
  q$mu <- family_kl$linkinv(q$eta + d$offset)
  q$dis <- family_kl$dis(p, d, q)
  kl <- family_kl$kl(p, d, q)

  alpha <- 0.1
  beta <- 0.5

  while(it < max_it) {

    dkl <- family_kl$dkl(p, d, q, family_kl$dme(q))
    db <- solve(dkl$s, dkl$f)
    decr <- crossprod(dkl$f, db)

    # check convergence
    if(decr*0.5 < eps) break

    step <- 1

    q$eta <- d$x%*%(q$b - step*db)
    q$mu <- family_kl$linkinv(q$eta + d$offset)
    q$dis <- family_kl$dis(p, d, q)
    kl_new <- family_kl$kl(p, d, q)

    # find stepsize
    while(kl_new > kl - alpha*step*decr) {
      step <- step*beta
      q$eta <- d$x%*%(q$b - step*db)
      q$mu <- family_kl$linkinv(q$eta + d$offset)
      q$dis <- family_kl$dis(p, d, q)
      kl_new <- family_kl$kl(p, d, q)
    }

    q$b <- q$b - step*db
    kl <- kl_new

    # check convergence
    if(step < 1e-8) break

    it <- it + 1
  }

  if(it >= max_it)
    warning(paste0("Maximum number of iterations reached."))

  res <- list(kl = kl, b = q$b)
  if(family_kl$family %in% c('gaussian','Gamma')) res$dis <- q$dis
  res
}
