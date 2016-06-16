#' Newton-Rhapson
#'
#' Minimizes the KL divergence between the full model and the submodel
#' using Newton-Rhapson algorithm.
#'
#' @param mu_p Fitted values of the full model.
#' @param x A model matrix of the selected variables.
#' @param b_p Sampled estimates of the coefficient of the full model.
#' @param w Observation weights.
#' @param dis_p dispersion parameter of the full model.
#' @param funs List of family-specific functions for the NR.
#' @param max_it Maximum number of iterations for the algorithm. Defaults to 50.
#' @param eps Tolerance, when derivative to any direction is at most eps, algorithm stops. Defaults to 1e-10.

NR <- function(mu_p, x, b_p, w, dis_p, funs, max_it = 50, eps = 1e-10) {

  # initial estimate for b_q
  b_q <- b_p

  it <- 1

  eta_q <- x%*%b_q
  mu_q <- funs$linkinv(x%*%b_q)
  dis_q <- funs$dis(mu_p, x, mu_q, dis_p)
  kl <- funs$kl(mu_p, x, mu_q, w, dis_p, dis_q)

  alpha <- 0.1
  beta <- 0.5

  while(it < max_it) {

    dme <- funs$dme(mu_q, eta_q)
    dkl <- funs$dkl(mu_p, x, mu_q, eta_q, w, dme)
    db_q <- solve(dkl$s, dkl$f)
    incr <- crossprod(dkl$f,db_q)

    # check convergence
    if(incr*0.5 < eps) break

    step <- 1

    eta_q <- x%*%(b_q - step*db_q)
    mu_q <- funs$linkinv(eta_q)
    dis_q <- funs$dis(mu_p, x, mu_q, dis_p)
    kl_new <- funs$kl(mu_p, x, mu_q, w, dis_p, dis_q)

    # find stepsize
    while(kl_new > kl + alpha*step*incr) {
      step <- step*beta
      eta_q <- x%*%(b_q - step*db_q)
      mu_q <- funs$linkinv(eta_q)
      dis_q <- funs$dis(mu_p, x, mu_q, dis_p)
      kl_new <- funs$kl(mu_p, x, mu_q, w, dis_p, dis_q)
    }

    b_q <- b_q - step*db_q
    kl <- kl_new

    # check convergence
    if(step < 1e-8) break

    it <- it + 1
  }

  if(it >= max_it) warning(paste0('Maximum number of iterations reached.'))

  list(kl = kl, b = b_q, dis = dis_q)
}
