#' Model-specific helper functions.
#'
#' \code{kl_helpers(fam)} returns a family object augmented with functions for
#' KL divergence, log predictive density and dispersion.
#'
#' Naming convention described in fsel.R.
#'
#' Missing: Quasi-families not implemented. If dis_gamma is the correct shape
#' parameter for projected Gamma regression, everything should be OK for gamma.

kl_helpers <- function(fam) {

  # define the functions for all families but
  # return only the ones that are needed.

  # kl-divergences
  # for binomial and poisson it is the mean of the dev.resids divided by 2
  kl_dev <- function(p_full, data, p_sub) {
    if(NCOL(p_full$mu)>1) {
      w <- rep(data$weights, NCOL(p_full$mu))
      colMeans(fam$dev.resids(p_full$mu, p_sub$mu, w))/2
    } else {
      mean(fam$dev.resids(p_full$mu, p_sub$mu, data$weights))/2
    }
  }
  kl_gauss <- function(p_full, data, p_sub) log(p_sub$dis) - log(p_full$dis)
  kl_gamma <- function(p_full, data, p_sub) {
    mean(data$weights*(
      p_sub$dis*(log(p_full$dis)-log(p_sub$dis)+log(p_sub$mu)-log(p_full$mu)) +
        digamma(p_full$dis)*(p_full$dis - p_sub$dis) - lgamma(p_full$dis) +
        lgamma(p_sub$dis) + p_full$mu*p_sub$dis/p_sub$mu - p_full$dis))
  }

  # for gaussian dispersion is sigma and for gamma it is the shape param
  dis_na <- function(p_full, data, p_sub) rep(1, length(p_full$dis))
  dis_gauss <- function(p_full, data, p_sub) {
    sqrt(mean(data$weights*(p_full$mu - p_sub$mu)^2) + p_full$dis^2)
  }

  dis_gamma <- function(p_full, data, p_sub) {
    mean(data$weights*((p_full$mu - p_sub$mu)/
                            fam$mu.eta(fam$linkfun(p_sub$mu))^2))
  }

  # log likelihoods
  ll_binom <- function(mu, dis, y, weights=1) dbinom(weights*y, weights, mu, log=T)
  ll_poiss <- function(mu, dis, y, weights=1) weights*dpois(y, mu, log=T)
  ll_gauss <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*dnorm(y, mu, dis, log=T)
  }
  ll_gamma <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*dgamma(y, dis, dis/matrix(mu), log=T)
  }
  
  # functions to sample from posterior predictive distribution
  ppd_gauss <- function(mu, dis, weights = 1) rnorm(length(mu), mu, sig)
  ppd_binom <- function(mu, dis, weights = 1) rbinom(length(mu), weights, mu)
  ppd_poiss <- function(mu, dis, weights = 1) rpois(length(mu), mu)
  ppd_gamma <- function(mu, dis, weights = 1) rgamma(length(mu), dis, dis/mu)
  

  # function for computing mu = E(y)
  mu_fun <- function(x, alpha, beta, offset) {
    if (!is.matrix(x)) stop('x must be a matrix.')
    if (!is.matrix(beta)) stop('beta must be a matrix')
    fam$linkinv(cbind(1, x) %*% rbind(alpha, beta) + offset)
  }

  # return the family object with the correct function handles
  c(switch(fam$family,
           'binomial' = list(kl = kl_dev, ll_fun = ll_binom, dis_fun = dis_na,
                             ppd_fun = ppd_binom),
           'poisson' = list(kl = kl_dev, ll_fun = ll_poiss, dis_fun = dis_na,
                            ppd_fun = ppd_poisson),
           'gaussian' = list(kl = kl_gauss, ll_fun = ll_gauss, dis_fun = dis_gauss,
                             ppd_fun = ppd_gauss),
           'Gamma' = list(kl = kl_gamma, ll_fun = ll_gamma, dis_fun = dis_gamma,
                          ppd_fun = ppd_gamma)),
    list(mu_fun = mu_fun), fam)

}



