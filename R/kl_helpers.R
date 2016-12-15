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
  kl_binom <- function(p_full, d_train, p_sub) {
    mean(fam$dev.resids(p_full$mu, p_sub$mu, d_train$weights))/2
  }
  kl_poiss <- function(p_full, d_train, p_sub, family) {
    mean(fam$dev.resids(p_full$mu, p_sub$mu, d_train$weights))/2
  }
  kl_gauss <- function(p_full, d_train, p_sub) log(p_sub$dis) - log(p_full$dis)
  kl_gamma <- function(p_full, d_train, p_sub) {
    mean(d_train$weights*(
      p_sub$dis*(log(p_full$dis)-log(p_sub$dis)+log(p_sub$mu)-log(p_full$mu)) +
      digamma(p_full$dis)*(p_full$dis - p_sub$dis) - lgamma(p_full$dis) +
        lgamma(p_sub$dis) + p_full$mu*p_sub$dis/p_sub$mu - p_full$dis))
  }

  # for gaussian dispersion is sigma and for gamma it is the shape param
  dis_na <- function(p_full, d_train, p_sub) rep(1, length(p_full$dis))
  dis_gauss <- function(p_full, d_train, p_sub) {
    sqrt(mean(d_train$weights*(p_full$mu - p_sub$mu)^2) + p_full$dis^2)
  }

  dis_gamma <- function(p_full, d_train, p_sub) {
    mean(d_train$weights*((p_full$mu - p_sub$mu)/
                            fam$mu.eta(fam$linkfun(p_sub$mu))^2))
  }

  # log likelihoods
  ll_binom <- function(mu, dis, y, weights) dbinom(weights*y, weights,mu,log=T)
  ll_poiss <- function(mu, dis, y, weights) weights*dpois(y, mu, log=T)
  ll_gauss <- function(mu, dis, y, weights) weights*dnorm(y, mu, dis, log=T)
  ll_gamma <- function(mu, dis, y, weights) weights*dgamma(y, dis,dis/mu,log=T)

  # return the family object with the correct function handles
  c(switch(fam$family,
      'binomial' = list(kl = kl_binom, ll_fun = ll_binom, dis_fun = dis_na),
      'poisson' = list(kl = kl_poiss, ll_fun = ll_poiss, dis_fun = dis_na),
      'gaussian' = list(kl = kl_gauss, ll_fun = ll_gauss, dis_fun = dis_gauss),
      'Gamma' = list(kl = kl_gamma, ll_fun = ll_gamma, dis_fun = dis_gamma))
  , list(mu_fun = function(x, alpha, beta, offset, intercept) {
      fam$linkinv(cbind(1*intercept, x)%*%rbind(alpha, beta) + offset)
    }), fam)

}



