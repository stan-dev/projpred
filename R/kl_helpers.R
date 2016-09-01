#' Model-specific helper functions.
#'
#' Contains functions for calculating KL divergence and its first and second
#' derivatives for different families and link functions.
#'
#' The NR-algorithm is general in the sense that in order to add a new model /
#' link function one simply needs to implement corresponding KL-divergence,
#' it's first and second derivatives + first and second derivatives of the
#' mu-parameter w.r.t eta.
#'
#' Naming convention described in fsel.R.


# add the functions to the family-object
kl_helpers <- function(family) {
  famfs <- switch(family$family,
      'gaussian' = list(kl = .kl_gauss, dkl = .dkl_gauss,
                        dis = .disp_ga, ll_fun = .ll_gauss),
      'binomial' = list(kl = .kl_binom, dkl = .dkl_binom,
                        dis = .disp_na, ll_fun = .ll_binom),
      'poisson' = list(kl = .kl_poiss, dkl = .dkl_poiss,
                       dis = .disp_na, ll_fun = .ll_poiss)
  )
  linkfs <- switch(family$link,
      'logit' = list(dme = .dme_logit),
      'probit' = list(dme = .dme_probit),
      'log' = list(dme = .dme_log),
      'identity' = list(dme = .dme_id)
  )
  c(famfs, linkfs, family)
}

# kl-divergences
.kl_gauss <- function(p_full, d_train, p_sub) {
  log(p_sub$dis) - log(p_full$dis)
}

.kl_binom <- function(p_full, d_train, p_sub) {
  sum(d_train$w*(p_full$mu*(log(p_full$mu)-log(p_sub$mu)) +
                 (1-p_full$mu)*(log1p(-1*p_full$mu) - log1p(-1*p_sub$mu))
                 ))/length(p_sub$mu)
}

.kl_poiss <- function(p_full, d_train, p_sub) {
  sum(p_sub$mu - p_full$mu + p_full$mu*(log(p_full$mu)-log(p_sub$mu)))/length(p_sub$mu)
}

# the gradients and the hessians of the kl divergences
.dkl_gauss <- function(p_full, d_train, p_sub, dme) {
  list(
    g = crossprod(d_train$x, p_sub$mu - p_full$mu)/(length(p_sub$mu)*p_sub$dis^2),
    h = crossprod(d_train$x)/(length(p_sub$mu)*p_sub$dis^2))
}

.dkl_binom <- function(p_full, d_train, p_sub, dme) {
  list(
    g = crossprod(d_train$x, d_train$w*((1-p_full$mu)/(1-p_sub$mu) - p_full$mu/p_sub$mu)*dme$g) / length(p_sub$mu),
    h = crossprod(d_train$x*d_train$w*drop(((1-p_full$mu)/(1-p_sub$mu) - p_full$mu/p_sub$mu)*dme$h +
                                             (p_full$mu/p_sub$mu^2 + (1-p_full$mu)/(1-p_sub$mu)^2)*dme$g^2),
                  d_train$x)/length(p_sub$mu))
}

.dkl_poiss <- function(p_full, d_train, p_sub, dme) {
  list(
    g = crossprod(d_train$x, (1-p_full$mu/p_sub$mu)*dme$g)/length(p_sub$mu),
    h = crossprod(d_train$x*drop((p_full$mu/p_sub$mu^2*dme$g^2 + (1-p_full$mu/p_sub$mu)*dme$h)),
                  d_train$x)/length(p_sub$mu))
}

# d mu / d eta
.dme_logit <- function(p_sub) {
  list(g = p_sub$mu*(1-p_sub$mu),
       h = p_sub$mu*(1-p_sub$mu)*(1-2*p_sub$mu))
}
.dme_probit <- function(p_sub) {
  tmp <- dnorm(p_sub$eta)
  list(g = tmp,
       h = -1*tmp*(p_sub$eta))
}
.dme_log <- function(p_sub) {
  list(g = p_sub$mu,
       h = p_sub$mu)
}
.dme_id <- function(p_sub) {
  list(g = rep(1, length(p_sub$mu)),
       h = rep(1, length(p_sub$mu)))
}

# dispersions
.disp_na <- function(p_full, d_train, p_sub) rep(1, length(p_sub$mu))

.disp_ga <- function(p_full, d_train, p_sub) sqrt(p_full$dis^2 + sum((p_full$mu-p_sub$mu)^2)/length(p_sub$mu))

# log likelihoods
.ll_gauss <- function(mu, dis, y, w) {
  dnorm(y, mean = mu, sd = rep(dis, each = length(y)), log = T)
}

.ll_binom <- function(mu, dis, y, w) {
  dbinom(w*y, size = w, prob = mu, log = T)
}

.ll_poiss <- function(mu, dis, y, w) {
  dpois(y, mu, log = T)
}

