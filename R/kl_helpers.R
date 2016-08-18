#' Model-specific helper functions.
#'
#' Contains functions for calculating KL divergence and its first and second
#' derivatives for different families and link functions.
#'
#' The NR-algorithm is general in the sense that in order to add a new model /
#' link function one simply needs to implement corresponding KL-divergence,
#' it's first and second derivatives + first and second derivatives of the
#' mu-parameter w.r.t eta.

kl_gauss <- function(p, d, q) {
  log(q$dis) - log(p$dis)
}

dkl_gauss <- function(p, d, q, dme) {
  invdenom <- 1/(length(q$mu)*q$dis^2)
  list(f = crossprod(d$x, q$mu - p$mu)*invdenom, s = crossprod(d$x)*invdenom)
}

kl_binom <- function(p, d, q) {
  sum(d$w*(p$mu*(log(p$mu)-log(q$mu)) +
           (1-p$mu)*(log1p(-1*p$mu)-log1p(-1*q$mu))))/length(q$mu)
}

dkl_binom <- function(p, d, q, dme) {
  ninv <- 1/length(q$mu)
  pq <- p$mu/q$mu
  pq1 <- (1-p$mu)/(1-q$mu)

  list(f = ninv*crossprod(d$x,d$w*(pq1 - pq)*dme$f),
       s = ninv*crossprod(d$x*drop(d$w*((pq1 - pq)*dme$s) +
                                    (pq/q$mu + pq1/(1-q$mu))*dme$f^2), d$x))
}

kl_poiss <- function(p, d, q) {
  sum(q$mu - p$mu + p$mu*(log(p$mu)-log(q$mu)))/length(q$mu)
}

dkl_poiss <- function(p, d, q, dme) {
  ninv <- 1/length(q$mu)
  pq <- p$mu/q$mu

  list(f = ninv*crossprod(d$x,(1-pq)*dme$f),
       s = ninv*crossprod(d$x*drop((pq/q$mu*dme$f^2 + (1-pq)*dme$s)),d$x))
}

# d mu / d eta
dme_logit <- function(q) {
  tmp <- q$mu*(1-q$mu)
  list(f = tmp, s = tmp*(1-2*q$mu))
}
dme_probit <- function(q) {
  tmp <- dnorm(q$eta)
  list(f = tmp, s = -1*tmp*(q$eta))
}
dme_log <- function(q) {
  list(f = q$mu, s = q$mu)
}
dme_id <- function(q) {
  tmp <- 1
  list(f = tmp, s = tmp)
}

# dispersions
disp_na <- function(...) NA

disp_ga <- function(p, d, q) {
  sqrt(p$dis^2 + sum((p$mu-q$mu)^2)/length(q$mu))
}

# log likelihoods for test data
ll_gauss <- function(mu, dis, y, w) {
  dnorm(y, mean = mu, sd = rep(dis, each = length(y)), log = T)
}

ll_binom <- function(mu, dis, y, w) {
  dbinom(y, size = w, prob = mu, log = T)
}

ll_poiss <- function(mu, dis, y, w) {
  dpois(y, mu, log = T)
}

kl_helpers <- function(family) {
  famfs <- switch(family$family,
      'gaussian' = list(kl = kl_gauss, dkl = dkl_gauss,
                        dis = disp_ga, ll_fun = ll_gauss),
      'binomial' = list(kl = kl_binom, dkl = dkl_binom,
                        dis = disp_na, ll_fun = ll_binom),
      'poisson' = list(kl = kl_poiss, dkl = dkl_poiss,
                       dis = disp_na, ll_fun = ll_poiss)
  )
  linkfs <- switch(family$link,
      'logit' = list(dme = dme_logit),
      'probit' = list(dme = dme_probit),
      'log' = list(dme = dme_log),
      'identity' = list(dme = dme_id)
  )
  c(famfs, linkfs, family)
}
