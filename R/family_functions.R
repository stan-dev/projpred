#' Helper functions for the NR-algorithm
#'
#' Contains functions for calculating KL divergence and its first and second
#' derivatives for different families and link functions.
#'
#' The NR-algorithm is general in the sense that in order to add a new model /
#' link function one simply needs to implement corresponding KL-divergence,
#' it's first and second derivatives + first and second derivatives of the
#' mu-parameter w.r.t eta.

kl_gauss <- function(mu_p, x, mu_q, w, dis_p, dis_q) {
  log(dis_q)-log(dis_p)
}

dkl_gauss <- function(mu_p, x, mu_q, eta_q, w, dme) {
  list(f = crossprod(x, mu_q - mu_p), s = crossprod(x))
}

kl_bin <- function(mu_p, x, mu_q, w, dis_p, dis_q) {
  ninv <- 1/length(mu_p)
  sum(w*(mu_p*(log(mu_p)-log(mu_q)) + (1-mu_p)*(log1p(-mu_p)-log1p(-mu_q))))*ninv
}

dkl_bin <- function(mu_p, x, mu_q, eta_q, w, dme) {
  pq <- mu_p/mu_q
  pq1 <- (1-mu_p)/(1-mu_q)
  ninv <- 1/length(mu_p)
  tmp1 <- -w*(pq - pq1)*dme$f*ninv
  tmp2 <- w*((pq/mu_q + pq1/(1-mu_q))*dme$f^2 - (pq - pq1)*dme$s)*ninv

  list(f = crossprod(x,tmp1), s = crossprod(x*drop(tmp2),x))
}

kl_poiss <- function(mu_p, x, mu_q, w, dis_p, dis_q) {
  ninv <- 1/length(mu_p)
  sum(- mu_p + mu_q + mu_p*(log(mu_p)-log(mu_q)))*ninv
}

dkl_poiss <- function(mu_p, x, mu_q, eta_q, w, dme) {
  pq <- mu_p/mu_q
  ninv <- 1/length(mu_p)
  tmp1 <- (1-pq)*dme$f*ninv
  tmp2 <- (pq/mu_q*dme$f^2 + (1-pq)*dme$s)*ninv

  list(f = crossprod(x,tmp1), s = crossprod(x*drop(tmp2),x))
}

# d mu / d eta
dme_logit <- function(mu_q, eta_q) {
  tmp <- mu_q*(1-mu_q)
  list(f = tmp, s = tmp*(1-2*mu_q))
}
dme_probit <- function(mu_q, eta_q) {
  tmp <- dnorm(eta_q)
  list(f = tmp, s = tmp*(-eta_q))
}
dme_log <- function(mu_q, eta_q) {
  list(f = mu_q, s = mu_q)
}
dme_id <- function(mu_q, eta_q) {
  tmp <- rep(1, length(mu_q))
  list(f = tmp, s = tmp)
}

# dispersions
disp_na <- function(...) NA

disp_ga <- function(mu_p, x, mu_q, dis_p) {
  ninv <- 1/length(mu_p)

  sqrt(dis_p^2 + sum((mu_p-mu_q)^2)*ninv)
}

family_kls <- function(family) {
  famfs <- switch(family$family,
      'gaussian' = list(kl = kl_gauss, dkl = dkl_gauss, dis = disp_ga),
      'binomial' = list(kl = kl_bin, dkl = dkl_bin, dis = disp_na),
      'poisson' = list(kl = kl_poiss, dkl = dkl_poiss, dis = disp_na)
  )
  linkfs <- switch(family$link,
      'logit' = list(dme = dme_logit),
      'log' = list(dme = dme_log),
      'identity' = list(dme = dme_id)
  )
  c(famfs, linkfs, family)
}
