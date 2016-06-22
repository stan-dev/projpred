#' KL divergence over samples
#'
#' Calculates the KL divergence over ns samples. For gaussian likelihood,
#' uses analytical solution to make computation a lot faster.
#'
#' @param \code{mu_p} Fitted values of the full model.
#' @param \code{x} A model matrix of the selected variables.
#' @param \code{b_p} Sampled estimates of the coefficients of the selected variables.
#' @param \code{w} Observation weights.
#' @param \code{dis_p} dispersion parameter of the full model.
#' @param \code{funs} Model-specific helper functions.

kl_over_samples <- function(mu_p, x, b_p, w, dis_p, funs) {

  if(funs$family == 'gaussian') {
    b_q <- solve(crossprod(x), crossprod(x, mu_p))
    dis_q <- sqrt(dis_p^2 + colSums((mu_p - x%*%b_q)^2)/nrow(x))
    kl <- sum(log(dis_q) - log(dis_p))/length(dis_q)
    res <- list(kl = kl, b = b_q, dis = dis_q)

  } else {
    helperf <- function(ind) NR(mu_p[,ind], x, b_p[,ind], w, dis_p[ind], funs)
    inds <- 1:ncol(mu_p)

    l_res <- lapply(inds,helperf)

    ul_res <- unlist(l_res, recursive = F)

    res <- list(kl = mean(unlist(ul_res[names(ul_res) == 'kl'])),
                b = unname(do.call('cbind',ul_res[names(ul_res) == 'b'])))
    if(funs$family %in% c('Gamma'))
      res$dis <- unname(unlist(ul_res[names(ul_res) == 'dis']))
  }

  res
}
