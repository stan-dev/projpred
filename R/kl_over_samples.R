#' KL divergence over samples
#'
#' Calculates the KL divergence over ns samples.
#'
#' @param mu_p Fitted values of the full model.
#' @param x A model matrix of the selected variables.
#' @param b_p Sampled estimates of the coefficients of the selected variables.
#' @param w Observation weights.
#' @param dis_p dispersion parameter of the full model.
#' @param funs List of family-specific functions for the NR.

kl_over_samples <- function(mu_p, x, b_p, w, dis_p, funs) {

  helperf <- function(ind) NR(mu_p[,ind], x, b_p[,ind], w, dis_p[ind], funs)
  inds <- 1:ncol(mu_p)

  l_res <- lapply(inds,helperf)

  ul_res <- unlist(l_res, recursive = F)

  res <- list(kl = mean(unlist(ul_res[names(ul_res) == 'kl'])),
              b = unname(do.call('cbind',ul_res[names(ul_res) == 'b'])))

  if(funs$family %in% c('gaussian','Gamma'))
    res$dis <- unname(unlist(ul_res[names(ul_res) == 'dis']))

  res
}
