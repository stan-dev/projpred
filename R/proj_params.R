#' Project parameters
#'
#' Find parameters that minimize the KL divergence between the full model
#' and the submodel.
#'
#' @param mu_p Fitted values of the full model.
#' @param x A model matrix.
#' @param b_p Sampled estimates of the coefficient of the full model.
#' @param w Observation weights.
#' @param dis_p Dispersion parameter of the full model.
#' @param funs List of family-specific functions for the NR.
#' @param chosen Indices for features used in the submodel.
#' @param cores Number of cores used.

proj_params <- function(mu_p, x, b_p, w, dis_p, funs, chosen, mc.cores) {
  d <- length(chosen)
  inds <- 1:d

  helperf <- function(ind) kl_over_samples(mu_p, x[, chosen[1:ind], drop = F], b_p[chosen[1:ind], , drop = F], w, dis_p, funs)

  if (mc.cores == 1 ) {
    l_res <- lapply(inds,helperf)
  } else {
    if (.Platform$OS.type == "windows") {
      cl <- makePSOCKcluster(mc.cores)
      on.exit(stopCluster(cl))
      l_res <- parLapply(cl, inds, helperf)
    } else {
      l_res <- mclapply(inds, helperf, mc.cores = mc.cores)
    }
  }

  ul_res <- unlist(l_res, recursive = F)

  res <- list(kl = unname(unlist(ul_res[names(ul_res) == 'kl'])),
              b = unname(ul_res[names(ul_res) == 'b']))

  if(funs$family %in% c('gaussian','Gamma'))
    res$dis <- unname(ul_res[names(ul_res) == 'dis'])

  res
}
