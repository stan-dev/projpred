#' Forward selection
#'
#' Performs forward selection on the data and returns
#' the 'optimal' order of the parameters.
#'
#' @param mu_p Fitted values of the full model.
#' @param x A model matrix.
#' @param b_p Sampled estimates of the coefficient of the full model.
#' @param w Observation weights.
#' @param dis_p dispersion parameter of the full model.
#' @param funs List of family-specific functions for the NR.
#' @param avg If TRUE, KL divergence is calculated using average of the full model, rather than individual samples.
#' @param d Maximum number of features in the projection (incl. intercept).
#' @param mc.cores Number of cores used.
#'
#' @importFrom parallel mclapply makePSOCKcluster stopCluster parLapply

fsel <- function(mu_p, x, b_p, w, dis_p, funs, avg, d, mc.cores) {

  chosen <- 1
  cols <- 1:ncol(x)
  notchosen <- setdiff(cols, chosen)

  if(avg) {
    mu_p_mean <- rowMeans(mu_p)
    b_p_mean <- rowMeans(b_p)
    dis_p_mean <- sqrt(mean(dis_p^2))
    helperf <- function(ind) NR(mu_p_mean, x[, c(chosen, ind), drop = F], b_p_mean[c(chosen, ind)], w, dis_p_mean, funs)$kl
  } else {
    helperf <- function(ind) kl_over_samples(mu_p, x[, c(chosen, ind), drop = F], b_p[c(chosen, ind), , drop = F], w, dis_p, funs)$kl
  }

  # start adding variables one at a time
  for (k in 2:d) {

    # with avg there the overhead of mclapply is too large
    # (this might not be the case if d >> 100)
    if (mc.cores == 1 || avg) {
      l_res <- lapply(notchosen, helperf)
    } else {
      if (.Platform$OS.type == 'windows') {
        cl <- makePSOCKcluster(mc.cores)
        on.exit(stopCluster(cl))
        l_res <- parLapply(cl, notchosen, helperf)
      } else {
        l_res <- mclapply(notchosen, helperf, mc.cores = mc.cores)
      }
    }

    imin <- which.min(unlist(l_res))

    chosen <- c(chosen, notchosen[imin])
    notchosen <- setdiff(cols, chosen)
  }

  chosen
}
