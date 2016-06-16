#' GLM projection to submodel using forward selection
#'
#' A helper function that calls the forward selection
#' and the projection functions.
#'
#' @param x A model matrix.
#' @param b_p Sampled estimates of the coefficient of the full model.
#' @param w Observation weights.
#' @param dis_p dispersion parameter of the full model.
#' @param n_out Number of samples used in the projection.
#' @param n_sel Number of samples used in the variable selection.
#' @param d Maximum number of features in the projection (incl. intercept).
#' @param family A \code{\link{family}}-object.
#' @param mc.cores Number of cores used.

glm_proj_fsel <- function(x, b_p, w, dis_p, n_out, n_sel, d, family, mc.cores) {

  if(ncol(x) < 2)
    stop('data must have at least 2 features.')
  if(!(family$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(family$family, 'family not yet supported'))

  if(missing(d)) d <- ncol(x) - 1
  if(missing(mc.cores)) mc.cores <- getOption('mc.cores', parallel::detectCores())

  # helper functions (kl, derivatives etc.)
  funs <- family_kls(family)

  # fitted values of the full model
  mu_p <- family$linkinv(x%*%b_p)

  # Set sample indices used in forward selection and final projection
  ns <- ncol(b_p)
  avg <- missing(n_sel)
  if(avg) n_sel <- ns
  sind <- round(seq(1, ns, length.out  = n_sel))

  # find 'optimal' sequence with forward selection
  chosen <- fsel(mu_p[, sind, drop = F], x, b_p[, sind, drop = F], w, dis_p[sind], funs, avg, d, mc.cores)

  if(missing(n_out)) n_out <- min(ns, 400)
  oind <- round(seq(1, ns, length.out  = n_out))

  # project parameters to submodels
  proj <- proj_params(mu_p[, oind, drop = F], x, b_p[, oind, drop = F], w, dis_p[oind], funs, chosen, mc.cores)

  c(chosen = list(chosen), proj)
}
