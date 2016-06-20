#' GLM projection to submodel using forward selection
#'
#' A helper function that calls the forward selection
#' and the projection functions. In order to use another
#' selection method, replace chosen <- fsel(...) with another
#' function.
#'
#' @param \code{x} A model matrix (incl. intercept).
#' @param \code{b_p} Sampled estimates of the coefficient of the full model.
#' @param \code{w} Observation weights.
#' @param \code{dis_p} dispersion parameter of the full model.
#' @param \code{family} A \code{\link{family}}-object.
#' @param \code{args} For possible additional arguments, see the
#'  documentation of \code{\link{glm_proj}}.
#'

glm_proj_fsel <- function(x, b_p, w, dis_p, family, args) {

  # number of samples
  ns <- ncol(b_p)

  # set additional argumets to defaults if they are missing etc.
  n_out <- min(ifelse(is.null(args$n_out), 800, args$n_out), ns)
  n_sel <- min(ifelse(is.null(args$n_sel), 400, args$n_sel), ns)
  d <- min(ncol(x)-1, args$d)
  avg <- ifelse(is.null(args$avg), F, args$avg)
  if(avg) n_sel <- ns
  glmproj.cores <- ifelse(is.null(args$glmproj.cores),
                          getOption('glmproj.cores', parallel::detectCores()),
                          args$glmproj.cores)
  if(ncol(x) < 3)
    stop('data must have at least 3 features.')
  if(!(family$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(family$family, 'family not yet supported'))

  # helper functions (kl, derivatives etc.)
  funs <- family_kls(family)

  # fitted values of the full model
  mu_p <- family$linkinv(x%*%b_p)

  # Set sample indices to be used with forward selection and final projection
  sind <- round(seq(1, ns, length.out  = n_sel))
  oind <- round(seq(1, ns, length.out  = n_out))

  # find 'optimal' sequence with forward selection
  chosen <- fsel(mu_p[, sind, drop = F], x, b_p[, sind, drop = F], w, dis_p[sind], funs, avg, d, glmproj.cores)

  # project parameters to submodels
  proj <- proj_params(mu_p[, oind, drop = F], x, b_p[, oind, drop = F], w, dis_p[oind], funs, chosen, glmproj.cores)
  proj$family <- family

  c(chosen = list(chosen), proj)
}
