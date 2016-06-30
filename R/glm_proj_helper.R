#' GLM projection to submodel using forward selection
#'
#' A helper function that calls the forward selection
#' and the projection functions. In order to use another
#' selection method, replace chosen <- fsel(...) with another
#' function.
#'
#' @param \code{x} A model matrix.
#' @param \code{b_p} Sampled estimates of the coefficient of the full model.
#' @param \code{w} Observation weights.
#' @param \code{dis_p} dispersion parameter of the full model.
#' @param \code{family} A \code{\link{family}}-object.
#' @param \code{args} For possible additional arguments, see the
#'  documentation of \code{\link{glm_proj}}.
#'
#' @importFrom Matrix rankMatrix
#'

glm_proj_helper <- function(x, b_p, w, dis_p, family, args) {

  x <- unname(x)
  b_p <- unname(b_p)
  w <- unname(w)
  dis_p <- unname(dis_p)

  # number of samples
  ns <- ncol(b_p)

  # set additional argumets to defaults if they are missing etc.
  n_out <- min(ifelse(is.null(args$n_out), 800, args$n_out), ns)
  n_sel <- min(ifelse(is.null(args$n_sel), 400, args$n_sel), ns)
  d <- min(ncol(x)-1, args$d, rankMatrix(x))
  avg <- ifelse(is.null(args$avg), FALSE, args$avg)
  intercept <- ifelse(is.null(args$intercept), FALSE, args$intercept)
  if(avg) n_sel <- ns
  glmproj.cores <- ifelse(is.null(args$glmproj.cores),
                          getOption('glmproj.cores', parallel::detectCores()),
                          args$glmproj.cores)
  verbose <- ifelse(is.null(args$verbose), F, args$verbose)
  if(ncol(x) < 2)
    stop('data must have at least 2 features.')
  if(!(family$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(family$family, 'family not yet supported'))

  # helper functions (kl, derivatives etc.)
  funs <- kl_helpers(family)

  # fitted values of the full model
  mu_p <- family$linkinv(x%*%b_p)

  # Set sample indices to be used with forward selection and final projection
  sind <- round(seq(1, ns, length.out  = n_sel))
  oind <- round(seq(1, ns, length.out  = n_out))

  # find 'optimal' sequence with forward selection
  if(verbose) print(paste0('Starting forward selection for up to ', d, ' variables...'))
  chosen <- fsel(mu_p[, sind, drop = F], x, b_p[, sind, drop = F], w,
                 dis_p[sind], funs, avg, d, glmproj.cores, verbose, intercept)
  if(verbose) print('Done.')

  if(verbose) print(paste0('Projecting parameters to submodels of size 1 to ', d, '.'))
  proj <- proj_params(mu_p[, oind, drop = F], x, b_p[, oind, drop = F], w,
                      dis_p[oind], funs, chosen, glmproj.cores, verbose)
  if(verbose) print('Done.')
  proj$family <- family

  c(chosen = list(chosen), proj)
}
