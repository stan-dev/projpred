#' GLM projection to submodel.
#'
#' @param fit Either a stanreg object or a list with the model matrix, the weight vector
#'  a \code{\link{family}}-object and the dispersion parameter, if applicable.
#' @param n_out Number of samples used in the projection.
#'  If not provided, \code{min(400, ncol(b))} is used.
#' @param n_sel Number of samples used in the variable selection.
#'  If not provided, average of mu is used instead of the samples.
#' @param d Maximum number of features in the projection (incl. intercept).
#'  If not provided, \code{ncol(x)} - 1 is used.
#' @param mc.cores Number of cores used when averaging KL divergence over samples.
#'  If not provided, defaults to \code{getOption("loo.cores", parallel::detectCores())}.
#'
#' @return A list with class \code{'glmproj'} containing
#' sequence of features added, distributions of the parameters
#' of the submodels, the corresponding KL divergences and the model family.
#'
#' @importFrom rstanarm get_x
#'
#' @export glm_proj

glm_proj <- function(fit, n_out, n_sel, d, mc.cores) {
  UseMethod('glm_proj')
}

#' @export

glm_proj.stanreg <- function(fit, n_out, n_sel, d, mc.cores) {

  x <- unname(get_x(fit))
  b_p <- unname(t(apply(as.array(fit$stanfit), 3, rbind)[, colnames(get_x(fit)), drop = FALSE]))
  w <- get_weights(fit)
  dis_p <- get_dispersion(fit)
  family <- fit$family

  res <- glm_proj_fsel(x, b_p, w, dis_p, n_out, n_sel, d, family, mc.cores)
  res$family <- family
  structure(res, class = 'glmproj')
}

#' @export

glm_proj.list <- function(fit, n_out, n_sel, d, mc.cores) {

  if(any(!(c('x','b','family') %in% names(fit))))
    stop('fit must contain elements x, b and family')

  x <- fit$x
  b_p <- fit$b
  w <- fit$w
  dis_p <- fit$dis
  family <- fit$family

  if(is.null(dis_p)) dis_p <- NA
  if(is.null(w)) w <- rep(1, nrow(x))
  if(family$family %in% c('gaussian','Gamma') && is.na(dis_p))
    stop('dispersion parameter not provided')

  res <- glm_proj_fsel(x, b_p, w, dis_p, n_out, n_sel, d, family, mc.cores)
  res$family <- family
  structure(res, class = 'glmproj')
}
