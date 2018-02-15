#' Projection to submodels of selected sizes
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param object The object returned by \link[=varsel]{varsel} or
#' \link[=cv_varsel]{cv_varsel}.
#' @param nv Number of variables in the submodel (the variable combination is taken from the
#' \code{varsel} information). If a list, then the projection is performed for each model size.
#' Default is all model sizes up to the maximum number of variables in \code{varsel}.
#'  Ignored if \code{vind} is specified.
#' @param vind Variable indices onto which the projection is done. If specified, \code{nv} is ignored.
#' @param ns Number of samples to be projected. Ignored if \code{nc} is specified. Default is 400.
#' @param nc Number of clusters in the clustered projection. 
#' @param intercept Whether to use intercept. Default is \code{TRUE}.
#' @param seed A seed used in the clustering (if \code{nc!=ns}). Can be used
#' to ensure same results every time.
#' @param regul Amount of regularization in the projection. Usually there is no need for
#' regularization, but sometimes for some models the projection can be ill-behaved and we
#' need to add some regularization to avoid numerical problems. Default is 1e-9.
#' @param ... Currently ignored.
#'
#' @return A list of submodels (or a single submodel if projection was performed onto
#' a single variable combination), each of which contain the following elements:
#' \describe{
#'  \item{\code{kl}}{The kl divergence from the full model to the submodel.}
#'  \item{\code{weights}}{Weights for each draw of the projected model.}
#'  \item{\code{dis}}{Draws from the projected dispersion parameter.}
#'  \item{\code{alpha}}{Draws from the projected intercept.}
#'  \item{\code{beta}}{Draws from the projected weight vector.}
#'  \item{\code{vind}}{The order in which the variables were added to the submodel.}
#'  \item{\code{intercept}}{Whether or not the model contains an intercept.}
#'  \item{\code{family_kl}}{A modified \code{\link{family}}-object.}
#' }
#'
#'
#' @examples
#' \donttest{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fit_v <- varsel(fit)
#' proj4 <- project(fit_v, nv = 4)
#' }
#'

#' @export
project <- function(object, nv = NULL, vind = NULL, ns = NULL, nc = NULL,
                    intercept = NULL, seed = NULL, regul=1e-6, ...) {

	if(!('varsel' %in% names(object)) && is.null(vind))
		stop(paste('The given object does not contain information about the ',
					'variable selection. Run the variable selection first, ',
                    'or provide the variable indices (vind).'))

    vars <- .extract_vars(object)

    if (is.null(vind)) {
        vind <- object$varsel$vind
    } else {
        nv <- length(vind) # if vind is given, nv is ignored (project only onto the given submodel)
    }

    # by default project with clusters
    if (is.null(ns) && is.null(nc))
      ns <- min(400, NCOL(vars$mu))
    # by default, run the projection up to the maximum number of variables
    # specified in the variable selection
    if (is.null(nv))
      nv <- c(0:length(vind))

	if (is.null(intercept))
	  intercept <- vars$intercept

	family_kl <- vars$fam

	if (max(nv) > length(vind))
	  stop(paste('Cannot perform the projection with', max(nv), 'variables,',
	             'because the variable selection has been run only up to',
	             length(object$varsel$vind), 'variables.'))

	# training data
	d_train <- .get_traindata(vars)

	# get the clustering or subsample
	p_full <- .get_refdist(vars, ns = ns, nc = nc, seed = seed)

	subm <- .get_submodels(vind, nv, family_kl, p_full, d_train, intercept, regul)

	# add family_kl
	proj <- lapply(subm, function(x) {
	  names(x$vind) <- sapply(x$vind, function(i, ch) names(ch)[which(ch == i)],
	                          object$varsel$vind)
	  x <- c(x, list(family_kl = family_kl), list(p_type = is.null(ns)))
	  class(x) <- 'projection'
	  return(x)
	 })

	# If only one model size, just return the proj instead of a list of projs
	.unlist_proj(proj)
}
