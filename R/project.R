#' Projection to a submodels of selected size.
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param object The object returned by \link[=varsel]{varsel} or
#' \link[=cv_varsel]{cv_varsel}.
#' @param nv Number of variables in the submodel. Can also be a list, in
#'  which the projection is performed for each of the sizes in the list. If
#'  \code{NULL}, the projection is performed for all possible submodel sizes.
#' @param ns Same as in \link[=varsel]{varsel}.
#' @param nc Same as in \link[=varsel]{varsel}.
#' @param intercept Same as in \link[=varsel]{varsel}.
#' @param seed A seed used in the clustering (if \code{nc!=ns}). Can be used
#' to ensure same results every time.
#' @param ... Currently ignored.
#'
#' @return The original fit-object object augmented with a field 'proj',
#' which is a list containing the following elements:
#' \describe{
#'  \item{\code{p_sub}}{A list of the projected parameters for each of the
#'  submodels.}
#'  \item{\code{intercept}}{A logical indicating whether the submodels contain
#'   intercept or not.}
#' }
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fit_v <- varsel(fit)
#' fit_p <- project(fit_v)
#' }
#'

#' @export
project <- function(object, nv = NULL, ns = NULL, nc = NULL, intercept = NULL, seed=NULL, ...) {

	# TODO, IMPLEMENT THE PROJECTION WITH AN ARBITRARY VARIABLE COMBINATION

	# .validate_for_varsel(object)
	if(!('varsel' %in% names(object)))
		stop(paste('The stanreg object doesn\'t contain information about the ',
					'variable selection. Run the variable selection first.'))

    vars <- .extract_vars(object)

	if (is.null(ns) && is.null(nc))
		# by default project with clusters
		nc <- min(50, NCOL(vars$mu))
	if (is.null(nv))
		# by default, run the projection up to the maximum number of variables
		# specified in the variable selection
		nv <- c(0:length(object$varsel$chosen))

	if(is.null(intercept)) intercept <- vars$intercept

	family_kl <- vars$fam

	if(max(nv) > length(object$varsel$chosen))
		stop(paste('Cannot perform the projection with', max(nv), 'variables,',
				'because the variable selection has been run only up to',
				length(object$varsel$chosen), 'variables.'))

	# training data
	d_train <- .get_traindata(object)

	# get the clustering or subsample
	p_full <- .get_refdist(object, ns=ns, nc=nc, seed=seed)

	object$proj <- list(
		p_sub = .get_submodels(object$varsel$chosen, nv, family_kl,
								p_full, d_train, intercept),
		intercept = intercept)

	object
}
