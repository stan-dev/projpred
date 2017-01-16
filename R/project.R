#' @export
project <- function(object, nv, ...) {
  UseMethod('project')
}

#' @export
project.stanreg <- function(object, nv, nc = NULL, ns = NULL, intercept = NULL, ...) {
	
	# TODO, IMPLEMENT THE PROJECTION WITH AN ARBITRARY VARIABLE COMBINATION
	
	if (is.null(ns) && is.null(nc))
		# by default project with clusters
		nc <- 50
	
	if(is.null(nv)) stop('nv not provided')
	.validate_for_varsel(object)
	if(!('varsel' %in% names(object)))
		stop(paste('The stanreg object doesn\'t contain information about the ',
					'variable selection. Run the variable selection first.'))

	vars <- .extract_vars(object)
	# if(ns > ncol(vars$beta)) {
	# 	warning(paste0('Setting the number of samples to ', ncol(vars$beta),'.'))
	# 	ns <- ncol(vars$beta)
	# }
	if(is.null(intercept)) intercept <- vars$intercept

	family_kl <- kl_helpers(family(object))

	if(max(nv) > length(object$varsel$chosen))
		stop(paste('Cannot perform the projection with', max(nv), 'variables,',
				'because the variable selection has been run only up to',
				length(object$varsel$chosen), 'variables.'))

	# training data
	d_train <- .get_traindata(object)
	
	# get the clustering or subsample
	p_full <- .get_refdist(object, ns=ns, nc=nc)
	
	object$proj <- list(
		p_sub = .get_submodels(object$varsel$chosen, nv, family_kl,
								p_full, d_train, intercept),
		intercept = intercept)

	object
}
