#' @export
project <- function(object, nv = NULL, nc = NULL, ns = NULL, intercept = NULL, seed=NULL, ...) {
	
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
