#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit Either a \link[=stanreg-objects]{stanreg}-object or an object returned
#' by \link[=init_refmodel]{init_refmodel}.
#' @param d_test A test dataset, which is used to evaluate model performance.
#' If not provided, training data is used. Currently this argument is for internal use only.
#' @param method The method used in the variable selection. Possible options are
#' \code{'L1'} for L1-search and \code{'forward'} for forward selection.
#' Default is 'forward' if the number of variables in the full data is at most 20, and
#' \code{'L1'} otherwise.
#' @param ns Number of posterior draws used in the variable selection.
#'    Cannot be larger than the number of draws in the full model.
#'    Ignored if nc is set.
#' @param nc Number of clusters to use in the clustered projection.
#'    Overrides the \code{ns} argument. Defaults to 1.
#' @param nspred Number of samples used for prediction (after selection). Ignored if ncpred is given.
#' @param ncpred Number of clusters used for prediction (after selection). Default is 5.
#' @param nv_max Maximum number of varibles until which the selection is continued.
#'    Defaults to min(D, floor(0.4*n)) where n is the number of observations and
#'    D the number of variables.
#' @param intercept Whether to use intercept in the submodels. Defaults to TRUE.
#' @param penalty Vector determining the relative penalties or costs for the variables.
#' Zero means that those variables have no cost and will therefore be selected first,
#' whereas Inf means that those variables will never be selected. Currently works only 
#' if method == 'L1'. By default 1 for each variable.
#' @param verbose If TRUE, may print out some information during the selection.
#'    Defaults to FALSE.
#' @param lambda_min_ratio Ratio between the smallest and largest lambda in the L1-penalized search.
#' This parameter essentially determines how long the search is carried out, i.e., how large submodels
#' are explored. No need to change the default value unless the program gives a warning about this.
#' @param nlambda Number of values in the lambda grid for L1-penalized search. No need to change unless
#' the program gives a warning about this.
#' @param regul Amount of regularization in the projection. Usually there is no need for 
#' regularization, but sometimes for some models the projection can be ill-behaved and we
#' need to add some regularization to avoid numerical problems. Default is 1e-9.
#' @param ... Currently ignored.
#'
#'
#' @return The original fit-object object augmented with a field 'varsel',
#' which is a list containing the following elements:
#' \describe{
#'  \item{\code{vind}}{The order in which the variables were added to the submodel.}
#'  \item{\code{kl}}{KL-divergence for each submodel size.}
#'  \item{\code{summaries}}{Summary statistics computed during the selection.}
#'  \item{\code{d_test}}{The data used to evaluate the summaries.}
#'  \item{\code{family_kl}}{A modified \link{family}-object.}
#' }
#'
#' @examples
#' \donttest{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fit_v <- varsel(fit)
#' plot_varsel(fit_v)
#' }
#'

#' @export
varsel <- function(fit, d_test = NULL, method = NULL, ns = NULL, nc = NULL, 
                   nspred = NULL, ncpred = NULL, nv_max = NULL, 
                   intercept = NULL, penalty=NULL, verbose = F, 
                   lambda_min_ratio=1e-5, nlambda=500, regul=1e-6, ...) {


  .validate_for_varsel(fit)
	vars <- .extract_vars(fit)
	family_kl <- vars$fam
	
	if (is.null(method)) {
		if (dim(vars$x)[2] <= 20)
			method <- 'forward'
		else
			method <- 'L1'
	}

  if ((is.null(ns) && is.null(nc)) || tolower(method)=='l1')
  	# use one cluster for selection by default, and always with L1-search
  	nc <- 1
  if (is.null(nspred) && is.null(ncpred))
    # use 5 clusters for prediction by default
		ncpred <- min(ncol(vars$mu), 5)

  if(is.null(intercept))
    intercept <- vars$intercept
  if(is.null(nv_max) || nv_max > NCOL(vars$x)) {
  	nv_max_default <- floor(0.4*length(vars$y)) # a somewhat sensible default limit for nv_max
  	nv_max <- min(NCOL(vars$x), nv_max_default)
  }

  # training and test data
  d_train <- .get_traindata(fit)
  if (is.null(d_test)) {
  	d_test <- d_train
  	d_type <- 'train'
  } else {
  	d_test <- .check_data(d_test)
  	d_type <- 'test'
  }

  # reference distributions for selection and prediction after selection
  p_sel <- .get_refdist(fit, ns, nc)
  p_pred <- .get_refdist(fit, nspred, ncpred)

  # perform the selection
  opt <- list(lambda_min_ratio=lambda_min_ratio, nlambda=nlambda, regul=regul)
  vind <- select(method, p_sel, d_train, family_kl, intercept, nv_max, penalty, verbose, opt)

  # statistics for the selected submodels
  p_sub <- .get_submodels(vind, c(0, seq_along(vind)), family_kl, p_pred,
                          d_train, intercept, regul)
  sub <- .get_sub_summaries(p_sub, d_test, family_kl)

  # predictive statistics of the reference model on test data. if no test data are provided, 
  # simply fetch the statistics on the train data
  if ('datafit' %in% class(fit)) {
  	# no actual reference model, so we don't know how to predict test observations
  	full <- list(mu=rep(NA,vars$nobs), lppd=rep(NA,vars$nobs))
  } else {
  	if (d_type == 'train') {
  		full <- .weighted_summary_means(d_test, family_kl, vars$wsample, vars$mu, vars$dis)
  	} else {
  		mu_test <- vars$predfun(d_test$x, d_test$offset)
  		full <- .weighted_summary_means(d_test, family_kl, vars$wsample, mu_test, vars$dis)
  	}
  }
  
  # store the relevant fields into fit
  fit$varsel <- list(vind = setNames(vind, vars$coefnames[vind]),
                     kl = sapply(p_sub, function(x) x$kl),
                     d_test = c(d_test[c('y','weights')], type = d_type),
                     summaries = list(sub = sub, full = full),
                     family_kl = family_kl)

  # suggest model size
  ssize <- .suggest_size(fit$varsel)
  # if(is.na(ssize)) {
    # try a more relaxed value, if this does not work either, issue a warning
    # ssize <- .suggest_size(fit$varsel, cutoff_pct = 0.2)
    # if(is.na(ssize))
      # warning('Submodels too close to each other, cant suggest a submodel.')
  # }
  fit$varsel$ssize <- ssize
  
  fit
}


select <- function(method, p_sel, d_train, family_kl, intercept, nv_max,
                   penalty, verbose, opt) {
  #
  # Auxiliary function, performs variable selection with the given method,
  # and returns the variable ordering.
  #
  if (NCOL(d_train$x) == 1)
    # special case, only one variable, so no need for selection
    return(1)
  if (tolower(method) == 'l1') {
    vind <- search_L1(p_sel, d_train, family_kl, intercept, nv_max, penalty, opt)
  } else if (tolower(method) == 'forward') {
    if ( NCOL(p_sel$mu) == 1)
      # only one mu column (one cluster or one sample), so use the optimized version of the forward search
      vind <- search_forward1(p_sel, d_train, family_kl, intercept, nv_max, verbose, opt)
    else
      # routine that can be used with several clusters
      tryCatch(vind <- search_forward(p_sel, d_train, family_kl, intercept, nv_max, verbose, opt),
               'error' = .varsel_errors)
  } else {
    stop(sprintf('Unknown search method: %s.', method))
  }
  return(vind)
}

