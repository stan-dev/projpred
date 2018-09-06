#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for generalized linear models using
#' generic reference models.
#' 
#' @param object Either a \code{refmodel}-type object created by \link[=get_refmodel]{get_refmodel}
#' or \link[=init_refmodel]{init_refmodel}, or an object which can be converted to a reference model
#' using \link[=get_refmodel]{get_refmodel}.
#' @param d_test A test dataset, which is used to evaluate model performance.
#' If not provided, training data is used. Currently this argument is for internal use only.
#' @param method The method used in the variable selection. Possible options are
#' \code{'L1'} for L1-search and \code{'forward'} for forward selection.
#' Default is 'forward' if the number of variables in the full data is at most 20, and
#' \code{'L1'} otherwise.
#' @param relax If TRUE, then the projected coefficients after L1-selection are computed
#' without any penalization (or using only the regularization determined by \code{regul}). If FALSE, then
#' the coefficients are the solution from the L1-penalized projection. This option is relevant only
#' if \code{method}='L1'. Default is TRUE for genuine reference models and FALSE if \code{object} is
#' datafit (see \link[=init_refmodel]{init_refmodel}).  
#' @param ns Number of posterior draws used in the variable selection.
#'    Cannot be larger than the number of draws in the full model.
#'    Ignored if nc is set.
#' @param nc Number of clusters to use in the clustered projection.
#'    Overrides the \code{ns} argument. Defaults to 1.
#' @param nspred Number of samples used for prediction (after selection). Ignored if ncpred is given.
#' @param ncpred Number of clusters used for prediction (after selection). Default is 5.
#' @param nv_max Maximum number of varibles until which the selection is continued.
#'    Defaults to min(20, D, floor(0.4*n)) where n is the number of observations and
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
#' @param thresh Convergence threshold when computing L1-path. Usually no need to change this.
#' @param regul Amount of regularization in the projection. Usually there is no need for 
#' regularization, but sometimes for some models the projection can be ill-behaved and we
#' need to add some regularization to avoid numerical problems. Default is 1e-9.
#' @param ... Additional arguments to be passed to the \code{get_refmodel}-function.
#'
#'
#' @return An object of type \code{vsel} that contains information about the feature selection. The fields are not 
#' meant to be accessed directly by the user but instead via the helper functions (see the vignettes or type ?projpred
#' to see the main functions in the package.)
#'
#' @examples
#' \donttest{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' vs <- varsel(fit)
#' varsel_plot(vs)
#' }
#'

#' @export
varsel <- function(object, d_test = NULL, method = NULL, ns = NULL, nc = NULL, 
                   nspred = NULL, ncpred = NULL, relax=NULL, nv_max = NULL, 
                   intercept = NULL, penalty=NULL, verbose = F, 
                   lambda_min_ratio=1e-5, nlambda=150, thresh=1e-6, regul=1e-6, ...) {

	refmodel <- get_refmodel(object, ...)
	family_kl <- refmodel$fam
	
	if (is.null(method)) {
		if (dim(refmodel$x)[2] <= 20)
			method <- 'forward'
		else
			method <- 'L1'
	}

	if (is.null(relax)) {
	  if ('datafit' %in% class(refmodel))
	    relax <- F
	  else
	    relax <- T 
	}
	
  if ((is.null(ns) && is.null(nc)) || tolower(method)=='l1')
  	# use one cluster for selection by default, and always with L1-search
  	nc <- 1
  if (is.null(nspred) && is.null(ncpred))
    # use 5 clusters for prediction by default
		ncpred <- min(ncol(refmodel$mu), 5)

  if(is.null(intercept))
    intercept <- refmodel$intercept
  if(is.null(nv_max) || nv_max > NCOL(refmodel$x)) {
  	nv_max_default <- floor(0.4*length(refmodel$y)) # a somewhat sensible default limit for nv_max
  	nv_max <- min(NCOL(refmodel$x), nv_max_default, 20)
  }

  # training and test data
  d_train <- .get_traindata(refmodel)
  if (is.null(d_test)) {
  	d_test <- d_train
  	d_type <- 'train'
  } else {
  	d_test <- .check_data(d_test)
  	d_type <- 'test'
  }

  # reference distributions for selection and prediction after selection
  p_sel <- .get_refdist(refmodel, ns, nc)
  p_pred <- .get_refdist(refmodel, nspred, ncpred)

  # perform the selection
  opt <- list(lambda_min_ratio=lambda_min_ratio, nlambda=nlambda, thresh=thresh, regul=regul)
  searchpath <- select(method, p_sel, d_train, family_kl, intercept, nv_max, penalty, verbose, opt)
  vind <- searchpath$vind
  
  # statistics for the selected submodels
  as.search <- !relax && !is.null(searchpath$beta) && !is.null(searchpath$alpha)
  p_sub <- .get_submodels(searchpath, c(0, seq_along(vind)), family_kl, p_pred,
                          d_train, intercept, regul, as.search=as.search)
  sub <- .get_sub_summaries(p_sub, d_test, family_kl)

  # predictive statistics of the reference model on test data. if no test data are provided, 
  # simply fetch the statistics on the train data
  if ('datafit' %in% class(refmodel)) {
  	# no actual reference model, so we don't know how to predict test observations
    ntest <- nrow(d_test$z)
  	full <- list(mu=rep(NA,ntest), lppd=rep(NA,ntest))
  } else {
  	if (d_type == 'train') {
  		full <- .weighted_summary_means(d_test, family_kl, refmodel$wsample, refmodel$mu, refmodel$dis)
  	} else {
  		mu_test <- refmodel$predfun(d_test$z, d_test$offset)
  		full <- .weighted_summary_means(d_test, family_kl, refmodel$wsample, mu_test, refmodel$dis)
  	}
  }
  
  # store the relevant fields into the object to be returned
  vs <- list(refmodel=refmodel,
  					 spath=searchpath,
             d_test = c(d_test[c('y','weights')], type = d_type),
             summaries = list(sub = sub, full = full),
             family_kl = family_kl,
  					 vind = setNames(vind, refmodel$coefnames[vind]),
  					 kl = sapply(p_sub, function(x) x$kl) )
  class(vs) <- 'vsel'

  # suggest model size
  vs$nv_max <- nv_max
  vs$nv_all <- ncol(refmodel$x)
  vs$ssize <- suggest_size(vs, warnings = F)
  
  vs
}


select <- function(method, p_sel, d_train, family_kl, intercept, nv_max,
                   penalty, verbose, opt) {
  #
  # Auxiliary function, performs variable selection with the given method,
  # and returns the searchpath, i.e., a list with the followint entries (the last three
  # are returned only if one cluster projection is used for selection):
  #   vind: the variable ordering
  #   beta: coefficients along the search path 
  #   alpha: intercepts along the search path 
  #   p_sel: the reference distribution used in the selection (the input argument p_sel)
  #
  if (tolower(method) == 'l1') {
    searchpath <- search_L1(p_sel, d_train, family_kl, intercept, nv_max, penalty, opt)
    searchpath$p_sel <- p_sel
    return(searchpath)
  } else if (tolower(method) == 'forward') {
    if ( NCOL(p_sel$mu) == 1) {
      # only one mu column (one cluster or one sample), so use the optimized version of the forward search
      searchpath <- search_forward1(p_sel, d_train, family_kl, intercept, nv_max, verbose, opt)
      searchpath$p_sel <- p_sel
      return(searchpath)
    } else {
      # routine that can be used with several clusters
      tryCatch(vind <- search_forward(p_sel, d_train, family_kl, intercept, nv_max, verbose, opt),
               'error' = .varsel_errors)
      searchpath <- list(vind=vind, p_sel=p_sel)
      return(searchpath)
    }
  } else {
    stop(sprintf('Unknown search method: %s.', method))
  }
}




# parse_varsel_args <- function(n, d, method = NULL, cv_method = NULL, 
#                               ns = NULL, nc = NULL, nspred = NULL, ncpred = NULL, relax = NULL,
#                               nv_max = NULL, intercept = NULL, penalty = NULL, verbose = NULL,
#                               nloo = NULL, K = NULL, k_fold = NULL, lambda_min_ratio = NULL, 
#                               nlambda = NULL, regul = NULL, validate_search = NULL, seed = NULL, ...) {
#   #
#   # Auxiliary function for figuring out the parameters for varsel and cv_varsel. The arguments
#   # specified by the user (or the function calling this function) are treated as they are, but if 
#   # some are not given, then this function fills them in with the default values (by default, use
#   # same values for both varsel and cv_varsel). The purpose of this function is to avoid repeating
#   # the same (longish) code both in varsel and cv_varsel.
#   #
#   if (is.null(seed))
#     seed <- 134654
#   
#   if (is.null(method)) {
#     if (dim(vars$x)[2] <= 20)
#       method <- 'forward'
#     else
#       method <- 'L1'
#   }
#   
#   if (is.null(relax)) {
#     if ('datafit' %in% class(refmodel))
#       relax <- F
#     else
#       relax <- T 
#   }
#   
#   if (is.null(cv_method)) {
#     if ('datafit' %in% class(refmodel))
#       # only data given, no actual reference model
#       cv_method <- 'kfold'
#     else
#       cv_method <- 'LOO'
#   }
#   if (cv_method == 'kfold' && is.null(K)) {
#     if ('datafit' %in% class(refmodel))
#       K <- 10
#     else 
#       K <- 4
#   }
#   
#   if ((is.null(ns) && is.null(nc)) || tolower(method)=='l1')
#     # use one cluster for selection by default, and always with L1-search
#     nc <- 1
#   if (is.null(nspred) && is.null(ncpred))
#     # use 5 clusters for prediction by default
#     ncpred <- min(ncol(vars$mu), 5)
#   
#   if (is.null(intercept))
#     intercept <- vars$intercept
#   if (is.null(nv_max) || nv_max > NCOL(vars$x)) {
#     nv_max_default <- floor(0.4*length(vars$y)) # a somewhat sensible default limit for nv_max
#     nv_max <- min(NCOL(vars$x), nv_max_default, 20)
#   }
#   
#   args <- list(method=method, cv_method=cv_method, ns=ns, nc=nc, nspred=nspred, ncpred=ncpred, 
#                relax=relax, nv_max=nv_max, intercept=intercept, penalty=penalty, verbose=verbose,
#                nloo=nloo, K=K, k_fold=k_fold, lambda_min_ratio=lambda_min_ratio, nlambda=nlambda, 
#                regul=regul, validate_search=validate_search, seed=seed)
#   
# }



