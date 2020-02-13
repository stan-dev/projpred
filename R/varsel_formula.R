varsel_poc <- function(object, d_test = NULL, method = NULL, ns = NULL, nc = NULL,
                       nspred = NULL, ncpred = NULL, relax=FALSE, nv_max = NULL,
                       intercept = NULL, penalty=NULL, verbose = F,
                       lambda_min_ratio=1e-5, nlambda=150, thresh=1e-6, regul=1e-4,
                       groups=NULL, ...) {

  refmodel <- get_refmodel_poc(object, ...)
	family_kl <- refmodel$fam

	# fetch the default arguments or replace them by the user defined values
	args <- parseargs_varsel_poc(refmodel, method, relax, intercept, nv_max, nc, ns, ncpred, nspred, groups)
	method <- args$method
	relax <- args$relax
	intercept <- args$intercept
	nv_max <- args$nv_max
	nc <- args$nc
	ns <- args$ns
	ncpred <- args$ncpred
	nspred <- args$nspred
  groups <- args$groups
  group_features <- formula_contains_group_terms(refmodel$formula)

  if (tolower(method) == 'l1' && group_features) {
    warning('l1 search is not supported for multilevel models, switching to forward search')
    method <- "forward"
  }

  if (is.null(d_test)) {
    d_type <- 'train'
    d_test <- list(y=refmodel$y, test_points=seq_along(refmodel$y))
  }

  # reference distributions for selection and prediction after selection
  p_sel <- .get_refdist(refmodel, ns, nc)
  p_pred <- .get_refdist(refmodel, nspred, ncpred)

  # perform the selection
  opt <- list(lambda_min_ratio=lambda_min_ratio, nlambda=nlambda, thresh=thresh, regul=regul)
  searchpath <- select_poc(method, p_sel, refmodel, family_kl, intercept, nv_max,
                           penalty, verbose, opt, groups=groups)
  vind <- searchpath$vind

  # statistics for the selected submodels
  as.search <- relax
  p_sub <- .get_submodels_poc(searchpath, c(0, seq_along(vind)), family_kl, p_pred,
                              refmodel, intercept, regul, as.search=as.search)
  sub <- .get_sub_summaries_poc(p_sub, seq_along(refmodel$y), refmodel, family_kl)

  # predictive statistics of the reference model on test data. if no test data are provided,
  # simply fetch the statistics on the train data
  if ('datafit' %in% class(refmodel)) {
  	# no actual reference model, so we don't know how to predict test observations
    ntest <- nrow(refmodel$y)
  	ref <- list(mu=rep(NA,ntest), lppd=rep(NA,ntest))
  } else {
    d_test$weights <- refmodel$wobs[d_test$test_points]
  	if (d_type == 'train') {
  		ref <- .weighted_summary_means_poc(d_test, family_kl, refmodel$wsample, refmodel$mu, refmodel$dis)
  	} else {
  		mu_test <- refmodel$predfun(refmodel$fit, newdata=d_test$data)
  		ref <- .weighted_summary_means_poc(d_test, family_kl, refmodel$wsample, mu_test, refmodel$dis)
  	}
  }

  # store the relevant fields into the object to be returned
  vs <- list(refmodel=refmodel,
  					 spath=searchpath,
             d_test = c(d_test["y"], type = d_type),
             summaries = list(sub = sub, ref = ref),
             family_kl = family_kl,
  					 vind = searchpath$vind,
  					 kl = sapply(p_sub, function(x) x$kl) )

  class(vs) <- 'vsel'
  # suggest model size
  vs$nv_max <- nv_max
  vs$nv_all <- count_terms_in_submodel(refmodel$formula)
  vs$ssize <- suggest_size(vs, warnings = F, group_features = group_features, groups = groups)

  vs
}


select_poc <- function(method, p_sel, refmodel, family_kl, intercept, nv_max,
                       penalty, verbose, opt, groups=NULL) {
  ##
  ## Auxiliary function, performs variable selection with the given method,
  ## and returns the searchpath, i.e., a list with the followint entries (the last three
  ## are returned only if one cluster projection is used for selection):
  ##   vind: the variable ordering
  ##   beta: coefficients along the search path
  ##   alpha: intercepts along the search path
  ##   p_sel: the reference distribution used in the selection (the input argument p_sel)
  ##
  ## routine that can be used with several clusters
  if (tolower(method) == 'l1') {
    searchpath <- search_L1_poc(p_sel, refmodel, family_kl, intercept, nv_max, penalty, opt)
    searchpath$p_sel <- p_sel
    return(searchpath)
  } else if (tolower(method) == 'forward') {
    tryCatch(searchpath <- search_forward_poc(p_sel, refmodel, family_kl, intercept, nv_max, verbose, opt, groups=groups),
             'error' = .varsel_errors)

    searchpath$p_sel <- p_sel
    return(searchpath)
  }
}


parseargs_varsel_poc <- function(refmodel, method, relax, intercept, nv_max, nc, ns, ncpred, nspred, groups) {

  #
  # Auxiliary function for parsing the input arguments for varsel. The arguments
  # specified by the user (or the function calling this function) are treated as they are, but if
  # some are not given, then this function fills them in with the default values. The purpose of this
  # function is to avoid repeating the same code both in varsel and cv_varsel.
  #

  if (is.null(groups))
    groups <- lapply(break_formula(refmodel$formula), function(t) t)

  group_features <- formula_contains_group_terms(refmodel$formula)

  if (is.null(method)) {
    if (group_features)
      if (length(groups) <= 20)
        method <- 'forward'
      else
        method <- 'l1'
    else
      ## if we are doing a grouped search we don't usually have that many groups
      method <- 'forward'
  } else
    if (group_features)
      method <- 'forward'

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

  max_nv_possible <- count_terms_in_submodel(refmodel$formula)
  if(is.null(intercept))
    intercept <- refmodel$intercept
  if(is.null(nv_max))
    nv_max <- min(max_nv_possible, 20)
  else
    nv_max <- min(max_nv_possible, nv_max)


  list(method=method, relax=relax, intercept=intercept, nv_max=nv_max,
       nc=nc, ns=ns, ncpred=ncpred, nspred=nspred, groups=groups)
}
