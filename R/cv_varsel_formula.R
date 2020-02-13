cv_varsel_poc <- function(fit,  method = NULL, cv_method = NULL,
                          ns = NULL, nc = NULL, nspred = NULL, ncpred = NULL, relax=FALSE,
                          nv_max = NULL, intercept = NULL, penalty = NULL, verbose = T,
                          nloo=NULL, K = NULL, lambda_min_ratio=1e-5, nlambda=150,
                          thresh=1e-6, regul=1e-4, validate_search=T, seed=NULL,
                          groups=NULL, ...) {

  refmodel <- get_refmodel_poc(fit, ...)

	# resolve the arguments similar to varsel
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
  group_features <- !is.null(groups)

	# arguments specific to this function
	args <- parseargs_cv_varsel(refmodel, cv_method, K)
	cv_method <- args$cv_method
	K <- args$K

	# search options
	opt <- list(lambda_min_ratio=lambda_min_ratio, nlambda=nlambda, thresh=thresh, regul=regul)

  if (tolower(cv_method) == 'loo')  {
    if (!(is.null(K))) warning('K provided, but cv_method is LOO.')
    sel_cv <- loo_varsel_poc(refmodel, method, nv_max, ns, nc, nspred, ncpred, relax, intercept, penalty,
                             verbose, opt, nloo = nloo, validate_search = validate_search, seed = seed,
                             groups=groups)
  } else if (tolower(cv_method) == 'kfold')  {
    sel_cv <- kfold_varsel_poc(refmodel, method, nv_max, ns, nc, nspred, ncpred, relax, intercept, penalty,
                               verbose, opt, K, seed = seed, groups=groups)
  } else {
    stop(sprintf('Unknown cross-validation method: %s.', method))
  }

  ## run the selection using the full dataset
  if (verbose)
    print(paste('Performing the selection using all the data..'))
  sel <- varsel_poc(refmodel, method=method, ns=ns, nc=nc, nspred=nspred, ncpred=ncpred,
                    relax=relax, nv_max=nv_max, intercept=intercept, penalty=penalty, verbose=verbose,
                    lambda_min_ratio=lambda_min_ratio, nlambda=nlambda, regul=regul,
                    groups=groups)


	# find out how many of cross-validated iterations select
	# the same variables as the selection with all the data.
	ch <- as.matrix(unname(as.data.frame(sel_cv$vind_cv)))
  w <- sel_cv$summaries$sub[[1]]$w # these weights might be non-constant in case of subsampling LOO
  selvind <- sel$vind
  if (is.null(w)) w <- rep(1/ncol(ch), ncol(ch)) # if weights are not set, then all validation folds have equal weight
  vars <- unlist(selvind)
  pctch <- t(sapply(seq_along(selvind), function(size) {
    c(size = size, sapply(vars, function(var) {
      sum(t(ch[1:size, ] == var) * w, na.rm = T)
    }))
  }))
  colnames(pctch)[-1] <- vars

  # create the object to be returned
	vs <- list()
	vs$refmodel <- refmodel
	vs$spath <- sel$spath
	vs <- c(vs, c(sel_cv[c('d_test', 'summaries')],
	              sel[c('family_kl', 'vind', 'kl')],
	              list(pctch = pctch)))
  class(vs) <- 'cvsel'
	vs$nv_max <- nv_max
	vs$nv_all <- count_terms_in_submodel(refmodel$formula)
	vs$ssize <- suggest_size(vs, warnings = F, group_features = group_features, groups = groups)

	if (verbose)
	  print('Done.')

	vs
}

parseargs_cv_varsel <- function(refmodel, cv_method, K) {

  ##
  ## Auxiliary function for parsing the input arguments for specific cv_varsel.
  ## This is similar in spirit to parseargs_varsel, that is, to avoid the main function to become
  ## too long and complicated to maintain.
  ##

  if (is.null(cv_method)) {
    if ('datafit' %in% class(refmodel))
      ## only data given, no actual reference model
      cv_method <- 'kfold'
    else
      cv_method <- 'LOO'
  }
  if (cv_method == 'kfold' && is.null(K)) {
    if ('datafit' %in% class(refmodel))
      K <- 10
    else
      K <- 5
  }

  list(cv_method=cv_method, K=K)
}

loo_varsel_poc <- function(refmodel, method, nv_max, ns, nc, nspred, ncpred, relax, intercept,
                       penalty, verbose, opt, nloo = NULL, validate_search = T, seed = NULL,
                       groups=NULL) {
  ##
  ## Performs the validation of the searching process using LOO.
  ## validate_search indicates whether the selection is performed separately for each
  ## fold (for each data point)
  ##

  group_features <- !is.null(groups)
  fam <- refmodel$fam
  mu <- refmodel$mu
  dis <- refmodel$dis
  # the clustering/subsampling used for selection
	p_sel <- .get_refdist(refmodel, ns=ns, nc=nc)
	cl_sel <- p_sel$cl # clustering information

	# the clustering/subsampling used for prediction
	p_pred <- .get_refdist(refmodel, ns=nspred, nc=ncpred)
	cl_pred <- p_pred$cl

	# fetch the log-likelihood for the reference model to obtain the LOO weights
	if (is.null(refmodel$loglik))
		# case where log-likelihood not available, i.e., the reference model is not a genuine model
		# => cannot compute LOO
		stop('LOO can be performed only if the reference model is a genuine probabilistic model for
          which the log-likelihood can be evaluated.')
	else
		# log-likelihood available
		loglik <- refmodel$loglik
	psisloo <- loo::psis(-loglik, cores = 1, r_eff = rep(1,ncol(loglik))) # TODO: should take r_eff:s into account
	lw <- weights(psisloo)
	pareto_k <- loo::pareto_k_values(psisloo)
	n <- length(pareto_k)
	nloo <- ifelse(is.null(nloo), n, nloo) # by default use all observations
	nloo <- min(nloo,n)

	# compute loo summaries for the reference model
	loo_ref <- apply(loglik+lw, 2, 'log_sum_exp')
	mu_ref <- rep(0,n)
	for (i in 1:n)
    mu_ref[i] <- mu[i,] %*% exp(lw[,i])

	# decide which points form the validation set based on the k-values
	validset <- .loo_subsample(n, nloo, pareto_k, seed)
	inds <- validset$inds

	# initialize matrices where to store the results
	vind_mat <- matrix(nrow=n, ncol=nv_max)
	loo_sub <- matrix(nrow=n, ncol=nv_max+1)
	mu_sub <- matrix(nrow=n, ncol=nv_max+1)

	if (verbose) {
    print('Computing LOOs...')
    pb <- utils::txtProgressBar(min = 0, max = nloo, style = 3, initial=0)
	}

	if (!validate_search) {
	  # perform selection only once using all the data (not separately for each fold),
	  # and perform the projection then for each submodel size
	  # vind <- select(method, p_sel, refmodel, fam, intercept, nv_max, penalty, verbose=F, opt)
	  spath <- select_poc(method, p_sel, refmodel, fam, intercept, nv_max, penalty,
                        verbose=F, opt, groups=groups)
    vind <- spath$vind
  }

	for (run_index in seq_along(inds)) {

	  # observation index
	  i <- inds[run_index]

	  # reweight the clusters/samples according to the is-loo weights
	  p_sel <- .get_p_clust(fam, mu, dis, wsample=exp(lw[,i]), cl=cl_sel)
	  p_pred <- .get_p_clust(fam, mu, dis, wsample=exp(lw[,i]), cl=cl_pred)

		if (validate_search) {
		  # perform selection with the reweighted clusters/samples
		  # vind <- select(method, p_sel, refmodel, fam, intercept, nv_max, penalty, verbose=F, opt)
		  spath <- select_poc(method, p_sel, refmodel, fam, intercept, nv_max, penalty,
                          verbose=F, opt, groups=groups)
      vind <- spath$vind
    }

		# project onto the selected models and compute the prediction accuracy for the left-out point
	  as.search <- relax
    submodels <- .get_submodels_poc(spath, c(0, seq_along(vind)), fam, p_pred,
                                    refmodel, intercept, opt$regul, as.search=as.search)
		summaries_sub <- .get_sub_summaries_poc(submodels, c(i), refmodel, fam)

    for (k in seq_along(summaries_sub)) {
      loo_sub[i,k] <- summaries_sub[[k]]$lppd
      mu_sub[i,k] <- summaries_sub[[k]]$mu
    }

    ## we are always doing group selection
    ## with `match` we get the indices of the variables as they enter the
    ## solution path in vind
    vind_mat[i,] <- match(vind, groups)

    if (verbose) {
		  utils::setTxtProgressBar(pb, run_index)
		}
	}

	if (verbose)
	  # close the progress bar object
	  close(pb)

	# put all the results together in the form required by cv_varsel
	summ_sub <-	lapply(0:nv_max, function(k){
	    list(lppd=loo_sub[,k+1], mu=mu_sub[,k+1], w=validset$w)
	})
	summ_ref <- list(lppd=loo_ref, mu=mu_ref)
	summaries <- list(sub=summ_sub, ref=summ_ref)

  vind_cv <- lapply(1:n, function(i){ vind_mat[i,] })

  d_test <- list(y=refmodel$y, type='loo',
                 test_points=seq_along(refmodel$y))

	return(list(vind_cv=vind_cv, summaries=summaries, d_test=d_test))

}

kfold_varsel_poc <- function(refmodel, method, nv_max, ns, nc, nspred, ncpred, relax,
                             intercept, penalty, verbose, opt, K, seed=NULL, groups=NULL) {

  group_features <- !is.null(groups)
	# fetch the k_fold list (or compute it now if not already computed)
	k_fold <- .get_kfold_poc(refmodel, K, verbose, seed)

	# check that k_fold has the correct form
	# .validate_kfold(refmodel, k_fold, refmodel$nobs)

	K <- length(k_fold$refmodel)
	family_kl <- refmodel$fam

  # extract variables from each fit-object (samples, x, y, etc.)
  # to a list of size K
	refmodels_cv <- k_fold$refmodel

  # List of K elements, each containing d_train, p_pred, etc. corresponding
  # to each fold.
  msgs <- paste0(method, ' search for fold ', 1:K, '/', K, '.')
  list_cv <- mapply(function(refmod, test_points, msg) {
  	p_sel <- .get_refdist(refmod, ns, nc)
  	p_pred <- .get_refdist(refmod, nspred, ncpred)
    newdata <- refmod$fetch_data(data_points=test_points)
  	mu_test <- refmod$predfun(refmod$fit, newdata=newdata)
  	list(refmodel = refmod, p_sel = p_sel, p_pred = p_pred,
         mu_test = mu_test, dis = refmod$dis, w_test = refmod$wsample,
         msg = msg, y_test = refmod$y[test_points], test_points = test_points)
  }, refmodels_cv, k_fold$test_points, msgs, SIMPLIFY = F)

  ## Perform the selection for each of the K folds
  if (verbose) {
    print('Performing selection for each fold..')
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial=0)
  }
  spath_cv <- lapply(seq_along(list_cv), function(fold_index) {
    fold <- list_cv[[fold_index]]
    family_kl <- fold$refmodel$family
    out <- select_poc(method, fold$p_sel, fold$refmodel, family_kl, intercept,
                      nv_max, penalty, verbose, opt, groups=groups)
    if (verbose)
      utils::setTxtProgressBar(pb, fold_index)
    out
  })

  vind_cv <- lapply(spath_cv, function(e) e$vind)
  if (verbose)
    close(pb)

  ## Construct submodel projections for each fold
  as.search <- relax
  if (verbose && !as.search) {
    print('Computing projections..')
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial=0)
  }
  p_sub_cv <- mapply(function(spath, fold_index) {
    fold <- list_cv[[fold_index]]
    family_kl <- fold$refmodel$family
    vind <- spath$vind
    p_sub <- .get_submodels_poc(spath, c(0, seq_along(vind)), family_kl,
                                fold$p_pred, fold$refmodel, intercept,
                                opt$regul, as.search=as.search)
    if (verbose && !as.search)
      utils::setTxtProgressBar(pb, fold_index)
    return(p_sub)
  }, spath_cv, seq_along(list_cv), SIMPLIFY = F)
  if (verbose && !as.search)
    close(pb)


  # Helper function extract and combine mu and lppd from K lists with each
  # n/K of the elements to one list with n elements
  hf <- function(x) as.list(do.call(rbind, x))

  ## Apply some magic to manipulate the structure of the list so that instead of
  ## list with K sub_summaries each containing n/K mu:s and lppd:s, we have only
  ## one sub_summary-list that contains with all n mu:s and lppd:s.
  sub <- apply(
    mapply(function(p_sub, fold) {
      family_kl <- fold$refmodel$family
      lapply(.get_sub_summaries_poc(p_sub, fold$test_points, fold$refmodel, family_kl),
             data.frame)
    }, p_sub_cv, list_cv),
    1, hf)

  ref <- hf(lapply(list_cv, function(fold) {
    family_kl <- fold$refmodel$family
    test <- fold$test_points
    d <- list(y=fold$y_test)
    data.frame(.weighted_summary_means_poc(d, family_kl, fold$mu_test, fold$dis))}))

  # Combine also the K separate test data sets into one list
  # with n y's and weights's.
  d_cv <- hf(lapply(list_cv, function(fold) {
    test <- fold$test_points
    data.frame(y=fold$y_test)}))

  list(vind_cv = vind_cv,
       summaries = list(sub = sub, ref = ref),
       d_test = c(d_cv, type = 'kfold'))
}


.get_kfold_poc <- function(refmodel, K, verbose, seed) {
  # Fetch the k_fold list or compute it now if not already computed. This function will
  # return a list of length K, where each element is a list with fields 'refmodel' (object
  # of type refmodel computed by init_refmodel) and index list 'test_points' that denotes which
  # of the data points were left out for the corresponding fold.

  ## TODO: pass cvfits for refmodel so that each fold has a different fit with
  ## the proper subset of the data
  indices <- seq_along(refmodel$y)
  test_points <- split(indices, sort(rank(indices) %% K) + 1)
  train_points <- lapply(test_points, function(test) setdiff(indices, test))

  if (is.null(refmodel$cvfits))
    stop("You need to provide cvfits!")

  refmodels <- lapply(1:K, function(i) {
    test <- test_points[[i]]
    train <- train_points[[i]]
    default_data <- refmodel$fetch_data(data_points=train)
    predfun <- function(fit, newdata=default_data)
      refmodel$predfun(fit, newdata=newdata)
    proj_predfun <- function(fit, newdata=default_data)
      refmodel$proj_predfun(fit, newdata=newdata)
    refmod <- get_refmodel_poc(refmodel$cvfits[[i]],
                               fetch_data(data=refmodel$fetch_data()),
                               refmodel$y, refmodel$formula, predfun,
                               proj_predfun, refmodel$mle, fetch_data,
                               family=refmodel$family,
                               folds=train)
    refmod
  })

  k_fold <- list(refmodel=refmodels, test_points=test_points, train_points=train_points)
  return(k_fold)
}
