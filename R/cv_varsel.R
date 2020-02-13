#' Cross-validate the variable selection (varsel)
#'
#' Perform cross-validation for the projective variable selection for a generalized
#' linear model.
#' @param fit Same as in \link[=varsel]{varsel}.
#' @param method Same as in \link[=varsel]{varsel}.
#' @param ns Number of samples used for selection. Ignored if nc is provided or if method='L1'.
#' @param nc Number of clusters used for selection. Default is 1 and ignored if method='L1' 
#' (L1-search uses always one cluster).
#' @param nspred Number of samples used for prediction (after selection). Ignored if ncpred is given.
#' @param ncpred Number of clusters used for prediction (after selection). Default is 5.
#' @param relax Same as in \link[=varsel]{varsel}.
#' @param nv_max Same as in \link[=varsel]{varsel}.
#' @param intercept Same as in \link[=varsel]{varsel}.
#' @param penalty Same as in \link[=varsel]{varsel}.
#' @param verbose Whether to print out some information during the validation, Default is TRUE.
#' @param cv_method The cross-validation method, either 'LOO' or 'kfold'. Default is 'LOO'.
#' @param nloo Number of observations used to compute the LOO validation (anything between 1 and the 
#' total number of observations). Smaller values lead to
#' faster computation but higher uncertainty (larger errorbars) in the accuracy estimation.
#' Default is to use all observations, but for faster experimentation, one can set this to a small value such as 100.
#' Only applicable if \code{cv_method = 'LOO'}. 
#' @param K Number of folds in the k-fold cross validation. Default is 5 for genuine
#' reference models and 10 for datafits (that is, for penalized maximum likelihood estimation).
#' @param lambda_min_ratio Same as in \link[=varsel]{varsel}.
#' @param nlambda Same as in \link[=varsel]{varsel}.
#' @param thresh Same as in \link[=varsel]{varsel}.
#' @param regul Amount of regularization in the projection. Usually there is no need for 
#' regularization, but sometimes for some models the projection can be ill-behaved and we
#' need to add some regularization to avoid numerical problems. 
#' @param validate_search Whether to cross-validate also the selection process, that is, whether to perform
#' selection separately for each fold. Default is TRUE and we strongly recommend not setting this
#' to FALSE, because this is known to bias the accuracy estimates for the selected submodels.
#' However, setting this to FALSE can sometimes be useful because comparing the results to the case
#' where this parameter is TRUE gives idea how strongly the feature selection is (over)fitted to the
#' data (the difference corresponds to the search degrees of freedom or the effective number 
#' of parameters introduced by the selectin process).
#' @param B number of bootstrap samples to generate from the Dirichlet distribution to compute pseudo-BMA+.
#' Only applies to loo.
#' @param seed Random seed used in the subsampling LOO. By default uses a fixed seed.
#' @param ... Additional arguments to be passed to the \code{get_refmodel}-function.
#'
#' @return An object of type \code{cvsel} that contains information about the feature selection. The fields are not 
#' meant to be accessed directly by the user but instead via the helper functions (see the vignettes or type ?projpred
#' to see the main functions in the package.)
#'
#' @examples
#' \donttest{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' cvs <- cv_varsel(fit)
#' varsel_plot(cvs)
#' }
#'

#' @export
cv_varsel <- function(fit,  method = NULL, cv_method = NULL, 
                      ns = NULL, nc = NULL, nspred = NULL, ncpred = NULL, relax=NULL,
                      nv_max = NULL, intercept = NULL, penalty = NULL, verbose = T,
                      nloo=NULL, K = NULL, lambda_min_ratio=1e-5, nlambda=150,
                      thresh=1e-6, regul=1e-4, validate_search=T, seed=NULL,
                      B = 20, ...) {

	refmodel <- get_refmodel(fit, ...)
	
	# resolve the arguments similar to varsel
	args <- parseargs_varsel(refmodel, method, relax, intercept, nv_max, nc, ns, ncpred, nspred)
	method <- args$method
	relax <- args$relax
	intercept <- args$intercept
	nv_max <- args$nv_max
	nc <- args$nc
	ns <- args$ns
	ncpred <- args$ncpred
	nspred <- args$nspred
	
	# arguments specific to this function
	args <- parseargs_cv_varsel(refmodel, cv_method, K)
	cv_method <- args$cv_method
	K <- args$K

	# search options
	opt <- list(lambda_min_ratio=lambda_min_ratio, nlambda=nlambda, thresh=thresh, regul=regul)
	
	if (tolower(cv_method) == 'kfold') {
	  # TODO: should we save the cvfits object to the reference model so that it need not be computed again
	  # if the user wants to compute the search again?
		sel_cv <- kfold_varsel(refmodel, method, nv_max, ns, nc, nspred, ncpred, relax, intercept, penalty,
		                       verbose, opt, K, seed=seed)
	} else if (tolower(cv_method) == 'loo')  {
	  if (!(is.null(K))) warning('K provided, but cv_method is LOO.')
		sel_cv <- loo_varsel(refmodel, method, nv_max, ns, nc, nspred, ncpred, relax, intercept, penalty, 
		                     verbose, opt, nloo = nloo, validate_search = validate_search, B = B, seed = seed)
	} else {
               stop(sprintf('Unknown cross-validation method: %s.', cv_method))
	}
	
	# run the selection using the full dataset
	if (verbose)
		cat('Performing variable selection using all data...\n')
	sel <- varsel(refmodel, method=method, ns=ns, nc=nc, nspred=nspred, ncpred=ncpred,
	              relax=relax, nv_max=nv_max, intercept=intercept, penalty=penalty, verbose=verbose, 
	              lambda_min_ratio=lambda_min_ratio, nlambda=nlambda, regul=regul)


	# find out how many of cross-validated iterations select
	# the same variables as the selection with all the data.
	ch <- as.matrix(unname(as.data.frame(sel_cv$vind_cv)))
	w <- sel_cv$summaries$sub[[1]]$w # these weights might be non-constant in case of subsampling LOO
	if (is.null(w)) w <- rep(1/ncol(ch), ncol(ch)) # if weights are not set, then all validation folds have equal weight
	pctch <- t(sapply(seq_along(sel$vind), function(size) {
	  c(size = size, sapply(sel$vind, function(var) {
	    sum(t(ch[1:size, ] == var) * w, na.rm = T)
	  }))
	}))
	colnames(pctch)[-1] <- names(sel$vind)

	
	# create the object to be returned
	vs <- list()
	vs$refmodel <- refmodel
	vs$spath <- sel$spath
	vs$method <- method
	vs$cv_method <- cv_method
	vs <- c(vs, c(sel_cv[c('d_test', 'summaries', 'pseudo_bma')],
	              sel[c('family_kl', 'vind', 'kl')],
	              list(pctch = pctch)))
	class(vs) <- 'cvsel'
	vs$nv_max <- nv_max
	vs$nv_all <- ncol(refmodel$x)
	vs$ssize <- suggest_size(vs, warnings = F)

	vs
}

parseargs_cv_varsel <- function(refmodel, cv_method, K) {
  
  #
  # Auxiliary function for parsing the input arguments for specific cv_varsel.
  # This is similar in spirit to parseargs_varsel, that is, to avoid the main function to become
  # too long and complicated to maintain.
  #
  
  if (is.null(cv_method)) {
    if ('datafit' %in% class(refmodel))
      # only data given, no actual reference model
      cv_method <- 'kfold'
    else
      cv_method <- 'LOO'
  }
  if (cv_method == 'kfold') {
    if (is.null(K))
      K <- if ('datafit' %in% class(refmodel)) 10 else 5
    else
      .validate_num_folds(K, refmodel$nobs)
  }
  
  list(cv_method=cv_method, K=K)
}


kfold_varsel <- function(refmodel, method, nv_max, ns, nc, nspred, ncpred, relax,
                         intercept, penalty, verbose, opt, K, seed=NULL) {
	
	# fetch the k_fold list (or compute it now if not already computed)
	k_fold <- .get_kfold(refmodel, K, verbose, seed)
	
	# check that k_fold has the correct form
	# .validate_kfold(refmodel, k_fold, refmodel$nobs)
	
	K <- length(k_fold) 
	family_kl <- refmodel$fam
	
  # extract variables from each fit-object (samples, x, y, etc.)
  # to a list of size K
	refmodels_cv <- lapply(k_fold, function(fold) fold$refmodel)

  # List of size K with test data for each fold
	d_test_cv <- lapply(k_fold, function(fold) {
	  list(z = refmodel$z[fold$omitted,,drop=F], 
	  		 x = refmodel$x[fold$omitted,,drop=F], y = refmodel$y[fold$omitted],
	       weights = refmodel$wobs[fold$omitted], offset = refmodel$offset[fold$omitted])
	})

  # List of K elements, each containing d_train, p_pred, etc. corresponding
  # to each fold.
  msgs <- paste0(method, ' search for fold ', 1:K, '/', K, '.')
  list_cv <- mapply(function(refmod, d_test, msg) {
  	d_train <- .get_traindata(refmod)
  	p_sel <- .get_refdist(refmod, ns, nc)
  	p_pred <- .get_refdist(refmod, nspred, ncpred)
  	mu_test <- refmod$predfun(d_test$z, d_test$offset)
  	list(d_train = d_train, d_test = d_test, p_sel = p_sel, p_pred = p_pred,
  	     mu_test = mu_test, dis = refmod$dis, w_test = refmod$wsample, msg = msg)
  }, refmodels_cv, d_test_cv, msgs, SIMPLIFY = F)
  
  # Perform the selection for each of the K folds
  if (verbose) {
    cat('Performing selection for each fold...\n')
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial=0)
  }
  spath_cv <- lapply(seq_along(list_cv), function(fold_index) {
    fold <- list_cv[[fold_index]]
    out <- select(method, fold$p_sel, fold$d_train, family_kl, intercept, nv_max, penalty, verbose, opt)
    if (verbose)
      utils::setTxtProgressBar(pb, fold_index)
    out
  })
  vind_cv <- lapply(spath_cv, function(e) e$vind)
  if (verbose)
    close(pb)
  

  # Construct submodel projections for each fold
  as.search <- !relax && !is.null(spath_cv[[1]]$beta) && !is.null(spath_cv[[1]]$alpha)
  if (verbose && !as.search) {
    cat('Computing projections...\n')
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial=0)
  }
  p_sub_cv <- mapply(function(spath, fold_index) {
    fold <- list_cv[[fold_index]]
    vind <- spath$vind
    p_sub <- .get_submodels(spath, c(0, seq_along(vind)), family_kl, fold$p_pred,
                            fold$d_train, intercept, opt$regul, as.search=as.search)
    if (verbose && !as.search)
      utils::setTxtProgressBar(pb, fold_index)
    return(p_sub)
  }, spath_cv, seq_along(list_cv), SIMPLIFY = F)
  if (verbose && !as.search)
    close(pb)
  
  
  # Helper function extract and combine mu and lppd from K lists with each
  # n/K of the elements to one list with n elements
  hf <- function(x) as.list(do.call(rbind, x))

  # Apply some magic to manipulate the structure of the list so that instead of
  # list with K sub_summaries each containing n/K mu:s and lppd:s, we have only
  # one sub_summary-list that contains with all n mu:s and lppd:s.
  sub <- apply(
    mapply(function(p_sub, fold) {
      lapply(.get_sub_summaries(p_sub, fold$d_test, family_kl), data.frame)
    }, p_sub_cv, list_cv),
    1, hf)

  ref <- hf(lapply(list_cv, function(fold) {
    data.frame(.weighted_summary_means(fold$d_test, family_kl, fold$w_test,
                                       fold$mu_test, fold$dis))
  }))

  # Combine also the K separate test data sets into one list
  # with n y's and weights's.
  d_cv <- hf(lapply(d_test_cv, function(d) {
    data.frame(d[c('y', 'weights')])
  }))

  list(vind_cv = vind_cv,
       summaries = list(sub = sub, ref = ref),
       pseudo_bma = NULL,
       d_test = c(d_cv, type = 'kfold'))
}


.get_kfold <- function(refmodel, K, verbose, seed) {
  # Fetch the k_fold list or compute it now if not already computed. This function will
  # return a list of length K, where each element is a list with fields 'refmodel' (object
  # of type refmodel computed by init_refmodel) and index list 'omitted' that denotes which
  # of the data points were left out for the corresponding fold.
	
  if (is.null(refmodel$cvfits)) {
    
    if (!is.null(refmodel$cvfun)) {
      
      # cv-function provided so perform the cross-validation now. In case refmodel
      # is datafit, cvfun will return an empty list and this will lead to normal cross-validation
      # for the submodels although we don't have an actual reference model
      if (verbose && !('datafit' %in% class(refmodel)))
        cat('Performing cross-validation for the reference model...\n')
      folds <- cvfolds(refmodel$nobs, k=K, seed=seed)
      cvfits <- refmodel$cvfun(folds)
      cvfits <- lapply(seq_along(cvfits), function(k) {
        # add the 'omitted' indices for the cvfits
        cvfit <- cvfits[[k]]
        cvfit$omitted <- which(folds==k)
        cvfit
      })
      
    } else
			# genuine probabilistic model but no k-fold fits nor cvfun provided, so raise an error
			stop('For a generic reference model, you must provide either cvfits or cvfun for k-fold cross-validation.
			      See function init_refmodel.')
	}	else
	  cvfits <- refmodel$cvfits
	
	# transform the cvfits-list to k_fold list, that is, initialize the reference models for each
	# fold given the prediction function and dispersion draws from the cvfits-list
	k_fold <- lapply(cvfits, function(cvfit) {
	  ref <- init_refmodel(z=refmodel$z[-cvfit$omitted,,drop=F], y=refmodel$y[-cvfit$omitted], x=refmodel$x[-cvfit$omitted,,drop=F],
												 family=refmodel$fam, predfun=cvfit$predfun, dis=cvfit$dis, offset=refmodel$offset[-cvfit$omitted],
												 wobs=refmodel$wobs[-cvfit$omitted], intercept=refmodel$intercept)
	  list(refmodel=ref, omitted=cvfit$omitted)
	})
	
	return(k_fold)
}



loo_varsel <- function(refmodel, method, nv_max, ns, nc, nspred, ncpred, relax, intercept, 
                       penalty, verbose, opt, nloo = NULL, validate_search = T, seed = NULL,
                       B = 20) {
	#
	# Performs the validation of the searching process using LOO.
	# validate_search indicates whether the selection is performed separately for each
  # fold (for each data point)
  #
  
	fam <- refmodel$fam
	mu <- refmodel$mu
	dis <- refmodel$dis
	n <- nrow(mu)

	# by default use all observations
	nloo <- ifelse(is.null(nloo), n, min(nloo, n))
	if (nloo < 1)
		stop('Value of \'nloo\' must be at least 1')

	# training data
	d_train <- .get_traindata(refmodel)

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

	# compute loo summaries for the reference model
	d_test <- d_train
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
		cat('Computing LOOs...\n')
		pb <- utils::txtProgressBar(min = 0, max = nloo, style = 3, initial=0)
	}

	if (!validate_search) {
	  # perform selection only once using all the data (not separately for each fold),
	  # and perform the projection then for each submodel size
	  # vind <- select(method, p_sel, d_train, fam, intercept, nv_max, penalty, verbose=F, opt)
	  spath <- select(method, p_sel, d_train, fam, intercept, nv_max, penalty, verbose=F, opt)
	  vind <- spath$vind
	}
	  
	for (run_index in seq_along(inds)) {
	  
	  # observation index
	  i <- inds[run_index]

	  # reweight the clusters/samples according to the is-loo weights
	  p_sel <- .get_p_clust(fam, mu, dis, wobs=refmodel$wobs, wsample=exp(lw[,i]), cl=cl_sel)
	  p_pred <- .get_p_clust(fam, mu, dis, wobs=refmodel$wobs, wsample=exp(lw[,i]), cl=cl_pred) 
	  
		if (validate_search) {
		  # perform selection with the reweighted clusters/samples
		  # vind <- select(method, p_sel, d_train, fam, intercept, nv_max, penalty, verbose=F, opt)
		  spath <- select(method, p_sel, d_train, fam, intercept, nv_max, penalty, verbose=F, opt)
		  vind <- spath$vind
		} 
	  
		# project onto the selected models and compute the prediction accuracy for the left-out point
	  as.search <- !relax && !is.null(spath$beta) && !is.null(spath$alpha)
	  submodels <- .get_submodels(spath, 0:nv_max, fam, p_pred,
	                              d_train, intercept, opt$regul, as.search=as.search)
		d_test <- list(x=matrix(refmodel$x[i,],nrow=1), y=refmodel$y[i], offset=d_train$offset[i], weights=d_train$weights[i])
		summaries_sub <- .get_sub_summaries(submodels, d_test, fam)

		for (k in seq_along(summaries_sub)) {
			loo_sub[i,k] <- summaries_sub[[k]]$lppd
			mu_sub[i,k] <- summaries_sub[[k]]$mu
		}
		vind_mat[i,] <- vind

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

  d_test <- list(y=d_train$y, weights=d_train$weights, type='loo')

  ## pseudo-BMA+ weights based on
  ## Y. Yao et al, "Using stacking to average Bayesian predictive distributions", 2017
  n <- length(inds)
  alpha <- rdirichlet(B, rep(1, n))
  z <- loo_sub
  z_b <- alpha %*% z * n
  ref <- exp(sum(loo_ref) - max(z_b))
  z_b <- exp(z_b - max(z_b))
  w <- rowMeans(sapply(1:nrow(z_b), function(i)
    z_b[i,] / (z_b[i,] + ref)))

	return(list(vind_cv=vind_cv, summaries=summaries, d_test=d_test, pseudo_bma=w))

}





.loo_subsample <- function(n, nloo, pareto_k, seed) {
  
  # decide which points to go through in the validation (i.e., which points
  # belong to the semi random subsample of validation points)
  
  # set random seed but ensure the old RNG state is restored on exit
  if (exists('.Random.seed')) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)
  
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  if (nloo < n) {
    
    bad <- which(pareto_k > 0.7)
    ok <- which(pareto_k <= 0.7 & pareto_k > 0.5)
    good <- which(pareto_k <= 0.5)
    inds <- resample(bad, min(length(bad), floor(nloo/3)) )
    inds <- c(inds, resample(ok, min(length(ok), floor(nloo/3))))
    inds <- c(inds, resample(good, min(length(good), floor(nloo/3))))
    if (length(inds) < nloo) {
      # not enough points selected, so choose randomly among the rest
      inds <- c(inds, resample(setdiff(1:n, inds), nloo-length(inds)))
    } 
    
    # assign the weights corresponding to this stratification (for example, the
    # 'bad' values are likely to be overpresented in the sample)
    w <- rep(0,n)
    w[inds[inds %in% bad]] <- length(bad) / sum(inds %in% bad)
    w[inds[inds %in% ok]] <- length(ok) / sum(inds %in% ok)
    w[inds[inds %in% good]] <- length(good) / sum(inds %in% good)
    
  } else {
    
    # all points used
    inds <- c(1:n)
    w <- rep(1,n)
  }
  
  # ensure weights are normalized
  w <- w/sum(w)
  
  return(list(inds=inds, w=w))
  
}
