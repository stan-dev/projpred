#' Cross-validate the variable selection (varsel)
#'
#' Perform cross-validation for the projective variable selection for a generalized
#' linear model.
#' @param fit Same as in \link[=varsel]{varsel}.
#' @param method Same as in \link[=varsel]{varsel}.
#' @param ns Same as in \link[=varsel]{varsel}.
#' @param nc Same as in \link[=varsel]{varsel}.
#' @param nv_max Same as in \link[=varsel]{varsel}.
#' @param intercept Same as in \link[=varsel]{varsel}.
#' @param verbose Whether to print out some information during the validation, Default is TRUE.
#' @param cv_method The cross-validation method, either 'LOO' or 'kfold'. Default is 'LOO'.
#' @param K Number of folds in the k-fold cross validation. Only applicable
#' if \code{cv_method = TRUE} and \code{k_fold = NULL}.
#' @param k_fold An array with cross-validated stanfits and the respective
#' test datasets returned by \link[=kfold]{kfold} with \code{save_fits=TRUE}.
#' If not provided, \link[=kfold]{kfold} is called inside the function.
#'
#' @return The original \link[=stanreg-objects]{stanreg} object augmented with an element 'varsel',
#' which is a list containing the following elements:
#' \describe{
#'  \item{\code{chosen}}{The order in which the variables were added to the submodel.}
#'  \item{\code{pctch}}{Percentage of cross-validation runs that included the given
#'    variable to a model of given size.}
#'  \item{\code{kl}}{KL-divergence for each submodel size.}
#'  \item{\code{summaries}}{CV-summary statistics computed during the selection.}
#'  \item{\code{d_test}}{The data used to evaluate the summaries.}
#'  \item{\code{family_kl}}{A modified \code{\link{family}}-object.}
#' }
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fit_cv <- cv_varsel(fit)
#' }
#'

#' @export
cv_varsel <- function(fit,  method = 'L1', cv_method = 'LOO', ns = NULL, nc = NULL,
                      nv_max = NULL, intercept = NULL, verbose = T,
                      K = NULL, k_fold = NULL, ...) {

	if ((is.null(ns) && is.null(nc)) || tolower(method)=='l1')
		# use one cluster for selection by default, and always with L1-search
		nc <- 1

	# .validate_for_varsel(fit)
	vars <- .extract_vars(fit)
	if(is.null(intercept))
		intercept <- vars$intercept
	if(is.null(nv_max) || nv_max > NCOL(vars$x)) {
		nv_max_default <- floor(0.4*length(vars$y)) # a somewhat sensible default limit for nv_max
		nv_max <- min(NCOL(vars$x), nv_max_default)
	}

	if (verbose)
		print(paste('Performing', method, 'search for the full model.'))
	sel <- varsel(fit, d_test=NULL, method=method, ns=ns, nv_max=nv_max, intercept=intercept, verbose=verbose)$varsel

	if(tolower(cv_method) == 'kfold') {
		sel_cv <- kfold_varsel(fit, method, nv_max, ns, nc, intercept, verbose, vars, K, k_fold)
	} else if (tolower(cv_method) == 'loo')  {
		sel_cv <- loo_varsel(fit, method, nv_max, ns, nc, intercept, verbose)
	} else {
		stop(sprintf('Unknown cross-validation method: %s.', method))
	}

	# find out how many of cross-validated iterations select
	# the same variables as the selection with all the data.
	pctch <- sapply(seq_along(sel$chosen), function(ind, chosen_array) {
		sum(chosen_array[, 1:ind] == sel$chosen[ind])/NROW(chosen_array)
	}, do.call(rbind, sel_cv$chosen_cv))

	fit$proj <- NULL
	fit$varsel <- c(sel[c('chosen', 'chosen_names', 'kl', 'family_kl')],
                  sel_cv[c('d_test', 'summaries')],
                  list(pctch = pctch))
	
	ssize <- .suggest_size(fit$varsel)
	fit$varsel$ssize <- ssize
	
	fit
}

kfold_varsel <- function(fit, method, nv_max, ns, nc, intercept, verbose, vars,
                         K, k_fold) {

  if (!('stanfit' %in% names(fit)))
  	stop('k-fold cross-validation not yet implemented for other than rstanarm reference models.')

  # KFOLD TEMPORARILY DISABLED
  stop(paste0('k-fold cross-validation is currently unavailable due to upcoming changes',
              ' in rstanarm. Please set cv_method to "LOO" to perform cross-validation.'))
  
  # Construct the kfold-objects. The resulting list contains an element 'fits',
  # which is a K x 2 dimensional array. Each row corresponds to one of the K
  # folds. First column contains the rstanarm-objects and the second column
  # the indices of the omitted observations (aka test data).
  if(is.null(k_fold)) {
    if(is.null(K)) K <- 10
    print(paste0('k_fold not provided, performing ', K,
                 '-fold cross-validation for the stan model.'))
    k_fold <- kfold(fit, K = K, save_fits = T)
  }
  family_kl <- kl_helpers(family(fit))

  # check that the fit-objects are valid for variable selection
  if(!all(apply(k_fold$fits, 1, function(fits, fit) {
    .validate_for_varsel(fits$fit)
    is.vector(fits$omitted) && max(fits$omitted) <= nobs(fit) && all(fits$omitted > 0)
  }, fit))) stop('k_fold does not have the correct form.')
  K <- attr(k_fold, 'K')

  # extract variables from each fit-object (stan-samples, x, y, etc.)
  # to a list of size K
  vars_cv <- lapply(k_fold$fits[,'fit'], .extract_vars)

  # List of size K with test data for each fold (note that vars is from
  # the full model, not from the cross-validated models).
  d_test <- lapply(k_fold$fits[,'omitted'], function(omitted) {
    list(x = vars$x[omitted,], y = vars$y[omitted],
         weights = vars$wobs[omitted], offset = vars$offset[omitted])
  })

  # List of K elements, each containing d_train, p_full, etc. corresponding
  # to each fold.
  msgs <- paste0(method, ' search for the fold number ', 1:K, '/', K, '.')
  list_cv <- mapply(function(vars, d_test, msg) {
    # .get_data_and_parameters(vars, d_test, intercept, ns, family_kl)
  	d_train <- .get_traindata(vars)
  	d_test <- d_test
  	p_full <- .get_refdist(vars, ns, nc)
  	mu <- family_kl$mu_fun(d_test$x, vars$alpha, vars$beta, d_test$offset)
  	dis <- vars$dis
  	list(d_train = d_train, d_test = d_test, p_full = p_full,
  	     mu_test = mu, dis_test = dis, w_test = vars$wsample, msg = msg)
  }, vars_cv, d_test, msgs, SIMPLIFY = F)

  # List of K elements, each a list of the variables selected for the
  # corresponding fold
  chosen_cv <- lapply(list_cv, function(x) {
    print(x$msg)
    select(method, x$p_full, x$d_train, family_kl, intercept, nv_max, verbose)
  })

  # Construct p_sub for each fold using .get_submodels.
  p_sub_cv <- mapply(function(chosen, x) {
    .get_submodels(chosen, c(0, seq_along(chosen)), family_kl, x$p_full,
                   x$d_train, intercept)
  }, chosen_cv, list_cv, SIMPLIFY = F)

  # Helper function extract and combine mu and lppd from K lists with each
  # n/K of the elements to one list with n elements
  hf <- function(x) as.list(do.call(rbind, x))

  # Apply some magic to manipulate the structure of the list so that instead of
  # list with K sub_summaries each containing n/K mu:s and lppd:s, we have only
  # one sub_summary-list that contains with all n mu:s and lppd:s.
  sub <- apply(
    mapply(function(p_sub, x, chosen) {
      lapply(.get_sub_summaries(chosen, x$d_test, p_sub, family_kl), data.frame)
    }, p_sub_cv, list_cv, chosen_cv),
    1, hf)

  full <- hf(lapply(list_cv, function(x) {
    data.frame(.weighted_summary_means(x$d_test, family_kl, x$w_test,
                                       x$mu_test, x$dis_test))
  }))

  # Combine also the K separate test data sets into one list
  # with n y's and weights's.
  d_cv <- hf(lapply(d_test, function(d) {
    data.frame(d[c('y', 'weights')])}))

  list(chosen_cv = chosen_cv, d_test = c(d_cv, type = 'kfold'),
       summaries = list(sub = sub, full = full))
}




loo_varsel <- function(fit, method, nv_max, ns, nc, intercept, verbose) {
	#
	# Performs the validation of the searching process using LOO.
	#
	#
	vars <- .extract_vars(fit)
	fam <- vars$fam
	mu <- vars$mu
	dis <- vars$dis

	# training data
	d_train <- .get_traindata(fit)

	# the reference distribution used for selection
	p_full <- .get_refdist(fit, ns=ns, nc=nc)
	cl <- p_full$cl # clustering information

	#### DEBUGGING ##
	# p_middle <- .get_refdist(fit, nc=3)
	# cl_middle <- p_middle$cl
	#################

	# fetch the log-likelihood for the full model to obtain the LOO weights
	if ('stanfit' %in% names(fit))
	    # stanreg-objects have a function log_lik
	    loglik <- log_lik(fit)
	else if (!is.null(fit$loglik))
	    # loglik given in the fit object (generic reference model)
	    loglik <- fit$loglik
	else
	    stop('To perform LOO for generic reference models, you must provide log-likelihood matrix to init_refmodel.')
	lw <- psislw(-loglik)$lw_smooth
	n <- dim(lw)[2]

	# compute loo summaries for the full model
	d_test <- d_train
	loo_full <- apply(loglik+lw, 2, 'log_sum_exp')
	mu_full <- rep(0,n)
	for (i in 1:n)
        mu_full[i] <- mu[i,] %*% exp(lw[,i])

	# initialize matrices where to store the results
	chosen_mat <- matrix(rep(0, n*nv_max), nrow=n)
	loo_sub <- matrix(nrow=n, ncol=nv_max+1)
	mu_sub <- matrix(nrow=n, ncol=nv_max+1)

	if (verbose) {
	    print('Start computing LOOs...')
	    nprints <- 10 # how many prints during the computation
	    print_at <- round( c(1:nprints)*(n/nprints) )
	    iprint <- 1
	}
		

	for (i in 1:n) {

		# reweight the clusters/samples according to the is-loo weights
		p_sel <- get_p_clust(fam, mu, dis, wobs=vars$wobs, wsample=exp(lw[,i]), cl=cl)
		# p_middle <- get_p_clust(fam, mu, dis, wobs=vars$wobs, wsample=exp(lw[,i]), cl=cl_middle)

		# perform selection
		chosen <- select(method, p_sel, d_train, fam, intercept, nv_max, verbose=F)
		chosen_mat[i,] <- chosen

		# project onto the selected models and compute the difference between
		# training and loo density for the left-out point
		p_sub <- .get_submodels(chosen, 0:nv_max, fam, p_sel, d_train, intercept) # replace p_sel by p_full here?
		# p_sub <- .get_submodels(chosen, 0:nv_max, fam, p_middle, d_train, intercept) # replace p_sel by p_full here?
		d_test <- list(x=matrix(vars$x[i,],nrow=1), y=vars$y[i], offset=d_train$offset[i], weights=d_train$weights[i])
		summaries_sub <- .get_sub_summaries(chosen, d_test, p_sub, fam)


		for (k in 0:nv_max) {
			loo_sub[i,k+1] <- summaries_sub[[k+1]]$lppd
			mu_sub[i,k+1] <- summaries_sub[[k+1]]$mu
		}

		if (verbose && (i %in% print_at)) {
		    print(sprintf('%d%% of LOOs done.', 10*iprint))
		    iprint <- iprint+1
		}
		    
	}

	###############
	## DEBUGGING ##
	# p_sel <- .get_refdist(fit, nc=1)
	# p_final <- .get_refdist(fit, nc=10)
	# chosen <- select(method, p_sel, d_train, fam, intercept, nv_max, verbose=F)
	# submod1 <- .get_submodels(chosen, 0:nv_max, fam, p_sel, d_train, intercept)
	# submod2 <- .get_submodels(chosen, 0:nv_max, fam, p_final, d_train, intercept)
	# summ1 <- .get_sub_summaries(chosen, d_train, submod1, fam)
	# summ2 <- .get_sub_summaries(chosen, d_train, submod2, fam)
	# peff <- matrix(0, nrow=n, ncol=length(summ1))
	# for (k in 1:length(summ1)) {
	# 	peff[,k] <- summ1[[k]]$lppd - loo_sub[,k]
	# 	print(mean(abs(summ2[[k]]$lppd - summ1[[k]]$lppd)))
	# }
	# print(apply(peff,2,sum))
	###############


	# put all the results together in the form required by cv_varsel
	summ_sub <-	lapply(0:nv_max, function(k){
	    list(lppd=loo_sub[,k+1], mu=mu_sub[,k+1])
	    # list(lppd= summ2[[k+1]]$lppd - peff[,k+1], mu=mu_sub[,k+1])
	})
	summ_full <- list(lppd=loo_full, mu=mu_full)
	summaries <- list(sub=summ_sub, full=summ_full)

    chosen_cv <- lapply(1:n, function(i){ chosen_mat[i,] })

    d_test <- list(y=d_train$y, weights=d_train$weights, type='loo')

	return(list(chosen_cv=chosen_cv, summaries=summaries, d_test=d_test))

}
