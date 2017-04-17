.onAttach <- function(...) {
  ver <- utils::packageVersion("glmproj")
  packageStartupMessage("This is glmproj version ", ver)
}

log_weighted_mean_exp <- function(x, w) {
  x <- x + log(w)
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

log_sum_exp <- function(x) {
	max_x <- max(x)
	max_x + log(sum(exp(x - max_x)))
}

# check if the fit object is suitable for variable selection
.validate_for_varsel <- function(fit) {

  if(!('stanreg' %in% class(fit)))
    stop('Object is not a stanreg object')

  if(!(gsub('rstanarm::', '', fit$call[1]) %in% c("stan_glm", "stan_lm")))
    stop('Only \'stan_lm\' and \'stan_glm\' are supported.')

  families <- c('gaussian','binomial','poisson')
  if(!(family(fit)$family %in% families))
    stop(paste0('Only the following families are supported:\n',
                paste(families, collapse = ', '), '.'))

  if(NCOL(get_x(fit)) < 4)
    stop('Not enough explanatory variables for variable selection')
}


# from rstanarm
`%ORifNULL%` <- function(a, b) if (is.null(a)) b else a


# extract all important information from the fit object for variable selection
.extract_vars <- function(fit) {

	if (!is.null(fit$stanfit)) {

		# the fit is an rstanarm-object
		e <- extract(fit$stanfit)

		# family and the predictor matrix x
		fam <- kl_helpers(family(fit))
		x <- unname(get_x(fit))
		coefnames <- names(coef(fit))[as.logical(attr(x, 'assign'))]
		x <- x[, as.logical(attr(x, 'assign'))]
		attr(x, 'assign') <- NULL

		# undo the random permutation to make results reproducible
		perm_inv <- c(mapply(function(p, i) order(p) + i*length(p),
							 fit$stanfit@sim$permutation,1:fit$stanfit@sim$chains-1))
		dis_name <- ifelse(grepl(fit$call[1], 'stan_lm'), 'sigma', 'aux')

		res <- list(
			fam = fam,
			x = x,
			alpha = unname(drop(e$alpha %ORifNULL% rep(0, NROW(e$beta))))[perm_inv], # EVENTUALLY NEED TO GET RID OFF THIS
			beta = t(unname(drop(e$beta)))[, perm_inv],                              # EVENTUALLY NEED TO GET RID OFF THIS
			dis = unname(e[[dis_name]]) %ORifNULL% rep(NA, nrow(e$beta))[perm_inv],
			offset = fit$offset %ORifNULL% rep(0, nobs(fit)),
			coefnames = coefnames,
			intercept = attr(fit$terms,'intercept') %ORifNULL% 0)

		res$mu <- fam$mu_fun(x, res$alpha, res$beta, res$offset)
		res$wsample <- rep(1/NCOL(res$mu), NCOL(res$mu)) # equal sample weights by default

		# y and the observation weights in a standard form
		temp <- .get_standard_y(unname(get_y(fit)), weights(fit))
		res$wobs <- temp$weights
		res$y <- temp$y

		# y <- unname(get_y(fit))
		# if(NCOL(y) == 1) {
		# 	res$wobs <- if(length(weights(fit))) unname(weights(fit)) else rep(1, nobs(fit))
		# 	if (is.factor(y))
		# 	    res$y <- as.vector(y, mode='integer') - 1 # zero-one vector
		# 	else
		# 	    res$y <- y
		# } else {
		# 	res$wobs <- rowSums(y)
		# 	res$y <- y[, 1] / res$wobs
		# }

		return(res)

	} else {

		# not and rstanarm-object, so look for the relevant fields

	    # DUMMY, simply return the object itself and assume it has all the relevant fiels
	    # (i.e. it was created by init_refmodel)
	    return(fit)
	    # if (!is.null(fit$x)) stop('')
		# stop('Other than rstanarm-fits are currently not supported, but will be in the near future.')
	}
}


.get_standard_y <- function(y, weights) {
    # return y and the corresponding observation weights into the 'standard' form:
    # for binomial family, y is transformed into a vector with values between 0 and 1,
    # and weights give the number of observations at each x.
    # for all other families, y and weights are kept as they are (unless weights is
    # a vector with length zero in which case it is replaced by a vector of ones).
    if(NCOL(y) == 1) {
        weights <- if(length(weights) > 0) unname(weights) else rep(1, length(y))
        if (is.factor(y))
            y <- as.vector(y, mode='integer') - 1 # zero-one vector
    } else {
        weights <- rowSums(y)
        y <- y[, 1] / weights
    }
    return(list(y=y,weights=weights))
}


.get_data_and_parameters <- function(vars, d_test, intercept, ns, family_kl) {

    # THIS FUNCTION SEEMS TO BE DEPRECATED, AND COULD THUS BE REMOVED

  # - Returns d_train, d_test, p_full, coef_full.
  # - If d_test is NA, it is set to d_train.

  mu <- family_kl$mu_fun(vars$x, vars$alpha, vars$beta, vars$offset)

  d_train <- list(x = vars$x, weights = vars$wobs, offset = vars$offset)

  # if test data doesn't exist, use training data to evaluate mse, r2, mlpd
  if(is.null(d_test)) {
    # check that d_test is of the correct form?
    d_test <- vars[c('x', 'weights', 'offset', 'y')]
  } else {
    if(is.null(d_test$weights)) d_test$weights <- rep(1, nrow(d_test$x))
    if(is.null(d_test$offset)) d_test$offset <- rep(0, nrow(d_test$x))
  }
  # indices of samples that are used in the projection
  s_ind <- round(seq(1, ncol(vars$beta), length.out  = ns))
  p_full <- list(mu = mu[, s_ind], dis = vars$dis[s_ind],
                 weights = rep(1/ns, ns))
  coef_full <- list(alpha = vars$alpha[s_ind], beta = vars$beta[, s_ind])

  # cluster the samples of the full model if clustering is wanted
  # for the variable selection
  do_clust <- F
  clust <- if(do_clust) get_p_clust(family_kl, mu, vars$dis, ns) else NULL
  if(do_clust) p_full <- clust

  list(d_train = d_train, d_test = d_test, p_full=p_full,
       coef_full = coef_full)
}


.get_refdist <- function(fit, ns=NULL, nc=NULL, seed=NULL) {
	#
	# Creates the reference distribution based on the fit-object, and the
	# desired number of clusters (nc) or number of subsamples (ns). Returns
	# a list with fields mu, dis and weights. If nc is specified, then
	# clustering is used and ns is ignored.
	#
	# It is possible to use this function by passing .extract_vars(fit) as
	# an argument in place of fit.
	#

	# save the old seed and initialize with the new one
	seed_old <- .Random.seed
	set.seed(seed)

	if ( all(c('fam', 'x', 'mu', 'dis') %in% names(fit)) )
		# all the relevant fields contained in the given structure
		vars <- fit
	else
		# fetch the relevant info from the fit object
		vars <- .extract_vars(fit)

	fam <- vars$fam
	n <- NROW(vars$x) # number of data points
	S <- NCOL(vars$mu) # sample size in the full model

	if (!is.null(nc)) {
		# use clustering (ignore ns argument)
		if (nc == 1) {
			# special case, only one cluster
			cl <- rep(1, S)
			p_ref <- get_p_clust(fam, vars$mu, vars$dis, wobs=vars$wobs, cl=cl)
		} else if (nc == NCOL(vars$mu)) {
		    # number of clusters equal to the number of samples, so return the samples
		    return(.get_refdist(fit, ns=nc))
		} else {
			# several clusters
		    if (nc > NCOL(vars$mu))
		        stop('The number of clusters nc cannot exceed the number of columns in mu.')
			p_ref <- get_p_clust(fam, vars$mu, vars$dis, wobs=vars$wobs, nc=nc)
		}
	} else if (!is.null(ns)) {
		# subsample from the full model
		# would it be safer to actually randomly draw the subsample?
	    if (ns > NCOL(vars$mu))
	        stop('The number of subsamples ns cannot exceed the number of columns in mu.')
		s_ind <- round(seq(1, S, length.out  = ns))
		cl <- rep(NA, S)
		cl[s_ind] <- c(1:ns)
		p_ref <- list(mu = vars$mu[, s_ind, drop=F], dis = vars$dis[s_ind], weights = rep(1/ns, ns), cl=cl)
	} else {
		# use all the samples from the full model
		p_ref <- list(mu = vars$mu, dis = vars$dis, weights = rep(1/S, S), cl=c(1:S))
	}

	# restore the old seed
	.Random.seed <- seed_old

	return(p_ref)
}

.get_traindata <- function(fit) {
	#
	# Returns the training data fetched from the fit object.
	# It is possible to use this function by passing .extract_vars(fit) as
	# an argument in place of fit.
	#
	if ( all(c('x', 'y', 'weights', 'offset') %in% names(fit)) )
		# all the relevant fields contained in the given structure
		vars <- fit
	else
		# fetch the relevant info from the fit object
		vars <- .extract_vars(fit)

	return(list(x = vars$x, y = vars$y, weights = vars$wobs, offset = vars$offset))
}

.fill_offset_and_weights <- function(data) {
	#
	# Simply checks whether the offset and weight fields exist in data structure,
	# fills them in if needed.
	#
	if(is.null(data$weights)) data$weights <- rep(1, nrow(data$x))
	if(is.null(data$offset)) data$offset <- rep(0, nrow(data$x))
	return(data)
}

# .validate_ns_nc <- function(nc, ns, nc_max, ns_max) {
# 	if (is.null(nc) && is.null(ns))
# 		stop('Either nc or ns must be non-null.')
# 	if (!is.null(nc)) {
# 		if (nc < 0)
# 			stop('nc must be > 0.')
# 		if (nc > 100) {
#
# 		}
# 	} else {
# 		if(ns > NCOL(vars$mu)) {
# warning(paste0('Setting the number of samples to ', NCOL(vars$mu),'.'))
# ns <- NCOL(vars$mu)
# }
# 	}
# }


.split_coef <- function(b, intercept) {
  if(intercept) {
    list(alpha = b[1, ], beta = b[-1, , drop = F])
  } else {
    list(alpha = rep(0, NCOL(b)), beta = b)
  }
}


.varsel_errors <- function(e) {
  if(grepl('computationally singular', e$message)) {
    stop(paste(
      'Numerical problems with inverting the covariance matrix. Possibly a',
      'problem with the convergence of the stan model?, If not, consider',
      'stopping the selection early by setting the variable nv_max accordingly.'
    ))
  } else {
    stop(e$message)
  }
}

.suggest_size <- function(varsel, alpha = 0.1, cutoff_pct = 0.1) {
  # suggest a model size. Currently finds the smallest model for which
  # the lower alpha/2-quantile of mlpd is at least cutoff_pct from the full model
  stats <- subset(.bootstrap_stats(varsel, alpha = alpha), statistic == 'mlpd'
                  & delta == TRUE & data %in% c('loo', 'kfold'))
  mlpd_null <- subset(stats, size == 0, 'value')
  mlpd_cutoff <- cutoff_pct*mlpd_null
  res <- subset(stats, lq >= mlpd_cutoff$value, 'size')

  if(nrow(res) == 0) {
    NA
  } else {
    min(subset(stats, lq >= mlpd_cutoff$value, 'size'))
  }
}

.is_proj_list <- function(proj) { !( 'family_kl' %in% names(proj) ) }

.unlist_proj <- function(p) if(length(p) == 1) p[[1]] else p
