.onAttach <- function(...) {
  ver <- utils::packageVersion("projpred")
  packageStartupMessage("This is projpred version ", ver)
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

    if ('stanreg' %in% class(fit)) {
        
        # a stanreg object
        
        if(!(gsub('rstanarm::', '', fit$call[1]) %in% c("stan_glm", "stan_lm")))
            stop('Only \'stan_lm\' and \'stan_glm\' are supported.')
        
        families <- c('gaussian','binomial','poisson')
        if(!(family(fit)$family %in% families))
            stop(paste0('Only the following families are supported:\n',
                        paste(families, collapse = ', '), '.'))
    
    } else if ('refmodel' %in% class(fit)) {
        # a fit object constructed by init_refmodel, so everything should be fine
        return()
    } else {
        stop('The class for the provided object is not recognized.')
    }
    
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
		x <- x[, as.logical(attr(x, 'assign')), drop=F]
		attr(x, 'assign') <- NULL

		# undo the random permutation to make results reproducible
		perm_inv <- c(mapply(function(p, i) order(p) + i*length(p),
							 fit$stanfit@sim$permutation,1:fit$stanfit@sim$chains-1))
		dis_name <- ifelse(grepl(fit$call[1], 'stan_lm'), 'sigma', 'aux')

		alpha <- unname(drop(e$alpha %ORifNULL% rep(0, NROW(e$beta))))[perm_inv]
		beta <- t(unname(as.matrix(drop(e$beta))))[, perm_inv, drop=F]
		
		res <- list(
			fam = fam,
			x = x,
			dis = unname(e[[dis_name]])[perm_inv] %ORifNULL% rep(NA, NROW(e$beta)),
			offset = fit$offset %ORifNULL% rep(0, nobs(fit)),
			coefnames = coefnames,
			intercept = as.logical(attr(fit$terms,'intercept') %ORifNULL% 0)
			)

		res$predfun <- function(x, offset) fam$mu_fun(x, alpha, beta, offset) #
		res$mu <- res$predfun(x, res$offset)
		res$wsample <- rep(1/NCOL(res$mu), NCOL(res$mu)) # equal sample weights by default

		# y and the observation weights in a standard form
		target <- .get_standard_y(unname(get_y(fit)), weights(fit), fam)
		res$wobs <- target$weights
		res$y <- target$y
		res$nobs <- length(res$y) # this assumes a single output model
		res$loglik <- t(fam$ll_fun(res$mu,res$dis,res$y,res$wobs))

		return(res)

	} else if ('refmodel' %in% class(fit)) {
		# an object constructed by init_refmodel so all the relavant fields should be there
		return(fit)
	} else {
		stop('The class for the provided object is not recognized.')
	}
}


.get_standard_y <- function(y, weights, fam) {
  # return y and the corresponding observation weights into the 'standard' form:
  # for binomial family, y is transformed into a vector with values between 0 and 1,
  # and weights give the number of observations at each x.
  # for all other families, y and weights are kept as they are (unless weights is
  # a vector with length zero in which case it is replaced by a vector of ones).
  if(NCOL(y) == 1) {
    # weights <- if(length(weights) > 0) unname(weights) else rep(1, length(y))
    if(length(weights) > 0) 
      weights <- unname(weights)
    else 
      weights <- rep(1, length(y))
    if (fam$family == 'binomial') {
      if (is.factor(y))
        y <- as.vector(y, mode='integer') - 1 # zero-one vector
      else {
        if (any(y < 0 | y > 1))
          stop("y values must be 0 <= y <= 1 for the binomial model.")
      }
    }
  } else if (NCOL(y) == 2) {
    weights <- rowSums(y)
    y <- y[, 1] / weights
  } else {
    stop('y cannot have more than two columns.')
  }
  return(list(y=y,weights=weights))
}



.get_refdist <- function(fit, ns=NULL, nc=NULL, seed=NULL) {
	#
	# Creates the reference distribution based on the fit-object, and the
	# desired number of clusters (nc) or number of subsamples (ns). Returns
	# a list with fields 
	#       mu, var, weights, cl
	# TODO EXPLAIN THESE
	#
	# If nc is specified, then clustering is used and ns is ignored.
	# It is possible to use this function by passing .extract_vars(fit) as
	# an argument in place of fit.
	#
  if (is.null(seed))
    seed <- 1

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
			p_ref <- .get_p_clust(fam, vars$mu, vars$dis, wobs=vars$wobs, cl=cl)
		} else if (nc == NCOL(vars$mu)) {
		    # number of clusters equal to the number of samples, so return the samples
		    return(.get_refdist(fit, ns=nc))
		} else {
			# several clusters
		    if (nc > NCOL(vars$mu))
		        stop('The number of clusters nc cannot exceed the number of columns in mu.')
			p_ref <- .get_p_clust(fam, vars$mu, vars$dis, wobs=vars$wobs, nc=nc)
		}
	} else if (!is.null(ns)) {
		# subsample from the full model
		# would it be safer to actually randomly draw the subsample?
		if (ns > NCOL(vars$mu))
			stop('The number of subsamples ns cannot exceed the number of columns in mu.')
		s_ind <- round(seq(1, S, length.out  = ns))
		cl <- rep(NA, S)
		cl[s_ind] <- c(1:ns)
		predvar <- sapply(s_ind, function(j) { fam$predvar(vars$mu[,j,drop=F], vars$dis[j]) })
		p_ref <- list(mu = vars$mu[, s_ind, drop=F], var = predvar, dis = vars$dis[s_ind], weights = rep(1/ns, ns), cl=cl)
	} else {
		# use all the samples from the full model
		predvar <- sapply(1:S, function(j) { fam$predvar(vars$mu[,j,drop=F], vars$dis[j])	})
		p_ref <- list(mu = vars$mu, var = predvar, dis = vars$dis, weights = vars$wsample, cl=c(1:S))
	}

	# restore the old seed
	.Random.seed <- seed_old

	return(p_ref)
}



.get_p_clust <- function(family_kl, mu, dis, nc=10, wobs=rep(1,dim(mu)[1]), wsample=rep(1,dim(mu)[2]), cl = NULL) {
	# Function for perfoming the clustering over the samples.
	#
	# cluster the samples in the latent space if no clustering provided
	if (is.null(cl)) {
		f <- family_kl$linkfun(mu)
		out <- kmeans(t(f), nc, iter.max = 50)
		cl <- out$cluster # cluster indices for each sample
	} else if (typeof(cl)=='list') {
		# old clustering solution provided, so fetch the cluster indices
		if (is.null(cl$cluster))
			stop('argument cl must be a vector of cluster indices or a clustering object returned by k-means.')
		cl <- cl$cluster
	}
	
	# (re)compute the cluster centers, because they may be different from the ones
	# returned by kmeans if the samples have differing weights
	nc <- max(cl, na.rm=T) # number of clusters (assumes labeling 1,...,nc)
	centers <- matrix(0, nrow=nc, ncol=dim(mu)[1])
	wcluster <- rep(0,nc) # cluster weights
	eps <- 1e-10
	for (j in 1:nc) {
		# compute normalized weights within the cluster, 1-eps is for numerical stability
		ind <- which(cl==j)
		ws <- wsample[ind]/sum(wsample[ind])*(1-eps)
		
		# cluster centers and their weights
		centers[j,] <- mu[,ind,drop=F] %*% ws
		wcluster[j] <- sum(wsample[ind]) # unnormalized weight for the jth cluster
	}
	wcluster <- wcluster/sum(wcluster)
	
	# predictive variances
	predvar <- sapply(1:nc, function(j) {
		# compute normalized weights within the cluster, 1-eps is for numerical stability
		ind <- which(cl == j)
		ws <- wsample[ind]/sum(wsample[ind])*(1-eps)
		family_kl$predvar( mu[,ind,drop=F], dis[ind], ws )
	})
	
	# combine the results
	p <- list(mu = unname(t(centers)),
						var = predvar,
						weights = wcluster,
						cl = cl)
	return(p)
}



.get_traindata <- function(fit) {
	#
	# Returns the training data fetched from the fit object.
	# It is possible to use this function by passing .extract_vars(fit) as
	# an argument in place of fit.
	#
	if ( all(c('x', 'y', 'wobs', 'offset') %in% names(fit)) )
		# all the relevant fields contained in the given structure
		vars <- fit
	else
		# fetch the relevant info from the fit object
		vars <- .extract_vars(fit)
	
	return(list(x = vars$x, y = vars$y, weights = vars$wobs, offset = vars$offset))
}

.check_data <- function(data) {
	#
	# Check that data object has the correct form for internal use. The object must
	# be a list with with fields 'x', 'y', 'weights' and 'offset'.
	# Raises error if x or y is missing, but fills weights and offset with default
	# values if missing.
	#
	if (is.null(data$x)) stop('The data object must be a list with field x giving the predictor values.')
	if (is.null(data$y)) stop('The data object must be a list with field y giving the target values.')
	if (is.null(data$weights)) data$weights <- rep(1, nrow(data$x))
	if (is.null(data$offset)) data$offset <- rep(0, nrow(data$x))
	return(data)
}


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
                  & delta == TRUE & data %in% c('train', 'loo', 'kfold'))
  
  if (!all(is.na(stats[,'value']))) {
  	
  	mlpd_null <- subset(stats, size == 0, 'value')
  	mlpd_cutoff <- cutoff_pct*mlpd_null
  	res <- subset(stats, lq >= mlpd_cutoff$value, 'size')
  	if(nrow(res) == 0) {
  		ssize <- NA
  	} else {
  		ssize <- min(subset(stats, lq >= mlpd_cutoff$value, 'size'))
  	}
  } else {
  	# special case; all values compared to the reference model are NA indicating
  	# that the reference model is missing, so suggest the smallest model which
    # has its mlpd estimate within one standard deviation of the highest mlpd estimate,
    # i.e. is contained in the 68% central region
  	stats <- subset(.bootstrap_stats(varsel, alpha = 0.32), statistic == 'mlpd'
  				 					& delta == F & data %in% c('loo', 'kfold'))
  	imax <- which.max(unname(unlist(stats['value'])))
  	thresh <- stats[imax, 'lq']
  	ssize <- min(subset(stats, value >= thresh, 'size'))
  }
  ssize
}

.df_to_model_mat <- function(dfnew, var_names) {
  f <- formula(paste('~', paste(c('0', var_names), collapse = ' + ')))
  model.matrix(terms(f, keep.order = T), data = dfnew)
}

.is_proj_list <- function(proj) { !( 'family_kl' %in% names(proj) ) }

.unlist_proj <- function(p) if(length(p) == 1) p[[1]] else p
