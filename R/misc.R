.onAttach <- function(...) {
  ver <- utils::packageVersion("projpred")
  msg <- paste0('This is projpred version ', ver, '\n\n',
                 'Note: The type of the returned objects of varsel/cv_varsel have changed\n',
                 'since the latest release, although this does not affect how the functions\n',
                 'are used.')
  packageStartupMessage(msg)
}

weighted.sd <- function(x, w, na.rm=F) {
	if (na.rm) {
		ind <- !is.na(w) & !is.na(x)
		n <- sum(ind)
	} else {
		n <- length(x)
		ind <- rep(T,n)
	}
	w <- w/sum(w[ind])
  m <- sum(x[ind]*w[ind])
  sqrt(n/(n-1)*sum(w[ind]*(x[ind] - m)^2))
}

weighted.cov <- function(x,y, w, na.rm=F) {
	if (na.rm) {
		ind <- !is.na(w) & !is.na(x) & !is.na(y)
		n <- sum(ind)
	} else {
		n <- length(x)
		ind <- rep(T,n)
	}
	w <- w/sum(w[ind])
	mx <- sum(x[ind]*w[ind])
	my <- sum(y[ind]*w[ind])
	n/(n-1)*sum(w[ind]*(x[ind] - mx)*(x[ind] - my))
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

bootstrap <- function(x, fun=mean, b=1000, oobfun=NULL, seed=NULL, ...) {
  #
  # bootstrap an arbitrary quantity fun that takes the sample x
  # as the first input. other parameters to fun can be passed in as ...
  # example: boostrap(x,mean)
  #
  
  # set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- rngtools::RNGseed()
  on.exit(rngtools::RNGseed(rng_state_old))
  set.seed(seed)
  
  bsstat <- rep(NA, b)
  oobstat <- rep(NA, b)
  for (i in 1:b) {
    bsind <- sample(seq_along(x), replace=T)
    bsstat[i] <- fun(x[bsind], ...)
    if (!is.null(oobfun)) {
      oobind <- setdiff(seq_along(x), unique(bsind))
      oobstat[i] <- oobfun(x[oobind], ...)
    }
  }
  if (!is.null(oobfun)) {
    return(list(bs=bsstat, oob=oobstat))
  } else
    return(bsstat)
}


.bbweights <- function(N,B) {
  # generate Bayesian bootstrap weights, N = original sample size,
  # B = number of bootstrap samples
  bbw <- matrix(rgamma(N*B, 1), ncol = N)
  bbw <- bbw/rowSums(bbw)
  return(bbw)
}



# from rstanarm
`%ORifNULL%` <- function(a, b) if (is.null(a)) b else a




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



.get_refdist <- function(refmodel, ns=NULL, nc=NULL, seed=NULL) {
	#
	# Creates the reference distribution based on the refmodel-object, and the
	# desired number of clusters (nc) or number of subsamples (ns). If nc is specified,
  # then clustering is used and ns is ignored. Returns a list with fields:
  #
	#   mu: n-by-s matrix, vector of expected values for y for each draw/cluster. here s
  #       means either the number of draws ns or clusters nc used, depending on which one is used.
  #   var: n-by-s matrix, vector of predictive variances for y for each draw/cluster which
  #         which are needed for projecting the dispersion parameter (note that this can be
  #         unintuitively zero for those families that do not have dispersion)
  #   weights: s-element vector of weights for the draws/clusters
  #   cl: cluster assignment for each posterior draw, that is, a vector that has length equal to the
  #       number of posterior draws and each value is an integer between 1 and s
	# 
  if (is.null(seed))
    seed <- 17249420

  # set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- rngtools::RNGseed()
  on.exit(rngtools::RNGseed(rng_state_old))
  set.seed(seed)
  
	fam <- refmodel$fam
	S <- NCOL(refmodel$mu) # number of draws in the reference model

	if (!is.null(nc)) {
		# use clustering (ignore ns argument)
		if (nc == 1) {
			# special case, only one cluster
			cl <- rep(1, S)
			p_ref <- .get_p_clust(fam, refmodel$mu, refmodel$dis, wobs=refmodel$wobs, cl=cl)
		} else if (nc == NCOL(refmodel$mu)) {
		    # number of clusters equal to the number of samples, so return the samples
		    return(.get_refdist(refmodel, ns=nc))
		} else {
			# several clusters
		    if (nc > NCOL(refmodel$mu))
		        stop('The number of clusters nc cannot exceed the number of columns in mu.')
			p_ref <- .get_p_clust(fam, refmodel$mu, refmodel$dis, wobs=refmodel$wobs, nc=nc)
		}
	} else if (!is.null(ns)) {
		# subsample from the reference model
		# would it be safer to actually randomly draw the subsample?
		if (ns > NCOL(refmodel$mu))
			stop('The number of subsamples ns cannot exceed the number of columns in mu.')
		s_ind <- round(seq(1, S, length.out  = ns))
		cl <- rep(NA, S)
		cl[s_ind] <- c(1:ns)
		predvar <- sapply(s_ind, function(j) { fam$predvar(refmodel$mu[,j,drop=F], refmodel$dis[j]) })
		p_ref <- list(mu = refmodel$mu[, s_ind, drop=F], var = predvar, dis = refmodel$dis[s_ind], weights = rep(1/ns, ns), cl=cl)
	} else {
		# use all the draws from the reference model
		predvar <- sapply(1:S, function(j) { fam$predvar(refmodel$mu[,j,drop=F], refmodel$dis[j])	})
		p_ref <- list(mu = refmodel$mu, var = predvar, dis = refmodel$dis, weights = refmodel$wsample, cl=c(1:S))
	}

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



.get_traindata <- function(refmodel) {
	#
	# Returns the training data fetched from the reference model object.
	return(list(z = refmodel$z, x = refmodel$x, y = refmodel$y, weights = refmodel$wobs, offset = refmodel$offset))
}

.check_data <- function(data) {
	#
	# Check that data object has the correct form for internal use. The object must
	# be a list with with fields 'x', 'y', 'weights' and 'offset'.
	# Raises error if x or y is missing, but fills weights and offset with default
	# values if missing.
	#
	if (is.null(data$z)) stop('The data object must be a list with field z giving the reference model inputs.')
	if (is.null(data$x)) stop('The data object must be a list with field x giving the feature values.')
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

.augmented_x <- function(x, intercept) {
  if (intercept)
    return(cbind(1, x))
  else
    return(x)
}

.nonaugmented_x <- function(x, intercept) {
  if (intercept) {
    if (ncol(x) == 1)
      # there is only the column of ones in x, so return empty matrix
      return(matrix(nrow=nrow(x), ncol=0))
    else
      return(x[,2:ncol(x),drop=F])
  } else
    return(x)
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


.df_to_model_mat <- function(dfnew, var_names) {
  f <- formula(paste('~', paste(c('0', var_names), collapse = ' + ')))
  model.matrix(terms(f, keep.order = T), data = dfnew)
}

.is_proj_list <- function(proj) { !( 'family_kl' %in% names(proj) ) }

.unlist_proj <- function(p) if(length(p) == 1) p[[1]] else p
