#' Extract draws of the linear predictor and draw from the predictive
#' distribution of the projected submodel
#'
#' \code{proj_linpred} extracts draws of the linear predictor and
#' \code{proj_predict} draws from the predictive distribution of the projected
#' submodel or submodels. If the projection has not been performed, the
#' functions also perform the projection.
#'
#' @name proj-pred
#'
#' @param object The object returned by \link{varsel},
#' \link{cv_varsel} or \link{project}.
#' @param xnew The predictor values used in the prediction. If \code{vind} is
#' specified, then \code{xnew} should either be a dataframe containing column names
#' that correspond to \code{vind} or a matrix with the number and order of columns
#' corresponding to \code{vind}. If \code{vind} is unspecified, then \code{xnew} must
#' either be a dataframe containing all the column names as in the original data or a matrix
#' with the same columns at the same positions as in the original data.
#' @param ynew New (test) target variables. If given, then the log predictive density
#' for the new observations is computed.
#' @param offsetnew Offsets for the new observations. By default a vector of
#' zeros.
#' @param weightsnew Weights for the new observations. For binomial model,
#' corresponds to the number trials per observation. For \code{proj_linpred},
#' this argument matters only if \code{ynew} is specified. By default a vector
#' of ones.
#' @param transform Should the linear predictor be transformed using the
#' inverse-link function? Default is \code{FALSE}. For \code{proj_linpred} only.
#' @param integrated If \code{TRUE}, the output is averaged over the
#' parameters. Default is \code{FALSE}. For \code{proj_linpred} only.
#' @param nv Number of variables in the submodel (the variable combination is
#' taken from the \code{varsel} information). If a list, then results for all
#' specified model sizes are returned. Ignored if \code{vind} is specified.
#' @param draws Number of draws to return from the predictive distribution of
#' the projection. The default is 1000.
#' For \code{proj_predict} only.
#' @param seed_samp An optional seed to use for drawing from the projection.
#' For \code{proj_predict} only.
#' @param ... Additional argument passed to \link{project} if \code{object}
#' is an object returned by \link{varsel} or \link{cv_varsel}.
#'
#' @return If the prediction is done for one submodel only (\code{nv} has length one
#' or \code{vind} is specified) and ynew is unspecified, a matrix or vector of
#' predictions (depending on the value of \code{integrated}). If \code{ynew} is specified,
#' returns a list with elements pred (predictions) and lpd (log predictive densities).
#' If the predictions are done for several submodel sizes, returns a list with one element
#' for each submodel.
#'
NULL

# The 'helper' for proj_linpred and proj_predict, ie. does all the
# functionality that is common to them. It essentially checks all the arguments
# and sets them to their respective defaults and then loops over the
# projections. For each projection, it evaluates the fun-function, which
# calculates the linear predictor if called from proj_linpred and samples from
# the predictive distribution if called from proj_predict.
proj_helper <- function(object, xnew, offsetnew, weightsnew, nv, seed_samp,
                        fun, ...) {

  if (is.null(offsetnew)) offsetnew <- rep(0, nrow(xnew))
  if (is.null(weightsnew)) weightsnew <- rep(1, nrow(xnew))

  if ('projection' %in% class(object) ||
      (length(object)>0 && 'projection' %in% class(object[[1]]))) {
    proj <- object
  } else {
    # reference model obtained, so run the projection
    proj <- project(object = object, nv = nv, ...)
  }

  if (!.is_proj_list(proj)) {
    proj <- list(proj)
  } else {
    # proj is not a projection object
    if(any(sapply(proj, function(x) !('family_kl' %in% names(x)))))
      stop(paste('proj_linpred only works with objects returned by',
                 ' varsel, cv_varsel or project'))
  }

  projected_sizes <- sapply(proj, function(x) NROW(x$beta))
  nv <- list(...)$nv %ORifNULL% projected_sizes

  if (!all(nv %in% projected_sizes))
    stop(paste0('Linear prediction requested for nv = ',
                paste(nv, collapse = ', '),
                ', but projection performed only for nv = ',
                paste(projected_sizes, collapse = ', '), '.'))

  projs <- Filter(function(x) NROW(x$beta) %in% nv, proj)
  names(projs) <- nv

  xnew_df <- is.data.frame(xnew)
  if (xnew_df) {
    terms <- unique(unlist(lapply(projs, function(x) names(x$vind))))
    xnew <- .df_to_model_mat(xnew, terms)
  }

  if (!is.matrix(xnew))
    stop('xnew not provided in the correct format. See ?proj-pred.')

  vind <- list(...)$vind
  if (!is.null(vind) && NCOL(xnew) != length(vind))
    stop(paste('The number of columns in xnew does not match with the given',
               'number of variable indices (vind).'))

  # save the old seed and initialize with the new one
  seed_old <- .Random.seed
  set.seed(seed_samp)

  preds <- lapply(projs, function(proj) {
    if (xnew_df) {
      xtemp <- xnew[, min(1, length(proj$vind)):length(proj$vind), drop = F]
    } else if (!is.null(vind)) {
      # columns of xnew are assumed to match to the given variable indices
      xtemp <- xnew
    } else {
      xtemp <- xnew[, proj$vind, drop = F]
    }
    mu <- proj$family_kl$mu_fun(xtemp, proj$alpha, proj$beta, offsetnew)

    fun(proj, mu, offsetnew, weightsnew)
  })

  # restore the old seed
  .Random.seed <- seed_old

  .unlist_proj(preds)
}

#' @rdname proj-pred
#' @export
proj_linpred <- function(object, xnew, ynew = NULL, offsetnew = NULL,
                         weightsnew = NULL, nv = NULL, transform = FALSE,
                         integrated = FALSE, ...) {

  # function to perform to each projected submodel
  fun <- function(proj, mu, offset, weights) {
    pred <- t(mu)
    if (!transform) pred <- proj$family_kl$linkfun(pred)
    if (integrated) {
      # average over the parameters
      pred <- as.vector( proj$weights %*% pred )
    } else if (!is.null(dim(pred)) && dim(pred)[1]==1) {
      # return a vector if pred contains only one row
      pred <- as.vector(pred)
    }

    if (!is.null(ynew)) {
      # compute also the log-density
      target <- .get_standard_y(ynew, weights, proj$family_kl)
      ynew <- target$y
      weights <- target$weights
      lpd <- proj$family_kl$ll_fun(mu, proj$dis, ynew, weights)
      if (integrated && !is.null(dim(lpd))) {
        lpd <- as.vector(apply(lpd, 1, log_weighted_mean_exp, proj$weights))
      } else if (!is.null(dim(lpd))) {
        lpd <- t(lpd)
      }
      list(pred = pred, lpd = lpd)
    } else {
      pred
    }
  }

  # proj_helper lapplies fun to each projection in object
  proj_helper(object = object, xnew = xnew, offsetnew = offsetnew,
              weightsnew = weightsnew, nv = nv, seed_samp = NULL, fun = fun,
              ...)
}

#' @rdname proj-pred
#' @export
proj_predict <- function(object, xnew, offsetnew = NULL, weightsnew = NULL,
                         nv = NULL, draws = NULL, seed_samp = NULL, ...) {

  # function to perform to each projected submodel
  fun <- function(proj, mu, offset, weights) {
    if(is.null(draws)) draws <- 1000
    draw_inds <- sample(x = seq_along(proj$weights), size = draws,
                        replace = TRUE, prob = proj$weights)

    t(sapply(draw_inds, function(i) {
      proj$family_kl$ppd_fun(mu[,i], proj$dis[i], weights)
    }))
  }

  # proj_helper lapplies fun to each projection in object
  proj_helper(object = object, xnew = xnew, offsetnew = offsetnew,
              weightsnew = weightsnew, nv = nv, seed_samp = seed_samp,
              fun = fun, ...)
}

#' Plotting or printing summary statistics related to variable selection
#'
#' \code{varsel_stats} can be used to obtain summary statistics related to
#' variable selection. The same statistics can be plotted with
#' \code{varsel_plot}.
#'
#' @name varsel-statistics
#'
#' @param object The object returned by \link[=varsel]{varsel} or
#' \link[=cv_varsel]{cv_varsel}.
#' @param ... Currently ignored.
#' @param nv_max Maximum submodel size for which the statistics are calculated.
#' @param stats A list of strings of statistics to calculate. Available
#' options are: mlpd, kl, mse (gaussian only), pctcorr (binomial only).
#' If \code{NULL}, set to varsel_plot plots only mlpd, but varsel_stats
#' return all the statistics.
#' @param type One of 'mean', 'lower', 'upper' indicating whether to compute mean,
#' or either the lower or upper credible bound. Upper and lower bounds are determined so
#' that \code{1-alpha} percent of the mass falls between them.
#' @param deltas If \code{TRUE}, the difference between the full model and the
#' submodel is returned instead of the actual value of the statistic.
#' Defaults to \code{FALSE}.
#' @param alpha A number indicating the desired coverage of the credible
#' intervals. Eg. \code{alpha=0.1} corresponds to 90\% probability mass
#' within the intervals. Defaults to \code{0.1}.
NULL

#' @rdname varsel-statistics
#' @export
varsel_plot <- function(object, ..., nv_max = NULL, stats = NULL, deltas = F, alpha = 0.1) {
  
	if(!('varsel' %in% names(object)))
	  stop(paste('The provided object doesn\'t contain information about the',
	             'variable selection. Run the variable selection first.'))

	if (all(is.na(object$varsel$summaries$full$mu))) {
		# no reference model (or the results missing for some other reason),
		# so cannot compute differences between the reference model and submodels
		refstat_found <- F
		if (deltas==T) {
			deltas <- F
			warning('Cannot use deltas = TRUE when there is no reference model, setting deltas = FALSE.')
		}
	} else
		refstat_found <- T
	
	if(is.null(stats)) 
	  stats <- 'mlpd' 
	
	# compute all the statistics and fetch only those that were asked
	tab <- .tabulate_stats(object$varsel, alpha)
	stats_table <- subset(tab, (tab$delta==deltas | tab$statistic=='kl') & tab$statistic %in% stats)
	stats_ref <- subset(stats_table, stats_table$size==Inf)
	stats_sub <- subset(stats_table, stats_table$size!=Inf)
	
	
	if(NROW(stats_sub) == 0) {
	    stop(paste0(ifelse(length(stats)==1, 'Statistics ', 'Statistic '),
	                paste0(unique(stats), collapse=', '), ' not available.'))
	}
	
	if(is.null(nv_max)) 
	  nv_max <- max(stats_sub$size)
	ylab <- if(deltas) 'Difference to the full model' else 'value'
	
	# make sure that breaks on the x-axis are integers
	n_opts <- c(4,5,6)
	n_possible <- Filter(function(x) nv_max %% x == 0, n_opts)
	n_alt <- n_opts[which.min(n_opts - (nv_max %% n_opts))]
	nb <- ifelse(length(n_possible) > 0, min(n_possible), n_alt)
	by <- ceiling(nv_max/min(nv_max, nb))
	breaks <- seq(0, by*min(nv_max, nb), by)
	minor_breaks <- if(by%%2 == 0)
	    seq(by/2, by*min(nv_max, nb), by)
	else
	  NULL

	# plot submodel results
	pp <- ggplot(data = subset(stats_sub, stats_sub$size <= nv_max), mapping = aes(x = size)) +
		geom_linerange(aes(ymin = lq, ymax = uq, alpha=0.1)) +
    geom_line(aes(y = value)) +
    geom_point(aes(y = value))
	
	if (refstat_found)
		# add reference model results if they exist
		pp <- pp + geom_hline(aes(yintercept = value), data = stats_ref,
													color = 'darkred', linetype=2)
	pp <- pp + 
		scale_x_continuous(breaks = breaks, minor_breaks = minor_breaks,
	                       limits = c(min(breaks), max(breaks))) +
	  labs(x = 'Number of variables in the submodel', y = ylab) +
	  theme(legend.position = 'none') +
	  facet_grid(statistic ~ ., scales = 'free_y') 
	pp
}

#' @rdname varsel-statistics
#' @export
varsel_stats <- function(object, ..., nv_max = NULL, type = 'mean', deltas = F, alpha=0.1) {
  
	if(!('varsel' %in% names(object)))
      stop(paste('The provided object doesn\'t contain information about the',
                   'variable selection. Run the variable selection first.'))

	if (all(is.na(object$varsel$summaries$full$mu)) && deltas==T) {
		# no reference model (or the results missing for some other reason),
		# so cannot compute differences between the reference model and submodels
		warning('Cannot compute statistics for deltas = TRUE when there is no reference model.')
	}
	
  tab <- .tabulate_stats(object$varsel, alpha=alpha)
  stats_table <- subset(tab, (tab$delta == deltas | tab$statistic == 'kl') & tab$size != Inf)
  stats <- as.character(unique(stats_table$statistic))

  # transform type to the names that appear in the statistic table, and pick the
  # required values
  type <- switch(type, mean='value', upper='uq', lower='lq')
  arr <- data.frame(sapply(stats, function(sname) {
      unname(subset(stats_table, stats_table$statistic == sname, type))
  }))
  arr <- cbind(size = unique(stats_table$size), arr)

  if(is.null(nv_max)) 
    nv_max <- max(stats_table$size)

  arr$vind <- c(NA, object$varsel$vind)
  if('pctch' %in% names(object$varsel))
    arr$pctch <- c(NA, diag(object$varsel$pctch[,-1]))

  subset(arr, arr$size <= nv_max)
}

#' Generic reference model initialization
#'
#' Initializes a structure that can be used as a reference fit for the
#' projective variable selection.
#' This function is provided to allow construction of the reference fit
#' using also other tools than \code{rstanarm}, because only certain specific
#' information is needed for the actual projection and variable selection.
#'
#' @param x Predictor matrix of dimension \code{n}-by-\code{D} containing the candidate
#' variables for selection (i.e. variables from which to select the submodel). Rows denote
#' the observations and columns the different variables. 
#' @param y Vector of length \code{n} giving the target variable values.
#' @param family \link{family} object giving the model family
#' @param predfun Function that takes a \code{nt}-by-\code{D} test predictor matrix as an input
#' (\code{nt} = # test points, \code{D} = # predictors) and outputs
#' a \code{nt}-by-\code{S} matrix of expected values for the target variable y,
#' each column corresponding to one posterior draw for the parameters in the reference model
#' (the number of draws \code{S} can also be 1).
#' The output should be computed without any offsets, these are automatically taken into account
#' internally, e.g. in cross-validation.
#' @param dis Vector of length \code{S} giving the posterior draws for the dispersion parameter
#' in the reference model if there is such a parameter in the model family. For Gaussian
#' observation model this is the noise std \code{sigma}.
#' @param offset Offset to be added to the linear predictor in the projection. (Same as in
#' function \code{glm}.)
#' @param wobs Observation weights. If omitted, equal weights are assumed.
#' @param wsample vector of length \code{S} giving the weights for the posterior draws. 
#' If omitted, equal weights are assumed.
#' @param intercept Whether to use intercept. Default is \code{TRUE}.
#' @param cvfits A list with K elements, each of which is a list with fields including at least
#' variables: tr, ts and predfun giving the training and test indices and prediction function
#' for each fold. Additionally each element can have field dis (dispersion samples for each fold)
#' if the model has a dispersion parameter. Can be omitted but needed for K-fold cross validation
#' for genuine reference models.
#'
#' @return An object that can be passed to all the functions that
#' take the reference fit as the first argument, such as \link{varsel}, \link{cv_varsel},
#' \link[=proj-pred]{proj_predict} and \link[=proj-pred]{proj_linpred}.

#' @export
init_refmodel <- function(x, y, family, predfun=NULL, dis=NULL, offset=NULL, 
													wobs=NULL, wsample=NULL, intercept=TRUE, cvfits=NULL) {

  n <- length(y)
	family <- kl_helpers(family)
	
	if (is.null(offset))
		offset <- rep(0, n)	
	
	if (is.null(predfun)) {
		# no prediction function given, so the 'reference model' will simply contain the
		# observed data as the fitted values
		predmu <- function(x,offset=0) matrix(rep(NA, NROW(x)))
		mu <- y
		proper_model <- FALSE
	}	else {
		# genuine reference mdoel. add impact of offset to the prediction function
		predmu <- function(x,offset=0) family$linkinv( family$linkfun(predfun(x)) + offset )
		mu <- predmu(x,offset)
		proper_model <- TRUE
	}
	
	if (proper_model)
		if (.has.dispersion(family) && is.null(dis))
			stop(sprintf('Family %s needs a dispersion parameter so you must specify input argument \'dis\'.', family$family))
	
	mu <- unname(as.matrix(mu))
	S <- NCOL(mu) # number of samples in the reference model
	
	if (is.null(dis))
		dis <- rep(0, S)
	if (is.null(wobs))
		wobs <- rep(1, n)
	if (is.null(wsample))
		wsample <- rep(1, S)
	if (is.null(intercept))
		intercept <- TRUE
	wsample <- wsample/sum(wsample)
	
	# compute log-likelihood
	if (proper_model)
	  loglik <- t(family$ll_fun(mu,dis,y,wobs))
	else
	  loglik <- NULL

	# figure out column names for the variables
	if (!is.null(colnames(x)))
		coefnames <- colnames(x)
	else
		coefnames <- paste0('x',1:ncol(x))

	# y and the observation weights in a standard form
	target <- .get_standard_y(y, wobs, family)
	
	# fetch information from the cross-validated fits and create a data structure
	# that will be understood by cv_varsel (or actually kfold_varsel)
	if (!is.null(cvfits)) {
		cvfits <- sapply(cvfits, function(fold) {
			# fold must contain: tr,ts,predfun,(dis),(wsample)
			tr <- fold$tr
			ts <- fold$ts
			fit <- init_refmodel(x[tr,], y[tr], family, predfun=fold$predfun, dis=fold$dis,
													 offset=offset[tr], wobs=wobs[tr], wsample=fold$wsample, intercept=intercept, cvfits=NULL)
			list(fit=fit, omitted=ts)
		})
		k_fold <- list(fits=t(cvfits))
	} else
		k_fold <- cvfits
    
	fit <- list(x=x, y=target$y, fam=family, mu=mu, dis=dis, nobs=length(y), coefnames=coefnames,
							offset=offset, wobs=target$weights, wsample=wsample, intercept=intercept, 
							predfun=predmu, loglik=loglik, k_fold=k_fold)
	
	# define the class of the retuned object to be 'refmodel' and additionally 'datafit'
	# if only the observed data was provided and no actual function for predicting test data
	class(fit) <- 'refmodel'
	if (!proper_model) 
		class(fit) <- c(class(fit),'datafit')
	
	return(fit)
}


#' @method as.matrix projection
#' @export
as.matrix.projection <- function(x, ...) {
  if (x$p_type) {
    warning(paste0('Note, that projection was performed using',
                   'clustering and the clusters might have different weights.'))
  }
  res <- t(x$beta)
  if (ncol(res) > 0) colnames(res) <- names(x$vind)
  if (x$intercept) res <- cbind('(Intercept)' = x$alpha, res)
  if (x$family_kl$family == 'gaussian') res <- cbind(res, sigma = x$dis)
  res
}


#' Create cross-validation indices
#'
#' Divide indices from 1 to \code{n} into subsets for \code{k}-fold cross validation.
#' This function is potentially useful for creating the cross-validation objects for 
#' \link[=init_refmodel]{init_refmodel}.
#'
#' @param n Number of data points.
#' @param k Number of folds. 
#' @param out Format of the output, either 'foldwise' (default) or 'indices'. See below for details.
#' @param seed Random seed so that the same division could be obtained again if needed.
#'
#' @return If \code{out} is 'foldwise', the returned value is a list with \code{k} elements, 
#' each having fields \code{tr} and \code{ts} which give the training and test indices, respectively,
#' for each fold. If \code{out} is 'indices', the returned value is a list with fields \code{tr} and \code{ts}
#' each of which is a list with \code{k} elements giving the training and test indices for each fold.
#' @examples
#' \donttest{
#' ### compute sample means within each fold
#' n <- 100
#' y <- rnorm(n)
#' cv <- cvind(n, k=5)
#' cvmeans <- lapply(cv, function(fold) mean(y[fold$itr]))
#' }

#' @export
cvind <- function(n, k, out='foldwise', seed=NULL) {

	ind <- c(1:n)
	if (!is.null(seed)) {
		# save the old seed and initialize with the new one
		seed_old <- .Random.seed
		set.seed(seed)
	}
	
	# shuffle the indices
	ind <- sample(ind, n, replace=FALSE)
	
	if (out == 'foldwise') {
		cv <- lapply(1:k, function(i) {
			ts <- sort(ind[seq(i,n,k)])  # test set
			tr <- setdiff(1:n, ts) # training set
			list(tr=tr,ts=ts)
		})
	}	else if (out == 'indices') {
		cv <- list()
		cv$tr <- list()
		cv$ts <- list()
		for (i in 1:k) {
			ts <- sort(ind[seq(i,n,k)]) # test set
			tr <- setdiff(1:n, ts) # training set
			cv$tr[[i]] <- tr
			cv$ts[[i]] <- ts
		}
	} else
		stop(paste0('Unknown output format requested: ', out))
	
	if (!is.null(seed))
		# restore the old seed
		.Random.seed <- seed_old
	
	return(cv)
}



