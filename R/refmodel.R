
#' Get reference model structure
#'
#' Generic function that can be used to create and fetch the reference model structure
#' for all those objects that have this method. All these implementations are wrappers
#' to the \code{\link{init_refmodel}}-function so the returned object has the same type.
#' 
#' @name get-refmodel
#'
#' @param object Object based on which the reference model is created. See possible types below. 
#' @param ... Arguments passed to the methods.
#'
#' @return An object of type \code{refmodel} (the same type as returned by \link{init_refmodel}) 
#' that can be passed to all the functions that
#' take the reference fit as the first argument, such as \link{varsel}, \link{cv_varsel}, \link{project},
#' \link[=proj-pred]{proj_predict} and \link[=proj-pred]{proj_linpred}.
#' 
#' @examples
#' \donttest{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' ref <- get_refmodel(fit)
#' print(class(ref))
#' 
#' # variable selection, use the already constructed reference model
#' vs <- varsel(ref) 
#' # this will first construct the reference model and then execute 
#' # exactly the same way as the previous command (the result is identical)
#' vs <- varsel(fit) 
#' }
#'
NULL

#' @rdname get-refmodel
#' @export
get_refmodel <- function (object, ...) {
	UseMethod("get_refmodel", object)
}

#' @rdname get-refmodel
#' @export
get_refmodel.refmodel <- function(object, ...) {
	# if the object is reference model already, then simply return it as is
	object
}

#' @rdname get-refmodel
#' @export
get_refmodel.vsel <- function(object, ...) {
	# the reference model is stored in vsel-object
	object$refmodel
}

#' @rdname get-refmodel
#' @export
get_refmodel.cvsel <- function(object, ...) {
	# the reference model is stored in cvsel object
	object$refmodel
}

#' @rdname get-refmodel
#' @export
get_refmodel.stanreg <- function(object, ...) {
	
	# the fit is an rstanarm-object
  
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    stop("You need package \"rstanarm\". Please install it.",
         call. = FALSE)
  }
	
	if ('lmerMod' %in% class(object))
		stop('stan_lmer and stan_glmer are not yet supported.')
	
	families <- c('gaussian','binomial','poisson')
	if (!(family(object)$family %in% families))
		stop(paste0('Only the following families are currently supported:\n',
								paste(families, collapse = ', '), '.'))
	
	# fetch the draws
	samp <- as.data.frame(object)
	ndraws <- nrow(samp)
	
	# data, family and the predictor matrix x
	z <- object$data # inputs of the reference model (this contains also the target or a transformation of it, but that shouldn't hurt)
	fam <- kl_helpers(family(object))
	x <- rstanarm::get_x(object)
	rownames(x) <- NULL # ignore the rownames
	x <- x[, as.logical(attr(x, 'assign')), drop=F] # drop the column of ones
	attr(x, 'assign') <- NULL
	
	y <- unname(rstanarm::get_y(object))
	dis <- samp$sigma %ORifNULL% rep(0, ndraws) # TODO: handle other than gaussian likelihoods..
	offset <- object$offset %ORifNULL% rep(0, nobs(object))
	intercept <- as.logical(attr(object$terms,'intercept') %ORifNULL% 0)
	predfun <- function(zt) t(rstanarm::posterior_linpred(object, newdata=data.frame(zt), transform=T, offset=rep(0,nrow(zt))))
	wsample <- rep(1/ndraws, ndraws) # equal sample weights by default
	wobs <- unname(weights(object)) # observation weights
	if (length(wobs)==0) wobs <- rep(1,nrow(z))
	
	# cvfun for k-fold cross-validation
	cvfun <- function(folds) {
	  cvres <- rstanarm::kfold(object, K = max(folds), save_fits = T, folds = folds)
	  fits <- cvres$fits[,'fit']
	  lapply(fits, function (fit) {
	    dis <- as.data.frame(fit)$sigma # NOTE: this works only for Gaussian family
	    predfun <- function(zt) t(rstanarm::posterior_linpred(fit, newdata=data.frame(zt), transform=T, offset=rep(0,nrow(zt))))
	    list(predfun=predfun, dis=dis)
	  })
	}
  
	init_refmodel(z=z, y=y, family=fam, x=x, predfun=predfun, dis=dis, offset=offset,
	              wobs=wobs, wsample=wsample, intercept=intercept, cvfits=NULL, cvfun=cvfun) 
}










#' Custom reference model initialization
#'
#' Initializes a structure that can be used as a reference fit for the
#' projective variable selection. This function is provided to allow construction 
#' of the reference fit from arbitrary fitted models, because only limited
#' information is needed for the actual projection and variable selection.
#'
#' @param z Predictor matrix of dimension \code{n}-by-\code{dz} containing the training
#' features for the reference model. Rows denote the observations and columns the different features. 
#' @param y Vector of length \code{n} giving the target variable values.
#' @param family \link{family} object giving the model family
#' @param x Predictor matrix of dimension \code{n}-by-\code{dx} containing the candidate
#' features for selection (i.e. variables from which to select the submodel).  Rows denote
#' the observations and columns the different features. Notice that this can
#' different from \code{z}. If missing, same as \code{z} by default.
#' @param predfun Function that takes a \code{nt}-by-\code{dz} test predictor matrix \code{zt} as an input
#' (\code{nt} = # test points, \code{dz} = number of features in the reference model) and outputs
#' a \code{nt}-by-\code{S} matrix of expected values for the target variable \code{y},
#' each column corresponding to one posterior draw for the parameters in the reference model
#' (the number of draws \code{S} can also be 1). Notice that the output should be computed without
#' any offsets, these are automatically taken into account internally, e.g. in cross-validation.
#' If omitted, then the returned object will be 'data reference', that is, it can be used to compute
#' penalized maximum likelihood solutions such as Lasso (see examples below and in the quickstart vignette.)
#' @param dis Vector of length \code{S} giving the posterior draws for the dispersion parameter
#' in the reference model if there is such a parameter in the model family. For Gaussian
#' observation model this is the noise std \code{sigma}.
#' @param offset Offset to be added to the linear predictor in the projection. (Same as in
#' function \code{glm}.)
#' @param wobs Observation weights. If omitted, equal weights are assumed.
#' @param wsample vector of length \code{S} giving the weights for the posterior draws. 
#' If omitted, equal weights are assumed.
#' @param intercept Whether to use intercept. Default is \code{TRUE}.
#' @param cvfun Function for performing K-fold cross-validation. The input is an \code{n}-element
#' vector where each value is an integer between 1 and K denoting the fold for each observation.
#' Should return a list with K elements, each of which is a list with fields \code{predfun} and
#' \code{dis} (if the model has a dispersion parameter) which are defined the same way as the arguments 
#' \code{predfun} and \code{dis} above but are computed using only the corresponding subset of the data. 
#' More precisely, if \code{cvres} denotes
#' the list returned by \code{cvfun}, then \code{cvres[[k]]$predfun} and \code{cvres[[k]]$dis} must be computed
#' using only data from indices \code{folds != k}, where \code{folds} is the \code{n}-element input for
#' \code{cvfun}. Can be omitted but either \code{cvfun} or \code{cvfits} is needed for K-fold cross-validation
#' for genuine reference models. See example below.
#' @param cvfits A list with K elements, that has the same format as the value returned by \code{cvind} but 
#' each element of \code{cvfits} must also contain a field \code{omitted} which indicates the indices that
#' were left out for the corresponding fold. Usually it is easier to specify \code{cvfun} but this can be useful
#' if you have already computed the cross-validation for the reference model and would like to avoid 
#' recomputing it. Can be omitted but either \code{cvfun} or \code{cvfits} is needed for K-fold cross-validation
#' for genuine reference models.
#' @param ... Currently ignored.
#'
#' @return An object that can be passed to all the functions that
#' take the reference fit as the first argument, such as \link{varsel}, \link{cv_varsel},
#' \link[=proj-pred]{proj_predict} and \link[=proj-pred]{proj_linpred}.
#' 
#' @examples
#' \donttest{
#' 
#' # generate some toy data
#' set.seed(1)
#' n <- 100
#' d <- 10
#' x <- matrix(rnorm(n*d), nrow=n, ncol=d)
#' b <- c(c(1,1),rep(0,d-2)) # first two variables are relevant
#' y <- x %*% b + rnorm(n)
#' 
#' # fit the model (this uses rstanarm for posterior inference, 
#' # but any other tool could also be used)
#' fit <- stan_glm(y~x, family=gaussian(), data=data.frame(x=I(x),y=y))
#' draws <- as.matrix(fit)
#' a <- draws[,1] # intercept
#' b <- draws[,2:(ncol(draws)-1)] # regression coefficients
#' sigma <- draws[,ncol(draws)] # noise std
#' 
#' # initialize the reference model structure
#' predfun <- function(xt) t( b %*% t(xt) + a )
#' ref <- init_refmodel(x,y, gaussian(), predfun=predfun, dis=sigma)
#' 
#' # variable selection based on the reference model
#' vs <- cv_varsel(ref)
#' varsel_plot(vs)
#' 
#' 
#' # pass in the original data as 'reference'; this allows us to compute 
#' # traditional estimates like Lasso
#' dref <- init_refmodel(x,y,gaussian())
#' lasso <- cv_varsel(dref, method='l1') # lasso
#' varsel_plot(lasso, stat='rmse')
#' 
#' }
#'

#' @export
init_refmodel <- function(z, y, family, x=NULL, predfun=NULL, dis=NULL, offset=NULL,
                          wobs=NULL, wsample=NULL, intercept=TRUE, cvfun=NULL, cvfits=NULL,  ...) {
	
	n <- NROW(z)
	family <- kl_helpers(family)
	
	if (is.null(x))
		x <- z
	if (is.null(offset))
		offset <- rep(0, n)	
	
	# y and the observation weights in a standard form
	target <- .get_standard_y(y, wobs, family)
	y <- target$y
	wobs <- target$weights
	
	if (is.null(predfun)) {
		# no prediction function given, so the 'reference model' will simply contain the
		# observed data as the fitted values
		predmu <- function(z,offset=0) matrix(rep(NA, NROW(z)))
		mu <- y
		proper_model <- FALSE
	}	else {
		# genuine reference model. add impact of offset to the prediction function
		predmu <- function(z,offset=0) family$linkinv( family$linkfun(predfun(z)) + offset )
		mu <- predmu(z,offset)
		if (NROW(y)!=NROW(mu)) 
			stop(paste0('The number of rows in the output of predfun(z) does not match with the given y;',
									'predfun seems to be misspecified.'))
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
	
	if (!proper_model) {
	  # this is a dummy definition for cvfun, but it will lead to standard cross-validation
	  # for datafit reference; see cv_varsel and get_kfold
	  cvfun <- function(folds) lapply(1:max(folds), function(k) list())
	}
	
	refmodel <- list(z=z, x=x, y=y, fam=family, mu=mu, dis=dis, nobs=n, coefnames=coefnames,
	                 offset=offset, wobs=wobs, wsample=wsample, intercept=intercept,
	                 predfun=predmu, loglik=loglik, cvfits=cvfits, cvfun=cvfun)
	
	# define the class of the retuned object to be 'refmodel' and additionally 'datafit'
	# if only the observed data was provided and no actual function for predicting test data
	class(refmodel) <- 'refmodel'
	if (!proper_model) 
		class(refmodel) <- c(class(refmodel),'datafit')
	
	return(refmodel)
}






#' Predict method for reference model objects
#'
#' Compute the predictions using the reference model, that is, compute the
#' expected value for the next observation, or evaluate the log-predictive
#' density at a given point.
#'
#' @param object The object of class \code{refmodel}.
#' @param znew Matrix of predictor values used in the prediction. 
#' @param ynew New (test) target variables. If given, then the log predictive density
#' for the new observations is computed.
#' @param offsetnew Offsets for the new observations. By default a vector of
#' zeros.
#' @param weightsnew Weights for the new observations. For binomial model,
#' corresponds to the number trials per observation. Has effect only if \code{ynew} is specified.
#' By default a vector of ones.
#' @param type Scale on which the predictions are returned. Either 'link' (the latent function
#' value, from -inf to inf) or 'response' (the scale on which the target \code{y} is measured, 
#' obtained by taking the inverse-link from the latent value).
#' @param ... Currently ignored.
#'
#' @return Returns either a vector of predictions, or vector of log predictive densities evaluated
#' at \code{ynew} if \code{ynew} is not \code{NULL}.

#' @export
predict.refmodel <- function(object, znew, ynew = NULL, offsetnew = NULL,
                             weightsnew = NULL, type = 'response', ...) {
	
	if ('datafit' %in% class(object))
		stop('Cannot make predictions with data reference only.')
	
	if (is.null(offsetnew)) offsetnew <- rep(0, nrow(znew))
	if (is.null(weightsnew)) weightsnew <- rep(1, nrow(znew))
	
	mu <- object$predfun(znew, offsetnew)
	
	if (is.null(ynew)) {
		
		if (type == 'link')
			pred <- object$family$linkfun(mu)
		else
			pred <- mu
		
		# integrate over the samples
		if (NCOL(pred) > 1)
			pred <- rowMeans(pred)
		
		return(pred)
		
	} else {
		
		# evaluate the log predictive density at the given ynew values
		loglik <- object$fam$ll_fun(mu, object$dis, ynew, weightsnew)
		S <- ncol(loglik)
		lpd <- apply(loglik, 1, log_sum_exp) - log(S)
		return(lpd)
	}
	
}





















