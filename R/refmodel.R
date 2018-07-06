
# TODO: should we make these available for the user
# #' @export
get_refmodel <- function (object, ...) {
	UseMethod("get_refmodel", object)
}


get_refmodel.refmodel <- function(object, ...) {
	# if the object is reference model already, then simply return it as is
	object
}

get_refmodel.vsel <- function(object, ...) {
	# the reference model is stored in vsel-object
	object$refmodel
}

get_refmodel.cvsel <- function(object, ...) {
	# the reference model is stored in cvsel object
	object$refmodel
}

# #' @export
get_refmodel.stanreg <- function(object, ...) {
	
	# the fit is an rstanarm-object
	
	# fetch the draws
	samp <- as.data.frame(object)
	ndraws <- nrow(samp)
	
	# family and the predictor matrix x
	fam <- kl_helpers(family(object))
	x <- get_x(object)
	rownames(x) <- NULL # ignore the rownames
	x <- x[, as.logical(attr(x, 'assign')), drop=F] # drop the column of ones
	attr(x, 'assign') <- NULL
	
	dis <- samp$sigma %ORifNULL% rep(0, ndraws) # TODO: handle other than gaussian likelihoods..
	offset <- object$offset %ORifNULL% rep(0, nobs(object))
	intercept <- as.logical(attr(object$terms,'intercept') %ORifNULL% 0)
	predfun <- function(xt) t(posterior_linpred(object, newdata=data.frame(xt), transform=T, offset=rep(0,nrow(xt))))
	wsample <- rep(1/ndraws, ndraws) # equal sample weights by default
	
	# y and the observation weights in a standard form
	target <- .get_standard_y(unname(get_y(object)), weights(object), fam)
	wobs <- target$weights
	y <- target$y

	# TODO: how to handle cvfits?
	# one idea: add new argument cvfun for init_refmodel, so that user has to specify one of these
	# and then we could also specify that here instead of cvfits...
	init_refmodel(x=x, y=y, family=fam, predfun=predfun, dis=dis, offset=offset, 
								wobs=wobs, wsample=wsample, intercept=intercept, cvfits=NULL) 
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
#' @param ... Currently ignored.
#'
#' @return An object that can be passed to all the functions that
#' take the reference fit as the first argument, such as \link{varsel}, \link{cv_varsel},
#' \link[=proj-pred]{proj_predict} and \link[=proj-pred]{proj_linpred}.

#' @export
init_refmodel <- function(x, y, family, predfun=NULL, dis=NULL, offset=NULL, 
													wobs=NULL, wsample=NULL, intercept=TRUE, cvfits=NULL, ...) {
	
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
		# genuine reference model. add impact of offset to the prediction function
		predmu <- function(x,offset=0) family$linkinv( family$linkfun(predfun(x)) + offset )
		mu <- predmu(x,offset)
		if (NROW(y)!=NROW(mu)) 
			stop(paste0('The number of rows in the output of predfun(x) does not match with the given y;',
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






#' Predictions using the reference model
#'
#' Compute the predictions using the reference model, that is, compute the
#' expected value for the next observation, or evaluate the log-predictive
#' density at a given point.
#'
#' @param object The object of class \code{refmodel}.
#' @param xnew Matrix of predictor values used in the prediction. 
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
predict.refmodel <- function(object, xnew, ynew = NULL, offsetnew = NULL,
														 weightsnew = NULL, type = 'response', ...) {
	
	if ('datafit' %in% class(object))
		stop('Cannot make predictions with data reference only.')
	
	if (is.null(offsetnew)) offsetnew <- rep(0, nrow(xnew))
	if (is.null(weightsnew)) weightsnew <- rep(1, nrow(xnew))
	
	mu <- object$predfun(xnew, offsetnew)
	
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





















