
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
## Usage with stanreg objects
#' dat <- data.frame(y = rnorm(100), x = rnorm(100))
#' fit <- stan_glm(y ~ x, family = gaussian(), data = dat)
#' ref <- get_refmodel(fit)
#' print(class(ref))
#'
## variable selection, use the already constructed reference model
#' vs <- varsel(ref)
## this will first construct the reference model and then execute
## exactly the same way as the previous command (the result is identical)
#' vs <- varsel(fit)
#' }
#'
NULL

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

  ## predfun returns link(mu)
  mu <- object$predfun(object$fit, znew)

  if (is.null(ynew)) {

    if (type == 'link')
      pred <- mu
    else
      pred <- object$family$linkinv(mu)

    ## integrate over the samples
    if (NCOL(pred) > 1)
      pred <- rowMeans(pred)

    return(pred)

  } else {

    ## evaluate the log predictive density at the given ynew values
    loglik <- object$fam$ll_fun(object$family$linkinv(mu), object$dis, ynew,
                                weightsnew)
    S <- ncol(loglik)
    lpd <- apply(loglik, 1, log_sum_exp) - log(S)
    return(lpd)
  }

}

#' @export
get_refmodel_poc <- function(fit, ...) {
  UseMethod("get_refmodel_poc", fit)
}

#' @export
get_refmodel_poc.refmodel <- function(object, ...) {
  ## if the object is reference model already, then simply return it as is
  object
}

#' @export
get_refmodel_poc.vsel <- function(object, ...) {
  ## the reference model is stored in vsel-object
  object$refmodel
}

#' @export
get_refmodel_poc.cvsel <- function(object, ...) {
  ## the reference model is stored in cvsel object
  object$refmodel
}

#' @export
get_refmodel_poc.default <- function(fit, data, y, formula, predfun, proj_predfun,
                                     mle, fetch_data, family=NULL, wobs=NULL,
                                     folds=NULL, cvfits=NULL, penalized=FALSE) {
  fetch_data_wrapper <- function(obs=folds, newdata=NULL)
    fetch_data(data, obs, newdata)

  if (is.null(family))
    family <- extend_family(family(fit))
  else
    family <- extend_family(family)

  family$mu_fun <- function(fit, obs=folds, xnew=NULL) {
    newdata <- fetch_data_wrapper(obs=obs, newdata=xnew)
    family$linkinv(proj_predfun(fit, newdata=newdata))
  }

  refmodel <- init_refmodel_poc(fit, data, y, formula, family, predfun, mle,
                                proj_predfun, fetch_data, fetch_data_wrapper,
                                penalized=penalized)
  refmodel$folds <- folds
  return(refmodel)
}

#' export
get_refmodel_poc.brmsfit <- function(fit, data=NULL, y=NULL, formula=NULL,
                                     predfun=NULL, proj_predfun=NULL, mle=NULL,
                                     folds=NULL, ...) {
  family_name <- family(fit)$family
  fam <- ifelse(family_name == "bernoulli", "binomial", family_name)
  fam <- get(fam, mode = "function")()
  brms_family <- family(fit)
  brms_family$mu.eta <- fam$mu.eta
  brms_family$variance <- fam$variance
  brms_family$dev.resids <- fam$dev.resids
  family <- extend_family(brms_family)

  brms_formula <- formula(fit)
  if (is.null(formula))
    formula <- brms_formula$formula

  terms <- extract_terms_response(formula)
  response_name <- brms_formula$resp

  formula <- update(formula, paste(response_name, " ~ ."))
  p <- parse_bf(brms_formula)

  if (is.null(data))
    data <- fit$data
  if (is.null(y))
    y <- fit$data[, response_name]

  if ("trials" %in% family$ad &&
      inherits(p$adforms$trials, "formula")) {
    trials_form <- p$adforms$trials
    trials_var <- eval(trials_form[[2]], data,
                       environment(trials_form))$vars$trials
    weights <- data[, trials_var]
    y <- y / weights
  } else {
    weights <- NULL
  }

  ## TODO: return y, weights and offset as right hand side formulas
  refmodel <- init_refmodel_poc(fit, data, y, formula, family, predfun, mle,
                                proj_predfun, folds, weights=weights)
  return(refmodel)
}

#' @export
get_refmodel_poc.stanreg <- function(fit, data=NULL, y=NULL, formula=NULL,
                                     predfun=NULL, proj_predfun=NULL, mle=NULL,
                                     folds=NULL, penalized=FALSE, ...) {
  family <- family(fit)
  family <- extend_family(family)

  if (is.null(formula))
    formula <- formula(fit)
  terms <- extract_terms_response(formula)
  response_name <- terms$response

  if (is.null(data))
    data <- fit$data
  if (is.null(y))
    y <- fit$data[, colnames(fit$data) == response_name]

  refmodel <- init_refmodel_poc(fit, data, y, formula, family, predfun, mle,
                                proj_predfun, folds, penalized)
  return(refmodel)
}

#' @export
init_refmodel_poc <- function(fit, data, y, formula, family, predfun=NULL, mle=NULL,
                              proj_predfun=NULL, folds=NULL, penalized=FALSE,
                              weights=NULL, offset=NULL, cvfun=NULL,
                              cvfits=NULL) {
  terms <- extract_terms_response(formula)
  if (is.null(predfun))
    predfun <- function(fit, newdata=NULL)
      t(posterior_linpred(fit, transform = FALSE, newdata = newdata))

  if (is.null(mle) && is.null(proj_predfun))
    if (length(terms$group_terms) != 0) {
      mle <- linear_multilevel_mle
      proj_predfun <- linear_multilevel_proj_predfun
    } else {
      if (!penalized) {
        mle <- linear_mle
        proj_predfun <- linear_proj_predfun
      } else {
        mle <- penalized_linear_mle
        proj_predfun <- penalized_linear_proj_predfun
      }
    }
  else if (is.null(mle) && !is.null(proj_predfun))
    if (length(terms$group_terms) != 0)
      mle <- linear_multilevel_mle
    else
      if (!penalized)
        mle <- linear_mle
    else
      mle <- penalized_linear_mle
  else if (!is.null(mle) && is.null(proj_predfun))
    if (length(terms$group_terms) != 0)
      proj_predfun <- linear_multilevel_proj_predfun
    else
      if (!penalized)
        proj_predfun <- linear_proj_predfun
    else
      proj_predfun <- penalized_linear_proj_predfun

  fetch_data_wrapper <- function(obs=folds, newdata=NULL)
    as.data.frame(fetch_data(data, obs, newdata))

  if (!.has_family_extras (family))
    family <- extend_family(family)

  ## TODO: ideally remove this, have to think about it
  family$mu_fun <- function(fit, obs=folds, xnew=NULL) {
    newdata <- fetch_data_wrapper(obs = obs, newdata = xnew)
    family$linkinv(proj_predfun(fit, newdata = newdata))
  }

  proper_model <- !is.null(fit)

  ## predfun should already take into account the family of the model
  ## we leave this here just in case
  if (proper_model) {
    mu <- predfun(fit)
    mu <- unname(as.matrix(mu))
    mu <- family$linkinv(mu)
  } else
    mu <- matrix(y, NROW(y), 1)

  ndraws <- ncol(mu)

  if (proper_model) {
    ## TODO: eventually this will be a function provided by the user
    dis <- rep(0, ndraws)
    tryCatch({
        dis <- as.data.frame(fit)$sigma %ORifNULL% rep(0, ndraws)
      },
      error = function(e) e
    )
  }
  target <- .get_standard_y(y, weights, family)
  y <- target$y

  ## equal sample weights by default
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  if (proper_model)
    loglik <- t(family$ll_fun(mu, dis, y, weights = weights))
  else
    loglik <- NULL

  if (!proper_model) {
    # this is a dummy definition for cvfun, but it will lead to standard cross-validation
    # for datafit reference; see cv_varsel and get_kfold
    cvfun <- function(folds) lapply(1:max(folds), function(k) list())
  }

  wsample <- rep(1 / ndraws, ndraws) # equal sample weights by default
  if (is.null(offset))
    offset <- rep(0, NROW(y))

  intercept <- as.logical(attr(terms(formula), 'intercept'))
  refmodel <- list(fit=fit, formula=formula, predfun=predfun, mle=mle,
                   family=family, mu=mu, dis=dis, y=y, loglik=loglik,
                   intercept=intercept, proj_predfun=proj_predfun,
                   fetch_data=fetch_data_wrapper, wobs=weights,
                   wsample=wsample, offset=offset, folds=folds,
                   cvfun=cvfun, cvfits=cvfits)
  if (proper_model)
    class(refmodel) <- "refmodel"
  else
    class(refmodel) <- c("datafit", "refmodel")

  return(refmodel)
}
