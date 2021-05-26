#' Get reference model structure
#'
#' Generic function that can be used to create and fetch the reference model
#' structure for all those objects that have this method. All these
#' implementations are wrappers to the \code{\link{init_refmodel}}-function so
#' the returned object has the same type.
#'
#' @name get-refmodel
#'
#' @param object Object from which the reference model is created. See possible
#'   types below.
#' @param data Data on which the reference model was fitted.
#' @param formula Reference model's lme4-like formula.
#' @param ref_predfun Prediction function for the linear predictor of the
#'   reference model.
#' @param proj_predfun Prediction function for the linear predictor of the
#'   projections.
#' @param div_minimizer Maximum likelihood estimator for the underlying
#'   projection.
#' @param fetch_data Wrapper function for fetching the data without directly
#'   accessing it. It should have a prototype fetch_data(data, data_points,
#'   newdata = NULL), where data_points is a vector of data indices and newdata,
#'   if not NULL, is a data frame with new data for testing.
#' @param extract_model_data A function with prototype
#'   extract_model_data(object, newdata, wrhs, orhs), where object is a
#'   reference model fit, newdata is either NULL or a data frame with new
#'   observations, wrhs is a right hand side formula to recover the weights from
#'   the data frame and orhs is a right hand side formula to recover the offset
#'   from the data frame.
#' @param family A family object that represents the observation model for the
#'   reference model.
#' @param wobs A weights vector for the observations in the data. The default is
#'   a vector of ones.
#' @param folds Only used for K-fold variable selection. It is a vector of fold
#'   indices for each data point in data.
#' @param cvfits Only used for K-fold variable selection. A list with one
#'   sublist called \code{"fits"} containing K-fold fitted objects from which
#'   reference models are created. The \code{cvfits} list (i.e., the superlist)
#'   needs to have attributes \code{"K"} and \code{"folds"}: \code{"K"} has to
#'   be a single integer giving the number of folds and \code{"folds"} has to be
#'   an integer vector giving the fold indices (one fold index per observation).
#'   Note that \code{cvfits} takes precedence over \code{cvfun}, i.e., if both
#'   are provided, \code{cvfits} is used.
#' @param cvfun Only used for K-fold variable selection. A function that, given
#'   a folds vector, fits a reference model per fold and returns the fitted
#'   object. Note that \code{cvfits} takes precedence over \code{cvfun}, i.e.,
#'   if both are provided, \code{cvfits} is used.
#' @param offset A vector of offsets per observation to add to the linear
#'   predictor.
#' @param dis A dispersion vector for each observation.
#' @param ... Arguments passed to the methods.
#'
#' @return An object of class \code{"refmodel"} that can be passed to all the
#'   functions that take the reference model fit as the first argument, such as
#'   \link{varsel}, \link{cv_varsel}, \link{project},
#'   \link[=proj-pred]{proj_predict}, and \link[=proj-pred]{proj_linpred}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace('rstanarm', quietly=TRUE)) {
#'   ### Usage with stanreg objects
#'   dat <- data.frame(y = rnorm(100), x = rnorm(100))
#'   fit <- rstanarm::stan_glm(y ~ x, family = gaussian(), data = dat)
#'   ref <- get_refmodel(fit)
#'   print(class(ref))
#'
#'   # variable selection, use the already constructed reference model
#'   vs <- varsel(ref)
#'   # this will first construct the reference model and then execute
#'   # exactly the same way as the previous command (the result is identical)
#'   vs <- varsel(fit)
#' }
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
#' @param newdata Matrix of predictor values used in the prediction.
#' @param ynew New (test) target variables. If given, then the log predictive
#'   density for the new observations is computed.
#' @param offsetnew Offsets for the new observations. By default a vector of
#'   zeros. By default we take the weights from newdata as in the original
#'   model. Either NULL or right hand side formulas.
#' @param weightsnew Weights for the new observations. For binomial model,
#'   corresponds to the number trials per observation. Has effect only if
#'   \code{ynew} is specified. By default a vector of ones. By default we take
#'   the weights from newdata as in the original model. Either NULL or right
#'   hand side formulas.
#' @param type Scale on which the predictions are returned. Either 'link' (the
#'   latent function value, from -inf to inf) or 'response' (the scale on which
#'   the target \code{y} is measured, obtained by taking the inverse-link from
#'   the latent value).
#' @param ... Currently ignored.
#'
#' @return Returns either a vector of predictions, or vector of log predictive
#'   densities evaluated at \code{ynew} if \code{ynew} is not \code{NULL}.

#' @export
predict.refmodel <- function(object, newdata, ynew = NULL, offsetnew = NULL,
                             weightsnew = NULL, type = "response", ...) {
  if (!(type %in% c("response", "link"))) {
    stop("type should be one of ('response', 'link')")
  }
  if ("datafit" %in% class(object)) {
    stop("Cannot make predictions with data reference only.")
  }
  if (!is.null(ynew)) {
    if (!(inherits(ynew, "numeric")) || NCOL(ynew) != 1) {
      stop("ynew must be a numerical vector")
    }
  }

  if (!is.null(offsetnew) && !inherits(offsetnew, "formula")) {
    stop("offsetnew specified but it's not a right hand side formula")
  }

  if (!is.null(weightsnew) && !inherits(weightsnew, "formula")) {
    stop("weightsnew specified but it's not a right hand side formula")
  }

  w_o <- object$extract_model_data(object$fit,
                                   newdata = newdata, weightsnew,
                                   offsetnew
  )

  weightsnew <- w_o$weights
  offsetnew <- w_o$offset

  ## ref_predfun returns link(mu)
  mu <- object$ref_predfun(object$fit, newdata)

  if (is.null(ynew)) {
    if (type == "link") {
      pred <- mu
    } else {
      pred <- object$family$linkinv(mu + offsetnew)
    }

    ## integrate over the samples
    if (NCOL(pred) > 1) {
      pred <- rowMeans(pred)
    }

    return(pred)
  } else {

    ## evaluate the log predictive density at the given ynew values
    loglik <- object$fam$ll_fun(
      object$family$linkinv(mu), object$dis, ynew,
      weightsnew
    )
    S <- ncol(loglik)
    lpd <- apply(loglik, 1, log_sum_exp) - log(S)
    return(lpd)
  }
}

.extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                orhs = NULL, resp_form = NULL) {
  if (is.null(newdata)) {
    newdata <- object$data
  }

  if (inherits(wrhs, "formula")) {
    weights <- eval_rhs(wrhs, newdata)
  } else if (is.null(wrhs)) {
    weights <- rep(1, NROW(newdata))
  } else {
    weights <- wrhs
  }

  if (inherits(orhs, "formula")) {
    offset <- eval_rhs(orhs, newdata)
  } else if (is.null(orhs)) {
    offset <- rep(0, NROW(newdata))
  } else {
    offset <- orhs
  }

  if (inherits(resp_form, "formula")) {
    y <- eval_rhs(resp_form, newdata)
  } else {
    y <- NULL
  }

  return(nlist(y, weights, offset))
}

#' @rdname get-refmodel
#' @export
get_refmodel <- function(object, ...) {
  UseMethod("get_refmodel", object)
}

#' @rdname get-refmodel
#' @export
get_refmodel.refmodel <- function(object, ...) {
  ## if the object is reference model already, then simply return it as is
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
get_refmodel.default <- function(object, data, formula, ref_predfun,
                                 proj_predfun, div_minimizer, fetch_data,
                                 family = NULL, wobs = NULL, folds = NULL,
                                 cvfits = NULL, offset = NULL, cvfun = NULL,
                                 dis = NULL, ...) {
  fetch_data_wrapper <- function(obs = folds, newdata = NULL) {
    fetch_data(data, obs, newdata)
  }

  if (is.null(family)) {
    family <- extend_family(family(object))
  } else {
    family <- extend_family(family)
  }

  extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                 orhs = NULL) {
    resp_form <- lhs(formula)
    args <- nlist(object, newdata, wrhs, orhs, resp_form)
    return(do_call(.extract_model_data, args))
  }

  refmodel <- init_refmodel(object, data, formula, family, ref_predfun,
                            div_minimizer, proj_predfun,
                            extract_model_data = extract_model_data,
                            cvfits = cvfits, folds = folds, cvfun = cvfun,
                            dis = dis
  )
  return(refmodel)
}

#' @rdname get-refmodel
#' @export
get_refmodel.stanreg <- function(object, data = NULL, ref_predfun = NULL,
                                 proj_predfun = NULL, div_minimizer = NULL,
                                 folds = NULL, ...) {
  family <- family(object)
  family <- extend_family(family)
  if (inherits(object, "gamm4")) {
    formula <- formula.gamm4(object)
  } else {
    formula <- object$formula
  }

  if (is.null(data)) {
    data <- object$data
  }

  stopifnot(inherits(formula, "formula"))
  formula <- expand_formula(formula, data)
  terms <- extract_terms_response(formula)
  response_name <- terms$response

  formula <- update(
    formula,
    as.formula(paste(response_name, "~ ."))
  )

  if (length(response_name) > 1) {
    resp_form <- as.formula(paste("~", response_name[[1]]))
    default_wrhs <- as.formula(paste(
      "~", response_name[[2]], "+",
      response_name[[1]]
    ))
  } else {
    resp_form <- as.formula(paste("~", response_name))
    default_wrhs <- NULL
  }

  extract_model_data <- function(object, newdata = NULL, wrhs = default_wrhs,
                                 orhs = NULL, extract_y = TRUE) {
    if (!extract_y) {
      resp_form <- NULL
    }

    if (is.null(newdata)) {
      newdata <- object$data
    }

    if (is.null(wrhs) && !is.null(object) &&
        !is.null(object$weights) && length(object$weights) != 0) {
      wrhs <- ~weights
      newdata <- cbind(newdata, weights = object$weights)
    }

    if (is.null(orhs) && !is.null(object) &&
        !is.null(object$offset) && length(object$offset) != 0) {
      orhs <- ~offset
      newdata <- cbind(newdata, offset = object$offset)
    }

    args <- nlist(object, newdata, wrhs, orhs, resp_form)
    return(do_call(.extract_model_data, args))
  }

  if (length(response_name) > 1) {
    response_name <- response_name[[1]]
  }

  if (.has_dispersion(family)) {
    dis <- data.frame(object)[, "sigma"]
  } else {
    dis <- NULL
  }

  cvfun <- function(folds) {
    cvres <- rstanarm::kfold(object,
                             K = max(folds), save_fits = TRUE,
                             folds = folds
    )
    fits <- cvres$fits[, "fit"]
    return(fits)
  }

  refmodel <- init_refmodel(
    object, data, formula, family,
    ref_predfun = ref_predfun, div_minimizer = div_minimizer,
    proj_predfun = proj_predfun, folds = folds,
    extract_model_data = extract_model_data, dis = dis,
    cvfun = cvfun, ...
  )
  return(refmodel)
}

#' @rdname get-refmodel
#' @importFrom rstantools posterior_linpred
#' @export
init_refmodel <- function(object, data, formula, family, ref_predfun = NULL,
                          div_minimizer = NULL, proj_predfun = NULL,
                          folds = NULL, extract_model_data = NULL, cvfun = NULL,
                          cvfits = NULL, dis = NULL, ...) {
  stopifnot(inherits(formula, "formula"))
  formula <- expand_formula(formula, data)
  terms <- extract_terms_response(formula)
  response_name <- terms$response
  if (is.null(ref_predfun)) {
    ref_predfun <- function(fit, newdata = NULL) {
      t(posterior_linpred(fit, transform = FALSE, newdata = newdata))
    }
  }

  ## remove parens from response
  response_name <- gsub("[()]", "", response_name)
  formula <- update(
    formula,
    paste(response_name, "~ .")
  )

  ## add (transformed) response with new name
  if (is.null(data)) {
    stop("Data was not provided.")
  }
  model_data <- extract_model_data(object, newdata = data)
  weights <- model_data$weights
  offset <- model_data$offset
  y <- model_data$y

  data[, response_name] <- y

  if (is.null(div_minimizer)) {
    if (length(terms$additive_terms) != 0) {
      div_minimizer <- additive_mle
    } else if (length(terms$group_terms) != 0) {
      div_minimizer <- linear_multilevel_mle
    } else {
      div_minimizer <- linear_mle
    }
  }

  if (is.null(proj_predfun)) {
    if (length(terms$additive_terms) != 0) {
      proj_predfun <- additive_proj_predfun
    } else if (length(terms$group_terms) != 0) {
      proj_predfun <- linear_multilevel_proj_predfun
    } else {
      proj_predfun <- linear_proj_predfun
    }
  }

  fetch_data_wrapper <- function(obs = folds, newdata = NULL) {
    as.data.frame(fetch_data(data, obs, newdata))
  }

  if (!.has_family_extras(family)) {
    family <- extend_family(family)
  }

  family$mu_fun <- function(fit, obs = folds, newdata = NULL, offset = NULL,
                            weights = NULL) {
    if (is.null(offset)) {
      offset <- rep(0, length(obs))
    }
    if (is.null(weights)) {
      weights <- rep(1, length(obs))
    }
    newdata <- fetch_data_wrapper(obs = obs, newdata = newdata)
    suppressWarnings(family$linkinv(proj_predfun(fit,
                                                 newdata = newdata,
                                                 weights = weights
    ) + offset))
  }

  proper_model <- !is.null(object)

  ## ref_predfun should already take into account the family of the model
  ## we leave this here just in case
  if (proper_model) {
    mu <- ref_predfun(object)
    mu <- unname(as.matrix(mu))
    mu <- family$linkinv(mu)
  } else {
    mu <- matrix(y / weights, NROW(y), 1)
    ref_predfun_datafit <- function(fit = NULL, newdata = NULL, offset = 0) {
      if (is.null(fit)) {
        if (is.null(newdata)) {
          matrix(rep(NA, NROW(y)))
        } else {
          matrix(rep(NA, NROW(newdata)))
        }
      } else {
        family$linkinv(ref_predfun(fit, newdata))
      }
    }
  }

  ndraws <- ncol(mu)
  if (is.null(dis)) {
    dis <- rep(0, ndraws)
  }

  if (is.null(offset)) {
    offset <- rep(0, NROW(y))
  }

  if (is.null(weights)) {
    weights <- rep(1, NROW(y))
  }

  target <- .get_standard_y(y, weights, family)
  y <- target$y

  if (proper_model) {
    loglik <- t(family$ll_fun(mu, dis, y, weights = weights))
  } else {
    loglik <- NULL
  }

  # this is a dummy definition for cvfun, but it will lead to standard
  # cross-validation for datafit reference; see cv_varsel and get_kfold
  if (is.null(cvfun)) {
    if (!proper_model) {
      cvfun <- function(folds, ...) lapply(1:max(folds), function(k) list())
    } else if (is.null(cvfits)) {
      stop("Please provide either 'cvfun' or 'cvfits'.")
    }
  }

  wsample <- rep(1 / ndraws, ndraws) # equal sample weights by default
  intercept <- as.logical(attr(terms(formula), "intercept"))
  refmodel <- nlist(
    fit = object, formula, div_minimizer, family, mu, dis, y,
    loglik, intercept, proj_predfun, fetch_data = fetch_data_wrapper,
    wobs = weights, wsample, offset, folds, cvfun, cvfits, extract_model_data
  )
  if (proper_model) {
    refmodel$ref_predfun <- ref_predfun
    class(refmodel) <- "refmodel"
  } else {
    refmodel$ref_predfun <- ref_predfun_datafit
    class(refmodel) <- c("datafit", "refmodel")
  }

  return(refmodel)
}
