#' Reference model structure
#'
#' Function \code{get_refmodel} is a generic function for creating the reference
#' model structure from a specific \code{object}. The \code{get_refmodel}
#' methods usually call \code{init_refmodel} which is the underlying workhorse
#' to create the reference model structure (and may also be used directly
#' without using \code{get_refmodel}).
#'
#' @name get-refmodel
#'
#' @param object Object from which the reference model is created. For
#'   \code{init_refmodel}, an object on which the functions from arguments
#'   \code{extract_model_data} and \code{ref_predfun} can be applied, with a
#'   \code{NULL} object being treated specially (see section "Value" below). For
#'   \code{get_refmodel.default}, an object on which function \code{family} can
#'   be applied to retrieve the family (if argument \code{family} is
#'   \code{NULL}), additionally to the properties required for
#'   \code{init_refmodel}. For non-default methods of \code{get_refmodel}, an
#'   object of the corresponding class.
#' @param data Data used for fitting the reference model.
#' @param formula Reference model's formula. For general information on formulas
#'   in \R, see \code{\link{formula}}. For multilevel formulas, see also package
#'   \pkg{lme4}, in particular \code{\link[lme4:lmer]{lme4::lmer}} and
#'   \code{\link[lme4:glmer]{lme4::glmer}}.
#' @param ref_predfun Prediction function for the linear predictor of the
#'   reference model. May be \code{NULL} for using an internal default.
#'   Otherwise, needs to have the prototype \code{ref_predfun(fit, newdata =
#'   NULL)} where:
#'   \itemize{
#'     \item{\code{fit} accepts the reference model fit as given in argument
#'     \code{object} (but possibly re-fitted to a subset of the observations, as
#'     done in K-fold cross-validation);}
#'     \item{\code{newdata} accepts either \code{NULL} (for using the original
#'     dataset, typically stored in \code{fit}) or data for new observations (at
#'     least in the form of a \code{data.frame}).}
#'   }
#'   Let \eqn{N} denote the number of observations in the original dataset (used
#'   for fitting the reference model) and \eqn{S} the number of posterior draws
#'   for the reference model's parameters. Then the return value of
#'   \code{ref_predfun} has to be a \eqn{N \times S}{N x S} matrix.
#' @param proj_predfun Prediction function for the linear predictor of a
#'   submodel onto which the reference model is projected. May be \code{NULL}
#'   for using an internal default. Otherwise, needs to have the prototype
#'   \code{proj_predfun(fit, newdata = NULL, weights = NULL)} where:
#'   \itemize{
#'     \item{\code{fit} accepts fit(s) of a submodel as returned by
#'     \code{\link{project}} in its output element \code{sub_fit} (which in turn
#'     is the same as the return value of \code{div_minimizer}, except if
#'     \code{\link{project}} was used with a \code{"vsel"} object from an L1
#'     search as well as \code{cv_search = FALSE});}
#'     \item{\code{newdata} accepts either \code{NULL} (for using the original
#'     dataset, typically stored in \code{fit}) or data for new observations (at
#'     least in the form of a \code{data.frame});}
#'     \item{\code{weights} accepts either \code{NULL} (for using a vector of
#'     ones) or weights for the new observations from \code{newdata} (at least
#'     in the form of a numeric vector).}
#'   }
#'   Let \eqn{N} denote the number of observations and
#'   \eqn{S_{\mbox{prj}}}{S_prj} the number of projected draws (corresponding to
#'   the number of fits in \code{fit}). Then the return value of
#'   \code{proj_predfun} has to be:
#'   \itemize{
#'     \item{a vector or a 1-column matrix of length \eqn{N} if
#'     \eqn{S_{\mbox{prj}} = 1}{S_prj = 1};}
#'     \item{a \eqn{N \times S_{\mbox{prj}}}{N x S_prj} matrix if
#'     \eqn{S_{\mbox{prj}} > 1}{S_prj > 1}.}
#'   }
#' @param div_minimizer A function for minimizing the Kullback-Leibler (KL)
#'   divergence from a submodel to the reference model (i.e., for performing the
#'   projection of the reference model onto a submodel). May be \code{NULL} for
#'   using an internal default. Otherwise, needs to have the prototype
#'   \code{div_minimizer(formula, data, family, weights = NULL, ...)} where
#'   (with \eqn{S_{\mbox{prj}}}{S_prj} denoting the number of resulting
#'   projected draws):
#'   \itemize{
#'     \item{\code{formula} accepts either a standard formula with a single
#'     response (if \eqn{S_{\mbox{prj}} = 1}{S_prj = 1}) or a formula with
#'     \eqn{S_{\mbox{prj}} > 1}{S_prj > 1} response variables
#'     \code{\link{cbind}}-ed on the left-hand side in which case the projection
#'     has to be performed for each of the response variables separately;}
#'     \item{\code{data} accepts a \code{data.frame} to be used for the
#'     projection;}
#'     \item{\code{family} accepts a \code{"family"} object (see argument
#'     \code{family});}
#'     \item{\code{weights} accepts either \code{NULL} (for using a vector of
#'     ones as weights) or weights for the observations from \code{data} (at
#'     least in the form of a numeric vector).}
#'   }
#'   The return value of \code{div_minimizer} has to be:
#'   \itemize{
#'     \item{a fitted model object if \eqn{S_{\mbox{prj}} = 1}{S_prj = 1};}
#'     \item{a list of \eqn{S_{\mbox{prj}}}{S_prj} fitted model objects if
#'     \eqn{S_{\mbox{prj}} > 1}{S_prj > 1}.}
#'   }
#'   This output of \code{div_minimizer} is used, e.g., by \code{proj_predfun}'s
#'   argument \code{fit}.
#' @param extract_model_data A function for fetching some variables (response,
#'   observation weights, offsets) from the original dataset (i.e., the dataset
#'   used for fitting the reference model) or from a new dataset. This function
#'   needs to have the prototype \code{extract_model_data(object, newdata,
#'   wrhs = NULL, orhs = NULL, extract_y = TRUE)}, where:
#'   \itemize{
#'     \item{\code{object} accepts the reference model fit as given in argument
#'     \code{object} (but possibly re-fitted to a subset of the observations, as
#'     done in K-fold cross-validation);}
#'     \item{\code{newdata} accepts data for new observations (at least in the
#'     form of a \code{data.frame});}
#'     \item{\code{wrhs} accepts at least either \code{NULL} (for using a vector
#'     of ones) or a right-hand side formula consisting only of the variable in
#'     \code{newdata} containing the weights;}
#'     \item{\code{orhs} accepts at least either \code{NULL} (for using a vector
#'     of zeros) or a right-hand side formula consisting only of the variable in
#'     \code{newdata} containing the offsets;}
#'     \item{\code{extract_y} accepts a single logical value indicating whether
#'     output element \code{y} (see below) shall be \code{NULL} (\code{TRUE}) or
#'     not (\code{FALSE}).}
#'   }
#'   The return value of \code{extract_model_data} needs to be a \code{list}
#'   with elements \code{y}, \code{weights}, and \code{offset}, each being a
#'   numeric vector containing the data for the response, the observation
#'   weights, and the offsets, respectively. An exception is that \code{y} may
#'   also be \code{NULL} (depending on argument \code{extract_y}).
#' @param family A \code{"family"} object representing the observational model
#'   (i.e., the distributional family for the response). For general information
#'   on \code{"family"} objects in \R, see \code{\link{family}}.
#' @param folds For K-fold cross-validation only. A vector of fold indices for
#'   each observation from \code{data}.
#' @param cvfits For K-fold cross-validation only. A list with one sublist
#'   called \code{"fits"} containing K-fold fitted objects from which reference
#'   models are created. The \code{cvfits} list (i.e., the superlist) needs to
#'   have attributes \code{"K"} and \code{"folds"}: \code{"K"} has to be a
#'   single integer giving the number of folds and \code{"folds"} has to be an
#'   integer vector giving the fold indices (one fold index per observation).
#'   Note that \code{cvfits} takes precedence over \code{cvfun}, i.e., if both
#'   are provided, \code{cvfits} is used.
#' @param cvfun For K-fold cross-validation only. A function that, given a folds
#'   vector, fits a reference model per fold and returns the fitted object. May
#'   be \code{NULL} if \code{object} is \code{NULL}. Note that \code{cvfits}
#'   takes precedence over \code{cvfun}, i.e., if both are provided,
#'   \code{cvfits} is used.
#' @param dis A vector of posterior draws for the dispersion parameter (if such
#'   a parameter exists; else \code{dis} may be \code{NULL}).
#' @param ... For \code{get_refmodel.default} and \code{get_refmodel.stanreg}:
#'   arguments passed to \code{init_refmodel}. For the \code{get_refmodel}
#'   generic: arguments passed to the appropriate method. Else: ignored.
#'
#' @return An object that can be passed to all the functions that take the
#'   reference model fit as the first argument, such as \link{varsel},
#'   \link{cv_varsel}, \link{project}, \link[=proj-pred]{proj_predict}, and
#'   \link[=proj-pred]{proj_linpred}. Usually, the returned object is of class
#'   \code{"refmodel"}. However, if \code{object} is \code{NULL}, the returned
#'   object is of class \code{c("datafit", "refmodel")} which is handled
#'   differently at several places throughout this package. In particular, for a
#'   \code{"datafit"}, argument \code{ref_predfun} is ignored and an internal
#'   function is used instead which always returns \code{NA}.
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
#' @template args-newdata
#' @param object The object of class \code{refmodel}.
#' @param ynew New (test) target variables. If given, then the log predictive
#'   density for the new observations is computed.
#' @param type Scale on which the predictions are returned. Either 'link' (the
#'   latent function value, from -inf to inf) or 'response' (the scale on which
#'   the target \code{y} is measured, obtained by taking the inverse-link from
#'   the latent value).
#' @param ... Currently ignored.
#'
#' @details Argument \code{weightsnew} is only relevant if
#'   \code{!is.null(ynew)}.
#'
#' @return Returns either a vector of predictions, or vector of log predictive
#'   densities evaluated at \code{ynew} if \code{ynew} is not \code{NULL}.

#' @export
predict.refmodel <- function(object, newdata, ynew = NULL, offsetnew = NULL,
                             weightsnew = NULL, type = "response", ...) {
  if (!(type %in% c("response", "link"))) {
    stop("type should be one of ('response', 'link')")
  }
  if (inherits(object, "datafit")) {
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
                                   offsetnew)

  weightsnew <- w_o$weights
  offsetnew <- w_o$offset

  ## ref_predfun returns link(mu)
  mu <- object$ref_predfun(object$fit, newdata) + offsetnew

  if (is.null(ynew)) {
    if (type == "link") {
      pred <- mu
    } else {
      pred <- object$family$linkinv(mu)
    }

    ## integrate over the samples
    if (NCOL(pred) > 1) {
      pred <- rowMeans(pred)
    }

    return(pred)
  } else {

    ## evaluate the log predictive density at the given ynew values
    loglik <- object$family$ll_fun(
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
get_refmodel.default <- function(object, formula, family = NULL, ...) {
  if (is.null(family)) {
    family <- extend_family(family(object))
  } else {
    family <- extend_family(family)
  }

  extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                 orhs = NULL, extract_y = TRUE) {
    resp_form <- if (!extract_y) NULL else lhs(formula)
    args <- nlist(object, newdata, wrhs, orhs, resp_form)
    return(do_call(.extract_model_data, args))
  }

  refmodel <- init_refmodel(
    object = object, formula = formula, family = family,
    extract_model_data = extract_model_data, ...
  )
  return(refmodel)
}

#' @rdname get-refmodel
#' @export
get_refmodel.stanreg <- function(object, ...) {
  family <- family(object)
  family <- extend_family(family)

  if (length(object$offset) > 0 &&
      is.null(attr(terms(object$formula), "offset"))) {
    # In this case, we would have to use argument `offset` of
    # posterior_linpred.stanreg() to allow for new offsets, requiring changes in
    # all ref_predfun() calls. Furthermore, there is rstanarm issue #541. Thus,
    # throw an error:
    stop("It looks like `object` was fitted with offsets specified via ",
         "argument `offset`. Currently, projpred does not support offsets ",
         "specified this way. Please use an `offset()` term in the model ",
         "formula instead.")
  }

  if (inherits(object, "gamm4")) {
    formula <- formula.gamm4(object)
  } else {
    formula <- object$formula
  }

  data <- object$data
  if (length(object$weights) != 0) {
    if ("projpred_internal_wobs_stanreg" %in% names(data)) {
      stop("Need to write to column `projpred_internal_wobs_stanreg` of ",
           "`data`, but that column already exists. Please rename this ",
           "column in `data` and try again.")
    }
    data$projpred_internal_wobs_stanreg <- object$weights
    default_wrhs <- ~ projpred_internal_wobs_stanreg
  } else {
    default_wrhs <- NULL
  }
  if (length(object$offset) != 0) {
    if ("projpred_internal_offs_stanreg" %in% names(data)) {
      stop("Need to write to column `projpred_internal_offs_stanreg` of ",
           "`data`, but that column already exists. Please rename this ",
           "column in `data` and try again.")
    }
    data$projpred_internal_offs_stanreg <- object$offset
    default_orhs <- ~ projpred_internal_offs_stanreg
  } else {
    default_orhs <- NULL
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
    if (family$family != "binomial") {
      stop("This case should not occur. Please notify the package maintainer.")
    }
    # This check is needed to be able to overwrite `default_wrhs` safely:
    if (length(object$weights) != 0) {
      stop("projpred cannot handle observation weights for a binomial family ",
           "with > 1 trials.")
    }
    resp_form <- as.formula(paste("~", response_name[[1]]))
    default_wrhs <- as.formula(paste(
      "~", response_name[[2]], "+",
      response_name[[1]]
    ))
  } else {
    resp_form <- as.formula(paste("~", response_name))
  }

  extract_model_data <- function(object, newdata = NULL, wrhs = default_wrhs,
                                 orhs = default_orhs, extract_y = TRUE) {
    if (!extract_y) {
      resp_form <- NULL
    }

    if (is.null(newdata)) {
      newdata <- data
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

  ref_predfun_stanreg <- function(fit, newdata = NULL) {
    linpred_out <- t(
      posterior_linpred(fit, transform = FALSE, newdata = newdata)
    )
    # Element `stan_function` is not documented in
    # `?rstanarm::`stanreg-objects``, so check at least its length:
    if (length(fit$stan_function) != 1) {
      stop("Unexpected length of `<stanreg_fit>$stan_function`. Please notify ",
           "the package maintainer.")
    }
    # Since posterior_linpred() is supposed to include the offsets in its
    # result, subtract them here and use a workaround for rstanarm issue #541
    # and rstanarm issue #542. This workaround consists of using `cond_no_offs`
    # which indicates whether posterior_linpred() excluded (`TRUE`) or included
    # (`FALSE`) the offsets:
    cond_no_offs <- (
      fit$stan_function %in% c("stan_lmer", "stan_glmer") &&
        !is.null(attr(terms(fit$formula), "offset"))
    ) || (
      fit$stan_function %in% c("stan_lm", "stan_glm") &&
        !is.null(newdata) && length(fit$offset) > 0
    )
    if (!cond_no_offs) {
      offs <- extract_model_data(fit, newdata = newdata)$offset
      stopifnot(identical(nrow(linpred_out), length(offs)))
      linpred_out <- linpred_out - offs
    }
    return(linpred_out)
  }

  cvfun <- function(folds) {
    cvres <- rstanarm::kfold(object,
                             K = max(folds), save_fits = TRUE,
                             folds = folds)
    fits <- cvres$fits[, "fit"]
    return(fits)
  }

  refmodel <- init_refmodel(
    object = object, data = data, formula = formula, family = family,
    ref_predfun = ref_predfun_stanreg, extract_model_data = extract_model_data,
    dis = dis, cvfun = cvfun, ...
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
      linpred_out <- t(
        posterior_linpred(fit, transform = FALSE, newdata = newdata)
      )
      # Since posterior_linpred() is supposed to include the offsets in its
      # result, subtract them here:
      offs <- extract_model_data(fit, newdata = newdata)$offset
      if (length(offs) > 0) {
        stopifnot(length(offs) %in% c(1L, nrow(linpred_out)))
        linpred_out <- linpred_out - offs
      }
      return(linpred_out)
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
    stop("Please provide argument `data`.")
  }
  if (is.null(extract_model_data)) {
    stop("Please provide argument `extract_model_data`.")
  }
  model_data <- extract_model_data(object, newdata = data)
  weights <- model_data$weights
  offset <- model_data$offset
  y <- model_data$y

  data[, response_name] <- y

  target <- .get_standard_y(y, weights, family)
  y <- target$y
  weights <- target$weights

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
                                                 weights = weights) +
                                      offset))
  }

  proper_model <- !is.null(object)

  ## ref_predfun should already take into account the family of the model
  ## we leave this here just in case
  if (proper_model) {
    mu <- ref_predfun(object)
    mu <- unname(as.matrix(mu))
    mu <- family$linkinv(mu)
  } else {
    if (family$family != "binomial") {
      mu <- y
    } else {
      mu <- y / weights
    }
    mu <- matrix(mu)
    ref_predfun_datafit <- function(fit = NULL, newdata = NULL) {
      stopifnot(is.null(fit))
      if (is.null(newdata)) {
        matrix(rep(NA, NROW(y)))
      } else {
        matrix(rep(NA, NROW(newdata)))
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
