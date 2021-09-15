#' Reference model structure
#'
#' Function [get_refmodel()] is a generic function for creating the reference
#' model structure from a specific `object`. The [get_refmodel()] methods
#' usually call [init_refmodel()] which is the underlying workhorse to create
#' the reference model structure (and may also be used directly without using
#' [get_refmodel()]). Some arguments are for K-fold cross-validation (K-fold CV)
#' only; see [cv_varsel()] for the use of K-fold CV in \pkg{projpred}.
#'
#' @name refmodel-init-get
#'
#' @param object Object from which the reference model is created. For
#'   [init_refmodel()], an object on which the functions from arguments
#'   `extract_model_data` and `ref_predfun` can be applied, with a `NULL` object
#'   being treated specially (see section "Value" below). For
#'   [get_refmodel.default()], an object on which function [family()] can be
#'   applied to retrieve the family (if argument `family` is `NULL`),
#'   additionally to the properties required for [init_refmodel()]. For
#'   non-default methods of [get_refmodel()], an object of the corresponding
#'   class.
#' @param data Data used for fitting the reference model.
#' @param formula Reference model's formula. For general information on formulas
#'   in \R, see [`formula`]. For multilevel formulas, see also package
#'   \pkg{lme4}, in particular [lme4::lmer()] and [lme4::glmer()].
#' @param ref_predfun Prediction function for the linear predictor of the
#'   reference model. See section "Details" below.
#' @param proj_predfun Prediction function for the linear predictor of a
#'   submodel onto which the reference model is projected. See section "Details"
#'   below.
#' @param div_minimizer A function for minimizing the Kullback-Leibler (KL)
#'   divergence from a submodel to the reference model (i.e., for performing the
#'   projection of the reference model onto a submodel). The output of
#'   `div_minimizer` is used, e.g., by `proj_predfun`'s argument `fit`. See
#'   section "Details" below.
#' @param extract_model_data A function for fetching some variables (response,
#'   observation weights, offsets) from the original dataset (i.e., the dataset
#'   used for fitting the reference model) or from a new dataset. See section
#'   "Details" below.
#' @param family A [`family`] object representing the observational model (i.e.,
#'   the distributional family for the response). For general information on
#'   [`family`] objects in \R, see [`family`].
#' @param cvfits For K-fold CV only. A list with one sublist called `fits`
#'   containing K-fold fitted objects from which reference models are created.
#'   The `cvfits` list (i.e., the superlist) needs to have attributes `K` and
#'   `folds`: `K` has to be a single integer giving the number of folds and
#'   `folds` has to be an integer vector giving the fold indices (one fold index
#'   per observation). Note that `cvfits` takes precedence over `cvfun`, i.e.,
#'   if both are provided, `cvfits` is used.
#' @param cvfun For K-fold CV only. A function that, given a folds vector, fits
#'   a reference model per fold and returns the fitted object. May be `NULL` if
#'   `object` is `NULL`. Note that `cvfits` takes precedence over `cvfun`, i.e.,
#'   if both are provided, `cvfits` is used.
#' @param dis A vector of posterior draws for the dispersion parameter (if such
#'   a parameter exists; else `dis` may be `NULL`).
#' @param ... For [get_refmodel.default()] and [get_refmodel.stanreg()]:
#'   arguments passed to [init_refmodel()]. For the [get_refmodel()] generic:
#'   arguments passed to the appropriate method. Else: ignored.
#'
#' @details
#'
#' # `ref_predfun`, `proj_predfun`, `div_minimizer`
#'
#' Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer` may be `NULL`
#' for using an internal default. Otherwise, let \eqn{N} denote the number of
#' observations in the original dataset (used for fitting the reference model),
#' \eqn{S} the number of posterior draws for the reference model's parameters,
#' and \eqn{S_{\mbox{prj}}}{S_prj} the number of resulting projected draws. Then
#' the functions supplied to these arguments need to have the following
#' prototypes:
#' * `ref_predfun(fit, newdata = NULL)` where:
#'     + `fit` accepts the reference model fit as given in argument `object`
#'     (but possibly re-fitted to a subset of the observations, as done in
#'     K-fold CV);
#'     + `newdata` accepts either `NULL` (for using the original dataset,
#'     typically stored in `fit`) or data for new observations (at least in the
#'     form of a `data.frame`).
#' * `proj_predfun(fit, newdata = NULL)` where:
#'     + `fit` accepts a list of length \eqn{S_{\mbox{prj}}}{S_prj} containing
#'     this number of submodel fits. This list is the same as that returned by
#'     [project()] in its output element `sub_fit` (which in turn is the same as
#'     the return value of `div_minimizer`, except if [project()] was used with
#'     an object of class `vsel` based on an L1 search as well as `cv_search =
#'     FALSE`);
#'     + `newdata` accepts either `NULL` (for using the original dataset,
#'     typically stored in `fit`) or data for new observations (at least in the
#'     form of a `data.frame`);
#' * `div_minimizer` does not need to have a specific prototype, but it needs to
#' be able to be called with the following arguments:
#'     + `formula` accepts either a standard [`formula`] with a single response
#'     (if \eqn{S_{\mbox{prj}} = 1}{S_prj = 1}) or a [`formula`] with
#'     \eqn{S_{\mbox{prj}} > 1}{S_prj > 1} response variables [cbind()]-ed on
#'     the left-hand side in which case the projection has to be performed for
#'     each of the response variables separately;
#'     + `data` accepts a `data.frame` to be used for the projection;
#'     + `family` accepts a [`family`] object (see argument `family`);
#'     + `weights` accepts either `NULL` (for using a vector of ones as weights)
#'     or observation weights (at least in the form of a numeric vector);
#'     + `projpred_var` accepts a numeric vector of length \eqn{N} containing
#'     predictive variances (necessary for \pkg{projpred}'s internal (G)LM
#'     fitter);
#'     + `projpred_regul` accepts a single numeric value as supplied to argument
#'     `regul` of [project()], for example.
#'
#' The return value of those functions needs to be:
#' * `ref_predfun`: a \eqn{N \times S}{N x S} matrix.
#' * `proj_predfun`: a \eqn{N \times S_{\mbox{prj}}}{N x S_prj} matrix.
#' * `div_minimizer`: a `list` of length \eqn{S_{\mbox{prj}}}{S_prj} containing
#' this number of submodel fits.
#'
#' # `extract_model_data`
#'
#' The function supplied to argument `extract_model_data` needs to have the
#' prototype `extract_model_data(object, newdata, wrhs = NULL, orhs = NULL,
#' extract_y = TRUE)`, where:
#' * `object` accepts the reference model fit as given in argument `object` (but
#' possibly re-fitted to a subset of the observations, as done in K-fold CV);
#' * `newdata` accepts data for new observations (at least in the form of a
#' `data.frame`);
#' * `wrhs` accepts at least either `NULL` (for using a vector of ones) or a
#' right-hand side formula consisting only of the variable in `newdata`
#' containing the weights;
#' * `orhs` accepts at least either `NULL` (for using a vector of zeros) or a
#' right-hand side formula consisting only of the variable in `newdata`
#' containing the offsets;
#' * `extract_y` accepts a single logical value indicating whether output
#' element `y` (see below) shall be `NULL` (`TRUE`) or not (`FALSE`).
#'
#' The return value of `extract_model_data` needs to be a `list` with elements
#' `y`, `weights`, and `offset`, each being a numeric vector containing the data
#' for the response, the observation weights, and the offsets, respectively. An
#' exception is that `y` may also be `NULL` (depending on argument `extract_y`).
#'
#' @return An object that can be passed to all the functions that take the
#'   reference model fit as the first argument, such as [varsel()],
#'   [cv_varsel()], [project()], [proj_linpred()], and [proj_predict()].
#'   Usually, the returned object is of class `refmodel`. However, if `object`
#'   is `NULL`, the returned object is of class `c("datafit", "refmodel")` which
#'   is handled differently at several places throughout this package. In
#'   particular, for a `datafit`, argument `ref_predfun` is ignored and an
#'   internal function is used instead which always returns `NA`.
#'
#' @examples
#' if (requireNamespace("rstanarm", quietly = TRUE)) {
#'   # Data:
#'   dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)
#'
#'   # The "stanreg" fit which will be used as the reference model:
#'   fit <- rstanarm::stan_glm(
#'     y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
#'     QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
#'   )
#'
#'   # Define the reference model formally:
#'   ref <- get_refmodel(fit)
#'   print(class(ref)) # gives `"refmodel"`
#'   # Now see, for example, `?varsel.refmodel`, `?cv_varsel.refmodel`,
#'   # `?project`, and `?predict.refmodel` for possible post-processing
#'   # functions. Most of them call get_refmodel() internally at the beginning
#'   # (or are methods whose generics do so), so you will rarely need to call
#'   # get_refmodel() yourself.
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
#' @param object The object of class `refmodel`.
#' @param ynew New (test) target variables. If given, then the log predictive
#'   density for the new observations is computed.
#' @param type Scale on which the predictions are returned. Either 'link' (the
#'   latent function value, from -inf to inf) or 'response' (the scale on which
#'   the target `y` is measured, obtained by taking the inverse-link from the
#'   latent value).
#' @param ... Currently ignored.
#'
#' @details Argument `weightsnew` is only relevant if `!is.null(ynew)`.
#'
#' @return Returns either a vector of predictions, or vector of log predictive
#'   densities evaluated at `ynew` if `ynew` is not `NULL`.
#'
#' @export
predict.refmodel <- function(object, newdata = NULL, ynew = NULL,
                             offsetnew = NULL, weightsnew = NULL,
                             type = "response", ...) {
  if (!type %in% c("response", "link")) {
    stop("type should be one of ('response', 'link')")
  }
  if (inherits(object, "datafit")) {
    stop("Cannot make predictions for an `object` of class \"datafit\".")
  }
  if (!is.null(ynew) && (!is.numeric(ynew) || NCOL(ynew) != 1)) {
    stop("Argument `ynew` must be a numeric vector.")
  }

  w_o <- object$extract_model_data(object$fit, newdata = newdata,
                                   wrhs = weightsnew, orhs = offsetnew)
  weightsnew <- w_o$weights
  offsetnew <- w_o$offset
  if (length(weightsnew) == 0) {
    weightsnew <- rep(1, length(w_o$y))
  }
  if (length(offsetnew) == 0) {
    offsetnew <- rep(0, length(w_o$y))
  }

  ## ref_predfun returns link(mu)
  mu <- object$ref_predfun(object$fit, newdata) + offsetnew

  if (is.null(ynew)) {
    pred <- if (type == "link") mu else object$family$linkinv(mu)
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

fetch_data <- function(data, obs = NULL, newdata = NULL) {
  if (is.null(obs)) {
    if (is.null(newdata)) {
      data_out <- data
    } else {
      data_out <- newdata
    }
  } else if (is.null(newdata)) {
    data_out <- data[obs, , drop = FALSE]
  } else {
    data_out <- newdata[obs, , drop = FALSE]
  }
  return(as.data.frame(data_out))
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

#' @rdname refmodel-init-get
#' @export
get_refmodel <- function(object, ...) {
  UseMethod("get_refmodel", object)
}

#' @rdname refmodel-init-get
#' @export
get_refmodel.refmodel <- function(object, ...) {
  # If the object is already of class "refmodel", then simply return it as is:
  object
}

#' @rdname refmodel-init-get
#' @export
get_refmodel.vsel <- function(object, ...) {
  # The reference model is stored in the `object` of class "vsel":
  object$refmodel
}

#' @rdname refmodel-init-get
#' @export
get_refmodel.default <- function(object, formula, family = NULL, ...) {
  if (is.null(family)) {
    family <- family(object)
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

#' @rdname refmodel-init-get
#' @export
get_refmodel.stanreg <- function(object, ...) {
  # Family ------------------------------------------------------------------

  family <- family(object)

  # Data --------------------------------------------------------------------

  data <- object$data
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`object$data` must be a `data.frame` or a `matrix` (but a ",
         "`data.frame` is recommended).")
  }

  # Weights (for the observations):
  if (family$family == "binomial" && length(object$weights) > 0) {
    stop("In case of the binomial family, projpred cannot handle observation ",
         "weights (apart from the numbers of trials).")
  }
  if (length(object$weights) > 0) {
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

  # Offsets:
  # Element `stan_function` (needed for handling the offsets) is not documented
  # in `?rstanarm::`stanreg-objects``, so check at least its length and type:
  if (length(object$stan_function) != 1 ||
      !is.vector(object$stan_function, mode = "character")) {
    stop("Unexpected value of `object$stan_function`. Please notify the ",
         "package maintainer.")
  }
  if (length(object$offset) > 0) {
    if (is.null(attr(terms(formula(object)), "offset"))) {
      # In this case, we would have to use argument `offset` of
      # posterior_linpred.stanreg() to allow for new offsets, requiring changes
      # in all ref_predfun() calls. Furthermore, there is rstanarm issue #541.
      # Thus, throw an error:
      stop("It looks like `object` was fitted with offsets specified via ",
           "argument `offset`. Currently, projpred does not support offsets ",
           "specified this way. Please use an `offset()` term in the model ",
           "formula instead.")
    }
    if (object$stan_function == "stan_gamm4") {
      stop("Because of rstanarm issue #546 (see GitHub), projpred cannot ",
           "allow offsets for additive models (fit with ",
           "rstanarm::stan_gamm4()).")
    }
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

  # Formula -----------------------------------------------------------------

  if (inherits(object, "gamm4")) {
    formula <- formula.gamm4(object)
  } else {
    formula <- formula(object)
  }

  stopifnot(inherits(formula, "formula"))
  formula <- expand_formula(formula, data)
  response_name <- extract_terms_response(formula)$response
  if (length(response_name) == 2) {
    if (family$family != "binomial") {
      stop("For non-binomial families, a two-column response is not allowed.")
    }
    default_wrhs <- as.formula(paste(
      "~", response_name[2], "+", response_name[1]
    ))
    response_name <- response_name[1]
  } else if (length(response_name) > 2) {
    stop("The response is not allowed to have more than two columns.")
  }
  resp_form <- as.formula(paste("~", response_name))
  formula <- update(formula, as.formula(paste(response_name, "~ .")))

  # Functions ---------------------------------------------------------------

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

  ref_predfun <- function(fit, newdata = NULL) {
    linpred_out <- t(
      posterior_linpred(fit, transform = FALSE, newdata = newdata)
    )
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
      # Observation weights are not needed here, so use `wrhs = NULL` to avoid
      # potential conflicts for a non-`NULL` default `wrhs`:
      offs <- extract_model_data(fit, newdata = newdata, wrhs = NULL)$offset
      stopifnot(identical(nrow(linpred_out), length(offs)))
      linpred_out <- linpred_out - offs
    }
    return(linpred_out)
  }

  cvfun <- function(folds, ...) {
    cvres <- rstanarm::kfold(object, K = max(folds), save_fits = TRUE,
                             folds = folds, ...)
    fits <- cvres$fits[, "fit"]
    return(fits)
  }

  # Miscellaneous -----------------------------------------------------------

  if (.has_dispersion(family)) {
    dis <- data.frame(object)[, "sigma"]
  } else {
    dis <- NULL
  }

  # Output ------------------------------------------------------------------

  return(init_refmodel(
    object = object, data = data, formula = formula, family = family,
    ref_predfun = ref_predfun, extract_model_data = extract_model_data,
    dis = dis, cvfun = cvfun, ...
  ))
}

#' @rdname refmodel-init-get
#' @importFrom rstantools posterior_linpred
#' @export
init_refmodel <- function(object, data, formula, family, ref_predfun = NULL,
                          div_minimizer = NULL, proj_predfun = NULL,
                          extract_model_data, cvfun = NULL,
                          cvfits = NULL, dis = NULL, ...) {
  proper_model <- !is.null(object)

  # Formula -----------------------------------------------------------------

  stopifnot(inherits(formula, "formula"))
  formula <- expand_formula(formula, data)
  response_name <- extract_terms_response(formula)$response
  if (length(response_name) == 2) {
    if (family$family != "binomial") {
      stop("For non-binomial families, a two-column response is not allowed.")
    }
  } else if (length(response_name) > 2) {
    stop("The response is not allowed to have more than two columns.")
  }
  # Remove parentheses from the response:
  response_name <- gsub("[()]", "", response_name)
  formula <- update(formula, paste(response_name[1], "~ ."))

  # Data --------------------------------------------------------------------

  model_data <- extract_model_data(object, newdata = data)
  weights <- model_data$weights
  offset <- model_data$offset
  y <- model_data$y

  # Add (transformed) response under the (possibly) new name:
  data[, response_name] <- y

  target <- .get_standard_y(y, weights, family)
  y <- target$y
  weights <- target$weights

  if (family$family == "binomial") {
    if (!all(.is.wholenumber(y))) {
      stop("In projpred, the response must contain numbers of successes (not ",
           "proportions of successes), in contrast to glm() where this is ",
           "the convention for a 1-column response.")
    } else if (all(y %in% c(0, 1)) &&
               length(response_name) == 1 &&
               !all(weights == 1)) {
      warning(
        "Assuming that the response contains numbers of successes (not ",
        "proportions of successes), in contrast to glm()."
      )
    }
  }

  if (is.null(offset)) {
    offset <- rep(0, NROW(y))
  }

  # Functions ---------------------------------------------------------------

  if (proper_model && is.null(ref_predfun)) {
    ref_predfun <- function(fit, newdata = NULL) {
      linpred_out <- t(
        posterior_linpred(fit, transform = FALSE, newdata = newdata)
      )
      # Since posterior_linpred() is supposed to include the offsets in its
      # result, subtract them here:
      # Observation weights are not needed here, so use `wrhs = NULL` to avoid
      # potential conflicts for a non-`NULL` default `wrhs`:
      offs <- extract_model_data(fit, newdata = newdata, wrhs = NULL)$offset
      if (length(offs) > 0) {
        stopifnot(length(offs) %in% c(1L, nrow(linpred_out)))
        linpred_out <- linpred_out - offs
      }
      return(linpred_out)
    }
  } else if (!proper_model) {
    if (!is.null(ref_predfun)) {
      warning("Ignoring argument `ref_predfun` because `object` is `NULL`.")
    }
    ref_predfun <- function(fit, newdata = NULL) {
      stopifnot(is.null(fit))
      if (is.null(newdata)) {
        return(matrix(rep(NA, NROW(y))))
      } else {
        return(matrix(rep(NA, NROW(newdata))))
      }
    }
  }

  if (is.null(div_minimizer)) {
    div_minimizer <- divmin
  }

  if (is.null(proj_predfun)) {
    proj_predfun <- subprd
  }

  fetch_data_wrapper <- function(obs = NULL) {
    fetch_data(data, obs, newdata = NULL)
  }

  # Family ------------------------------------------------------------------

  if (!.has_family_extras(family)) {
    family <- extend_family(family)
  }

  family$mu_fun <- function(fit, obs = NULL, newdata = NULL, offset = NULL) {
    newdata <- fetch_data(data, obs = obs, newdata = newdata)
    if (is.null(offset)) {
      offset <- rep(0, nrow(newdata))
    }
    family$linkinv(proj_predfun(fit, newdata = newdata) + offset)
  }

  # mu ----------------------------------------------------------------------

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
  }

  # Miscellaneous -----------------------------------------------------------

  ndraws <- ncol(mu)
  if (is.null(dis)) {
    dis <- rep(0, ndraws)
  }

  if (proper_model) {
    loglik <- t(family$ll_fun(mu, dis, y, weights = weights))
  } else {
    loglik <- NULL
  }

  if (is.null(cvfun)) {
    if (!proper_model) {
      # This is a dummy definition for cvfun(), but it will lead to standard CV
      # for `datafit`s; see cv_varsel() and .get_kfold():
      cvfun <- function(folds, ...) {
        lapply(seq_len(max(folds)), function(k) list())
      }
    } else if (is.null(cvfits)) {
      stop("Please provide either argument `cvfun` or argument `cvfits`.")
    }
  }

  # Equal sample (draws) weights by default:
  wsample <- rep(1 / ndraws, ndraws)

  intercept <- as.logical(attr(terms(formula), "intercept"))

  # Output ------------------------------------------------------------------

  refmodel <- nlist(
    fit = object, formula, div_minimizer, family, mu, dis, y, loglik, intercept,
    proj_predfun, fetch_data = fetch_data_wrapper, wobs = weights, wsample,
    offset, cvfun, cvfits, extract_model_data, ref_predfun
  )
  if (proper_model) {
    class(refmodel) <- "refmodel"
  } else {
    class(refmodel) <- c("datafit", "refmodel")
  }

  return(refmodel)
}
