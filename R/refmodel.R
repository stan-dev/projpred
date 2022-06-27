#' Reference model structure
#'
#' Function [get_refmodel()] is a generic function for creating the reference
#' model structure from a specific `object`. The methods for [get_refmodel()]
#' usually call [init_refmodel()] which is the underlying workhorse (and may
#' also be used directly without a call to [get_refmodel()]). Some arguments are
#' for \eqn{K}-fold cross-validation (\eqn{K}-fold CV) only; see [cv_varsel()]
#' for the use of \eqn{K}-fold CV in \pkg{projpred}.
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
#' @param data Data used for fitting the reference model. Any `contrasts`
#'   attributes of the dataset's columns are silently removed.
#' @param formula Reference model's formula. For general information on formulas
#'   in \R, see [`formula`]. For multilevel formulas, see also package
#'   \pkg{lme4} (in particular, functions [lme4::lmer()] and [lme4::glmer()]).
#'   For additive formulas, see also packages \pkg{mgcv} (in particular,
#'   function [mgcv::gam()]) and \pkg{gamm4} (in particular, function
#'   [gamm4::gamm4()]) as well as the notes in section "Formula terms" below.
#' @param ref_predfun Prediction function for the linear predictor of the
#'   reference model, including offsets (if existing). See also section
#'   "Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`" below. If
#'   `object` is `NULL`, `ref_predfun` is ignored and an internal default is
#'   used instead.
#' @param proj_predfun Prediction function for the linear predictor of a
#'   submodel onto which the reference model is projected. See also section
#'   "Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`" below.
#' @param div_minimizer A function for minimizing the Kullback-Leibler (KL)
#'   divergence from a submodel to the reference model (i.e., for performing the
#'   projection of the reference model onto a submodel). The output of
#'   `div_minimizer` is used, e.g., by `proj_predfun`'s argument `fits`. See
#'   also section "Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`"
#'   below.
#' @param extract_model_data A function for fetching some variables (response,
#'   observation weights, offsets) from the original dataset (i.e., the dataset
#'   used for fitting the reference model) or from a new dataset. See also
#'   section "Argument `extract_model_data`" below.
#' @param family A [`family`] object representing the observational model (i.e.,
#'   the distributional family for the response). May be `NULL` for
#'   [get_refmodel.default()] in which case the family is retrieved from
#'   `object`.
#' @param cvfits For \eqn{K}-fold CV only. A `list` containing a sub-`list`
#'   called `fits` containing the \eqn{K} model fits from which reference model
#'   structures are created. The `cvfits` `list` (i.e., the super-`list`) needs
#'   to have attributes `K` and `folds`: `K` has to be a single integer giving
#'   the number of folds and `folds` has to be an integer vector giving the fold
#'   indices (one fold index per observation). Each element of `cvfits$fits`
#'   (i.e., each of the \eqn{K} model fits) needs to be a list. Only one of
#'   `cvfits` and `cvfun` needs to be provided (for \eqn{K}-fold CV). Note that
#'   `cvfits` takes precedence over `cvfun`, i.e., if both are provided,
#'   `cvfits` is used.
#' @param cvfun For \eqn{K}-fold CV only. A function that, given a fold indices
#'   vector, fits the reference model separately for each fold and returns the
#'   \eqn{K} model fits as a `list`. Each of the \eqn{K} model fits needs to be
#'   a `list`. If `object` is `NULL`, `cvfun` may be `NULL` for using an
#'   internal default. Only one of `cvfits` and `cvfun` needs to be provided
#'   (for \eqn{K}-fold CV). Note that `cvfits` takes precedence over `cvfun`,
#'   i.e., if both are provided, `cvfits` is used.
#' @param cvrefbuilder For \eqn{K}-fold CV only. A function that, given a
#'   reference model fit for fold \eqn{k \in \{1, ..., K\}}{k = 1, ..., K} (this
#'   model fit is the \eqn{k}-th element of the return value of `cvfun` or the
#'   \eqn{k}-th element of `cvfits$fits`, extended by elements `omitted`
#'   (containing the indices of the left-out observations in that fold) and
#'   `projpred_k` (containing the integer \eqn{k})), returns an object of the
#'   same type as [init_refmodel()] does. Argument `cvrefbuilder` may be `NULL`
#'   for using an internal default: [get_refmodel()] if `object` is not `NULL`
#'   and a function calling [init_refmodel()] appropriately (with the assumption
#'   `dis = 0`) if `object` is `NULL`.
#' @param dis A vector of posterior draws for the dispersion parameter (if
#'   existing). May be `NULL` if the model has no dispersion parameter or if the
#'   model does have a dispersion parameter, but `object` is `NULL` (in which
#'   case `0` is used for `dis`). Note that for the [gaussian()] `family`, `dis`
#'   is the standard deviation, not the variance.
#' @param latent_proj A single logical value indicating whether to use the
#'   latent projection (`TRUE`) or not (`FALSE`).
#' @param ... For [get_refmodel.default()] and [get_refmodel.stanreg()]:
#'   arguments passed to [init_refmodel()]. For the [get_refmodel()] generic:
#'   arguments passed to the appropriate method. For [init_refmodel()]:
#'   arguments passed to [extend_family()].
#'
#' @details
#'
#' # Formula terms
#'
#' For additive models (still an experimental feature), only [mgcv::s()] and
#' [mgcv::t2()] are currently supported as smooth terms. Furthermore, these need
#' to be called without any arguments apart from the predictor names (symbols).
#' For example, for smoothing the effect of a predictor `x`, only `s(x)` or
#' `t2(x)` are allowed. As another example, for smoothing the joint effect of
#' two predictors `x` and `z`, only `s(x, z)` or `t2(x, z)` are allowed (and
#' analogously for higher-order joint effects, e.g., of three predictors).
#'
#' # Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`
#'
#' Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer` may be `NULL`
#' for using an internal default. Otherwise, let \eqn{N} denote the number of
#' observations (in case of CV, these may be reduced to each fold),
#' \eqn{S_{\mbox{ref}}}{S_ref} the number of posterior draws for the reference
#' model's parameters, and \eqn{S_{\mbox{prj}}}{S_prj} the number of (possibly
#' clustered) parameter draws for projection (short: the number of projected
#' draws). For the augmented-data projection, let \eqn{C_{\mbox{cat}}}{C_cat}
#' denote the number of response categories, \eqn{C_{\mbox{lat}}}{C_lat} the
#' number of latent response categories (which typically equals
#' \eqn{C_{\mbox{cat}} - 1}{C_cat - 1}), and define \eqn{N_{\mbox{augcat}} := N
#' \cdot C_{\mbox{cat}}}{N_augcat := N * C_cat} as well as
#' \eqn{N_{\mbox{auglat}} := N \cdot C_{\mbox{lat}}}{N_auglat := N * C_lat}.
#' Then the functions supplied to these arguments need to have the following
#' prototypes:
#' * `ref_predfun`: `ref_predfun(fit, newdata = NULL)` where:
#'     + `fit` accepts the reference model fit as given in argument `object`
#'     (but possibly re-fitted to a subset of the observations, as done in
#'     \eqn{K}-fold CV).
#'     + `newdata` accepts either `NULL` (for using the original dataset,
#'     typically stored in `fit`) or data for new observations (at least in the
#'     form of a `data.frame`).
#' * `proj_predfun`: `proj_predfun(fits, newdata)` where:
#'     + `fits` accepts a `list` of length \eqn{S_{\mbox{prj}}}{S_prj}
#'     containing this number of submodel fits. This `list` is the same as that
#'     returned by [project()] in its output element `submodl` (which in turn is
#'     the same as the return value of `div_minimizer`, except if [project()]
#'     was used with an `object` of class `vsel` based on an L1 search as well
#'     as with `refit_prj = FALSE`).
#'     + `newdata` accepts data for new observations (at least in the form of a
#'     `data.frame`).
#' * `div_minimizer` does not need to have a specific prototype, but it needs to
#' be able to be called with the following arguments:
#'     + `formula` accepts either a standard [`formula`] with a single response
#'     (if \eqn{S_{\mbox{prj}} = 1}{S_prj = 1} or in case of the augmented-data
#'     projection) or a [`formula`] with \eqn{S_{\mbox{prj}} > 1}{S_prj > 1}
#'     response variables [cbind()]-ed on the left-hand side in which case the
#'     projection has to be performed for each of the response variables
#'     separately.
#'     + `data` accepts a `data.frame` to be used for the projection. In case of
#'     the traditional (non-augmented-data) projection, this dataset has \eqn{N}
#'     rows. In case of the augmented-data projection, this dataset has
#'     \eqn{N_{\mbox{augcat}}}{N_augcat} rows.
#'     + `family` accepts a [`family`] object.
#'     + `weights` accepts either observation weights (at least in the form of a
#'     numeric vector) or `NULL` (for using a vector of ones as weights).
#'     + `projpred_var` accepts an \eqn{N \times S_{\mbox{prj}}}{N x S_prj}
#'     matrix of predictive variances (necessary for \pkg{projpred}'s internal
#'     GLM fitter) in case of the traditional projection and an
#'     \eqn{N_{\mbox{augcat}} \times S_{\mbox{prj}}}{N_augcat x S_prj} matrix of
#'     `NA`s in case of the augmented-data projection.
#'     + `projpred_regul` accepts a single numeric value as supplied to argument
#'     `regul` of [project()], for example.
#'     + `projpred_ws_aug` accepts an \eqn{N \times S_{\mbox{prj}}}{N x S_prj}
#'     matrix of expected values for the response in case of the traditional
#'     projection and an \eqn{N_{\mbox{augcat}} \times S_{\mbox{prj}}}{N_augcat
#'     x S_prj} matrix of probabilities for the response categories in case of
#'     the augmented-data projection.
#'     + `...` accepts further arguments specified by the user.
#'
#' The return value of these functions needs to be:
#' * `ref_predfun`: for the traditional projection, an \eqn{N \times
#' S_{\mbox{ref}}}{N x S_ref} matrix; for the augmented-data projection, an
#' \eqn{S_{\mbox{ref}} \times N \times C_{\mbox{lat}}}{S_ref x N x C_lat} array
#' (the only exception is the augmented-data projection for the [binomial()]
#' family in which case `ref_predfun` needs to return an \eqn{N \times
#' S_{\mbox{ref}}}{N x S_ref} matrix just like for the traditional projection
#' because the array is constructed by an internal wrapper function).
#' * `proj_predfun`: for the traditional projection, an \eqn{N \times
#' S_{\mbox{prj}}}{N x S_prj} matrix; for the augmented-data projection, an
#' \eqn{N \times C_{\mbox{lat}} \times S_{\mbox{prj}}}{N x C_lat x S_prj} array.
#' * `div_minimizer`: a `list` of length \eqn{S_{\mbox{prj}}}{S_prj} containing
#' this number of submodel fits.
#'
#' # Argument `extract_model_data`
#'
#' The function supplied to argument `extract_model_data` needs to have the
#' prototype
#' ```{r, eval = FALSE}
#' extract_model_data(object, newdata, wrhs = NULL, orhs = NULL, extract_y = TRUE)
#' ```
#' where:
#' * `object` accepts the reference model fit as given in argument `object` (but
#' possibly re-fitted to a subset of the observations, as done in \eqn{K}-fold
#' CV).
#' * `newdata` accepts either `NULL` (for using the original dataset, typically
#' stored in `object`) or data for new observations (at least in the form of a
#' `data.frame`).
#' * `wrhs` accepts at least either `NULL` (for using a vector of ones) or a
#' right-hand side formula consisting only of the variable in `newdata`
#' containing the weights.
#' * `orhs` accepts at least either `NULL` (for using a vector of zeros) or a
#' right-hand side formula consisting only of the variable in `newdata`
#' containing the offsets.
#' * `extract_y` accepts a single logical value indicating whether output
#' element `y` (see below) shall be `NULL` (`TRUE`) or not (`FALSE`).
#'
#' The return value of `extract_model_data` needs to be a `list` with elements
#' `y`, `weights`, and `offset`, each being a numeric vector containing the data
#' for the response, the observation weights, and the offsets, respectively. An
#' exception is that `y` may also be `NULL` (depending on argument `extract_y`).
#'
#' # Augmented-data projection
#'
#' If a custom reference model for an augmented-data projection is needed, see
#' also [extend_family()].
#'
#' For the augmented-data projection, the response vector resulting from
#' `extract_model_data` is internally coerced to a `factor` (using
#' [as.factor()]). The levels of this `factor` have to be identical to
#' `family$cats` (see [extend_family()]'s argument `augdat_y_unqs`).
#'
#' Note that response-specific offsets (i.e., one length-\eqn{N} offset vector
#' per response category) are not supported by \pkg{projpred} yet. So far, only
#' offsets which are the same across all response categories are supported.
#'
#' # Latent projection
#'
#' In case of the latent projection (see argument `latent_proj`), the
#' [rstantools::log_lik()] generic is applied to `object` (inside of
#' [init_refmodel()]), so there needs to be an appropriate method (which does
#' exist in packages \pkg{rstanarm} and \pkg{brms} for `stanreg`s and
#' `brmsfit`s, respectively).
#'
#' @return An object that can be passed to all the functions that take the
#'   reference model fit as the first argument, such as [varsel()],
#'   [cv_varsel()], [project()], [proj_linpred()], and [proj_predict()].
#'   Usually, the returned object is of class `refmodel`. However, if `object`
#'   is `NULL`, the returned object is of class `datafit` as well as of class
#'   `refmodel` (with `datafit` being first). Objects of class `datafit` are
#'   handled differently at several places throughout this package.
#'
#' @examples
#' if (requireNamespace("rstanarm", quietly = TRUE)) {
#'   # Data:
#'   dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)
#'
#'   # The "stanreg" fit which will be used as the reference model (with small
#'   # values for `chains` and `iter`, but only for technical reasons in this
#'   # example; this is not recommended in general):
#'   fit <- rstanarm::stan_glm(
#'     y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
#'     QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
#'   )
#'
#'   # Define the reference model explicitly:
#'   ref <- get_refmodel(fit)
#'   print(class(ref)) # gives `"refmodel"`
#'   # Now see, for example, `?varsel`, `?cv_varsel`, and `?project` for
#'   # possible post-processing functions. Most of the post-processing functions
#'   # call get_refmodel() internally at the beginning, so you will rarely need
#'   # to call get_refmodel() yourself.
#'
#'   # A custom reference model which may be used in a variable selection where
#'   # the candidate predictors are not a subset of those used for the reference
#'   # model's predictions:
#'   ref_cust <- init_refmodel(
#'     fit,
#'     data = dat_gauss,
#'     formula = y ~ X6 + X7,
#'     family = gaussian(),
#'     extract_model_data = function(object, newdata = NULL, wrhs = NULL,
#'                                   orhs = NULL, extract_y = TRUE) {
#'       if (!extract_y) {
#'         resp_form <- NULL
#'       } else {
#'         resp_form <- ~ y
#'       }
#'
#'       if (is.null(newdata)) {
#'         newdata <- dat_gauss
#'       }
#'
#'       args <- projpred:::nlist(object, newdata, wrhs, orhs, resp_form)
#'       return(projpred::do_call(projpred:::.extract_model_data, args))
#'     },
#'     cvfun = function(folds) {
#'       kfold(
#'         fit, K = max(folds), save_fits = TRUE, folds = folds, cores = 1
#'       )$fits[, "fit"]
#'     },
#'     dis = as.matrix(fit)[, "sigma"]
#'   )
#'   # Now, the post-processing functions mentioned above (for example,
#'   # varsel(), cv_varsel(), and project()) may be applied to `ref_cust`.
#' }
#'
NULL

#' Predictions or log predictive densities from a reference model
#'
#' This is the [predict()] method for `refmodel` objects (returned by
#' [get_refmodel()] or [init_refmodel()]). It offers three types of output which
#' are all based on the reference model and new (or old) observations: Either
#' the linear predictor on link scale, the linear predictor transformed to
#' response scale, or the log predictive density.
#'
#' @template args-newdata
#' @param object An object of class `refmodel` (returned by [get_refmodel()] or
#'   [init_refmodel()]).
#' @param ynew If not `NULL`, then this needs to be a vector of new (or old)
#'   response values. See also section "Value" below. In case of the
#'   augmented-data projection, `ynew` is internally coerced to a `factor`
#'   (using [as.factor()]). The levels of this `factor` have to be a subset of
#'   `object$family$cats` (see [extend_family()]'s argument `augdat_y_unqs`).
#' @param type Only relevant if `is.null(ynew)`. The scale on which the
#'   predictions are returned, either `"link"` or `"response"` (see
#'   [predict.glm()] but note that [predict.refmodel()] does not adhere to the
#'   typical \R convention of a default prediction on link scale). For both
#'   scales, the predictions are averaged across the posterior draws.
#' @param ... Currently ignored.
#'
#' @details Argument `weightsnew` is only relevant if `!is.null(ynew)`.
#'
#' @return In the following, \eqn{N}, \eqn{C_{\mbox{cat}}}{C_cat}, and
#'   \eqn{C_{\mbox{lat}}}{C_lat} from help topic [refmodel-init-get] are used.
#'   Furthermore, let \eqn{C} denote either \eqn{C_{\mbox{cat}}}{C_cat} (if
#'   `type = "response"`) or \eqn{C_{\mbox{lat}}}{C_lat} (if `type = "link"`).
#'   Then, if `is.null(ynew)`, the returned object contains the reference
#'   model's predictions (with the scale depending on argument `type`) as a
#'   length-\eqn{N} vector in case of the traditional projection and as an
#'   \eqn{N \times C}{N x C} matrix in case of the augmented-data projection. If
#'   `!is.null(ynew)`, the returned object is a length-\eqn{N} vector of log
#'   predictive densities evaluated at `ynew`.
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
  if (!is.null(ynew) && (!is.numeric(ynew) || NCOL(ynew) != 1) &&
      !object$family$for_augdat) {
    stop("Argument `ynew` must be a numeric vector.")
  }
  if (!is.null(ynew) && object$family$for_augdat) {
    ynew <- as.factor(ynew)
    if (!all(levels(ynew) %in% object$family$cats)) {
      stop("The levels of the response variable (after coercing it to a ",
           "`factor`) have to be a subset of `family$cats`. Either modify ",
           "`ynew` accordingly or see the documentation for extend_family()'s ",
           "argument `augdat_y_unqs` to solve this.")
    }
    # Re-assign the original levels because some levels might be missing:
    ynew <- factor(ynew, levels = object$family$cats)
  }

  if (is.null(newdata)) {
    newdata <- object$fetch_data()
  } else {
    newdata <- na.fail(newdata)
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
  if (object$family$for_augdat && !all(weightsnew == 1)) {
    stop("Currently, the augmented-data projection may not be combined with ",
         "observation weights (other than 1).")
  }
  if (inherits(object$fit, "stanreg") && length(object$offset) > 0) {
    if ("projpred_internal_offs_stanreg" %in% names(newdata)) {
      stop("Need to write to column `projpred_internal_offs_stanreg` of ",
           "`newdata`, but that column already exists. Please rename this ",
           "column in `newdata` and try again.")
    }
    newdata$projpred_internal_offs_stanreg <- offsetnew
  }

  ## ref_predfun returns eta = link(mu)
  eta <- object$ref_predfun(object$fit, newdata = newdata) + offsetnew

  if (is.null(ynew)) {
    pred <- if (type == "link") eta else object$family$linkinv(eta)
    ## integrate over the samples
    if (NCOL(pred) > 1) {
      pred <- rowMeans(pred)
    }
    if (object$family$for_augdat) {
      pred <- structure(pred,
                        nobs_orig = nrow(newdata),
                        class = "augvec")
      pred <- augmat2arr(augvec2augmat(pred))
      pred <- matrix(pred, nrow = dim(pred)[1], ncol = dim(pred)[2])
    }
    return(pred)
  } else {
    ## evaluate the log predictive density at the given ynew values
    loglik <- object$family$ll_fun(
      object$family$linkinv(eta), object$dis, ynew, weightsnew
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

refprd <- function(fit, newdata = NULL) {
  # For safety reasons, keep `transform = FALSE` even though this should
  # be the default in all posterior_linpred() methods (but we cannot be
  # sure with regard to user-defined posterior_linpred() methods):
  linpred_out <- posterior_linpred(
    fit, transform = FALSE, newdata = newdata
  )
  if (length(dim(linpred_out)) == 2) {
    linpred_out <- t(linpred_out)
  } else if (length(dim(linpred_out)) != 3) {
    # A 3-dimensional array would be ok for the augmented-data projection
    # (and doesn't need any transposition or permutation of its
    # dimensions). Everything else is unexpected.
    stop("Unexpected structure for posterior_linpred()'s output. Please ",
         "notify the package maintainer.")
  }
  return(linpred_out)
}

.extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                orhs = NULL, resp_form = NULL) {
  if (is.null(newdata)) {
    newdata <- object$data
  }

  if (inherits(wrhs, "formula")) {
    weights <- eval_rhs(wrhs, newdata)
  } else if (is.null(wrhs)) {
    weights <- rep(1, nrow(newdata))
  } else {
    weights <- wrhs
  }

  if (inherits(orhs, "formula")) {
    offset <- eval_rhs(orhs, newdata)
  } else if (is.null(orhs)) {
    offset <- rep(0, nrow(newdata))
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
  UseMethod("get_refmodel")
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
get_refmodel.stanreg <- function(object, latent_proj = FALSE, ...) {
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    stop("Please install the 'rstanarm' package.")
  }

  # Family ------------------------------------------------------------------

  family <- family(object)
  if (object$stan_function == "stan_polr") {
    # Create a custom family (in particular, to have `family$family`):
    if (family == "logistic") {
      family <- "logit"
    } else if (family == "loglog") {
      stop("Currently, the \"", family, "\" link is not supported by ",
           "projpred.")
    }
    family <- structure(list(family = "cumulative_rstanarm",
                             link = family,
                             cats = levels(object$y)),
                        class = "family")
  }
  aug_data <- object$stan_function == "stan_polr" && !latent_proj
  if (aug_data) {
    # Currently, we need brms for the special link and inverse link function.
    # It shouldn't be hard to implement these separately so that brms is not
    # needed here, but that would introduce redundancies and for now, relying
    # on brms (>= 2.16.3) is the quickest solution and not too demanding.
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package \"brms\" needed. Please install it.",
           call. = FALSE)
    }
    stopifnot(utils::packageVersion("brms") >= "2.16.3")
  }

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
  if (length(object$offset) > 0) {
    # Element `stan_function` (needed here for handling rstanarm issue #546) is
    # not documented in `?rstanarm::`stanreg-objects``, so check at least its
    # length and type:
    if (length(object$stan_function) != 1 ||
        !is.vector(object$stan_function, mode = "character")) {
      stop("Unexpected value of `object$stan_function`. Please notify the ",
           "package maintainer.")
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
    # The easiest way to deal with rstanarm issue #541 and rstanarm issue #542,
    # changes between rstanarm versions 2.21.2 and 2.21.3 with respect to these
    # issues, and the fact that offsets may be specified via argument `offset`
    # of the respective model-fitting function (e.g., rstanarm::stan_glm()) is
    # to include offsets explicitly in the call to
    # rstanarm:::posterior_linpred.stanreg().

    # Observation weights are not needed here, so use `wrhs = NULL` to avoid
    # potential conflicts for a non-`NULL` default `wrhs`:
    offs <- extract_model_data(fit, newdata = newdata, wrhs = NULL)$offset
    n_obs <- nrow(newdata %||% data)
    if (length(offs) == 0) {
      offs <- rep(0, n_obs)
    } else if (length(offs) == 1) {
      offs <- rep(offs, n_obs)
    } else if (length(offs) != n_obs) {
      stop("Unexpected length of element `offset` returned by ",
           "extract_model_data() (see `?init_refmodel`).")
    }
    linpred_out <- posterior_linpred(fit, newdata = newdata, offset = offs)
    stopifnot(length(dim(linpred_out)) == 2)
    aug_data <- fit$stan_function == "stan_polr" && !latent_proj
    if (aug_data) {
      # Since rstanarm::posterior_linpred.stanreg() doesn't offer an argument
      # like `incl_thres` of brms::posterior_linpred.brmsfit(), we need to
      # incorporate the thresholds into the linear predictors by hand:
      linpred_out <- apply(
        as.matrix(fit)[, names(fit$zeta), drop = FALSE],
        2,
        function(x) {
          x - linpred_out
        },
        simplify = FALSE
      )
      linpred_out <- abind::abind(linpred_out, rev.along = 0)
    } else {
      linpred_out <- t(linpred_out)
    }
    return(linpred_out)
  }

  cvfun <- function(folds) {
    # Use `cores = 1` because of rstanarm issue #551. In fact, this issue only
    # affects Windows systems, but since `cores = 1` leads to an *inner*
    # parallelization (i.e., across chains, not across CV folds) with
    # `stan_cores <- getOption("mc.cores", 1)` cores, this should also be
    # suitable for other systems:
    kfold(
      object, K = max(folds), save_fits = TRUE, folds = folds, cores = 1
    )$fits[, "fit"]
  }

  cvrefbuilder <- function(cvfit) {
    get_refmodel(cvfit, latent_proj = latent_proj, ...)
  }

  # Miscellaneous -----------------------------------------------------------

  if (.has_dispersion(family)) {
    dis <- data.frame(object)[, "sigma"]
  } else {
    dis <- NULL
  }

  # Augmented-data projection -----------------------------------------------

  if (aug_data) {
    args_augdat <- list(
      augdat_link = "brms" %:::% "link_cumulative",
      augdat_ilink = "brms" %:::% "inv_link_cumulative",
      augdat_args_link = list(link = family$link),
      augdat_args_ilink = list(link = family$link)
    )
  } else {
    args_augdat <- list()
  }

  # Output ------------------------------------------------------------------

  args_basic <- list(
    object = object, data = data, formula = formula, family = family,
    ref_predfun = ref_predfun, extract_model_data = extract_model_data,
    dis = dis, cvfun = cvfun, cvrefbuilder = cvrefbuilder,
    latent_proj = latent_proj
  )
  return(do.call(init_refmodel, args = c(args_basic, args_augdat, list(...))))
}

#' @rdname refmodel-init-get
#' @export
init_refmodel <- function(object, data, formula, family, ref_predfun = NULL,
                          div_minimizer = NULL, proj_predfun = NULL,
                          extract_model_data, cvfun = NULL,
                          cvfits = NULL, dis = NULL, cvrefbuilder = NULL,
                          latent_proj = FALSE, ...) {
  # Family ------------------------------------------------------------------

  if (latent_proj) {
    family <- extend_family(gaussian())
  }

  if (family$family == "Student_t") {
    warning("Support for the `Student_t` family is still experimental.")
  } else if (family$family == "Gamma") {
    warning("Support for the `Gamma` family is still experimental.")
  }

  if (!.has_family_extras(family)) {
    family <- extend_family(family, ...)
  }
  aug_data <- family$for_augdat
  if (aug_data &&
      isTRUE(getOption("projpred.warn_augdat_experimental", TRUE))) {
    warning("The augmented-data projection is still experimental.")
  }

  family$mu_fun <- function(fits, obs = NULL, newdata = NULL, offset = NULL) {
    newdata <- fetch_data(data, obs = obs, newdata = newdata)
    if (is.null(offset)) {
      offset <- rep(0, nrow(newdata))
    } else {
      stopifnot(length(offset) %in% c(1L, nrow(newdata)))
    }
    family$linkinv(proj_predfun(fits, newdata = newdata) + offset)
  }

  if (family$family == "categorical" && family$link != "logit") {
    stop("For the brms::categorical() family, projpred only supports the ",
         "logit link.")
  }

  # Special case: `datafit` -------------------------------------------------

  proper_model <- !is.null(object)
  if (!proper_model && aug_data) {
    stop("Currently, the augmented-data projection may not be combined with ",
         "`object = NULL` (i.e., a `datafit`).")
  }

  # Formula -----------------------------------------------------------------

  stopifnot(inherits(formula, "formula"))
  data <- na.fail(data)
  formula <- expand_formula(formula, data)
  response_name <- extract_terms_response(formula)$response
  if (length(response_name) == 2) {
    if (family$family != "binomial") {
      stop("For non-binomial families, a two-column response is not allowed.")
    } else if (aug_data) {
      stop("Currently, the augmented-data projection may not be combined with ",
           "a 2-column response.")
    }
  } else if (length(response_name) > 2) {
    stop("The response is not allowed to have more than two columns.")
  }
  # Remove parentheses from the response:
  response_name <- gsub("[()]", "", response_name)
  if (latent_proj) {
    response_name <- paste0(".", response_name[1])
  }
  formula <- update(formula, paste(response_name[1], "~ ."))
  if (formula_contains_additive_terms(formula)) {
    if (aug_data) {
      stop("Currently, the augmented-data projection may not be combined with ",
           "additive models.")
    } else if (isTRUE(getOption("projpred.warn_additive_experimental", TRUE))) {
      warning("Support for additive models is still experimental.")
    }
  }

  # Functions ---------------------------------------------------------------

  if (proper_model) {
    if (is.null(ref_predfun)) {
      ref_predfun <- refprd
    }
    if (aug_data && family$family == "binomial") {
      ref_predfun_mat <- ref_predfun
      ref_predfun <- function(fit, newdata = NULL) {
        linpred1 <- ref_predfun_mat(fit = fit, newdata = newdata)
        linpred1 <- t(linpred1)
        return(array(linpred1, dim = c(dim(linpred1), 1L)))
      }
    }
    # Since posterior_linpred() is supposed to include any offsets but (at least
    # currently) projpred expects the final ref_predfun() to exclude any offsets
    # (see issue #186), the offsets have to be subtracted here by a wrapper
    # function. This wrapper function also performs some preparations for the
    # augmented-data projection:
    ref_predfun_usr <- ref_predfun
    ref_predfun <- function(fit, newdata = NULL) {
      linpred_out <- ref_predfun_usr(fit = fit, newdata = newdata)
      if (length(dim(linpred_out)) == 2) {
        n_obs <- nrow(linpred_out)
      } else if (length(dim(linpred_out)) == 3) {
        # For the augmented-data projection, `linpred_out` is expected to be a
        # 3-dimensional array with dimensions S_ref x N x C_lat (see
        # `?init_refmodel` for a definition of these dimensions). Therefore, it
        # is converted to an augmented-rows matrix (see file "augdat.R" for a
        # definition):
        linpred_out <- arr2augmat(linpred_out, margin_draws = 1)
        n_obs <- attr(linpred_out, "nobs_orig")
      } else {
        stop("Unexpected structure for `linpred_out`. Does the return value ",
             "of `ref_predfun` have the correct structure?")
      }
      linpred_out <- unname(linpred_out)

      # Observation weights are not needed here, so use `wrhs = NULL` to avoid
      # potential conflicts for a non-`NULL` default `wrhs`:
      offs <- extract_model_data(fit, newdata = newdata, wrhs = NULL)$offset
      if (length(offs) > 0) {
        stopifnot(length(offs) %in% c(1L, n_obs))
        linpred_out <- linpred_out - offs
      }
      return(linpred_out)
    }
  } else {
    if (!is.null(ref_predfun)) {
      warning("Ignoring argument `ref_predfun` because `object` is `NULL`.")
    }
    ref_predfun <- function(fit, newdata = NULL) {
      stopifnot(is.null(fit))
      if (is.null(newdata)) {
        return(matrix(rep(NA, nrow(data))))
      } else {
        return(matrix(rep(NA, nrow(newdata))))
      }
    }
  }

  if (is.null(div_minimizer)) {
    if (!aug_data) {
      div_minimizer <- divmin
    } else {
      div_minimizer <- divmin_augdat
    }
  }

  if (is.null(proj_predfun)) {
    if (!aug_data) {
      proj_predfun <- subprd
    } else if (family$family == "binomial") {
      proj_predfun <- subprd_augdat_binom
    } else {
      proj_predfun <- subprd_augdat
    }
  }
  if (aug_data) {
    proj_predfun_usr <- proj_predfun
    proj_predfun <- function(fits, newdata) {
      augprd_arr <- proj_predfun_usr(fits, newdata = newdata)
      return(arr2augmat(augprd_arr))
    }
  }

  fetch_data_wrapper <- function(obs = NULL) {
    fetch_data(data, obs = obs)
  }

  if (is.null(cvfun)) {
    if (!proper_model) {
      # This is a dummy definition for cvfun(), but it will lead to standard CV
      # for `datafit`s; see cv_varsel() and .get_kfold():
      cvfun <- function(folds) {
        lapply(seq_len(max(folds)), function(k) list())
      }
    }
  }

  if (is.null(cvrefbuilder)) {
    if (proper_model) {
      cvrefbuilder <- get_refmodel
    } else {
      cvrefbuilder <- function(cvfit) {
        init_refmodel(
          object = NULL,
          data = fetch_data_wrapper(obs = setdiff(seq_len(nrow(data)),
                                                  cvfit$omitted)),
          formula = formula,
          family = family,
          div_minimizer = div_minimizer,
          proj_predfun = proj_predfun,
          extract_model_data = extract_model_data
        )
      }
    }
  }

  # Data --------------------------------------------------------------------

  model_data <- extract_model_data(object, newdata = data)
  weights <- model_data$weights
  offset <- model_data$offset
  if (latent_proj) {
    y <- rowMeans(ref_predfun(object))
  } else {
    y <- model_data$y
  }

  # Add (transformed) response under the (possibly) new name:
  data[, response_name] <- y

  target <- .get_standard_y(y, weights, family)
  y <- target$y
  weights <- target$weights

  if (aug_data) {
    y <- as.factor(y)
    stopifnot(nlevels(y) >= 2)
    if (!identical(levels(y), family$cats)) {
      stop("The levels of the response variable (after coercing it to a ",
           "`factor`) have to be identical to `family$cats`. See the ",
           "documentation for extend_family()'s argument `augdat_y_unqs` to ",
           "solve this.")
    }
  } else if (family$family == "binomial") {
    if (!all(.is.wholenumber(y))) {
      stop("In projpred, the response must contain numbers of successes (not ",
           "proportions of successes), in contrast to glm() where this is ",
           "possible for a 1-column response if the multiplication with the ",
           "weights gives whole numbers.")
    } else if (all(y %in% c(0, 1)) &&
               length(response_name) == 1 &&
               !all(weights == 1)) {
      warning("Assuming that the response contains numbers of successes (not ",
              "proportions of successes), in contrast to glm().")
    }
  }

  if (aug_data && !all(weights == 1)) {
    stop("Currently, the augmented-data projection may not be combined with ",
         "observation weights (other than 1).")
  }
  if (latent_proj && !all(weights == 1)) {
    stop("Currently, the latent projection may not be combined with ",
         "observation weights (other than 1).")
  }

  if (is.null(offset)) {
    offset <- rep(0, length(y))
  }

  # For avoiding the warning "contrasts dropped from factor <factor_name>" when
  # predicting for each projected draw, e.g., for submodels fit with lm()/glm():
  has_contr <- sapply(data, function(data_col) {
    !is.null(attr(data_col, "contrasts"))
  })
  for (idx_col in which(has_contr)) {
    attr(data[[idx_col]], "contrasts") <- NULL
  }

  # mu ----------------------------------------------------------------------

  # Note: For the augmented-data projection, in particular for nominal and
  # ordinal families with more than 2 categories, the final matrix `mu` is an
  # augmented-rows matrix containing the probabilities for each of the response
  # categories (at each observation and each posterior draw).
  if (proper_model) {
    eta <- ref_predfun(object)
    mu <- family$linkinv(eta)
  } else {
    if (family$family != "binomial") {
      mu <- y
    } else {
      mu <- y / weights
    }
    mu <- matrix(mu)
    eta <- family$linkfun(mu)
  }

  # Miscellaneous -----------------------------------------------------------

  ndraws <- ncol(mu)
  if (latent_proj) {
    ## latent noise is fixed
    dis <- rep(1, ndraws)
  } else if (is.null(dis)) {
    if (!.has_dispersion(family)) {
      dis <- rep(NA, ndraws)
    } else {
      if (proper_model) {
        stop("Please supply argument `dis`.")
      } else {
        dis <- 0
      }
    }
  } else {
    stopifnot(length(dis) == ndraws)
  }

  if (proper_model) {
    loglik <- t(family$ll_fun(
      family$linkinv(eta + offset), dis, y, weights = weights
    ))
  } else {
    loglik <- NULL
  }
  if (proper_model && latent_proj) {
    loglik_forPSIS <- log_lik(object)
  } else {
    loglik_forPSIS <- NULL
  }

  # Equal sample (draws) weights by default:
  wsample <- rep(1 / ndraws, ndraws)

  intercept <- as.logical(attr(terms(formula), "intercept"))
  if (!intercept) {
    stop("Reference models without an intercept are currently not supported.")
  }

  # Output ------------------------------------------------------------------

  refmodel <- nlist(
    fit = object, formula, div_minimizer, family, mu, eta, dis, y, loglik,
    loglik_forPSIS, intercept, proj_predfun, fetch_data = fetch_data_wrapper,
    wobs = weights, wsample, offset, cvfun, cvfits, extract_model_data,
    ref_predfun, cvrefbuilder
  )
  if (proper_model) {
    class(refmodel) <- "refmodel"
  } else {
    class(refmodel) <- c("datafit", "refmodel")
  }

  return(refmodel)
}
