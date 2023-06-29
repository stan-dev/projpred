# Common documentation ----------------------------------------------------

#' Reference model and more general information
#'
#' @description
#'
#' Function [get_refmodel()] is a generic function whose methods usually call
#' [init_refmodel()] which is the underlying workhorse (and may also be used
#' directly without a call to [get_refmodel()]).
#'
#' Both, [get_refmodel()] and [init_refmodel()], create an object containing
#' information needed for the projection predictive variable selection, namely
#' about the reference model, the submodels, and how the projection should be
#' carried out. For the sake of simplicity, the documentation may refer to the
#' resulting object also as "reference model" or "reference model object", even
#' though it also contains information about the submodels and the projection.
#'
#' A "typical" reference model object is created by [get_refmodel.stanreg()] and
#' [brms::get_refmodel.brmsfit()], either implicitly by a call to a top-level
#' function such as [project()], [varsel()], and [cv_varsel()] or explicitly by
#' a call to [get_refmodel()]. All non-"typical" reference model objects will be
#' called "custom" reference model objects.
#'
#' Some arguments are for \eqn{K}-fold cross-validation (\eqn{K}-fold CV) only;
#' see [cv_varsel()] for the use of \eqn{K}-fold CV in \pkg{projpred}.
#'
#' @name refmodel-init-get
#'
#' @inheritParams extend_family
#' @param object For [init_refmodel()], an object that the functions from
#'   arguments `extract_model_data` and `ref_predfun` can be applied to, with a
#'   `NULL` object being treated specially (see section "Value" below). For
#'   [get_refmodel.default()], an object of type `list` that (i) function
#'   [family()] can be applied to in order to retrieve the family (if argument
#'   `family` is `NULL`) and (ii) has an element called `data` containing the
#'   original dataset (see argument `data` of [init_refmodel()]), additionally
#'   to the properties required for [init_refmodel()]. For non-default methods
#'   of [get_refmodel()], an object of the corresponding class.
#' @param data A `data.frame` containing the data to use for the projection
#'   predictive variable selection. Any `contrasts` attributes of the dataset's
#'   columns are silently removed. For custom reference models, the columns of
#'   `data` do not necessarily have to coincide with those of the dataset used
#'   for fitting the reference model, but keep in mind that a row-subset of
#'   `data` is used for argument `newdata` of `ref_predfun` during \eqn{K}-fold
#'   CV.
#' @param formula The full formula to use for the search procedure. For custom
#'   reference models, this does not necessarily coincide with the reference
#'   model's formula. For general information about formulas in \R, see
#'   [`formula`]. For information about possible right-hand side (i.e.,
#'   predictor) terms in `formula` here, see the main vignette and section
#'   "Formula terms" below. For multilevel formulas, see also package \pkg{lme4}
#'   (in particular, functions [lme4::lmer()] and [lme4::glmer()]). For additive
#'   formulas, see also packages \pkg{mgcv} (in particular, function
#'   [mgcv::gam()]) and \pkg{gamm4} (in particular, function [gamm4::gamm4()]).
#' @param ref_predfun Prediction function for the linear predictor of the
#'   reference model, including offsets (if existing). See also section
#'   "Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`" below. If
#'   `object` is `NULL`, `ref_predfun` is ignored and an internal default is
#'   used instead.
#' @param proj_predfun Prediction function for the linear predictor of a
#'   submodel onto which the reference model is projected. See also section
#'   "Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`" below.
#' @param div_minimizer A function for minimizing the Kullback-Leibler (KL)
#'   divergence from the reference model to a submodel (i.e., for performing the
#'   projection of the reference model onto a submodel). The output of
#'   `div_minimizer` is used, e.g., by `proj_predfun`'s argument `fits`. See
#'   also section "Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`"
#'   below.
#' @param extract_model_data A function for fetching some variables (response,
#'   observation weights, offsets) from the original dataset (supplied to
#'   argument `data`) or from a new dataset. See also section "Argument
#'   `extract_model_data`" below.
#' @param family An object of class `family` representing the observation model
#'   (i.e., the distributional family for the response) of the *submodels*.
#'   (However, the link and the inverse-link function of this `family` are also
#'   used for quantities like predictions and fitted values related to the
#'   *reference model*.) May be `NULL` for [get_refmodel.default()] in which
#'   case the family is retrieved from `object`. For custom reference models,
#'   `family` does not have to coincide with the family of the reference model
#'   (if the reference model possesses a formal `family` at all). In typical
#'   reference models, however, these families do coincide. Furthermore, the
#'   latent projection is an exception where `family` is not the family of the
#'   submodels (in that case, the family of the submodels is the [gaussian()]
#'   family).
#' @param cvfits For \eqn{K}-fold CV only. A `list` containing a sub-`list`
#'   called `fits` containing the \eqn{K} model fits from which reference model
#'   structures are created. The `cvfits` `list` (i.e., the super-`list`) needs
#'   to have an attribute called `folds`, consisting of an integer vector giving
#'   the fold indices (one fold index per observation). Each element of
#'   `cvfits$fits` (i.e., each of the \eqn{K} model fits) needs to be a list.
#'   Only one of `cvfits` and `cvfun` needs to be provided (for \eqn{K}-fold
#'   CV). Note that `cvfits` takes precedence over `cvfun`, i.e., if both are
#'   provided, `cvfits` is used.
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
#' @param called_from_cvrefbuilder A single logical value indicating whether
#'   [init_refmodel()] is called from a `cvrefbuilder` function (`TRUE`) or not
#'   (`FALSE`). Currently, `TRUE` only causes some warnings to be suppressed
#'   (warnings which don't need to be thrown for each of the \eqn{K} reference
#'   model objects because it is sufficient to throw them for the original
#'   reference model object only). This argument is mainly for internal use, but
#'   may also be helpful for users with a custom `cvrefbuilder` function.
#' @param dis A vector of posterior draws for the reference model's dispersion
#'   parameter or---more precisely---the posterior values for the reference
#'   model's parameter-conditional predictive variance (assuming that this
#'   variance is the same for all observations). May be `NULL` if the submodels
#'   have no dispersion parameter or if the submodels do have a dispersion
#'   parameter, but `object` is `NULL` (in which case `0` is used for `dis`).
#'   Note that for the [gaussian()] `family`, `dis` is the standard deviation,
#'   not the variance.
#' @param ... For [get_refmodel.default()] and [get_refmodel.stanreg()]:
#'   arguments passed to [init_refmodel()]. For the [get_refmodel()] generic:
#'   arguments passed to the appropriate method. For [init_refmodel()]:
#'   arguments passed to [extend_family()] (apart from `family`).
#'
#' @details
#'
#' # Formula terms
#'
#' Although bad practice (in general), a reference model lacking an intercept
#' can be used within \pkg{projpred}. However, it will always be projected onto
#' submodels which *include* an intercept. The reason is that even if the true
#' intercept in the reference model is zero, this does not need to hold for the
#' submodels.
#'
#' In multilevel (group-level) terms, function calls on the right-hand side of
#' the `|` character (e.g., `(1 | gr(group_variable))`, which is possible in
#' \pkg{brms}) are currently not allowed in \pkg{projpred}.
#'
#' For additive models (still an experimental feature), only [mgcv::s()] and
#' [mgcv::t2()] are currently supported as smooth terms. Furthermore, these need
#' to be called without any arguments apart from the predictor names (symbols).
#' For example, for smoothing the effect of a predictor `x`, only `s(x)` or
#' `t2(x)` are allowed. As another example, for smoothing the joint effect of
#' two predictors `x` and `z`, only `s(x, z)` or `t2(x, z)` are allowed (and
#' analogously for higher-order joint effects, e.g., of three predictors). Note
#' that all smooth terms need to be included in `formula` (there is no `random`
#' argument as in [rstanarm::stan_gamm4()], for example).
#'
#' # Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer`
#'
#' Arguments `ref_predfun`, `proj_predfun`, and `div_minimizer` may be `NULL`
#' for using an internal default (see [projpred-package] for the functions used
#' by the default divergence minimizers). Otherwise, let \eqn{N} denote the
#' number of observations (in case of CV, these may be reduced to each fold),
#' \eqn{S_{\mathrm{ref}}}{S_ref} the number of posterior draws for the reference
#' model's parameters, and \eqn{S_{\mathrm{prj}}}{S_prj} the number of draws for
#' the parameters of a submodel that the reference model has been projected onto
#' (short: the number of projected draws). For the augmented-data projection,
#' let \eqn{C_{\mathrm{cat}}}{C_cat} denote the number of response categories,
#' \eqn{C_{\mathrm{lat}}}{C_lat} the number of latent response categories (which
#' typically equals \eqn{C_{\mathrm{cat}} - 1}{C_cat - 1}), and define
#' \eqn{N_{\mathrm{augcat}} := N \cdot C_{\mathrm{cat}}}{N_augcat := N * C_cat}
#' as well as \eqn{N_{\mathrm{auglat}} := N \cdot C_{\mathrm{lat}}}{N_auglat :=
#' N * C_lat}. Then the functions supplied to these arguments need to have the
#' following prototypes:
#' * `ref_predfun`: `ref_predfun(fit, newdata = NULL)` where:
#'     + `fit` accepts the reference model fit as given in argument `object`
#'     (but possibly re-fitted to a subset of the observations, as done in
#'     \eqn{K}-fold CV).
#'     + `newdata` accepts either `NULL` (for using the original dataset,
#'     typically stored in `fit`) or data for new observations (at least in the
#'     form of a `data.frame`).
#' * `proj_predfun`: `proj_predfun(fits, newdata)` where:
#'     + `fits` accepts a `list` of length \eqn{S_{\mathrm{prj}}}{S_prj}
#'     containing this number of submodel fits. This `list` is the same as that
#'     returned by [project()] in its output element `outdmin` (which in turn is
#'     the same as the return value of `div_minimizer`, except if [project()]
#'     was used with an `object` of class `vsel` based on an L1 search as well
#'     as with `refit_prj = FALSE`).
#'     + `newdata` accepts data for new observations (at least in the form of a
#'     `data.frame`).
#' * `div_minimizer` does not need to have a specific prototype, but it needs to
#' be able to be called with the following arguments:
#'     + `formula` accepts either a standard [`formula`] with a single response
#'     (if \eqn{S_{\mathrm{prj}} = 1}{S_prj = 1} or in case of the
#'     augmented-data projection) or a [`formula`] with \eqn{S_{\mathrm{prj}} >
#'     1}{S_prj > 1} response variables [cbind()]-ed on the left-hand side in
#'     which case the projection has to be performed for each of the response
#'     variables separately.
#'     + `data` accepts a `data.frame` to be used for the projection. In case of
#'     the traditional or the latent projection, this dataset has \eqn{N} rows.
#'     In case of the augmented-data projection, this dataset has
#'     \eqn{N_{\mathrm{augcat}}}{N_augcat} rows.
#'     + `family` accepts an object of class `family`.
#'     + `weights` accepts either observation weights (at least in the form of a
#'     numeric vector) or `NULL` (for using a vector of ones as weights).
#'     + `projpred_var` accepts an \eqn{N \times S_{\mathrm{prj}}}{N x S_prj}
#'     matrix of predictive variances (necessary for \pkg{projpred}'s internal
#'     GLM fitter) in case of the traditional or the latent projection and an
#'     \eqn{N_{\mathrm{augcat}} \times S_{\mathrm{prj}}}{N_augcat x S_prj}
#'     matrix (containing only `NA`s) in case of the augmented-data projection.
#'     + `projpred_regul` accepts a single numeric value as supplied to argument
#'     `regul` of [project()], for example.
#'     + `projpred_ws_aug` accepts an \eqn{N \times S_{\mathrm{prj}}}{N x S_prj}
#'     matrix of expected values for the response in case of the traditional or
#'     the latent projection and an \eqn{N_{\mathrm{augcat}} \times
#'     S_{\mathrm{prj}}}{N_augcat x S_prj} matrix of probabilities for the
#'     response categories in case of the augmented-data projection.
#'     + `...` accepts further arguments specified by the user.
#'
#' The return value of these functions needs to be:
#' * `ref_predfun`: for the traditional or the latent projection, an \eqn{N
#' \times S_{\mathrm{ref}}}{N x S_ref} matrix; for the augmented-data
#' projection, an \eqn{S_{\mathrm{ref}} \times N \times C_{\mathrm{lat}}}{S_ref
#' x N x C_lat} array (the only exception is the augmented-data projection for
#' the [binomial()] family in which case `ref_predfun` needs to return an \eqn{N
#' \times S_{\mathrm{ref}}}{N x S_ref} matrix just like for the traditional
#' projection because the array is constructed by an internal wrapper function).
#' * `proj_predfun`: for the traditional or the latent projection, an \eqn{N
#' \times S_{\mathrm{prj}}}{N x S_prj} matrix; for the augmented-data
#' projection, an \eqn{N \times C_{\mathrm{lat}} \times S_{\mathrm{prj}}}{N x
#' C_lat x S_prj} array.
#' * `div_minimizer`: a `list` of length \eqn{S_{\mathrm{prj}}}{S_prj}
#' containing this number of submodel fits.
#'
#' # Argument `extract_model_data`
#'
#' The function supplied to argument `extract_model_data` needs to have the
#' prototype
#' ```{r, eval = FALSE}
#' extract_model_data(object, newdata, wrhs = NULL, orhs = NULL,
#'                    extract_y = TRUE)
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
#' exception is that `y` may also be `NULL` (depending on argument `extract_y`),
#' a non-numeric vector, or a `factor`.
#'
#' The weights and offsets returned by `extract_model_data` will be assumed to
#' hold for the reference model as well as for the submodels.
#'
#' # Augmented-data projection
#'
#' If a custom reference model for an augmented-data projection is needed, see
#' also [extend_family()].
#'
#' For the augmented-data projection, the response vector resulting from
#' `extract_model_data` is internally coerced to a `factor` (using
#' [as.factor()]). The levels of this `factor` have to be identical to
#' `family$cats` (*after* applying [extend_family()] internally; see
#' [extend_family()]'s argument `augdat_y_unqs`).
#'
#' Note that response-specific offsets (i.e., one length-\eqn{N} offset vector
#' per response category) are not supported by \pkg{projpred} yet. So far, only
#' offsets which are the same across all response categories are supported. This
#' is why in case of the [brms::categorical()] family, offsets are currently not
#' supported at all.
#'
#' Currently, `object = NULL` (i.e., a `datafit`; see section "Value") is not
#' supported in case of the augmented-data projection.
#'
#' # Latent projection
#'
#' If a custom reference model for a latent projection is needed, see also
#' [extend_family()].
#'
#' For the latent projection, `family$cats` (*after* applying [extend_family()]
#' internally; see [extend_family()]'s argument `latent_y_unqs`) currently must
#' not be `NULL` if the original (i.e., non-latent) response is a `factor`.
#' Conversely, if `family$cats` (*after* applying [extend_family()]) is
#' non-`NULL`, the response vector resulting from `extract_model_data` is
#' internally coerced to a `factor` (using [as.factor()]). The levels of this
#' `factor` have to be identical to that non-`NULL` element `family$cats`.
#'
#' Currently, `object = NULL` (i.e., a `datafit`; see section "Value") is not
#' supported in case of the latent projection.
#'
#' @return An object that can be passed to all the functions that take the
#'   reference model fit as the first argument, such as [varsel()],
#'   [cv_varsel()], [project()], [proj_linpred()], and [proj_predict()].
#'   Usually, the returned object is of class `refmodel`. However, if `object`
#'   is `NULL`, the returned object is of class `datafit` as well as of class
#'   `refmodel` (with `datafit` being first). Objects of class `datafit` are
#'   handled differently at several places throughout this package.
#'
#'   The elements of the returned object are not meant to be accessed directly
#'   but instead via downstream functions (see the functions mentioned above as
#'   well as [predict.refmodel()]).
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

# Function definitions ----------------------------------------------------

#' Predictions or log posterior predictive densities from a reference model
#'
#' This is the [predict()] method for `refmodel` objects (returned by
#' [get_refmodel()] or [init_refmodel()]). It offers three types of output which
#' are all based on the reference model and new (or old) observations: Either
#' the linear predictor on link scale, the linear predictor transformed to
#' response scale, or the log posterior predictive density.
#'
#' @template args-newdata
#' @param object An object of class `refmodel` (returned by [get_refmodel()] or
#'   [init_refmodel()]).
#' @param ynew If not `NULL`, then this needs to be a vector of new (or old)
#'   response values. See also section "Value" below. In case of (i) the
#'   augmented-data projection or (ii) the latent projection with `type =
#'   "response"` and `object$family$cats` being not `NULL`, `ynew` is internally
#'   coerced to a `factor` (using [as.factor()]). The levels of this `factor`
#'   have to be a subset of `object$family$cats` (see [extend_family()]'s
#'   arguments `augdat_y_unqs` and `latent_y_unqs`, respectively).
#' @param type Usually only relevant if `is.null(ynew)`, but for the latent
#'   projection, this also affects the `!is.null(ynew)` case (see below). The
#'   scale on which the predictions are returned, either `"link"` or
#'   `"response"` (see [predict.glm()] but note that [predict.refmodel()] does
#'   not adhere to the typical \R convention of a default prediction on link
#'   scale). For both scales, the predictions are averaged across the posterior
#'   draws. In case of the latent projection, argument `type` is similar in
#'   spirit to argument `resp_oscale` from other functions: If (i)
#'   `is.null(ynew)`, then argument `type` affects the predictions as described
#'   above. In that case, note that `type = "link"` yields the linear predictors
#'   without any modifications that may be due to the original response
#'   distribution (e.g., for a [brms::cumulative()] model, the ordered
#'   thresholds are not taken into account). If (ii) `!is.null(ynew)`, then
#'   argument `type` also affects the scale of the log posterior predictive
#'   densities (`type = "response"` for the original response scale, `type =
#'   "link"` for the latent Gaussian scale).
#' @param ... Currently ignored.
#'
#' @details Argument `weightsnew` is only relevant if `!is.null(ynew)`.
#'
#'   In case of a multilevel reference model, group-level effects for new group
#'   levels are drawn randomly from a (multivariate) Gaussian distribution. When
#'   setting `projpred.mlvl_pred_new` to `TRUE`, all group levels from `newdata`
#'   (even those that already exist in the original dataset) are treated as new
#'   group levels (if `is.null(newdata)`, all group levels from the original
#'   dataset are considered as new group levels in that case).
#'
#' @return In the following, \eqn{N}, \eqn{C_{\mathrm{cat}}}{C_cat}, and
#'   \eqn{C_{\mathrm{lat}}}{C_lat} from help topic [refmodel-init-get] are used.
#'   Furthermore, let \eqn{C} denote either \eqn{C_{\mathrm{cat}}}{C_cat} (if
#'   `type = "response"`) or \eqn{C_{\mathrm{lat}}}{C_lat} (if `type = "link"`).
#'   Then, if `is.null(ynew)`, the returned object contains the reference
#'   model's predictions (with the scale depending on argument `type`) as:
#'   * a length-\eqn{N} vector in case of (i) the traditional projection, (ii)
#'   the latent projection with `type = "link"`, or (iii) the latent projection
#'   with `type = "response"` and `object$family$cats` being `NULL`;
#'   * an \eqn{N \times C}{N x C} matrix in case of (i) the augmented-data
#'   projection or (ii) the latent projection with `type = "response"` and
#'   `object$family$cats` being not `NULL`.
#'
#'   If `!is.null(ynew)`, the returned object is a length-\eqn{N} vector of log
#'   posterior predictive densities evaluated at `ynew`.
#'
#' @export
predict.refmodel <- function(object, newdata = NULL, ynew = NULL,
                             offsetnew = NULL, weightsnew = NULL,
                             type = "response", ...) {
  if (inherits(object, "datafit")) {
    stop("Cannot make predictions for an `object` of class \"datafit\".")
  }
  refmodel <- object
  if (!type %in% c("response", "link")) {
    stop("type should be one of ('response', 'link')")
  }
  if (!is.null(ynew) && (!is.numeric(ynew) || NCOL(ynew) != 1) &&
      is.null(refmodel$family$cats)) {
    stop("Argument `ynew` must be a numeric vector.")
  }
  if (!is.null(ynew) && !is.null(refmodel$family$cats) &&
      (!refmodel$family$for_latent || type == "response")) {
    ynew <- as.factor(ynew)
    if (!all(levels(ynew) %in% refmodel$family$cats)) {
      if (refmodel$family$for_augdat) {
        y_unqs_str <- "augdat_y_unqs"
      } else {
        y_unqs_str <- "latent_y_unqs"
      }
      stop("The levels of the response variable (after coercing it to a ",
           "`factor`) have to be a subset of `family$cats`. Either modify ",
           "`ynew` accordingly or see the documentation for extend_family()'s ",
           "argument `", y_unqs_str, "` to solve this.")
    }
    # Re-assign the original levels because some levels might be missing:
    ynew <- factor(ynew, levels = refmodel$family$cats)
  } else if (!is.null(ynew) &&
             refmodel$family$for_latent &&
             is.null(refmodel$family$cats) &&
             (is.factor(ynew) || is.character(ynew) || is.logical(ynew))) {
    stop("If the original (i.e., non-latent) response is `factor`-like, ",
         "`family$cats` must not be `NULL`. See the documentation for ",
         "extend_family()'s argument `latent_y_unqs` to solve this.")
  }

  if (!is.null(newdata)) {
    newdata <- na.fail(newdata)
  }
  nobs_new <- nrow(newdata) %||% refmodel$nobs
  w_o <- refmodel$extract_model_data(refmodel$fit, newdata = newdata,
                                     wrhs = weightsnew, orhs = offsetnew,
                                     extract_y = FALSE)
  weightsnew <- w_o$weights
  offsetnew <- w_o$offset
  if (length(weightsnew) == 0) {
    weightsnew <- rep(1, nobs_new)
  }
  if (length(offsetnew) == 0) {
    offsetnew <- rep(0, nobs_new)
  }
  if (refmodel$family$for_augdat && !all(weightsnew == 1)) {
    stop("Currently, the augmented-data projection may not be combined with ",
         "observation weights (other than 1).")
  }
  if (refmodel$family$for_latent && !all(weightsnew == 1)) {
    stop("Currently, the latent projection may not be combined with ",
         "observation weights (other than 1).")
  }
  if (!is.null(newdata) && inherits(refmodel$fit, "stanreg") &&
      length(refmodel$fit$offset) > 0) {
    if ("projpred_internal_offs_stanreg" %in% names(newdata)) {
      stop("Need to write to column `projpred_internal_offs_stanreg` of ",
           "`newdata`, but that column already exists. Please rename this ",
           "column in `newdata` and try again.")
    }
    newdata$projpred_internal_offs_stanreg <- offsetnew
  }

  ## ref_predfun returns eta = link(mu)
  eta <- refmodel$ref_predfun(refmodel$fit, newdata = newdata,
                              excl_offs = FALSE)

  if (is.null(ynew)) {
    if (type == "link") {
      pred <- eta
    } else {
      if (refmodel$family$for_latent) {
        pred <- refmodel$family$latent_ilink(
          t(eta), cl_ref = seq_along(refmodel$wdraws_ref),
          wdraws_ref = rep(1, length(refmodel$wdraws_ref))
        )
        if (length(dim(pred)) < 2) {
          stop("Unexpected structure for the output of `latent_ilink`.")
        }
        if (length(dim(pred)) == 3) {
          pred <- arr2augmat(pred, margin_draws = 1)
        }
        if (all(is.na(pred))) {
          message(
            "`latent_ilink` returned only `NA`s, so the output will also be ",
            "`NA` as long as `type = \"response\"`."
          )
        }
      } else {
        pred <- refmodel$family$linkinv(eta)
      }
    }
    was_augmat <- inherits(pred, "augmat")
    ## integrate over the draws
    if (type == "link" || !refmodel$family$for_latent || was_augmat) {
      if (ncol(pred) > 1) {
        pred <- rowMeans(pred)
      }
    } else {
      if (nrow(pred) > 1) {
        pred <- colMeans(pred)
      }
    }
    if (was_augmat) {
      pred <- structure(pred, nobs_orig = nobs_new, class = "augvec")
      pred <- augmat2arr(augvec2augmat(pred))
      pred <- matrix(pred, nrow = dim(pred)[1], ncol = dim(pred)[2])
    }
    return(pred)
  } else {
    ## evaluate the log posterior predictive density at the given ynew values
    if (refmodel$family$for_latent && type == "response") {
      mu_oscale <- refmodel$family$latent_ilink(
        t(eta), cl_ref = seq_along(refmodel$wdraws_ref),
        wdraws_ref = rep(1, length(refmodel$wdraws_ref))
      )
      if (length(dim(mu_oscale)) < 2) {
        stop("Unexpected structure for the output of `latent_ilink`.")
      }
      loglik <- refmodel$family$latent_ll_oscale(
        mu_oscale, y_oscale = ynew, wobs = weightsnew,
        cl_ref = seq_along(refmodel$wdraws_ref),
        wdraws_ref = rep(1, length(refmodel$wdraws_ref))
      )
      if (!is.matrix(loglik)) {
        stop("Unexpected structure for the output of `latent_ll_oscale`.")
      }
      if (all(is.na(mu_oscale))) {
        message(
          "`latent_ilink` returned only `NA`s, so the output will also be ",
          "`NA` as long as `type = \"response\"`."
        )
      } else if (all(is.na(loglik))) {
        message(
          "`latent_ll_oscale` returned only `NA`s, so the output will also be ",
          "`NA` as long as `type = \"response\"`."
        )
      }
      S <- nrow(loglik)
      marg_obs <- 2
    } else {
      if (refmodel$family$for_latent) {
        if (all(is.na(refmodel$dis))) {
          message(
            "Cannot calculate LPD values if `type = \"link\"` and ",
            "`<refmodel>$dis` consists of only `NA`s. If it's not possible to ",
            "supply a suitable argument `dis` to init_refmodel(), consider ",
            "switching to `type = \"response\"` (which might require the ",
            "specification of functions needed by extend_family())."
          )
        }
        if (is.null(newdata)) {
          newdata_lat <- newdata
          if (inherits(refmodel$fit, "stanreg") &&
              length(refmodel$fit$offset) > 0) {
            newdata_lat <- refmodel$fetch_data()
            newdata_lat$projpred_internal_offs_stanreg <- offsetnew
          }
          ynew <- rowMeans(refmodel$ref_predfun(
            fit = refmodel$fit,
            newdata = newdata_lat,
            excl_offs = FALSE,
            mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)
          ))
        }
      }
      loglik <- refmodel$family$ll_fun(
        refmodel$family$linkinv(eta), refmodel$dis, ynew, weightsnew
      )
      S <- ncol(loglik)
      marg_obs <- 1
    }
    lpd <- apply(loglik, marg_obs, log_sum_exp) - log(S)
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
    y <- eval_el2(resp_form, newdata)
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
get_refmodel.stanreg <- function(object, latent = FALSE, dis = NULL, ...) {
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

  # Data --------------------------------------------------------------------

  data <- object$data
  stopifnot(is.data.frame(data))

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
    offs <- extract_model_data(fit, newdata = newdata, wrhs = NULL,
                               extract_y = FALSE)$offset
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
    aug_data <- fit$stan_function == "stan_polr" && !latent
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
    get_refmodel(cvfit, latent = latent, dis = dis,
                 called_from_cvrefbuilder = TRUE, ...)
  }

  # Miscellaneous -----------------------------------------------------------

  if (is.null(dis) && !latent && has_dispersion(family)) {
    dis <- data.frame(object)[, "sigma"]
  }

  # Augmented-data projection -----------------------------------------------

  aug_data <- object$stan_function == "stan_polr" && !latent
  if (aug_data) {
    args_augdat <- list(
      augdat_link = augdat_link_cumul,
      augdat_ilink = augdat_ilink_cumul,
      augdat_args_link = list(link = family$link),
      augdat_args_ilink = list(link = family$link)
    )
  } else {
    args_augdat <- list()
  }

  # Latent projection -------------------------------------------------------

  args_latent <- list(latent = latent)
  if (latent) {
    if (object$stan_function == "stan_polr") {
      draws_mat <- as.matrix(object)
      thres_nms <- names(object$zeta)
      thres_draws <- draws_mat[, thres_nms, drop = FALSE]
      latent_ilink_tmp <- function(lpreds, cl_ref,
                                   wdraws_ref = rep(1, length(cl_ref))) {
        thres_agg <- cl_agg(thres_draws, cl = cl_ref,
                            wdraws = wdraws_ref)
        lpreds_thres <- apply(thres_agg, 2, function(thres_agg_c) {
          # Notes on dimensionalities (with S_agg = `nrow(lpreds)`):
          # * `thres_agg` is S_agg x C_lat (with C_lat = `ncats - 1L` =
          #   `nthres`) and thus `thres_agg_c` is a vector of length S_agg,
          # * `lpreds` is S_agg x N (with N denoting the number of (possibly
          #   new) observations (not necessarily the original number of
          #   observations)).
          thres_agg_c - lpreds
        }, simplify = FALSE)
        # Coerce to an S_agg x N x C_lat array:
        lpreds_thres <- do.call(abind::abind, c(lpreds_thres, rev.along = 0))
        # Transform to response space, yielding an S_agg x N x C_cat array:
        return(augdat_ilink_cumul(lpreds_thres, link = family$link))
      }
      args_latent <- c(args_latent, list(latent_ilink = latent_ilink_tmp))
      # Free up some memory:
      rm(draws_mat)
    }
    # TODO (latent): Add response-scale support for more families: For
    # response-scale support, they all need a specific `latent_ilink` function;
    # some families (those for which the response can be numeric) also require
    # specific `latent_ll_oscale` and `latent_ppd_oscale` functions. The
    # binomial family has response-scale support implemented natively in
    # projpred.
  }

  # Output ------------------------------------------------------------------

  args_basic <- list(
    object = object, data = data, formula = formula, family = family,
    ref_predfun = ref_predfun, extract_model_data = extract_model_data,
    dis = dis, cvfun = cvfun, cvrefbuilder = cvrefbuilder
  )
  return(do.call(init_refmodel, args = c(args_basic, args_augdat, args_latent,
                                         list(...))))
}

#' @rdname refmodel-init-get
#' @export
init_refmodel <- function(object, data, formula, family, ref_predfun = NULL,
                          div_minimizer = NULL, proj_predfun = NULL,
                          extract_model_data, cvfun = NULL,
                          cvfits = NULL, dis = NULL, cvrefbuilder = NULL,
                          called_from_cvrefbuilder = FALSE, ...) {
  # Family ------------------------------------------------------------------

  family <- extend_family(family, ...)

  if (!called_from_cvrefbuilder) {
    if (family$family == "Student_t") {
      warning("Support for the `Student_t` family is still experimental.")
    } else if (family$family == "Gamma") {
      warning("Support for the `Gamma` family is still experimental.")
    }
  }

  family$mu_fun <- function(fits, obs = NULL, newdata = NULL, offset = NULL,
                            transform = TRUE) {
    newdata <- fetch_data(data, obs = obs, newdata = newdata)
    if (is.null(offset)) {
      offset <- rep(0, nrow(newdata))
    } else {
      stopifnot(length(offset) %in% c(1L, nrow(newdata)))
    }
    pred_sub <- proj_predfun(fits, newdata = newdata)
    if (family$family %in% fams_neg_linpred()) {
      pred_sub <- pred_sub - offset
    } else {
      pred_sub <- pred_sub + offset
    }
    if (transform) {
      pred_sub <- family$linkinv(pred_sub)
    }
    return(pred_sub)
  }

  if (family$family == "categorical" && family$link != "logit") {
    stop("For the brms::categorical() family, projpred only supports the ",
         "logit link.")
  }

  # Special case: `datafit` -------------------------------------------------

  proper_model <- !is.null(object)
  if (!proper_model && family$for_augdat) {
    stop("Currently, the augmented-data projection may not be combined with ",
         "`object = NULL` (i.e., a `datafit`).")
  } else if (!proper_model && family$for_latent) {
    stop("Currently, the latent projection may not be combined with ",
         "`object = NULL` (i.e., a `datafit`).")
  }

  # Formula -----------------------------------------------------------------

  stopifnot(inherits(formula, "formula"))
  data <- na.fail(data)
  stopifnot(is.data.frame(data))
  formula <- expand_formula(formula, data)
  if (!as.logical(attr(terms(formula), "intercept"))) {
    # Add an intercept to `formula` so that we always project onto submodels
    # *including* an intercept (see the discussion at #96):
    message("Adding an intercept to `formula` (the full-model formula used ",
            "for the search) so that the projection is always performed onto ",
            "submodels *including* an intercept.")
    formula <- update(formula, . ~ . + 1)
  }
  fml_extractions <- extract_terms_response(formula)
  response_name <- fml_extractions$response
  if (length(response_name) == 2) {
    if (family$family != "binomial") {
      stop("For non-binomial families, a two-column response is not allowed.")
    } else if (family$for_augdat) {
      stop("Currently, the augmented-data projection may not be combined with ",
           "a 2-column response.")
    }
  } else if (length(response_name) > 2) {
    stop("The response is not allowed to have more than two columns.")
  }
  # Remove parentheses from the response:
  response_name <- gsub("[()]", "", response_name)
  if (family$for_latent) {
    response_name <- paste0(".", response_name[1])
  }
  formula <- update(formula, paste(response_name[1], "~ ."))
  if (formula_contains_additive_terms(formula)) {
    if (family$for_augdat) {
      stop("Currently, the augmented-data projection may not be combined with ",
           "additive models.")
    } else if (getOption("projpred.warn_additive_experimental", TRUE) &&
               !called_from_cvrefbuilder) {
      warning("Support for additive models is still experimental.")
    }
  }
  if (formula_contains_group_terms(formula) &&
      getOption("projpred.warn_instable_projections", TRUE) &&
      !called_from_cvrefbuilder) {
    if (family$for_augdat) {
      warning(
        "For multilevel models, the augmented-data projection may not work ",
        "properly. The latent projection may be a remedy. See section ",
        "\"Troubleshooting\" of the main vignette for more information."
      )
    } else if (family$family == "binomial") {
      warning(
        "For multilevel binomial models, the traditional projection may not ",
        "work properly. The latent projection may be a remedy. See section ",
        "\"Troubleshooting\" of the main vignette for more information."
      )
    } else if (family$family == "poisson") {
      warning(
        "For multilevel Poisson models, the traditional projection may take ",
        "very long. The latent projection may be a remedy. See section ",
        "\"Troubleshooting\" of the main vignette for more information."
      )
    }
  }
  if (family$family == "categorical" &&
      length(fml_extractions$offset_terms) > 0) {
    stop("Currently, offsets are not supported in case of the ",
         "brms::categorical() family.")
  }

  # Functions ---------------------------------------------------------------

  if (proper_model) {
    if (is.null(ref_predfun)) {
      ref_predfun <- refprd
    }
    if (family$for_augdat && family$family == "binomial") {
      ref_predfun_mat <- ref_predfun
      # The assignment to a dummy object is just needed to avoid a `NOTE` in `R
      # CMD check`, namely "init_refmodel: multiple local function definitions
      # for 'ref_predfun' with different formal arguments":
      ref_predfun_dummy <- function(fit, newdata = NULL) {
        linpred1 <- ref_predfun_mat(fit = fit, newdata = newdata)
        linpred1 <- t(linpred1)
        return(array(linpred1, dim = c(dim(linpred1), 1L)))
      }
      ref_predfun <- ref_predfun_dummy
    }
    # Since posterior_linpred() is supposed to include any offsets, but in
    # general (i.e., in the default case `excl_offs = TRUE`, see below),
    # projpred currently expects the final ref_predfun() to exclude any offsets
    # (see issue #186), the offsets have to be subtracted (or added, in case of
    # some ordinal families). This is done here by defining the final
    # ref_predfun() as a wrapper function around the user-supplied (or
    # automatically derived) preliminary ref_predfun(). This wrapper function
    # also ensures that in the case `mlvl_allrandom = TRUE`, we draw new
    # group-level effects for *all* group levels (existing and new ones) and
    # performs some preparations for the augmented-data projection:
    ref_predfun_usr <- ref_predfun
    ref_predfun <- function(fit, newdata = NULL, excl_offs = TRUE,
                            mlvl_allrandom = getOption("projpred.mlvl_pred_new",
                                                       FALSE)) {
      if (length(fml_extractions$group_terms) > 0 && mlvl_allrandom) {
        # Need to replace existing group levels by dummy ones to ensure that we
        # draw new group-level effects for *all* group levels (existing and new
        # ones):
        if (is.null(newdata)) newdata <- data
        vnms <- flatten_group_terms(fml_extractions$group_terms)
        vnms <- sub("^.*\\|[[:blank:]]*", "", vnms)
        vnms <- sub("[[:blank:]]*\\)$", "", vnms)
        lvls_list <- lapply(setNames(nm = vnms), function(vnm) {
          if (!vnm %in% names(data)) {
            stop("Could not find column `", vnm, "` in `data`.")
          }
          if (!vnm %in% names(newdata)) {
            stop("Could not find column `", vnm, "` in `newdata`.")
          }
          from_fit <- unique(data[, vnm])
          from_new <- unique(newdata[, vnm])

          # Strictly speaking, this is not necessary (currently), but include it
          # for safety reasons, in case downstream code is changed in the future
          # (or in case the behavior of `factor`s in R is changed in general):
          if (is.factor(from_fit)) {
            from_fit <- as.character(from_fit)
          }
          if (is.factor(from_new)) {
            from_new <- as.character(from_new)
          }

          list(comb = union(from_fit, from_new),
               exist = from_fit,
               new = from_new)
        })
        for (vnm in vnms) {
          ex_lvl <- newdata[[vnm]] %in% lvls_list[[vnm]]$exist
          if (is.numeric(newdata[[vnm]])) {
            stopifnot(is.numeric(data[[vnm]]))
            if (!all(lvls_list[[vnm]]$exist >= 0)) {
              stop("In case of a numeric group variable, projpred requires ",
                   "this to have values >= 0.")
            }
            newdata[[vnm]][ex_lvl] <- max(lvls_list[[vnm]]$comb) + 1L +
              newdata[[vnm]][ex_lvl]
          } else if (is.character(newdata[[vnm]]) ||
                     is.factor(newdata[[vnm]])) {
            timestamp <- gsub("\\.", "", as.character(as.numeric(Sys.time())))
            dummy_lvls_ex <- paste("projpred_DUMMY", timestamp,
                                   newdata[[vnm]][ex_lvl], sep = "_")
            if (is.factor(newdata[[vnm]])) {
              orig_lvls <- levels(newdata[[vnm]])
              orig_ord <- is.ordered(newdata[[vnm]])
              newdata[[vnm]] <- as.character(newdata[[vnm]])
            } else {
              orig_lvls <- NULL
              orig_ord <- NULL
            }
            dummy_lvls <- unique(dummy_lvls_ex)
            if (any(dummy_lvls %in% lvls_list[[vnm]]$comb)) {
              stop("Need to assign dummy levels to existing group levels of ",
                   "variable `", vnm, "`, but encountered a conflict. Please ",
                   "try again or rename the group levels.")
            }
            newdata[[vnm]][ex_lvl] <- dummy_lvls_ex
            if (!is.null(orig_lvls) && !is.null(orig_ord)) {
              newdata[[vnm]] <- factor(newdata[[vnm]],
                                       levels = c(orig_lvls, dummy_lvls),
                                       ordered = orig_ord)
            }
          } else {
            stop("Unknown type of group variable. Please use factor, ",
                 "character, or numeric.")
          }
        }
      }

      linpred_out <- ref_predfun_usr(fit = fit, newdata = newdata)
      if (length(dim(linpred_out)) == 2) {
        n_obs <- nrow(linpred_out)
      } else if (length(dim(linpred_out)) == 3) {
        # For the augmented-data projection, `linpred_out` is expected to be a
        # 3-dimensional array with dimensions S_ref x N x C_lat (see
        # `?init_refmodel` for a definition of these dimensions). Therefore, it
        # is converted to an augmented-rows matrix (see `?`augdat-internals``
        # for a definition):
        linpred_out <- arr2augmat(linpred_out, margin_draws = 1)
        n_obs <- attr(linpred_out, "nobs_orig")
      } else {
        stop("Unexpected structure for `linpred_out`. Does the return value ",
             "of `ref_predfun` have the correct structure?")
      }
      linpred_out <- unname(linpred_out)

      if (excl_offs) {
        # Observation weights are not needed here, so use `wrhs = NULL` to avoid
        # potential conflicts for a non-`NULL` default `wrhs`:
        offs <- extract_model_data(fit, newdata = newdata, wrhs = NULL,
                                   extract_y = FALSE)$offset
        if (length(offs) > 0) {
          stopifnot(length(offs) %in% c(1L, n_obs))
          if (family$family %in% fams_neg_linpred()) {
            linpred_out <- linpred_out + offs
          } else {
            linpred_out <- linpred_out - offs
          }
        }
      }
      return(linpred_out)
    }
  } else {
    if (!is.null(ref_predfun)) {
      warning("Ignoring argument `ref_predfun` because `object` is `NULL`.")
    }
    ref_predfun <- function(fit, newdata = NULL, excl_offs = TRUE,
                            mlvl_allrandom = getOption("projpred.mlvl_pred_new",
                                                       FALSE)) {
      stopifnot(is.null(fit))
      if (is.null(newdata)) {
        return(matrix(rep(NA_real_, nrow(data))))
      } else {
        return(matrix(rep(NA_real_, nrow(newdata))))
      }
    }
  }

  if (is.null(div_minimizer)) {
    if (!family$for_augdat) {
      div_minimizer <- divmin
    } else {
      div_minimizer <- divmin_augdat
    }
  }

  if (is.null(proj_predfun)) {
    if (!family$for_augdat) {
      proj_predfun <- subprd
    } else if (family$family == "binomial") {
      proj_predfun <- subprd_augdat_binom
    } else {
      proj_predfun <- subprd_augdat
    }
  }
  if (family$for_augdat) {
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
      # This is a dummy definition for cvfun(), but the cvrefbuilder() function
      # defined below will lead to standard CV nonetheless (at least for the
      # submodels; for the reference model, we don't have an actual reference
      # model, only a `datafit`):
      cvfun <- function(folds) {
        lapply(seq_len(max(folds)), function(k) list())
      }
    }
  }

  if (is.null(cvrefbuilder)) {
    if (proper_model) {
      cvrefbuilder <- function(cvfit) {
        get_refmodel(cvfit, called_from_cvrefbuilder = TRUE)
      }
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
          extract_model_data = extract_model_data,
          called_from_cvrefbuilder = TRUE
        )
      }
    }
  }

  # Data --------------------------------------------------------------------

  model_data <- extract_model_data(object, newdata = data, extract_y = TRUE)
  weights <- model_data$weights
  offset <- model_data$offset
  if (family$for_latent) {
    y <- rowMeans(ref_predfun(
      object, excl_offs = FALSE,
      mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)
    ))
    y_oscale <- model_data$y
    if (is.null(family$cats) &&
        (is.factor(y_oscale) || is.character(y_oscale) ||
         is.logical(y_oscale))) {
      stop("If the original (i.e., non-latent) response is `factor`-like, ",
           "`family$cats` must not be `NULL`. See the documentation for ",
           "extend_family()'s argument `latent_y_unqs` to solve this.")
      # Alternatively, we could think about `family$cats <- levels(y_oscale)`.
      # But the error message is conceptually more desirable because it avoids
      # the retrospective modification of extend_family() output.
    }
    if (!is.null(family$cats)) {
      y_oscale <- as.factor(y_oscale)
      stopifnot(nlevels(y_oscale) >= 2)
      if (!identical(levels(y_oscale), family$cats)) {
        stop("The levels of the response variable (after coercing it to a ",
             "`factor`) have to be identical to `family$cats`. See the ",
             "documentation for extend_family()'s argument `latent_y_unqs` to ",
             "solve this.")
      }
    } else if (family$family_oscale == "binomial") {
      if (!all(is_wholenumber(y_oscale))) {
        stop(
          "In projpred, the response must contain numbers of successes (not ",
          "proportions of successes), in contrast to glm() where this is ",
          "possible for a 1-column response if the multiplication with the ",
          "weights gives whole numbers."
        )
      } else if (all(y_oscale %in% c(0, 1)) &&
                 length(response_name) == 1 &&
                 !all(weights == 1)) {
        warning(
          "Assuming that the response contains numbers of successes (not ",
          "proportions of successes), in contrast to glm()."
        )
      }
    }
  } else {
    y <- model_data$y
    y_oscale <- NULL
  }

  # Add (transformed) response under the (possibly) new name:
  data[, response_name] <- y

  target <- get_standard_y(y, weights, family)
  y <- target$y
  weights <- target$weights

  if (family$for_augdat) {
    y <- as.factor(y)
    stopifnot(nlevels(y) >= 2)
    if (!identical(levels(y), family$cats)) {
      stop("The levels of the response variable (after coercing it to a ",
           "`factor`) have to be identical to `family$cats`. See the ",
           "documentation for extend_family()'s argument `augdat_y_unqs` to ",
           "solve this.")
    }
  } else if (family$family == "binomial") {
    if (!all(is_wholenumber(y))) {
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

  if (family$for_augdat && !all(weights == 1)) {
    stop("Currently, the augmented-data projection may not be combined with ",
         "observation weights (other than 1).")
  }
  if (family$for_latent && !all(weights == 1)) {
    stop("Currently, the latent projection may not be combined with ",
         "observation weights (other than 1).")
  }

  if (is.null(offset)) {
    offset <- rep(0, length(y))
  }

  if (!proper_model && !all(offset == 0)) {
    # Disallow offsets for `datafit`s because the submodel fitting does not take
    # offsets into account (but `<refmodel>$mu` contains the observed response
    # values which inevitably "include" the offsets):
    stop("For a `datafit`, offsets are not allowed.")
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
    eta <- ref_predfun(
      object, mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)
    )
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

  # Same as `mu`, but taking offsets into account:
  if (!all(offset == 0)) {
    if (family$family %in% fams_neg_linpred()) {
      eta_offs <- eta - offset
    } else {
      eta_offs <- eta + offset
    }
    mu_offs <- family$linkinv(eta_offs)
    rm(eta_offs)
  } else {
    mu_offs <- mu
  }

  # Miscellaneous -----------------------------------------------------------

  ndraws <- ncol(mu)
  warn_allrandom_dis <- getOption("projpred.warn_allrandom_dis", TRUE)
  if (is.null(dis)) {
    if (family$for_latent && proper_model) {
      if (!is.null(family$link_oscale)) {
        if (family$link_oscale %in% c("probit", "probit_approx")) {
          dis <- rep(1, ndraws)
          warn_allrandom_dis <- FALSE
        } else if (family$link_oscale %in% c("logit", "logistic")) {
          dis <- rep(1.6, ndraws)
          warn_allrandom_dis <- FALSE
        } else {
          dis <- rep(NA, ndraws)
        }
      } else {
        dis <- rep(NA, ndraws)
      }
      if (all(is.na(dis))) {
        message(
          "Since `<refmodel>$dis` will consist of only `NA`s, downstream ",
          "analyses based on this reference model object won't be able to use ",
          "log predictive density (LPD) values on latent scale. Furthermore, ",
          "proj_predict() won't be able to draw from the latent Gaussian ",
          "distribution."
        )
      }
    } else if (!has_dispersion(family)) {
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
  if (getOption("projpred.mlvl_pred_new", FALSE) && warn_allrandom_dis &&
      !all(is.na(dis))) {
    warning("Option `projpred.mlvl_pred_new` has been set to `TRUE`, but the ",
            "reference model includes non-trivial dispersion parameter ",
            "values. Since option `projpred.mlvl_pred_new` also affects the ",
            "projected dispersion parameter values, you need to ensure ",
            "yourself that the reference model's dispersion parameter values ",
            "are the correct ones in the sense that they should typically ",
            "result from integrating out group-level effects. In case of the ",
            "latent projection, a remedy is to switch to response-scale ",
            "analyses as they do not make use of the latent projected ",
            "dispersion parameter values.")
  }

  # Equal weights for the posterior draws by default:
  wdraws_ref <- rep(1 / ndraws, ndraws)

  # Output ------------------------------------------------------------------

  refmodel <- nlist(
    fit = object, formula, div_minimizer, family, eta, mu, mu_offs, dis, y,
    proj_predfun, fetch_data = fetch_data_wrapper, wobs = weights, wdraws_ref,
    offset, cvfun, cvfits, extract_model_data, ref_predfun, cvrefbuilder,
    y_oscale = y_oscale %||% y, nobs = nrow(data)
  )
  if (proper_model) {
    class(refmodel) <- "refmodel"
  } else {
    class(refmodel) <- c("datafit", "refmodel")
  }

  return(refmodel)
}
