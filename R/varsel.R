#' Run search and performance evaluation without cross-validation
#'
#' Run the *search* part and the *evaluation* part for a projection predictive
#' variable selection. The search part determines the predictor ranking (also
#' known as solution path), i.e., the best submodel for each submodel size
#' (number of predictor terms). The evaluation part determines the predictive
#' performance of the submodels along the predictor ranking. A special method is
#' [varsel.vsel()] which re-uses the search results from an earlier
#' [varsel()] (or [cv_varsel()]) run, as illustrated in the main vignette.
#'
#' @param object An object of class `refmodel` (returned by [get_refmodel()] or
#'   [init_refmodel()]) or an object that can be passed to argument `object` of
#'   [get_refmodel()].
#' @param d_test A `list` of the structure outlined in section "Argument
#'   `d_test`" below, providing test data for evaluating the predictive
#'   performance of the submodels as well as of the reference model. If `NULL`,
#'   the training data is used.
#' @param method The method for the search part. Possible options are
#'   `"forward"` for forward search and `"L1"` for L1 search. See also section
#'   "Details" below.
#' @param refit_prj For the evaluation part, should the projections onto the
#' submodels along the predictor ranking be performed again using `ndraws_pred`
#' draws or `nclusters_pred` clusters (`TRUE`) or should their projections from
#' the search part, which used `ndraws` draws or `nclusters` clusters, be re-used
#' (`FALSE`)?
#' @param ndraws Number of posterior draws used in the search part. Ignored if
#'   `nclusters` is not `NULL` or in case of L1 search (because L1 search always
#'   uses a single cluster). If both (`nclusters` and `ndraws`) are `NULL`, the
#'   number of posterior draws from the reference model is used for `ndraws`.
#'   See also section "Details" below.
#' @param nclusters Number of clusters of posterior draws used in the search
#'   part. Ignored in case of L1 search (because L1 search always uses a single
#'   cluster). For the meaning of `NULL`, see argument `ndraws`. See also
#'   section "Details" below.
#' @param ndraws_pred Only relevant if `refit_prj` is `TRUE`. Number of
#'   posterior draws used in the evaluation part. Ignored if `nclusters_pred` is
#'   not `NULL`. If both (`nclusters_pred` and `ndraws_pred`) are `NULL`, the
#'   number of posterior draws from the reference model is used for
#'   `ndraws_pred`. See also section "Details" below.
#' @param nclusters_pred Only relevant if `refit_prj` is `TRUE`. Number of
#'   clusters of posterior draws used in the evaluation part. For the meaning of
#'   `NULL`, see argument `ndraws_pred`. See also section "Details" below.
#' @param nterms_max Maximum submodel size (number of predictor terms) up to
#'   which the search is continued. If `NULL`, then `min(19, D)` is used where
#'   `D` is the number of terms in the reference model (or in `search_terms`, if
#'   supplied). Note that `nterms_max` does not count the intercept, so use
#'   `nterms_max = 0` for the intercept-only model. (Correspondingly, `D` above
#'   does not count the intercept.)
#' @param penalty Only relevant for L1 search. A numeric vector determining the
#'   relative penalties or costs for the predictors. A value of `0` means that
#'   those predictors have no cost and will therefore be selected first, whereas
#'   `Inf` means those predictors will never be selected. If `NULL`, then `1` is
#'   used for each predictor.
#' @param search_control A `list` of "control" arguments (i.e., tuning
#'   parameters) for the search. In case of forward search, these arguments are
#'   passed to the divergence minimizer (see argument `div_minimizer` of
#'   [init_refmodel()] as well as section "Draw-wise divergence minimizers" of
#'   [projpred-package]). In case of forward search, `NULL` causes `...` to be
#'   used not only for the performance evaluation, but also for the search. In
#'   case of L1 search, possible arguments are:
#'   * `lambda_min_ratio`: Ratio between the smallest and largest lambda in the
#'   L1-penalized search (default: `1e-5`). This parameter essentially
#'   determines how long the search is carried out, i.e., how large submodels
#'   are explored. No need to change this unless the program gives a warning
#'   about this.
#'   * `nlambda`: Number of values in the lambda grid for L1-penalized search
#'   (default: `150`). No need to change this unless the program gives a warning
#'   about this.
#'   * `thresh`: Convergence threshold when computing the L1 path (default:
#'   `1e-6`). Usually, there is no need to change this.
#' @param lambda_min_ratio Deprecated (please use `search_control` instead).
#'   Only relevant for L1 search. Ratio between the smallest and largest lambda
#'   in the L1-penalized search. This parameter essentially determines how long
#'   the search is carried out, i.e., how large submodels are explored. No need
#'   to change this unless the program gives a warning about this.
#' @param nlambda Deprecated (please use `search_control` instead). Only
#'   relevant for L1 search. Number of values in the lambda grid for
#'   L1-penalized search. No need to change this unless the program gives a
#'   warning about this.
#' @param thresh Deprecated (please use `search_control` instead). Only relevant
#'   for L1 search. Convergence threshold when computing the L1 path. Usually,
#'   there is no need to change this.
#' @param search_terms Only relevant for forward search. A custom character
#'   vector of predictor term blocks to consider for the search. Section
#'   "Details" below describes more precisely what "predictor term block" means.
#'   The intercept (`"1"`) is always included internally via `union()`, so
#'   there's no difference between including it explicitly or omitting it. The
#'   default `search_terms` considers all the terms in the reference model's
#'   formula.
#' @param search_out Intended for internal use.
#' @param verbose A single logical value indicating whether to print out
#'   additional information during the computations.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. If not
#'   `NA`, then the PRNG state is reset (to the state before calling [varsel()])
#'   upon exiting [varsel()]. Here, `seed` is used for clustering the reference
#'   model's posterior draws (if `!is.null(nclusters)` or
#'   `!is.null(nclusters_pred)`) and for drawing new group-level effects when
#'   predicting from a multilevel submodel (however, not yet in case of a GAMM).
#' @param ... For [varsel.default()]: Arguments passed to [get_refmodel()] as
#'   well as to [varsel.refmodel()]. For [varsel.vsel()]: Arguments passed to
#'   [varsel.refmodel()]. For [varsel.refmodel()]: Arguments passed to the
#'   divergence minimizer (see argument `div_minimizer` of [init_refmodel()] as
#'   well as section "Draw-wise divergence minimizers" of [projpred-package])
#'   when refitting the submodels for the performance evaluation (if `refit_prj`
#'   is `TRUE`).
#'
#' @details
#'
#' # Argument `d_test`
#'
#' If not `NULL`, then `d_test` needs to be a `list` with the following
#' elements:
#' * `data`: a `data.frame` containing the predictor variables for the test set.
#' * `offset`: a numeric vector containing the offset values for the test set
#' (if there is no offset, use a vector of zeros).
#' * `weights`: a numeric vector containing the observation weights for the test
#' set (if there are no observation weights, use a vector of ones).
#' * `y`: a vector or a `factor` containing the response values for the test
#' set. In case of the latent projection, this has to be a vector containing the
#' *latent* response values, but it can also be a vector full of `NA`s if
#' latent-scale post-processing is not needed.
#' * `y_oscale`: Only needs to be provided in case of the latent projection
#' where this needs to be a vector or a `factor` containing the *original*
#' (i.e., non-latent) response values for the test set.
#'
#' @details Arguments `ndraws`, `nclusters`, `nclusters_pred`, and `ndraws_pred`
#'   are automatically truncated at the number of posterior draws in the
#'   reference model (which is `1` for `datafit`s). Using less draws or clusters
#'   in `ndraws`, `nclusters`, `nclusters_pred`, or `ndraws_pred` than posterior
#'   draws in the reference model may result in slightly inaccurate projection
#'   performance. Increasing these arguments affects the computation time
#'   linearly.
#'
#'   For argument `method`, there are some restrictions: For a reference model
#'   with multilevel or additive formula terms or a reference model set up for
#'   the augmented-data projection, only the forward search is available.
#'   Furthermore, argument `search_terms` requires a forward search to take
#'   effect.
#'
#'   L1 search is faster than forward search, but forward search may be more
#'   accurate. Furthermore, forward search may find a sparser model with
#'   comparable performance to that found by L1 search, but it may also
#'   overfit when more predictors are added. This overfit can be detected
#'   by running search validation (see [cv_varsel()]).
#'
#'   An L1 search may select an interaction term before all involved lower-order
#'   interaction terms (including main-effect terms) have been selected. In
#'   \pkg{projpred} versions > 2.6.0, the resulting predictor ranking is
#'   automatically modified so that the lower-order interaction terms come
#'   before this interaction term, but if this is conceptually undesired, choose
#'   the forward search instead.
#'
#'   The elements of the `search_terms` character vector don't need to be
#'   individual predictor terms. Instead, they can be building blocks consisting
#'   of several predictor terms connected by the `+` symbol. To understand how
#'   these building blocks work, it is important to know how \pkg{projpred}'s
#'   forward search works: It starts with an empty vector `chosen` which will
#'   later contain already selected predictor terms. Then, the search iterates
#'   over model sizes \eqn{j \in \{0, ..., J\}}{j = 0, ..., J} (with \eqn{J}
#'   denoting the maximum submodel size, not counting the intercept). The
#'   candidate models at model size \eqn{j} are constructed from those elements
#'   from `search_terms` which yield model size \eqn{j} when combined with the
#'   `chosen` predictor terms. Note that sometimes, there may be no candidate
#'   models for model size \eqn{j}. Also note that internally, `search_terms` is
#'   expanded to include the intercept (`"1"`), so the first step of the search
#'   (model size 0) always consists of the intercept-only model as the only
#'   candidate.
#'
#'   As a `search_terms` example, consider a reference model with formula `y ~
#'   x1 + x2 + x3`. Then, to ensure that `x1` is always included in the
#'   candidate models, specify `search_terms = c("x1", "x1 + x2", "x1 + x3",
#'   "x1 + x2 + x3")` (or, in a simpler way that leads to the same results,
#'   `search_terms = c("x1", "x1 + x2", "x1 + x3")`, for which helper function
#'   [force_search_terms()] exists). This search would start with `y ~ 1` as the
#'   only candidate at model size 0. At model size 1, `y ~ x1` would be the only
#'   candidate. At model size 2, `y ~ x1 + x2` and `y ~ x1 + x3` would be the
#'   two candidates. At the last model size of 3, `y ~ x1 + x2 + x3` would be
#'   the only candidate. As another example, to exclude `x1` from the search,
#'   specify `search_terms = c("x2", "x3", "x2 + x3")` (or, in a simpler way
#'   that leads to the same results, `search_terms = c("x2", "x3")`).
#'
#' @return An object of class `vsel`. The elements of this object are not meant
#'   to be accessed directly but instead via helper functions (see the main
#'   vignette and [projpred-package]).
#'
#' @seealso [cv_varsel()]
#'
#' @examplesIf requireNamespace("rstanarm", quietly = TRUE)
#' # Data:
#' dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)
#'
#' # The `stanreg` fit which will be used as the reference model (with small
#' # values for `chains` and `iter`, but only for technical reasons in this
#' # example; this is not recommended in general):
#' fit <- rstanarm::stan_glm(
#'   y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
#'   QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
#' )
#'
#' # Run varsel() (here without cross-validation, with L1 search, and with small
#' # values for `nterms_max` and `nclusters_pred`, but only for the sake of
#' # speed in this example; this is not recommended in general):
#' vs <- varsel(fit, method = "L1", nterms_max = 3, nclusters_pred = 10,
#'              seed = 5555)
#' # Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
#' # and `?ranking` for possible post-processing functions.
#'
#' @export
varsel <- function(object, ...) {
  UseMethod("varsel")
}

#' @rdname varsel
#' @export
varsel.default <- function(object, ...) {
  refmodel <- get_refmodel(object, ...)
  return(varsel(refmodel, ...))
}

#' @rdname varsel
#' @export
varsel.vsel <- function(object, ...) {
  arg_nms_internal <- c("method", "ndraws", "nclusters", "nterms_max",
                        "search_control", "penalty", "search_terms")
  arg_nms_internal_used <- intersect(arg_nms_internal, ...names())
  n_arg_nms_internal_used <- length(arg_nms_internal_used)
  if (n_arg_nms_internal_used > 0) {
    stop("Argument", if (n_arg_nms_internal_used > 1) "s" else "", " ",
         paste(paste0("`", arg_nms_internal_used, "`"), collapse = ", "), " ",
         "cannot be specified in this case because varsel.vsel() specifies ",
         if (n_arg_nms_internal_used > 1) "them" else "it", " ", "internally.")
  }
  return(varsel(
    object = get_refmodel(object),
    method = object[["args_search"]][["method"]],
    ndraws = object[["args_search"]][["ndraws"]],
    nclusters = object[["args_search"]][["nclusters"]],
    nterms_max = object[["args_search"]][["nterms_max"]],
    search_control = object[["args_search"]][["search_control"]],
    penalty = object[["args_search"]][["penalty"]],
    search_terms = object[["args_search"]][["search_terms"]],
    search_out = list(search_path = object[["search_path"]]),
    ...
  ))
}

#' @rdname varsel
#' @export
varsel.refmodel <- function(
    object,
    d_test = NULL,
    method = "forward",
    ndraws = NULL,
    nclusters = 20,
    ndraws_pred = 400,
    nclusters_pred = NULL,
    refit_prj = !inherits(object, "datafit"),
    nterms_max = NULL,
    verbose = getOption("projpred.verbose", interactive()),
    search_control = NULL,
    lambda_min_ratio = 1e-5,
    nlambda = 150,
    thresh = 1e-6,
    penalty = NULL,
    search_terms = NULL,
    search_out = NULL,
    seed = NA,
    ...
) {
  if (!missing(lambda_min_ratio)) {
    warning("Argument `lambda_min_ratio` is deprecated. Please specify ",
            "control arguments for the search via argument `search_control`. ",
            "Now using `lambda_min_ratio` as element `lambda_min_ratio` of ",
            "`search_control`.")
    search_control$lambda_min_ratio <- lambda_min_ratio
  }
  if (!missing(nlambda)) {
    warning("Argument `nlambda` is deprecated. Please specify control ",
            "arguments for the search via argument `search_control`. ",
            "Now using `nlambda` as element `nlambda` of `search_control`.")
    search_control$nlambda <- nlambda
  }
  if (!missing(thresh)) {
    warning("Argument `thresh` is deprecated. Please specify control ",
            "arguments for the search via argument `search_control`. ",
            "Now using `thresh` as element `thresh` of `search_control`.")
    search_control$thresh <- thresh
  }

  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  if (!is.na(seed)) {
    # Set seed, but ensure the old RNG state is restored on exit:
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(seed)
  }

  refmodel <- object
  nterms_all <- count_terms_in_formula(refmodel$formula) - 1L

  # Parse arguments:
  args <- parse_args_varsel(
    refmodel = refmodel, method = method, refit_prj = refit_prj,
    nterms_max = nterms_max, nclusters = nclusters, search_terms = search_terms,
    nterms_all = nterms_all
  )
  method <- args$method
  refit_prj <- args$refit_prj
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  search_terms <- args$search_terms
  search_terms_was_null <- args$search_terms_was_null

  # Pre-process `d_test`:
  if (is.null(d_test)) {
    d_test <- list(type = "train", data = NULL, offset = refmodel$offset,
                   weights = refmodel$wobs, y = refmodel$y,
                   y_oscale = refmodel$y_oscale)
  } else {
    d_test$type <- "test_hold-out"
    if (!refmodel$family$for_latent) {
      d_test$y_oscale <- d_test$y
    }
    d_test <- d_test[nms_d_test()]
    invisible(lapply(setNames(nm = setdiff(nms_d_test(), c("y", "y_oscale"))),
                     function(d_nm) na.fail(d_test[[d_nm]])))
    hasNA_y_test <- is.na(d_test[["y"]])
    if (any(hasNA_y_test)) {
      stopifnot(all(hasNA_y_test))
    }
    hasNA_y_oscale_test <- is.na(d_test[["y_oscale"]])
    if (any(hasNA_y_oscale_test)) {
      stopifnot(all(hasNA_y_oscale_test))
    }
    if (length(d_test[["weights"]]) != nrow(d_test[["data"]])) {
      stop("Element `d_test$weights` needs to have length equal to the number ",
           "of test observations.")
    }
    if (length(d_test[["offset"]]) != nrow(d_test[["data"]])) {
      stop("Element `d_test$offset` needs to have length equal to the number ",
           "of test observations.")
    }
    if (refmodel$family$for_augdat) {
      d_test$y <- as.factor(d_test$y)
      if (!all(levels(d_test$y) %in% refmodel$family$cats)) {
        stop("The levels of the response variable (after coercing it to a ",
             "`factor`) have to be a subset of `family$cats`. Either modify ",
             "`d_test$y` accordingly or see the documentation for ",
             "extend_family()'s argument `augdat_y_unqs` to solve this.")
      }
      # Re-assign the original levels because some levels might be missing:
      d_test$y <- factor(d_test$y, levels = refmodel$family$cats)
    } else if (refmodel$family$for_latent) {
      if (is.null(refmodel$family$cats) &&
          (is.factor(d_test$y_oscale) || is.character(d_test$y_oscale) ||
           is.logical(d_test$y_oscale))) {
        stop("If the original (i.e., non-latent) response is `factor`-like, ",
             "`family$cats` must not be `NULL`. See the documentation for ",
             "extend_family()'s argument `latent_y_unqs` to solve this.")
      }
      if (!is.null(refmodel$family$cats)) {
        d_test$y_oscale <- as.factor(d_test$y_oscale)
        if (!all(levels(d_test$y_oscale) %in% refmodel$family$cats)) {
          stop("The levels of the response variable (after coercing it to a ",
               "`factor`) have to be a subset of `family$cats`. Either modify ",
               "`d_test$y_oscale` accordingly or see the documentation for ",
               "extend_family()'s argument `latent_y_unqs` to solve this.")
        }
        # Re-assign the original levels because some levels might be missing:
        d_test$y_oscale <- factor(d_test$y_oscale,
                                  levels = refmodel$family$cats)
      }
    }
  }
  y_wobs_test <- setNames(
    as.data.frame(d_test[nms_y_wobs_test(wobs_nm = "weights")]),
    nms_y_wobs_test()
  )
  nobs_test <- nrow(y_wobs_test)

  # Run the search:
  if (!is.null(search_out)) {
    search_path <- search_out[["search_path"]]
  } else {
    search_path <- .select(
      refmodel = refmodel, ndraws = ndraws, nclusters = nclusters,
      method = method, nterms_max = nterms_max, penalty = penalty,
      verbose = verbose, search_control = search_control,
      search_terms = search_terms,
      search_terms_was_null = search_terms_was_null, ...
    )
  }

  # "Run" the performance evaluation for the submodels along the predictor
  # ranking (in fact, we only prepare the performance evaluation by computing
  # precursor quantities, but for users, this difference is not perceivable):
  perf_eval_out <- perf_eval(
    search_path = search_path, refmodel = refmodel, refit_prj = refit_prj,
    ndraws = ndraws_pred, nclusters = nclusters_pred, indices_test = NULL,
    newdata_test = d_test$data, offset_test = d_test$offset,
    wobs_test = d_test$weights, y_test = d_test$y,
    y_oscale_test = d_test$y_oscale, verbose = verbose, ...
  )

  # Predictive performance of the reference model:
  if (inherits(refmodel, "datafit")) {
    # In this case, there is no actual reference model, so we don't know how to
    # predict for actual new data.
    ref <- list(mu = rep(NA, nobs_test), lppd = rep(NA, nobs_test))
    if (refmodel$family$for_latent) {
      # In general, we could use `ref$oscale <- ref` here, but the case where
      # refmodel$family$latent_ilink() returns a 3-dimensional array (S x N x C)
      # needs special care.
      if (!is.null(refmodel$family$cats)) {
        mu_oscale <- structure(rep(NA,
                                   nobs_test * length(refmodel$family$cats)),
                               ndiscrete = length(refmodel$family$cats),
                               class = "augvec")
      } else {
        mu_oscale <- ref$mu
      }
      ref$oscale <- list(mu = mu_oscale, lppd = ref$lppd)
    }
  } else {
    if (d_test$type == "train") {
      if (formula_contains_group_terms(refmodel$formula) &&
          getOption("projpred.mlvl_pred_new", FALSE)) {
        # Need to use `mlvl_allrandom = TRUE` (`refmodel$mu_offs` is based on
        # `mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)`):
        eta_test <- refmodel$ref_predfun(refmodel$fit, excl_offs = FALSE)
        mu_test <- refmodel$family$linkinv(eta_test)
      } else {
        mu_test <- refmodel$mu_offs
      }
    } else {
      newdata_for_ref <- d_test$data
      if (inherits(refmodel$fit, "stanreg") &&
          length(refmodel$fit$offset) > 0) {
        if ("projpred_internal_offs_stanreg" %in% names(newdata_for_ref)) {
          stop("Need to write to column `projpred_internal_offs_stanreg` of ",
               "`d_test$data`, but that column already exists. Please rename ",
               "this column in `d_test$data` and try again.")
        }
        newdata_for_ref$projpred_internal_offs_stanreg <- d_test$offset
      }
      eta_test <- refmodel$ref_predfun(refmodel$fit, newdata = newdata_for_ref,
                                       excl_offs = FALSE)
      mu_test <- refmodel$family$linkinv(eta_test)
    }
    ref <- weighted_summary_means(
      y_wobs_test = y_wobs_test, family = refmodel$family,
      wdraws = refmodel$wdraws_ref, mu = mu_test, dis = refmodel$dis,
      cl_ref = seq_along(refmodel$wdraws_ref)
    )
  }

  # The object to be returned:
  vs <- nlist(refmodel,
              nobs_train = refmodel$nobs,
              search_path,
              predictor_ranking = search_path$predictor_ranking,
              predictor_ranking_cv = NULL,
              ce = perf_eval_out[["ce"]],
              type_test = d_test$type,
              y_wobs_test,
              nobs_test,
              summaries = nlist(sub = perf_eval_out[["sub_summaries"]], ref),
              summaries_fast = NULL,
              nterms_all,
              nterms_max,
              method,
              cv_method = NULL,
              nloo = NULL,
              loo_inds = NULL,
              K = NULL,
              validate_search = NULL,
              ### Not set to `NULL` because in K-fold CV (relevant when using
              ### cv_varsel.vsel() based on this `vsel` object), `cvfits = NULL`
              ### would mean to use `cvfun`; instead, the default should be
              ### element `cvfits` from `refmodel`, so:
              cvfits = refmodel$cvfits,
              ###
              args_search = nlist(
                method, ndraws, nclusters, nterms_max,
                search_control = if (
                  method == "forward" && is.null(search_control)
                ) list(...) else search_control,
                penalty,
                search_terms = if (search_terms_was_null) NULL else search_terms
              ),
              clust_used_search = search_path$p_sel$clust_used,
              clust_used_eval = perf_eval_out[["clust_used"]],
              nprjdraws_search = search_path$p_sel$nprjdraws,
              nprjdraws_eval = perf_eval_out[["nprjdraws"]],
              refit_prj,
              projpred_version = utils::packageVersion("projpred"))
  class(vs) <- "vsel"

  return(vs)
}

# Workhorse function for the search
#
# @param reweighting_args If the projected draws (i.e., those obtained from
#   clustering or thinning) should be re-weighted (usually according to PSIS
#   weights), then this needs to be a `list` with elements `wdraws_ref` and
#   `cl_ref`. For these two elements, see the (internal) documentation of
#   weighted_summary_means().
# @param verbose_txt_obs Passed to `...` of verb_out(), so character string(s)
#   to be included in the verbose message indicating the start of the search.
#   May also be `NULL` to omit that verbose message completely.
# For all other arguments, see the documentation of varsel().
#
# @return A list with elements `predictor_ranking` (the predictor ranking
#   resulting from the search, see `?ranking` for a formal definition),
#   `outdmins` (the submodel fits along the predictor ranking, with the number
#   of fits per model size being equal to the number of projected draws), and
#   `p_sel` (the output from get_refdist() for the search).
.select <- function(refmodel, ndraws, nclusters, reweighting_args = NULL,
                    method, nterms_max, penalty, verbose, verbose_txt_obs = "",
                    search_control, ...) {
  if (is.null(reweighting_args)) {
    p_sel <- get_refdist(refmodel, ndraws = ndraws, nclusters = nclusters)
  } else {
    # Reweight the clusters (or thinned draws) according to the PSIS weights:
    p_sel <- get_p_clust(
      family = refmodel$family, eta = refmodel$eta, mu = refmodel$mu,
      mu_offs = refmodel$mu_offs, dis = refmodel$dis,
      wdraws = reweighting_args$wdraws_ref, cl = reweighting_args$cl_ref,
      clust_used_forced = clust_info(
        ndraws = ndraws,
        nclusters = nclusters,
        S = length(refmodel$wdraws_ref)
      )[["clust_used"]]
    )
  }

  if (!is.null(verbose_txt_obs)) {
    verb_out("-----\nRunning ", method, " search ", verbose_txt_obs, "with ",
             txt_clust_draws(p_sel[["clust_used"]], p_sel[["nprjdraws"]]),
             " ...", verbose = verbose)
  }
  if (method == "L1") {
    search_path <- search_L1(
      p_ref = p_sel, refmodel = refmodel, nterms_max = nterms_max,
      penalty = penalty, search_control = search_control
    )
  } else if (method == "forward") {
    search_path <- search_forward(
      p_ref = p_sel, refmodel = refmodel, nterms_max = nterms_max,
      verbose = verbose, search_control = search_control, ...
    )
  }
  if (!is.null(verbose_txt_obs)) {
    verb_out("-----", verbose = verbose)
  }

  search_path$p_sel <- p_sel
  return(search_path)
}

# Auxiliary function for parsing the arguments of varsel()
#
# The arguments specified by the user (or the function calling this function)
# are treated as they are, but if some are not given, then this function fills
# them in with the default values. The purpose of this function is to avoid
# repeating the same code both in varsel() and cv_varsel().
parse_args_varsel <- function(refmodel, method, refit_prj, nterms_max,
                              nclusters, search_terms, nterms_all) {
  search_terms_was_null <- is.null(search_terms)
  if (search_terms_was_null) {
    search_terms <- split_formula(refmodel$formula,
                                  data = refmodel$fetch_data())
  }
  search_terms <- union("1", search_terms)
  has_group_features <- formula_contains_group_terms(refmodel$formula)
  has_additive_features <- formula_contains_additive_terms(refmodel$formula)

  stopifnot(!is.null(method))
  if (method == "l1") {
    method <- toupper(method)
  }
  if (method == "L1") {
    if (has_group_features || has_additive_features) {
      stop("L1 search is only supported for reference models without ",
           "multilevel and without additive (\"smoothing\") terms.")
    }
    if (!search_terms_was_null) {
      warning("Argument `search_terms` only takes effect if ",
              "`method = \"forward\"`.")
    }
    if (refmodel$family$for_augdat) {
      stop("Currently, the augmented-data projection may not be combined ",
           "with an L1 search.")
    }
  }
  if (!method %in% c("L1", "forward")) {
    stop("Unexpected value for argument `method`.")
  }

  stopifnot(!is.null(refit_prj))
  if (refit_prj && inherits(refmodel, "datafit")) {
    warning("For an `object` of class `datafit`, `refit_prj` is automatically ",
            "set to `FALSE`.")
    refit_prj <- FALSE
  }

  if (method == "L1") {
    nclusters <- 1
  }

  if (!has_group_features &&
      getOption("projpred.warn_instable_projections", TRUE) &&
      method == "L1" && refmodel$family$family == "poisson") {
    warning(
      "For non-multilevel Poisson models, an L1 search based on the ",
      "traditional projection may be instable. The latent projection may be a ",
      "remedy. See section \"Troubleshooting\" of the main vignette for more ",
      "information."
    )
  }

  search_terms_unq <- unique(unlist(
    strsplit(search_terms, split = "+", fixed = TRUE)
  ))
  max_nv_possible <- count_terms_chosen(search_terms_unq) - 1L
  if (is.null(nterms_max)) {
    nterms_max <- 19
    if (nterms_max < max_nv_possible &&
        getOption("projpred.mssg_cut_search", TRUE)) {
      message("Cutting off the search at size ", nterms_max, ". Use argument ",
              "`nterms_max` to customize this.")
    }
  }
  nterms_max <- min(max_nv_possible, nterms_max)

  if (nterms_max == nterms_all && has_group_features &&
      getOption("projpred.warn_instable_projections", TRUE) &&
      (refmodel$family$family == "gaussian" || refmodel$family$for_latent)) {
    warning(
      "In case of the Gaussian family (also in case of the latent projection) ",
      "and multilevel terms, the projection onto the full model can be ",
      "instable and even lead to an error, see GitHub issue #323. If you ",
      "experience this and may refrain from the projection onto the full ",
      "model, set `nterms_max` to the number of predictor terms in the full ",
      "model minus 1 (possibly accounting for submodel sizes skipped by ",
      "custom `search_terms`)."
    )
  }

  return(nlist(method, refit_prj, nterms_max, nclusters, search_terms,
               search_terms_was_null))
}
