# General functions for CV ------------------------------------------------

#' Run search and performance evaluation with cross-validation
#'
#' Run the *search* part and the *evaluation* part for a projection predictive
#' variable selection. The search part determines the predictor ranking (also
#' known as solution path), i.e., the best submodel for each submodel size
#' (number of predictor terms). The evaluation part determines the predictive
#' performance of the submodels along the predictor ranking. In contrast to
#' [varsel()], [cv_varsel()] performs a cross-validation (CV) by running the
#' search part with the training data of each CV fold separately (an exception
#' is explained in section "Note" below) and by running the evaluation part on
#' the corresponding test set of each CV fold. A special method is
#' [cv_varsel.vsel()] because it re-uses the search results from an earlier
#' [cv_varsel()] (or [varsel()]) run, as illustrated in the main vignette.
#'
#' @inheritParams varsel
#' @param cv_method The CV method, either `"LOO"` or `"kfold"`. In the `"LOO"`
#'   case, a Pareto-smoothed importance sampling leave-one-out CV (PSIS-LOO-CV)
#'   is performed, which avoids refitting the reference model `nloo` times (in
#'   contrast to a standard LOO-CV). In the `"kfold"` case, a \eqn{K}-fold-CV is
#'   performed. See also section "Note" below.
#' @param nloo Only relevant if `cv_method = "LOO"` and `validate_search =
#'   TRUE`. If `nloo > 0` is smaller than the number of all observations, full
#'   LOO-CV (i.e., PSIS-LOO CV with `validate_search = TRUE` and with `nloo = n`
#'   where `n` denotes the number of all observations) is approximated by
#'   subsampled LOO-CV, i.e., by combining the fast (i.e., `validate_search =
#'   FALSE`) LOO result for the selected models and `nloo` leave-one-out
#'   searches using the difference estimator with simple random sampling (SRS)
#'   without replacement (WOR) (Magnusson et al., 2020). Smaller `nloo` values
#'   lead to faster computation, but higher uncertainty in the evaluation part.
#'   If `NULL`, all observations are used (as by default). Note that performance
#'   statistic `"auc"` (see argument `stats` of [summary.vsel()] and
#'   [plot.vsel()]) is not supported in case of subsampled LOO-CV. Furthermore,
#'   option `"best"` for argument `baseline` of [summary.vsel()] and
#'   [plot.vsel()] is not supported in case of subsampled LOO-CV.
#' @param K Only relevant if `cv_method = "kfold"` and if `cvfits` is `NULL`
#'   (which is the case for reference model objects created by
#'   [get_refmodel.stanreg()] or [brms::get_refmodel.brmsfit()]). Number of
#'   folds in \eqn{K}-fold-CV.
#' @param cvfits Only relevant if `cv_method = "kfold"`. The same as argument
#'   `cvfits` of [init_refmodel()], but repeated here so that output from
#'   [run_cvfun()] can be inserted here straightforwardly.
#' @param validate_search A single logical value indicating whether to
#'   cross-validate also the search part, i.e., whether to run the search
#'   separately for each CV-fold (`TRUE`) or not (`FALSE`). With `FALSE`
#'   the computation is faster, but the predictive performance estimates
#'   of the selected submodels are optimistically biased. However, these fast
#'   biased estimated can be useful to obtain initial information on the
#'   usefulness of projection predictive variable selection.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. If not
#'   `NA`, then the PRNG state is reset (to the state before calling
#'   [cv_varsel()]) upon exiting [cv_varsel()]. Here, `seed` is used for
#'   clustering the reference model's posterior draws (if `!is.null(nclusters)`
#'   or `!is.null(nclusters_pred)`), for subsampling PSIS-LOO-CV folds (if
#'   `nloo` is smaller than the number of observations), for sampling the folds
#'   in \eqn{K}-fold-CV, and for drawing new group-level effects when predicting
#'   from a multilevel submodel (however, not yet in case of a GAMM).
#' @param parallel A single logical value indicating whether to run costly parts
#'   of the CV in parallel (`TRUE`) or not (`FALSE`). See also section "Note"
#'   below as well as section "Parallelization" in [projpred-package].
#' @param ... For [cv_varsel.default()]: Arguments passed to [get_refmodel()] as
#'   well as to [cv_varsel.refmodel()]. For [cv_varsel.vsel()]: Arguments passed
#'   to [cv_varsel.refmodel()]. For [cv_varsel.refmodel()]: Arguments passed to
#'   the divergence minimizer (see argument `div_minimizer` of [init_refmodel()]
#'   as well as section "Draw-wise divergence minimizers" of [projpred-package])
#'   when refitting the submodels for the performance evaluation (if `refit_prj`
#'   is `TRUE`).
#'
#' @inherit varsel details return
#'
#' @note If `validate_search` is `FALSE`, the search is not included in the CV
#'   so that only a single full-data search is run. If the number of
#'   observations is big, the fast PSIS-LOO-CV along the full-data search path
#'   is likely to be accurate. If the number of observations is small or
#'   moderate, the fast PSIS-LOO-CV along the full-data search path is likely to
#'   have optimistic bias in the middle of the search path. This result can be
#'   used to guide further actions and the optimistic bias can be greatly
#'   reduced by using `validate_search = TRUE`.
#'
#'   PSIS uses Pareto-\eqn{\hat{k}} diagnostic to assess the reliability of
#'   PSIS-LOO-CV. Global option `projpred.warn_psis` (default `TRUE`) controls
#'   whether the Pareto-\eqn{\hat{k}} diagnostics may result in warnings. See
#'   [loo::loo-glossary] for how to interpret the Pareto-\eqn{\hat{k}} values
#'   and the warning thresholds. \pkg{projpred} does not support the usually
#'   recommended moment-matching (see [loo::loo_moment_match()] and
#'   [brms::loo_moment_match()]), mixture importance sampling
#'   (`vignette("loo2-mixis", package="loo")`), or `reloo`-ing
#'   ([brms::reloo()]). If the reference model PSIS-LOO-CV Pareto-\eqn{\hat{k}}
#'   values are good, but there are high Pareto-\eqn{\hat{k}} values for the
#'   projected models, you can try increasing the number of draws used for the
#'   PSIS-LOO-CV (`ndraws` in case of `refit_prj = FALSE`; `ndraws_pred` in case
#'   of `refit_prj = TRUE`). If increasing the number of draws does not help and
#'   if the reference model PSIS-LOO-CV Pareto-\eqn{\hat{k}} values are high,
#'   and the reference model PSIS-LOO-CV results change substantially when using
#'   moment-matching, mixture importance sampling, or `reloo`-ing, we recommend
#'   to use \eqn{K}-fold-CV within `projpred`.
#'
#'   For PSIS-LOO-CV, \pkg{projpred} calls [loo::psis()] (or, exceptionally,
#'   [loo::sis()], see below) with `r_eff = NA`. This is only a problem if there
#'   was extreme autocorrelation between the MCMC iterations when the reference
#'   model was built. In those cases however, the reference model should not
#'   have been used anyway, so we don't expect \pkg{projpred}'s `r_eff = NA` to
#'   be a problem.
#'
#'   PSIS cannot be used if the number of draws or clusters is too small. In
#'   such cases, \pkg{projpred} resorts to standard importance sampling (SIS)
#'   and shows a message about this. Throughout the documentation, the term
#'   "PSIS" is used even though in fact, \pkg{projpred} resorts to SIS in these
#'   special cases. If SIS is used, check that the reference model PSIS-LOO-CV
#'   Pareto-\eqn{\hat{k}} values are good.
#'
#'   With `parallel = TRUE`, costly parts of \pkg{projpred}'s CV can be run in
#'   parallel. Costly parts are the fold-wise searches and performance
#'   evaluations in case of `validate_search = TRUE`. (Note that in case of
#'   \eqn{K}-fold CV, the \eqn{K} reference model refits are not affected by
#'   argument `parallel`; only \pkg{projpred}'s CV is affected.) The
#'   parallelization is powered by the \pkg{foreach} package. Thus, any parallel
#'   (or sequential) backend compatible with \pkg{foreach} can be used, e.g.,
#'   the backends from packages \pkg{doParallel}, \pkg{doMPI}, or
#'   \pkg{doFuture}. For GLMs, this CV parallelization should work reliably, but
#'   for other models (such as GLMMs), it may lead to excessive memory usage
#'   which in turn may crash the R session (on Unix systems, setting an
#'   appropriate memory limit via [unix::rlimit_as()] may avoid crashing the
#'   whole machine). However, the problem of excessive memory usage is less
#'   pronounced for the CV parallelization than for the projection
#'   parallelization described in [projpred-package]. In that regard, the CV
#'   parallelization is recommended over the projection parallelization.
#'
#' @references
#'
#' Måns Magnusson, Michael Riis Andersen, Johan Jonasson, Aki Vehtari
#' (2020). Leave-one-out cross-validation for Bayesian model
#' comparison in large data. Proceedings of the 23rd International
#' Conference on Artificial Intelligence and Statistics (AISTATS),
#' PMLR 108:341-351.
#'
#' Aki Vehtari, Andrew Gelman, and Jonah Gabry (2017). Practical Bayesian Model
#' Evaluation Using Leave-One-Out Cross-Validation and WAIC. Statistics and
#' Computing, 27(5):1413--32. \doi{10.1007/s11222-016-9696-4}.
#'
#' Aki Vehtari, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah
#' Gabry (2024). Pareto smoothed importance sampling. Journal of
#' Machine Learning Research, 25(72):1-58.
#'
#' @seealso [varsel()]
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
#'   QR = TRUE, chains = 2, iter = 1000, refresh = 0, seed = 9876
#' )
#'
#' # Run cv_varsel() (with L1 search and small values for `K`, `nterms_max`, and
#' # `nclusters_pred`, but only for the sake of speed in this example; this is
#' # not recommended in general):
#' cvvs <- cv_varsel(fit, method = "L1", cv_method = "kfold", K = 2,
#'                   nterms_max = 3, nclusters_pred = 10, seed = 5555)
#' # Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
#' # and `?ranking` for possible post-processing functions.
#'
#' @export
cv_varsel <- function(object, ...) {
  UseMethod("cv_varsel")
}

#' @rdname cv_varsel
#' @export
cv_varsel.default <- function(object, ...) {
  refmodel <- get_refmodel(object, ...)
  return(cv_varsel(refmodel, ...))
}

#' @rdname cv_varsel
#' @export
cv_varsel.vsel <- function(
    object,
    cv_method = object$cv_method %||% "LOO",
    nloo = object$nloo,
    K = object$K %||% if (!inherits(object, "datafit")) 5 else 10,
    cvfits = object$cvfits,
    validate_search = object$validate_search %||% TRUE,
    ...
) {
  ## the following arguments should not change
  arg_nms_internal <- c("method", "ndraws", "nclusters", "nterms_max",
                        "search_control", "penalty", "search_terms")
  dots <- list(...)
  arg_nms_internal_used <- intersect(arg_nms_internal, names(dots))
  for (arg in arg_nms_internal_used) {
    if (!identical(object[["args_search"]][[arg]], dots[[arg]])) {
      message("Argument `", arg, "` ignored. Using the argument value stored ",
              "in the `vsel` object.")
    }
    ## remove duplicate arguments
    dots[[arg]] <- NULL
  }

  refmodel <- get_refmodel(object)
  rk_foldwise <- ranking(object)[["foldwise"]]
  if (validate_search && !is.null(rk_foldwise)) {
    if (!identical(cv_method, object[["cv_method"]]) ||
        (identical(cv_method, object[["cv_method"]]) &&
         identical(cv_method, "kfold") &&
         (is.null(cvfits) || !identical(cvfits, object[["cvfits"]]))) ||
        (identical(cv_method, object[["cv_method"]]) &&
         (identical(cv_method, "LOO") || identical(cv_method, "loo")) &&
         !identical(nloo, refmodel[["nobs"]]))) {
      # In these cases, previous fold-wise predictor rankings cannot be re-used
      # for the `validate_search = TRUE` run requested here:
      message("In this case, the previous fold-wise search results cannot be ",
              "re-used, so the fold-wise searches are run again.")
      rk_foldwise <- NULL
    }
    if (identical(cv_method, object[["cv_method"]]) &&
        identical(cv_method, "kfold") &&
        identical(cvfits, object[["cvfits"]]) &&
        inherits(refmodel[["fit"]], "brmsfit") &&
        getOption("projpred.mlvl_proj_ref_new", FALSE) &&
        formula_contains_group_terms(refmodel[["formula"]])) {
      # In this case, the call(s) to ref_predfun() that is/are performed when
      # initializing the fold-wise reference model objects via init_refmodel()
      # (within cvrefbuilder()) involve(s) using the PRNG, so in order to be
      # able to re-use previous fold-wise predictor rankings, argument
      # `brms_seed` of brms:::get_refmodel.brmsfit() needs to be set:
      warning("Please make sure that you have set argument `brms_seed` of ",
              "brms:::get_refmodel.brmsfit() to some non-`NULL` value.")
    }
  }

  return(do_call(cv_varsel, c(
    list(
      object = refmodel,
      method = object[["args_search"]][["method"]],
      ndraws = object[["args_search"]][["ndraws"]],
      nclusters = object[["args_search"]][["nclusters"]],
      nterms_max = object[["args_search"]][["nterms_max"]],
      search_control = object[["args_search"]][["search_control"]],
      penalty = object[["args_search"]][["penalty"]],
      search_terms = object[["args_search"]][["search_terms"]],
      cv_method = cv_method,
      nloo = nloo,
      K = K,
      cvfits = cvfits,
      validate_search = validate_search,
      search_out = nlist(search_path = object[["search_path"]], rk_foldwise)
    ),
    dots
  )))
}

#' @rdname cv_varsel
#' @export
cv_varsel.refmodel <- function(
    object,
    method = "forward",
    cv_method = if (!inherits(object, "datafit")) "LOO" else "kfold",
    ndraws = NULL,
    nclusters = 20,
    ndraws_pred = 400,
    nclusters_pred = NULL,
    refit_prj = !inherits(object, "datafit"),
    nterms_max = NULL,
    penalty = NULL,
    verbose = getOption("projpred.verbose", as.integer(interactive())),
    nloo = if (cv_method == "LOO") object$nobs else NULL,
    K = if (!inherits(object, "datafit")) 5 else 10,
    cvfits = object$cvfits,
    search_control = NULL,
    lambda_min_ratio = 1e-5,
    nlambda = 150,
    thresh = 1e-6,
    validate_search = TRUE,
    seed = NA,
    search_terms = NULL,
    search_out = NULL,
    parallel = getOption("projpred.parallel_cv", FALSE),
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
  if (missing(parallel) && is.null(getOption("projpred.parallel_cv")) &&
      !is.null(getOption("projpred.prll_cv"))) {
    warning(
      "Global option `projpred.prll_cv` is deprecated. Please use global ",
      "option `projpred.parallel_cv` instead (or use argument `parallel` ",
      "directly). Now using the value from global option `projpred.prll_cv` ",
      "for argument `parallel`."
    )
    parallel <- getOption("projpred.prll_cv")
  }
  verbose <- verbose_from_deprecated_options(verbose, with_cv = TRUE)

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

  # Parse arguments which also exist in varsel():
  args <- parse_args_varsel(
    refmodel = refmodel, method = method, refit_prj = refit_prj,
    nterms_max = nterms_max, nclusters = nclusters, search_terms = search_terms,
    nterms_all = nterms_all, verbose = verbose
  )
  method <- args$method
  refit_prj <- args$refit_prj
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  search_terms <- args$search_terms
  search_terms_was_null <- args$search_terms_was_null
  verbose <- args$verbose
  # Parse arguments specific to cv_varsel():
  args <- parse_args_cv_varsel(
    refmodel = refmodel, cv_method = cv_method, nloo = nloo, K = K,
    cvfits = cvfits, validate_search = validate_search, refit_prj = refit_prj,
    search_out = search_out
  )
  cv_method <- args$cv_method
  nloo <- args$nloo
  K <- args$K
  cvfits <- args$cvfits

  # Full-data search:
  if (!is.null(search_out)) {
    search_path_fulldata <- search_out[["search_path"]]
  } else {
    search_path_fulldata <- .select(
      refmodel = refmodel, ndraws = ndraws, nclusters = nclusters,
      method = method, nterms_max = nterms_max, penalty = penalty,
      verbose = if ((validate_search || cv_method == "kfold") && verbose >= 2L) {
        verbose - 1L
      } else {
        verbose
      },
      # NOTE: If `!validate_search`, this is still a full-data search, but in
      # that case, there are no fold-wise searches, so declaring this as a
      # full-data search could be confusing:
      verbose_txt_add = if (validate_search) "using the full dataset " else "",
      search_control = search_control, search_terms = search_terms,
      search_terms_was_null = search_terms_was_null, ...
    )
  }

  if (!is.null(search_out) && validate_search) {
    # Extract the fold-wise predictor rankings (to avoid passing the large
    # object `search_out` itself) and coerce them to a `list` (in a row-wise
    # manner) which is needed for the K-fold-CV parallelization:
    search_out_rks <- search_out[["rk_foldwise"]]
    if (!is.null(search_out_rks)) {
      n_folds <- nrow(search_out_rks)
      search_out_rks <- lapply(seq_len(n_folds), function(row_idx) {
        search_out_rks[row_idx, ]
      })
    }
  } else {
    search_out_rks <- NULL
  }

  summaries_fast <- NULL
  if (cv_method == "LOO") {
    sel_cv <- loo_varsel(
      refmodel = refmodel, method = method, nterms_max = nterms_max,
      ndraws = ndraws, nclusters = nclusters, ndraws_pred = ndraws_pred,
      nclusters_pred = nclusters_pred, refit_prj = refit_prj, penalty = penalty,
      verbose = verbose, search_control = search_control, nloo = nloo,
      validate_search = validate_search,
      search_path_fulldata = if (validate_search) {
        # Not needed in this case, so for computational efficiency, avoiding
        # passing the large object `search_path_fulldata` to loo_varsel():
        NULL
      } else {
        search_path_fulldata
      },
      search_terms = search_terms,
      search_terms_was_null = search_terms_was_null,
      search_out_rks = search_out_rks, parallel = parallel, ...
    )
    if (validate_search && nloo < refmodel$nobs) {
      # Run fast LOO-CV to be used in subsampling difference estimator
      summaries_fast <- loo_varsel(
        refmodel = refmodel, method = method, nterms_max = nterms_max,
        ndraws = ndraws, nclusters = nclusters, ndraws_pred = ndraws_pred,
        nclusters_pred = nclusters_pred, refit_prj = refit_prj, penalty = penalty,
        verbose = if (verbose >= 2L) verbose - 1L else verbose,
        verbose_txt_add = "using the full dataset ",
        search_control = search_control,
        nloo = refmodel$nobs,    # fast LOO-CV (using all observations)
        validate_search = FALSE, # fast LOO-CV (using all observations)
        search_path_fulldata = search_path_fulldata,
        search_terms = search_terms,
        search_terms_was_null = search_terms_was_null,
        search_out_rks = search_out_rks, parallel = parallel, ...
      )[["summaries"]]
    }
  } else if (cv_method == "kfold") {
    sel_cv <- kfold_varsel(
      refmodel = refmodel, method = method, nterms_max = nterms_max,
      ndraws = ndraws, nclusters = nclusters, ndraws_pred = ndraws_pred,
      nclusters_pred = nclusters_pred, refit_prj = refit_prj, penalty = penalty,
      verbose = verbose, search_control = search_control, K = K,
      cvfits = cvfits, validate_search = validate_search,
      search_path_fulldata = if (validate_search) {
        # Not needed in this case, so for computational efficiency, avoiding
        # passing the large object `search_path_fulldata` to loo_varsel():
        NULL
      } else {
        # For K-fold-CV, `validate_search = FALSE` may not be combined with
        # `refit_prj = FALSE`, so element `predictor_ranking` is all we need:
        search_path_fulldata["predictor_ranking"]
      },
      search_terms = search_terms, search_out_rks = search_out_rks,
      parallel = parallel, ...
    )
  }

  if (!validate_search && cv_method == "LOO") {
    ce_out <- sel_cv$ce
  } else {
    ce_out <- rep(NA_real_, length(search_path_fulldata$predictor_ranking) + 1L)
  }

  # Defined here for `nobs_test` later:
  y_wobs_test <- sel_cv$y_wobs_test

  # Information about the clustering/thinning used for the search:
  refdist_info_search <- search_path_fulldata$p_sel[c("clust_used",
                                                      "nprjdraws")]
  # Information about the clustering/thinning used for the performance
  # evaluation:
  if (refit_prj) {
    refdist_info_eval <- sel_cv[c("clust_used_eval", "nprjdraws_eval")]
  } else {
    refdist_info_eval <- refdist_info_search
  }

  # The object to be returned:
  vs <- nlist(refmodel,
              nobs_train = refmodel$nobs,
              search_path = search_path_fulldata,
              predictor_ranking = search_path_fulldata$predictor_ranking,
              predictor_ranking_cv = sel_cv$predictor_ranking_cv,
              ce = ce_out,
              type_test = cv_method,
              y_wobs_test,
              nobs_test = nrow(y_wobs_test),
              summaries = sel_cv$summaries,
              summaries_fast,
              nterms_all,
              nterms_max,
              method,
              cv_method,
              nloo,
              loo_inds = sel_cv$inds,
              K,
              validate_search,
              cvfits,
              args_search = nlist(
                method, ndraws, nclusters, nterms_max,
                search_control = if (
                  method == "forward" && is.null(search_control)
                ) list(...) else search_control,
                penalty,
                search_terms = if (search_terms_was_null) NULL else search_terms
              ),
              clust_used_search = refdist_info_search$clust_used,
              clust_used_eval = refdist_info_eval$clust_used,
              nprjdraws_search = refdist_info_search$nprjdraws,
              nprjdraws_eval = refdist_info_eval$nprjdraws,
              refit_prj,
              projpred_version = utils::packageVersion("projpred"))
  class(vs) <- "vsel"
  return(vs)
}

# Auxiliary function for parsing the arguments specific to cv_varsel()
#
# This is similar in spirit to parse_args_varsel(), in that it prevents the main
# function from becoming too long and complicated to maintain.
#
# @param refmodel See argument `object` of cv_varsel().
# @param cv_method See argument `cv_method` of cv_varsel().
# @param nloo See argument `nloo` of cv_varsel().
# @param K See argument `K` of cv_varsel().
# @param cvfits See argument `cvfits` of cv_varsel().
# @param validate_search See argument `validate_search` of cv_varsel().
#
# @return A list with the processed elements `cv_method`, `nloo`, `K`, and
#   `cvfits`.
parse_args_cv_varsel <- function(refmodel, cv_method, nloo, K, cvfits,
                                 validate_search, refit_prj, search_out) {
  stopifnot(!is.null(cv_method))
  if (cv_method == "loo") {
    cv_method <- toupper(cv_method)
  }
  if (!cv_method %in% c("kfold", "LOO")) {
    stop("Unknown `cv_method`.")
  }
  if (cv_method == "LOO" && inherits(refmodel, "datafit")) {
    warning("For an `object` of class `datafit`, `cv_method` is automatically ",
            "set to \"kfold\".")
    cv_method <- "kfold"
  }

  if (cv_method == "kfold") {
    if (!is.null(cvfits)) {
      if (identical(names(cvfits), "fits")) {
        warning(
          "The content of `cvfits`'s sub-list called `fits` should be moved ",
          "one level up (and element `fits` removed). The old structure will ",
          "continue to work for a while, but is deprecated."
        )
        cvfits <- structure(cvfits$fits, folds = attr(cvfits, "folds"))
      }
      K <- length(cvfits)
    }
    stopifnot(!is.null(K))
    if (length(K) > 1 || !is.numeric(K) || !is_wholenumber(K)) {
      stop("`K` must be a single integer value.")
    }
    if (K < 2) {
      stop("`K` must be at least 2.")
    }
    if (K > NROW(refmodel$y)) {
      stop("`K` cannot exceed the number of observations.")
    }
    if (!validate_search && !refit_prj) {
      # Not allowed because this would induce a dependency between training and
      # test data:
      stop("For K-fold-CV, `validate_search = FALSE` may not be combined with ",
           "`refit_prj = FALSE`.")
    }
    nloo <- NULL
  } else {
    stopifnot(!is.null(refmodel[["nobs"]]))
    nloo <- min(nloo, refmodel[["nobs"]])
    if (nloo < 1) {
      stop("nloo must be at least 1")
    }
    if (!validate_search && nloo < refmodel[["nobs"]]) {
      stop("Subsampled PSIS-LOO-CV is not supported for ",
           "`validate_search = FALSE`.")
    }
  }

  # Restrictions in case of previous search results which should be re-used:
  if (!is.null(search_out)) {
    if (validate_search && !is.null(search_out[["rk_foldwise"]]) &&
        !refit_prj) {
      # In this case, we would need the fold-wise submodel fits (along the
      # fold-wise predictor rankings), which are currently not available:
      stop("If `validate_search = TRUE`, then in general, `refit_prj = FALSE` ",
           "cannot be combined with the re-use of previous search results.")
    }
  }

  return(nlist(cv_method, nloo, K, cvfits))
}

# PSIS-LOO-CV -------------------------------------------------------------

# Workhorse function for a variable selection with PSIS-LOO-CV
#
# Argument `validate_search` indicates whether the search is performed
# separately for each LOO-CV fold (i.e., separately for each observation). For
# all other arguments, see the documentation of cv_varsel().
loo_varsel <- function(refmodel, method, nterms_max, ndraws,
                       nclusters, ndraws_pred, nclusters_pred, refit_prj,
                       penalty, verbose, search_control, nloo, validate_search,
                       search_path_fulldata, search_terms,
                       search_terms_was_null, search_out_rks, parallel, ...) {
  ## Pre-processing ---------------------------------------------------------

  has_grp <- formula_contains_group_terms(refmodel$formula)

  if (inherits(refmodel, "datafit")) {
    stop("LOO can be performed only if the reference model is a genuine ",
         "probabilistic model for which the log-likelihood can be evaluated.")
  }

  # Log-likelihood values for the reference model (necessary for the PSIS-LOO-CV
  # weights, but also for performance statistics like ELPD, MLPD, and GMPD):
  if (refmodel$family$for_latent) {
    mu_offs_oscale <- refmodel$family$latent_ilink(
      t(refmodel$mu_offs), cl_ref = seq_along(refmodel$wdraws_ref),
      wdraws_ref = refmodel$wdraws_ref
    )
    if (length(dim(mu_offs_oscale)) < 2) {
      stop("Unexpected structure for the output of `latent_ilink`.")
    }
    loglik_forPSIS <- refmodel$family$latent_ll_oscale(
      mu_offs_oscale, dis = refmodel$dis, y_oscale = refmodel$y_oscale,
      wobs = refmodel$wobs, cl_ref = seq_along(refmodel$wdraws_ref),
      wdraws_ref = refmodel$wdraws_ref
    )
    if (!is.matrix(loglik_forPSIS)) {
      stop("Unexpected structure for the output of `latent_ll_oscale`.")
    }
    if (all(is.na(loglik_forPSIS))) {
      stop("In case of the latent projection, `cv_method = \"LOO\"` requires ",
           "a function `latent_ll_oscale` that does not return only `NA`s.")
    }
    if (length(dim(mu_offs_oscale)) == 3) {
      # In this case, `mu_offs_oscale` is a 3-dimensional array (S x N x C), so
      # coerce it to an augmented-rows matrix:
      mu_offs_oscale <- arr2augmat(mu_offs_oscale, margin_draws = 1)
      # In the corresponding `else` case, `mu_offs_oscale` is a matrix (S x N).
      # Transposing it to an N x S matrix would be more consistent with
      # projpred's internal convention, but avoiding the transposition is
      # computationally more efficient.
    }
  } else {
    loglik_forPSIS <- t(refmodel$family$ll_fun(
      refmodel$mu_offs, refmodel$dis, refmodel$y, refmodel$wobs
    ))
  }
  n <- ncol(loglik_forPSIS)

  # PSIS-LOO-CV weights:
  if (length(unique(refmodel$wdraws_ref)) != 1) {
    stop("Currently, projpred requires the reference model's posterior draws ",
         "to have constant weights.")
  }
  if (nrow(loglik_forPSIS) <= 1) {
    stop("Currently, more than one posterior draw from the reference model is ",
         "needed (because projpred relies on loo::psis() for PSIS-LOO-CV).")
  }
  # Call loo::psis() and while doing so, catch messages and warnings via
  # capt_mssgs_warns() to filter out some of them.
  mssgs_warns_capt <- capt_mssgs_warns(
    psisloo <- loo::psis(-loglik_forPSIS, cores = 1, r_eff = NA)
  )
  # Filter out the Pareto k-value warning (we throw a customized one instead):
  mssgs_warns_capt <- grep(
    "Some Pareto k diagnostic values are (too|slightly) high", mssgs_warns_capt,
    value = TRUE, invert = TRUE
  )
  if (length(mssgs_warns_capt) > 0) {
    warning(mssgs_warns_capt)
  }
  pareto_k <- loo::pareto_k_values(psisloo)
  # Within projpred, moment matching and mixture importance sampling (as well
  # as reference model refits leaving out each problematic observation in
  # turn, i.e., brms's `reloo` argument) currently cannot be used because all
  # these techniques result in new MCMC draws for the reference model, meaning
  # that the projection would have to be adapted. Therefore, it is easier to
  # recommend K-fold-CV (for the reference model refits, i.e., brms's `reloo`
  # argument, another reason is that they can quickly become as costly as
  # K-fold-CV).
  warn_pareto(
    n07 = sum(pareto_k > .ps_khat_threshold(dim(psisloo)[1])), n = n,
    khat_threshold = .ps_khat_threshold(dim(psisloo)[1]),
    warn_txt = paste0("Some (%d / %d) Pareto k's for the reference model's ",
                      "PSIS-LOO weights are > %s.")
  )
  lw <- weights(psisloo)

  if (refmodel$family$for_latent) {
    # Need to re-calculate the latent response values in `refmodel$y` by
    # incorporating the PSIS weights because `refmodel$y` resulted from applying
    # `colMeans(posterior_linpred())` to the original (full-data) reference
    # model fit, so using `refmodel$y` would induce a dependency between
    # training and test data:
    y_lat_E <- loo::E_loo(
      t(refmodel$ref_predfun(
        refmodel$fit, excl_offs = FALSE,
        mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)
      )),
      psis_object = psisloo,
      log_ratios = -loglik_forPSIS
    )
    # The k-values are h-specific (expectation-specific) here (see Vehtari et
    # al., 2024, <https://jmlr.org/papers/v25/19-556.html>, beginning of
    # section 3, section 3.2.8, appendix D, and appendix E).
    warn_pareto(
      n07 = sum(y_lat_E$pareto_k > .ps_khat_threshold(dim(psisloo)[1])), n = n,
      khat_threshold = .ps_khat_threshold(dim(psisloo)[1]),
      warn_txt = paste0(
        "In the recalculation of the latent response values, some (%d / %d) ",
        "expectation-specific Pareto k-values are > %s.\n",
        "In general, we recommend K-fold-CV in this case."
      )
    )
    refmodel$y <- y_lat_E$value
  }

  loo_ref_oscale <- apply(loglik_forPSIS + lw, 2, log_sum_exp)

  if (validate_search && nloo < n) {
    # Select which LOO-folds get more accurate computation using simple
    # random sampling without resampling (Magnusson et al., 2020)
    inds <- sample.int(n, size = nloo, replace = FALSE)
  } else {
    inds <- seq_len(n)
  }

  # Initialize objects where to store the results:
  loo_sub <- replicate(nterms_max + 1L, rep(NA, n), simplify = FALSE)
  mu_sub <- replicate(
    nterms_max + 1L,
    structure(rep(NA, nrow(refmodel$mu_offs)),
              ndiscrete = attr(refmodel$mu_offs, "ndiscrete"),
              class = sub("augmat", "augvec", oldClass(refmodel$mu_offs),
                          fixed = TRUE)),
    simplify = FALSE
  )
  if (refmodel$family$for_latent) {
    loo_sub_oscale <- loo_sub
    # In general, we could use `mu_sub_oscale <- mu_sub` here, but the case
    # where refmodel$family$latent_ilink() returns a 3-dimensional array (S x N
    # x C) needs special care.
    if (!is.null(refmodel$family$cats)) {
      mu_sub_oscale <- replicate(
        nterms_max + 1L,
        structure(rep(NA, n * length(refmodel$family$cats)),
                  ndiscrete = length(refmodel$family$cats),
                  class = "augvec"),
        simplify = FALSE
      )
    } else {
      mu_sub_oscale <- mu_sub
    }
  }

  if (!validate_search) {
    ## Case `validate_search = FALSE` -----------------------------------------

    # NOTE: The case where `inds` is an actual subset of the set of all
    # observation indices should never occur here in the
    # `validate_search = FALSE` case. Thus, in principle, the code could be
    # simplified here, but keeping `inds` in case this might be helpful in the
    # future.
    if (nloo < n) {
      stop("`nloo < n` is unexpected if `validate_search = FALSE`")
    }

    # "Run" the performance evaluation for the submodels along the predictor
    # ranking (in fact, we only prepare the performance evaluation by computing
    # precursor quantities, but for users, this difference is not perceivable):
    # * Step 1: Re-project (using the full dataset) onto the submodels along the
    #   full-data predictor ranking and evaluate their predictive performance.
    perf_eval_out <- perf_eval(
      search_path = search_path_fulldata, refmodel = refmodel,
      refit_prj = refit_prj, ndraws = ndraws_pred, nclusters = nclusters_pred,
      return_preds = TRUE, return_p_ref = TRUE, return_dis_sub = TRUE,
      indices_test = inds, verbose = verbose, ...
    )
    clust_used_eval <- perf_eval_out[["clust_used"]]
    nprjdraws_eval <- perf_eval_out[["nprjdraws"]]
    refdist_eval <- perf_eval_out[["p_ref"]]
    # * Step 2: Weight the full-data performance evaluation results according to
    #   the PSIS-LOO-CV weights.
    verb_out("-----\nWeighting the full-data performance evaluation results ",
             "using the PSIS-LOO-CV weights ...", verbose = verbose)
    if (refmodel$family$for_latent) {
      refdist_eval_mu_offs_oscale <- refmodel$family$latent_ilink(
        t(refdist_eval$mu_offs), cl_ref = refdist_eval$cl,
        wdraws_ref = refdist_eval$wdraws_orig
      )
      if (length(dim(refdist_eval_mu_offs_oscale)) == 3) {
        refdist_eval_mu_offs_oscale <- refdist_eval_mu_offs_oscale[, inds, ,
                                                                   drop = FALSE]
      } else {
        refdist_eval_mu_offs_oscale <- refdist_eval_mu_offs_oscale[, inds,
                                                                   drop = FALSE]
      }
      log_lik_ref <- refmodel$family$latent_ll_oscale(
        refdist_eval_mu_offs_oscale, dis = refdist_eval$dis,
        y_oscale = refmodel$y_oscale[inds], wobs = refmodel$wobs[inds],
        cl_ref = refdist_eval$cl, wdraws_ref = refdist_eval$wdraws_orig
      )
      if (all(is.na(log_lik_ref))) {
        stop("In case of the latent projection, `validate_search = FALSE` ",
             "requires a function `latent_ll_oscale` that does not return ",
             "only `NA`s.")
      }
    } else {
      inds_aug <- inds
      if (refmodel$family$for_augdat) {
        inds_aug <- inds_aug + rep(
          (seq_along(refmodel$family$cats) - 1L) * n,
          each = length(inds_aug)
        )
      }
      log_lik_ref <- t(refmodel$family$ll_fun(
        refdist_eval$mu_offs[inds_aug, , drop = FALSE], refdist_eval$dis,
        refmodel$y[inds], refmodel$wobs[inds]
      ))
    }
    if (nrow(log_lik_ref) > 1) {
      # Take into account that clustered draws usually have different weights:
      lwref_lwdraws <- -log_lik_ref + log(refdist_eval$wdraws_prj)
      # This re-weighting requires a re-normalization (as.array() is applied to
      # have stricter consistency checks, see `?sweep`):
      lwref_lwdraws <- sweep(lwref_lwdraws, 2, as.array(apply(lwref_lwdraws, 2,
                                                              log_sum_exp)))
      # Internally, loo::psis() doesn't perform the Pareto smoothing if the
      # number of draws is small (as indicated by object `no_psis_eval`, see
      # below). In projpred, this can occur, e.g., if users request a number
      # of projected draws (for performance evaluation, either after
      # clustering or thinning the reference model's posterior draws) that is
      # much smaller than the default of 400. In order to throw a customized
      # warning message (and to avoid the calculation of Pareto k-values, see
      # loo issue stan-dev/loo#227), object `no_psis_eval` indicates whether
      # loo::psis() would perform the Pareto smoothing or not (for the
      # decision rule, see loo:::n_pareto() and loo:::enough_tail_samples(),
      # keeping in mind that we have `r_eff = 1` for all observations here).
      S_for_psis_eval <- nrow(log_lik_ref)
      no_psis_eval <- ceiling(min(0.2 * S_for_psis_eval,
                                  3 * sqrt(S_for_psis_eval))) < 5
      if (no_psis_eval) {
        if (getOption("projpred.warn_psis", TRUE)) {
          message(
            "Using standard importance sampling (SIS) due to a small number ",
            "of ", if (clust_used_eval) "clusters" else "draws", "."
          )
        }
        # Use loo::sis().
        # In principle, we could rely on loo::psis() here (because in such a
        # case, it would internally switch to SIS automatically), but using
        # loo::sis() explicitly is safer because if the loo package changes
        # its decision rule, we would get a mismatch between our customized
        # warning here and the IS method used by loo. See also loo issue
        # stan-dev/loo#227.
        importance_sampling_nm <- "sis"
      } else {
        # Use loo::psis().
        # Usually, we have a small number of projected draws here (400 by
        # default), which means that the 'loo' package will automatically
        # perform the regularization from Vehtari et al. (2024,
        # <https://jmlr.org/papers/v25/19-556.html>, appendix G).
        importance_sampling_nm <- "psis"
      }
      importance_sampling_func <- get(importance_sampling_nm,
                                      asNamespace("loo"))
      mssgs_warns_capt <- capt_mssgs_warns(
        sub_psisloo <- importance_sampling_func(lwref_lwdraws, cores = 1,
                                                r_eff = NA)
      )
      # Filter out Pareto k-value warnings (we throw a customized one instead):
      mssgs_warns_capt <- grep(
        "Some Pareto k diagnostic values are (too|slightly) high",
        mssgs_warns_capt, value = TRUE, invert = TRUE
      )
      if (length(mssgs_warns_capt) > 0) {
        warning(mssgs_warns_capt)
      }

      if (importance_sampling_nm == "psis") {
        pareto_k_eval <- loo::pareto_k_values(sub_psisloo)
        warn_pareto(
          n07 = sum(pareto_k_eval > .ps_khat_threshold(dim(psisloo)[1])), n = n,
          khat_threshold = .ps_khat_threshold(dim(sub_psisloo)[1]),
          warn_txt = paste0(
            "Some (%d / %d) Pareto k's for the reference model's PSIS-LOO ",
            "weights given ", txt_clust_draws(clust_used_eval, nprjdraws_eval),
            " are > %s."
          )
        )
      }
      lw_sub <- weights(sub_psisloo)
    } else {
      lw_sub <- matrix(0, nrow = nrow(log_lik_ref), ncol = ncol(log_lik_ref))
    }
    for (k in seq_len(1 + length(search_path_fulldata$predictor_ranking))) {
      # TODO: For consistency, replace `k` in this `for` loop by `j`.
      mu_k <- perf_eval_out[["mu_by_size"]][[k]]
      log_lik_sub <- perf_eval_out[["lppd_by_size"]][[k]]
      loo_sub[[k]][inds] <- apply(log_lik_sub + lw_sub, 2, log_sum_exp)
      if (refmodel$family$for_latent) {
        mu_k_oscale <- refmodel$family$latent_ilink(
          t(mu_k), cl_ref = refdist_eval$cl,
          wdraws_ref = refdist_eval$wdraws_orig
        )
        log_lik_sub_oscale <- refmodel$family$latent_ll_oscale(
          mu_k_oscale, dis = perf_eval_out[["dis_sub"]][[k]],
          y_oscale = refmodel$y_oscale[inds], wobs = refmodel$wobs[inds],
          cl_ref = refdist_eval$cl, wdraws_ref = refdist_eval$wdraws_orig
        )
        loo_sub_oscale[[k]][inds] <- apply(log_lik_sub_oscale + lw_sub, 2,
                                           log_sum_exp)
        if (length(dim(mu_k_oscale)) == 3) {
          mu_k_oscale <- arr2augmat(mu_k_oscale, margin_draws = 1)
        }
      }
      for (run_index in seq_along(inds)) {
        i_aug <- inds[run_index]
        run_index_aug <- run_index
        if (!is.null(refmodel$family$cats)) {
          i_aug <- i_aug + (seq_along(refmodel$family$cats) - 1L) * n
          run_index_aug <- run_index_aug +
            (seq_along(refmodel$family$cats) - 1L) * nloo
        }
        i_flx <- i_aug
        run_index_flx <- run_index_aug
        if (refmodel$family$for_latent && !is.null(refmodel$family$cats)) {
          i_flx <- inds[run_index]
          run_index_flx <- run_index
        }
        mu_sub[[k]][i_flx] <- mu_k[run_index_flx, , drop = FALSE] %*%
          exp(lw_sub[, run_index])
        if (refmodel$family$for_latent) {
          if (inherits(mu_k_oscale, "augmat")) {
            mu_sub_oscale[[k]][i_aug] <- mu_k_oscale[run_index_aug, ,
                                                     drop = FALSE] %*%
              exp(lw_sub[, run_index])
          } else {
            # In principle, we could use the same code for averaging across the
            # draws as above in the `augmat` case. However, that would require
            # `mu_k_oscale <- t(mu_k_oscale)` beforehand, so the following
            # should be more efficient:
            mu_sub_oscale[[k]][i_aug] <- exp(lw_sub[, run_index]) %*%
              mu_k_oscale[, run_index_aug, drop = FALSE]
          }
        }
      }
    }
    verb_out("-----", verbose = verbose)
    # Needed for cutting off post-processed results later:
    prv_len_rk <- length(search_path_fulldata$predictor_ranking)
  } else {
    ## Case `validate_search = TRUE` ------------------------------------------

    search_out_rks_was_null <- is.null(search_out_rks)
    if (search_out_rks_was_null) {
      refdist_sel <- get_refdist(refmodel, ndraws = ndraws,
                                 nclusters = nclusters)
    } else {
      refdist_sel <- NULL
    }
    cl_sel <- refdist_sel$cl
    if (refit_prj) {
      refdist_pred <- get_refdist(refmodel, ndraws = ndraws_pred,
                                  nclusters = nclusters_pred)
    } else {
      refdist_pred <- NULL
    }
    cl_pred <- refdist_pred$cl

    if (verbose) {
      if (refit_prj) {
        vtxt_clust_used_eval <- refdist_pred[["clust_used"]]
        vtxt_nprjdraws_eval <- refdist_pred[["nprjdraws"]]
      } else {
        # NOTE: `!refit_prj` cannot occur in combination with
        # `!search_out_rks_was_null`, so it is correct and safe to use
        # `refdist_sel` here.
        vtxt_clust_used_eval <- refdist_sel[["clust_used"]]
        vtxt_nprjdraws_eval <- refdist_sel[["nprjdraws"]]
      }
      verb_out("-----\nRunning ",
               if (!search_out_rks_was_null) {
                 ""
               } else {
                 paste0(method, " search with ",
                        txt_clust_draws(refdist_sel[["clust_used"]],
                                        refdist_sel[["nprjdraws"]]),
                        " and ")
               },
               "the performance evaluation with ",
               txt_clust_draws(vtxt_clust_used_eval, vtxt_nprjdraws_eval),
               " (`refit_prj = ", refit_prj, "`) for each of the `nloo = ",
               nloo, "` LOO-CV folds separately ...")
    }
    one_obs <- function(run_index,
                        verbose_obs = max(verbose - 1L, 0L),
                        ...) {
      # Observation index:
      i <- inds[run_index]

      # For (extra-)verbose mode:
      vtxt_obs_i <- paste0(
        "for LOO-CV fold ", run_index, " out of `nloo = ", nloo, "` ",
        "(", round(100 * run_index / nloo), " %) ",
        if (nloo < n) {
          paste0("(this is observation ", i, ") ")
        } else {
          ""
        }
      )

      # Run the search with the reweighted clusters (or thinned draws) (so the
      # *reweighted* fitted response values from the reference model act as
      # artifical response values in the projection (or L1-penalized
      # projection)):
      if (!search_out_rks_was_null) {
        search_path <- list(predictor_ranking = search_out_rks[[run_index]])
      } else {
        search_path <- .select(
          refmodel = refmodel, ndraws = ndraws, nclusters = nclusters,
          reweighting_args = list(cl_ref = cl_sel, wdraws_ref = exp(lw[, i])),
          method = method, nterms_max = nterms_max, penalty = penalty,
          verbose = verbose_obs, verbose_line_length = 3,
          verbose_txt_add = vtxt_obs_i, search_control = search_control,
          search_terms = search_terms, est_runtime = FALSE, ...
        )
      }

      # Run the performance evaluation for the submodels along the predictor
      # ranking:
      perf_eval_out <- perf_eval(
        search_path = search_path, refmodel = refmodel, refit_prj = refit_prj,
        ndraws = ndraws_pred, nclusters = nclusters_pred,
        reweighting_args = list(cl_ref = cl_pred, wdraws_ref = exp(lw[, i])),
        indices_test = i, verbose = verbose_obs, verbose_line_length = 3,
        verbose_txt_add = vtxt_obs_i, ...
      )

      return(nlist(predictor_ranking = search_path[["predictor_ranking"]],
                   summaries_sub = perf_eval_out[["sub_summaries"]],
                   clust_used_eval = perf_eval_out[["clust_used"]],
                   nprjdraws_eval = perf_eval_out[["nprjdraws"]]))
    }
    if (!parallel) {
      # Sequential case. Actually, we could simply use ``%do_projpred%` <-
      # foreach::`%do%`` here and then proceed as in the parallel case, but that
      # would require adding more "hard" dependencies (because packages
      # 'foreach' and 'doRNG' would have to be moved from `Suggests:` to
      # `Imports:`).
      if (verbose == 1L) {
        pb <- utils::txtProgressBar(min = 0, max = nloo, style = 3, initial = 0,
                                    file = stderr())
      }
      res_cv <- lapply(seq_along(inds), function(run_index) {
        if (verbose == 1L) {
          on.exit(utils::setTxtProgressBar(pb, run_index))
        }
        one_obs(run_index, ...)
      })
      if (verbose == 1L) {
        close(pb)
      }
    } else {
      # Parallel case.
      if (!requireNamespace("foreach", quietly = TRUE)) {
        stop("Please install the 'foreach' package.")
      }
      if (!requireNamespace("doRNG", quietly = TRUE)) {
        stop("Please install the 'doRNG' package.")
      }
      if (verbose && use_progressr()) {
        progressor_obj <- progressr::progressor(length(inds))
      } else {
        progressor_obj <- NULL
      }
      dot_args <- list(...)
      `%do_projpred%` <- doRNG::`%dorng%`
      res_cv <- foreach::foreach(
        run_index = seq_along(inds),
        .packages = c("projpred"),
        .export = c("one_obs", "dot_args", "progressor_obj",
                    getOption("projpred.export_to_workers", character())),
        .noexport = c("mu_offs_oscale", "loglik_forPSIS", "psisloo", "y_lat_E",
                      "loo_ref_oscale", "validset", "loo_sub", "mu_sub",
                      "loo_sub_oscale", "mu_sub_oscale"),
        .errorhandling = getOption("projpred.foreach_errorhandling", "stop"),
        .verbose = getOption("projpred.foreach_verbose", FALSE)
      ) %do_projpred% {
        out_one_obs <- do.call(one_obs, c(list(run_index = run_index,
                                               verbose_obs = FALSE),
                                          dot_args))
        if (!is.null(progressor_obj)) progressor_obj()
        return(out_one_obs)
      }
    }
    # For storing the fold-wise predictor rankings:
    predictor_ranking_mat <- matrix(nrow = n, ncol = nterms_max)
    # Needed for checking that the length of the predictor ranking is the same
    # across all CV folds and for cutting off post-processed results later:
    prv_len_rk <- NULL
    # For checking that `clust_used_eval` is the same across all CV folds (and
    # also for storing it):
    clust_used_eval <- NULL
    # For checking that `nprjdraws_eval` is the same across all CV folds (and
    # also for storing it):
    nprjdraws_eval <- NULL
    for (run_index in seq_along(inds)) {
      i <- inds[run_index]

      summaries_sub <- res_cv[[run_index]][["summaries_sub"]]
      i_aug <- i
      if (!is.null(refmodel$family$cats)) {
        i_aug <- i_aug + (seq_along(refmodel$family$cats) - 1L) * n
      }
      i_flx <- i_aug
      if (refmodel$family$for_latent && !is.null(refmodel$family$cats)) {
        i_flx <- i
      }
      for (k in seq_along(summaries_sub)) {
        loo_sub[[k]][i] <- summaries_sub[[k]]$lppd
        mu_sub[[k]][i_flx] <- summaries_sub[[k]]$mu
        if (!is.null(summaries_sub[[k]]$oscale)) {
          loo_sub_oscale[[k]][i] <- summaries_sub[[k]]$oscale$lppd
          mu_sub_oscale[[k]][i_aug] <- summaries_sub[[k]]$oscale$mu
        }
      }

      rk_i <- res_cv[[run_index]][["predictor_ranking"]]
      if (is.null(prv_len_rk)) {
        prv_len_rk <- length(rk_i)
      } else if (getOption("projpred.additional_checks", FALSE)) {
        stopifnot(identical(length(rk_i), prv_len_rk))
      }
      predictor_ranking_mat[i, seq_along(rk_i)] <- rk_i

      if (is.null(clust_used_eval)) {
        clust_used_eval <- res_cv[[run_index]][["clust_used_eval"]]
      } else if (getOption("projpred.additional_checks", FALSE)) {
        stopifnot(identical(res_cv[[run_index]][["clust_used_eval"]],
                            clust_used_eval))
      }
      if (is.null(nprjdraws_eval)) {
        nprjdraws_eval <- res_cv[[run_index]][["nprjdraws_eval"]]
      } else if (getOption("projpred.additional_checks", FALSE)) {
        stopifnot(identical(res_cv[[run_index]][["nprjdraws_eval"]],
                            nprjdraws_eval))
      }
    }
    verb_out("-----", verbose = verbose)
  }

  ## Post-processing --------------------------------------------------------

  # Submodel predictive performance:
  summ_sub <- lapply(seq_len(prv_len_rk + 1L), function(k) {
    summ_k <- list(lppd = loo_sub[[k]], mu = mu_sub[[k]])
    if (refmodel$family$for_latent) {
      summ_k$oscale <- list(lppd = loo_sub_oscale[[k]], mu = mu_sub_oscale[[k]])
    }
    return(summ_k)
  })

  # Reference model predictive performance:
  if (has_grp && getOption("projpred.mlvl_pred_new", FALSE)) {
    # Need to use `mlvl_allrandom = TRUE` (`refmodel$mu_offs` is based on
    # `mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)`):
    eta_offs_mlvlRan <- refmodel$ref_predfun(refmodel$fit, excl_offs = FALSE)
    mu_offs_mlvlRan <- refmodel$family$linkinv(eta_offs_mlvlRan)
  } else {
    mu_offs_mlvlRan <- refmodel$mu_offs
  }
  mu_ref <- as.vector(do.call(rbind, lapply(seq_len(n), function(i) {
    # For the augmented-data projection, `mu_offs_mlvlRan` is an augmented-rows
    # matrix whereas the columns of `lw` refer to the original (non-augmented)
    # observations. Since `i` refers to the columns of `lw` (we have
    # `n == ncol(lw)`), the indices for `mu_offs_mlvlRan` need to be adapted:
    i_aug <- i
    if (!is.null(refmodel$family$cats)) {
      i_aug <- i_aug + (seq_along(refmodel$family$cats) - 1L) * n
    }
    i_flx <- i_aug
    if (refmodel$family$for_latent && !is.null(refmodel$family$cats)) {
      i_flx <- i
    }
    return(as.vector(mu_offs_mlvlRan[i_flx, , drop = FALSE] %*% exp(lw[, i])))
  })))
  mu_ref <- structure(
    mu_ref,
    ndiscrete = attr(mu_offs_mlvlRan, "ndiscrete"),
    class = sub("augmat", "augvec", oldClass(mu_offs_mlvlRan), fixed = TRUE)
  )
  if (refmodel$family$for_latent) {
    loglik_lat <- t(refmodel$family$ll_fun(
      mu_offs_mlvlRan, refmodel$dis, refmodel$y, refmodel$wobs
    ))
    lppd_ref <- apply(loglik_lat + lw, 2, log_sum_exp)
  } else {
    if (has_grp && getOption("projpred.mlvl_pred_new", FALSE)) {
      # Need to use `mlvl_allrandom = TRUE` (`loo_ref_oscale` is based on
      # `mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)`):
      loglik_mlvlRan <- t(refmodel$family$ll_fun(
        mu_offs_mlvlRan, refmodel$dis, refmodel$y, refmodel$wobs
      ))
      lppd_ref <- apply(loglik_mlvlRan + lw, 2, log_sum_exp)
    } else {
      lppd_ref <- loo_ref_oscale
    }
  }
  summ_ref <- list(lppd = lppd_ref, mu = mu_ref)
  if (refmodel$family$for_latent) {
    if (has_grp && getOption("projpred.mlvl_pred_new", FALSE)) {
      # Need to use `mlvl_allrandom = TRUE` (`mu_offs_oscale` is based on
      # `mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)`):
      mu_offs_mlvlRan_oscale <- refmodel$family$latent_ilink(
        t(mu_offs_mlvlRan), cl_ref = seq_along(refmodel$wdraws_ref),
        wdraws_ref = refmodel$wdraws_ref
      )
      mu_offs_mlvlRan_oscale_odim <- mu_offs_mlvlRan_oscale
      if (length(dim(mu_offs_mlvlRan_oscale)) == 3) {
        mu_offs_mlvlRan_oscale <- arr2augmat(mu_offs_mlvlRan_oscale,
                                             margin_draws = 1)
      }
    } else {
      mu_offs_mlvlRan_oscale <- mu_offs_oscale
    }
    mu_ref_oscale <- as.vector(do.call(rbind, lapply(seq_len(n), function(i) {
      i_aug <- i
      if (!is.null(refmodel$family$cats)) {
        i_aug <- i_aug + (seq_along(refmodel$family$cats) - 1L) * n
      }
      if (inherits(mu_offs_mlvlRan_oscale, "augmat")) {
        return(as.vector(mu_offs_mlvlRan_oscale[i_aug, , drop = FALSE] %*%
                           exp(lw[, i])))
      } else {
        # In principle, we could use the same code for averaging across the
        # draws as above in the `augmat` case. However, that would require
        # `mu_offs_mlvlRan_oscale <- t(mu_offs_mlvlRan_oscale)` beforehand, so
        # the following should be more efficient:
        return(exp(lw[, i]) %*% mu_offs_mlvlRan_oscale[, i_aug, drop = FALSE])
      }
    })))
    mu_ref_oscale <- structure(
      mu_ref_oscale,
      ndiscrete = attr(mu_offs_mlvlRan_oscale, "ndiscrete"),
      class = sub("augmat", "augvec", oldClass(mu_offs_mlvlRan_oscale),
                  fixed = TRUE)
    )
    if (has_grp && getOption("projpred.mlvl_pred_new", FALSE)) {
      # Need to use `mlvl_allrandom = TRUE` (`loo_ref_oscale` is based on
      # `mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)`):
      loglik_mlvlRan <- refmodel$family$latent_ll_oscale(
        mu_offs_mlvlRan_oscale_odim, dis = refmodel$dis,
        y_oscale = refmodel$y_oscale, wobs = refmodel$wobs,
        cl_ref = seq_along(refmodel$wdraws_ref),
        wdraws_ref = refmodel$wdraws_ref
      )
      lppd_ref_oscale <- apply(loglik_mlvlRan + lw, 2, log_sum_exp)
    } else {
      lppd_ref_oscale <- loo_ref_oscale
    }
    summ_ref$oscale <- list(lppd = lppd_ref_oscale, mu = mu_ref_oscale)
  }

  # Combined submodel and reference model predictive performance:
  summaries <- list(sub = summ_sub, ref = summ_ref)

  if (!validate_search) {
    out_list <- nlist(ce = perf_eval_out[["ce"]])
  } else {
    out_list <- nlist(predictor_ranking_cv = predictor_ranking_mat[
      , seq_len(prv_len_rk), drop = FALSE
    ])
  }
  out_list <- c(out_list,
                nlist(summaries,
                      y_wobs_test = as.data.frame(refmodel[nms_y_wobs_test()]),
                      clust_used_eval, nprjdraws_eval, inds))
  return(out_list)
}

warn_pareto <- function(n07, n, khat_threshold = 0.7, warn_txt) {
  if (!getOption("projpred.warn_psis", TRUE) || (n07 == 0)) return()
  warning(sprintf(warn_txt, n07, n, as.character(round(khat_threshold, 2))),
          call. = FALSE)
  return()
}

#' Pareto-smoothing k-hat threshold
#'
#' Copied from loo package. Remove after loo package exposes this.
#'
#' Given sample size S computes khat threshold for reliable Pareto
#' smoothed estimate (to have small probability of large error). See
#' section 3.2.4, equation (13). Sample sizes 100, 320, 1000, 2200,
#' 10000 correspond to thresholds 0.5, 0.6, 0.67, 0.7, 0.75. Although
#' with bigger sample size S we can achieve estimates with small
#' probability of large error, it is difficult to get accurate MCSE
#' estimates as the bias starts to dominate when k > 0.7 (see Section 3.2.3).
#' Thus the sample size dependend k-ht threshold is capped at 0.7.
#' @param S sample size
#' @param ... unused
#' @return threshold
#' @noRd
.ps_khat_threshold <- function(S, ...) {
  min(1 - 1 / log10(S), 0.7)
}

# K-fold-CV ---------------------------------------------------------------

# Needed to avoid a NOTE in `R CMD check`:
if (getRversion() >= package_version("2.15.1")) {
  utils::globalVariables("list_cv_k")
  utils::globalVariables("search_out_rks_k")
  utils::globalVariables("ks_k")
}

kfold_varsel <- function(refmodel, method, nterms_max, ndraws, nclusters,
                         ndraws_pred, nclusters_pred, refit_prj, penalty,
                         verbose, search_control, K, cvfits, validate_search,
                         search_path_fulldata, search_terms, search_out_rks,
                         parallel, ...) {
  # Fetch the K reference model fits (or fit them now if not already done) and
  # create objects of class `refmodel` from them (and also store the `omitted`
  # indices):
  list_cv <- get_kfold(refmodel, K = K, cvfits = cvfits, verbose = verbose)
  K <- length(list_cv)

  search_out_rks_was_null <- is.null(search_out_rks)
  if (search_out_rks_was_null) {
    search_out_rks <- replicate(K, NULL, simplify = FALSE)
  }

  if (refmodel$family$for_latent) {
    # Need to set the latent response values in `refmodel$y` to `NA`s because
    # `refmodel$y` resulted from applying `colMeans(posterior_linpred())` to the
    # original (full-data) reference model fit, so using the `fold$omitted`
    # subset of `refmodel$y` as (latent) response values in fold k of K would
    # induce a dependency between training and test data:
    refmodel$y <- rep(NA, refmodel$nobs)
  }
  y_wobs_test <- as.data.frame(refmodel[nms_y_wobs_test()])

  if (verbose) {
    # Here in kfold_varsel(), we have no get_refdist() (or get_p_clust()) output
    # available, so use clust_info() as a workaround:
    clust_info_sel <- clust_info(
      ndraws = ndraws,
      nclusters = nclusters,
      S = length(refmodel$wdraws_ref)
    )
    vtxt_clust_used_sel <- clust_info_sel[["clust_used"]]
    vtxt_nprjdraws_sel <- clust_info_sel[["nprjdraws"]]
    if (refit_prj) {
      clust_info_eval <- clust_info(
        ndraws = ndraws_pred,
        nclusters = nclusters_pred,
        S = length(refmodel$wdraws_ref)
      )
      vtxt_clust_used_eval <- clust_info_eval[["clust_used"]]
      vtxt_nprjdraws_eval <- clust_info_eval[["nprjdraws"]]
    } else {
      # NOTE: `!refit_prj` cannot occur in combination with
      # `!search_out_rks_was_null || !validate_search`, so it is correct and
      # safe to use `vtxt_clust_used_sel` and `vtxt_nprjdraws_sel` here.
      vtxt_clust_used_eval <- vtxt_clust_used_sel
      vtxt_nprjdraws_eval <- vtxt_nprjdraws_sel
    }
    verb_out("-----\nRunning ",
             if (!search_out_rks_was_null || !validate_search) {
               ""
             } else {
               paste0(method, " search with ",
                      txt_clust_draws(vtxt_clust_used_sel, vtxt_nprjdraws_sel),
                      " and ")
             },
             "the performance evaluation with ",
             txt_clust_draws(vtxt_clust_used_eval, vtxt_nprjdraws_eval),
             " (`refit_prj = ", refit_prj, "`) for each of the K = ", K,
             " CV folds separately ...")
  }
  one_fold <- function(fold,
                       rk,
                       k,
                       verbose_fold = max(verbose - 1L, 0L),
                       ...) {
    # For (extra-)verbose mode:
    vtxt_fold_k <- paste0("for CV fold ", k, " out of K = ", K, " (",
                          round(100 * k / K), " %) ")

    # Run the search for the current fold:
    if (!validate_search) {
      search_path <- search_path_fulldata
    } else if (!search_out_rks_was_null) {
      search_path <- list(predictor_ranking = rk)
    } else {
      search_path <- .select(
        refmodel = fold$refmodel, ndraws = ndraws, nclusters = nclusters,
        method = method, nterms_max = nterms_max, penalty = penalty,
        verbose = verbose_fold, verbose_line_length = 3,
        verbose_txt_add = vtxt_fold_k, search_control = search_control,
        search_terms = search_terms, est_runtime = FALSE, ...
      )
    }

    # Run the performance evaluation for the submodels along the predictor
    # ranking:
    perf_eval_out <- perf_eval(
      search_path = search_path, refmodel = fold$refmodel,
      refit_prj = refit_prj, ndraws = ndraws_pred, nclusters = nclusters_pred,
      refmodel_fulldata = refmodel, indices_test = fold$omitted,
      verbose = verbose_fold, verbose_line_length = 3,
      verbose_txt_add = vtxt_fold_k, ...
    )

    # Performance evaluation for the reference model of the current fold:
    eta_test <- fold$refmodel$ref_predfun(
      fold$refmodel$fit,
      newdata = refmodel$fetch_data(obs = fold$omitted),
      excl_offs = FALSE
    )
    mu_test <- fold$refmodel$family$linkinv(eta_test)
    summaries_ref <- weighted_summary_means(
      y_wobs_test = y_wobs_test[fold$omitted, , drop = FALSE],
      family = fold$refmodel$family,
      wdraws = fold$refmodel$wdraws_ref,
      mu = mu_test,
      dis = fold$refmodel$dis,
      cl_ref = seq_along(fold$refmodel$wdraws_ref)
    )

    return(nlist(predictor_ranking = search_path[["predictor_ranking"]],
                 summaries_sub = perf_eval_out[["sub_summaries"]],
                 summaries_ref, clust_used_eval = perf_eval_out[["clust_used"]],
                 nprjdraws_eval = perf_eval_out[["nprjdraws"]]))
  }
  if (!parallel) {
    # Sequential case. Actually, we could simply use ``%do_projpred%` <-
    # foreach::`%do%`` here and then proceed as in the parallel case, but that
    # would require adding more "hard" dependencies (because packages 'foreach'
    # and 'doRNG' would have to be moved from `Suggests:` to `Imports:`).
    if (verbose == 1L) {
      pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial = 0,
                                  file = stderr())
    }
    res_cv <- lapply(seq_len(K), function(k) {
      if (verbose == 1L) {
        on.exit(utils::setTxtProgressBar(pb, k))
      }
      one_fold(fold = list_cv[[k]], rk = search_out_rks[[k]], k = k, ...)
    })
    if (verbose == 1L) {
      close(pb)
    }
  } else {
    # Parallel case.
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Please install the 'foreach' package.")
    }
    if (!requireNamespace("doRNG", quietly = TRUE)) {
      stop("Please install the 'doRNG' package.")
    }
    if (verbose && use_progressr()) {
      progressor_obj <- progressr::progressor(K)
    } else {
      progressor_obj <- NULL
    }
    dot_args <- list(...)
    `%do_projpred%` <- doRNG::`%dorng%`
    res_cv <- foreach::foreach(
      list_cv_k = list_cv,
      search_out_rks_k = search_out_rks,
      ks_k = seq_len(K),
      .packages = c("projpred"),
      .export = c("one_fold", "dot_args", "progressor_obj",
                  getOption("projpred.export_to_workers", character())),
      .noexport = c("list_cv", "search_out_rks"),
      .errorhandling = getOption("projpred.foreach_errorhandling", "stop"),
      .verbose = getOption("projpred.foreach_verbose", FALSE)
    ) %do_projpred% {
      out_one_fold <- do_call(one_fold, c(list(fold = list_cv_k,
                                               rk = search_out_rks_k,
                                               k = ks_k,
                                               verbose_fold = FALSE),
                                          dot_args))
      if (!is.null(progressor_obj)) progressor_obj()
      return(out_one_fold)
    }
  }
  verb_out("-----", verbose = verbose)
  predictor_ranking_cv <- do.call(rbind,
                                  lapply(res_cv, "[[", "predictor_ranking"))
  clust_used_eval <- element_unq(res_cv, nm = "clust_used_eval")
  nprjdraws_eval <- element_unq(res_cv, nm = "nprjdraws_eval")

  # Handle the submodels' performance evaluation results:
  sub_foldwise <- lapply(res_cv, "[[", "summaries_sub")
  if (getRversion() >= package_version("4.2.0")) {
    sub_foldwise <- simplify2array(sub_foldwise, higher = FALSE, except = NULL)
  } else {
    sub_foldwise <- simplify2array(sub_foldwise, higher = FALSE)
    if (is.null(dim(sub_foldwise))) {
      sub_dim <- dim(predictor_ranking_cv)
      sub_dim[2] <- sub_dim[2] + 1L # +1 is for the empty model
      dim(sub_foldwise) <- rev(sub_dim)
    }
  }
  sub <- apply(sub_foldwise, 1, rbind2list)
  idxs_sorted_by_fold <- unlist(lapply(list_cv, function(fold) {
    fold$omitted
  }))
  idxs_sorted_by_fold_aug <- idxs_sorted_by_fold
  if (!is.null(refmodel$family$cats)) {
    idxs_sorted_by_fold_aug <- idxs_sorted_by_fold_aug + rep(
      (seq_along(refmodel$family$cats) - 1L) * refmodel$nobs,
      each = length(idxs_sorted_by_fold_aug)
    )
  }
  idxs_sorted_by_fold_flx <- idxs_sorted_by_fold_aug
  if (refmodel$family$for_latent && !is.null(refmodel$family$cats)) {
    idxs_sorted_by_fold_flx <- idxs_sorted_by_fold
  }
  sub <- lapply(sub, function(summ) {
    summ$mu <- summ$mu[order(idxs_sorted_by_fold_flx)]
    summ$lppd <- summ$lppd[order(idxs_sorted_by_fold)]

    if (!is.null(summ$oscale)) {
      summ$oscale$mu <- summ$oscale$mu[order(idxs_sorted_by_fold_aug)]
      summ$oscale$lppd <- summ$oscale$lppd[order(idxs_sorted_by_fold)]
    }
    return(summ)
  })

  # Handle the reference model's performance evaluation results:
  ref <- rbind2list(lapply(res_cv, "[[", "summaries_ref"))
  ref$mu <- ref$mu[order(idxs_sorted_by_fold_flx)]
  ref$lppd <- ref$lppd[order(idxs_sorted_by_fold)]
  if (!is.null(ref$oscale)) {
    ref$oscale$mu <- ref$oscale$mu[order(idxs_sorted_by_fold_aug)]
    ref$oscale$lppd <- ref$oscale$lppd[order(idxs_sorted_by_fold)]
  }

  if (!validate_search) {
    out_list <- list()
  } else {
    out_list <- nlist(predictor_ranking_cv)
  }
  out_list <- c(out_list,
                nlist(summaries = nlist(sub, ref), y_wobs_test, clust_used_eval,
                      nprjdraws_eval))
  return(out_list)
}

# Refit the reference model K times (once for each fold; `cvfun` case) or fetch
# the K reference model fits if already computed (`cvfits` case). This function
# will return a list of length K, where each element is a list with elements
# `refmodel` (output of init_refmodel()) and `omitted` (vector of indices of
# those observations which were left out for the corresponding fold).
get_kfold <- function(refmodel, K, cvfits, verbose) {
  if (is.null(cvfits)) {
    if (!is.null(refmodel$cvfun)) {
      # In this case, cvfun() provided (and `cvfits` not), so run cvfun() now.
      if (verbose && !inherits(refmodel, "datafit")) {
        verb_out("-----\nRefitting the reference model K = ", K, " times ",
                 "(using the fold-wise training data) ...")
      }
      folds <- cv_folds(refmodel$nobs, K = K,
                        seed = sample.int(.Machine$integer.max, 1))
      if (getOption("projpred.warn_kfold_refits", TRUE)) {
        cvfits <- refmodel$cvfun(folds)
      } else {
        cvfits <- suppressWarnings(refmodel$cvfun(folds))
      }
      verb_out("-----", verbose = verbose)
    } else {
      stop("For a reference model which is not of class `datafit`, either ",
           "`cvfits` or `cvfun` needs to be provided for K-fold-CV (see ",
           "`?init_refmodel`).")
    }
  } else {
    folds <- attr(cvfits, "folds")
  }
  stopifnot(!is.null(folds))
  return(lapply(seq_len(K), function(k) {
    cvfit <- cvfits[[k]]
    # Add the omitted observation indices for this fold (and the fold index `k`
    # itself):
    omitted_idxs <- which(folds == k)
    if (is.list(cvfit)) {
      cvfit$omitted <- omitted_idxs
      cvfit$projpred_k <- k
    } else {
      attr(cvfit, "omitted") <- omitted_idxs
      attr(cvfit, "projpred_k") <- k
    }
    return(list(refmodel = refmodel$cvrefbuilder(cvfit),
                omitted = omitted_idxs))
  }))
}

#' Create `cvfits` from `cvfun`
#'
#' A helper function that can be used to create input for
#' [cv_varsel.refmodel()]'s argument `cvfits` by running first [cv_folds()] and
#' then the reference model object's `cvfun` (see [init_refmodel()]). This is
#' helpful if \eqn{K}-fold CV is run multiple times based on the same \eqn{K}
#' reference model refits.
#'
#' @param object An object of class `refmodel` (returned by [get_refmodel()] or
#'   [init_refmodel()]) or an object that can be passed to argument `object` of
#'   [get_refmodel()].
#' @param K Number of folds. Must be at least 2 and not exceed the number of
#'   observations. Ignored if `folds` is not `NULL`.
#' @param folds Either `NULL` for determining the CV folds automatically via
#'   [cv_folds()] (using argument `K`) or a numeric (in fact, integer) vector
#'   giving the fold index for each observation. In the latter case, argument
#'   `K` is ignored.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. If not
#'   `NA`, then the PRNG state is reset (to the state before calling
#'   [run_cvfun()]) upon exiting [run_cvfun()].
#' @param ... For [run_cvfun.default()]: Arguments passed to [get_refmodel()].
#'   For [run_cvfun.refmodel()]: Currently ignored.
#'
#' @return An object that can be used as input for [cv_varsel.refmodel()]'s
#'   argument `cvfits`.
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
#' # Define the reference model object explicitly (not really necessary here
#' # because the get_refmodel() call is quite fast in this example, but in
#' # general, this approach is faster than defining the reference model object
#' # multiple times implicitly):
#' ref <- get_refmodel(fit)
#'
#' # Run the reference model object's `cvfun` (with a small value for `K`, but
#' # only for the sake of speed in this example; this is not recommended in
#' # general):
#' cv_fits <- run_cvfun(ref, K = 2, seed = 184)
#'
#' # Run cv_varsel() (with L1 search and small values for `nterms_max` and
#' # `nclusters_pred`, but only for the sake of speed in this example; this is
#' # not recommended in general) and use `cv_fits` there:
#' cvvs_L1 <- cv_varsel(ref, method = "L1", cv_method = "kfold",
#'                      cvfits = cv_fits, nterms_max = 3, nclusters_pred = 10,
#'                      seed = 5555)
#' # Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
#' # and `?ranking` for possible post-processing functions.
#'
#' # The purpose of run_cvfun() is to create an object that can be used in
#' # multiple cv_varsel() calls, e.g., to check the sensitivity to the search
#' # method (L1 or forward):
#' cvvs_fw <- cv_varsel(ref, method = "forward", cv_method = "kfold",
#'                      cvfits = cv_fits, nterms_max = 3, nclusters = 5,
#'                      nclusters_pred = 10, seed = 5555)
#'
#' # Stratified K-fold-CV is straightforward:
#' n_strat <- 3L
#' set.seed(692)
#' # Some example strata:
#' strat_fac <- sample(paste0("lvl", seq_len(n_strat)), size = nrow(dat_gauss),
#'                     replace = TRUE,
#'                     prob = diff(c(0, pnorm(seq_len(n_strat - 1L) - 0.5), 1)))
#' table(strat_fac)
#' # Use loo::kfold_split_stratified() to create the folds vector:
#' folds_strat <- loo::kfold_split_stratified(K = 2, x = strat_fac)
#' table(folds_strat, strat_fac)
#' # Call run_cvfun(), but this time with argument `folds` instead of `K` (here,
#' # specifying argument `seed` would not be necessary because of the set.seed()
#' # call above, but we specify it nonetheless for the sake of generality):
#' cv_fits_strat <- run_cvfun(ref, folds = folds_strat, seed = 391)
#' # Now use `cv_fits_strat` analogously to `cv_fits` from above.
#'
#' @export
run_cvfun <- function(object, ...) {
  UseMethod("run_cvfun")
}

#' @rdname run_cvfun
#' @export
run_cvfun.default <- function(object, ...) {
  refmodel <- get_refmodel(object, ...)
  return(run_cvfun(refmodel, ...))
}

#' @rdname run_cvfun
#' @export
run_cvfun.refmodel <- function(object,
                               K = if (!inherits(object, "datafit")) 5 else 10,
                               folds = NULL, seed = NA, ...) {
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
  stopifnot(!is.null(refmodel$cvfun))

  if (is.null(folds)) {
    folds <- cv_folds(refmodel$nobs, K = K)
  }
  if (getOption("projpred.warn_kfold_refits", TRUE)) {
    cvfits <- refmodel$cvfun(folds)
  } else {
    cvfits <- suppressWarnings(refmodel$cvfun(folds))
  }
  return(structure(cvfits, folds = folds))
}
