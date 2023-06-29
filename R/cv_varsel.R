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
#' the corresponding test set of each CV fold.
#'
#' @inheritParams varsel
#' @param cv_method The CV method, either `"LOO"` or `"kfold"`. In the `"LOO"`
#'   case, a Pareto-smoothed importance sampling leave-one-out CV (PSIS-LOO CV)
#'   is performed, which avoids refitting the reference model `nloo` times (in
#'   contrast to a standard LOO CV). In the `"kfold"` case, a \eqn{K}-fold CV is
#'   performed.
#' @param nloo **Caution:** Still experimental. Only relevant if `cv_method =
#'   "LOO"`. Number of subsampled LOO CV folds, i.e., number of observations
#'   used for the LOO CV (anything between 1 and the original number of
#'   observations). Smaller values lead to faster computation but higher
#'   uncertainty in the evaluation part. If `NULL`, all observations are used,
#'   but for faster experimentation, one can set this to a smaller value.
#' @param K Only relevant if `cv_method = "kfold"` and if the reference model
#'   was created with `cvfits` being `NULL` (which is the case for
#'   [get_refmodel.stanreg()] and [brms::get_refmodel.brmsfit()]). Number of
#'   folds in \eqn{K}-fold CV.
#' @param validate_search Only relevant if `cv_method = "LOO"`. A single logical
#'   value indicating whether to cross-validate also the search part, i.e.,
#'   whether to run the search separately for each CV fold (`TRUE`) or not
#'   (`FALSE`). We strongly do not recommend setting this to `FALSE`, because
#'   this is known to bias the predictive performance estimates of the selected
#'   submodels. However, setting this to `FALSE` can sometimes be useful because
#'   comparing the results to the case where this argument is `TRUE` gives an
#'   idea of how strongly the search is (over-)fitted to the data (the
#'   difference corresponds to the search degrees of freedom or the effective
#'   number of parameters introduced by the search).
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. If not
#'   `NA`, then the PRNG state is reset (to the state before calling
#'   [cv_varsel()]) upon exiting [cv_varsel()]. Here, `seed` is used for
#'   clustering the reference model's posterior draws (if `!is.null(nclusters)`
#'   or `!is.null(nclusters_pred)`), for subsampling LOO CV folds (if `nloo` is
#'   smaller than the number of observations), for sampling the folds in
#'   \eqn{K}-fold CV, and for drawing new group-level effects when predicting
#'   from a multilevel submodel (however, not yet in case of a GAMM).
#' @param parallel A single logical value indicating whether to run costly parts
#'   of the CV in parallel (`TRUE`) or not (`FALSE`). See also section "Note"
#'   below.
#'
#' @inherit varsel details return
#'
#' @note If `validate_search` is `FALSE`, the search is not included in the CV
#'   so that only a single full-data search is run.
#'
#'   For PSIS-LOO CV, \pkg{projpred} calls [loo::psis()] with `r_eff = NA`. This
#'   is only a problem if there was extreme autocorrelation between the MCMC
#'   iterations when the reference model was built. In those cases however, the
#'   reference model should not have been used anyway, so we don't expect
#'   \pkg{projpred}'s `r_eff = NA` to be a problem.
#'
#'   With `parallel = TRUE`, costly parts of \pkg{projpred}'s CV are run in
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
#' Magnusson, Måns, Michael Andersen, Johan Jonasson, and Aki Vehtari. 2019.
#' "Bayesian Leave-One-Out Cross-Validation for Large Data." In *Proceedings of
#' the 36th International Conference on Machine Learning*, edited by Kamalika
#' Chaudhuri and Ruslan Salakhutdinov, 97:4244--53. Proceedings of Machine
#' Learning Research. PMLR.
#' <https://proceedings.mlr.press/v97/magnusson19a.html>.
#'
#' Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. "Practical Bayesian Model
#' Evaluation Using Leave-One-Out Cross-Validation and WAIC." *Statistics and
#' Computing* 27 (5): 1413--32. \doi{10.1007/s11222-016-9696-4}.
#'
#' Vehtari, Aki, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah Gabry.
#' 2022. "Pareto Smoothed Importance Sampling." arXiv.
#' \doi{10.48550/arXiv.1507.02646}.
#'
#' @seealso [varsel()]
#'
#' @examplesIf identical(Sys.getenv("RUN_EX"), "true")
#' # Note: The code from this example is not executed when called via example().
#' # To execute it, you have to copy and paste it manually to the console.
#' if (requireNamespace("rstanarm", quietly = TRUE)) {
#'   # Data:
#'   dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)
#'
#'   # The "stanreg" fit which will be used as the reference model (with small
#'   # values for `chains` and `iter`, but only for technical reasons in this
#'   # example; this is not recommended in general):
#'   fit <- rstanarm::stan_glm(
#'     y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
#'     QR = TRUE, chains = 2, iter = 1000, refresh = 0, seed = 9876
#'   )
#'
#'   # Run cv_varsel() (with small values for `K`, `nterms_max`, `nclusters`,
#'   # and `nclusters_pred`, but only for the sake of speed in this example;
#'   # this is not recommended in general):
#'   cvvs <- cv_varsel(fit, cv_method = "kfold", K = 2, nterms_max = 3,
#'                     nclusters = 5, nclusters_pred = 10, seed = 5555)
#'   # Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
#'   # and `?ranking` for possible post-processing functions.
#' }
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
cv_varsel.refmodel <- function(
    object,
    method = NULL,
    cv_method = if (!inherits(object, "datafit")) "LOO" else "kfold",
    ndraws = NULL,
    nclusters = 20,
    ndraws_pred = 400,
    nclusters_pred = NULL,
    refit_prj = !inherits(object, "datafit"),
    nterms_max = NULL,
    penalty = NULL,
    verbose = TRUE,
    nloo = NULL,
    K = if (!inherits(object, "datafit")) 5 else 10,
    lambda_min_ratio = 1e-5,
    nlambda = 150,
    thresh = 1e-6,
    regul = 1e-4,
    validate_search = TRUE,
    seed = NA,
    search_terms = NULL,
    parallel = getOption("projpred.prll_cv", FALSE),
    ...
) {
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
    nterms_all = nterms_all
  )
  method <- args$method
  refit_prj <- args$refit_prj
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  search_terms <- args$search_terms
  # Parse arguments specific to cv_varsel():
  args <- parse_args_cv_varsel(
    refmodel = refmodel, cv_method = cv_method, K = K,
    validate_search = validate_search
  )
  cv_method <- args$cv_method
  K <- args$K
  # Arguments specific to the search:
  opt <- nlist(lambda_min_ratio, nlambda, thresh, regul)

  if (validate_search) {
    # Clustering or thinning for the final full-data search (already clustering
    # or thinning here for consistent PRNG states between the full-data search
    # in the `validate_search == FALSE` case and the full-data search in the
    # `validate_search == TRUE` case we are in here):
    p_sel <- get_refdist(refmodel, ndraws, nclusters)
  }

  if (cv_method == "LOO") {
    sel_cv <- loo_varsel(
      refmodel = refmodel, method = method, nterms_max = nterms_max,
      ndraws = ndraws, nclusters = nclusters, ndraws_pred = ndraws_pred,
      nclusters_pred = nclusters_pred, refit_prj = refit_prj, penalty = penalty,
      verbose = verbose, opt = opt, nloo = nloo,
      validate_search = validate_search, search_terms = search_terms,
      parallel = parallel, ...
    )
  } else if (cv_method == "kfold") {
    sel_cv <- kfold_varsel(
      refmodel = refmodel, method = method, nterms_max = nterms_max,
      ndraws = ndraws, nclusters = nclusters, ndraws_pred = ndraws_pred,
      nclusters_pred = nclusters_pred, refit_prj = refit_prj, penalty = penalty,
      verbose = verbose, opt = opt, K = K, search_terms = search_terms,
      parallel = parallel, ...
    )
  }

  if (validate_search) {
    verb_out("-----\nRunning a final search using the full dataset ...",
             verbose = verbose)
    search_path_full_data <- select(
      method = method, p_sel = p_sel, refmodel = refmodel,
      nterms_max = nterms_max, penalty = penalty, verbose = verbose, opt = opt,
      search_terms = search_terms, ...
    )
    verb_out("-----", verbose = verbose)
    ce_out <- rep(NA_real_, length(search_path_full_data$solution_terms) + 1L)
  } else {
    search_path_full_data <- sel_cv$search_path
    ce_out <- sel_cv$ce
  }

  # Defined here for `nobs_test` later:
  y_wobs_test <- sel_cv$y_wobs_test

  # Just a dummy object which is not used as usual, but only for inferring the
  # output elements `clust_used_eval` and `nprjdraws_eval` (this get_refdist()
  # call is much cheaper than calling varsel() with its re-projections (if
  # `refit_prj = TRUE`) instead of select() above in the case
  # `if (validate_search)`, see GitHub PR #385):
  if (refit_prj) {
    refdist_eval_dummy <- get_refdist(refmodel, ndraws_pred, nclusters_pred)
  } else {
    refdist_eval_dummy <- search_path_full_data$p_sel
  }

  # The object to be returned:
  vs <- nlist(refmodel,
              nobs_train = refmodel$nobs,
              search_path = search_path_full_data,
              solution_terms = search_path_full_data$solution_terms,
              solution_terms_cv = sel_cv$solution_terms_cv,
              ce = ce_out,
              type_test = cv_method,
              y_wobs_test,
              nobs_test = nrow(y_wobs_test),
              summaries = sel_cv$summaries,
              nterms_all,
              nterms_max,
              method,
              cv_method,
              K = K,
              validate_search,
              clust_used_search = search_path_full_data$p_sel$clust_used,
              clust_used_eval = refdist_eval_dummy$clust_used,
              nprjdraws_search = NCOL(search_path_full_data$p_sel$mu),
              nprjdraws_eval = NCOL(refdist_eval_dummy$mu),
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
# @param K See argument `K` of cv_varsel().
# @param validate_search See argument `validate_search` of cv_varsel().
#
# @return A list with elements `cv_method` and `K`.
parse_args_cv_varsel <- function(refmodel, cv_method, K, validate_search) {
  stopifnot(!is.null(cv_method))
  if (cv_method == "loo") {
    cv_method <- toupper(cv_method)
  }
  if (!cv_method %in% c("kfold", "LOO")) {
    stop("Unknown `cv_method`.")
  }
  if (cv_method == "LOO" && inherits(refmodel, "datafit")) {
    warning("For an `object` of class \"datafit\", `cv_method` is ",
            "automatically set to \"kfold\".")
    cv_method <- "kfold"
  }

  if (cv_method == "kfold") {
    if (!is.null(refmodel$cvfits)) {
      K <- length(refmodel$cvfits$fits)
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
    if (!validate_search) {
      stop("`cv_method = \"kfold\"` cannot be used with ",
           "`validate_search = FALSE`.")
    }
  } else {
    K <- NULL
  }

  return(nlist(cv_method, K))
}

# PSIS-LOO CV -------------------------------------------------------------

# Workhorse function for a variable selection with PSIS-LOO CV
#
# Argument `validate_search` indicates whether the search is performed
# separately for each LOO CV fold (i.e., separately for each observation). For
# all other arguments, see the documentation of cv_varsel().
loo_varsel <- function(refmodel, method, nterms_max, ndraws,
                       nclusters, ndraws_pred, nclusters_pred, refit_prj,
                       penalty, verbose, opt, nloo, validate_search,
                       search_terms, parallel, ...) {
  ## Pre-processing ---------------------------------------------------------

  # Clustering or thinning for the search (note that in case of
  # `validate_search = TRUE`, only `cl_sel` is used later, not `p_sel` itself):
  p_sel <- get_refdist(refmodel, ndraws = ndraws, nclusters = nclusters)
  cl_sel <- p_sel$cl
  # Clustering or thinning for the performance evaluation (note that in case of
  # `validate_search = TRUE`, only `cl_pred` is used later, not `p_pred`
  # itself):
  p_pred <- get_refdist(refmodel, ndraws = ndraws_pred,
                        nclusters = nclusters_pred)
  cl_pred <- p_pred$cl

  if (inherits(refmodel, "datafit")) {
    stop("LOO can be performed only if the reference model is a genuine ",
         "probabilistic model for which the log-likelihood can be evaluated.")
  }

  # Log-likelihood values for the reference model (necessary for the PSIS-LOO CV
  # weights, but also for performance statistics like ELPD and MLPD):
  if (refmodel$family$for_latent) {
    mu_offs_oscale <- refmodel$family$latent_ilink(
      t(refmodel$mu_offs), cl_ref = seq_along(refmodel$wdraws_ref),
      wdraws_ref = refmodel$wdraws_ref
    )
    if (length(dim(mu_offs_oscale)) < 2) {
      stop("Unexpected structure for the output of `latent_ilink`.")
    }
    loglik_forPSIS <- refmodel$family$latent_ll_oscale(
      mu_offs_oscale, y_oscale = refmodel$y_oscale, wobs = refmodel$wobs,
      cl_ref = seq_along(refmodel$wdraws_ref), wdraws_ref = refmodel$wdraws_ref
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
      n_aug <- nrow(mu_offs_oscale)
    } else {
      # In this case, `mu_offs_oscale` is a matrix (S x N). Transposing it to an
      # N x S matrix would be more consistent with projpred's internal
      # convention, but avoiding the transposition is computationally more
      # efficient:
      n_aug <- ncol(mu_offs_oscale)
    }
  } else {
    loglik_forPSIS <- t(refmodel$family$ll_fun(
      refmodel$mu_offs, refmodel$dis, refmodel$y, refmodel$wobs
    ))
  }
  n <- ncol(loglik_forPSIS)

  # PSIS-LOO CV weights:
  psisloo <- loo::psis(-loglik_forPSIS, cores = 1, r_eff = NA)
  lw <- weights(psisloo)
  # pareto_k <- loo::pareto_k_values(psisloo)

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
    if (any(y_lat_E$pareto_k > 0.7)) {
      warning("In the recalculation of the latent response values, ",
              sum(y_lat_E$pareto_k > 0.7), " (of ", n, ") Pareto k-value(s) ",
              "exceeded the threshold of 0.7.")
    }
    refmodel$y <- y_lat_E$value
  }

  # LOO subsampling (by default, don't subsample, but use all observations):
  nloo <- min(nloo, n)
  if (nloo < 1) {
    stop("nloo must be at least 1")
  } else if (nloo < n) {
    warning("Subsampled LOO CV is still experimental.")
  }
  # validset <- loo_subsample(n, nloo, pareto_k)
  loo_ref_oscale <- apply(loglik_forPSIS + lw, 2, log_sum_exp)
  validset <- loo_subsample_pps(nloo, loo_ref_oscale)
  inds <- validset$inds

  # Initialize objects where to store the results:
  loo_sub <- replicate(nterms_max + 1L, rep(NA, n), simplify = FALSE)
  mu_sub <- replicate(
    nterms_max + 1L,
    structure(rep(NA, nrow(refmodel$mu_offs)),
              nobs_orig = attr(refmodel$mu_offs, "nobs_orig"),
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
                  nobs_orig = n,
                  class = "augvec"),
        simplify = FALSE
      )
    } else {
      mu_sub_oscale <- mu_sub
    }
  }

  if (!validate_search) {
    ## Case `validate_search = FALSE` -----------------------------------------

    verb_out("-----\nRunning the search using the full dataset ...",
             verbose = verbose)
    search_path <- select(
      method = method, p_sel = p_sel, refmodel = refmodel,
      nterms_max = nterms_max, penalty = penalty, verbose = verbose, opt = opt,
      search_terms = search_terms, ...
    )
    verb_out("-----", verbose = verbose)

    verb_out("-----\nFor performance evaluation: Re-projecting (using the ",
             "full dataset) onto the submodels along the full-data solution ",
             "path ...", verbose = verbose && refit_prj)
    submodls <- get_submodls(
      search_path = search_path,
      nterms = c(0, seq_along(search_path$solution_terms)), p_ref = p_pred,
      refmodel = refmodel, regul = opt$regul, refit_prj = refit_prj, ...
    )
    verb_out("-----", verbose = verbose && refit_prj)

    verb_out("-----\nCalculating the full-data performance evaluation ",
             "results and weighting those results according to the PSIS-LOO ",
             "CV weights ...", verbose = verbose)
    if (refit_prj) {
      refdist_eval <- p_pred
    } else {
      refdist_eval <- p_sel
    }
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
        refdist_eval_mu_offs_oscale, y_oscale = refmodel$y_oscale[inds],
        wobs = refmodel$wobs[inds], cl_ref = refdist_eval$cl,
        wdraws_ref = refdist_eval$wdraws_orig
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
    sub_psisloo <- suppressWarnings(
      loo::psis(-log_lik_ref, cores = 1, r_eff = NA)
    )
    lw_sub <- suppressWarnings(weights(sub_psisloo))
    # Take into account that clustered draws usually have different weights:
    lw_sub <- lw_sub + log(refdist_eval$wdraws_prj)
    # This re-weighting requires a re-normalization (as.array() is applied to
    # have stricter consistency checks, see `?sweep`):
    lw_sub <- sweep(lw_sub, 2, as.array(apply(lw_sub, 2, log_sum_exp)))
    for (k in seq_along(submodls)) {
      mu_k <- refmodel$family$mu_fun(submodls[[k]]$outdmin,
                                     obs = inds,
                                     offset = refmodel$offset[inds])
      log_lik_sub <- t(refmodel$family$ll_fun(
        mu_k, submodls[[k]]$dis, refmodel$y[inds], refmodel$wobs[inds]
      ))
      loo_sub[[k]][inds] <- apply(log_lik_sub + lw_sub, 2, log_sum_exp)
      if (refmodel$family$for_latent) {
        mu_k_oscale <- refmodel$family$latent_ilink(
          t(mu_k), cl_ref = refdist_eval$cl,
          wdraws_ref = refdist_eval$wdraws_orig
        )
        log_lik_sub_oscale <- refmodel$family$latent_ll_oscale(
          mu_k_oscale, y_oscale = refmodel$y_oscale[inds],
          wobs = refmodel$wobs[inds], cl_ref = refdist_eval$cl,
          wdraws_ref = refdist_eval$wdraws_orig
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
        mu_sub[[k]][i_flx] <- mu_k[run_index_flx, ] %*% exp(lw_sub[, run_index])
        if (refmodel$family$for_latent) {
          if (inherits(mu_k_oscale, "augmat")) {
            mu_sub_oscale[[k]][i_aug] <- mu_k_oscale[run_index_aug, ] %*%
              exp(lw_sub[, run_index])
          } else {
            # In principle, we could use the same code for averaging across the
            # draws as above in the `"augmat"` case. However, that would require
            # `mu_k_oscale <- t(mu_k_oscale)` beforehand, so the following
            # should be more efficient:
            mu_sub_oscale[[k]][i_aug] <- exp(lw_sub[, run_index]) %*%
              mu_k_oscale[, run_index_aug]
          }
        }
      }
    }
    verb_out("-----", verbose = verbose)
  } else {
    ## Case `validate_search = TRUE` ------------------------------------------

    verb_out("-----\nRunning the search and the performance evaluation for ",
             "each of the N = ", nloo, " LOO CV folds separately ...",
             verbose = verbose)
    one_obs <- function(run_index,
                        verbose_search = verbose &&
                          getOption("projpred.extra_verbose", FALSE),
                        ...) {
      # Observation index:
      i <- inds[run_index]

      # Reweight the clusters (or thinned draws) according to the PSIS weights:
      p_sel <- get_p_clust(
        family = refmodel$family, eta = refmodel$eta, mu = refmodel$mu,
        mu_offs = refmodel$mu_offs, dis = refmodel$dis, wdraws = exp(lw[, i]),
        cl = cl_sel
      )
      p_pred <- get_p_clust(
        family = refmodel$family, eta = refmodel$eta, mu = refmodel$mu,
        mu_offs = refmodel$mu_offs, dis = refmodel$dis, wdraws = exp(lw[, i]),
        cl = cl_pred
      )

      # Run the search with the reweighted clusters (or thinned draws) (so the
      # *reweighted* fitted response values from the reference model act as
      # artifical response values in the projection (or L1-penalized
      # projection)):
      search_path <- select(
        method = method, p_sel = p_sel, refmodel = refmodel,
        nterms_max = nterms_max, penalty = penalty, verbose = verbose_search,
        opt = opt, search_terms = search_terms, ...
      )

      # Re-project along the solution path (or fetch the projections from the
      # search results) of the current fold:
      submodls <- get_submodls(
        search_path = search_path,
        nterms = c(0, seq_along(search_path$solution_terms)),
        p_ref = p_pred, refmodel = refmodel, regul = opt$regul,
        refit_prj = refit_prj, ...
      )
      # Predictive performance at the omitted observation:
      summaries_sub <- get_sub_summaries(submodls = submodls,
                                         refmodel = refmodel,
                                         test_points = i)

      return(nlist(predictor_ranking = search_path[["solution_terms"]],
                   summaries_sub))
    }
    if (!parallel) {
      # Sequential case. Actually, we could simply use ``%do_projpred%` <-
      # foreach::`%do%`` here and then proceed as in the parallel case, but that
      # would require adding more "hard" dependencies (because packages
      # 'foreach' and 'doRNG' would have to be moved from `Suggests:` to
      # `Imports:`).
      if (verbose) {
        pb <- utils::txtProgressBar(min = 0, max = nloo, style = 3, initial = 0)
      }
      res_cv <- lapply(seq_along(inds), function(run_index) {
        if (verbose) {
          on.exit(utils::setTxtProgressBar(pb, run_index))
        }
        one_obs(run_index, ...)
      })
      if (verbose) {
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
      dot_args <- list(...)
      `%do_projpred%` <- doRNG::`%dorng%`
      res_cv <- foreach::foreach(
        run_index = seq_along(inds),
        .export = c("one_obs", "dot_args"),
        .noexport = c("p_sel", "p_pred", "mu_offs_oscale", "loglik_forPSIS",
                      "psisloo", "y_lat_E", "loo_ref_oscale", "validset",
                      "loo_sub", "mu_sub", "loo_sub_oscale", "mu_sub_oscale")
      ) %do_projpred% {
        do.call(one_obs, c(list(run_index = run_index, verbose_search = FALSE),
                           dot_args))
      }
    }
    # For storing the fold-wise solution paths:
    solution_terms_mat <- matrix(nrow = n, ncol = nterms_max)
    # For checking that the length of the predictor ranking is the same across
    # all CV folds (and also for cutting off `solution_terms_mat` later):
    prv_len_soltrms <- NULL
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
      if (is.null(prv_len_soltrms)) {
        prv_len_soltrms <- length(rk_i)
      } else {
        stopifnot(identical(length(rk_i), prv_len_soltrms))
      }
      solution_terms_mat[i, seq_along(rk_i)] <- rk_i
    }
    verb_out("-----", verbose = verbose)
  }

  ## Post-processing --------------------------------------------------------

  # Submodel predictive performance:
  summ_sub <- lapply(seq_len(nterms_max + 1L), function(k) {
    summ_k <- list(lppd = loo_sub[[k]], mu = mu_sub[[k]], wcv = validset$wcv)
    if (refmodel$family$for_latent) {
      summ_k$oscale <- list(lppd = loo_sub_oscale[[k]], mu = mu_sub_oscale[[k]],
                            wcv = validset$wcv)
    }
    return(summ_k)
  })

  # Reference model predictive performance:
  if (formula_contains_group_terms(refmodel$formula) &&
      getOption("projpred.mlvl_pred_new", FALSE)) {
    # Need to use `mlvl_allrandom = TRUE` (`refmodel$mu_offs` is based on
    # `mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)`):
    eta_offs_mlvlRan <- refmodel$ref_predfun(refmodel$fit, excl_offs = FALSE)
    mu_offs_mlvlRan <- refmodel$family$linkinv(eta_offs_mlvlRan)
  } else {
    mu_offs_mlvlRan <- refmodel$mu_offs
  }
  mu_ref <- do.call(c, lapply(seq_len(nrow(mu_offs_mlvlRan)), function(i) {
    # For the augmented-data projection, `mu_offs_mlvlRan` is an augmented-rows
    # matrix whereas the columns of `lw` refer to the original (non-augmented)
    # observations. Since `i` refers to the rows of `mu_offs_mlvlRan`, the index
    # for `lw` needs to be adapted:
    i_nonaug <- i %% n
    if (i_nonaug == 0) {
      i_nonaug <- n
    }
    mu_offs_mlvlRan[i, ] %*% exp(lw[, i_nonaug])
  }))
  mu_ref <- structure(
    mu_ref,
    nobs_orig = attr(mu_offs_mlvlRan, "nobs_orig"),
    class = sub("augmat", "augvec", oldClass(mu_offs_mlvlRan), fixed = TRUE)
  )
  if (refmodel$family$for_latent) {
    loglik_lat <- t(refmodel$family$ll_fun(
      mu_offs_mlvlRan, refmodel$dis, refmodel$y, refmodel$wobs
    ))
    lppd_ref <- apply(loglik_lat + lw, 2, log_sum_exp)
  } else {
    if (formula_contains_group_terms(refmodel$formula) &&
        getOption("projpred.mlvl_pred_new", FALSE)) {
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
    if (formula_contains_group_terms(refmodel$formula) &&
        getOption("projpred.mlvl_pred_new", FALSE)) {
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
    mu_ref_oscale <- do.call(c, lapply(seq_len(n_aug), function(i) {
      i_nonaug <- i %% n
      if (i_nonaug == 0) {
        i_nonaug <- n
      }
      if (inherits(mu_offs_mlvlRan_oscale, "augmat")) {
        return(mu_offs_mlvlRan_oscale[i, ] %*% exp(lw[, i_nonaug]))
      } else {
        # In principle, we could use the same code for averaging across the
        # draws as above in the `"augmat"` case. However, that would require
        # `mu_offs_mlvlRan_oscale <- t(mu_offs_mlvlRan_oscale)` beforehand, so
        # the following should be more efficient:
        return(exp(lw[, i_nonaug]) %*% mu_offs_mlvlRan_oscale[, i])
      }
    }))
    mu_ref_oscale <- structure(
      mu_ref_oscale,
      nobs_orig = attr(mu_offs_mlvlRan_oscale, "nobs_orig"),
      class = sub("augmat", "augvec", oldClass(mu_offs_mlvlRan_oscale),
                  fixed = TRUE)
    )
    if (formula_contains_group_terms(refmodel$formula) &&
        getOption("projpred.mlvl_pred_new", FALSE)) {
      # Need to use `mlvl_allrandom = TRUE` (`loo_ref_oscale` is based on
      # `mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)`):
      loglik_mlvlRan <- refmodel$family$latent_ll_oscale(
        mu_offs_mlvlRan_oscale_odim, y_oscale = refmodel$y_oscale,
        wobs = refmodel$wobs, cl_ref = seq_along(refmodel$wdraws_ref),
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
    out_list <- nlist(search_path, ce = sapply(submodls, "[[", "ce"))
  } else {
    out_list <- nlist(solution_terms_cv = solution_terms_mat[
      , seq_len(prv_len_soltrms), drop = FALSE
    ])
  }
  out_list <- c(out_list,
                nlist(summaries,
                      y_wobs_test = as.data.frame(refmodel[nms_y_wobs_test()])))
  return(out_list)
}

# K-fold CV ---------------------------------------------------------------

# Needed to avoid a NOTE in `R CMD check`:
if (getRversion() >= package_version("2.15.1")) {
  utils::globalVariables("list_cv_k")
}

kfold_varsel <- function(refmodel, method, nterms_max, ndraws,
                         nclusters, ndraws_pred, nclusters_pred,
                         refit_prj, penalty, verbose, opt, K,
                         search_terms, parallel, ...) {
  # Fetch the K reference model fits (or fit them now if not already done) and
  # create objects of class `refmodel` from them (and also store the `omitted`
  # indices):
  list_cv <- get_kfold(refmodel, K, verbose)
  K <- length(list_cv)

  if (refmodel$family$for_latent) {
    # Need to set the latent response values in `refmodel$y` to `NA`s because
    # `refmodel$y` resulted from applying `colMeans(posterior_linpred())` to the
    # original (full-data) reference model fit, so using the `fold$omitted`
    # subset of `refmodel$y` as (latent) response values in fold k of K would
    # induce a dependency between training and test data:
    refmodel$y <- rep(NA, refmodel$nobs)
  }
  y_wobs_test <- as.data.frame(refmodel[nms_y_wobs_test()])

  verb_out("-----\nRunning the search and the performance evaluation for ",
           "each of the K = ", K, " CV folds separately ...", verbose = verbose)
  one_fold <- function(fold,
                       verbose_search = verbose &&
                         getOption("projpred.extra_verbose", FALSE),
                       ...) {
    # Run the search for the current fold:
    p_sel <- get_refdist(fold$refmodel, ndraws, nclusters)
    search_path <- select(
      method = method, p_sel = p_sel, refmodel = fold$refmodel,
      nterms_max = nterms_max, penalty = penalty, verbose = verbose_search,
      opt = opt, search_terms = search_terms, ...
    )

    # For performance evaluation: Re-project (using the training data of the
    # current fold) along the predictor ranking (or fetch the projections from
    # the search output) of the current fold:
    p_pred <- get_refdist(fold$refmodel, ndraws_pred, nclusters_pred)
    submodls <- get_submodls(
      search_path = search_path,
      nterms = c(0, seq_along(search_path$solution_terms)),
      p_ref = p_pred, refmodel = fold$refmodel, regul = opt$regul,
      refit_prj = refit_prj, ...
    )

    # Performance evaluation for the re-projected or fetched submodels of the
    # current fold:
    summaries_sub <- get_sub_summaries(submodls = submodls,
                                       refmodel = refmodel,
                                       test_points = fold$omitted)

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

    return(nlist(predictor_ranking = search_path[["solution_terms"]],
                 summaries_sub, summaries_ref))
  }
  if (!parallel) {
    # Sequential case. Actually, we could simply use ``%do_projpred%` <-
    # foreach::`%do%`` here and then proceed as in the parallel case, but that
    # would require adding more "hard" dependencies (because packages 'foreach'
    # and 'doRNG' would have to be moved from `Suggests:` to `Imports:`).
    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial = 0)
    }
    res_cv <- lapply(seq_along(list_cv), function(k) {
      if (verbose) {
        on.exit(utils::setTxtProgressBar(pb, k))
      }
      one_fold(list_cv[[k]], ...)
    })
    if (verbose) {
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
    dot_args <- list(...)
    `%do_projpred%` <- doRNG::`%dorng%`
    res_cv <- foreach::foreach(
      list_cv_k = list_cv,
      .export = c("one_fold", "dot_args"),
      .noexport = c("list_cv")
    ) %do_projpred% {
      do.call(one_fold, c(list(fold = list_cv_k, verbose_search = FALSE),
                          dot_args))
    }
  }
  verb_out("-----", verbose = verbose)
  solution_terms_cv <- do.call(rbind, lapply(res_cv, "[[", "predictor_ranking"))

  # Handle the submodels' performance evaluation results:
  sub_foldwise <- lapply(res_cv, "[[", "summaries_sub")
  if (getRversion() >= package_version("4.2.0")) {
    sub_foldwise <- simplify2array(sub_foldwise, higher = FALSE, except = NULL)
  } else {
    sub_foldwise <- simplify2array(sub_foldwise, higher = FALSE)
    if (is.null(dim(sub_foldwise))) {
      sub_dim <- dim(solution_terms_cv)
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

    # Add fold-specific weights (see the discussion at GitHub issue #94 for why
    # this might have to be changed):
    summ$wcv <- rep(1, length(summ$lppd))
    summ$wcv <- summ$wcv / sum(summ$wcv)

    if (!is.null(summ$oscale)) {
      summ$oscale$mu <- summ$oscale$mu[order(idxs_sorted_by_fold_aug)]
      summ$oscale$lppd <- summ$oscale$lppd[order(idxs_sorted_by_fold)]
      summ$oscale$wcv <- summ$wcv
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

  return(nlist(solution_terms_cv, summaries = nlist(sub, ref), y_wobs_test))
}

# Re-fit the reference model K times (once for each fold; `cvfun` case) or fetch
# the K reference model fits if already computed (`cvfits` case). This function
# will return a list of length K, where each element is a list with elements
# `refmodel` (output of init_refmodel()) and `omitted` (vector of indices of
# those observations which were left out for the corresponding fold).
get_kfold <- function(refmodel, K, verbose) {
  if (is.null(refmodel$cvfits)) {
    if (!is.null(refmodel$cvfun)) {
      # In this case, cvfun() provided (and `cvfits` not), so run cvfun() now.
      if (verbose && !inherits(refmodel, "datafit")) {
        verb_out("-----\nRefitting the reference model K = ", K, " times ",
                 "(using the fold-wise training data) ...")
      }
      folds <- cv_folds(refmodel$nobs, K = K,
                        seed = sample.int(.Machine$integer.max, 1))
      cvfits <- refmodel$cvfun(folds)
      verb_out("-----", verbose = verbose)
    } else {
      stop("For a reference model which is not of class `datafit`, either ",
           "`cvfits` or `cvfun` needs to be provided for K-fold CV (see ",
           "`?init_refmodel`).")
    }
  } else {
    folds <- attr(refmodel$cvfits, "folds")
    cvfits <- refmodel$cvfits$fits
  }
  return(lapply(seq_len(K), function(k) {
    cvfit <- cvfits[[k]]
    # Add the omitted observation indices for this fold:
    cvfit$omitted <- which(folds == k)
    # Add the fold index:
    cvfit$projpred_k <- k
    return(list(refmodel = refmodel$cvrefbuilder(cvfit),
                omitted = cvfit$omitted))
  }))
}

# ## decide which points to go through in the validation (i.e., which points
# ## belong to the semi random subsample of validation points)
# loo_subsample <- function(n, nloo, pareto_k) {
#   # Note: A seed is not set here because this function is not exported and has
#   # a calling stack at the beginning of which a seed is set.
#
#   resample <- function(x, ...) x[sample.int(length(x), ...)]
#
#   if (nloo < n) {
#     bad <- which(pareto_k > 0.7)
#     ok <- which(pareto_k <= 0.7 & pareto_k > 0.5)
#     good <- which(pareto_k <= 0.5)
#     inds <- resample(bad, min(length(bad), floor(nloo / 3)))
#     inds <- c(inds, resample(ok, min(length(ok), floor(nloo / 3))))
#     inds <- c(inds, resample(good, min(length(good), floor(nloo / 3))))
#     if (length(inds) < nloo) {
#       ## not enough points selected, so choose randomly among the rest
#       inds <- c(inds, resample(setdiff(seq_len(n), inds), nloo - length(inds)))
#     }
#
#     ## assign the weights corresponding to this stratification (for example,
#     ## the 'bad' values are likely to be overpresented in the sample)
#     wcv <- rep(0, n)
#     wcv[inds[inds %in% bad]] <- length(bad) / sum(inds %in% bad)
#     wcv[inds[inds %in% ok]] <- length(ok) / sum(inds %in% ok)
#     wcv[inds[inds %in% good]] <- length(good) / sum(inds %in% good)
#   } else {
#     ## all points used
#     inds <- seq_len(n)
#     wcv <- rep(1, n)
#   }
#
#   ## ensure weights are normalized
#   wcv <- wcv / sum(wcv)
#
#   return(nlist(inds, wcv))
# }

## decide which points to go through in the validation based on
## proportional-to-size subsampling as implemented in Magnusson, M., Riis
## Andersen, M., Jonasson, J. and Vehtari, A. (2019). Leave-One-Out
## Cross-Validation for Large Data. In International Conference on Machine
## Learning.
loo_subsample_pps <- function(nloo, lppd) {
  # Note: A seed is not set here because this function is not exported and has a
  # calling stack at the beginning of which a seed is set.

  if (nloo > length(lppd)) {
    stop("Argument `nloo` must not be larger than the number of observations.")
  } else if (nloo == length(lppd)) {
    inds <- seq_len(nloo)
    wcv <- rep(1, nloo)
  } else if (nloo < length(lppd)) {
    wcv <- exp(lppd - max(lppd))
    inds <- sample(seq_along(lppd), size = nloo, prob = wcv)
  }
  wcv <- wcv / sum(wcv)

  return(nlist(inds, wcv))
}
