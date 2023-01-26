#' Variable selection with cross-validation
#'
#' Run the *search* part and the *evaluation* part for a projection predictive
#' variable selection. The search part determines the solution path, i.e., the
#' best submodel for each submodel size (number of predictor terms). The
#' evaluation part determines the predictive performance of the submodels along
#' the solution path. In contrast to [varsel()], [cv_varsel()] performs a
#' cross-validation (CV) by running the search part with the training data of
#' each CV fold separately (an exception is explained in section "Note" below)
#' and running the evaluation part on the corresponding test set of each CV
#' fold.
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
#'   idea of how strongly the variable selection is (over-)fitted to the data
#'   (the difference corresponds to the search degrees of freedom or the
#'   effective number of parameters introduced by the search).
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. Here,
#'   this seed is used for clustering the reference model's posterior draws (if
#'   `!is.null(nclusters)` or `!is.null(nclusters_pred)`), for subsampling LOO
#'   CV folds (if `nloo` is smaller than the number of observations), for
#'   sampling the folds in K-fold CV, and for drawing new group-level effects
#'   when predicting from a multilevel submodel (however, not yet in case of a
#'   GAMM).
#'
#' @inherit varsel details return
#'
#' @note The case `cv_method == "LOO" && !validate_search` constitutes an
#'   exception where the search part is not cross-validated. In that case, the
#'   evaluation part is based on a PSIS-LOO CV also for the submodels.
#'
#'   For all PSIS-LOO CVs, \pkg{projpred} calls [loo::psis()] with `r_eff = NA`.
#'   This is only a problem if there was extreme autocorrelation between the
#'   MCMC iterations when the reference model was built. In those cases however,
#'   the reference model should not have been used anyway, so we don't expect
#'   \pkg{projpred}'s `r_eff = NA` to be a problem.
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
#'     QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
#'   )
#'
#'   # Variable selection with cross-validation (with small values
#'   # for `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the
#'   # sake of speed in this example; this is not recommended in general):
#'   cvvs <- cv_varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                     seed = 5555)
#'   # Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
#'   # and `?solution_terms.vsel` for possible post-processing functions.
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
    seed = sample.int(.Machine$integer.max, 1),
    search_terms = NULL,
    ...
) {
  # Set seed, but ensure the old RNG state is restored on exit:
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  if (!is.na(seed)) set.seed(seed)

  refmodel <- object
  # Needed to avoid a warning when calling varsel() later:
  search_terms_usr <- search_terms
  ## resolve the arguments similar to varsel
  args <- parse_args_varsel(
    refmodel = refmodel, method = method, refit_prj = refit_prj,
    nterms_max = nterms_max, nclusters = nclusters, search_terms = search_terms
  )
  method <- args$method
  refit_prj <- args$refit_prj
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  search_terms <- args$search_terms

  ## arguments specific to this function
  args <- parse_args_cv_varsel(
    refmodel = refmodel, cv_method = cv_method, K = K,
    validate_search = validate_search
  )
  cv_method <- args$cv_method
  K <- args$K

  ## search options
  opt <- nlist(lambda_min_ratio, nlambda, thresh, regul)

  if (cv_method == "LOO") {
    sel_cv <- loo_varsel(
      refmodel = refmodel, method = method, nterms_max = nterms_max,
      ndraws = ndraws, nclusters = nclusters, ndraws_pred = ndraws_pred,
      nclusters_pred = nclusters_pred, refit_prj = refit_prj, penalty = penalty,
      verbose = verbose, opt = opt, nloo = nloo,
      validate_search = validate_search, search_terms = search_terms, ...
    )
  } else if (cv_method == "kfold") {
    sel_cv <- kfold_varsel(
      refmodel = refmodel, method = method, nterms_max = nterms_max,
      ndraws = ndraws, nclusters = nclusters, ndraws_pred = ndraws_pred,
      nclusters_pred = nclusters_pred, refit_prj = refit_prj, penalty = penalty,
      verbose = verbose, opt = opt, K = K, search_terms = search_terms, ...
    )
  }

  if (validate_search || cv_method == "kfold") {
    ## run the selection using the full dataset
    if (verbose) {
      print(paste("Performing the selection using all the data.."))
    }
    sel <- varsel(refmodel,
                  method = method, ndraws = ndraws, nclusters = nclusters,
                  ndraws_pred = ndraws_pred, nclusters_pred = nclusters_pred,
                  refit_prj = refit_prj, nterms_max = nterms_max - 1,
                  penalty = penalty, verbose = verbose,
                  lambda_min_ratio = lambda_min_ratio, nlambda = nlambda,
                  regul = regul, search_terms = search_terms_usr, seed = seed,
                  ...)
  } else {
    sel <- sel_cv$sel
  }

  # Create `pct_solution_terms_cv`, a summary table of the fold-wise solution
  # paths. For the column names (and therefore the order of the solution terms
  # in the columns), the solution path from the full-data search is used. Note
  # that the following code assumes that all CV folds have equal weight.
  pct_solution_terms_cv <- cbind(
    size = seq_len(ncol(sel_cv$solution_terms_cv)),
    do.call(cbind, lapply(
      setNames(nm = sel$solution_terms),
      function(soltrm_k) {
        colMeans(sel_cv$solution_terms_cv == soltrm_k, na.rm = TRUE)
      }
    ))
  )

  ## create the object to be returned
  vs <- nlist(refmodel,
              search_path = sel$search_path,
              d_test = sel_cv$d_test,
              summaries = sel_cv$summaries,
              ce = sel$ce,
              solution_terms = sel$solution_terms,
              pct_solution_terms_cv,
              nterms_all = count_terms_in_formula(refmodel$formula),
              nterms_max,
              method,
              cv_method,
              validate_search,
              clust_used_search = sel$clust_used_search,
              clust_used_eval = sel$clust_used_eval,
              nprjdraws_search = sel$nprjdraws_search,
              nprjdraws_eval = sel$nprjdraws_eval)
  class(vs) <- "vsel"
  if (verbose) {
    print("Done.")
  }
  return(vs)
}

# Auxiliary function for parsing the input arguments for specific cv_varsel.
# This is similar in spirit to parse_args_varsel(), that is, to avoid the main
# function to become too long and complicated to maintain.
#
# @param refmodel Reference model as extracted by get_refmodel
# @param cv_method The cross-validation method, either 'LOO' or 'kfold'.
#   Default is 'LOO'.
# @param K Number of folds in the K-fold cross validation. Default is 5 for
#   genuine reference models and 10 for datafits (that is, for penalized
#   maximum likelihood estimation).
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
    stopifnot(!is.null(K))
    if (length(K) > 1 || !is.numeric(K) || !.is.wholenumber(K)) {
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
  }

  return(nlist(cv_method, K))
}

loo_varsel <- function(refmodel, method, nterms_max, ndraws,
                       nclusters, ndraws_pred, nclusters_pred, refit_prj,
                       penalty, verbose, opt, nloo = NULL,
                       validate_search = TRUE, search_terms = NULL, ...) {
  ##
  ## Perform the validation of the searching process using LOO. validate_search
  ## indicates whether the selection is performed separately for each fold (for
  ## each data point)
  ##

  # Pre-processing ----------------------------------------------------------

  eta <- refmodel$eta
  mu <- refmodel$mu
  mu_offs <- refmodel$mu_offs
  dis <- refmodel$dis
  ## the clustering/subsampling used for selection
  p_sel <- .get_refdist(refmodel, ndraws = ndraws, nclusters = nclusters)
  cl_sel <- p_sel$cl # clustering information

  ## the clustering/subsampling used for prediction
  p_pred <- .get_refdist(refmodel, ndraws = ndraws_pred,
                         nclusters = nclusters_pred)
  cl_pred <- p_pred$cl

  ## fetch the log-likelihood for the reference model to obtain the LOO weights
  if (inherits(refmodel, "datafit")) {
    ## case where log-likelihood not available, i.e., the reference model is not
    ## a genuine model => cannot compute LOO
    stop("LOO can be performed only if the reference model is a genuine ",
         "probabilistic model for which the log-likelihood can be evaluated.")
  }

  if (refmodel$family$for_latent) {
    mu_offs_oscale <- refmodel$family$latent_ilink(
      t(mu_offs), cl_ref = seq_along(refmodel$wsample),
      wdraws_ref = refmodel$wsample
    )
    if (length(dim(mu_offs_oscale)) < 2) {
      stop("Unexpected structure for the output of `latent_ilink`.")
    }
    loglik_forPSIS <- refmodel$family$latent_ll_oscale(
      mu_offs_oscale, y_oscale = refmodel$y_oscale, wobs = refmodel$wobs,
      cl_ref = seq_along(refmodel$wsample), wdraws_ref = refmodel$wsample
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
      mu_offs, dis, refmodel$y, refmodel$wobs
    ))
  }
  n <- ncol(loglik_forPSIS)
  psisloo <- loo::psis(-loglik_forPSIS, cores = 1, r_eff = NA)
  lw <- weights(psisloo)
  # pareto_k <- loo::pareto_k_values(psisloo)

  if (refmodel$family$for_latent) {
    # Need to re-calculate the latent response values in `refmodel$y` by
    # incorporating the PSIS weights because `refmodel$y` resulted from applying
    # `colMeans(posterior_linpred())` to the original (full-data) reference
    # model fit, so using `refmodel$y` would induce a dependency between
    # training and test data:
    y_lat_E <- loo::E_loo(t(refmodel$ref_predfun(refmodel$fit,
                                                 excl_offs = FALSE)),
                          psis_object = psisloo,
                          log_ratios = -loglik_forPSIS)
    if (any(y_lat_E$pareto_k > 0.7)) {
      warning("In the recalculation of the latent response values, ",
              sum(y_lat_E$pareto_k > 0.7), " (of ", n, ") Pareto k-value(s) ",
              "exceeded the threshold of 0.7.")
    }
    refmodel$y <- y_lat_E$value
  }

  ## by default use all observations
  nloo <- min(nloo, n)
  if (nloo < 1) {
    stop("nloo must be at least 1")
  } else if (nloo < n) {
    warning("Subsampled LOO CV is still experimental.")
  }

  ## decide which points form the validation set based on the k-values
  ## validset <- .loo_subsample(n, nloo, pareto_k)
  loo_ref_oscale <- apply(loglik_forPSIS + lw, 2, log_sum_exp)
  validset <- .loo_subsample_pps(nloo, loo_ref_oscale)
  inds <- validset$inds

  ## initialize objects where to store the results
  solution_terms_mat <- matrix(nrow = n, ncol = nterms_max - refmodel$intercept)
  loo_sub <- replicate(nterms_max, rep(NA, n), simplify = FALSE)
  mu_sub <- replicate(
    nterms_max,
    structure(rep(NA, nrow(mu_offs)),
              nobs_orig = attr(mu_offs, "nobs_orig"),
              class = sub("augmat", "augvec", oldClass(mu_offs), fixed = TRUE)),
    simplify = FALSE
  )
  if (refmodel$family$for_latent) {
    loo_sub_oscale <- loo_sub
    # In general, we could use `mu_sub_oscale <- mu_sub` here, but the case
    # where refmodel$family$latent_ilink() returns a 3-dimensional array (S x N
    # x C) needs special care.
    if (!is.null(refmodel$family$cats)) {
      mu_sub_oscale <- replicate(
        nterms_max,
        structure(rep(NA, n * length(refmodel$family$cats)),
                  nobs_orig = n,
                  class = "augvec"),
        simplify = FALSE
      )
    } else {
      mu_sub_oscale <- mu_sub
    }
  }

  if (verbose) {
    if (validate_search) {
      msg <- paste("Repeating", method, "search for", nloo, "LOO folds...")
    } else {
      msg <- paste("Computing LOO for", nterms_max, "models...")
    }
  }

  if (!validate_search) {
    # Case `validate_search = FALSE` ------------------------------------------

    if (verbose) {
      print(paste("Performing the selection using all the data.."))
    }
    ## perform selection only once using all the data (not separately for each
    ## fold), and perform the projection then for each submodel size
    search_path <- select(
      method = method, p_sel = p_sel, refmodel = refmodel,
      nterms_max = nterms_max, penalty = penalty, verbose = FALSE, opt = opt,
      search_terms = search_terms, ...
    )

    ## project onto the selected models and compute the prediction accuracy for
    ## the full data
    submodels <- .get_submodels(
      search_path = search_path,
      nterms = c(0, seq_along(search_path$solution_terms)),
      p_ref = p_pred, refmodel = refmodel, regul = opt$regul,
      refit_prj = refit_prj, ...
    )

    if (verbose) {
      print(msg)
      pb <- utils::txtProgressBar(
        min = 0, max = nterms_max, style = 3,
        initial = 0
      )
    }

    ## compute approximate LOO with PSIS weights
    if (refit_prj) {
      refdist_eval <- p_pred
    } else {
      refdist_eval <- p_sel
    }
    if (refmodel$family$for_latent) {
      refdist_eval_mu_offs_oscale <- refmodel$family$latent_ilink(
        t(refdist_eval$mu_offs), cl_ref = refdist_eval$cl,
        wdraws_ref = refdist_eval$wsample_orig
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
        wdraws_ref = refdist_eval$wsample_orig
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
    lw_sub <- lw_sub + log(refdist_eval$weights)
    # This re-weighting requires a re-normalization (as.array() is applied to
    # have stricter consistency checks, see `?sweep`):
    lw_sub <- sweep(lw_sub, 2, as.array(apply(lw_sub, 2, log_sum_exp)))
    for (k in seq_along(submodels)) {
      mu_k <- refmodel$family$mu_fun(submodels[[k]]$submodl,
                                     obs = inds,
                                     offset = refmodel$offset[inds])
      log_lik_sub <- t(refmodel$family$ll_fun(
        mu_k, submodels[[k]]$dis, refmodel$y[inds], refmodel$wobs[inds]
      ))
      loo_sub[[k]][inds] <- apply(log_lik_sub + lw_sub, 2, log_sum_exp)
      if (refmodel$family$for_latent) {
        mu_k_oscale <- refmodel$family$latent_ilink(
          t(mu_k), cl_ref = refdist_eval$cl,
          wdraws_ref = refdist_eval$wsample_orig
        )
        log_lik_sub_oscale <- refmodel$family$latent_ll_oscale(
          mu_k_oscale, y_oscale = refmodel$y_oscale[inds],
          wobs = refmodel$wobs[inds], cl_ref = refdist_eval$cl,
          wdraws_ref = refdist_eval$wsample_orig
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

      if (verbose) {
        utils::setTxtProgressBar(pb, k)
      }
    }

    for (i in seq_len(n)) {
      solution_terms_mat[
        i, seq_along(search_path$solution_terms)
      ] <- search_path$solution_terms
    }
    sel <- nlist(search_path, ce = sapply(submodels, "[[", "ce"),
                 solution_terms = search_path$solution_terms,
                 clust_used_search = p_sel$clust_used,
                 clust_used_eval = refdist_eval$clust_used,
                 nprjdraws_search = NCOL(p_sel$mu),
                 nprjdraws_eval = NCOL(refdist_eval$mu))
  } else {
    # Case `validate_search = TRUE` -------------------------------------------

    # For checking that the number of solution terms is the same across all CV
    # folds:
    prv_len_soltrms <- NULL

    if (verbose) {
      print(msg)
      pb <- utils::txtProgressBar(min = 0, max = nloo, style = 3, initial = 0)
    }

    for (run_index in seq_along(inds)) {
      ## observation index
      i <- inds[run_index]

      ## reweight the clusters/samples according to the psis-loo weights
      p_sel <- .get_p_clust(
        family = refmodel$family, eta = eta, mu = mu, mu_offs = mu_offs,
        dis = dis, wsample = exp(lw[, i]), cl = cl_sel
      )
      p_pred <- .get_p_clust(
        family = refmodel$family, eta = eta, mu = mu, mu_offs = mu_offs,
        dis = dis, wsample = exp(lw[, i]), cl = cl_pred
      )

      ## perform selection with the reweighted clusters/samples
      search_path <- select(
        method = method, p_sel = p_sel, refmodel = refmodel,
        nterms_max = nterms_max, penalty = penalty, verbose = FALSE, opt = opt,
        search_terms = search_terms, ...
      )

      ## project onto the selected models and compute the prediction accuracy
      ## for the left-out point
      submodels <- .get_submodels(
        search_path = search_path,
        nterms = c(0, seq_along(search_path$solution_terms)),
        p_ref = p_pred, refmodel = refmodel, regul = opt$regul,
        refit_prj = refit_prj, ...
      )
      summaries_sub <- .get_sub_summaries(submodels = submodels,
                                          refmodel = refmodel,
                                          test_points = i)
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

      if (is.null(prv_len_soltrms)) {
        prv_len_soltrms <- length(search_path$solution_terms)
      } else {
        stopifnot(identical(length(search_path$solution_terms),
                            prv_len_soltrms))
      }
      solution_terms_mat[
        i, seq_along(search_path$solution_terms)
      ] <- search_path$solution_terms

      if (verbose) {
        utils::setTxtProgressBar(pb, run_index)
      }
    }
  }

  # Post-processing ---------------------------------------------------------

  if (verbose) {
    ## close the progress bar object
    close(pb)
  }

  ## put all the results together in the form required by cv_varsel
  summ_sub <- lapply(seq_len(nterms_max), function(k) {
    summ_k <- list(lppd = loo_sub[[k]], mu = mu_sub[[k]], wcv = validset$wcv)
    if (refmodel$family$for_latent) {
      summ_k$oscale <- list(lppd = loo_sub_oscale[[k]], mu = mu_sub_oscale[[k]],
                            wcv = validset$wcv)
    }
    return(summ_k)
  })
  mu_ref <- do.call(c, lapply(seq_len(nrow(mu_offs)), function(i) {
    # For the augmented-data projection, `mu_offs` is an augmented-rows matrix
    # whereas the columns of `lw` refer to the original (non-augmented)
    # observations. Since `i` refers to the rows of `mu_offs`, the index for
    # `lw` needs to be adapted:
    i_nonaug <- i %% n
    if (i_nonaug == 0) {
      i_nonaug <- n
    }
    mu_offs[i, ] %*% exp(lw[, i_nonaug])
  }))
  mu_ref <- structure(
    mu_ref,
    nobs_orig = attr(mu_offs, "nobs_orig"),
    class = sub("augmat", "augvec", oldClass(mu_offs), fixed = TRUE)
  )
  if (refmodel$family$for_latent) {
    loglik_lat <- t(refmodel$family$ll_fun(
      mu_offs, dis, refmodel$y, refmodel$wobs
    ))
    lppd_ref <- apply(loglik_lat + lw, 2, log_sum_exp)
  } else {
    lppd_ref <- loo_ref_oscale
  }
  summ_ref <- list(lppd = lppd_ref, mu = mu_ref)
  if (refmodel$family$for_latent) {
    mu_ref_oscale <- do.call(c, lapply(seq_len(n_aug), function(i) {
      i_nonaug <- i %% n
      if (i_nonaug == 0) {
        i_nonaug <- n
      }
      if (inherits(mu_offs_oscale, "augmat")) {
        return(mu_offs_oscale[i, ] %*% exp(lw[, i_nonaug]))
      } else {
        # In principle, we could use the same code for averaging across the
        # draws as above in the `"augmat"` case. However, that would require
        # `mu_offs_oscale <- t(mu_offs_oscale)` beforehand, so the following
        # should be more efficient:
        return(exp(lw[, i_nonaug]) %*% mu_offs_oscale[, i])
      }
    }))
    mu_ref_oscale <- structure(
      mu_ref_oscale,
      nobs_orig = attr(mu_offs_oscale, "nobs_orig"),
      class = sub("augmat", "augvec", oldClass(mu_offs_oscale), fixed = TRUE)
    )
    summ_ref$oscale <- list(lppd = loo_ref_oscale, mu = mu_ref_oscale)
  }
  summaries <- list(sub = summ_sub, ref = summ_ref)

  d_test <- list(type = "LOO", data = NULL, offset = refmodel$offset,
                 weights = refmodel$wobs, y = refmodel$y,
                 y_oscale = refmodel$y_oscale)

  solution_terms_mat <- solution_terms_mat[
    , seq_along(search_path$solution_terms), drop = FALSE
  ]
  out_list <- nlist(solution_terms_cv = solution_terms_mat, summaries, d_test)
  if (!validate_search) {
    out_list <- c(out_list, nlist(sel))
  }
  return(out_list)
}

kfold_varsel <- function(refmodel, method, nterms_max, ndraws,
                         nclusters, ndraws_pred, nclusters_pred,
                         refit_prj, penalty, verbose, opt, K,
                         search_terms = NULL, ...) {
  # Fetch the K reference model fits (or fit them now if not already done) and
  # create objects of class `refmodel` from them (and also store the `omitted`
  # indices):
  list_cv <- .get_kfold(refmodel, K, verbose)
  K <- length(list_cv)

  if (refmodel$family$for_latent) {
    # Need to set the latent response values in `refmodel$y` to `NA`s because
    # `refmodel$y` resulted from applying `colMeans(posterior_linpred())` to the
    # original (full-data) reference model fit, so using the `fold$omitted`
    # subset of `refmodel$y` as (latent) response values in fold k of K would
    # induce a dependency between training and test data:
    refmodel$y <- rep(NA, length(refmodel$y))
  }

  # Extend `list_cv` to also contain `y`, `weights`, and `offset`:
  extend_list_cv <- function(fold) {
    d_test <- list(
      y = refmodel$y[fold$omitted],
      y_oscale = refmodel$y_oscale[fold$omitted],
      weights = refmodel$wobs[fold$omitted],
      offset = refmodel$offset[fold$omitted],
      omitted = fold$omitted
    )
    return(nlist(refmodel = fold$refmodel, d_test))
  }
  list_cv <- mapply(extend_list_cv, list_cv, SIMPLIFY = FALSE)

  # Perform the search for each fold:
  if (verbose) {
    print("Performing selection for each fold..")
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial = 0)
  }
  search_path_cv <- lapply(seq_along(list_cv), function(fold_index) {
    fold <- list_cv[[fold_index]]
    p_sel <- .get_refdist(fold$refmodel, ndraws, nclusters)
    out <- select(
      method = method, p_sel = p_sel, refmodel = fold$refmodel,
      nterms_max = nterms_max, penalty = penalty, verbose = FALSE, opt = opt,
      search_terms = search_terms, ...
    )
    if (verbose) {
      utils::setTxtProgressBar(pb, fold_index)
    }
    return(out)
  })
  solution_terms_cv <- do.call(rbind, lapply(search_path_cv, function(e) {
    e$solution_terms
  }))
  if (verbose) {
    close(pb)
  }

  # Re-project along the solution path (or fetch the projections from the search
  # results) for each fold:
  if (verbose && refit_prj) {
    print("Computing projections..")
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial = 0)
  }
  get_submodels_cv <- function(search_path, fold_index) {
    fold <- list_cv[[fold_index]]
    p_pred <- .get_refdist(fold$refmodel, ndraws_pred, nclusters_pred)
    submodels <- .get_submodels(
      search_path = search_path,
      nterms = c(0, seq_along(search_path$solution_terms)),
      p_ref = p_pred, refmodel = fold$refmodel, regul = opt$regul,
      refit_prj = refit_prj, ...
    )
    if (verbose && refit_prj) {
      utils::setTxtProgressBar(pb, fold_index)
    }
    return(submodels)
  }
  submodels_cv <- mapply(get_submodels_cv, search_path_cv, seq_along(list_cv),
                         SIMPLIFY = FALSE)
  if (verbose && refit_prj) {
    close(pb)
  }

  # Perform the evaluation of the submodels for each fold (and make sure to
  # combine the results from the K folds into a single results list):
  get_summaries_submodels_cv <- function(submodels, fold) {
    .get_sub_summaries(submodels = submodels,
                       refmodel = refmodel,
                       test_points = fold$d_test$omitted)
  }
  sub_cv_summaries <- mapply(get_summaries_submodels_cv, submodels_cv, list_cv)
  if (is.null(dim(sub_cv_summaries))) {
    summ_dim <- dim(solution_terms_cv)
    summ_dim[2] <- summ_dim[2] + 1L # +1 is for the empty model
    dim(sub_cv_summaries) <- rev(summ_dim)
  }
  sub <- apply(sub_cv_summaries, 1, rbind2list)
  idxs_sorted_by_fold <- unlist(lapply(list_cv, function(fold) {
    fold$d_test$omitted
  }))
  idxs_sorted_by_fold_aug <- idxs_sorted_by_fold
  if (!is.null(refmodel$family$cats)) {
    idxs_sorted_by_fold_aug <- idxs_sorted_by_fold_aug + rep(
      (seq_along(refmodel$family$cats) - 1L) * length(refmodel$y),
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

  # Perform the evaluation of the reference model for each fold:
  ref <- rbind2list(lapply(list_cv, function(fold) {
    eta_test <- fold$refmodel$ref_predfun(
      fold$refmodel$fit,
      newdata = refmodel$fetch_data(obs = fold$d_test$omitted),
      excl_offs = FALSE
    )
    mu_test <- fold$refmodel$family$linkinv(eta_test)
    .weighted_summary_means(
      y_test = fold$d_test, family = fold$refmodel$family,
      wsample = fold$refmodel$wsample, mu = mu_test,
      dis = fold$refmodel$dis, cl_ref = seq_along(fold$refmodel$wsample)
    )
  }))
  ref$mu <- ref$mu[order(idxs_sorted_by_fold_flx)]
  ref$lppd <- ref$lppd[order(idxs_sorted_by_fold)]
  if (!is.null(ref$oscale)) {
    ref$oscale$mu <- ref$oscale$mu[order(idxs_sorted_by_fold_aug)]
    ref$oscale$lppd <- ref$oscale$lppd[order(idxs_sorted_by_fold)]
  }

  # Combine the K separate test "datasets" (rather "information objects") into a
  # single list:
  d_cv <- rbind2list(lapply(list_cv, function(fold) {
    list(offset = fold$d_test$offset,
         weights = fold$d_test$weights,
         y = fold$d_test$y,
         y_oscale = fold$d_test$y_oscale)
  }))
  d_cv <- as.list(
    as.data.frame(d_cv)[order(idxs_sorted_by_fold), , drop = FALSE]
  )

  return(nlist(solution_terms_cv,
               summaries = nlist(sub, ref),
               d_test = c(list(type = "kfold", data = NULL), d_cv)))
}

# Re-fit the reference model K times (once for each fold; `cvfun` case) or fetch
# the K reference model fits if already computed (`cvfits` case). This function
# will return a list of length K, where each element is a list with elements
# `refmodel` (output of init_refmodel()) and `omitted` (vector of indices of
# those observations which were left out for the corresponding fold).
.get_kfold <- function(refmodel, K, verbose, approximate = FALSE) {
  if (is.null(refmodel$cvfits)) {
    if (!is.null(refmodel$cvfun)) {
      # cv-function provided so perform the cross-validation now. In case
      # refmodel is datafit, cvfun will return an empty list and this will lead
      # to normal cross-validation for the submodels although we don't have an
      # actual reference model
      if (verbose && !inherits(refmodel, "datafit")) {
        print("Performing cross-validation for the reference model..")
      }
      nobs <- NROW(refmodel$y)
      folds <- cvfolds(nobs, K = K)
      cvfits <- refmodel$cvfun(folds)
    } else {
      ## genuine probabilistic model but no K-fold fits nor cvfun provided,
      ## this only works for approximate kfold computation
      if (approximate) {
        nobs <- NROW(refmodel$y)
        folds <- cvfolds(nobs, K = K)
        cvfits <- lapply(seq_len(K), function(k) {
          refmodel$fit
        })
      } else {
        stop("For a reference model which is not of class `datafit`, either ",
             "`cvfits` or `cvfun` needs to be provided for K-fold CV (see ",
             "`?init_refmodel`).")
      }
    }
  } else {
    cvfits <- refmodel$cvfits
    K <- attr(cvfits, "K")
    folds <- attr(cvfits, "folds")
    cvfits <- cvfits$fits
  }
  cvfits <- lapply(seq_len(K), function(k) {
    cvfit <- cvfits[[k]]
    # Add the omitted observation indices for this fold:
    cvfit$omitted <- which(folds == k)
    # Add the fold index:
    cvfit$projpred_k <- k
    return(cvfit)
  })
  return(lapply(cvfits, .init_kfold_refmodel, refmodel = refmodel))
}

.init_kfold_refmodel <- function(cvfit, refmodel) {
  return(list(refmodel = refmodel$cvrefbuilder(cvfit),
              omitted = cvfit$omitted))
}

# ## decide which points to go through in the validation (i.e., which points
# ## belong to the semi random subsample of validation points)
# .loo_subsample <- function(n, nloo, pareto_k) {
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
.loo_subsample_pps <- function(nloo, lppd) {
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
