#' Variable selection with cross-validation
#'
#' Perform the projection predictive variable selection for (G)LMs, (G)LMMs,
#' (G)AMs, and (G)AMMs. This variable selection consists of a *search* part and
#' an *evaluation* part. The search part determines the solution path, i.e., the
#' best submodel for each number of predictor terms (model size). The evaluation
#' part determines the predictive performance of the submodels along the
#' solution path. In contrast to [varsel()], [cv_varsel()] performs a
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
#' @param nloo Only relevant if `cv_method == "LOO"`. Number of subsampled LOO
#'   CV folds, i.e., number of observations used for the LOO CV (anything
#'   between 1 and the original number of observations). Smaller values lead to
#'   faster computation but higher uncertainty in the evaluation part. If
#'   `NULL`, all observations are used, but for faster experimentation, one can
#'   set this to a smaller value.
#' @param K Only relevant if `cv_method == "kfold"`. Number of folds in the
#'   \eqn{K}-fold CV.
#' @param validate_search Only relevant if `cv_method == "LOO"`. A single
#'   logical value indicating whether to cross-validate also the search part,
#'   i.e., whether to run the search separately for each CV fold (`TRUE`) or not
#'   (`FALSE`). We strongly do not recommend setting this to `FALSE`, because
#'   this is known to bias the predictive performance estimates of the selected
#'   submodels. However, setting this to `FALSE` can sometimes be useful because
#'   comparing the results to the case where this argument is `TRUE` gives an
#'   idea of how strongly the variable selection is (over-)fitted to the data
#'   (the difference corresponds to the search degrees of freedom or the
#'   effective number of parameters introduced by the search).
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. If `NULL`, no seed is set and
#'   therefore, the results are not reproducible. See [set.seed()] for details.
#'   Here, this seed is used for clustering the reference model's posterior
#'   draws (if `!is.null(nclusters)`), for subsampling LOO CV folds (if `nloo`
#'   is smaller than the number of observations), and for sampling the folds in
#'   K-fold CV.
#'
#' @inherit varsel details return
#'
#' @note The case `cv_method == "LOO" && !validate_search` constitutes an
#'   exception where the search part is not cross-validated. In that case, the
#'   evaluation part is based on a PSIS-LOO CV.
#'
#' @references
#'
#' Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics and
#' Computing*, **27**(5), 1413-1432. DOI:
#' [10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4).
#'
#' Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2021). Pareto
#' smoothed importance sampling. *arXiv:1507.02646*. URL:
#' <https://arxiv.org/abs/1507.02646>.
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
  ndraws = 20,
  nclusters = NULL,
  ndraws_pred = 400,
  nclusters_pred = NULL,
  cv_search = !inherits(object, "datafit"),
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
  seed = NULL,
  search_terms = NULL,
  ...
) {
  refmodel <- object
  ## resolve the arguments similar to varsel
  args <- parse_args_varsel(
    refmodel = refmodel, method = method, cv_search = cv_search,
    intercept = NULL, nterms_max = nterms_max, nclusters = nclusters,
    ndraws = ndraws, nclusters_pred = nclusters_pred,
    ndraws_pred = ndraws_pred, search_terms = search_terms
  )
  method <- args$method
  cv_search <- args$cv_search
  intercept <- args$intercept
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  ndraws <- args$ndraws
  nclusters_pred <- args$nclusters_pred
  ndraws_pred <- args$ndraws_pred
  search_terms <- args$search_terms
  has_group_features <- formula_contains_group_terms(refmodel$formula)
  has_additive_features <- formula_contains_additive_terms(refmodel$formula)

  if (method == "l1" && (has_group_features || has_additive_features)) {
    stop("L1 search is only supported for GLMs.")
  }

  ## arguments specific to this function
  args <- parse_args_cv_varsel(
    refmodel, cv_method, K, nclusters,
    nclusters_pred
  )
  cv_method <- args$cv_method
  K <- args$K
  nclusters <- args$nclusters
  nclusters_pred <- args$nclusters_pred

  ## search options
  opt <- nlist(lambda_min_ratio, nlambda, thresh, regul)

  if (cv_method == "LOO") {
    sel_cv <- loo_varsel(
      refmodel = refmodel, method = method, nterms_max = nterms_max,
      ndraws = ndraws, nclusters = nclusters,
      ndraws_pred = ndraws_pred,
      nclusters_pred = nclusters_pred,
      cv_search = cv_search, intercept = intercept, penalty = penalty,
      verbose = verbose, opt = opt, nloo = nloo,
      validate_search = validate_search, seed = seed,
      search_terms = search_terms
    )
  } else if (cv_method == "kfold") {
    sel_cv <- kfold_varsel(
      refmodel = refmodel, method = method, nterms_max = nterms_max,
      ndraws = ndraws, nclusters = nclusters,
      ndraws_pred = ndraws_pred,
      nclusters_pred = nclusters_pred,
      cv_search = cv_search, intercept = intercept,
      penalty = penalty, verbose = verbose, opt = opt, K = K,
      seed = seed, search_terms = search_terms
    )
  } else {
    stop(sprintf("Unknown `cv_method`: %s.", method))
  }

  if (validate_search || cv_method == "kfold") {
    ## run the selection using the full dataset
    if (verbose) {
      print(paste("Performing the selection using all the data.."))
    }
    sel <- varsel(refmodel,
                  method = method, ndraws = ndraws, nclusters = nclusters,
                  ndraws_pred = ndraws_pred, nclusters_pred = nclusters_pred,
                  cv_search = cv_search, nterms_max = nterms_max - 1,
                  intercept = intercept, penalty = penalty, verbose = verbose,
                  lambda_min_ratio = lambda_min_ratio, nlambda = nlambda,
                  regul = regul, search_terms = search_terms, seed = seed)
  } else if (cv_method == "LOO") {
    sel <- sel_cv$sel
  }

  # Find out how many CV folds select the same variables as the selection with
  # all the data (assuming all CV folds have equal weight):
  candidate_terms <- split_formula(refmodel$formula,
                                   data = refmodel$fetch_data(),
                                   add_main_effects = FALSE)
  candidate_terms <- unlist(utils::tail(candidate_terms, -1))
  solution_terms_cv_ch <- sapply(
    seq_len(NROW(sel_cv$solution_terms_cv)),
    function(i) {
      if (!is.character(sel_cv$solution_terms_cv[i, ])) {
        return(candidate_terms[sel_cv$solution_terms_cv[i, ]])
      } else {
        return(sel_cv$solution_terms_cv[i, ])
      }
    }
  )
  sel_solution_terms <- unlist(sel$solution_terms)
  if (!is.matrix(solution_terms_cv_ch)) {
    stop("Unexpected `solution_terms_cv_ch`. Please notify the package ",
         "maintainer.")
  }
  if (!identical(nrow(solution_terms_cv_ch), length(sel_solution_terms))) {
    stop("Unexpected number of rows in `solution_terms_cv_ch`. Please notify ",
         "the package maintainer.")
  }
  pct_solution_terms_cv <- cbind(
    size = seq_len(nrow(solution_terms_cv_ch)),
    sapply(sel_solution_terms, function(var_nm) {
      rowMeans(solution_terms_cv_ch == var_nm, na.rm = TRUE)
    })
  )

  ## create the object to be returned
  vs <- nlist(refmodel,
              search_path = sel$search_path,
              d_test = sel_cv$d_test,
              summaries = sel_cv$summaries,
              family = refmodel$family,
              kl = sel$kl,
              solution_terms = sel$solution_terms,
              pct_solution_terms_cv,
              nterms_all = count_terms_in_formula(refmodel$formula),
              nterms_max,
              method,
              cv_method,
              validate_search,
              nclusters,
              nclusters_pred,
              ndraws,
              ndraws_pred)
  class(vs) <- "vsel"
  vs$suggested_size <- suggest_size(vs, warnings = FALSE)
  summary <- summary(vs)
  vs$summary <- summary$selection
  if (verbose) {
    print("Done.")
  }
  return(vs)
}

#
# Auxiliary function for parsing the input arguments for specific cv_varsel.
# This is similar in spirit to parse_args_varsel, that is, to avoid the main
# function to become too long and complicated to maintain.
#
# @param refmodel Reference model as extracted by get_refmodel
# @param cv_method The cross-validation method, either 'LOO' or 'kfold'.
#   Default is 'LOO'.
# @param K Number of folds in the K-fold cross validation. Default is 5 for
#   genuine reference models and 10 for datafits (that is, for penalized
#   maximum likelihood estimation).
parse_args_cv_varsel <- function(refmodel, cv_method, K,
                                 nclusters, nclusters_pred) {
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
  }

  return(nlist(cv_method, K, nclusters, nclusters_pred))
}

loo_varsel <- function(refmodel, method, nterms_max, ndraws,
                       nclusters, ndraws_pred, nclusters_pred, cv_search,
                       intercept, penalty, verbose, opt, nloo = NULL,
                       validate_search = TRUE, seed = NULL,
                       search_terms = NULL) {
  ##
  ## Perform the validation of the searching process using LOO. validate_search
  ## indicates whether the selection is performed separately for each fold (for
  ## each data point)
  ##

  family <- refmodel$family
  mu <- refmodel$mu
  dis <- refmodel$dis
  ## the clustering/subsampling used for selection
  p_sel <- .get_refdist(refmodel,
                        ndraws = ndraws,
                        nclusters = nclusters,
                        seed = seed)
  cl_sel <- p_sel$cl # clustering information

  ## the clustering/subsampling used for prediction
  p_pred <- .get_refdist(refmodel,
                         ndraws = ndraws_pred,
                         nclusters = nclusters_pred,
                         seed = seed)
  cl_pred <- p_pred$cl

  ## fetch the log-likelihood for the reference model to obtain the LOO weights
  if (is.null(refmodel$loglik)) {
    ## case where log-likelihood not available, i.e., the reference model is not
    ## a genuine model => cannot compute LOO
    stop("LOO can be performed only if the reference model is a genuine ",
         "probabilistic model for which the log-likelihood can be evaluated.")
  } else {
    ## log-likelihood available
    loglik <- refmodel$loglik
  }
  n <- ncol(loglik)
  ## TODO: should take r_eff:s into account
  psisloo <- loo::psis(-loglik, cores = 1, r_eff = rep(1, n))
  lw <- weights(psisloo)
  # pareto_k <- loo::pareto_k_values(psisloo)
  ## by default use all observations
  nloo <- min(nloo, n)

  if (nloo < 1) {
    stop("nloo must be at least 1")
  }

  ## compute loo summaries for the reference model
  loo_ref <- apply(loglik + lw, 2, log_sum_exp)
  mu_ref <- rep(0, n)
  for (i in seq_len(n)) {
    mu_ref[i] <- mu[i, ] %*% exp(lw[, i])
  }

  ## decide which points form the validation set based on the k-values
  ## validset <- .loo_subsample(n, nloo, pareto_k, seed)
  validset <- .loo_subsample_pps(nloo, loo_ref, seed)
  inds <- validset$inds

  ## initialize matrices where to store the results
  solution_terms_mat <- matrix(nrow = n, ncol = nterms_max - 1)
  loo_sub <- matrix(nrow = n, ncol = nterms_max)
  mu_sub <- matrix(nrow = n, ncol = nterms_max)

  if (verbose) {
    if (validate_search) {
      msg <- paste("Repeating", method, "search for", nloo, "LOO folds...")
    } else {
      msg <- paste("Computing LOO for", nterms_max, "models...")
    }
  }

  if (!validate_search) {
    if (verbose) {
      print(paste("Performing the selection using all the data.."))
    }
    ## perform selection only once using all the data (not separately for each
    ## fold), and perform the projection then for each submodel size
    search_path <- select(
      method = method, p_sel = p_sel, refmodel = refmodel, family = family,
      intercept = intercept, nterms_max = nterms_max, penalty = penalty,
      verbose = FALSE, opt = opt, search_terms = search_terms
    )
    solution_terms <- search_path$solution_terms

    ## project onto the selected models and compute the prediction accuracy for
    ## the full data
    submodels <- .get_submodels(
      search_path = search_path, nterms = c(0, seq_along(solution_terms)),
      family = family, p_ref = p_pred, refmodel = refmodel,
      intercept = intercept, regul = opt$regul, cv_search = cv_search
    )

    if (verbose) {
      print(msg)
      pb <- utils::txtProgressBar(
        min = 0, max = nterms_max, style = 3,
        initial = 0
      )
    }

    ## compute approximate LOO with PSIS weights
    y <- matrix(refmodel$y, nrow = n)
    for (k in seq_along(submodels)) {
      mu_k <- family$mu_fun(submodels[[k]]$sub_fit,
                            obs = inds,
                            offset = refmodel$offset)
      log_lik_sub <- t(family$ll_fun(
        mu_k, submodels[[k]]$dis,
        y[inds], refmodel$wobs[inds]
      ))
      sub_psisloo <- suppressWarnings(
        loo::psis(-log_lik_sub,
                  cores = 1,
                  r_eff = rep(1, ncol(log_lik_sub)))
      )
      lw_sub <- suppressWarnings(loo::weights.importance_sampling(sub_psisloo))
      loo_sub[inds, k] <- apply(
        log_lik_sub[,] + lw_sub[,], 2,
        log_sum_exp
      )
      for (i in seq_along(inds)) {
        mu_sub[inds[i], k] <- mu_k[i, ] %*% exp(lw_sub[, i])
      }

      if (verbose) {
        utils::setTxtProgressBar(pb, k)
      }
    }

    candidate_terms <- split_formula(refmodel$formula,
                                     data = refmodel$fetch_data(),
                                     add_main_effects = FALSE)
    ## with `match` we get the indices of the variables as they enter the
    ## solution path in solution_terms
    solution <- match(solution_terms, candidate_terms[-1])
    for (i in seq_len(n)) {
      solution_terms_mat[i, seq_along(solution)] <- solution
    }
    sel <- nlist(search_path, kl = sapply(submodels, function(x) x$kl),
                 solution_terms)
  } else {
    if (verbose) {
      print(msg)
      pb <- utils::txtProgressBar(min = 0, max = nloo, style = 3, initial = 0)
    }

    for (run_index in seq_along(inds)) {
      ## observation index
      i <- inds[run_index]

      ## reweight the clusters/samples according to the psis-loo weights
      p_sel <- .get_p_clust(
        family = family, mu = mu, dis = dis, wsample = exp(lw[, i]),
        cl = cl_sel
      )
      p_pred <- .get_p_clust(
        family = family, mu = mu, dis = dis, wsample = exp(lw[, i]),
        cl = cl_pred
      )

      ## perform selection with the reweighted clusters/samples
      search_path <- select(
        method = method, p_sel = p_sel, refmodel = refmodel,
        family = family, intercept = intercept, nterms_max = nterms_max,
        penalty = penalty, verbose = FALSE, opt = opt,
        search_terms = search_terms
      )
      solution_terms <- search_path$solution_terms

      ## project onto the selected models and compute the prediction accuracy
      ## for the left-out point
      submodels <- .get_submodels(
        search_path = search_path, nterms = c(0, seq_along(solution_terms)),
        family = family, p_ref = p_pred, refmodel = refmodel,
        intercept = intercept, regul = opt$regul, cv_search = cv_search
      )
      summaries_sub <- .get_sub_summaries(
        submodels = submodels, test_points = c(i), refmodel = refmodel,
        family = family
      )
      for (k in seq_along(summaries_sub)) {
        loo_sub[i, k] <- summaries_sub[[k]]$lppd
        mu_sub[i, k] <- summaries_sub[[k]]$mu
      }

      candidate_terms <- split_formula(refmodel$formula,
                                       data = refmodel$fetch_data(),
                                       add_main_effects = FALSE)
      ## with `match` we get the indices of the variables as they enter the
      ## solution path in solution_terms
      solution <- match(solution_terms, candidate_terms[-1])
      solution_terms_mat[i, seq_along(solution)] <- solution

      if (verbose) {
        utils::setTxtProgressBar(pb, run_index)
      }
    }
  }

  if (verbose) {
    ## close the progress bar object
    close(pb)
  }

  ## put all the results together in the form required by cv_varsel
  summ_sub <- lapply(seq_len(nterms_max), function(k) {
    list(lppd = loo_sub[, k], mu = mu_sub[, k], w = validset$w)
  })
  summ_ref <- list(lppd = loo_ref, mu = mu_ref)
  summaries <- list(sub = summ_sub, ref = summ_ref)

  d_test <- list(
    y = refmodel$y, type = "LOO",
    test_points = seq_along(refmodel$y),
    weights = refmodel$wobs,
    data = NULL, offset = refmodel$offset
  )

  if (!validate_search) {
    return(nlist(
      solution_terms_cv = solution_terms_mat,
      summaries, d_test, sel
    ))
  } else {
    return(nlist(solution_terms_cv = solution_terms_mat, summaries, d_test))
  }
}

kfold_varsel <- function(refmodel, method, nterms_max, ndraws,
                         nclusters, ndraws_pred, nclusters_pred,
                         cv_search, intercept, penalty, verbose, opt,
                         K, seed = NULL, search_terms = NULL) {
  ## fetch the k_fold list (or compute it now if not already computed)
  k_fold <- .get_kfold(refmodel, K, verbose, seed)

  ## check that k_fold has the correct form
  ## .validate_kfold(refmodel, k_fold, refmodel$nobs)

  K <- length(k_fold)
  family <- refmodel$family

  ## extract variables from each fit-object (samples, x, y, etc.)
  ## to a list of size K
  refmodels_cv <- lapply(k_fold, function(fold) fold$refmodel)

  # List of size K with test data for each fold
  d_test_cv <- lapply(k_fold, function(fold) {
    list(
      newdata = refmodel$fetch_data(obs = fold$omitted),
      y = refmodel$y[fold$omitted],
      weights = refmodel$wobs[fold$omitted],
      offset = refmodel$offset[fold$omitted],
      omitted = fold$omitted
    )
  })

  ## List of K elements, each containing d_test, p_pred, etc. corresponding
  ## to each fold.
  make_list_cv <- function(refmodel, d_test, msg) {
    if (!is.null(nclusters_pred) || !is.null(refmodel$nclusters_pred)) {
      nclusters_pred <- min(
        refmodel$nclusters_pred,
        nclusters_pred
      )
    }
    p_sel <- .get_refdist(refmodel, ndraws, nclusters, seed = seed)
    p_pred <- .get_refdist(refmodel, ndraws_pred, nclusters_pred, seed = seed)
    pred <- refmodel$ref_predfun(refmodel$fit, newdata = d_test$newdata) +
      d_test$offset
    pred <- matrix(
      as.numeric(pred), nrow = NROW(pred), ncol = NCOL(pred)
    )
    mu_test <- family$linkinv(pred)
    nlist(refmodel, p_sel, p_pred, mu_test,
          dis = refmodel$dis, w_test = refmodel$wsample, d_test, msg)
  }

  msgs <- paste0(method, " search for fold ", 1:K, "/", K, ".")
  list_cv <- mapply(make_list_cv, refmodels_cv, d_test_cv, msgs,
                    SIMPLIFY = FALSE)

  ## Perform the selection for each of the K folds
  if (verbose) {
    print("Performing selection for each fold..")
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial = 0)
  }
  search_path_cv <- lapply(seq_along(list_cv), function(fold_index) {
    fold <- list_cv[[fold_index]]
    family <- fold$refmodel$family
    out <- select(
      method = method, p_sel = fold$p_sel, refmodel = fold$refmodel,
      family = family, intercept = intercept, nterms_max = nterms_max,
      penalty = penalty, verbose = FALSE, opt = opt,
      search_terms = search_terms
    )
    if (verbose) {
      utils::setTxtProgressBar(pb, fold_index)
    }
    out
  })

  solution_terms_cv <- t(sapply(search_path_cv, function(e) e$solution_terms))
  if (verbose) {
    close(pb)
  }

  ## Construct submodel projections for each fold
  if (verbose && cv_search) {
    print("Computing projections..")
    pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial = 0)
  }

  get_submodels_cv <- function(search_path, fold_index) {
    fold <- list_cv[[fold_index]]
    family <- fold$refmodel$family
    solution_terms <- search_path$solution_terms
    submodels <- .get_submodels(
      search_path = search_path, nterms = c(0, seq_along(solution_terms)),
      family = family, p_ref = fold$p_pred, refmodel = fold$refmodel,
      intercept = intercept, regul = opt$regul, cv_search = FALSE
    )
    if (verbose && cv_search) {
      utils::setTxtProgressBar(pb, fold_index)
    }
    return(submodels)
  }

  submodels_cv <- mapply(get_submodels_cv, search_path_cv, seq_along(list_cv),
                         SIMPLIFY = FALSE)
  if (verbose && cv_search) {
    close(pb)
  }

  ## Apply some magic to manipulate the structure of the list so that instead of
  ## list with K sub_summaries each containing n/K mu:s and lppd:s, we have only
  ## one sub_summary-list that contains with all n mu:s and lppd:s.
  get_summaries_submodel_cv <- function(submodels, fold) {
    omitted <- fold$d_test$omitted
    fold_summaries <- .get_sub_summaries(
      submodels = submodels, test_points = omitted, refmodel = refmodel,
      family = family
    )
    summ <- lapply(fold_summaries, data.frame)
    return(summ)
  }
  sub_cv_summaries <- mapply(get_summaries_submodel_cv, submodels_cv, list_cv)
  sub <- apply(sub_cv_summaries, 1, hf)
  sub <- lapply(sub, function(summ) {
    summ$w <- rep(1, length(summ$mu))
    summ$w <- summ$w / sum(summ$w)
    summ
  })

  ref <- hf(lapply(list_cv, function(fold) {
    data.frame(.weighted_summary_means(
      y_test = fold$d_test, family = family, wsample = fold$refmodel$wsample,
      mu = fold$mu_test, dis = fold$refmodel$dis
    ))
  }))

  ## Combine also the K separate test data sets into one list
  ## with n y's and weights's.
  d_cv <- hf(lapply(d_test_cv, function(fold) {
    data.frame(
      y = fold$y, weights = fold$weights,
      test_points = fold$omitted,
      offset = fold$offset
    )
  }))

  return(nlist(solution_terms_cv,
               summaries = list(sub = sub, ref = ref),
               d_test = c(d_cv, type = "kfold")))
}


.get_kfold <- function(refmodel, K, verbose, seed, approximate = FALSE) {
  ## Fetch the k_fold list or compute it now if not already computed. This
  ## function will return a list of length K, where each element is a list
  ## with fields 'refmodel' (object of type refmodel computed by init_refmodel)
  ## and index list 'test_points' that denotes which of the data points were
  ## left out for the corresponding fold.

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
      folds <- cvfolds(nobs, K = K, seed = seed)
      cvfits <- refmodel$cvfun(folds)
      cvfits <- lapply(seq_along(cvfits), function(k) {
        # add the 'omitted' indices for the cvfits
        cvfit <- cvfits[[k]]
        cvfit$omitted <- which(folds == k)
        cvfit
      })
    } else {
      ## genuine probabilistic model but no K-fold fits nor cvfun provided,
      ## this only works for approximate kfold computation
      if (approximate) {
        nobs <- NROW(refmodel$y)
        folds <- cvfolds(nobs, K = K, seed = seed)
        cvfits <- lapply(seq_len(K), function(k) {
          ## add the 'omitted' indices for the cvfits
          cvfit <- refmodel$fit
          cvfit$omitted <- which(folds == k)
          cvfit
        })
      } else {
        stop(
          "For a generic reference model, you must provide either `cvfits` or ",
          "`cvfun` for K-fold cross-validation. See function init_refmodel()."
        )
      }
    }
  } else {
    cvfits <- refmodel$cvfits
    K <- attr(cvfits, "K")
    folds <- attr(cvfits, "folds")
    cvfits <- lapply(seq_len(K), function(k) {
      cvfit <- cvfits$fits[[k]]
      cvfit$omitted <- which(folds == k)
      cvfit
    })
  }

  train <- seq_along(refmodel$y)

  k_fold <- lapply(cvfits, .init_kfold_refmodel, refmodel, train)

  return(k_fold)
}

.init_kfold_refmodel <- function(cvfit, refmodel, train) {
  fold <- setdiff(
    train,
    cvfit$omitted
  )
  default_data <- refmodel$fetch_data(obs = fold)
  ref_predfun <- function(fit, newdata = default_data) {
    refmodel$ref_predfun(fit, newdata = newdata)
  }
  proj_predfun <- function(fit, newdata = default_data) {
    refmodel$proj_predfun(fit, newdata = newdata)
  }
  extract_model_data <- function(object, newdata = default_data, ...) {
    refmodel$extract_model_data(object = object, newdata = newdata, ...)
  }
  if (!inherits(refmodel, "datafit")) {
    k_refmodel <- get_refmodel(cvfit)
  } else {
    k_refmodel <- init_refmodel(
      object = NULL, data = default_data,
      formula = refmodel$formula, family = refmodel$family,
      div_minimizer = refmodel$div_minimizer,
      proj_predfun = proj_predfun,
      extract_model_data = extract_model_data
    )
  }
  ## k_refmodel$nclusters_pred <- min(NCOL(k_refmodel$mu), 5)
  return(nlist(refmodel = k_refmodel, omitted = cvfit$omitted))
}

# .loo_subsample <- function(n, nloo, pareto_k, seed) {
#   ## decide which points to go through in the validation (i.e., which points
#   ## belong to the semi random subsample of validation points)
#
#   ## set random seed but ensure the old RNG state is restored on exit
#   if (exists(".Random.seed")) {
#     rng_state_old <- .Random.seed
#     on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
#   }
#   set.seed(seed)
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
#     ## assign the weights corresponding to this stratification (for example, the
#     ## 'bad' values are likely to be overpresented in the sample)
#     w <- rep(0, n)
#     w[inds[inds %in% bad]] <- length(bad) / sum(inds %in% bad)
#     w[inds[inds %in% ok]] <- length(ok) / sum(inds %in% ok)
#     w[inds[inds %in% good]] <- length(good) / sum(inds %in% good)
#   } else {
#     ## all points used
#     inds <- seq_len(n)
#     w <- rep(1, n)
#   }
#
#   ## ensure weights are normalized
#   w <- w / sum(w)
#
#   return(nlist(inds, w))
# }

.loo_subsample_pps <- function(nloo, lppd, seed) {
  ## decide which points to go through in the validation based on
  ## proportional-to-size subsampling as implemented in Magnusson, M., Riis
  ## Andersen, M., Jonasson, J. and Vehtari, A. (2019). Leave-One-Out
  ## Cross-Validation for Large Data. In International Conference on Machine
  ## Learning.

  ## set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)
  if (nloo > length(lppd)) {
    stop("Argument `nloo` must not be larger than the number of observations.")
  } else if (nloo == length(lppd)) {
    inds <- seq_len(nloo)
    w <- rep(1, nloo)
  } else if (nloo < length(lppd)) {
    w <- exp(lppd - max(lppd))
    inds <- sample(seq_along(lppd), size = nloo, prob = w)
  }
  w <- w / sum(w)

  return(nlist(inds, w))
}
