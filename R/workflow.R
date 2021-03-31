#' Variable and structure selection workflow
#'
#' Perform variable and structure selection in a cascade workflow. First, it
#' searches through model space to determine the ordering of model terms for
#' projection. Second, we compute the projections with approximate cross
#' validation. Third, if selection diagnostics indicate that more posterior
#' draws are required, we recommend running a full search with more draws.
#'
#' @name workflow
#'
NULL

#' @rdname workflow
#' @export
varsel_search <- function(object, ...) {
    UseMethod("varsel_search")
}

#' @rdname workflow
#' @export
varsel_search.default <- function(object, ...) {
    refmodel <- get_refmodel(object, ...)
    return(varsel_search(refmodel, ...))
}


#' @rdname workflow
#' @export
varsel_search.refmodel <- function(object, method = NULL,
                                   ndraws = NULL, nclusters = NULL,
                                   nterms_max = NULL, verbose = TRUE,
                                   lambda_min_ratio = 1e-5, nlambda = 150,
                                   thresh = 1e-6, regul = 1e-4,
                                   penalty = NULL, search_terms = NULL,
                                   ...) {
  refmodel <- object
  family <- refmodel$family

  ## fetch the default arguments or replace them by the user defined values
  args <- parse_args_varsel(
    refmodel, method, NULL, TRUE, nterms_max,
    nclusters, ndraws, NULL, NULL, search_terms
  )
  method <- args$method
  intercept <- args$intercept
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  ndraws <- args$ndraws
  search_terms <- args$search_terms
  has_group_features <- formula_contains_group_terms(refmodel$formula)

  if (method == "l1" && has_group_features) {
    stop(
      "l1 search is not supported for multilevel models",
      )
  }

  ## reference distributions for selection and prediction after selection
  p_sel <- .get_refdist(refmodel, ndraws, nclusters)

  ## perform the selection
  opt <- nlist(lambda_min_ratio, nlambda, thresh, regul)
  start <- Sys.time()
  search_path <- select(
    method = method, p_sel = p_sel, refmodel = refmodel,
    family = family, intercept = intercept, nterms_max = nterms_max,
    penalty = penalty, verbose = verbose, opt = opt,
    search_terms = search_terms
  )
  class(search_path) <- "vselsearch"
  search_path$control <- nlist(
    method,
    time = difftime(Sys.time(), start, units = "secs"),
    opt,
    nterms_max,
    ndraws,
    nclusters
  )
  search_path$refmodel <- refmodel
  return(search_path)
}

#' @export
summary.vselsearch <- function(object, ...) {
  out <- list(
    formula = object$refmodel$formula,
    fit = object$refmodel$fit,
    family = object$refmodel$family,
    nobs = NROW(object$refmodel$fetch_data()),
    method = object$control$method,
    ndraws = object$control$ndraws,
    nclusters = object$control$nclusters,
    nterms_max = object$control$nterms_max,
    time = object$control$time,
    nsubmodels = object$nsubmodels,
    solution_terms = object$solution_terms
  )
  class(out) <- "vselsearchsummary"
  return(out)
}

#' @export
#' @method print vselsearch
print.vselsearch <- function(x, digits = 1, ...) {
  summary <- summary.vselsearch(x, digits = digits, ...)
  print(summary)
  return(invisible(summary))
}

#' @export
print.vselsearchsummary <- function(x, digits = 1, ...) {
  print(x$family)
  cat("Formula: ")
  print(x$formula)
  cat(paste0("Observations: ", x$nobs, "\n"))
  nterms_max <- x$nterms_max
  cat(paste0(
      "Search method: ", x$method, ", maximum number of terms ",
      nterms_max, "\n"
  ))
  cat(paste0(
      "\nSolution search path: ", paste(x$solution_terms, collapse = ", "),
      "\n"
  ))
  cat(paste0(
    "Draws used for selection: ", x$ndraws, ", in ",
    x$nclusters, " clusters\n"
  ))
  cat(paste0("Number of submodels visited during search: ", x$nsubmodels, "\n"))
  cat(paste0("\nThe search took ", round(x$time, digits), " seconds.\n"))
  if (x$method != "l1") {
    one_proj <- round(x$time / x$nsubmodels / x$nclusters, digits)
  } else {
    one_proj <- round(x$time, digits)
  }
  cat(paste0("Projecting one draw takes roughly ", one_proj, " seconds\n"))
  approx_cv <- length(x$solution_terms)
  cat(paste0(
      "\nApproximate LOO would roughly take ndraws_pred * ", approx_cv *
      one_proj, " + ", approx_cv, " times one LOO seconds\n"
  ))
  cat(paste0(
    "Approximate KFold would roughly take K * reference fit time",
    " + ndraws_pred * ", approx_cv * one_proj, " + K * ", approx_cv,
    " times one kfold seconds\n"
  ))
  cat(paste0(
    "\nFull LOO would roughly take n * ", round(x$time, digits),
    " + nloo * ndraws_pred * ",
    approx_cv * one_proj, " seconds\n"
  ))
  cat(paste0(
    "Full KFold would roughly take K * reference fit time + K * ",
    round(x$time, digits), " + K * ndraws_pred * ",
    approx_cv * one_proj, " seconds\n"
  ))
  return(invisible(x))
}

#' @rdname workflow
#' @export
approximate_loo <- function(object, ...) {
    UseMethod("approximate_loo")
}

#' @rdname workflow
#' @export
approximate_loo.default <- function(object, ...) {
    return(approximate_loo(varsel_search(object, ...)))
}

#' @rdname workflow
#' @export
approximate_loo.refmodel <- function(object, ...) {
    return(approximate_loo(varsel_search(object, ...)))
}

#' @rdname workflow
#' @export
approximate_loo.vselsearch <- function(object,
                                       ndraws_pred = NULL,
                                       nclusters_pred = NULL,
                                       verbose = TRUE,
                                       penalty = NULL,
                                       nloo = NULL,
                                       refit_proj = TRUE,
                                       seed = NULL,
                                       ...) {
  refmodel <- object$refmodel
  family <- refmodel$family
  nterms_max <- object$control$nterms_max

  ## fetch the default arguments or replace them by the user defined values
  args <- parse_args_varsel(
    refmodel, NULL, NULL, NULL, NULL,
    NULL, NULL, ndraws_pred, nclusters_pred, NULL
  )
  nclusters_pred <- args$nclusters_pred
  ndraws_pred <- args$ndraws_pred

  mu <- refmodel$mu
  dis <- refmodel$dis
  ## the clustering/subsampling used for selection

  ## the clustering/subsampling used for prediction
  p_pred <- .get_refdist(refmodel,
    ndraws = ndraws_pred,
    nclusters = nclusters_pred
  )
  cl_pred <- p_pred$cl
  ## fetch the log-likelihood for the reference model to obtain the LOO
  ## weights
  if (is.null(refmodel$loglik)) {
    ## case where log-likelihood not available, i.e., the reference model is
    ## not a genuine model => cannot compute LOO
    stop(
      "LOO can be performed only if the reference model is a genuine ",
      "probabilistic model for which the log-likelihood can be evaluated."
    )
  } else {
    ## log-likelihood available
    loglik <- refmodel$loglik
  }
  psisloo <- loo::psis(-loglik, cores = 1, r_eff = rep(1, ncol(loglik)))
  lw <- weights(psisloo)
  pareto_k <- loo::pareto_k_values(psisloo)
  n <- length(pareto_k)
  ## by default use all observations
  nloo <- min(nloo, n)

  if (nloo < 0) {
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
  solution_terms_cv <- matrix(nrow = n, ncol = nterms_max - 1)
  loo_sub <- matrix(nrow = n, ncol = nterms_max)
  mu_sub <- matrix(nrow = n, ncol = nterms_max)

  start <- Sys.time()
  if (verbose) {
    msg <- paste("Computing LOO for", nterms_max, "models...")
  }
  solution_terms <- object$solution_terms

  ## project onto the selected models and compute the prediction accuracy for
  ## the full data
  submodels <- .get_submodels(
    search_path = object, nterms = c(0, seq_along(solution_terms)),
    family = family, p_ref = p_pred, refmodel = refmodel,
    intercept = TRUE, regul = object$control$opt$regul, cv_search = refit_proj
  )
  summaries_sub <- .get_sub_summaries(
    submodels = submodels, test_points = seq_len(n),
    refmodel = refmodel, family = family
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
  for (k in seq_along(summaries_sub)) {
    mu_k <- family$mu_fun(submodels[[k]]$sub_fit,
        obs = inds,
        offset = refmodel$offset,
        weights = 1
    )
    log_lik_sub <- t(family$ll_fun(
        mu_k, submodels[[k]]$dis,
        y[inds], refmodel$wobs
    ))
    sub_psisloo <- suppressWarnings(loo::psis(-log_lik_sub,
        cores = 1,
        r_eff = rep(1, ncol(log_lik_sub))
    ))
    lw_sub <- suppressWarnings(
        loo::weights.importance_sampling(sub_psisloo)
    )
    loo_sub[inds, k] <- apply(
        log_lik_sub[, ] + lw_sub[, ], 2,
        log_sum_exp
    )
    for (i in seq_along(inds)) {
        mu_sub[inds[i], k] <- mu_k[i, ] %*% exp(lw_sub[, i])
    }

    if (verbose) {
        utils::setTxtProgressBar(pb, k)
    }
  }
  end <- Sys.time()
  candidate_terms <- split_formula(refmodel$formula,
      data = refmodel$fetch_data(),
      add_main_effects = FALSE
  )
  ## with `match` we get the indices of the variables as they enter the
  ## solution path in solution_terms
  solution <- match(solution_terms, candidate_terms[-1])
  solution_terms_cv[, seq_along(solution)] <- solution
  if (length(solution) < (nterms_max - 1)) {
    not_in_solution <- setdiff(seq_len(nterms_max - 1), seq_along(solution))
    solution_terms_cv[, not_in_solution] <- NA
  }

  ## find out how many of cross-validated iterations select
  ## the same variables as the selection with all the data.
  solution_terms_cv_ch <- sapply(
    seq_len(NROW(solution_terms_cv)),
    function(i) {
      if (!is.character(solution_terms_cv[i, ])) {
        unlist(candidate_terms[-1])[solution_terms_cv[i, ]]
      } else {
        solution_terms_cv[i, ]
      }
    }
  )

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
    data = NULL
  )
  sel_loo <- nlist(
    search_path = object,
    solution_terms,
    summaries,
    d_test,
    kl = sapply(submodels, function(x) x$kl),
    control = nlist(
      time = difftime(end, start, units = "secs"),
      ndraws = object$control$ndraws,
      nclusters = object$control$nclusters,
      ndraws_pred,
      nclusters_pred,
      cv_method = "LOO"
    )
  )
  class(sel_loo) <- c("vselapproxcv", "vsel")
  return(sel_loo)
}

#' @rdname workflow
#' @export
approximate_kfold <- function(object, ...) {
    UseMethod("approximate_kfold")
}

#' @rdname workflow
#' @export
approximate_kfold.default <- function(object, ...) {
    return(approximate_kfold(varsel_search(object, ...)))
}

#' @rdname workflow
#' @export
approximate_kfold.refmodel <- function(object, ...) {
    return(approximate_kfold(varsel_search(object, ...)))
}

#' @rdname workflow
#' @export
approximate_kfold.vselsearch <- function(object,
                                         ndraws_pred = NULL,
                                         nclusters_pred = NULL,
                                         K = NULL,
                                         verbose = TRUE,
                                         penalty = NULL,
                                         refit_proj = TRUE,
                                         seed = NULL,
                                         ...) {
  refmodel <- object$refmodel
  family <- refmodel$family
  nterms_max <- object$object$nterms_max

  solution_terms <- object$solution_terms
  p_sel <- object$p_sel

  ## fetch the default arguments or replace them by the user defined values
  args <- parse_args_varsel(
      refmodel, NULL, NULL, NULL, NULL,
      NULL, NULL, ndraws_pred, nclusters_pred, NULL
  )
  nclusters_pred <- args$nclusters_pred
  ndraws_pred <- args$ndraws_pred

  args <- parse_args_cv_varsel(
    refmodel, "kfold", K, NULL,
    nclusters_pred
  )
  K <- args$K

  ## the clustering/subsampling used for prediction
  p_pred <- .get_refdist(refmodel,
      ndraws = ndraws_pred,
      nclusters = nclusters_pred
  )

  start <- Sys.time()
  ## fetch the k_fold list (or compute it now if not already computed)
  k_fold <- .get_kfold(refmodel, K, verbose, seed, approximate = TRUE)

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

  ## List of K elements, each containing d_train, p_pred, etc. corresponding
  ## to each fold.
  make_list_cv <- function(refmodel, d_test, msg, K) {
      nclusters_pred <- min(
          refmodel$nclusters_pred,
          nclusters_pred
      )
      newdata <- d_test$newdata
      pred <- refmodel$ref_predfun(refmodel$fit, newdata = newdata)
      pred <- matrix(
          as.numeric(pred),
          nrow = NROW(pred), ncol = NCOL(pred)
      )
      mu_test <- family$linkinv(pred)
      nlist(refmodel, mu_test,
          dis = refmodel$dis,
          w_test = refmodel$wsample, d_test, msg,
          fold_index = K
      )
  }

  msgs <- paste0(object$method, " search for fold ", seq_len(K), "/", K, ".")
  list_cv <- mapply(make_list_cv, refmodels_cv, d_test_cv, msgs,
      SIMPLIFY = FALSE, seq_len(K)
  )

  ## Construct submodel projections for each fold
  if (verbose) {
      print(paste0(
        "Computing kfold for ", K, " folds for ",
        length(solution_terms), " models."
      ))
      pb <- utils::txtProgressBar(
          min = 0, max = K,
          style = 3, initial = 0
      )
  }

  p_sub <- .get_submodels(
      search_path = object, nterms = c(0, seq_along(solution_terms)),
      family = family, p_ref = p_pred, refmodel = refmodel,
      intercept = TRUE, regul = object$control$opt$regul, cv_search = refit_proj
  )

  ## Helper function extract and combine mu and lppd from K lists with each
  ## n/K of the elements to one list with n elements
  hf <- function(x) as.list(do.call(rbind, x))

  ## Apply some magic to manipulate the structure of the list so that instead
  ## of list with K sub_summaries each containing n/K mu:s and lppd:s, we have
  ## only one sub_summary-list that contains with all n mu:s and lppd:s.
  get_summaries_submodel_cv <- function(fold) {
      omitted <- fold$d_test$omitted
      fold_summaries <- .get_sub_summaries(
          submodels = p_sub, test_points = omitted, refmodel = refmodel,
          family = family
      )
      summ <- lapply(fold_summaries, data.frame)
      if (verbose) {
          utils::setTxtProgressBar(pb, fold$fold_index)
      }
      return(summ)
  }
  sub_cv_summaries <- mapply(get_summaries_submodel_cv, list_cv)
  sub <- apply(sub_cv_summaries, 1, hf)
  sub <- lapply(sub, function(summ) {
      summ$w <- rep(1, length(summ$mu))
      summ$w <- summ$w / sum(summ$w)
      summ
  })

  if (verbose) {
      close(pb)
  }

  ref <- hf(lapply(list_cv, function(fold) {
      data.frame(.weighted_summary_means(
          y_test = fold$d_test, family = family,
          wsample = fold$refmodel$wsample,
          mu = fold$mu_test, dis = fold$refmodel$dis
      ))
  }))

  ## Combine also the K separate test data sets into one list
  ## with n y's and weights's.
  d_cv <- hf(lapply(d_test_cv, function(fold) {
      data.frame(
          y = fold$y, weights = fold$weights,
          test_points = fold$omitted
      )
  }))

  end <- Sys.time()
  sel_kfold <- nlist(
    search_path = object,
    solution_terms,
    summaries = list(sub = sub, ref = ref),
    d_test = c(d_cv, type = "kfold"),
    kl = sapply(p_sub, function(x) x$kl),
    control = nlist(
      time = difftime(end, start, units = "secs"),
      ndraws = object$control$ndraws,
      nclusters = object$control$nclusters,
      ndraws_pred,
      nclusters_pred,
      cv_method = "KFold"
    )
  )
  class(sel_kfold) <- c("vselapproxcv", "vsel")
  return(sel_kfold)
}

#' @export
diagnostic <- function(x, ...) {
  diff <- x$diff[NROW(x)]
  diff_se <- x$diff_se[NROW(x)]
  if (diff > diff_se) {
    return(paste0(
      "The projections' ELPDs seems overoptimistic, we recommend ",
      "increasing `ndraws_pred`."
    ))
  } else {
    return(paste0(
      "The projections' ELPDs match the reference model's."
    ))
  }
}

#' @export
summary.vselapproxcv <- function(object, stats = "elpd",
                                 type = c("mean", "se", "diff", "diff_se"),
                                 alpha = 0.32, baseline = NULL, deltas = FALSE,
                                 digits = 1, ...) {
  search_path <- object$search_path
  refmodel <- object$search_path$refmodel
  out <- list(
    formula = refmodel$formula,
    fit = refmodel$fit,
    family = refmodel$family,
    nobs = NROW(refmodel$fetch_data()),
    ndraws = object$control$ndraws,
    nclusters = object$control$nclusters,
    ndraws_pred = object$control$ndraws_pred,
    nclusters_pred = object$control$nclusters_pred,
    time = object$control$time,
    cv_method = object$control$cv_method,
    solution_terms = search_path$solution_terms,
    search_included = "search not included"
  )
  class(out) <- "vselapproxcvsummary"

  stats_table <- summary_stats_table(
    object, refmodel, stats, type, alpha,
    baseline, deltas
  )

  out$stats_table <- stats_table
  out$diagnostic <- diagnostic(stats_table)
  return(out)
}

#' @export
#' @method print vselapproxcv
print.vselapproxcv <- function(x, digits = 1, ...) {
  summary <- summary.vselapproxcv(x, digits = digits, ...)
  print(summary)
  return(invisible(summary))
}

#' @export
print.vselapproxcvsummary <- function(x, digits = 1, ...) {
  cat(paste0(
      "Approximate ", x$cv_method, " CV selection took ",
      round(x$time, digits), " seconds.\n"
  ))
  print(x$family)
  cat("Formula: ")
  print(x$formula)
  cat(paste0("Observations: ", x$nobs, "\n"))
  if (!is.null(x$cv_method)) {
    cat(paste("CV method:", x$cv_method, x$search_included, "\n"))
  }
  nterms_max <- max(x$stats_table$size)
  cat(paste0(
    "Draws used for selection: ", x$ndraws, ", in ",
    x$nclusters, " clusters\n"
  ))
  cat(paste0(
    "Draws used for prediction: ", x$ndraws_pred, ", in ",
    x$nclusters_pred, " clusters\n"
  ))
  cat(paste0("\nDiagnostics:\n", x$diagnostic, "\n"))
  cat("\nSelection Summary:\n")
  print(x$stats_table %>% dplyr::mutate(dplyr::across(
    where(is.numeric),
    ~ round(., digits)
  )),
  row.names = FALSE
  )
  return(invisible(x))
}


#' @rdname workflow
#' @export
cv_loo <- function(object, ...) {
    UseMethod("cv_loo")
}

#' @rdname workflow
#' @export
cv_loo.default <- function(object, ...) {
    cv_loo(varsel_search(object, ...))
}

#' @rdname workflow
#' @export
cv_loo.refmodel <- function(object, ...) {
    cv_loo(varsel_search(object, ...))
}

#' @rdname workflow
#' @export
cv_loo.vselsearch <- function(object, ...) {
    sel_cv <- nlist(search_path = object)
    class(sel_cv) <- "vselapproxcv"
    return(cv_loo.vselapproxcv(sel_cv, ...))
}

#' @rdname workflow
#' @importFrom doRNG %dorng%
#' @importFrom doParallel %dopar%
#' @export
cv_loo.vselapproxcv <- function(object,
                                method = NULL,
                                ndraws = NULL,
                                nclusters = NULL,
                                ndraws_pred = NULL,
                                nclusters_pred = NULL,
                                nterms_max = NULL,
                                penalty = NULL,
                                verbose = TRUE,
                                nloo = NULL,
                                lambda_min_ratio = 1e-5,
                                nlambda = 150,
                                thresh = 1e-6,
                                regul = 1e-4,
                                seed = NULL,
                                cores = parallel::detectCores(),
                                search_terms = NULL,
                                ...) {
  search_path_prev <- object$search_path
  refmodel <- search_path_prev$refmodel
  family <- refmodel$family
  nterms_max <- length(search_path_prev$solution_terms) + 1

  ## fetch the default arguments or replace them by the user defined values
  args <- parse_args_varsel(
    refmodel, method, NULL, NULL, NULL,
    ndraws, nclusters, ndraws_pred, nclusters_pred, NULL
  )
  nclusters <- args$nclusters
  ndraws <- args$ndraws
  nclusters_pred <- args$nclusters_pred
  ndraws_pred <- args$ndraws_pred
  search_terms <- args$search_terms
  method <- args$method

  opt <- nlist(lambda_min_ratio, nlambda, thresh, regul)
  mu <- refmodel$mu
  dis <- refmodel$dis

  ## the clustering/subsampling used for selection
  p_sel <- .get_refdist(refmodel,
    ndraws = ndraws,
    nclusters = nclusters
  )
  cl_sel <- p_sel$cl # clustering information

  ## the clustering/subsampling used for prediction
  p_pred <- .get_refdist(refmodel,
    ndraws = ndraws_pred,
    nclusters = nclusters_pred
  )
  cl_pred <- p_pred$cl

  ## fetch the log-likelihood for the reference model to obtain the LOO
  ## weights
  if (is.null(refmodel$loglik)) {
    ## case where log-likelihood not available, i.e., the reference model is
    ## not a genuine model => cannot compute LOO
    stop(
      "LOO can be performed only if the reference model is a genuine ",
      "probabilistic model for which the log-likelihood can be evaluated."
    )
  } else {
    ## log-likelihood available
    loglik <- refmodel$loglik
  }
  psisloo <- loo::psis(-loglik, cores = 1, r_eff = rep(1, ncol(loglik)))
  lw <- weights(psisloo)
  pareto_k <- loo::pareto_k_values(psisloo)
  n <- length(pareto_k)
  ## by default use all observations
  nloo <- min(nloo, n)

  if (nloo < 0) {
    stop("nloo must be at least 1")
  }

  ## compute loo summaries for the reference model
  loo_ref <- apply(loglik + lw, 2, log_sum_exp)
  mu_ref <- rep(0, n)
  for (i in seq_len(n)) {
    mu_ref[i] <- mu[i, ] %*% exp(lw[, i])
  }

  doFuture::registerDoFuture()
  if (cores > 1) {
    future::plan(future::multicore, workers = cores)
  } else {
    future::plan(future::sequential, workers = cores)
  }

  if (verbose) {
    msg <- paste("Repeating", method, "search for", nloo, "LOO folds...")
    print(msg)
    pb <- utils::txtProgressBar(min = 0, max = nloo, style = 3, initial = 0)
  }

  start <- Sys.time()
  validset <- .loo_subsample_pps(nloo, loo_ref, seed)
  inds <- validset$inds

  ## solution_terms_cv <- matrix(nrow = n, ncol = nterms_max - 1)
  ## loo_sub <- matrix(nrow = n, ncol = nterms_max)
  ## mu_sub <- matrix(nrow = n, ncol = nterms_max)

  comb <- function(r, ...) {
    for (pr in list(...)) {
      i <- pr$i
      r$mu_sub[i, ] <- pr$mu_sub
      r$loo_sub[i, ] <- pr$loo_sub
      r$solution_terms_cv[i, ] <- pr$solution_terms_cv
    }
    return(r)
  }

  candidate_terms <- split_formula(refmodel$formula,
    data = refmodel$fetch_data(),
    add_main_effects = FALSE
  )

  ## for (run_index in seq_along(inds)) {
  r <- foreach::foreach(
    run_index = seq_along(inds),
    .combine = "comb",
    .multicombine = TRUE,
    .init = list(
      mu_sub = matrix(nrow = n, ncol = nterms_max),
      loo_sub = matrix(nrow = n, ncol = nterms_max),
      solution_terms_cv = matrix(nrow = n, ncol = nterms_max - 1)
    )
  ) %dorng% {
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
      family = family, intercept = TRUE, nterms_max = nterms_max,
      penalty = penalty, verbose = FALSE, opt = opt,
      search_terms = search_terms
    )
    solution_terms <- search_path$solution_terms

    ## project onto the selected models and compute the prediction accuracy
    ## for the left-out point
    submodels <- .get_submodels(
      search_path = search_path, nterms = c(0, seq_along(solution_terms)),
      family = family, p_ref = p_pred, refmodel = refmodel,
      intercept = TRUE, regul = opt$regul, cv_search = FALSE
    )
    summaries_sub <- .get_sub_summaries(
      submodels = submodels, test_points = c(i), refmodel = refmodel,
      family = family
    )

    loo_sub <- sapply(seq_along(summaries_sub), function(k) {
      summaries_sub[[k]]$lppd
    })
    mu_sub <- sapply(seq_along(summaries_sub), function(k) {
      summaries_sub[[k]]$mu
    })

    ## with `match` we get the indices of the variables as they enter the
    ## solution path in solution_terms
    solution <- match(solution_terms, candidate_terms[-1])
    solution_terms_cv <- rep(NA, nterms_max - 1)
    solution_terms_cv[seq_along(solution)] <- solution
    if (length(solution) < (nterms_max - 1)) {
      not_in_solution <- setdiff(
        seq_len(nterms_max - 1),
        seq_along(solution)
      )
      solution_terms_cv[not_in_solution] <- NA
    }

    if (verbose) {
      utils::setTxtProgressBar(pb, run_index)
    }
    return(nlist(i, mu_sub, loo_sub, solution_terms_cv))
  }
  end <- Sys.time()

  mu_sub <- r$mu_sub
  loo_sub <- r$loo_sub
  solution_terms_cv <- r$solution_terms_cv

  if (verbose) {
    ## close the progress bar object
    close(pb)
  }

  ## find out how many of cross-validated iterations select
  ## the same variables as the selection with all the data.
  solution_terms_cv_ch <- sapply(
    seq_len(NROW(solution_terms_cv)),
    function(i) {
      if (!is.character(solution_terms_cv[i, ])) {
        unlist(candidate_terms[-1])[solution_terms_cv[i, ]]
      } else {
        solution_terms_cv[i, ]
      }
    }
  )

  solution_terms <- search_path_prev$solution_terms
  ## these weights might be non-constant in case of subsampling LOO
  sel_solution_terms <- solution_terms_cv
  ## if weights are not set, then all validation folds have equal weight
  w <- rep(1, NROW(sel_solution_terms))
  w <- w / sum(w)
  vars <- unlist(solution_terms)
  pct_solution_terms_cv <- t(sapply(
    seq_along(solution_terms),
    function(size) {
      c(
        size = size,
        sapply(vars, function(var) {
          sum((solution_terms_cv_ch[seq_len(size), ,
            drop = FALSE
          ] == var) * w,
          na.rm = TRUE
          )
        })
      )
    }
  ))

  ## put all the results together in the form required by cv_varsel
  summ_sub <- lapply(seq_len(nterms_max), function(k) {
    list(lppd = loo_sub[, k], mu = mu_sub[, k], w = validset$w)
  })
  summ_ref <- list(lppd = loo_ref, mu = mu_ref)
  summaries <- list(sub = summ_sub, ref = summ_ref)

  search_path_prev$refmodel <- refmodel
  d_test <- list(
    y = refmodel$y, type = "LOO",
    test_points = seq_along(refmodel$y),
    weights = refmodel$wobs,
    data = NULL
  )
  sel_cv_loo <- nlist(
    solution_terms_cv,
    pct_solution_terms_cv,
    solution_terms,
    summaries,
    d_test,
    search_path = search_path_prev,
    kl = object$kl,
    control = nlist(
      nterms_max,
      method,
      time = difftime(end, start, units = "secs"),
      ndraws,
      nclusters,
      ndraws_pred,
      nclusters_pred,
      cv_method = "LOO"
    )
  )
  class(sel_cv_loo) <- c("vselcv", "vsel")
  return(sel_cv_loo)
}

#' @rdname workflow
#' @export
cv_kfold <- function(object, ...) {
    UseMethod("cv_kfold")
}

#' @rdname workflow
#' @export
cv_kfold.default <- function(object, ...) {
    cv_kfold(varsel_search(object, ...), ...)
}

#' @rdname workflow
#' @export
cv_kfold.refmodel <- function(object, ...) {
    cv_kfold(varsel_search(object, ...), ...)
}

#' @rdname workflow
#' @export
cv_kfold.vselsearch <- function(object, ...) {
    sel_cv <- nlist(search_path = object)
    class(sel_cv) <- "vselapproxcv"
    return(cv_kfold.vselapproxcv(sel_cv, ...))
}

cv_kfold.vselapproxcv <- function(object,
                                  method = NULL,
                                  ndraws = NULL,
                                  nclusters = NULL,
                                  K = NULL,
                                  ndraws_pred = NULL,
                                  nclusters_pred = NULL,
                                  nterms_max = NULL,
                                  penalty = NULL,
                                  verbose = TRUE,
                                  nloo = NULL,
                                  lambda_min_ratio = 1e-5,
                                  nlambda = 150,
                                  thresh = 1e-6,
                                  regul = 1e-4,
                                  cores = parallel::detectCores(),
                                  seed = NULL,
                                  search_terms = NULL,
                                  ...) {
  search_path <- object$search_path
  refmodel <- search_path$refmodel
  family <- refmodel$family
  solution_terms <- search_path$solution_terms
  nterms_max <- length(solution_terms) + 1

  ## fetch the default arguments or replace them by the user defined values
  args <- parse_args_varsel(
      refmodel, method, NULL, NULL, NULL, ndraws,
      nclusters, ndraws_pred, nclusters_pred, NULL
  )
  nclusters <- args$nclusters
  ndraws <- args$ndraws
  nclusters_pred <- args$nclusters_pred
  ndraws_pred <- args$ndraws_pred
  search_terms <- args$search_terms
  method <- args$method

  opt <- nlist(lambda_min_ratio, nlambda, thresh, regul)

  ## fetch the k_fold list (or compute it now if not already computed)
  k_fold <- .get_kfold(refmodel, K, verbose, seed)

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

  ## List of K elements, each containing d_train, p_pred, etc. corresponding
  ## to each fold.
  make_list_cv <- function(refmodel, d_test, msg) {
      nclusters_pred <- min(
          refmodel$nclusters_pred,
          nclusters_pred
      )
      p_sel <- .get_refdist(refmodel, ndraws, nclusters)
      p_pred <- .get_refdist(refmodel, ndraws_pred, nclusters_pred)
      newdata <- d_test$newdata
      pred <- refmodel$ref_predfun(refmodel$fit, newdata = newdata)
      pred <- matrix(
          as.numeric(pred),
          nrow = NROW(pred), ncol = NCOL(pred)
      )
      mu_test <- family$linkinv(pred)
      nlist(refmodel, p_sel, p_pred, mu_test,
          dis = refmodel$dis, w_test = refmodel$wsample, d_test, msg
      )
  }

  msgs <- paste0(method, " search for fold ", 1:K, "/", K, ".")
  list_cv <- mapply(make_list_cv, refmodels_cv, d_test_cv, msgs,
      SIMPLIFY = FALSE
  )

  if (cores > 1) {
    future::plan(future::multicore, workers = cores)
  } else {
    future::plan(future::sequential, workers = cores)
  }

  ## Perform the selection for each of the K folds
  start <- Sys.time()
  if (verbose) {
      print(paste0("Repeating ", method, " search for ", K, " folds.."))
      pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial = 0)
  }
  search_path_cv <- future.apply::future_lapply(
    seq_along(list_cv),
    function(fold_index) {
      fold <- list_cv[[fold_index]]
      family <- fold$refmodel$family
      out <- select(
        method = method, p_sel = fold$p_sel, refmodel = fold$refmodel,
        family = family, intercept = TRUE, nterms_max = nterms_max,
        penalty = penalty, verbose = FALSE, opt = opt,
        search_terms = search_terms
      )
      if (verbose) {
        utils::setTxtProgressBar(pb, fold_index)
      }
      out
    }
  )

  solution_terms_cv <- t(sapply(search_path_cv, function(e) e$solution_terms))
  if (verbose) {
      close(pb)
  }

  ## Construct submodel projections for each fold
  if (verbose) {
      print("Computing projections..")
      pb <- utils::txtProgressBar(min = 0, max = K, style = 3, initial = 0)
  }

  get_submodels_cv <- function(search_path, fold_index) {
      fold <- list_cv[[fold_index]]
      family <- fold$refmodel$family
      solution_terms <- search_path$solution_terms
      p_sub <- .get_submodels(
          search_path = search_path, nterms = c(0, seq_along(solution_terms)),
          family = family, p_ref = fold$p_pred, refmodel = fold$refmodel,
          intercept = TRUE, regul = opt$regul, cv_search = TRUE
      )
      if (verbose) {
          utils::setTxtProgressBar(pb, fold_index)
      }
      return(p_sub)
  }

  p_sub_cv <- future.apply::future_mapply(get_submodels_cv,
    search_path_cv,
    seq_along(list_cv),
    future.seed = seed,
    SIMPLIFY = FALSE
  )
  if (verbose) {
    close(pb)
  }

  ## Helper function extract and combine mu and lppd from K lists with each
  ## n/K of the elements to one list with n elements
  hf <- function(x) as.list(do.call(rbind, x))

  ## Apply some magic to manipulate the structure of the list so that instead of
  ## list with K sub_summaries each containing n/K mu:s and lppd:s, we have only
  ## one sub_summary-list that contains with all n mu:s and lppd:s.
  get_summaries_submodel_cv <- function(p_sub, fold) {
      omitted <- fold$d_test$omitted
      fold_summaries <- .get_sub_summaries(
          submodels = p_sub, test_points = omitted, refmodel = refmodel,
          family = family
      )
      summ <- lapply(fold_summaries, data.frame)
      return(summ)
  }
  sub_cv_summaries <- future.apply::future_mapply(
    get_summaries_submodel_cv,
    p_sub_cv, list_cv
  )
  sub <- apply(sub_cv_summaries, 1, hf)
  sub <- lapply(sub, function(summ) {
      summ$w <- rep(1, length(summ$mu))
      summ$w <- summ$w / sum(summ$w)
      summ
  })

  ref <- hf(lapply(list_cv, function(fold) {
      data.frame(.weighted_summary_means(
          y_test = fold$d_test, family = family,
          wsample = fold$refmodel$wsample,
          mu = fold$mu_test, dis = fold$refmodel$dis
      ))
  }))

  end <- Sys.time()
  ## Combine also the K separate test data sets into one list
  ## with n y's and weights's.
  d_cv <- hf(lapply(d_test_cv, function(fold) {
      data.frame(
          y = fold$y, weights = fold$weights,
          test_points = fold$omitted
      )
  }))

  ## these weights might be non-constant in case of subsampling LOO
  sel_solution_terms <- t(solution_terms_cv)
  ## if weights are not set, then all validation folds have equal weight
  w <- rep(1, NCOL(sel_solution_terms))
  w <- w / sum(w)
  vars <- unlist(solution_terms)
  pct_solution_terms_cv <- t(sapply(
      seq_along(solution_terms),
      function(size) {
          c(
              size = size,
              sapply(vars, function(var) {
                  sum((sel_solution_terms[seq_len(size), ,
                      drop = FALSE
                  ] == var) * w,
                  na.rm = TRUE
                  )
              })
          )
      }
  ))

  search_path$refmodel <- refmodel
  sel_cv_kfold <- nlist(
    solution_terms_cv,
    pct_solution_terms_cv,
    solution_terms,
    summaries = nlist(ref, sub),
    search_path,
    kl = object$kl,
    d_test = c(d_cv, type = "kfold"),
    control = nlist(
      nterms_max,
      method,
      time = difftime(end, start, units = "secs"),
      ndraws,
      nclusters,
      ndraws_pred,
      nclusters_pred,
      cv_method = "KFold"
    )
  )
  class(sel_cv_kfold) <- c("vselcv", "vsel")
  return(sel_cv_kfold)
}

#' @export
summary.vselcv <- function(object, stats = "elpd",
                           type = c("mean", "se", "diff", "diff_se"),
                           alpha = 0.32, baseline = NULL, deltas = FALSE,
                           digits = 1, ...) {
  search_path <- object$search_path
  refmodel <- object$search_path$refmodel
  out <- list(
    formula = refmodel$formula,
    fit = refmodel$fit,
    family = refmodel$family,
    nobs = NROW(refmodel$fetch_data()),
    ndraws = object$control$ndraws,
    nclusters = object$control$nclusters,
    ndraws_pred = object$control$ndraws_pred,
    nclusters_pred = object$control$nclusters_pred,
    time = object$control$time,
    cv_method = object$control$cv_method,
    method = object$control$method,
    nterms_max = object$control$nterms_max,
    solution_terms = search_path$solution_terms,
    search_included = "search included"
  )
  class(out) <- "vselcvsummary"

  stats_table <- summary_stats_table(
    object, refmodel, stats, type, alpha,
    baseline, deltas
  )
  out$stats_table <- stats_table
  out$diagnostic <- diagnostic(stats_table)
  
  return(out)
}

#' @export
#' @method print vselcv
print.vselcv <- function(x, digits = 1, ...) {
  summary <- summary.vselcv(x, digits = digits, ...)
  print(summary)
  return(invisible(summary))
}

#' @export
print.vselcvsummary <- function(x, digits = 1, ...) {
  cat(paste0(
      "Full ", x$cv_method, " CV selection took ",
      round(x$time, digits), " seconds.\n"
  ))
  print(x$family)
  cat("Formula: ")
  print(x$formula)
  cat(paste0("Observations: ", x$nobs, "\n"))
  cat(paste0(
      "Search method: ", x$method, ", maximum number of terms ",
      x$nterms_max, "\n"
  ))
  if (!is.null(x$cv_method)) {
    cat(paste("CV method:", x$cv_method, x$search_included, "\n"))
  }
  nterms_max <- max(x$stats_table$size)
  cat(paste0(
    "Draws used for selection: ", x$ndraws, ", in ",
    x$nclusters, " clusters\n"
  ))
  cat(paste0(
    "Draws used for prediction: ", x$ndraws_pred, ", in ",
    x$nclusters_pred, " clusters\n"
  ))
  cat(paste0("\nDiagnostics:\n", x$diagnostic, "\n"))
  cat("\nSelection Summary:\n")
  print(x$stats_table %>% dplyr::mutate(dplyr::across(
    where(is.numeric),
    ~ round(., digits)
  )),
  row.names = FALSE
  )
  return(invisible(x))
}

#' @rdname workflow
#' @export
varsel_cv <- function(object, ...) {
    UseMethod("varsel_cv")
}

#' @rdname workflow
#' @export
varsel_cv.default <- function(object, ...) {
    varsel_cv(varsel_search(object, ...), ...)
}

#' @rdname workflow
#' @export
varsel_cv.refmodel <- function(object, ...) {
    varsel_cv(varsel_search(object, ...), ...)
}

#' @rdname workflow
#' @export
varsel_cv.vselsearch <- function(object, ...) {
    sel_cv <- nlist(search_path = object)
    class(sel_cv) <- "vselapproxcv"
    return(varsel_cv.vselapproxcv(sel_cv, ...))
}

#' @rdname workflow
#' @export
varsel_cv.vselapproxcv <- function(object,
                                   method = NULL,
                                   cv_method = NULL,
                                   ndraws = NULL,
                                   nclusters = NULL,
                                   ndraws_pred = NULL,
                                   nclusters_pred = NULL,
                                   nterms_max = NULL,
                                   K = NULL,
                                   penalty = NULL,
                                   verbose = TRUE,
                                   nloo = NULL,
                                   lambda_min_ratio = 1e-5,
                                   nlambda = 150,
                                   thresh = 1e-6,
                                   regul = 1e-4,
                                   seed = NULL,
                                   search_terms = NULL,
                                   cv_search = TRUE,
                                   ...) {
  search_path <- object$search_path
  refmodel <- search_path$refmodel
  ## resolve the arguments similar to varsel
  args <- parse_args_varsel(
    refmodel = refmodel, method = method, cv_search = cv_search,
    intercept = TRUE, nterms_max = nterms_max, nclusters = nclusters,
    ndraws = ndraws, nclusters_pred = nclusters_pred,
    ndraws_pred = ndraws_pred, search_terms = search_terms
  )
  method <- args$method
  cv_search <- args$cv_search
  nterms_max <- args$nterms_max
  nclusters <- args$nclusters
  ndraws <- args$ndraws
  nclusters_pred <- args$nclusters_pred
  ndraws_pred <- args$ndraws_pred
  search_terms <- args$search_terms

  ## arguments specific to this function
  args <- parse_args_cv_varsel(
    refmodel, cv_method, K, nclusters,
    nclusters_pred
  )
  cv_method <- args$cv_method
  K <- args$K
  nclusters <- args$nclusters
  nclusters_pred <- args$nclusters_pred

  if (tolower(cv_method) == "loo") {
    if (cv_search == TRUE) {
      sel_cv <- cv_loo(object,
          method = method, ndraws = ndraws, nclusters = nclusters,
          ndraws_pred = ndraws_pred, nclusters_pred = nclusters_pred,
          nterms_max = nterms_max, nloo = nloo, penalty = penalty,
          lambda_min_ratio = lambda_min_ratio, thresh = thresh,
          regul = regul, seed = seed, search_terms = search_terms
      )
    } else {
      sel_cv <- approximate_loo(search_path,
          ndraws_pred = ndraws_pred, penalty = penalty,
          nclusters_pred = nclusters_pred, nloo = nloo
      )
    }
  } else if (tolower(cv_method) == "kfold") {
    if (cv_search == TRUE) {
      sel_cv <- cv_kfold(object,
          method = method, ndraws = ndraws, nclusters = nclusters,
          ndraws_pred = ndraws_pred, nclusters_pred = nclusters_pred,
          nterms_max = nterms_max, nloo = nloo, penalty = penalty,
          lambda_min_ratio = lambda_min_ratio, thresh = thresh,
          regul = regul, seed = seed, search_terms = search_terms,
          K = K
      )
    } else {
      sel_cv <- approximate_kfold(search_path,
          ndraws_pred = ndraws_pred, penalty = penalty,
          nclusters_pred = nclusters_pred, nloo = nloo, K = K
      )
    }
  }

  return(sel_cv)
}

summary_stats_table <- function(object, refmodel, stats = "elpd",
                                type = c("mean", "se", "diff", "diff_se"),
                                alpha = 0.32, baseline = NULL, deltas = FALSE) {
  if (!inherits(object, "vselapproxcv") && !inherits(object, "vselcv")
      && !inherits(object, "vsel")) {
    stop("Can only get summary stats table for vsel objects with summaries.")
  }
  baseline <- .validate_baseline(refmodel, baseline, deltas)
  if (deltas) {
    nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1])
    tab <- .tabulate_stats(object, stats,
                           alpha = alpha, nfeat_baseline = nfeat_baseline
                           )
  } else {
    tab <- .tabulate_stats(object, stats, alpha = alpha)
  }
  stats_table <- subset(tab, tab$size != Inf) %>%
    dplyr::group_by(statistic) %>%
    dplyr::slice_head(n = length(object$solution_terms) + 1)

  if (deltas) {
    type <- setdiff(type, c("diff", "diff_se"))
  }
  ## these are the corresponding names for mean, se, upper and lower in the
  ## stats_table, and their suffices in the table to be returned
  qty <- unname(sapply(type, function(t) {
    switch(t, mean = "value", upper = "uq", lower = "lq", se = "se",
           diff = "diff", diff_se = "diff_se")
  }))
  if (!is.null(object$cv_method)) {
    cv_suffix <- unname(switch(object$cv_method,
      LOO = "_loo", kfold = "_kfold"
    ))
  } else {
    cv_suffix <- NULL
  }

  if (length(stats) > 1) {
    suffix <- lapply(stats, function(s) {
      paste0(
        s,
        unname(sapply(type, function(t) {
          switch(t, mean = cv_suffix, upper = "_upper", lower = "_lower",
            se = "_se", diff = "_diff", diff_se = "_diff_se"
          )
        }))
      )
    })
  } else {
    suffix <- list(unname(sapply(type, function(t) {
      switch(t, mean = paste0(stats, cv_suffix), upper = "upper",
        lower = "lower", se = "se",
        diff = "diff", diff_se = "diff_se"
      )
    })))
  }

  ## loop through all the required statistics
  arr <- data.frame(
    size = unique(stats_table$size),
    solution_terms = c(NA, object$search_path$solution_terms)
  )
  for (i in seq_along(stats)) {
    temp <- subset(stats_table, stats_table$statistic == stats[i], qty)
    newnames <- suffix[[i]]
    colnames(temp) <- newnames
    arr <- cbind(arr, temp)
  }

  return(arr)
}
