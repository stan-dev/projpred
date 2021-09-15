#' Predictions from a projected submodel
#'
#' [proj_linpred()] gives draws of the linear predictor (possibly transformed to
#' response scale) of a projected submodel (i.e., a submodel resulting from
#' projecting the reference model onto it). [proj_predict()] draws from the
#' predictive distribution of a projected submodel. If the projection has not
#' been performed, both functions also perform the projection. Both functions
#' can also handle multiple projected submodels at once (if the input object is
#' of class `vsel`).
#'
#' @name projection-linpred-predict
#'
#' @template args-newdata
#' @param object Either an object returned by [project()] or alternatively any
#'   object that can be passed to argument `object` of [project()].
#' @param filter_nterms Only applies if `object` is an object returned by
#'   [project()]. In that case, `filter_nterms` can be used to filter `object`
#'   for only those elements (submodels) with a number of solution terms in
#'   `filter_nterms`. Therefore, needs to be a numeric vector or `NULL`. If
#'   `NULL`, use all submodels.
#' @param transform For [proj_linpred()] only. A single logical value indicating
#'   whether the linear predictor should be transformed using the inverse-link
#'   function (`TRUE`) or not (`FALSE`).
#' @param integrated For [proj_linpred()] only. A single logical value
#'   indicating whether the output should be averaged over the projected
#'   posterior draws (`TRUE`) or not (`FALSE`).
#' @param nresample_clusters For [proj_predict()] with clustered projection
#'   only. Number of draws to return from the predictive distribution of the
#'   projection. Not to be confused with argument `nclusters` of [project()]:
#'   `nresample_clusters` gives the number of draws (*with* replacement) from
#'   the set of clustered posterior draws after projection (as determined by
#'   argument `nclusters` of [project()]).
#' @param .seed For [proj_predict()] only. Pseudorandom number generation (PRNG)
#'   seed by which the same results can be obtained again if needed. If `NULL`,
#'   no seed is set and therefore, the results are not reproducible. See
#'   [set.seed()] for details. Here, this seed is used for drawing from the
#'   predictive distribution of the submodel(s) onto which the reference model
#'   was projected. If a clustered projection was performed, `.seed` is also
#'   used for drawing from the set of the projected clusters of posterior draws
#'   (see argument `nresample_clusters`).
#' @param ... Additional arguments passed to [project()] if `object` is not
#'   already an object returned by [project()].
#'
#' @return If the prediction is done for one submodel only (i.e., `nterms` has
#'   length one or `solution_terms` is specified):
#'   * [proj_linpred()] returns a `list` with elements `pred` (predictions) and
#'   `lpd` (log predictive densities). Each of these two elements is a \eqn{S
#'   \times N}{S x N} matrix.
#'   * [proj_predict()] returns a \eqn{S \times N}{S x N} matrix of predictions.
#'
#'   Thereby, \eqn{S} denotes the number of (possibly clustered) projected
#'   posterior draws and \eqn{N} denotes the number of observations.
#'
#'   If the prediction is done for more than one submodel, the output from above
#'   is returned for each submodel, giving a named `list` with one element for
#'   each submodel (the names of this `list` being the numbers of solutions
#'   terms of the submodels when counting the intercept, too).
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
#'   # Projection onto an arbitrary combination of predictor terms (with a small
#'   # value for `nclusters`, but only for the sake of speed in this example;
#'   # this is not recommended in general):
#'   prj <- project(fit, solution_terms = c("X1", "X3", "X5"), nclusters = 10,
#'                  seed = 9182)
#'
#'   # Predictions (at the training points) from the projected submodel:
#'   prjl <- proj_linpred(prj)
#'   prjp <- proj_predict(prj, .seed = 7364)
#' }
#'
NULL

## The 'helper' for proj_linpred and proj_predict, ie. does all the
## functionality that is common to them. It essentially checks all the arguments
## and sets them to their respective defaults and then loops over the
## projections. For each projection, it evaluates the fun-function, which
## calculates the linear predictor if called from proj_linpred and samples from
## the predictive distribution if called from proj_predict.
proj_helper <- function(object, newdata,
                        offsetnew, weightsnew,
                        onesub_fun, filter_nterms = NULL,
                        transform = NULL, integrated = NULL,
                        nresample_clusters = NULL, ...) {
  if (inherits(object, "projection") ||
      (length(object) > 0 && inherits(object[[1]], "projection"))) {
    if (!is.null(filter_nterms)) {
      if (!.is_proj_list(object)) {
        object <- list(object)
      }
      projs <- Filter(function(x) {
        count_terms_chosen(x$solution_terms, add_icpt = TRUE) %in%
          (filter_nterms + 1)
      }, object)
    } else {
      projs <- object
    }
  } else {
    ## reference model or varsel object obtained, so run the projection
    projs <- project(object = object, ...)
  }

  if (!.is_proj_list(projs)) {
    projs <- list(projs)
  } else {
    ## projs is some other object, not containing an element called "family" (so
    ## it could be a `proj_list` but must not necessarily)
    if (any(sapply(projs, function(x) !("family" %in% names(x))))) {
      stop("Invalid object supplied to argument `object`.")
    }
  }

  if (is.null(newdata)) {
    ## pick first projection's function
    newdata <- projs[[1]]$refmodel$fetch_data()
    extract_y_ind <- TRUE
  } else {
    if (!inherits(newdata, c("matrix", "data.frame"))) {
      stop("newdata must be a data.frame or a matrix")
    }
    y_nm <- extract_terms_response(projs[[1]]$refmodel$formula)$response
    # Note: At this point, even for the binomial family with > 1 trials, we
    # expect only one response column name (the one for the successes), as
    # handled by get_refmodel.stanreg(), for example. Therefore, perform the
    # following check (needed for `extract_y_ind` later):
    stopifnot(length(y_nm) == 1)
    ### Might be helpful as a starting point in the future, but commented
    ### because some prediction functions might require only those columns from
    ### the original dataset which are needed for the corresponding submodel:
    # newdata_dummy <- projs[[1]]$refmodel$fetch_data()
    # if (is.data.frame(newdata) ||
    #     (is.matrix(newdata) && !is.null(colnames(newdata)))) {
    #   if (!setequal(setdiff(colnames(newdata), y_nm),
    #                 setdiff(colnames(newdata_dummy), y_nm))) {
    #     stop("`newdata` has to contain the same columns as the original ",
    #          "dataset (apart from ", paste(y_nm, collapse = ", "), ").")
    #   }
    # } else {
    #   warning("It seems like `newdata` is a matrix without column names. ",
    #           "It is safer to provide column names.")
    # }
    ###
    extract_y_ind <- ifelse(y_nm %in% colnames(newdata), TRUE, FALSE)
  }

  names(projs) <- sapply(projs, function(proj) {
    count_terms_chosen(proj$solution_terms, add_icpt = TRUE)
  })

  solution_terms <- list(...)$solution_terms
  if (!is.null(solution_terms) &&
      length(solution_terms) > NCOL(newdata)) {
    stop("The number of solution terms is greater than the number of columns ",
         "in `newdata`.")
  }

  preds <- lapply(projs, function(proj) {
    w_o <- proj$extract_model_data(proj$refmodel$fit, newdata = newdata,
                                   wrhs = weightsnew, orhs = offsetnew,
                                   extract_y = FALSE)
    weightsnew <- w_o$weights
    offsetnew <- w_o$offset
    if (length(weightsnew) == 0) {
      weightsnew <- rep(1, NROW(newdata))
    }
    if (length(offsetnew) == 0) {
      offsetnew <- rep(0, NROW(newdata))
    }
    mu <- proj$family$mu_fun(proj$sub_fit,
                             newdata = newdata, offset = offsetnew)

    onesub_fun(proj, mu, weightsnew,
               offset = offsetnew, newdata = newdata,
               extract_y_ind = extract_y_ind,
               transform = transform, integrated = integrated,
               nresample_clusters = nresample_clusters)
  })

  return(.unlist_proj(preds))
}

#' @rdname projection-linpred-predict
#' @export
proj_linpred <- function(object, newdata = NULL,
                         offsetnew = NULL, weightsnew = NULL,
                         filter_nterms = NULL,
                         transform = FALSE, integrated = FALSE, ...) {
  ## proj_helper lapplies fun to each projection in object
  proj_helper(
    object = object, newdata = newdata,
    offsetnew = offsetnew, weightsnew = weightsnew,
    onesub_fun = proj_linpred_aux, filter_nterms = filter_nterms,
    transform = transform, integrated = integrated, ...
  )
}

## function applied to each projected submodel in case of proj_linpred()
proj_linpred_aux <- function(proj, mu, weights, ...) {
  dot_args <- list(...)
  stopifnot(!is.null(dot_args$transform))
  stopifnot(!is.null(dot_args$integrated))
  stopifnot(!is.null(dot_args$newdata))
  stopifnot(!is.null(dot_args$offset))
  stopifnot(!is.null(dot_args$extract_y_ind))
  w_o <- proj$extract_model_data(proj$refmodel$fit, newdata = dot_args$newdata,
                                 wrhs = weights, orhs = dot_args$offset,
                                 extract_y = dot_args$extract_y_ind)
  ynew <- w_o$y
  lpd_out <- compute_lpd(
    ynew = ynew, mu = mu, proj = proj, weights = weights
  )
  pred_out <- if (!dot_args$transform) proj$family$linkfun(mu) else mu
  if (dot_args$integrated) {
    ## average over the posterior draws
    pred_out <- pred_out %*% proj$weights
    if (!is.null(lpd_out)) {
      lpd_out <- as.matrix(
        apply(lpd_out, 1, log_weighted_mean_exp, proj$weights)
      )
    }
  }
  return(nlist(pred = t(pred_out),
               lpd = if (is.null(lpd_out)) lpd_out else t(lpd_out)))
}

compute_lpd <- function(ynew, mu, proj, weights) {
  if (!is.null(ynew)) {
    ## compute also the log-density
    target <- .get_standard_y(ynew, weights, proj$family)
    ynew <- target$y
    weights <- target$weights
    return(as.matrix(proj$family$ll_fun(mu, proj$dis, ynew, weights)))
  } else {
    return(NULL)
  }
}

#' @rdname projection-linpred-predict
#' @export
proj_predict <- function(object, newdata = NULL,
                         offsetnew = NULL, weightsnew = NULL,
                         filter_nterms = NULL,
                         nresample_clusters = 1000, .seed = NULL, ...) {
  ## set random seed but ensure the old RNG state is restored on exit
  rng_state_old <- .Random.seed
  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  set.seed(.seed)

  ## proj_helper lapplies fun to each projection in object
  proj_helper(
    object = object, newdata = newdata,
    offsetnew = offsetnew, weightsnew = weightsnew,
    onesub_fun = proj_predict_aux, filter_nterms = filter_nterms,
    nresample_clusters = nresample_clusters, ...
  )
}

## function applied to each projected submodel in case of proj_predict()
proj_predict_aux <- function(proj, mu, weights, ...) {
  dot_args <- list(...)
  if (proj$p_type) {
    # In this case, the posterior draws have been clustered.
    stopifnot(!is.null(dot_args$nresample_clusters))
    draw_inds <- sample(
      x = seq_along(proj$weights), size = dot_args$nresample_clusters,
      replace = TRUE, prob = proj$weights
    )
  } else {
    draw_inds <- seq_along(proj$weights)
  }

  do.call(rbind, lapply(draw_inds, function(i) {
    proj$family$ppd(mu[, i], proj$dis[i], weights)
  }))
}

#' Plot summary statistics of a variable selection
#'
#' This is the [plot()] method for `vsel` objects (returned by [varsel()] or
#' [cv_varsel()]).
#'
#' @inheritParams summary.vsel
#' @param x An object of class `vsel` (returned by [varsel()] or [cv_varsel()]).
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
#'   # Variable selection (here without cross-validation and with small values
#'   # for `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the
#'   # sake of speed in this example; this is not recommended in general):
#'   vs <- varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                seed = 5555)
#'   print(plot(vs))
#' }
#'
#' @export
plot.vsel <- function(x, nterms_max = NULL, stats = "elpd",
                      deltas = FALSE, alpha = 0.32, baseline = NULL,
                      ...) {
  object <- x
  .validate_vsel_object_stats(object, stats)
  baseline <- .validate_baseline(object$refmodel, baseline, deltas)

  ## compute all the statistics and fetch only those that were asked
  nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1])
  tab <- rbind(
    .tabulate_stats(object, stats,
                    alpha = alpha,
                    nfeat_baseline = nfeat_baseline
    ),
    .tabulate_stats(object, stats, alpha = alpha)
  )
  stats_table <- subset(tab, tab$delta == deltas)
  stats_ref <- subset(stats_table, stats_table$size == Inf)
  stats_sub <- subset(stats_table, stats_table$size != Inf)
  stats_bs <- subset(stats_table, stats_table$size == nfeat_baseline)


  if (NROW(stats_sub) == 0) {
    stop(paste0(
      ifelse(length(stats) == 1, "Statistics ", "Statistic "),
      paste0(unique(stats), collapse = ", "), " not available."
    ))
  }

  max_size <- max(stats_sub$size)
  if (max_size == 0) {
    stop("plot.vsel() cannot be used if there is just the intercept-only ",
         "submodel.")
  }
  if (is.null(nterms_max)) {
    nterms_max <- max_size
  } else {
    # don't exceed the maximum submodel size
    nterms_max <- min(nterms_max, max_size)
  }
  if (nterms_max < 1) {
    stop("nterms_max must be at least 1")
  }
  ylab <- if (deltas) "Difference to the baseline" else "Value"

  # make sure that breaks on the x-axis are integers
  n_opts <- c(4, 5, 6)
  n_possible <- Filter(function(x) nterms_max %% x == 0, n_opts)
  n_alt <- n_opts[which.min(n_opts - (nterms_max %% n_opts))]
  nb <- ifelse(length(n_possible) > 0, min(n_possible), n_alt)
  by <- ceiling(nterms_max / min(nterms_max, nb))
  breaks <- seq(0, by * min(nterms_max, nb), by)
  minor_breaks <- if (by %% 2 == 0) {
    seq(by / 2, by * min(nterms_max, nb), by)
  } else {
    NULL
  }

  # plot submodel results
  pp <- ggplot(
    data = subset(stats_sub, stats_sub$size <= nterms_max),
    mapping = aes_string(x = "size")
  ) +
    geom_linerange(aes_string(ymin = "lq", ymax = "uq", alpha = 0.1)) +
    geom_line(aes_string(y = "value")) +
    geom_point(aes_string(y = "value"))

  if (!all(is.na(stats_ref$se))) {
    # add reference model results if they exist
    pp <- pp + geom_hline(aes_string(yintercept = "value"),
                          data = stats_ref,
                          color = "darkred", linetype = 2
    )
  }
  if (baseline != "ref") {
    # add the baseline result (if different from the reference model)
    pp <- pp + geom_hline(aes_string(yintercept = "value"),
                          data = stats_bs,
                          color = "black", linetype = 3
    )
  }
  pp <- pp +
    scale_x_continuous(
      breaks = breaks, minor_breaks = minor_breaks,
      limits = c(min(breaks), max(breaks))
    ) +
    labs(x = "Number of terms in the submodel", y = ylab) +
    theme(legend.position = "none") +
    facet_grid(statistic ~ ., scales = "free_y")
  return(pp)
}

#' Summary statistics of a variable selection
#'
#' This is the [summary()] method for `vsel` objects (returned by [varsel()] or
#' [cv_varsel()]).
#'
#' @param object An object of class `vsel` (returned by [varsel()] or
#'   [cv_varsel()]).
#' @param nterms_max Maximum submodel size for which the statistics are
#'   calculated. Note that `nterms_max` does not count the intercept, so use
#'   `nterms_max = 0` for the intercept-only model. For [plot.vsel()],
#'   `nterms_max` must be at least `1`.
#' @param stats One or several strings determining which statistics to
#'   calculate. Available statistics are:
#'   * `elpd`: (expected) sum of log predictive densities;
#'   * `mlpd`: mean log predictive density, that is, `elpd` divided by the
#'   number of observations (data points);
#'   * `mse`: mean squared error ([gaussian()] family only);
#'   * `rmse`: root mean squared error ([gaussian()] family only);
#'   * `acc`/`pctcorr`: classification accuracy ([binomial()] family only);
#'   * `auc`: area under the ROC curve ([binomial()] family only).
#' @param type One or more items from `"mean"`, `"se"`, `"lower"`, `"upper"`,
#'   `"diff"`, and `"diff.se"` indicating which of these to compute (mean,
#'   standard error, lower and upper credible bounds, difference to the
#'   corresponding statistic of the reference model, and standard error of this
#'   difference). The credible bounds are determined so that `1 - alpha` percent
#'   of the probability mass falls between them. Items `"diff"` and `"diff.se"`
#'   are only supported if `!deltas`.
#' @param deltas If `TRUE`, the submodel statistics are estimated relative to
#'   the baseline model (see argument `baseline`) instead of estimating the
#'   actual values of the statistics.
#' @param alpha A number giving the desired coverage of the credible intervals.
#'   For example, `alpha = 0.32` corresponds to 68% probability mass within the
#'   intervals, that is, one-standard-error intervals.
#' @param baseline Either `"ref"` or `"best"` indicating whether the baseline is
#'   the reference model or the best submodel found, respectively. If `NULL`,
#'   then `"ref"` is used, except for `datafit`s for which `"best"` is used.
#' @param ... Currently ignored.
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
#'   # Variable selection (here without cross-validation and with small values
#'   # for `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the
#'   # sake of speed in this example; this is not recommended in general):
#'   vs <- varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                seed = 5555)
#'   print(summary(vs))
#' }
#'
#' @export
summary.vsel <- function(object, nterms_max = NULL, stats = "elpd",
                         type = c("mean", "se", "diff", "diff.se"),
                         deltas = FALSE, alpha = 0.32, baseline = NULL, ...) {
  .validate_vsel_object_stats(object, stats)
  baseline <- .validate_baseline(object$refmodel, baseline, deltas)

  out <- list(
    formula = object$refmodel$formula,
    fit = object$fit,
    family = object$family,
    nobs = NROW(object$refmodel$fetch_data()),
    method = object$method,
    cv_method = object$cv_method,
    validate_search = object$validate_search,
    ndraws = object$ndraws,
    ndraws_pred = object$ndraws_pred,
    nclusters = object$nclusters,
    nclusters_pred = object$nclusters_pred
  )

  if (!is.null(out$validate_search)) {
    if (out$validate_search == TRUE) {
      out$search_included <- "search included"
    } else {
      out$search_included <- "search not included"
    }
  } else {
    out$search_included <- "search not included"
  }

  class(out) <- "vselsummary"
  ## fetch statistics
  if (deltas) {
    nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1])
    tab <- .tabulate_stats(object, stats,
                           alpha = alpha, nfeat_baseline = nfeat_baseline
    )
  } else {
    tab <- .tabulate_stats(object, stats, alpha = alpha)
  }
  stats_table <- subset(tab, tab$size != Inf) %>%
    dplyr::group_by(.data$statistic) %>%
    dplyr::slice_head(n = length(object$solution_terms) + 1)

  if (deltas) {
    type <- setdiff(type, c("diff", "diff.se"))
  }
  ## these are the corresponding names for mean, se, upper and lower in the
  ## stats_table, and their suffices in the table to be returned
  qty <- unname(sapply(type, function(t) {
    switch(t, mean = "value", upper = "uq", lower = "lq", se = "se",
           diff = "diff", diff.se = "diff.se")
  }))
  if (!is.null(object$cv_method)) {
    cv_suffix <- unname(switch(object$cv_method,
                               LOO = ".loo", kfold = ".kfold"
    ))
  } else {
    cv_suffix <- NULL
  }

  if (length(stats) > 1) {
    suffix <- lapply(stats, function(s) {
      unname(sapply(type, function(t) {
        paste0(
          s,
          switch(t, mean = cv_suffix, upper = ".upper", lower = ".lower",
                 se = ".se", diff = ".diff", diff.se = ".diff.se")
        )
      }))
    })
  } else {
    suffix <- list(unname(sapply(type, function(t) {
      switch(t, mean = paste0(stats, cv_suffix), upper = "upper",
             lower = "lower", se = "se",
             diff = "diff", diff.se = "diff.se"
      )
    })))
  }

  ## loop through all the required statistics
  arr <- data.frame(
    size = unique(stats_table$size),
    solution_terms = c(NA, object$solution_terms)
  )
  for (i in seq_along(stats)) {
    temp <- subset(stats_table, stats_table$statistic == stats[i], qty)
    newnames <- suffix[[i]]
    colnames(temp) <- newnames
    arr <- cbind(arr, temp)
  }

  if (is.null(nterms_max)) {
    nterms_max <- max(stats_table$size)
  }

  out$nterms <- nterms_max
  if ("pct_solution_terms_cv" %in% names(object)) {
    out$pct_solution_terms_cv <- object$pct_solution_terms_cv
  }

  out$suggested_size <- object$suggested_size
  out$selection <- subset(arr, arr$size <= nterms_max)
  return(out)
}

#' Print summary of variable selection
#'
#' This is the [print()] method for summary objects created by [summary.vsel()].
#' It displays a summary of the results of the projection predictive variable
#' selection.
#'
#' @param x An object of class `vselsummary`.
#' @param digits Number of decimal places to be reported.
#' @param ... Currently ignored.
#'
#' @return The output of [summary.vsel()] (invisible).
#'
#' @export
print.vselsummary <- function(x, digits = 1, ...) {
  print(x$family)
  cat("Formula: ")
  print(x$formula, showEnv = FALSE)
  cat(paste0("Observations: ", x$nobs, "\n"))
  if (!is.null(x$cv_method)) {
    cat(paste("CV method:", x$cv_method, x$search_included, "\n"))
  }
  cat(paste0("Search method: ", x$method, ", maximum number of terms ",
             max(x$selection$size), "\n"))
  if (!is.null(x$nclusters)) {
    cat(paste0(
      "Number of clusters used for selection: ", x$nclusters, "\n"
    ))
  } else {
    cat(paste0(
      "Number of draws used for selection: ", x$ndraws, "\n"
    ))
  }
  if (!is.null(x$nclusters_pred)) {
    cat(paste0(
      "Number of clusters used for prediction: ", x$nclusters_pred, "\n"
    ))
  } else {
    cat(paste0(
      "Number of draws used for prediction: ", x$ndraws_pred, "\n"
    ))
  }
  cat(paste0("Suggested Projection Size: ", x$suggested_size, "\n"))
  cat("\n")
  cat("Selection Summary:\n")
  print(
    x$selection %>% dplyr::mutate(dplyr::across(
      where(is.numeric),
      ~ round(., digits)
    )),
    row.names = FALSE
  )
  return(invisible(x))
}

#' Print results (summary) of variable selection
#'
#' This is the [print()] method for `vsel` objects (returned by [varsel()] or
#' [cv_varsel()]). It displays a summary of the results of the projection
#' predictive variable selection by first calling [summary.vsel()] and then
#' [print.vselsummary()].
#'
#' @param x An object of class `vsel` (returned by [varsel()] or [cv_varsel()]).
#' @param ... Further arguments passed to [summary.vsel()] (apart from
#'   argument `digits` which is passed to [print.vselsummary()]).
#'
#' @return The `data.frame` returned by [summary.vsel()] (invisible).
#'
#' @export
print.vsel <- function(x, ...) {
  stats <- summary.vsel(x, ...)
  print(stats, ...)
  return(invisible(stats))
}

#' Suggest model size
#'
#' This function can suggest an appropriate model size based on a decision rule
#' described in section "Details". Note that the decision rule is heuristic and
#' should only be interpreted as a guideline. It is recommended to examine the
#' results via [plot.vsel()] and/or [summary.vsel()] and make the final decision
#' based on what is most appropriate for the given problem.
#'
#' @inheritParams summary.vsel
#' @param object An object of class `vsel` (returned by [varsel()] or
#'   [cv_varsel()]).
#' @param stat Statistic used for the decision. See [summary.vsel()] for
#'   possible choices.
#' @param alpha A number giving the desired coverage of the credible intervals
#'   based on which the decision is made. For example, `alpha = 0.32`
#'   corresponds to 68% probability mass within the intervals, that is,
#'   one-standard-error intervals. See section "Details" for more information.
#' @param pct A number giving the relative proportion (*not* percents) between
#'   baseline model and null model utilities one is willing to sacrifice. See
#'   section "Details" for more information.
#' @param type Either `"upper"` or `"lower"` determining whether the decisions
#'   are based on the upper or lower credible bounds, respectively. See section
#'   "Details" for more information.
#' @param warnings Whether to give warnings if automatic suggestion fails,
#'   mainly for internal use. Usually there is no reason to set this to `FALSE`.
#' @param ... Currently ignored.
#'
#' @details The suggested model size is the smallest model size for which either
#'   the lower or upper bound (depending on argument `type`) of the credible
#'   interval (with coverage `1 - alpha`) for \eqn{u_k - u_{\mbox{base}}}{u_k -
#'   u_base} (with \eqn{u_k} denoting the \eqn{k}-th submodel's utility and
#'   \eqn{u_{\mbox{base}}}{u_base} denoting the baseline model's utility) falls
#'   above (or is equal to) \deqn{\mbox{pct} * (u_0 - u_{\mbox{base}})}{pct *
#'   (u_0 - u_base)} where \eqn{u_0} denotes the null model utility. The
#'   baseline is either the reference model or the best submodel found (see
#'   argument `baseline`).
#'
#'   For example, `alpha = 0.32`, `pct = 0`, and `type = "upper"` means that we
#'   select the smallest model size for which the upper bound of the
#'   corresponding credible interval exceeds (or is equal to) the baseline model
#'   utility, that is, which is better than the baseline model with a
#'   probability of at least 0.16 (and consequently, worse with a probability of
#'   at most 0.84). In other words, the estimated difference between the
#'   baseline model utility and the submodel utility is at most one standard
#'   error away from zero, so the two utilities are considered to be close.
#'
#' @note Loss statistics like the root mean-squared error (RMSE) and the
#'   mean-squared error (MSE) are converted to utilities by multiplying them by
#'   `-1`, so a call such as `suggest_size(object, stat = "rmse", type =
#'   "upper")` finds the smallest model size whose upper credible bound of the
#'   *negative* RMSE or MSE exceeds the cutoff (or equivalently has the lower
#'   credible bound of RMSE or MSE below the cutoff). This is done to make the
#'   interpretation of argument `type` the same regardless of argument `stat`.
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
#'   # Variable selection (here without cross-validation and with small values
#'   # for `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the
#'   # sake of speed in this example; this is not recommended in general):
#'   vs <- varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                seed = 5555)
#'   print(suggest_size(vs))
#' }
#'
#' @export
suggest_size <- function(object, ...) {
  UseMethod("suggest_size")
}

#' @rdname suggest_size
#' @export
suggest_size.vsel <- function(object, stat = "elpd", alpha = 0.32, pct = 0,
                              type = "upper", baseline = NULL, warnings = TRUE,
                              ...) {
  .validate_vsel_object_stats(object, stat)
  if (length(stat) > 1) {
    stop("Only one statistic can be specified to suggest_size")
  }

  if (.is_util(stat)) {
    sgn <- 1
  } else {
    sgn <- -1
    if (type == "upper") {
      type <- "lower"
    } else {
      type <- "upper"
    }
  }
  if (!is.null(object$cv_method)) {
    suffix <- paste0(".", tolower(object$cv_method))
  } else {
    suffix <- ""
  }
  bound <- type
  stats <- summary.vsel(object,
                        stats = stat, alpha = alpha,
                        type = c("mean", "upper", "lower"),
                        baseline = baseline, deltas = TRUE
  )$selection
  util_null <- sgn * unlist(unname(subset(
    stats, stats$size == 0,
    paste0(stat, suffix)
  )))
  util_cutoff <- pct * util_null
  res <- subset(stats, sgn * stats[, bound] >= util_cutoff, "size")

  if (nrow(res) == 0) {
    ## no submodel satisfying the criterion found
    if (object$nterms_max == object$nterms_all) {
      suggested_size <- object$nterms_max
    } else {
      suggested_size <- NA
      if (warnings) {
        warning(paste(
          "Could not suggest model size. Investigate plot.vsel to identify",
          "if the search was terminated too early. If this is the case,",
          "run variable selection with larger value for nterms_max."
        ))
      }
    }
  } else {
    suggested_size <- min(res) + 1
  }

  return(suggested_size - 1) ## substract the intercept
}

replace_intercept_name <- function(names) {
  return(gsub(
    "\\(Intercept\\)",
    "Intercept",
    names
  ))
}

replace_population_names <- function(population_effects) {
  # Use brms's naming convention:
  names(population_effects) <- replace_intercept_name(names(population_effects))
  names(population_effects) <- paste0("b_", names(population_effects))
  return(population_effects)
}

#' @keywords internal
#' @export
coef.subfit <- function(object, ...) {
  return(with(object, c(
    "Intercept" = alpha,
    setNames(beta, colnames(x))
  )))
}

#' @method as.matrix lm
#' @keywords internal
#' @export
as.matrix.lm <- function(x, ...) {
  return(coef(x) %>%
           replace_population_names())
}

#' @method as.matrix ridgelm
#' @keywords internal
#' @export
as.matrix.ridgelm <- function(x, ...) {
  return(as.matrix.lm(x))
}

#' @method as.matrix subfit
#' @keywords internal
#' @export
as.matrix.subfit <- function(x, ...) {
  return(as.matrix.lm(x, ...))
}

#' @method as.matrix glm
#' @keywords internal
#' @export
as.matrix.glm <- function(x, ...) {
  return(as.matrix.lm(x, ...))
}

#' @method as.matrix lmerMod
#' @keywords internal
#' @export
as.matrix.lmerMod <- function(x, ...) {
  population_effects <- lme4::fixef(x) %>%
    replace_population_names()

  # Extract variance components:
  group_vc_raw <- lme4::VarCorr(x)
  group_vc <- unlist(lapply(group_vc_raw, function(vc_obj) {
    # The vector of standard deviations:
    vc_out <- c("sd" = attr(vc_obj, "stddev"))
    # The correlation matrix:
    cor_mat <- attr(vc_obj, "correlation")
    if (!is.null(cor_mat)) {
      # Auxiliary object: A matrix of the same dimension as cor_mat, but
      # containing the paste()-d dimnames:
      cor_mat_nms <- matrix(apply(expand.grid(
        rownames(cor_mat),
        colnames(cor_mat)
      ),
      1, paste,
      collapse = "."
      ),
      nrow = nrow(cor_mat), ncol = ncol(cor_mat)
      )
      # Note: With upper.tri() (and also with lower.tri()), the indexed matrix
      # is coerced to a vector in column-major order:
      vc_out <- c(
        vc_out,
        "cor" = setNames(
          cor_mat[upper.tri(cor_mat)],
          cor_mat_nms[upper.tri(cor_mat_nms)]
        )
      )
    }
    return(vc_out)
  }))

  # Use brms's naming convention:
  names(group_vc) <- replace_intercept_name(names(group_vc))

  # We will have to move the substrings "sd\\." and "cor\\." up front (i.e. in
  # front of the group name), so make sure that they don't occur in the group
  # names:
  stopifnot(!any(grepl("sd\\.|cor\\.", names(group_vc_raw))))
  # Move the substrings "sd\\." and "cor\\." up front and replace the dot
  # following the group name by double underscores:
  names(group_vc) <- sub(
    paste0(
      "(",
      paste(
        gsub("\\.", "\\\\.", names(group_vc_raw)),
        collapse = "|"
      ),
      ")\\.(sd|cor)\\."
    ),
    "\\2_\\1__",
    names(group_vc)
  )
  # Replace dots between coefficient names by double underscores:
  coef_nms <- lapply(group_vc_raw, rownames)
  for (coef_nms_i in coef_nms) {
    coef_nms_i <- replace_intercept_name(coef_nms_i)
    names(group_vc) <- gsub(
      paste0(
        "(",
        paste(
          gsub("\\.", "\\\\.", coef_nms_i),
          collapse = "|"
        ),
        ")\\."
      ),
      "\\1__",
      names(group_vc)
    )
  }

  # Extract the group-level effects themselves:
  group_ef <- unlist(lapply(lme4::ranef(x), function(ranef_df) {
    ranef_mat <- as.matrix(ranef_df)
    setNames(
      as.vector(ranef_mat),
      apply(
        expand.grid(rownames(ranef_mat), colnames(ranef_mat)),
        1, function(row_col_nm) {
          paste(rev(row_col_nm), collapse = ".")
        }
      )
    )
  }))

  # Use brms's naming convention:
  names(group_ef) <- replace_intercept_name(names(group_ef))
  names(group_ef) <- paste0("r_", names(group_ef))
  for (coef_nms_idx in seq_along(coef_nms)) {
    group_nm_i <- names(coef_nms)[coef_nms_idx]
    coef_nms_i <- coef_nms[[coef_nms_idx]]
    coef_nms_i <- replace_intercept_name(coef_nms_i)
    # Put the part following the group name in square brackets, reorder its two
    # subparts (coefficient name and group level) and separate them by comma:
    names(group_ef) <- sub(
      paste0(
        "(",
        gsub("\\.", "\\\\.", group_nm_i),
        ")\\.(",
        paste(
          gsub("\\.", "\\\\.", coef_nms_i),
          collapse = "|"
        ),
        ")\\.(.*)$"
      ),
      "\\1[\\3,\\2]",
      names(group_ef)
    )
  }

  return(c(population_effects, group_vc, group_ef))
}

#' @method as.matrix glmerMod
#' @keywords internal
#' @export
as.matrix.glmerMod <- function(x, ...) {
  return(as.matrix.lmerMod(x, ...))
}

#' @method as.matrix gamm4
#' @keywords internal
#' @export
as.matrix.gamm4 <- function(x, ...) {
  return(as.matrix.lm(x, ...))
}

#' @method as.matrix list
#' @keywords internal
#' @export
as.matrix.list <- function(x, ...) {
  return(do.call(cbind, lapply(x, as.matrix.glm, ...)))
}

#' @keywords internal
#' @export
t.glm <- function(x, ...) {
  return(t(as.matrix(x), ...))
}

#' @keywords internal
#' @export
t.lm <- function(x, ...) {
  return(t(as.matrix(x), ...))
}

#' @keywords internal
#' @export
t.ridgelm <- function(x, ...) {
  return(t(as.matrix(x), ...))
}

#' @keywords internal
#' @export
t.list <- function(x, ...) {
  return(t(as.matrix(x), ...))
}

#' Extract projected parameter draws
#'
#' This is the [as.matrix()] method for `projection` objects (returned by
#' [project()], possibly as elements of a `list`). It extracts the projected
#' parameter draws and returns them as a matrix.
#'
#' @param x An object of class `projection` (returned by [project()], possibly
#'   as elements of a `list`).
#' @param ... Currently ignored.
#'
#' @return An \eqn{S_{\mbox{prj}} \times Q}{S_prj x Q} matrix of projected
#'   draws, with \eqn{S_{\mbox{prj}}}{S_prj} denoting the number of projected
#'   draws and \eqn{Q} the number of parameters.
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
#'   # Projection onto an arbitrary combination of predictor terms (with a small
#'   # value for `nclusters`, but only for the sake of speed in this example;
#'   # this is not recommended in general):
#'   prj <- project(fit, solution_terms = c("X1", "X3", "X5"), nclusters = 10,
#'                  seed = 9182)
#'   prjmat <- as.matrix(prj)
#'   ### For further post-processing (e.g., via packages `bayesplot` and
#'   ### `posterior`), we will here ignore the fact that clustering was used
#'   ### (due to argument `nclusters` above). CAUTION: Ignoring the clustering
#'   ### is not recommended and only shown here for demonstration purposes. A
#'   ### better solution for the clustering case is explained below.
#'   # If the `bayesplot` package is installed, the output from
#'   # as.matrix.projection() can be used there, for example:
#'   if (requireNamespace("bayesplot", quietly = TRUE)) {
#'     print(bayesplot::mcmc_intervals(prjmat))
#'   }
#'   # If the `posterior` package is installed, the output from
#'   # as.matrix.projection() can be used there, for example:
#'   if (requireNamespace("posterior", quietly = TRUE)) {
#'     prjdrws <- posterior::as_draws_matrix(prjmat)
#'     print(posterior::summarize_draws(
#'       prjdrws,
#'       "median", "mad", function(x) quantile(x, probs = c(0.025, 0.975))
#'     ))
#'   }
#'   ### Better solution for post-processing clustered draws (e.g., via
#'   ### `bayesplot` or `posterior`): Don't ignore the fact that clustering was
#'   ### used. Instead, resample the clusters according to their weights (e.g.,
#'   ### via posterior::resample_draws()). However, this requires access to the
#'   ### cluster weights which is not implemented in `projpred` yet. This
#'   ### example will be extended as soon as those weights are accessible.
#' }
#'
#' @method as.matrix projection
#' @export
as.matrix.projection <- function(x, ...) {
  if (inherits(x$refmodel, "datafit")) {
    stop("as.matrix.projection() does not work for objects based on ",
         "\"datafit\"s.")
  }
  if (x$p_type) {
    warning(paste(
      "Note that projection was performed using",
      "clustering and the clusters might have different weights."
    ))
  }
  if (!all(sapply(x$sub_fit, inherits, what = get_as.matrix_cls_projpred()))) {
    # Throw an error because in this case, we probably need a new
    # as.matrix.<class_name>() method.
    stop("This case should not occur. Please notify the package maintainer.")
  }
  res <- t(do.call(cbind, lapply(x$sub_fit, as.matrix)))
  if (x$family$family == "gaussian") res <- cbind(res, sigma = x$dis)
  return(res)
}

#' Create cross-validation folds
#'
#' These are helper functions to create cross-validation (CV) folds, i.e., to
#' split up the indices from 1 to `n` into `K` subsets ("folds") for K-fold CV.
#' These functions are potentially useful when creating the `cvfits` and `cvfun`
#' arguments for [init_refmodel()]. The return value is different for these two
#' methods, see below for details.
#'
#' @name cv-indices
#'
#' @param n Number of observations (data points).
#' @param K Number of folds. Must be at least 2 and not exceed `n`.
#' @param out Format of the output, either `"foldwise"` or `"indices"`. See
#'   below for details.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. If `NULL`, no seed is set and
#'   therefore, the results are not reproducible. See [set.seed()] for details.
#'
#' @return [cvfolds()] returns a vector of length `n` such that each element is
#'   an integer between 1 and `k` denoting which fold the corresponding data
#'   point belongs to. The return value of [cv_ids()] depends on the `out`
#'   argument. If `out = "foldwise"`, the return value is a `list` with `k`
#'   elements, each being a `list` with elements `tr` and `ts` giving the
#'   training and test indices, respectively, for the corresponding fold. If
#'   `out = "indices"`, the return value is a `list` with elements `tr` and `ts`
#'   each being a `list` with `k` elements giving the training and test indices,
#'   respectively, for each fold.
#'
#' @examples
#' n <- 100
#' set.seed(1234)
#' y <- rnorm(n)
#' cv <- cv_ids(n, K = 5, seed = 9876)
#' # Mean within the test set of each fold:
#' cvmeans <- sapply(cv, function(fold) mean(y[fold$ts]))
#'
NULL

#' @rdname cv-indices
#' @export
cvfolds <- function(n, K, seed = NULL) {
  .validate_num_folds(K, n)

  ## set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  ## create and shuffle the indices
  folds <- rep_len(seq_len(K), length.out = n)
  folds <- sample(folds, n, replace = FALSE)

  return(folds)
}

#' @rdname cv-indices
#' @export
cv_ids <- function(n, K, out = c("foldwise", "indices"), seed = NULL) {
  .validate_num_folds(K, n)
  out <- match.arg(out)

  # set random seed but ensure the old RNG state is restored on exit
  if (exists(".Random.seed")) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)

  # shuffle the indices
  ind <- sample(seq_len(n), n, replace = FALSE)

  if (out == "foldwise") {
    cv <- lapply(seq_len(K), function(i) {
      ts <- sort(ind[seq(i, n, K)]) # test set
      tr <- setdiff(seq_len(n), ts) # training set
      list(tr = tr, ts = ts)
    })
  } else if (out == "indices") {
    cv <- list()
    cv$tr <- list()
    cv$ts <- list()
    for (i in seq_len(K)) {
      ts <- sort(ind[seq(i, n, K)]) # test set
      tr <- setdiff(seq_len(n), ts) # training set
      cv$tr[[i]] <- tr
      cv$ts[[i]] <- ts
    }
  }

  return(cv)
}

#' Retrieve predictor solution path or predictor combination
#'
#' This function retrieves the "solution terms" from an object. For `vsel`
#' objects (returned by [varsel()] or [cv_varsel()]), this is the predictor
#' solution path of the variable selection. For `projection` objects (returned
#' by [project()], possibly as elements of a `list`), this is the predictor
#' combination upon which the projection was performed.
#'
#' @param object The object from which to retrieve the solution terms. Possible
#'   classes may be inferred from the names of the corresponding methods (see
#'   also the description).
#' @param ... Currently ignored.
#'
#' @return A character vector of solution terms.
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
#'   # Variable selection (here without cross-validation and with small values
#'   # for `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the
#'   # sake of speed in this example; this is not recommended in general):
#'   vs <- varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                seed = 5555)
#'   print(solution_terms(vs))
#'
#'   # Projection onto an arbitrary combination of predictor terms (with a small
#'   # value for `nclusters`, but only for the sake of speed in this example;
#'   # this is not recommended in general):
#'   prj <- project(fit, solution_terms = c("X1", "X3", "X5"), nclusters = 10,
#'                  seed = 9182)
#'   print(solution_terms(prj)) # gives `c("X1", "X3", "X5")`
#' }
#'
#' @export
solution_terms <- function(object, ...) {
  UseMethod("solution_terms")
}

#' @rdname solution_terms
#' @export
solution_terms.vsel <- function(object, ...) {
  return(object$solution_terms)
}

#' @rdname solution_terms
#' @export
solution_terms.projection <- function(object, ...) {
  return(object$solution_terms)
}
