#' Predictions from a submodel (after projection)
#'
#' After the projection of the reference model onto a submodel, the linear
#' predictors (for the original dataset or new data) based on that submodel can
#' be calculated by [proj_linpred()]. The linear predictors can also be
#' transformed to response scale. Furthermore, [proj_linpred()] returns the
#' corresponding log predictive density values if the new dataset contains
#' response values. The [proj_predict()] function draws from the predictive
#' distribution of the submodel that the reference model has been projected
#' onto. If the projection has not been performed yet, both functions call
#' [project()] internally to perform the projection. Both functions can also
#' handle multiple submodels at once (for `object`s of class `vsel` or `object`s
#' returned by a [project()] call to an object of class `vsel`; see
#' [project()]).
#'
#' @name pred-projection
#'
#' @template args-newdata
#' @param object An object returned by [project()] or an object that can be
#'   passed to argument `object` of [project()].
#' @param filter_nterms Only applies if `object` is an object returned by
#'   [project()]. In that case, `filter_nterms` can be used to filter `object`
#'   for only those elements (submodels) with a number of solution terms in
#'   `filter_nterms`. Therefore, needs to be a numeric vector or `NULL`. If
#'   `NULL`, use all submodels.
#' @param transform For [proj_linpred()] only. A single logical value indicating
#'   whether the linear predictor should be transformed to response scale using
#'   the inverse-link function (`TRUE`) or not (`FALSE`). In case of the latent
#'   projection, argument `transform` is similar in spirit to argument
#'   `respOrig` from other functions and affects the scale of both output
#'   elements `pred` and `lpd` (see section "Value" below). In that case, note
#'   that `transform = FALSE` yields the linear predictors (in output element
#'   `pred`) without any modifications that may be due to the original response
#'   distribution (e.g., for a [brms::cumulative()] model, the ordered
#'   thresholds are not taken into account).
#' @param integrated For [proj_linpred()] only. A single logical value
#'   indicating whether the output should be averaged across the projected
#'   posterior draws (`TRUE`) or not (`FALSE`).
#' @param nresample_clusters For [proj_predict()] with clustered projection
#'   only. Number of draws to return from the predictive distribution of the
#'   submodel. Not to be confused with argument `nclusters` of [project()]:
#'   `nresample_clusters` gives the number of draws (*with* replacement) from
#'   the set of clustered posterior draws after projection (with this set being
#'   determined by argument `nclusters` of [project()]).
#' @param .seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. Here,
#'   this seed is used for drawing new group-level effects in case of a
#'   multilevel submodel (however, not yet in case of a GAMM) and for drawing
#'   from the predictive distribution of the submodel(s) in case of
#'   [proj_predict()]. If a clustered projection was performed, then in
#'   [proj_predict()], `.seed` is also used for drawing from the set of the
#'   projected clusters of posterior draws (see argument `nresample_clusters`).
#' @param respOrig Only relevant for the latent projection. A single logical
#'   value indicating whether to draw from the predictive distribution on the
#'   original response scale (`TRUE`) or on latent scale (`FALSE`).
#' @param ... Arguments passed to [project()] if `object` is not already an
#'   object returned by [project()].
#'
#' @return In the following, \eqn{S_{\mathrm{prj}}}{S_prj}, \eqn{N},
#'   \eqn{C_{\mathrm{cat}}}{C_cat}, and \eqn{C_{\mathrm{lat}}}{C_lat} from help
#'   topic [refmodel-init-get] are used. Furthermore, let \eqn{C} denote either
#'   \eqn{C_{\mathrm{cat}}}{C_cat} (if `transform = TRUE`) or
#'   \eqn{C_{\mathrm{lat}}}{C_lat} (if `transform = FALSE`). Then, if the
#'   prediction is done for one submodel only (i.e., `length(nterms) == 1 ||
#'   !is.null(solution_terms)` in the call to [project()]):
#'   * [proj_linpred()] returns a `list` with elements `pred` (predictions,
#'   i.e., the linear predictors, possibly transformed to response scale) and
#'   `lpd` (log predictive densities; only calculated if `newdata` contains
#'   response values). In case of (i) the traditional projection or (ii) the
#'   latent projection with `<refmodel>$family$cats` (where `<refmodel>` is an
#'   object resulting from [init_refmodel()]; see also [extend_family()]'s
#'   argument `latent_y_unqs`) being `NULL`, both elements are
#'   \eqn{S_{\mathrm{prj}} \times N}{S_prj x N} matrices. In case of (i) the
#'   augmented-data projection or (ii) the latent projection with `transform =
#'   TRUE` and `<refmodel>$family$cats` being not `NULL`, `pred` is an
#'   \eqn{S_{\mathrm{prj}} \times N \times C}{S_prj x N x C} array and `lpd` is
#'   an \eqn{S_{\mathrm{prj}} \times N}{S_prj x N} matrix.
#'   * [proj_predict()] returns an \eqn{S_{\mathrm{prj}} \times N}{S_prj x N}
#'   matrix of predictions where \eqn{S_{\mathrm{prj}}}{S_prj} denotes
#'   `nresample_clusters` in case of clustered projection. In case of (i) the
#'   augmented-data projection or (ii) the latent projection with `respOrig =
#'   TRUE` and `<refmodel>$family$cats` being not `NULL`, this matrix has an
#'   attribute called `cats` (the character vector of response categories) and
#'   the values of the matrix are the predicted indices of the response
#'   categories (these indices refer to the order of the response categories
#'   from attribute `cats`).
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
#'   # The "stanreg" fit which will be used as the reference model (with small
#'   # values for `chains` and `iter`, but only for technical reasons in this
#'   # example; this is not recommended in general):
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
#'   # Predictions (at the training points) from the submodel onto which the
#'   # reference model was projected:
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
proj_helper <- function(object, newdata, offsetnew, weightsnew, onesub_fun,
                        filter_nterms = NULL, ...) {
  if (inherits(object, "projection") || .is_proj_list(object)) {
    if (!is.null(filter_nterms)) {
      if (!.is_proj_list(object)) {
        object <- list(object)
      }
      projs <- Filter(
        function(x) {
          count_terms_chosen(x$solution_terms, add_icpt = TRUE) %in%
            (filter_nterms + 1)
        },
        object
      )
      if (!length(projs)) {
        stop("Invalid `filter_nterms`.")
      }
    } else {
      projs <- object
    }
  } else {
    ## reference model or varsel object obtained, so run the projection
    projs <- project(object = object, ...)
  }

  if (!.is_proj_list(projs)) {
    projs <- list(projs)
  }

  if (is.null(newdata)) {
    ## pick first projection's function
    newdata <- projs[[1]]$refmodel$fetch_data()
    extract_y_ind <- TRUE
  } else {
    if (!inherits(newdata, c("matrix", "data.frame"))) {
      stop("newdata must be a data.frame or a matrix")
    }
    newdata <- na.fail(newdata)
    y_nm <- extract_terms_response(projs[[1]]$refmodel$formula)$response
    # Note: At this point, even for the binomial family with > 1 trials, we
    # expect only one response column name (the one for the successes), as
    # handled by get_refmodel.stanreg(), for example. Therefore, perform the
    # following check (needed for `extract_y_ind` later):
    stopifnot(length(y_nm) == 1)
    if (projs[[1]]$refmodel$family$for_latent) {
      # Remove the leading dot which was added in init_refmodel():
      y_nm <- sub("^\\.", "", y_nm)
    }
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
    extract_y_ind <- y_nm %in% colnames(newdata)
  }

  names(projs) <- sapply(projs, function(proj) {
    count_terms_chosen(proj$solution_terms, add_icpt = TRUE)
  })

  preds <- lapply(projs, function(proj) {
    w_o <- proj$refmodel$extract_model_data(
      proj$refmodel$fit, newdata = newdata, wrhs = weightsnew, orhs = offsetnew,
      extract_y = FALSE
    )
    weightsnew <- w_o$weights
    offsetnew <- w_o$offset
    if (length(weightsnew) == 0) {
      weightsnew <- rep(1, NROW(newdata))
    }
    if (length(offsetnew) == 0) {
      offsetnew <- rep(0, NROW(newdata))
    }
    if (proj$refmodel$family$for_augdat && !all(weightsnew == 1)) {
      stop("Currently, the augmented-data projection may not be combined with ",
           "observation weights (other than 1).")
    }
    if (proj$refmodel$family$for_latent && !all(weightsnew == 1)) {
      stop("Currently, the latent projection may not be combined with ",
           "observation weights (other than 1).")
    }
    onesub_fun(proj, newdata = newdata, offset = offsetnew,
               weights = weightsnew, extract_y_ind = extract_y_ind, ...)
  })

  return(.unlist_proj(preds))
}

#' @rdname pred-projection
#' @export
proj_linpred <- function(object, newdata = NULL, offsetnew = NULL,
                         weightsnew = NULL, filter_nterms = NULL,
                         transform = FALSE, integrated = FALSE,
                         .seed = sample.int(.Machine$integer.max, 1), ...) {
  # Set seed, but ensure the old RNG state is restored on exit:
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  if (!is.na(.seed)) set.seed(.seed)

  ## proj_helper lapplies fun to each projection in object
  proj_helper(
    object = object, newdata = newdata,
    offsetnew = offsetnew, weightsnew = weightsnew,
    onesub_fun = proj_linpred_aux, filter_nterms = filter_nterms,
    transform = transform, integrated = integrated, ...
  )
}

## function applied to each projected submodel in case of proj_linpred()
proj_linpred_aux <- function(proj, newdata, offset, weights, transform = FALSE,
                             integrated = FALSE, extract_y_ind = TRUE, ...) {
  pred_sub <- proj$refmodel$family$mu_fun(proj$submodl, newdata = newdata,
                                          offset = offset,
                                          transform = transform)
  if (proj$refmodel$family$for_latent && transform) {
    pred_sub <- proj$refmodel$family$latent_ilink(
      t(pred_sub), cl_ref = proj$cl_ref, wdraws_ref = proj$wdraws_ref
    )
    if (length(dim(pred_sub)) < 2) {
      stop("Unexpected structure for the output of `latent_ilink`.")
    }
    if (all(is.na(pred_sub))) {
      message(
        "`latent_ilink` returned only `NA`s, so the corresponding output will ",
        "also be `NA` as long as `transform = TRUE`."
      )
    }
  }
  w_o <- proj$refmodel$extract_model_data(
    proj$refmodel$fit, newdata = newdata, wrhs = weights,
    orhs = offset, extract_y = extract_y_ind
  )
  ynew <- w_o$y
  if (!is.null(ynew) && proj$refmodel$family$for_latent && !transform) {
    # Use `ref_predfun_usr` here (instead of `ref_predfun`) to include
    # offsets:
    refprd_with_offs <- get("ref_predfun_usr",
                            envir = environment(proj$refmodel$ref_predfun))
    ynew <- rowMeans(unname(
      refprd_with_offs(fit = proj$refmodel$fit, newdata = newdata)
    ))
  }
  lpd_out <- compute_lpd(ynew = ynew, pred_sub = pred_sub, proj = proj,
                         weights = weights, transformed = transform)
  if (proj$refmodel$family$for_latent && transform) {
    if (length(dim(pred_sub)) == 3 && integrated) {
      pred_sub <- arr2augmat(pred_sub, margin_draws = 1)
    }
  }
  if (integrated) {
    ## average over the projected draws
    if (proj$refmodel$family$for_latent && transform &&
        !inherits(pred_sub, "augmat")) {
      pred_sub <- proj$weights %*% pred_sub
    } else {
      pred_sub <- structure(pred_sub %*% proj$weights,
                            nobs_orig = attr(pred_sub, "nobs_orig"),
                            class = oldClass(pred_sub))
    }
    if (!is.null(lpd_out)) {
      if (!proj$refmodel$family$for_latent) {
        marg_obs <- 1
      } else {
        marg_obs <- 2
      }
      lpd_out <- as.matrix(
        apply(lpd_out, marg_obs, log_weighted_mean_exp, proj$weights)
      )
    }
  }
  if (inherits(pred_sub, "augmat")) {
    pred_sub <- augmat2arr(pred_sub, margin_draws = 1)
  } else if (!(proj$refmodel$family$for_latent && transform)) {
    pred_sub <- t(pred_sub)
  }
  if (!is.null(lpd_out) &&
      (!proj$refmodel$family$for_latent ||
       (proj$refmodel$family$for_latent && integrated))) {
    lpd_out <- t(lpd_out)
  }
  return(nlist(pred = pred_sub, lpd = lpd_out))
}

compute_lpd <- function(ynew, pred_sub, proj, weights, transformed) {
  if (!is.null(ynew)) {
    ## compute also the log-density
    target <- .get_standard_y(ynew, weights, proj$refmodel$family)
    ynew <- target$y
    weights <- target$weights
    if ((!proj$refmodel$family$for_latent ||
         (proj$refmodel$family$for_latent && transformed)) &&
        !is.null(proj$refmodel$family$cats)) {
      ynew <- as.factor(ynew)
      if (!all(levels(ynew) %in% proj$refmodel$family$cats)) {
        if (proj$refmodel$family$for_augdat) {
          y_unqs_str <- "augdat_y_unqs"
        } else {
          y_unqs_str <- "latent_y_unqs"
        }
        stop("The levels of the response variable (after coercing it to a ",
             "`factor`) have to be a subset of `family$cats`. Either modify ",
             "`newdata` or the function supplied to `extract_model_data` in ",
             "init_refmodel() accordingly or see the documentation for ",
             "extend_family()'s argument `", y_unqs_str, "` to solve this.")
      }
      # Re-assign the original levels because some levels might be missing:
      ynew <- factor(ynew, levels = proj$refmodel$family$cats)
    } else if (proj$refmodel$family$for_latent && transformed &&
               is.null(proj$refmodel$family$cats) &&
               (is.factor(ynew) || is.character(ynew) || is.logical(ynew))) {
      stop("If the original (i.e., non-latent) response is `factor`-like, ",
           "`family$cats` must not be `NULL`. See the documentation for ",
           "extend_family()'s argument `latent_y_unqs` to solve this.")
    }
    if (!transformed) {
      pred_sub <- proj$refmodel$family$linkinv(pred_sub)
    }
    if (proj$refmodel$family$for_latent && transformed) {
      llOrig_out <- proj$refmodel$family$latent_llOrig(
        pred_sub, yOrig = ynew, wobs = weights, cl_ref = proj$cl_ref,
        wdraws_ref = proj$wdraws_ref
      )
      if (!is.matrix(llOrig_out)) {
        stop("Unexpected structure for the output of `latent_llOrig`.")
      }
      if (all(is.na(llOrig_out))) {
        message(
          "`latent_llOrig` returned only `NA`s, so the corresponding output ",
          "will also be `NA` as long as `transform = TRUE`."
        )
      }
      return(llOrig_out)
    } else {
      if (proj$refmodel$family$for_latent && all(is.na(proj$refmodel$dis))) {
        # There's already a corresponding message thrown at the time when the
        # reference model was built, but users might have forgotten about it, so
        # throw another one here:
        message(
          "Cannot calculate LPD values if `transform = FALSE` and ",
          "`<refmodel>$dis` consists of only `NA`s. You should have received ",
          "a message describing possible remedies when the reference model ",
          "was built."
        )
      }
      return(proj$refmodel$family$ll_fun(pred_sub, proj$dis, ynew, weights))
    }
  } else {
    return(NULL)
  }
}

#' @rdname pred-projection
#' @export
proj_predict <- function(object, newdata = NULL, offsetnew = NULL,
                         weightsnew = NULL, filter_nterms = NULL,
                         nresample_clusters = 1000,
                         .seed = sample.int(.Machine$integer.max, 1),
                         respOrig = TRUE, ...) {
  # Set seed, but ensure the old RNG state is restored on exit:
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  if (!is.na(.seed)) set.seed(.seed)

  ## proj_helper lapplies fun to each projection in object
  proj_helper(
    object = object, newdata = newdata,
    offsetnew = offsetnew, weightsnew = weightsnew,
    onesub_fun = proj_predict_aux, filter_nterms = filter_nterms,
    nresample_clusters = nresample_clusters, respOrig = respOrig, ...
  )
}

## function applied to each projected submodel in case of proj_predict()
proj_predict_aux <- function(proj, newdata, offset, weights,
                             nresample_clusters = 1000, respOrig = TRUE, ...) {
  if (!proj$refmodel$family$for_latent && !respOrig) {
    stop("`respOrig = FALSE` can only be used in case of the latent ",
         "projection.")
  }
  mu <- proj$refmodel$family$mu_fun(proj$submodl,
                                    newdata = newdata,
                                    offset = offset)
  if (proj$p_type) {
    # In this case, the posterior draws have been clustered.
    draw_inds <- sample(x = seq_along(proj$weights), size = nresample_clusters,
                        replace = TRUE, prob = proj$weights)
  } else {
    draw_inds <- seq_along(proj$weights)
  }
  cats_aug <- proj$refmodel$family$cats
  if (proj$refmodel$family$for_latent && respOrig) {
    mu_Orig <- proj$refmodel$family$latent_ilink(t(mu), cl_ref = proj$cl_ref,
                                                 wdraws_ref = proj$wdraws_ref)
    if (length(dim(mu_Orig)) < 2) {
      stop("Unexpected structure for the output of `latent_ilink`.")
    }
    if (length(dim(mu_Orig)) == 3) {
      mu_Orig_resamp <- mu_Orig[draw_inds, , , drop = FALSE]
    } else {
      mu_Orig_resamp <- mu_Orig[draw_inds, , drop = FALSE]
    }
    pppd_out <- proj$refmodel$family$latent_ppdOrig(
      mu_Orig_resamp, wobs = weights, cl_ref = proj$cl_ref,
      wdraws_ref = proj$wdraws_ref, idxs_prjdraws = draw_inds
    )
    if (!is.matrix(pppd_out)) {
      stop("Unexpected structure for the output of `latent_ppdOrig`.")
    }
    if (all(is.na(pppd_out))) {
      message(
        "`latent_ppdOrig` returned only `NA`s, so the output will also be ",
        "`NA` as long as `respOrig = TRUE`."
      )
    }
  } else {
    if (proj$refmodel$family$for_latent) {
      # In this case, the PPPD will be on latent scale, so the response-scale
      # categories should not be appended as an attribute to the output:
      if (!is.null(cats_aug)) {
        cats_aug <- NULL
      }
      if (all(is.na(proj$refmodel$dis))) {
        # There's already a corresponding message thrown at the time when the
        # reference model was built, but users might have forgotten about it, so
        # throw another one here:
        message(
          "Cannot draw from the latent Gaussian distribution if ",
          "`<refmodel>$dis` consists of only `NA`s. You should have received ",
          "a message describing possible remedies when the reference model ",
          "was built."
        )
      }
    }
    pppd_out <- do.call(rbind, lapply(draw_inds, function(i) {
      proj$refmodel$family$ppd(mu[, i], proj$dis[i], weights)
    }))
  }
  return(structure(pppd_out, cats = cats_aug))
}

#' Plot summary statistics of a variable selection
#'
#' This is the [plot()] method for `vsel` objects (returned by [varsel()] or
#' [cv_varsel()]).
#'
#' @inheritParams summary.vsel
#' @param x An object of class `vsel` (returned by [varsel()] or [cv_varsel()]).
#' @param thres_elpd Only relevant if `any(stats %in% c("elpd", "mlpd"))`. The
#'   threshold for the ELPD difference (taking the submodel's ELPD minus the
#'   baseline model's ELPD) above which the submodel's ELPD is considered to be
#'   close enough to the baseline model's ELPD. An equivalent rule is applied in
#'   case of the MLPD. See [suggest_size()] for a formalization. Supplying `NA`
#'   deactivates this.
#'
#' @inherit summary.vsel details
#'
#' @details
#'
#' # Horizontal lines
#'
#' As long as the reference model's performance is computable, it is always
#' shown in the plot as a dashed red horizontal line. If `baseline = "best"`,
#' the baseline model's performance is shown as a dotted black horizontal line.
#' If `!is.na(thres_elpd)` and `any(stats %in% c("elpd", "mlpd"))`, the value
#' supplied to `thres_elpd` (which is automatically adapted internally in case
#' of the MLPD or `deltas = FALSE`) is shown as a dot-dashed gray horizontal
#' line for the reference model and, if `baseline = "best"`, as a long-dashed
#' green horizontal line for the baseline model.
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
#'   # Variable selection (here without cross-validation and with small values
#'   # for `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the
#'   # sake of speed in this example; this is not recommended in general):
#'   vs <- varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                seed = 5555)
#'   print(plot(vs))
#' }
#'
#' @export
plot.vsel <- function(
    x,
    nterms_max = NULL,
    stats = "elpd",
    deltas = FALSE,
    alpha = 0.32,
    baseline = if (!inherits(x$refmodel, "datafit")) "ref" else "best",
    thres_elpd = NA,
    respOrig = TRUE,
    ...
) {
  object <- x
  .validate_vsel_object_stats(object, stats, respOrig = respOrig)
  baseline <- .validate_baseline(object$refmodel, baseline, deltas)

  ## compute all the statistics and fetch only those that were asked
  nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1],
                                        respOrig = respOrig)
  tab <- rbind(
    .tabulate_stats(object, stats, alpha = alpha,
                    nfeat_baseline = nfeat_baseline, respOrig = respOrig, ...),
    .tabulate_stats(object, stats, alpha = alpha, respOrig = respOrig, ...)
  )
  stats_table <- subset(tab, tab$delta == deltas)
  stats_ref <- subset(stats_table, stats_table$size == Inf)
  stats_sub <- subset(stats_table, stats_table$size != Inf)
  stats_bs <- subset(stats_table, stats_table$size == nfeat_baseline)


  if (NROW(stats_sub) == 0) {
    stop(ifelse(length(stats) > 1, "Statistics ", "Statistic "),
         paste0(unique(stats), collapse = ", "), " not available.")
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
  if (baseline == "ref") {
    baseline_pretty <- "reference model"
  } else {
    baseline_pretty <- "best submodel"
  }
  if (deltas) {
    ylab <- paste0("Difference vs. ", baseline_pretty)
  } else {
    ylab <- "Value"
  }

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

  if (!is.na(thres_elpd)) {
    # Table of thresholds used in extended suggest_size() heuristics (only in
    # case of ELPD and MLPD):
    nobs_test <- nrow(object$d_test$data %||% object$refmodel$fetch_data())
    thres_tab_basic <- data.frame(statistic = c("elpd", "mlpd"),
                                  thres = c(thres_elpd, thres_elpd / nobs_test))
  }

  # plot submodel results
  pp <- ggplot(data = subset(stats_sub, stats_sub$size <= nterms_max),
               mapping = aes_string(x = "size"))
  if (!all(is.na(stats_ref$se))) {
    # add reference model results if they exist

    pp <- pp +
      # The reference model's dashed red horizontal line:
      geom_hline(aes_string(yintercept = "value"),
                 data = stats_ref,
                 color = "darkred", linetype = 2)

    if (!is.na(thres_elpd)) {
      # The thresholds used in extended suggest_size() heuristics:
      thres_tab_ref <- merge(thres_tab_basic,
                             stats_ref[, c("statistic", "value")],
                             by = "statistic")
      thres_tab_ref$thres <- thres_tab_ref$value + thres_tab_ref$thres
      pp <- pp +
        geom_hline(aes_string(yintercept = "thres"),
                   data = thres_tab_ref,
                   color = "gray50", linetype = "dotdash")
    }
  }
  if (baseline != "ref") {
    # add baseline model results (if different from the reference model)

    pp <- pp +
      # The baseline model's dotted black horizontal line:
      geom_hline(aes_string(yintercept = "value"),
                 data = stats_bs,
                 color = "black", linetype = 3)

    if (!is.na(thres_elpd)) {
      # The thresholds used in extended suggest_size() heuristics:
      thres_tab_bs <- merge(thres_tab_basic,
                            stats_bs[, c("statistic", "value")],
                            by = "statistic")
      thres_tab_bs$thres <- thres_tab_bs$value + thres_tab_bs$thres
      pp <- pp +
        geom_hline(aes_string(yintercept = "thres"),
                   data = thres_tab_bs,
                   color = "darkgreen", linetype = "longdash")
    }
  }
  pp <- pp +
    # The submodel-specific graphical elements:
    geom_linerange(aes_string(ymin = "lq", ymax = "uq", alpha = 0.1)) +
    geom_line(aes_string(y = "value")) +
    geom_point(aes_string(y = "value")) +
    # Miscellaneous stuff (axes, theming, faceting, etc.):
    scale_x_continuous(
      breaks = breaks, minor_breaks = minor_breaks,
      limits = c(min(breaks), max(breaks))
    ) +
    labs(x = "Submodel size (number of predictor terms)", y = ylab) +
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
#' @param stats One or more character strings determining which performance
#'   statistics (i.e., utilities or losses) to calculate. Available statistics
#'   are:
#'   * `"elpd"`: (expected) sum of log predictive densities.
#'   * `"mlpd"`: mean log predictive density, that is, `"elpd"` divided by the
#'   number of observations.
#'   * `"mse"`: mean squared error (only available in the situations mentioned
#'   in section "Details" below).
#'   * `"rmse"`: root mean squared error (only available in the situations
#'   mentioned in section "Details" below). For the corresponding standard error
#'   and lower and upper confidence interval bounds, bootstrapping is used.
#'   * `"acc"` (or its alias, `"pctcorr"`): classification accuracy (only
#'   available in the situations mentioned in section "Details" below).
#'   * `"auc"`: area under the ROC curve (only available in the situations
#'   mentioned in section "Details" below). For the corresponding standard error
#'   and lower and upper confidence interval bounds, bootstrapping is used.
#' @param type One or more items from `"mean"`, `"se"`, `"lower"`, `"upper"`,
#'   `"diff"`, and `"diff.se"` indicating which of these to compute for each
#'   item from `stats` (mean, standard error, lower and upper confidence
#'   interval bounds, mean difference to the corresponding statistic of the
#'   reference model, and standard error of this difference, respectively). The
#'   confidence interval bounds belong to normal-approximation (or bootstrap;
#'   see argument `stats`) confidence intervals with (nominal) coverage `1 -
#'   alpha`. Items `"diff"` and `"diff.se"` are only supported if `deltas` is
#'   `FALSE`.
#' @param deltas If `TRUE`, the submodel statistics are estimated as differences
#'   from the baseline model (see argument `baseline`) instead of estimating the
#'   actual values of the statistics.
#' @param alpha A number determining the (nominal) coverage `1 - alpha` of the
#'   normal-approximation confidence intervals. For example, `alpha = 0.32`
#'   corresponds to a coverage of 68%, i.e., one-standard-error intervals
#'   (because of the normal approximation).
#' @param baseline For [summary.vsel()]: Only relevant if `deltas` is `TRUE`.
#'   For [plot.vsel()]: Always relevant. Either `"ref"` or `"best"`, indicating
#'   whether the baseline is the reference model or the best submodel found (in
#'   terms of `stats[1]`), respectively.
#' @param respOrig Only relevant for the latent projection. A single logical
#'   value indicating whether to calculate the performance statistics on
#'   response scale (`TRUE`) or on latent scale (`FALSE`).
#' @param ... Arguments passed to the internal function which is used for
#'   bootstrapping (if applicable; see argument `stats`). Currently, relevant
#'   arguments are `B` (the number of bootstrap samples, defaulting to `2000`)
#'   and `seed` (see [set.seed()], defaulting to
#'   `sample.int(.Machine$integer.max, 1)`, but can also be `NA` to not call
#'   [set.seed()] at all).
#'
#' @details The `stats` options `"mse"` and `"rmse"` are only available for:
#'   * the traditional projection,
#'   * the latent projection with `respOrig = FALSE`,
#'   * the latent projection with `respOrig = TRUE` in combination with
#'   `<refmodel>$family$cats` being `NULL`.
#'
#'   The `stats` option `"acc"` (= `"pctcorr"`) is only available for:
#'   * the [binomial()] family in case of the traditional projection,
#'   * all families in case of the augmented-data projection,
#'   * the [binomial()] family (on the original response scale) in case of the
#'   latent projection with `respOrig = TRUE` in combination with
#'   `<refmodel>$family$cats` being `NULL`,
#'   * all families (on the original response scale) in case of the latent
#'   projection with `respOrig = TRUE` in combination with
#'   `<refmodel>$family$cats` being not `NULL`.
#'
#'   The `stats` option `"auc"` is only available for:
#'   * the [binomial()] family in case of the traditional projection,
#'   * the [binomial()] family (on the original response scale) in case of the
#'   latent projection with `respOrig = TRUE` in combination with
#'   `<refmodel>$family$cats` being `NULL`.
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
#'   # Variable selection (here without cross-validation and with small values
#'   # for `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the
#'   # sake of speed in this example; this is not recommended in general):
#'   vs <- varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                seed = 5555)
#'   print(summary(vs))
#' }
#'
#' @export
summary.vsel <- function(
    object,
    nterms_max = NULL,
    stats = "elpd",
    type = c("mean", "se", "diff", "diff.se"),
    deltas = FALSE,
    alpha = 0.32,
    baseline = if (!inherits(object$refmodel, "datafit")) "ref" else "best",
    respOrig = TRUE,
    ...
) {
  .validate_vsel_object_stats(object, stats, respOrig = respOrig)
  baseline <- .validate_baseline(object$refmodel, baseline, deltas)

  # Initialize output:
  out <- list(
    formula = object$refmodel$formula,
    family = object$refmodel$family,
    nobs_train = nrow(object$refmodel$fetch_data()),
    nobs_test = nrow(object$d_test$data),
    method = object$method,
    cv_method = object$cv_method,
    validate_search = object$validate_search,
    clust_used_search = object$clust_used_search,
    clust_used_eval = object$clust_used_eval,
    nprjdraws_search = object$nprjdraws_search,
    nprjdraws_eval = object$nprjdraws_eval
  )
  if (isTRUE(out$validate_search)) {
    out$search_included <- "search included"
  } else {
    out$search_included <- "search not included"
  }
  class(out) <- "vselsummary"

  # The full table of the performance statistics from `stats`:
  if (deltas) {
    nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1],
                                          respOrig = respOrig)
    tab <- .tabulate_stats(object, stats, alpha = alpha,
                           nfeat_baseline = nfeat_baseline, respOrig = respOrig,
                           ...)
  } else {
    tab <- .tabulate_stats(object, stats, alpha = alpha, respOrig = respOrig,
                           ...)
  }
  stats_table <- subset(tab, tab$size != Inf) %>%
    dplyr::group_by(.data$statistic) %>%
    dplyr::slice_head(n = length(object$solution_terms) + 1)

  # Get the names of `stats_table` corresponding to all items from `type`, and
  # set up their suffices in the table to be returned:
  if (deltas) {
    type <- setdiff(type, c("diff", "diff.se"))
  }
  qty <- unname(sapply(type, function(t) {
    switch(t, mean = "value", upper = "uq", lower = "lq", se = "se",
           diff = "diff", diff.se = "diff.se")
  }))
  if (!is.null(object$cv_method)) {
    cv_suffix <- unname(switch(object$cv_method,
                               LOO = ".loo", kfold = ".kfold"))
  } else {
    cv_suffix <- NULL
  }
  if (length(stats) > 1) {
    suffix <- lapply(stats, function(s) {
      unname(sapply(type, function(t) {
        paste0(s,
               switch(t, mean = cv_suffix, upper = ".upper", lower = ".lower",
                      se = ".se", diff = ".diff", diff.se = ".diff.se"))
      }))
    })
  } else {
    suffix <- list(unname(sapply(type, function(t) {
      switch(t, mean = paste0(stats, cv_suffix), upper = "upper",
             lower = "lower", se = "se",
             diff = "diff", diff.se = "diff.se")
    })))
  }

  # Construct the (almost) final output table by looping over all requested
  # statistics, reshaping the corresponding data in `stats_table`, and selecting
  # only the requested `type`s:
  arr <- data.frame(size = unique(stats_table$size),
                    solution_terms = c(NA, object$solution_terms))
  for (i in seq_along(stats)) {
    temp <- subset(stats_table, stats_table$statistic == stats[i], qty)
    newnames <- suffix[[i]]
    colnames(temp) <- newnames
    arr <- cbind(arr, temp)
  }

  # Output (and also cut `arr` at `nterms_max` (if provided)):
  if (is.null(nterms_max)) {
    nterms_max <- max(stats_table$size)
  }
  out$nterms <- nterms_max
  if ("pct_solution_terms_cv" %in% names(object)) {
    out$pct_solution_terms_cv <- object$pct_solution_terms_cv
  }
  out$selection <- subset(arr, arr$size <= nterms_max)
  out$respOrig <- respOrig
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
  if (is.null(x$nobs_test)) {
    cat(paste0("Observations: ", x$nobs_train, "\n"))
  } else {
    cat(paste0("Observations (training set): ", x$nobs_train, "\n"))
    cat(paste0("Observations (test set): ", x$nobs_test, "\n"))
  }
  if (x$family$for_augdat) {
    prj_meth <- "augmented-data"
  } else if (x$family$for_latent) {
    prj_meth <- "latent"
  } else {
    prj_meth <- "traditional"
  }
  cat("Projection method: ", prj_meth, "\n", sep = "")
  if (!is.null(x$cv_method)) {
    cat(paste("CV method:", x$cv_method, x$search_included, "\n"))
  }
  cat(paste0("Search method: ", x$method, ", maximum number of terms ",
             max(x$selection$size), "\n"))
  cat("Number of ", ifelse(x$clust_used_search, "clusters", "draws"),
      " used for selection: ", x$nprjdraws_search, "\n", sep = "")
  cat("Number of ", ifelse(x$clust_used_eval, "clusters", "draws"),
      " used for prediction: ", x$nprjdraws_eval, "\n", sep = "")
  cat("\n")
  if (x$family$for_latent) {
    if (x$respOrig) {
      scale_string <- " (response scale)"
    } else {
      scale_string <- " (latent scale)"
    }
  } else {
    scale_string <- ""
  }
  cat("Selection Summary", scale_string, ":\n", sep = "")
  where <- "tidyselect" %:::% "where"
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
#' @return The output of [summary.vsel()] (invisible).
#'
#' @export
print.vsel <- function(x, ...) {
  dot_args <- list(...)
  stats <- do.call(summary.vsel, c(list(object = x),
                                   dot_args[names(dot_args) != "digits"]))
  do.call(print, c(list(x = stats),
                   dot_args[names(dot_args) == "digits"]))
  return(invisible(stats))
}

#' Suggest submodel size
#'
#' This function can suggest an appropriate submodel size based on a decision
#' rule described in section "Details" below. Note that this decision is quite
#' heuristic and should be interpreted with caution. It is recommended to
#' examine the results via [plot.vsel()] and/or [summary.vsel()] and to make the
#' final decision based on what is most appropriate for the problem at hand.
#'
#' @param object An object of class `vsel` (returned by [varsel()] or
#'   [cv_varsel()]).
#' @param stat Performance statistic (i.e., utility or loss) used for the
#'   decision. See argument `stats` of [summary.vsel()] for possible choices.
#' @param pct A number giving the relative proportion (*not* percents) between
#'   baseline model and null model utilities one is willing to sacrifice. See
#'   section "Details" below for more information.
#' @param type Either `"upper"` or `"lower"` determining whether the decision is
#'   based on the upper or lower confidence interval bound, respectively. See
#'   section "Details" below for more information.
#' @param thres_elpd Only relevant if `stat %in% c("elpd", "mlpd")`. The
#'   threshold for the ELPD difference (taking the submodel's ELPD minus the
#'   baseline model's ELPD) above which the submodel's ELPD is considered to be
#'   close enough to the baseline model's ELPD. An equivalent rule is applied in
#'   case of the MLPD. See section "Details" for a formalization. Supplying `NA`
#'   deactivates this.
#' @param warnings Mainly for internal use. A single logical value indicating
#'   whether to throw warnings if automatic suggestion fails. Usually there is
#'   no reason to set this to `FALSE`.
#' @param ... Arguments passed to [summary.vsel()], except for `object`, `stats`
#'   (which is set to `stat`), `type`, and `deltas` (which is set to `TRUE`).
#'   See section "Details" below for some important arguments which may be
#'   passed here.
#'
#' @details In general (beware of special extensions below), the suggested model
#'   size is the smallest model size \eqn{k \in \{0, 1, ...,
#'   \texttt{nterms\_max}\}}{{k = 0, 1, ..., nterms_max}} for which either the
#'   lower or upper bound (depending on argument `type`) of the
#'   normal-approximation confidence interval (with nominal coverage `1 -
#'   alpha`; see argument `alpha` of [summary.vsel()]) for \eqn{U_k -
#'   U_{\mathrm{base}}}{U_k - U_base} (with \eqn{U_k} denoting the \eqn{k}-th
#'   submodel's true utility and \eqn{U_{\mathrm{base}}}{U_base} denoting the
#'   baseline model's true utility) falls above (or is equal to)
#'   \deqn{\texttt{pct} \cdot (u_0 - u_{\mathrm{base}})}{pct * (u_0 - u_base)}
#'   where \eqn{u_0} denotes the null model's estimated utility and
#'   \eqn{u_{\mathrm{base}}}{u_base} the baseline model's estimated utility. The
#'   baseline model is either the reference model or the best submodel found
#'   (see argument `baseline` of [summary.vsel()]).
#'
#'   If `!is.na(thres_elpd)` and `stat = "elpd"`, the decision rule above is
#'   extended: The suggested model size is then the smallest model size \eqn{k}
#'   fulfilling the rule above *or* \eqn{u_k - u_{\mathrm{base}} >
#'   \texttt{thres\_elpd}}{u_k - u_base > thres_elpd}. Correspondingly, in case
#'   of `stat = "mlpd"` (and `!is.na(thres_elpd)`), the suggested model size is
#'   the smallest model size \eqn{k} fulfilling the rule above *or* \eqn{u_k -
#'   u_{\mathrm{base}} > \frac{\texttt{thres\_elpd}}{N}}{u_k - u_base >
#'   thres_elpd / N} with \eqn{N} denoting the number of observations.
#'
#'   For example (disregarding the special extensions in case of `stat = "elpd"`
#'   or `stat = "mlpd"`), `alpha = 0.32`, `pct = 0`, and `type = "upper"` means
#'   that we select the smallest model size for which the upper bound of the 68%
#'   confidence interval for \eqn{U_k - U_{\mathrm{base}}}{U_k - U_base} exceeds
#'   (or is equal to) zero, that is, for which the submodel's utility estimate
#'   is at most one standard error smaller than the baseline model's utility
#'   estimate (with that standard error referring to the utility *difference*).
#'
#'   Apart from the two [summary.vsel()] arguments mentioned above (`alpha` and
#'   `baseline`), `respOrig` is another important [summary.vsel()] argument that
#'   may be passed via `...`.
#'
#' @note Loss statistics like the root mean-squared error (RMSE) and the
#'   mean-squared error (MSE) are converted to utilities by multiplying them by
#'   `-1`, so a call such as `suggest_size(object, stat = "rmse", type =
#'   "upper")` finds the smallest model size whose upper confidence interval
#'   bound for the *negative* RMSE or MSE exceeds the cutoff (or, equivalently,
#'   has the lower confidence interval bound for the RMSE or MSE below the
#'   cutoff). This is done to make the interpretation of argument `type` the
#'   same regardless of argument `stat`.
#'
#'   The intercept is not counted by [suggest_size()], so a suggested size of
#'   zero stands for the intercept-only model.
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
suggest_size.vsel <- function(
    object,
    stat = "elpd",
    pct = 0,
    type = "upper",
    thres_elpd = NA,
    warnings = TRUE,
    ...
) {
  if (length(stat) > 1) {
    stop("Only one statistic can be specified to suggest_size")
  }
  stats <- summary.vsel(object,
                        stats = stat,
                        type = c("mean", "upper", "lower"),
                        deltas = TRUE,
                        ...)
  nobs_test <- stats$nobs_test %||% stats$nobs_train
  stats <- stats$selection

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

  util_null <- sgn * unlist(unname(subset(
    stats, stats$size == 0,
    paste0(stat, suffix)
  )))
  util_cutoff <- pct * util_null
  if (is.na(thres_elpd)) {
    thres_elpd <- Inf
  }
  res <- stats[
    (sgn * stats[, bound] >= util_cutoff) |
      (stat == "elpd" & stats[, paste0(stat, suffix)] > thres_elpd) |
      (stat == "mlpd" & stats[, paste0(stat, suffix)] > thres_elpd / nobs_test),
    "size", drop = FALSE
  ]

  if (nrow(res) == 0) {
    ## no submodel satisfying the criterion found
    if (object$nterms_max == object$nterms_all) {
      suggested_size <- object$nterms_max
    } else {
      suggested_size <- NA
      if (warnings) {
        warning("Could not suggest model size. Investigate plot.vsel() to ",
                "identify if the search was terminated too early. If this is ",
                "the case, run variable selection with larger value for ",
                "`nterms_max`.")
      }
    }
  } else {
    # Above, `object$nterms_max` includes the intercept (if present), so we need
    # to include it here, too:
    suggested_size <- min(res) + object$refmodel$intercept
    # We don't use `na.rm = TRUE` in min() to be as cautious as possible. In
    # fact, we could refine this to remove `NA`s after the first non-`NA` value
    # (meaning that if there is no non-`NA` value, no `NA`s will be removed),
    # but this gets overly complicated and it's better to be as cautious as
    # possible (because `NA`s after the first non-`NA` value are also strange).
  }

  return(suggested_size - object$refmodel$intercept)
}

# Make the parameter name(s) for the intercept(s) adhere to the naming scheme
# `nm_scheme`:
mknms_icpt <- function(nms, nm_scheme) {
  if (nm_scheme == "brms") {
    nms <- gsub("\\(Intercept\\)", "Intercept", nms)
  }
  return(nms)
}

# Replace the names of an object containing population-level effects so that
# these names adhere to the naming scheme `nm_scheme`:
replace_population_names <- function(population_effects, nm_scheme) {
  if (nm_scheme == "brms") {
    # Use brms's naming convention:
    names(population_effects) <- mknms_icpt(
      names(population_effects),
      nm_scheme = nm_scheme
    )
    if (length(population_effects) > 0) {
      # We could also use `recycle0 = TRUE` here, but that would
      # require R >= 4.0.1.
      names(population_effects) <- paste0("b_", names(population_effects))
    }
  }
  return(population_effects)
}

# Make the parameter names for variance components adhere to the naming scheme
# `nm_scheme`:
mknms_VarCorr <- function(nms, nms_lats = NULL, nm_scheme, coef_nms) {
  if (!is.null(nms_lats)) {
    stopifnot(nm_scheme == "brms")
  }
  grp_nms <- names(coef_nms)
  # We will have to search for the substrings "\\sd\\." and "\\cor\\.", so make
  # sure that they don't occur in the coefficient or group names:
  stopifnot(!any(grepl("\\.sd\\.|\\.cor\\.", grp_nms)))
  stopifnot(!any(unlist(lapply(
    coef_nms, grepl, pattern = "\\.sd\\.|\\.cor\\."
  ))))
  if (nm_scheme == "brms") {
    nms <- mknms_icpt(nms, nm_scheme = nm_scheme)
    # Escape special characters in the group names and collapse them with "|":
    grp_nms_esc <- paste(gsub("\\)", "\\\\)",
                              gsub("\\(", "\\\\(",
                                   gsub("\\.", "\\\\.", grp_nms))),
                         collapse = "|")
    # Move the substrings "\\.sd\\." and "\\.cor\\." up front (i.e. in front of
    # the group name), replace their dots, and replace the dot following the
    # group name by double underscores:
    nms <- sub(paste0("(", grp_nms_esc, ")\\.(sd|cor)\\."),
               "\\2_\\1__",
               nms)
  }
  for (coef_nms_i in coef_nms) {
    if (nm_scheme == "brms") {
      coef_nms_i <- mknms_icpt(coef_nms_i, nm_scheme = nm_scheme)
    }
    # Escape special characters in the coefficient names and collapse them
    # with "|":
    coef_nms_i_esc <- paste(gsub("\\)", "\\\\)",
                                 gsub("\\(", "\\\\(",
                                      gsub("\\.", "\\\\.", coef_nms_i))),
                            collapse = "|")
    if (nm_scheme == "brms") {
      # Replace dots between coefficient names by double underscores:
      nms <- gsub(paste0("(", coef_nms_i_esc, ")\\."),
                  "\\1__",
                  nms)
    } else if (nm_scheme == "rstanarm") {
      # For the substring "\\.sd\\.":
      nms <- sub(paste0("\\.sd\\.(", coef_nms_i_esc, ")$"),
                 ":\\1,\\1",
                 nms)
      # For the substring "\\.cor\\.":
      nms <- sub(
        paste0("\\.cor\\.(", coef_nms_i_esc, ")\\.(", coef_nms_i_esc, ")$"),
        ":\\2,\\1",
        nms
      )
    }
  }
  if (nm_scheme == "rstanarm") {
    nms <- paste0("Sigma[", nms, "]")
  }
  if (!is.null(nms_lats)) {
    # Escape special characters in the latent category names and collapse them
    # with "|":
    nms_lats_esc <- paste(gsub("\\)", "\\\\)",
                               gsub("\\(", "\\\\(",
                                    gsub("\\.", "\\\\.", nms_lats))),
                          collapse = "|")
    # Put the string `mu` in front of the latent category names and replace the
    # following tilde by an underscore:
    nms <- gsub(paste0("(", nms_lats_esc, ")~"),
                "mu\\1_",
                nms)
  }
  return(nms)
}

# Make the parameter names for group-level effects adhere to the naming scheme
# `nm_scheme`:
mknms_ranef <- function(nms, nms_lats = NULL, nm_scheme, coef_nms) {
  if (!is.null(nms_lats)) {
    stopifnot(nm_scheme == "brms")
  }
  if (nm_scheme == "brms") {
    nms <- mknms_icpt(nms, nm_scheme = nm_scheme)
  }
  for (coef_nms_idx in seq_along(coef_nms)) {
    coef_nms_i <- coef_nms[[coef_nms_idx]]
    if (nm_scheme == "brms") {
      coef_nms_i <- mknms_icpt(coef_nms_i, nm_scheme = nm_scheme)
    }
    # Escape special characters in the coefficient names and collapse them with
    # "|":
    coef_nms_i_esc <- paste(gsub("\\)", "\\\\)",
                                 gsub("\\(", "\\\\(",
                                      gsub("\\.", "\\\\.", coef_nms_i))),
                            collapse = "|")
    if (nm_scheme == "brms") {
      # Put the part following the group name in square brackets, reorder its
      # two subparts (coefficient name and group level), and separate them by
      # comma:
      nms <- sub(paste0("\\.(", coef_nms_i_esc, ")\\.(.*)$"),
                 "[\\2,\\1]",
                 nms)
    } else if (nm_scheme == "rstanarm") {
      grp_nm_i <- names(coef_nms)[coef_nms_idx]
      # Escape special characters in the group name:
      grp_nm_i_esc <- gsub("\\)", "\\\\)",
                           gsub("\\(", "\\\\(",
                                gsub("\\.", "\\\\.", grp_nm_i)))
      # Re-arrange as required:
      nms <- sub(paste0("^(", grp_nm_i_esc, ")\\.(", coef_nms_i_esc, ")\\."),
                 "\\2 \\1:",
                 nms)
    }
  }
  if (nm_scheme == "brms") {
    nms <- paste0("r_", nms)
  } else if (nm_scheme == "rstanarm") {
    nms <- paste0("b[", nms, "]")
  }
  if (!is.null(nms_lats)) {
    # Escape special characters in the latent category names and collapse them
    # with "|":
    nms_lats_esc <- paste(gsub("\\)", "\\\\)",
                               gsub("\\(", "\\\\(",
                                    gsub("\\.", "\\\\.", nms_lats))),
                          collapse = "|")
    # Put the string `mu` in front of the latent category names, remove the
    # following tilde, and place all this in front of the first square bracket:
    nms <- gsub(paste0("\\[(.*),(", nms_lats_esc, ")~"),
                "__mu\\2[\\1,",
                nms)
  }
  return(nms)
}

# Make the parameter names for the thresholds of an ordinal model adhere to the
# naming scheme `nm_scheme`:
mknms_thres <- function(nms, nm_scheme) {
  if (nm_scheme == "brms") {
    nms <- paste0("b_Intercept[", seq_along(nms), "]")
  }
  return(nms)
}

# Make the non-multilevel parameter names of a categorical model adhere to the
# naming scheme `nm_scheme`:
mknms_categ <- function(dimnms, nm_scheme) {
  # rstanarm currently doesn't support categorical models:
  stopifnot(nm_scheme == "brms")
  nmsdf <- expand.grid(dimnms, stringsAsFactors = FALSE)
  nmsdf[, 1] <- paste0("mu", nmsdf[, 1])
  nmsdf[, 2] <- mknms_icpt(nmsdf[, 2], nm_scheme = nm_scheme)
  nmsdf <- cbind("b", nmsdf)
  return(apply(nmsdf, 1, paste, collapse = "_"))
}

#' @noRd
#' @export
coef.subfit <- function(object, ...) {
  return(with(object, c(
    "(Intercept)" = alpha,
    setNames(beta, colnames(x))
  )))
}

# To process the multilevel variance components (from a submodel fit):
proc_VarCorr <- function(group_vc_raw, nms_lats = NULL, ...) {
  group_vc <- unlist(lapply(group_vc_raw, function(vc_obj) {
    # The vector of standard deviations:
    if (is.null(nms_lats)) {
      vc_out <- c("sd" = attr(vc_obj, "stddev"))
    } else {
      vc_out <- c("sd" = sqrt(diag(vc_obj)))
    }
    # The correlation matrix:
    if (is.null(nms_lats)) {
      cor_mat <- attr(vc_obj, "correlation")
      has_cor <- !is.null(cor_mat)
    } else {
      cor_mat <- cov2cor(vc_obj)
      has_cor <- ncol(cor_mat) > 1
    }
    if (has_cor) {
      # Auxiliary object: A matrix of the same dimension as cor_mat, but
      # containing the paste()-d dimnames:
      cor_mat_nms <- matrix(
        apply(expand.grid(rownames(cor_mat),
                          colnames(cor_mat)),
              1,
              paste,
              collapse = "."),
        nrow = nrow(cor_mat),
        ncol = ncol(cor_mat)
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
  names(group_vc) <- mknms_VarCorr(names(group_vc), nms_lats = nms_lats, ...)
  return(group_vc)
}

# To process the raw group-level effects themselves (from a submodel fit):
proc_ranef <- function(group_ef_raw, nms_lats = NULL, ncoefs, grps_lvls, VarCov,
                       ...) {
  if (!is.null(nms_lats)) {
    coef_nms <- list(...)$coef_nms
    stopifnot(!is.null(coef_nms))
    nlats <- length(nms_lats)
    group_ef_raw <- lapply(setNames(nm = names(group_ef_raw)), function(vnm) {
      ranef_tmp <- group_ef_raw[[vnm]]
      if (utils::packageVersion("mclogit") < "0.9") {
        ncoefs_vnm <- ncoefs
      } else {
        ncoefs_vnm <- ncoefs[vnm]
      }
      # Coerce the random effects into the same format as the output of ranef()
      # from packages 'lme4' and 'ordinal':
      ranef_tmp <- matrix(ranef_tmp,
                          nrow = nlats * ncoefs_vnm,
                          ncol = length(grps_lvls[[vnm]]),
                          dimnames = list(coef_nms[[vnm]],
                                          grps_lvls[[vnm]]))
      return(as.data.frame(t(ranef_tmp)))
    })
  }
  group_ef <- unlist(lapply(group_ef_raw, function(ranef_df) {
    ranef_mat <- as.matrix(ranef_df)
    setNames(
      as.vector(ranef_mat),
      apply(expand.grid(rownames(ranef_mat),
                        colnames(ranef_mat)),
            1,
            function(row_col_nm) {
              paste(rev(row_col_nm), collapse = ".")
            })
    )
  }))
  names(group_ef) <- mknms_ranef(names(group_ef), nms_lats = nms_lats, ...)
  return(group_ef)
}

# An (internal) generic for extracting the coefficients and any other parameter
# estimates from a submodel fit.
get_subparams <- function(x, ...) {
  UseMethod("get_subparams")
}

#' @noRd
#' @export
get_subparams.lm <- function(x, ...) {
  return(coef(x) %>%
           replace_population_names(...))
}

#' @noRd
#' @export
get_subparams.subfit <- function(x, ...) {
  return(get_subparams.lm(x, ...))
}

#' @noRd
#' @export
get_subparams.glm <- function(x, ...) {
  return(get_subparams.lm(x, ...))
}

#' @noRd
#' @export
get_subparams.lmerMod <- function(x, ...) {
  population_effects <- lme4::fixef(x) %>%
    replace_population_names(...)

  group_vc_raw <- lme4::VarCorr(x)
  group_vc <- proc_VarCorr(group_vc_raw,
                           coef_nms = lapply(group_vc_raw, rownames), ...)

  group_ef <- proc_ranef(lme4::ranef(x, condVar = FALSE),
                         coef_nms = lapply(group_vc_raw, rownames), ...)

  return(c(population_effects, group_vc, group_ef))
}

#' @noRd
#' @export
get_subparams.glmerMod <- function(x, ...) {
  return(get_subparams.lmerMod(x, ...))
}

#' @noRd
#' @export
get_subparams.gamm4 <- function(x, ...) {
  return(get_subparams.lm(x, ...))
}

#' @noRd
#' @export
get_subparams.polr <- function(x, ...) {
  thres <- x$zeta
  names(thres) <- mknms_thres(names(thres), ...)
  return(c(thres, get_subparams.lm(x, ...)))
}

#' @noRd
#' @export
get_subparams.clmm <- function(x, ...) {
  thres <- x$alpha
  names(thres) <- mknms_thres(names(thres), ...)

  group_vc_raw <- ordinal::VarCorr(x)
  group_vc <- proc_VarCorr(group_vc_raw,
                           coef_nms = lapply(group_vc_raw, rownames), ...)

  group_ef <- proc_ranef(ordinal::ranef(x),
                         coef_nms = lapply(group_vc_raw, rownames), ...)

  return(c(thres, replace_population_names(x$beta, ...), group_vc, group_ef))
}

#' @noRd
#' @export
get_subparams.multinom <- function(x, ...) {
  coefs <- coef(x)
  nms <- mknms_categ(dimnames(coefs), ...)
  return(setNames(as.vector(coefs), nms))
}

#' @noRd
#' @export
get_subparams.mmblogit <- function(x, ...) {
  coefs <- x$coefmat
  group_vc_raw <- x$VarCov
  if (utils::packageVersion("mclogit") < "0.9") {
    group_vc_raw <- setNames(group_vc_raw, names(x$groups))
  }
  group_vc_raw <- lapply(group_vc_raw, function(vc_obj) {
    nms_lats_coefs <- rownames(vc_obj)
    stopifnot(identical(nms_lats_coefs, colnames(vc_obj)))
    tilde_check <- gregexpr("~", nms_lats_coefs)
    stopifnot(all(lengths(tilde_check) == 1))
    stopifnot(all(unlist(tilde_check) != -1))
    rownames(vc_obj) <- colnames(vc_obj) <-
      sub("~1$", "~(Intercept)", nms_lats_coefs)
    return(vc_obj)
  })
  group_vc <- proc_VarCorr(group_vc_raw, nms_lats = colnames(x$D),
                           coef_nms = lapply(group_vc_raw, rownames), ...)

  if (utils::packageVersion("mclogit") < "0.9") {
    ncoefs_all <- length(all.vars(x$random$formula)) + 1L
  } else {
    ncoefs_all <- sapply(
      setNames(x$random, names(x$groups)),
      function(re_info_i) {
        length(all.vars(re_info_i$formula)) + 1L
      }
    )
  }
  group_ef <- proc_ranef(setNames(x$random.effects, names(x$groups)),
                         nms_lats = colnames(x$D),
                         ncoefs = ncoefs_all,
                         grps_lvls = lapply(x$groups, levels),
                         coef_nms = lapply(group_vc_raw, rownames), ...)

  nms <- mknms_categ(dimnames(coefs), ...)
  return(c(setNames(as.vector(coefs), nms), group_vc, group_ef))
}

#' Extract projected parameter draws
#'
#' This is the [as.matrix()] method for `projection` objects (returned by
#' [project()], possibly as elements of a `list`). It extracts the projected
#' parameter draws and returns them as a matrix.
#'
#' @param x An object of class `projection` (returned by [project()], possibly
#'   as elements of a `list`).
#' @param nm_scheme The naming scheme for the columns of the output matrix.
#'   Either `"auto"`, `"rstanarm"`, or `"brms"`, where `"auto"` chooses
#'   `"rstanarm"` or `"brms"` based on the class of the reference model fit (and
#'   uses `"rstanarm"` if the reference model fit is of an unknown class).
#' @param ... Currently ignored.
#'
#' @details In case of the augmented-data projection for a multilevel submodel
#'   of a [brms::categorical()] reference model, the multilevel parameters (and
#'   therefore also their names) slightly differ from those in the \pkg{brms}
#'   reference model fit (see section "Augmented-data projection" in
#'   [extend_family()]'s documentation).
#'
#' @return An \eqn{S_{\mathrm{prj}} \times Q}{S_prj x Q} matrix of projected
#'   draws, with \eqn{S_{\mathrm{prj}}}{S_prj} denoting the number of projected
#'   draws and \eqn{Q} the number of parameters.
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
#'   # Projection onto an arbitrary combination of predictor terms (with a small
#'   # value for `nclusters`, but only for the sake of speed in this example;
#'   # this is not recommended in general):
#'   prj <- project(fit, solution_terms = c("X1", "X3", "X5"), nclusters = 10,
#'                  seed = 9182)
#'   prjmat <- as.matrix(prj)
#'   ### For further post-processing (e.g., via packages `bayesplot` and
#'   ### `posterior`), we will here ignore the fact that clustering was used
#'   ### (due to argument `nclusters` above). CAUTION: Ignoring the clustering
#'   ### is not recommended and only shown here for demonstrative purposes. A
#'   ### better solution for the clustering case is explained below.
#'   # If the `bayesplot` package is installed, the output from
#'   # as.matrix.projection() can be used there. For example:
#'   if (requireNamespace("bayesplot", quietly = TRUE)) {
#'     print(bayesplot::mcmc_intervals(prjmat))
#'   }
#'   # If the `posterior` package is installed, the output from
#'   # as.matrix.projection() can be used there. For example:
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
as.matrix.projection <- function(x, nm_scheme = "auto", ...) {
  if (inherits(x$refmodel, "datafit")) {
    stop("as.matrix.projection() does not work for objects based on ",
         "\"datafit\"s.")
  }
  if (x$p_type) {
    warning("Note that projection was performed using clustering and the ",
            "clusters might have different weights.")
  }
  if (identical(nm_scheme, "auto")) {
    if (inherits(x$refmodel$fit, "brmsfit")) {
      nm_scheme <- "brms"
    } else {
      nm_scheme <- "rstanarm"
    }
  }
  stopifnot(nm_scheme %in% c("rstanarm", "brms"))
  res <- do.call(rbind, lapply(x$submodl, get_subparams, nm_scheme = nm_scheme))
  if (x$refmodel$family$family == "gaussian") res <- cbind(res, sigma = x$dis)
  return(res)
}

#' Create cross-validation folds
#'
#' These are helper functions to create cross-validation (CV) folds, i.e., to
#' split up the indices from 1 to `n` into `K` subsets ("folds") for
#' \eqn{K}-fold CV. These functions are potentially useful when creating the
#' `cvfits` and `cvfun` arguments for [init_refmodel()]. The return value is
#' different for these two methods, see below for details.
#'
#' @name cv-indices
#'
#' @param n Number of observations.
#' @param K Number of folds. Must be at least 2 and not exceed `n`.
#' @param out Format of the output, either `"foldwise"` or `"indices"`. See
#'   below for details.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all.
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
cvfolds <- function(n, K, seed = sample.int(.Machine$integer.max, 1)) {
  .validate_num_folds(K, n)

  # Set seed, but ensure the old RNG state is restored on exit:
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  if (!is.na(seed)) set.seed(seed)

  ## create and shuffle the indices
  folds <- rep_len(seq_len(K), length.out = n)
  folds <- sample(folds, n, replace = FALSE)

  return(folds)
}

#' @rdname cv-indices
#' @export
cv_ids <- function(n, K, out = c("foldwise", "indices"),
                   seed = sample.int(.Machine$integer.max, 1)) {
  .validate_num_folds(K, n)
  out <- match.arg(out)

  # Set seed, but ensure the old RNG state is restored on exit:
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  if (!is.na(seed)) set.seed(seed)

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
#' combination onto which the projection was performed.
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
#'   # The "stanreg" fit which will be used as the reference model (with small
#'   # values for `chains` and `iter`, but only for technical reasons in this
#'   # example; this is not recommended in general):
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
