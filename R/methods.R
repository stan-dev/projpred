# Common documentation ----------------------------------------------------

#' Predictions from a submodel (after projection)
#'
#' After the projection of the reference model onto a submodel, the linear
#' predictors (for the original or a new dataset) based on that submodel can be
#' calculated by [proj_linpred()]. These linear predictors can also be
#' transformed to response scale and averaged across the projected parameter
#' draws. Furthermore, [proj_linpred()] returns the corresponding log predictive
#' density values if the (original or new) dataset contains response values. The
#' [proj_predict()] function draws from the predictive distributions (there is
#' one such distribution for each observation from the original or new dataset)
#' of the submodel that the reference model has been projected onto. If the
#' projection has not been performed yet, both functions call [project()]
#' internally to perform the projection. Both functions can also handle multiple
#' submodels at once (for `object`s of class `vsel` or `object`s returned by a
#' [project()] call to an object of class `vsel`; see [project()]).
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
#'   `resp_oscale` from other functions and affects the scale of both output
#'   elements `pred` and `lpd` (see sections "Details" and "Value" below).
#' @param integrated For [proj_linpred()] only. A single logical value
#'   indicating whether the output should be averaged across the projected
#'   posterior draws (`TRUE`) or not (`FALSE`).
#' @param nresample_clusters For [proj_predict()] with clustered projection (and
#'   nonconstant weights for the projected draws) only. Number of draws to
#'   return from the predictive distributions of the submodel(s). Not to be
#'   confused with argument `nclusters` of [project()]: `nresample_clusters`
#'   gives the number of draws (*with* replacement) from the set of clustered
#'   posterior draws after projection (with this set being determined by
#'   argument `nclusters` of [project()]).
#' @param .seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. If not
#'   `NA`, then the PRNG state is reset (to the state before calling
#'   [proj_linpred()] or [proj_predict()]) upon exiting [proj_linpred()] or
#'   [proj_predict()]. Here, `.seed` is used for drawing new group-level effects
#'   in case of a multilevel submodel (however, not yet in case of a GAMM) and
#'   for drawing from the predictive distributions of the submodel(s) in case of
#'   [proj_predict()]. If a clustered projection was performed, then in
#'   [proj_predict()], `.seed` is also used for drawing from the set of
#'   projected clusters of posterior draws (see argument `nresample_clusters`).
#'   If [project()] is called internally with `seed = NA` (or with `seed` being
#'   a lazily evaluated expression that uses the PRNG), then `.seed` also
#'   affects the PRNG usage there.
#' @param resp_oscale Only relevant for the latent projection. A single logical
#'   value indicating whether to draw from the posterior-projection predictive
#'   distributions on the original response scale (`TRUE`) or on latent scale
#'   (`FALSE`).
#' @param ... Arguments passed to [project()] if `object` is not already an
#'   object returned by [project()].
#'
#' @details Currently, [proj_predict()] ignores observation weights that are not
#'   equal to `1`. A corresponding warning is thrown if this is the case.
#'
#'   In case of the latent projection and `transform = FALSE`:
#'   * Output element `pred` contains the linear predictors without any
#'   modifications that may be due to the original response distribution (e.g.,
#'   for a [brms::cumulative()] model, the ordered thresholds are not taken into
#'   account).
#'   * Output element `lpd` contains the *latent* log predictive density values,
#'   i.e., those corresponding to the latent Gaussian distribution. If `newdata`
#'   is not `NULL`, this requires the latent response values to be supplied in a
#'   column called `.<response_name>` of `newdata` where `<response_name>` needs
#'   to be replaced by the name of the original response variable (if
#'   `<response_name>` contained parentheses, these have been stripped off by
#'   [init_refmodel()]; see the left-hand side of `formula(<refmodel>)`). For
#'   technical reasons, the existence of column `<response_name>` in `newdata`
#'   is another requirement (even though `.<response_name>` is actually used).
#'
#' @return In the following, \eqn{S_{\mathrm{prj}}}{S_prj}, \eqn{N},
#'   \eqn{C_{\mathrm{cat}}}{C_cat}, and \eqn{C_{\mathrm{lat}}}{C_lat} from help
#'   topic [refmodel-init-get] are used. (For [proj_linpred()] with `integrated
#'   = TRUE`, we have \eqn{S_{\mathrm{prj}} = 1}{S_prj = 1}.) Furthermore, let
#'   \eqn{C} denote either \eqn{C_{\mathrm{cat}}}{C_cat} (if `transform = TRUE`)
#'   or \eqn{C_{\mathrm{lat}}}{C_lat} (if `transform = FALSE`). Then, if the
#'   prediction is done for one submodel only (i.e., `length(nterms) == 1 ||
#'   !is.null(solution_terms)` in the call to [project()]):
#'   * [proj_linpred()] returns a `list` with the following elements:
#'       + Element `pred` contains the actual predictions, i.e., the linear
#'       predictors, possibly transformed to response scale (depending on
#'       argument `transform`).
#'       + Element `lpd` is non-`NULL` only if `newdata` is `NULL` or if
#'       `newdata` contains response values in the corresponding column. In that
#'       case, it contains the log predictive density values (conditional on
#'       each of the projected parameter draws if `integrated = FALSE` and
#'       averaged across the projected parameter draws if `integrated = TRUE`).
#'
#'       In case of (i) the traditional projection, (ii) the latent projection
#'       with `transform = FALSE`, or (iii) the latent projection with
#'       `transform = TRUE` and `<refmodel>$family$cats` (where `<refmodel>` is
#'       an object resulting from [init_refmodel()]; see also
#'       [extend_family()]'s argument `latent_y_unqs`) being `NULL`, both
#'       elements are \eqn{S_{\mathrm{prj}} \times N}{S_prj x N} matrices. In
#'       case of (i) the augmented-data projection or (ii) the latent projection
#'       with `transform = TRUE` and `<refmodel>$family$cats` being not `NULL`,
#'       `pred` is an \eqn{S_{\mathrm{prj}} \times N \times C}{S_prj x N x C}
#'       array and `lpd` is an \eqn{S_{\mathrm{prj}} \times N}{S_prj x N}
#'       matrix.
#'   * [proj_predict()] returns an \eqn{S_{\mathrm{prj}} \times N}{S_prj x N}
#'   matrix of predictions where \eqn{S_{\mathrm{prj}}}{S_prj} denotes
#'   `nresample_clusters` in case of clustered projection. In case of (i) the
#'   augmented-data projection or (ii) the latent projection with `resp_oscale =
#'   TRUE` and `<refmodel>$family$cats` being not `NULL`, this matrix has an
#'   attribute called `cats` (the character vector of response categories) and
#'   the values of the matrix are the predicted indices of the response
#'   categories (these indices refer to the order of the response categories
#'   from attribute `cats`).
#'
#'   If the prediction is done for more than one submodel, the output from above
#'   is returned for each submodel, giving a named `list` with one element for
#'   each submodel (the names of this `list` being the numbers of solution terms
#'   of the submodels when counting the intercept, too).
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

# Function definitions ----------------------------------------------------

## The 'helper' for proj_linpred and proj_predict, ie. does all the
## functionality that is common to them. It essentially checks all the arguments
## and sets them to their respective defaults and then loops over the
## projections. For each projection, it evaluates the fun-function, which
## calculates the linear predictor if called from proj_linpred and samples from
## the predictive distribution if called from proj_predict.
proj_helper <- function(object, newdata, offsetnew, weightsnew, onesub_fun,
                        filter_nterms = NULL, ...) {
  if (inherits(object, "projection") || is_proj_list(object)) {
    if (!is.null(filter_nterms)) {
      if (!is_proj_list(object)) {
        object <- list(object)
      }
      projs <- Filter(
        function(x) {
          count_terms_chosen(x$solution_terms) %in% (filter_nterms + 1)
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

  if (!is_proj_list(projs)) {
    projs <- list(projs)
  }

  if (is.null(newdata)) {
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
    count_terms_chosen(proj$solution_terms)
  })

  preds <- lapply(projs, function(proj) {
    w_o <- proj$refmodel$extract_model_data(
      proj$refmodel$fit, newdata = newdata, wrhs = weightsnew, orhs = offsetnew,
      extract_y = FALSE
    )
    weightsnew <- w_o$weights
    offsetnew <- w_o$offset
    if (length(weightsnew) == 0) {
      weightsnew <- rep(1, nrow(newdata) %||% proj$refmodel$nobs)
    }
    if (length(offsetnew) == 0) {
      offsetnew <- rep(0, nrow(newdata) %||% proj$refmodel$nobs)
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

  return(unlist_proj(preds))
}

#' @rdname pred-projection
#' @export
proj_linpred <- function(object, newdata = NULL, offsetnew = NULL,
                         weightsnew = NULL, filter_nterms = NULL,
                         transform = FALSE, integrated = FALSE, .seed = NA,
                         ...) {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  if (!is.na(.seed)) {
    # Set seed, but ensure the old RNG state is restored on exit:
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(.seed)
  }

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
  pred_sub <- proj$refmodel$family$mu_fun(proj$outdmin, newdata = newdata,
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
    if (is.null(newdata)) {
      newdata_lat <- newdata
      if (inherits(proj$refmodel$fit, "stanreg") &&
          length(proj$refmodel$fit$offset) > 0) {
        newdata_lat <- proj$refmodel$fetch_data()
        newdata_lat$projpred_internal_offs_stanreg <- offset
      }
      ynew <- rowMeans(proj$refmodel$ref_predfun(
        fit = proj$refmodel$fit,
        newdata = newdata_lat,
        excl_offs = FALSE,
        mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)
      ))
    } else {
      ynew <- eval_lhs(formula = proj$refmodel$formula, data = newdata)
    }
  }
  lpd_out <- compute_lpd(ynew = ynew, pred_sub = pred_sub, proj = proj,
                         weights = weights, transformed = transform)
  if (integrated) {
    if (proj$refmodel$family$for_latent && transform &&
        length(dim(pred_sub)) == 3) {
      pred_sub <- arr2augmat(pred_sub, margin_draws = 1)
    }
    ## average over the projected draws
    if (proj$refmodel$family$for_latent && transform &&
        !inherits(pred_sub, "augmat")) {
      pred_sub <- proj$wdraws_prj %*% pred_sub
    } else {
      pred_sub <- structure(pred_sub %*% proj$wdraws_prj,
                            nobs_orig = attr(pred_sub, "nobs_orig"),
                            class = oldClass(pred_sub))
    }
    if (!is.null(lpd_out)) {
      if (!(proj$refmodel$family$for_latent && transform)) {
        marg_obs <- 1
      } else {
        marg_obs <- 2
      }
      lpd_out <- as.matrix(
        apply(lpd_out, marg_obs, log_weighted_mean_exp, proj$wdraws_prj)
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
       (proj$refmodel$family$for_latent && integrated) ||
       (proj$refmodel$family$for_latent && !transform))) {
    lpd_out <- t(lpd_out)
  }
  return(nlist(pred = pred_sub, lpd = lpd_out))
}

compute_lpd <- function(ynew, pred_sub, proj, weights, transformed) {
  if (!is.null(ynew)) {
    ## compute also the log-density
    target <- get_standard_y(ynew, weights, proj$refmodel$family)
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
      ll_oscale_out <- proj$refmodel$family$latent_ll_oscale(
        pred_sub, y_oscale = ynew, wobs = weights, cl_ref = proj$cl_ref,
        wdraws_ref = proj$wdraws_ref
      )
      if (!is.matrix(ll_oscale_out)) {
        stop("Unexpected structure for the output of `latent_ll_oscale`.")
      }
      if (all(is.na(ll_oscale_out))) {
        message(
          "`latent_ll_oscale` returned only `NA`s, so the corresponding ",
          "output will also be `NA` as long as `transform = TRUE`."
        )
      }
      return(ll_oscale_out)
    } else {
      if (proj$refmodel$family$for_latent && all(is.na(proj$refmodel$dis))) {
        message(
          "Cannot calculate LPD values if `transform = FALSE` and ",
          "`<refmodel>$dis` consists of only `NA`s. If it's not possible to ",
          "supply a suitable argument `dis` to init_refmodel(), consider ",
          "switching to `transform = TRUE` (which might require the ",
          "specification of functions needed by extend_family())."
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
                         nresample_clusters = 1000, .seed = NA,
                         resp_oscale = TRUE, ...) {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  if (!is.na(.seed)) {
    # Set seed, but ensure the old RNG state is restored on exit:
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(.seed)
  }

  ## proj_helper lapplies fun to each projection in object
  proj_helper(
    object = object, newdata = newdata,
    offsetnew = offsetnew, weightsnew = weightsnew,
    onesub_fun = proj_predict_aux, filter_nterms = filter_nterms,
    nresample_clusters = nresample_clusters, resp_oscale = resp_oscale, ...
  )
}

## function applied to each projected submodel in case of proj_predict()
proj_predict_aux <- function(proj, newdata, offset, weights,
                             nresample_clusters = 1000, resp_oscale = TRUE,
                             ...) {
  if (!proj$refmodel$family$for_latent && !resp_oscale) {
    stop("`resp_oscale = FALSE` can only be used in case of the latent ",
         "projection.")
  }
  mu <- proj$refmodel$family$mu_fun(proj$outdmin,
                                    newdata = newdata,
                                    offset = offset)
  if (!proj$const_wdraws_prj) {
    # In this case, the posterior draws have nonconstant weights.
    draw_inds <- sample(x = seq_along(proj$wdraws_prj),
                        size = nresample_clusters, replace = TRUE,
                        prob = proj$wdraws_prj)
  } else {
    draw_inds <- seq_along(proj$wdraws_prj)
  }
  cats_aug <- proj$refmodel$family$cats
  if (proj$refmodel$family$for_latent && resp_oscale) {
    mu_oscale <- proj$refmodel$family$latent_ilink(t(mu), cl_ref = proj$cl_ref,
                                                   wdraws_ref = proj$wdraws_ref)
    if (length(dim(mu_oscale)) < 2) {
      stop("Unexpected structure for the output of `latent_ilink`.")
    }
    if (length(dim(mu_oscale)) == 3) {
      mu_oscale_resamp <- mu_oscale[draw_inds, , , drop = FALSE]
    } else {
      mu_oscale_resamp <- mu_oscale[draw_inds, , drop = FALSE]
    }
    pppd_out <- proj$refmodel$family$latent_ppd_oscale(
      mu_oscale_resamp, wobs = weights, cl_ref = proj$cl_ref,
      wdraws_ref = proj$wdraws_ref, idxs_prjdraws = draw_inds
    )
    if (!is.matrix(pppd_out)) {
      stop("Unexpected structure for the output of `latent_ppd_oscale`.")
    }
    if (all(is.na(mu_oscale))) {
      message(
        "`latent_ilink` returned only `NA`s, so the output will also be ",
        "`NA` as long as `resp_oscale = TRUE`."
      )
    } else if (all(is.na(pppd_out))) {
      message(
        "`latent_ppd_oscale` returned only `NA`s, so the output will also be ",
        "`NA` as long as `resp_oscale = TRUE`."
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
        message(
          "Cannot draw from the latent Gaussian distribution if ",
          "`<refmodel>$dis` consists of only `NA`s. If it's not possible to ",
          "supply a suitable argument `dis` to init_refmodel(), consider ",
          "switching to `resp_oscale = TRUE` (which might require the ",
          "specification of functions needed by extend_family())."
        )
      }
    }
    pppd_out <- do.call(rbind, lapply(draw_inds, function(i) {
      proj$refmodel$family$ppd(mu[, i], proj$dis[i], weights)
    }))
  }
  return(structure(pppd_out, cats = cats_aug))
}

#' Plot predictive performance
#'
#' This is the [plot()] method for `vsel` objects (returned by [varsel()] or
#' [cv_varsel()]). It visualizes the predictive performance of the reference
#' model (possibly also that of some other "baseline" model) and that of the
#' submodels along the full-data predictor ranking. Basic information about the
#' (CV) variability in the ranking of the predictors is included as well (if
#' available; inferred from [cv_proportions()]). For a tabular representation,
#' see [summary.vsel()].
#'
#' @inheritParams summary.vsel
#' @param x An object of class `vsel` (returned by [varsel()] or [cv_varsel()]).
#' @param thres_elpd Only relevant if `any(stats %in% c("elpd", "mlpd"))`. The
#'   threshold for the ELPD difference (taking the submodel's ELPD minus the
#'   baseline model's ELPD) above which the submodel's ELPD is considered to be
#'   close enough to the baseline model's ELPD. An equivalent rule is applied in
#'   case of the MLPD. See [suggest_size()] for a formalization. Supplying `NA`
#'   deactivates this.
#' @param ranking_nterms_max Maximum submodel size (number of predictor terms)
#'   for which the predictor names and the corresponding ranking proportions are
#'   added on the x-axis. Using `NULL` is effectively the same as using
#'   `nterms_max`. Using `NA` causes the predictor names and the corresponding
#'   ranking proportions to be omitted. Note that `ranking_nterms_max` does not
#'   count the intercept, so `ranking_nterms_max = 1` corresponds to the
#'   submodel consisting of the first (non-intercept) predictor term.
#' @param ranking_abbreviate A single logical value indicating whether the
#'   predictor names in the full-data predictor ranking should be abbreviated by
#'   [abbreviate()] (`TRUE`) or not (`FALSE`). See also argument
#'   `ranking_abbreviate_args` and section "Value".
#' @param ranking_abbreviate_args A `list` of arguments (except for `names.arg`)
#'   to be passed to [abbreviate()] in case of `ranking_abbreviate = TRUE`.
#' @param ranking_repel Either `NULL`, `"text"`, or `"label"`. By `NULL`, the
#'   full-data predictor ranking and the corresponding ranking proportions are
#'   placed below the x-axis. By `"text"` or `"label"`, they are placed within
#'   the plotting area, using [ggrepel::geom_text_repel()] or
#'   [ggrepel::geom_label_repel()], respectively. See also argument
#'   `ranking_repel_args`.
#' @param ranking_repel_args A `list` of arguments (except for `mapping`) to be
#'   passed to [ggrepel::geom_text_repel()] or [ggrepel::geom_label_repel()] in
#'   case of `ranking_repel = "text"` or `ranking_repel = "label"`,
#'   respectively.
#' @param ranking_colored A single logical value indicating whether the points
#'   and the uncertainty bars should be gradient-colored according to the CV
#'   ranking proportions (`TRUE`) or not (`FALSE`). The CV ranking proportions
#'   may be cumulated (see argument `cumulate`). Note that the point and the
#'   uncertainty bar at submodel size 0 (i.e., at the intercept-only model) are
#'   always colored in gray because the intercept is forced to be selected
#'   before any predictors are selected (in other words, the reason is that for
#'   submodel size 0, the question of variability across CV folds is not
#'   appropriate in the first place).
#' @param cumulate Passed to argument `cumulate` of [cv_proportions()]. Affects
#'   the ranking proportions given on the x-axis (below the full-data predictor
#'   ranking).
#' @param text_angle Passed to argument `angle` of [ggplot2::element_text()] for
#'   the x-axis tick labels. In case of long predictor names (and/or large
#'   `nterms_max`), `text_angle = 45` might be helpful (for example).
#'
#' @inherit summary.vsel details
#'
#' @return A \pkg{ggplot2} plotting object (of class `gg` and `ggplot`). If
#'   `ranking_abbreviate` is `TRUE`, the output of [abbreviate()] is stored in
#'   an attribute called `projpred_ranking_abbreviated` (to allow the
#'   abbreviations to be easily mapped back to the original predictor names).
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
#'   # Run varsel() (here without cross-validation and with small values for
#'   # `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the sake of
#'   # speed in this example; this is not recommended in general):
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
    alpha = 2 * pnorm(-1),
    baseline = if (!inherits(x$refmodel, "datafit")) "ref" else "best",
    thres_elpd = NA,
    resp_oscale = TRUE,
    ranking_nterms_max = NULL,
    ranking_abbreviate = FALSE,
    ranking_abbreviate_args = list(),
    ranking_repel = NULL,
    ranking_repel_args = list(),
    ranking_colored = FALSE,
    cumulate = FALSE,
    text_angle = NULL,
    ...
) {
  object <- x
  validate_vsel_object_stats(object, stats, resp_oscale = resp_oscale)
  baseline <- validate_baseline(object$refmodel, baseline, deltas)
  if (!is.null(ranking_repel) && !requireNamespace("ggrepel", quietly = TRUE)) {
    warning("Package 'ggrepel' is needed for a non-`NULL` argument ",
            "`ranking_repel`, but could not be found. Setting `ranking_repel` ",
            "to `NULL` now.")
    ranking_repel <- NULL
  } else if (!is.null(ranking_repel)) {
    stopifnot(isTRUE(ranking_repel %in% c("text", "label")))
  }

  ## compute all the statistics and fetch only those that were asked
  nfeat_baseline <- get_nfeat_baseline(object, baseline, stats[1],
                                       resp_oscale = resp_oscale)
  tab <- rbind(
    .tabulate_stats(object, stats, alpha = alpha,
                    nfeat_baseline = nfeat_baseline, resp_oscale = resp_oscale,
                    ...),
    .tabulate_stats(object, stats, alpha = alpha, resp_oscale = resp_oscale,
                    ...)
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
  if (!is_wholenumber(nterms_max)) {
    stop("`nterms_max` must be a whole number.")
  }
  nterms_max <- as.integer(nterms_max)
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
  if (object$refmodel$family$for_latent) {
    if (resp_oscale) {
      ylab <- paste(ylab, "(response scale)")
    } else {
      ylab <- paste(ylab, "(latent scale)")
    }
  }

  # make sure that breaks on the x-axis are integers
  n_opts <- 4:6
  n_possible <- Filter(function(x) nterms_max %% x == 0, n_opts)
  n_alt <- n_opts[which.min(n_opts - (nterms_max %% n_opts))]
  nb <- ifelse(length(n_possible) > 0, min(n_possible), n_alt)
  # Using as.integer() only to make it clear that this is an integer (just like
  # `breaks` and `minor_breaks`):
  by <- as.integer(ceiling(nterms_max / min(nterms_max, nb)))
  breaks <- seq(0L, by * min(nterms_max, nb), by)
  minor_breaks <- if (by %% 2 == 0) {
    seq(by %/% 2L, by * min(nterms_max, nb), by)
  } else {
    NULL
  }
  if (is.null(ranking_nterms_max)) {
    ranking_nterms_max <- nterms_max
  } else if (!is.na(ranking_nterms_max)) {
    ranking_nterms_max <- min(ranking_nterms_max, nterms_max)
    if (!is_wholenumber(ranking_nterms_max)) {
      stop("`ranking_nterms_max` must be a whole number.")
    }
    ranking_nterms_max <- as.integer(ranking_nterms_max)
  }
  if (!is.na(ranking_nterms_max)) {
    breaks <- sort(union(breaks, seq_len(ranking_nterms_max)))
    minor_breaks <- setdiff(minor_breaks, breaks)
  }

  if (!is.na(thres_elpd)) {
    # Table of thresholds used in extended suggest_size() heuristics (only in
    # case of ELPD and MLPD):
    thres_tab_basic <- data.frame(
      statistic = c("elpd", "mlpd"),
      thres = c(thres_elpd, thres_elpd / object$nobs_test)
    )
  }

  # Start x-axis label (title):
  xlab <- "Submodel size (number of predictor terms)"

  if (!is.na(ranking_nterms_max)) {
    # Predictor ranking(s):
    rk <- ranking(object, nterms_max = ranking_nterms_max)
    if (!is.null(rk[["foldwise"]])) {
      pr_rk <- diag(cv_proportions(rk, cumulate = cumulate))
    } else {
      pr_rk <- rep(NA, length(rk[["fulldata"]]))
    }
    rk_dfr <- data.frame(
      size = c(0L, seq_along(rk[["fulldata"]])),
      rk_fulldata = c("", rk[["fulldata"]]),
      cv_props_diag = c(NA, pr_rk)
    )
    rk_dfr[["cv_props_diag_num"]] <- rk_dfr[["cv_props_diag"]]
    rk_dfr[["cv_props_diag"]] <- paste(round(100 * rk_dfr[["cv_props_diag"]]),
                                       "%")
    rk_dfr[["cv_props_diag"]][1] <- "" # empty model
    rk_dfr_empty <- do.call(rbind, lapply(
      setdiff(breaks, rk_dfr[["size"]]),
      function(br_j) {
        data.frame(size = br_j, rk_fulldata = "", cv_props_diag = "",
                   cv_props_diag_num = NA)
      }
    ))
    rk_dfr <- rbind(rk_dfr, rk_dfr_empty)
    if (ranking_abbreviate) {
      rk_fulldata_abbv <- do.call(abbreviate, c(
        list(names.arg = rk_dfr[["rk_fulldata"]]),
        ranking_abbreviate_args
      ))
      rk_dfr[["rk_fulldata"]] <- rk_fulldata_abbv
    }
    rk_dfr[["rkfulldt_cvpropdiag"]] <- rk_dfr[["rk_fulldata"]]
    if (!is.null(rk[["foldwise"]])) {
      rk_dfr[["rkfulldt_cvpropdiag"]] <- paste(rk_dfr[["rkfulldt_cvpropdiag"]],
                                               rk_dfr[["cv_props_diag"]],
                                               sep = "\n")
    }
    rk_dfr[["size_rkfulldt_cvpropdiag"]] <- paste(
      rk_dfr[["size"]], rk_dfr[["rkfulldt_cvpropdiag"]], sep = "\n"
    )

    # Continue x-axis label (title):
    xlab_rk <- "Corresponding predictor from full-data predictor ranking"
    if (identical(ranking_repel, "text")) {
      xlab_rk <- paste("Text:", xlab_rk)
    } else if (identical(ranking_repel, "label")) {
      xlab_rk <- paste("Label:", xlab_rk)
    }
    xlab <- paste(xlab, xlab_rk, sep = "\n")
    if (!is.null(rk[["foldwise"]])) {
      if (cumulate) {
        cumul_pretty <- " cumulated "
      } else {
        cumul_pretty <- " "
      }
      xlab_cumul <- paste0("Corresponding main diagonal element from",
                           cumul_pretty, "CV ranking proportions matrix")
      if (identical(ranking_repel, "text")) {
        xlab_cumul <- paste("Text:", xlab_cumul)
      } else if (identical(ranking_repel, "label")) {
        xlab_cumul <- paste("Label:", xlab_cumul)
      }
      xlab <- paste(xlab, xlab_cumul, sep = "\n")
    }
  }

  # plot submodel results
  data_gg <- subset(stats_sub, stats_sub$size <= nterms_max)
  if (!is.na(ranking_nterms_max) &&
      (!is.null(ranking_repel) ||
       (ranking_colored && !is.null(rk[["foldwise"]])))) {
    colnms_orig <- names(data_gg)
    data_gg[["row_idx"]] <- seq_len(nrow(data_gg))
    cols_add <- c("cv_props_diag_num", "rkfulldt_cvpropdiag")
    data_gg <- merge(data_gg,
                     rk_dfr[, c("size", cols_add), drop = FALSE],
                     by = "size", all.x = TRUE, all.y = FALSE, sort = FALSE)
    data_gg <- data_gg[order(data_gg[["row_idx"]]), , drop = FALSE]
    data_gg[["row_idx"]] <- NULL
    data_gg <- data_gg[, c(colnms_orig, cols_add), drop = FALSE]
  }
  pp <- ggplot(data = data_gg,
               mapping = aes(x = .data[["size"]], y = .data[["value"]],
                             ymin = .data[["lq"]], ymax = .data[["uq"]]))
  if (!all(is.na(stats_ref$se))) {
    # add reference model results if they exist

    pp <- pp +
      # The reference model's dashed red horizontal line:
      geom_hline(aes(yintercept = .data[["value"]]),
                 data = stats_ref,
                 color = "darkred", linetype = 2)

    if (!is.na(thres_elpd)) {
      # The thresholds used in extended suggest_size() heuristics:
      thres_tab_ref <- merge(thres_tab_basic,
                             stats_ref[, c("statistic", "value")],
                             by = "statistic")
      thres_tab_ref$thres <- thres_tab_ref$value + thres_tab_ref$thres
      pp <- pp +
        geom_hline(aes(yintercept = .data[["thres"]]),
                   data = thres_tab_ref,
                   color = "gray50", linetype = "dotdash")
    }
  }
  if (baseline != "ref") {
    # add baseline model results (if different from the reference model)

    pp <- pp +
      # The baseline model's dotted black horizontal line:
      geom_hline(aes(yintercept = .data[["value"]]),
                 data = stats_bs,
                 color = "black", linetype = 3)

    if (!is.na(thres_elpd)) {
      # The thresholds used in extended suggest_size() heuristics:
      thres_tab_bs <- merge(thres_tab_basic,
                            stats_bs[, c("statistic", "value")],
                            by = "statistic")
      thres_tab_bs$thres <- thres_tab_bs$value + thres_tab_bs$thres
      pp <- pp +
        geom_hline(aes(yintercept = .data[["thres"]]),
                   data = thres_tab_bs,
                   color = "darkgreen", linetype = "longdash")
    }
  }
  if (!is.na(ranking_nterms_max) && ranking_colored &&
      !is.null(rk[["foldwise"]])) {
    aes_linerg_pt <- aes(color = .data[["cv_props_diag_num"]])
    alpha_linerg <- 1
  } else {
    aes_linerg_pt <- NULL
    alpha_linerg <- 0.55
  }
  if (!is.na(ranking_nterms_max) && is.null(ranking_repel)) {
    tick_labs_x <- rk_dfr[order(match(rk_dfr[["size"]], breaks), na.last = NA),
                          "size_rkfulldt_cvpropdiag"]
  } else {
    tick_labs_x <- waiver()
  }
  # The submodel-specific graphical elements:
  pp <- pp +
    geom_linerange(aes_linerg_pt, alpha = alpha_linerg, linewidth = 1) +
    geom_line() +
    geom_point(aes_linerg_pt, size = 3)
  # Miscellaneous stuff (axes, theming, faceting, etc.):
  if (!is.na(ranking_nterms_max) && ranking_colored &&
      !is.null(rk[["foldwise"]])) {
    ### Option 1:
    pp <- pp +
      scale_color_gradient(name = "Proportion\nof CV folds",
                           labels = scales::label_percent(suffix = " %"),
                           limits = c(0, 1),
                           low = "#ededed", high = "#0f365c")
    ###
    ### Option 2 (requires the 'RColorBrewer' package):
    # pp <- pp +
    #   scale_color_distiller(name = "Proportion\nof CV folds",
    #                        labels = scales::label_percent(suffix = " %"),
    #                        direction = 1)
    ###
  }
  pp <- pp +
    scale_x_continuous(breaks = breaks, minor_breaks = minor_breaks,
                       limits = c(min(breaks), max(breaks)),
                       labels = tick_labs_x) +
    labs(x = xlab, y = ylab) +
    theme(axis.text.x = element_text(angle = text_angle, hjust = 0.5,
                                     vjust = 0.5)) +
    facet_grid(statistic ~ ., scales = "free_y")
  if (!is.na(ranking_nterms_max) && !is.null(ranking_repel)) {
    if (identical(ranking_repel, "text")) {
      geom_repel_fun <- ggrepel::geom_text_repel
    } else if (identical(ranking_repel, "label")) {
      geom_repel_fun <- ggrepel::geom_label_repel
    }
    pp <- pp +
      do.call(geom_repel_fun, c(
        list(mapping = aes(label = .data[["rkfulldt_cvpropdiag"]])),
        ranking_repel_args
      ))
  }
  if (!is.na(ranking_nterms_max) && ranking_abbreviate) {
    attr(pp, "projpred_ranking_abbreviated") <- rk_fulldata_abbv[
      rk_fulldata_abbv != ""
    ]
  }
  return(pp)
}

#' Summary of a [varsel()] or [cv_varsel()] run
#'
#' This is the [summary()] method for `vsel` objects (returned by [varsel()] or
#' [cv_varsel()]). Apart from some general information about the [varsel()] or
#' [cv_varsel()] run, it shows the full-data predictor ranking, basic
#' information about the (CV) variability in the ranking of the predictors (if
#' available; inferred from [cv_proportions()]), and estimates for
#' user-specified predictive performance statistics. For a graphical
#' representation, see [plot.vsel()].
#'
#' @param object An object of class `vsel` (returned by [varsel()] or
#'   [cv_varsel()]).
#' @param nterms_max Maximum submodel size (number of predictor terms) for which
#'   the performance statistics are calculated. Using `NULL` is effectively the
#'   same as `length(ranking(object)[["fulldata"]])`. Note that `nterms_max`
#'   does not count the intercept, so use `nterms_max = 0` for the
#'   intercept-only model. For [plot.vsel()], `nterms_max` must be at least `1`.
#' @param stats One or more character strings determining which performance
#'   statistics (i.e., utilities or losses) to estimate based on the
#'   observations in the evaluation (or "test") set (in case of
#'   cross-validation, these are all observations because they are partitioned
#'   into multiple test sets; in case of [varsel()] with `d_test = NULL`, these
#'   are again all observations because the test set is the same as the training
#'   set). Available statistics are:
#'   * `"elpd"`: expected log (pointwise) predictive density (for a new
#'   dataset). Estimated by the sum of the observation-specific log predictive
#'   density values (with each of these predictive density values being
#'   a---possibly weighted---average across the parameter draws).
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
#'   from the baseline model (see argument `baseline`). With a "difference
#'   *from* the baseline model", we mean to take the submodel statistic minus
#'   the baseline model statistic (not the other way round).
#' @param alpha A number determining the (nominal) coverage `1 - alpha` of the
#'   normal-approximation (or bootstrap; see argument `stats`) confidence
#'   intervals. For example, in case of the normal approximation, `alpha = 2 *
#'   pnorm(-1)` corresponds to a confidence interval stretching by one standard
#'   error on either side of the point estimate.
#' @param baseline For [summary.vsel()]: Only relevant if `deltas` is `TRUE`.
#'   For [plot.vsel()]: Always relevant. Either `"ref"` or `"best"`, indicating
#'   whether the baseline is the reference model or the best submodel found (in
#'   terms of `stats[1]`), respectively.
#' @param resp_oscale Only relevant for the latent projection. A single logical
#'   value indicating whether to calculate the performance statistics on the
#'   original response scale (`TRUE`) or on latent scale (`FALSE`).
#' @param cumulate Passed to argument `cumulate` of [cv_proportions()]. Affects
#'   column `cv_proportions_diag` of the summary table.
#' @param ... Arguments passed to the internal function which is used for
#'   bootstrapping (if applicable; see argument `stats`). Currently, relevant
#'   arguments are `B` (the number of bootstrap samples, defaulting to `2000`)
#'   and `seed` (see [set.seed()], but defaulting to `NA` so that [set.seed()]
#'   is not called within that function at all).
#'
#' @details The `stats` options `"mse"` and `"rmse"` are only available for:
#'   * the traditional projection,
#'   * the latent projection with `resp_oscale = FALSE`,
#'   * the latent projection with `resp_oscale = TRUE` in combination with
#'   `<refmodel>$family$cats` being `NULL`.
#'
#'   The `stats` option `"acc"` (= `"pctcorr"`) is only available for:
#'   * the [binomial()] family in case of the traditional projection,
#'   * all families in case of the augmented-data projection,
#'   * the [binomial()] family (on the original response scale) in case of the
#'   latent projection with `resp_oscale = TRUE` in combination with
#'   `<refmodel>$family$cats` being `NULL`,
#'   * all families (on the original response scale) in case of the latent
#'   projection with `resp_oscale = TRUE` in combination with
#'   `<refmodel>$family$cats` being not `NULL`.
#'
#'   The `stats` option `"auc"` is only available for:
#'   * the [binomial()] family in case of the traditional projection,
#'   * the [binomial()] family (on the original response scale) in case of the
#'   latent projection with `resp_oscale = TRUE` in combination with
#'   `<refmodel>$family$cats` being `NULL`.
#'
#' @return An object of class `vselsummary`.
#'
#' @seealso [print.vselsummary()]
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
#'   # Run varsel() (here without cross-validation and with small values for
#'   # `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the sake of
#'   # speed in this example; this is not recommended in general):
#'   vs <- varsel(fit, nterms_max = 3, nclusters = 5, nclusters_pred = 10,
#'                seed = 5555)
#'   print(summary(vs), digits = 1)
#' }
#'
#' @export
summary.vsel <- function(
    object,
    nterms_max = NULL,
    stats = "elpd",
    type = c("mean", "se", "diff", "diff.se"),
    deltas = FALSE,
    alpha = 2 * pnorm(-1),
    baseline = if (!inherits(object$refmodel, "datafit")) "ref" else "best",
    resp_oscale = TRUE,
    cumulate = FALSE,
    ...
) {
  validate_vsel_object_stats(object, stats, resp_oscale = resp_oscale)
  baseline <- validate_baseline(object$refmodel, baseline, deltas)

  # Initialize output:
  out <- c(
    object$refmodel[c("formula", "family")],
    object[c("nobs_train", "type_test", "nobs_test", "method", "cv_method", "K",
             "validate_search", "clust_used_search", "clust_used_eval",
             "nprjdraws_search", "nprjdraws_eval")]
  )
  if (isTRUE(out$validate_search)) {
    out$search_included <- "search included (i.e., fold-wise searches)"
  } else {
    out$search_included <- "search not included (i.e., a full-data search only)"
  }
  class(out) <- "vselsummary"

  # The full table of the performance statistics from `stats`:
  if (deltas) {
    nfeat_baseline <- get_nfeat_baseline(object, baseline, stats[1],
                                         resp_oscale = resp_oscale)
    tab <- .tabulate_stats(object, stats, alpha = alpha,
                           nfeat_baseline = nfeat_baseline,
                           resp_oscale = resp_oscale, ...)
  } else {
    tab <- .tabulate_stats(object, stats, alpha = alpha,
                           resp_oscale = resp_oscale, ...)
  }
  stats_table <- subset(tab, tab$size != Inf)
  stats_table <- do.call(rbind,
                         lapply(split(stats_table, stats_table$statistic),
                                utils::head,
                                n = length(object$solution_terms) + 1))
  row.names(stats_table) <- NULL

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

  # Predictor ranking(s) and associated ranking proportions from fold-wise
  # predictor rankings (if existing):
  rk <- ranking(object)
  if (!is.null(rk[["foldwise"]]) && ncol(rk[["foldwise"]]) > 0) {
    pr_rk <- diag(cv_proportions(rk, cumulate = cumulate))
  } else {
    pr_rk <- rep(NA, length(rk[["fulldata"]]))
  }

  # Construct the (almost) final output table by looping over all requested
  # statistics, reshaping the corresponding data in `stats_table`, and selecting
  # only the requested `type`s:
  arr <- data.frame(size = unique(stats_table$size),
                    solution_terms = c(NA_character_, rk[["fulldata"]]),
                    cv_proportions_diag = c(NA, pr_rk))
  for (i in seq_along(stats)) {
    temp <- subset(stats_table, stats_table$statistic == stats[i], qty)
    newnames <- suffix[[i]]
    colnames(temp) <- newnames
    arr <- cbind(arr, temp)
  }
  row.names(arr) <- NULL

  # Output (and also cut `arr` at `nterms_max` (if provided)):
  if (is.null(nterms_max)) {
    nterms_max <- max(stats_table$size)
  }
  out$nterms <- nterms_max
  out$selection <- subset(arr, arr$size <= nterms_max)
  out$resp_oscale <- resp_oscale
  out$deltas <- deltas
  out$cumulate <- cumulate
  return(out)
}

#' Print summary of a [varsel()] or [cv_varsel()] run
#'
#' This is the [print()] method for summary objects created by [summary.vsel()].
#' It displays a summary of the results from a [varsel()] or [cv_varsel()] run.
#'
#' @param x An object of class `vselsummary`.
#' @param ... Arguments passed to [print.data.frame()].
#'
#' @details In the table printed at the bottom, column `solution_terms` contains
#'   the full-data predictor ranking and column `cv_proportions_diag` contains
#'   the main diagonal of the matrix returned by [cv_proportions()] (with
#'   `cumulate` as set in the [summary.vsel()] call that created `x`).
#'
#' @return The output of [summary.vsel()] (invisible).
#'
#' @export
print.vselsummary <- function(x, ...) {
  if (x$family$for_latent) {
    cat("------\nResponse-scale family:\n")
    print(structure(x$family[c("family_oscale", "link_oscale")],
                    class = "family"))
    cat("------\nLatent-scale family:\n")
  }
  print(x$family)
  if (x$family$for_latent) {
    cat("------\n")
  }
  cat("Formula: ")
  print(x$formula, showEnv = FALSE)
  if (x$type_test != "test_hold-out") {
    cat("Observations: ", x$nobs_train, "\n", sep = "")
  } else {
    cat("Observations (training set): ", x$nobs_train, "\n", sep = "")
    cat("Observations (test set): ", x$nobs_test, "\n", sep = "")
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
    cv_meth_pretty <- sub("^kfold$", "K-fold", x$cv_method)
    cv_meth_pretty <- sub("^LOO$", "PSIS-LOO", cv_meth_pretty)
    if (x$cv_method == "kfold") {
      K_pretty <- paste("K =", x$K, "and ")
    } else {
      K_pretty <- ""
    }
    cat("CV method: ", cv_meth_pretty, " CV with ", K_pretty, x$search_included,
        "\n", sep = "")
  }
  cat("Search method: ", x$method, "\n", sep = "")
  cat("Maximum submodel size for the search: ", x$nterms, "\n", sep = "")
  if (x$clust_used_search) {
    clust_search_pretty <- " (from clustered projection)"
  } else {
    clust_search_pretty <- ""
  }
  if (x$clust_used_eval) {
    clust_eval_pretty <- " (from clustered projection)"
  } else {
    clust_eval_pretty <- ""
  }
  cat("Number of projected draws in the search: ", x$nprjdraws_search,
      clust_search_pretty, "\n", sep = "")
  cat("Number of projected draws in the performance evaluation: ",
      x$nprjdraws_eval, clust_eval_pretty, "\n", sep = "")
  cat("\n")
  if (x$family$for_latent) {
    if (x$resp_oscale) {
      scale_string <- " (response scale)"
    } else {
      scale_string <- " (latent scale)"
    }
  } else {
    scale_string <- ""
  }
  cat("Performance evaluation summary", scale_string, " with `deltas = ",
      x$deltas, "` and `cumulate = ", x$cumulate, "`:\n", sep = "")
  print(x$selection, row.names = FALSE, ...)
  if (isTRUE(x$validate_search)) {
    message(
      "Column `solution_terms` contains the full-data predictor ranking. To ",
      "retrieve the fold-wise predictor rankings, use the ranking() function, ",
      "possibly followed by cv_proportions() for computing the ranking ",
      "proportions (which can be visualized by plot.cv_proportions()). The ",
      "main diagonal of the matrix returned by cv_proportions() (with ",
      "`cumulate = ", x$cumulate, "`) is contained in column ",
      "`cv_proportions_diag`."
    )
  }
  return(invisible(x))
}

#' Print results (summary) of a [varsel()] or [cv_varsel()] run
#'
#' This is the [print()] method for `vsel` objects (returned by [varsel()] or
#' [cv_varsel()]). It displays a summary of a [varsel()] or [cv_varsel()] run by
#' first calling [summary.vsel()] and then [print.vselsummary()].
#'
#' @param x An object of class `vsel` (returned by [varsel()] or [cv_varsel()]).
#' @param ... Arguments passed to [summary.vsel()] (apart from argument `digits`
#'   which is passed to [print.vselsummary()]).
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
#' @param pct A number giving the proportion (*not* percents) of the *relative*
#'   null model utility one is willing to sacrifice. See section "Details" below
#'   for more information.
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
#'   size is the smallest model size \eqn{j \in \{0, 1, ...,
#'   \texttt{nterms\_max}\}}{{j = 0, 1, ..., nterms_max}} for which either the
#'   lower or upper bound (depending on argument `type`) of the
#'   normal-approximation (or bootstrap; see argument `stat`) confidence
#'   interval (with nominal coverage `1 - alpha`; see argument `alpha` of
#'   [summary.vsel()]) for \eqn{U_j - U_{\mathrm{base}}}{U_j - U_base} (with
#'   \eqn{U_j} denoting the \eqn{j}-th submodel's true utility and
#'   \eqn{U_{\mathrm{base}}}{U_base} denoting the baseline model's true utility)
#'   falls above (or is equal to) \deqn{\texttt{pct} \cdot (u_0 -
#'   u_{\mathrm{base}})}{pct * (u_0 - u_base)} where \eqn{u_0} denotes the null
#'   model's estimated utility and \eqn{u_{\mathrm{base}}}{u_base} the baseline
#'   model's estimated utility. The baseline model is either the reference model
#'   or the best submodel found (see argument `baseline` of [summary.vsel()]).
#'
#'   If `!is.na(thres_elpd)` and `stat = "elpd"`, the decision rule above is
#'   extended: The suggested model size is then the smallest model size \eqn{j}
#'   fulfilling the rule above *or* \eqn{u_j - u_{\mathrm{base}} >
#'   \texttt{thres\_elpd}}{u_j - u_base > thres_elpd}. Correspondingly, in case
#'   of `stat = "mlpd"` (and `!is.na(thres_elpd)`), the suggested model size is
#'   the smallest model size \eqn{j} fulfilling the rule above *or* \eqn{u_j -
#'   u_{\mathrm{base}} > \frac{\texttt{thres\_elpd}}{N}}{u_j - u_base >
#'   thres_elpd / N} with \eqn{N} denoting the number of observations.
#'
#'   For example (disregarding the special extensions in case of
#'   `!is.na(thres_elpd)` with `stat = "elpd"` or `stat = "mlpd"`),
#'   `alpha = 2 * pnorm(-1)`, `pct = 0`, and `type = "upper"` means that we
#'   select the smallest model size for which the upper bound of the
#'   `1 - 2 * pnorm(-1)` (approximately 68.3%) confidence interval for
#'   \eqn{U_j - U_{\mathrm{base}}}{U_j - U_base} exceeds (or is equal to) zero,
#'   that is (if `stat` is a performance statistic for which the normal
#'   approximation is used, not the bootstrap), for which the submodel's utility
#'   estimate is at most one standard error smaller than the baseline model's
#'   utility estimate (with that standard error referring to the utility
#'   *difference*).
#'
#'   Apart from the two [summary.vsel()] arguments mentioned above (`alpha` and
#'   `baseline`), `resp_oscale` is another important [summary.vsel()] argument
#'   that may be passed via `...`.
#'
#' @note Loss statistics like the root mean squared error (RMSE) and the mean
#'   squared error (MSE) are converted to utilities by multiplying them by `-1`,
#'   so a call such as `suggest_size(object, stat = "rmse", type = "upper")`
#'   finds the smallest model size whose upper confidence interval bound for the
#'   *negative* RMSE or MSE exceeds the cutoff (or, equivalently, has the lower
#'   confidence interval bound for the RMSE or MSE below the cutoff). This is
#'   done to make the interpretation of argument `type` the same regardless of
#'   argument `stat`.
#'
#' @return A single numeric value, giving the suggested submodel size (or `NA`
#'   if the suggestion failed).
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
#'   # Run varsel() (here without cross-validation and with small values for
#'   # `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the sake of
#'   # speed in this example; this is not recommended in general):
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
  stats <- stats$selection

  if (is_util(stat)) {
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
  nobs_test <- object$nobs_test
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
        warning("Could not suggest submodel size. Investigate plot.vsel() to ",
                "identify if the search was terminated too early. If this is ",
                "the case, run variable selection with larger value for ",
                "`nterms_max`.")
      }
    }
  } else {
    suggested_size <- min(res)
    # We don't use `na.rm = TRUE` in min() to be as cautious as possible. In
    # fact, we could refine this to remove `NA`s after the first non-`NA` value
    # (meaning that if there is no non-`NA` value, no `NA`s will be removed),
    # but this gets overly complicated and it's better to be as cautious as
    # possible (because `NA`s after the first non-`NA` value are also strange).
  }

  return(suggested_size)
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

# Escape special characters in each element of a character vector, to give a
# character vector of the same length which may be used in regular expressions.
# Copied over from brms::escape_all() (GitHub commit
# e42e8da64fc48919085fabd6cba40b7b86668f4b) with Paul Bürkner's consent.
# Slightly refactored afterwards.
esc_chars <- function(chr_vec) {
  for (chr_spcl in c(".", "*", "+", "?", "^", "$", "(", ")", "[", "]", "|")) {
    chr_vec <- gsub(chr_spcl, paste0("\\", chr_spcl), chr_vec, fixed = TRUE)
  }
  return(chr_vec)
}

# Helper function for removing underscores in response category names (as done
# by brms) contained in a special character vector. Unfortunately, for these
# special character vectors, this replacement doesn't seem to be feasible with
# regular expressions, so we need to iterate over the elements of such a special
# character vector as well as over the category names and perform the
# replacement manually:
rm_underscore <- function(nms, nms_lats, preceding_char = ".") {
  preceding_char_esc <- esc_chars(preceding_char)
  unlist(lapply(strsplit(nms, "~"), function(nm_split) {
    paste(unlist(lapply(nm_split, function(nm_split_part) {
      for (nm_lat in grep("_", nms_lats, value = TRUE)) {
        nm_lat_esc <- esc_chars(nm_lat)
        nm_lat_regex <- paste0(preceding_char_esc, "(", nm_lat_esc, ")")
        if (grepl(nm_lat_regex, nm_split_part)) {
          nm_split_part <- paste0(sub(nm_lat_regex, "", nm_split_part),
                                  preceding_char, gsub("_", "", nm_lat))
        }
      }
      return(nm_split_part)
    })), collapse = "~")
  }))
}

# Make the parameter names for variance components adhere to the naming scheme
# `nm_scheme`:
mknms_VarCorr <- function(nms, nms_lats = NULL, nm_scheme, coef_nms) {
  if (!is.null(nms_lats)) {
    stopifnot(nm_scheme == "brms")
    # Remove underscores in the response category names (as done by brms):
    if (any(grepl("_", nms_lats))) {
      nms <- rm_underscore(nms, nms_lats = nms_lats)
      coef_nms <- lapply(coef_nms, rm_underscore, nms_lats = nms_lats,
                         preceding_char = "")
      nms_lats <- gsub("_", "", nms_lats)
    }
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
    grp_nms_esc <- paste(esc_chars(grp_nms), collapse = "|")
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
    coef_nms_i_esc <- paste(esc_chars(coef_nms_i), collapse = "|")
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
    nms_lats_esc <- paste(esc_chars(nms_lats), collapse = "|")
    # Put the string `mu` in front of the latent category names and replace the
    # following tilde by an underscore:
    nms <- gsub(paste0("(", nms_lats_esc, ")~"), "mu\\1_", nms)
  }
  return(nms)
}

# Make the parameter names for group-level effects adhere to the naming scheme
# `nm_scheme`:
mknms_ranef <- function(nms, nms_lats = NULL, nm_scheme, coef_nms) {
  if (!is.null(nms_lats)) {
    stopifnot(nm_scheme == "brms")
    # Remove underscores in the response category names (as done by brms):
    if (any(grepl("_", nms_lats))) {
      nms <- rm_underscore(nms, nms_lats = nms_lats)
      coef_nms <- lapply(coef_nms, rm_underscore, nms_lats = nms_lats,
                         preceding_char = "")
      nms_lats <- gsub("_", "", nms_lats)
    }
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
    coef_nms_i_esc <- paste(esc_chars(coef_nms_i), collapse = "|")
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
      grp_nm_i_esc <- esc_chars(grp_nm_i)
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
    nms_lats_esc <- paste(esc_chars(nms_lats), collapse = "|")
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
  nmsdf[, 1] <- paste0("mu", gsub("_", "", nmsdf[, 1]))
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
  return(replace_population_names(coef(x), ...))
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
get_subparams.glmmPQL <- function(x, ...) {
  ### TODO (glmmPQL): Implement the get_subparams.glmmPQL() method:
  stop("Under construction (the get_subparams.glmmPQL() method needs to be ",
       "implemented.")
  ###
}

#' @noRd
#' @export
get_subparams.lmerMod <- function(x, ...) {
  population_effects <- replace_population_names(lme4::fixef(x), ...)

  group_vc_raw <- lme4::VarCorr(x)
  group_vc <- proc_VarCorr(group_vc_raw,
                           coef_nms = lapply(group_vc_raw, rownames), ...)

  subparams <- c(population_effects, group_vc)

  if (!getOption("projpred.mlvl_pred_new", FALSE)) {
    group_ef <- proc_ranef(lme4::ranef(x, condVar = FALSE),
                           coef_nms = lapply(group_vc_raw, rownames), ...)
    subparams <- c(subparams, group_ef)
  }

  return(subparams)
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

  subparams <- c(thres, replace_population_names(x$beta, ...), group_vc)

  if (!getOption("projpred.mlvl_pred_new", FALSE)) {
    group_ef <- proc_ranef(ordinal::ranef(x),
                           coef_nms = lapply(group_vc_raw, rownames), ...)
    subparams <- c(subparams, group_ef)
  }

  return(subparams)
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

  nms <- mknms_categ(dimnames(coefs), ...)
  subparams <- c(setNames(as.vector(coefs), nms), group_vc)

  if (!getOption("projpred.mlvl_pred_new", FALSE)) {
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
    subparams <- c(subparams, group_ef)
  }

  return(subparams)
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
  if (!x$const_wdraws_prj) {
    warning("The projected draws have different (i.e., nonconstant) weights. ",
            "Therefore, the results from this as.matrix() method should not ",
            "be used without taking these weights into account.")
  }
  if (identical(nm_scheme, "auto")) {
    if (inherits(x$refmodel$fit, "brmsfit")) {
      nm_scheme <- "brms"
    } else {
      nm_scheme <- "rstanarm"
    }
  }
  stopifnot(nm_scheme %in% c("rstanarm", "brms"))
  res <- do.call(rbind, lapply(x$outdmin, get_subparams, nm_scheme = nm_scheme))
  if (x$refmodel$family$family == "gaussian") res <- cbind(res, sigma = x$dis)
  return(res)
}

#' Create cross-validation folds
#'
#' These are helper functions to create cross-validation (CV) folds, i.e., to
#' split up the indices from 1 to `n` into `K` subsets ("folds") for
#' \eqn{K}-fold CV. These functions are potentially useful when creating the
#' `cvfits` and `cvfun` arguments for [init_refmodel()]. Function [cvfolds()] is
#' deprecated; please use [cv_folds()] instead (apart from the name, they are
#' the same). The return value of [cv_folds()] and [cv_ids()] is different, see
#' below for details.
#'
#' @name cv-indices
#'
#' @param n Number of observations.
#' @param K Number of folds. Must be at least 2 and not exceed `n`.
#' @param out Format of the output, either `"foldwise"` or `"indices"`. See
#'   below for details.
#' @param seed Pseudorandom number generation (PRNG) seed by which the same
#'   results can be obtained again if needed. Passed to argument `seed` of
#'   [set.seed()], but can also be `NA` to not call [set.seed()] at all. If not
#'   `NA`, then the PRNG state is reset (to the state before calling
#'   [cv_folds()] or [cv_ids()]) upon exiting [cv_folds()] or [cv_ids()].
#'
#' @return [cv_folds()] returns a vector of length `n` such that each element is
#'   an integer between 1 and `K` denoting which fold the corresponding data
#'   point belongs to. The return value of [cv_ids()] depends on the `out`
#'   argument. If `out = "foldwise"`, the return value is a `list` with `K`
#'   elements, each being a `list` with elements `tr` and `ts` giving the
#'   training and test indices, respectively, for the corresponding fold. If
#'   `out = "indices"`, the return value is a `list` with elements `tr` and `ts`
#'   each being a `list` with `K` elements giving the training and test indices,
#'   respectively, for each fold.
#'
#' @examples
#' n <- 100
#' set.seed(1234)
#' y <- rnorm(n)
#' cv <- cv_ids(n, K = 5)
#' # Mean within the test set of each fold:
#' cvmeans <- sapply(cv, function(fold) mean(y[fold$ts]))
#'
NULL

#' @rdname cv-indices
#' @export
cv_folds <- function(n, K, seed = NA) {
  validate_num_folds(K, n)

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

  ## create and shuffle the indices
  folds <- rep_len(seq_len(K), length.out = n)
  folds <- sample(folds, n, replace = FALSE)

  return(folds)
}

#' @rdname cv-indices
#' @export
cvfolds <- function(n, K, seed = NA) {
  warning("cvfolds() is deprecated. Please use cv_folds() instead.")
  cv_folds(n = n, K = K, seed = seed)
}

#' @rdname cv-indices
#' @export
cv_ids <- function(n, K, out = c("foldwise", "indices"), seed = NA) {
  validate_num_folds(K, n)
  out <- match.arg(out)

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

#' Retrieve the full-data solution path from a [varsel()] or [cv_varsel()] run
#' or the predictor combination from a [project()] run
#'
#' The [solution_terms.vsel()] method retrieves the solution path from a
#' full-data search (`vsel` objects are returned by [varsel()] or
#' [cv_varsel()]). The [solution_terms.projection()] method retrieves the
#' predictor combination onto which a projection was performed (`projection`
#' objects are returned by [project()], possibly as elements of a `list`). Both
#' methods (and hence also the [solution_terms()] generic) are deprecated and
#' will be removed in a future release. Please use [ranking()] instead of
#' [solution_terms.vsel()] ([ranking()]'s output element `fulldata` contains the
#' full-data predictor ranking that is extracted by [solution_terms.vsel()];
#' [ranking()]'s output element `foldwise` contains the fold-wise predictor
#' rankings---if available---which were previously not accessible via a built-in
#' function) and [predictor_terms()] instead of [solution_terms.projection()].
#'
#' @param object The object from which to retrieve the predictor terms. Possible
#'   classes may be inferred from the names of the corresponding methods (see
#'   also the description).
#' @param ... Currently ignored.
#'
#' @return A character vector of predictor terms.
#'
#' @export
solution_terms <- function(object, ...) {
  UseMethod("solution_terms")
}

#' @rdname solution_terms
#' @export
solution_terms.vsel <- function(object, ...) {
  warning("solution_terms.vsel() is deprecated. Please use ranking() instead ",
          "(ranking()'s output element `fulldata` contains the full-data ",
          "predictor ranking that is also extracted by solution_terms.vsel(); ",
          "ranking()'s output element `foldwise` contains fold-wise predictor ",
          "rankings which were previously not accessible via a function).")
  return(ranking(object)[["fulldata"]])
}

#' @rdname solution_terms
#' @export
solution_terms.projection <- function(object, ...) {
  warning("solution_terms.projection() is deprecated. Please use ",
          "predictor_terms() instead.")
  return(predictor_terms(object))
}

#' Predictor terms used in a [project()] run
#'
#' For a `projection` object (returned by [project()], possibly as elements of a
#' `list`), this function extracts the combination of predictor terms onto which
#' the projection was performed.
#'
#' @param object An object of class `projection` (returned by [project()],
#'   possibly as elements of a `list`) from which to retrieve the predictor
#'   terms.
#' @param ... Currently ignored.
#'
#' @return A character vector of predictor terms.
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
#'   print(predictor_terms(prj)) # gives `c("X1", "X3", "X5")`
#' }
#'
#' @export
predictor_terms <- function(object, ...) {
  UseMethod("predictor_terms")
}

#' @rdname predictor_terms
#' @export
predictor_terms.projection <- function(object, ...) {
  return(object[["solution_terms"]])
}

#' Predictor ranking(s)
#'
#' Extracts the *predictor ranking(s)* from an object of class `vsel` (returned
#' by [varsel()] or [cv_varsel()]). A predictor ranking is simply a character
#' vector of predictor terms ranked by predictive relevance (with the most
#' relevant term first). In any case, objects of class `vsel` contain the
#' predictor ranking based on the *full-data* search. If an object of class
#' `vsel` is based on a cross-validation (CV) with fold-wise searches (i.e., if
#' it was created by [cv_varsel()] with `validate_search = TRUE`), then it also
#' contains *fold-wise* predictor rankings.
#'
#' @param object The object from which to retrieve the predictor ranking(s).
#'   Possible classes may be inferred from the names of the corresponding
#'   methods (see also the description).
#' @param nterms_max Maximum submodel size (number of predictor terms) for the
#'   predictor ranking(s), i.e., the submodel size at which to cut off the
#'   predictor ranking(s). Using `NULL` is effectively the same as setting
#'   `nterms_max` to the full model size, i.e., this means to not cut off the
#'   predictor ranking(s) at all. Note that `nterms_max` does not count the
#'   intercept, so `nterms_max = 1` corresponds to the submodel consisting of
#'   the first (non-intercept) predictor term.
#' @param ... Currently ignored.
#'
#' @return An object of class `ranking` which is a `list` with the following
#'   elements:
#'   * `fulldata`: The predictor ranking from the full-data search.
#'   * `foldwise`: The predictor rankings from the fold-wise
#'   searches in the form of a character matrix (only available if `object` is
#'   based on a CV with fold-wise searches, otherwise element `foldwise` is
#'   `NULL`). The rows of this matrix correspond to the CV folds and the columns
#'   to the submodel sizes. Each row contains the predictor ranking from the
#'   search of that CV fold.
#'
#' @seealso [cv_proportions()]
#'
#' @examples
#' # For an example, see `?plot.cv_proportions`.
#'
#' @export
ranking <- function(object, ...) {
  UseMethod("ranking")
}

#' @rdname ranking
#' @export
ranking.vsel <- function(object, nterms_max = NULL, ...) {
  if (is.null(object$projpred_version) && !is.null(object$cv_method)) {
    warning(
      "It seems like a projpred version <= 2.5.0 was used for creating the ",
      "`vsel` object. Thus, even if there are fold-wise searches, the ",
      "corresponding fold-wise predictor rankings cannot be extracted."
    )
  }
  out <- list(fulldata = object[["solution_terms"]],
              foldwise = object[["solution_terms_cv"]])
  if (!is.null(nterms_max)) {
    out[["fulldata"]] <- utils::head(out[["fulldata"]], nterms_max)
    if (!is.null(out[["foldwise"]])) {
      out[["foldwise"]] <- out[["foldwise"]][, seq_len(nterms_max),
                                             drop = FALSE]
    }
  }
  if (!is.null(out[["foldwise"]]) &&
      length(out[["fulldata"]]) != ncol(out[["foldwise"]])) {
    stop("Unexpected dimensions of ranking() output. Please notify the ",
         "package maintainer.")
  }
  class(out) <- "ranking"
  return(out)
}

#' Ranking proportions from fold-wise predictor rankings
#'
#' Calculates the *ranking proportions* from the fold-wise predictor rankings in
#' a cross-validation (CV) with fold-wise searches. For a given predictor
#' \eqn{x} and a given submodel size \eqn{j}, the ranking proportion is the
#' proportion of CV folds which have predictor \eqn{x} at position \eqn{j} of
#' their predictor ranking. While these ranking proportions are helpful for
#' investigating variability in the predictor ranking, they can also be
#' *cumulated* across submodel sizes. The cumulated ranking proportions are more
#' helpful when it comes to model selection.
#'
#' @param object For [cv_proportions.ranking()]: an object of class `ranking`
#'   (returned by [ranking()]). For [cv_proportions.vsel()]: an object of class
#'   `vsel` (returned by [varsel()] or [cv_varsel()]) that [ranking()] will be
#'   applied to internally before then calling [cv_proportions.ranking()].
#' @param cumulate A single logical value indicating whether the ranking
#'   proportions should be cumulated across increasing submodel sizes (`TRUE`)
#'   or not (`FALSE`).
#' @param ... For [cv_proportions.vsel()]: arguments passed to [ranking.vsel()]
#'   and [cv_proportions.ranking()]. For [cv_proportions.ranking()]: currently
#'   ignored.
#'
#' @return A numeric matrix containing the ranking proportions. This matrix has
#'   `nterms_max` rows and `nterms_max` columns, with `nterms_max` as specified
#'   in the (possibly implicit) [ranking()] call. The rows correspond to the
#'   submodel sizes and the columns to the predictor terms (sorted according to
#'   the full-data predictor ranking). If `cumulate` is `FALSE`, then the
#'   returned matrix is of class `cv_proportions`. If `cumulate` is `TRUE`, then
#'   the returned matrix is of classes `cv_proportions_cumul` and
#'   `cv_proportions` (in this order).
#'
#'   Note that if `cumulate` is `FALSE`, then the values in the returned matrix
#'   only need to sum to 1 (column-wise and row-wise) if `nterms_max` (see
#'   above) is equal to the full model size. Likewise, if `cumulate` is `TRUE`,
#'   then the value `1` only needs to occur in each column of the returned
#'   matrix if `nterms_max` is equal to the full model size.
#'
#'   The [cv_proportions()] function is only applicable if the `ranking` object
#'   includes fold-wise predictor rankings (i.e., if it is based on a `vsel`
#'   object created by [cv_varsel()] with `validate_search = TRUE`). If the
#'   `ranking` object contains only a full-data predictor ranking (i.e., if it
#'   is based on a `vsel` object created by [varsel()] or by [cv_varsel()], but
#'   the latter with `validate_search = FALSE`), then an error is thrown because
#'   in that case, there are no fold-wise predictor rankings from which to
#'   calculate ranking proportions.
#'
#' @seealso [plot.cv_proportions()]
#'
#' @examples
#' # For an example, see `?plot.cv_proportions`.
#'
#' @export
cv_proportions <- function(object, ...) {
  UseMethod("cv_proportions")
}

#' @rdname cv_proportions
#' @export
cv_proportions.ranking <- function(object, cumulate = FALSE, ...) {
  cv_paths <- object[["foldwise"]]
  if (is.null(cv_paths)) {
    stop("Could not find fold-wise predictor rankings from which to calculate ",
         "ranking proportions. The reason is probably that `object` is not ",
         "based on a cross-validation or that the search has been excluded ",
         "from the cross-validation.")
  }
  if (ncol(cv_paths) == 0) {
    stop("Needing `nterms_max >= 1` in the (possibly implicit) ranking() call.")
  }
  # Calculate the ranking proportions. Note that the following code assumes that
  # all CV folds have equal weight.
  cv_props <- do.call(cbind, lapply(
    setNames(nm = object[["fulldata"]]),
    function(predictor_j) {
      # We need `na.rm = TRUE` for subsampled LOO CV:
      colMeans(cv_paths == predictor_j, na.rm = TRUE)
    }
  ))
  rownames(cv_props) <- seq_len(nrow(cv_props))
  classes_out <- "cv_proportions"
  if (cumulate) {
    cv_props <- do.call(cbind, apply(cv_props, 2, cumsum, simplify = FALSE))
    rownames(cv_props) <- paste0("<=", rownames(cv_props))
    classes_out <- c("cv_proportions_cumul", classes_out)
  }
  # Setting the `dimnames` names here (not before the `if (cumulate)` part)
  # because `simplify = FALSE` in apply() makes it impossible to keep these:
  names(dimnames(cv_props)) <- c("size", "predictor")
  class(cv_props) <- classes_out
  return(cv_props)
}

#' @rdname cv_proportions
#' @export
cv_proportions.vsel <- function(object, ...) {
  cv_proportions(ranking(object, ...), ...)
}

#' Plot ranking proportions from fold-wise predictor rankings
#'
#' Plots the ranking proportions (see [cv_proportions()]) from the fold-wise
#' predictor rankings in a cross-validation with fold-wise searches. This is a
#' visualization of the *transposed* matrix returned by [cv_proportions()]. The
#' proportions printed as text inside of the colored tiles are rounded to whole
#' percentage points (the plotted proportions themselves are not rounded).
#'
#' @param x For [plot.cv_proportions()]: an object of class `cv_proportions`
#'   (returned by [cv_proportions()], possibly with `cumulate = TRUE`). For
#'   [plot.ranking()]: an object of class `ranking` (returned by [ranking()])
#'   that [cv_proportions()] will be applied to internally before then calling
#'   [plot.cv_proportions()].
#' @param text_angle Passed to argument `angle` of [ggplot2::element_text()] for
#'   the y-axis tick labels. In case of long predictor names, `text_angle = 45`
#'   might be helpful (for example).
#' @param ... For [plot.ranking()]: arguments passed to
#'   [cv_proportions.ranking()] and [plot.cv_proportions()]. For
#'   [plot.cv_proportions()]: currently ignored.
#'
#' @return A \pkg{ggplot2} plotting object (of class `gg` and `ggplot`).
#'
#' @author Idea and original code by Aki Vehtari. Slight modifications of the
#'   original code by Frank Weber, Yann McLatchie, and Sölvi Rögnvaldsson. Final
#'   implementation in \pkg{projpred} by Frank Weber.
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
#'
#'   # Extract predictor rankings:
#'   rk <- ranking(cvvs)
#'
#'   # Compute ranking proportions:
#'   pr_rk <- cv_proportions(rk)
#'
#'   # Visualize the ranking proportions:
#'   gg_pr_rk <- plot(pr_rk)
#'   print(gg_pr_rk)
#'
#'   # Since the object returned by plot.cv_proportions() is a standard ggplot2
#'   # plotting object, you can modify the plot easily, e.g., to remove the
#'   # legend:
#'   print(gg_pr_rk + ggplot2::theme(legend.position = "none"))
#' }
#'
#' @export
plot.cv_proportions <- function(x, text_angle = NULL, ...) {
  cv_props_long <- data.frame(
    msize = factor(rep(rownames(x), times = ncol(x)), levels = rownames(x)),
    pterm = factor(rep(colnames(x), each = nrow(x)), levels = colnames(x)),
    propcv = as.vector(x)
  )
  cv_props_long$txtcolor <- ifelse(cv_props_long$propcv > 0.5, "white", "black")
  gg_cv_props <- ggplot(data = cv_props_long,
                        mapping = aes(x = .data[["msize"]],
                                      y = .data[["pterm"]])) +
    geom_tile(mapping = aes(fill = .data[["propcv"]]),
              width = 1, height = 1, linewidth = 1, color = "white") +
    # Note: The original code for this function specified argument `fontface`
    # in the aes() call of geom_text(), but incorrectly (as constantly `1`):
    geom_text(mapping = aes(label = paste(round(100 * .data[["propcv"]]), "%"),
                            color = I(.data[["txtcolor"]])),
              size = 3) +
    scale_y_discrete(limits = rev(levels(cv_props_long$pterm))) +
    # Filling color:
    ### Option 1:
    scale_fill_gradient(name = "Proportion\nof CV folds",
                        labels = scales::label_percent(suffix = " %"),
                        limits = c(0, 1),
                        low = "#ededed", high = "#0f365c") +
    ###
    ### Option 2 (requires the 'RColorBrewer' package):
    # scale_fill_distiller(name = "Proportion\nof CV folds",
    #                      labels = scales::label_percent(suffix = " %"),
    #                      direction = 1) +
    ###
    labs(x = "Submodel size (number of predictor terms)", y = "Predictor") +
    coord_cartesian(expand = FALSE) +
    theme(axis.text.y = element_text(angle = text_angle))
  return(gg_cv_props)
}

#' @rdname plot.cv_proportions
#' @export
plot.ranking <- function(x, ...) {
  plot(cv_proportions(x, ...), ...)
}
