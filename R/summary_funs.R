# Calculate log posterior(-projection) predictive density values and average
# them across parameter draws (together with the corresponding expected response
# values).
#
# @param y_wobs_test A `list` (but we encourage to use a `data.frame`), at least
#   with elements (columns) `y` (response values) and `wobs` (observation
#   weights). In case of the latent projection, this `list` (or `data.frame`)
#   also needs to contain `y_oscale` (response values on the original response
#   scale, i.e., the non-latent response values).
# @param family A `family` object.
# @param wdraws A vector of weights for the parameter draws.
# @param mu A matrix of expected values for `y`.
# @param dis A vector of dispersion parameter draws.
# @param cl_ref A numeric vector of length \eqn{S} (with \eqn{S} denoting the
#   number of parameter draws in the reference model), giving the cluster
#   indices for the parameter draws in the reference model. Draws that should be
#   dropped (e.g., because of thinning by `ndraws` or `ndraws_pred`) need to
#   have an `NA` in `cl_ref`. Caution: This always refers to the reference
#   model's parameter draws, not necessarily to the columns of `mu`, the entries
#   of `wdraws`, or the entries of `dis`!
# @param wdraws_ref A numeric vector of length \eqn{S} (with \eqn{S} denoting
#   the number of parameter draws in the reference model), giving the weights of
#   the parameter draws in the reference model. It doesn't matter whether these
#   are normalized (i.e., sum to `1`) or not because the family$latent_ilink()
#   function that receives them should treat them as unnormalized. Draws that
#   should be dropped (e.g., because of thinning by `ndraws` or `ndraws_pred`)
#   can (but must not necessarily) have an `NA` in `wdraws_ref`.
#
# @return A `list` with elements `mu` and `lppd` which are both vectors
#   containing the values for the quantities from the description above.
weighted_summary_means <- function(y_wobs_test, family, wdraws, mu, dis, cl_ref,
                                   wdraws_ref = rep(1, length(cl_ref))) {
  if (!is.matrix(mu) || any(dim(mu) == 0)) {
    stop("Unexpected structure for `mu`. Do the return values of ",
         "`proj_predfun` and `ref_predfun` have the correct structure?")
  }
  loglik <- family$ll_fun(mu, dis, y_wobs_test$y, y_wobs_test$wobs)
  if (!is.matrix(loglik) || any(dim(loglik) == 0)) {
    stop("Unexpected structure for `loglik`. Please notify the package ",
         "maintainer.")
  }
  # Average over the draws, taking their weights into account:
  avg <- list(
    mu = structure(c(mu %*% wdraws),
                   ndiscrete = attr(mu, "ndiscrete"),
                   class = sub("augmat", "augvec", oldClass(mu), fixed = TRUE)),
    lppd = apply(loglik, 1, log_weighted_mean_exp, wdraws)
  )
  if (family$for_latent) {
    mu_oscale <- family$latent_ilink(t(mu), cl_ref = cl_ref,
                                     wdraws_ref = wdraws_ref)
    if (length(dim(mu_oscale)) < 2) {
      stop("Unexpected structure for the output of `latent_ilink`.")
    }
    loglik_oscale <- family$latent_ll_oscale(
      mu_oscale, y_oscale = y_wobs_test$y_oscale, wobs = y_wobs_test$wobs,
      cl_ref = cl_ref, wdraws_ref = wdraws_ref
    )
    if (!is.matrix(loglik_oscale)) {
      stop("Unexpected structure for the output of `latent_ll_oscale`.")
    }
    if (length(dim(mu_oscale)) == 3) {
      # In this case, `mu_oscale` is a 3-dimensional array (S x N x C), so
      # coerce it to an augmented-rows matrix:
      mu_oscale <- arr2augmat(mu_oscale, margin_draws = 1)
      mu_oscale_avg <- structure(
        c(mu_oscale %*% wdraws),
        ndiscrete = attr(mu_oscale, "ndiscrete"),
        class = sub("augmat", "augvec", oldClass(mu_oscale), fixed = TRUE)
      )
    } else {
      # In principle, we could use the same code for `mu_oscale_avg` as above.
      # However, that would require `mu_oscale <- t(mu_oscale)` beforehand, so
      # the following should be more efficient:
      mu_oscale_avg <- c(wdraws %*% mu_oscale)
    }
    avg$oscale <- list(
      mu = mu_oscale_avg,
      lppd = apply(loglik_oscale, 2, log_weighted_mean_exp, wdraws)
    )
  }
  return(avg)
}

# A function to calculate the desired performance statistics, their standard
# errors, and confidence intervals with coverage `1 - alpha` based on the
# variable selection output. If `nfeat_baseline` is given, then compute the
# statistics relative to the baseline model of that size (`nfeat_baseline = Inf`
# means that the baseline model is the reference model).
.tabulate_stats <- function(varsel, stats, alpha = 0.05,
                            nfeat_baseline = NULL, resp_oscale = TRUE, ...) {
  stat_tab <- data.frame()
  summaries_ref <- varsel$summaries$ref
  summaries_sub <- varsel$summaries$sub
  summaries_fast_sub <- varsel$summaries_fast$sub
  if (!is.null(summaries_fast_sub) && any(stats %in% c("auc"))) {
    stop("Subsampled LOO-CV with AUC not implemented. Alternatives using ",
         "`validate_search = TRUE` are full (i.e., non-subsampled) LOO-CV and ",
         "K-fold CV. Otherwise, results from `validate_search = FALSE` (which ",
         "often already exist at this point of the workflow) can be used, ",
         "with the downside that the search part is not cross-validated in ",
         "that case.")
  }

  if (!varsel$refmodel$family$for_latent && !resp_oscale) {
    stop("`resp_oscale = FALSE` can only be used in case of the latent ",
         "projection.")
  }
  if (varsel$refmodel$family$for_latent) {
    if (resp_oscale) {
      summaries_ref <- summaries_ref$oscale
      summaries_sub <- lapply(summaries_sub, "[[", "oscale")
      if (!is.null(summaries_fast_sub)) {
        summaries_fast_sub <- lapply(summaries_fast_sub, "[[", "oscale")
      }
      ref_lppd_NA <- all(is.na(summaries_ref$lppd))
      sub_lppd_NA <- any(sapply(summaries_sub, check_sub_NA, el_nm = "lppd"))
      if (!is.null(summaries_fast_sub)) {
        fast_sub_lppd_NA <- any(sapply(summaries_fast_sub, check_sub_NA, el_nm = "lppd"))
      } else {
        fast_sub_lppd_NA <- FALSE
      }
      ref_mu_NA <- all(is.na(summaries_ref$mu))
      sub_mu_NA <- any(sapply(summaries_sub, check_sub_NA, el_nm = "mu"))
      if (!is.null(summaries_fast_sub)) {
        fast_sub_mu_NA <- any(sapply(summaries_fast_sub, check_sub_NA, el_nm = "mu"))
      } else {
        fast_sub_mu_NA <- FALSE
      }
      if (all(is.na(varsel$y_wobs_test$y_oscale))) {
        message(
          "Cannot calculate performance statistics if `resp_oscale = TRUE` ",
          "and `<vsel>$y_wobs_test$y_oscale` consists of only `NA`s."
        )
      } else if (ref_mu_NA || sub_mu_NA || fast_sub_mu_NA) {
        message(
          "`latent_ilink` returned only `NA`s, so all performance statistics ",
          "will also be `NA` as long as `resp_oscale = TRUE`."
        )
      } else if (any(stats %in% c("elpd", "mlpd", "gmpd")) &&
                 (ref_lppd_NA || sub_lppd_NA || fast_sub_lppd_NA)) {
        message(
          "`latent_ll_oscale` returned only `NA`s, so ELPD, MLPD, and GMPD ",
          "will also be `NA` as long as `resp_oscale = TRUE`."
        )
      }
      varsel$y_wobs_test$y <- varsel$y_wobs_test$y_oscale
    } else {
      if (all(is.na(varsel$refmodel$dis)) &&
          any(stats %in% c("elpd", "mlpd", "gmpd"))) {
        message(
          "Cannot calculate ELPD, MLPD, or GMPD if `resp_oscale = FALSE` and ",
          "`<refmodel>$dis` consists of only `NA`s. If it is not possible to ",
          "supply values to argument `dis` of init_refmodel(), consider (i) ",
          "switching to `resp_oscale = TRUE` (which might require the ",
          "specification of functions needed by extend_family()) or (ii) ",
          "using a performance statistic other than ELPD, MLPD, or GMPD."
        )
      }
      if (all(is.na(varsel$y_wobs_test$y))) {
        mssg_y_NA <- paste0(
          "Cannot calculate performance statistics if `resp_oscale = FALSE` ",
          "and `<vsel>$y_wobs_test$y` consists of only `NA`s."
        )
        if (identical(varsel$cv_method, "kfold")) {
          mssg_y_NA <- paste0(
            mssg_y_NA, " The reason for these `NA`s is probably that `<vsel>` ",
            "was created by cv_varsel() with `cv_method = \"kfold\"`. (In ",
            "case of K-fold cross-validation, the latent response values for ",
            "the test datasets cannot be defined in a straightforward manner ",
            "without inducing dependencies between training and test datasets.)"
          )
        }
        message(mssg_y_NA)
      }
    }
  }
  # Just to avoid that `$y` gets expanded to `$y_oscale` if element `y` does not
  # exist (for whatever reason; actually, it should always exist):
  varsel$y_wobs_test$y_oscale <- NULL

  if (resp_oscale && !is.null(varsel$refmodel$family$cats) &&
      any(stats %in% c("acc", "pctcorr"))) {
    summaries_ref$mu <- catmaxprb(summaries_ref$mu, lvls = varsel$refmodel$family$cats)
    summaries_sub <- lapply(summaries_sub, function(summaries_sub_k) {
      summaries_sub_k$mu <- catmaxprb(summaries_sub_k$mu,
                                      lvls = varsel$refmodel$family$cats)
      return(summaries_sub_k)
    })
    if (!is.null(summaries_fast_sub)) {
      summaries_fast_sub <- lapply(summaries_fast_sub, function(summaries_fast_sub_k) {
        summaries_fast_sub_k$mu <- catmaxprb(summaries_fast_sub_k$mu,
                                             lvls = varsel$refmodel$family$cats)
        return(summaries_fast_sub_k)
      })
    }
    # Since `mu` is an unordered factor, `y` needs to be unordered, too (or both
    # would need to be ordered; however, unordered is the simpler type):
    varsel$y_wobs_test$y <- factor(varsel$y_wobs_test$y, ordered = FALSE)
  }

  if (varsel$refmodel$family$family == "binomial" &&
      !all(varsel$y_wobs_test$wobs == 1)) {
    # This case should not occur (yet) for the augmented-data or the latent
    # projection:
    stopifnot(!varsel$refmodel$family$for_augdat)
    stopifnot(!varsel$refmodel$family$for_latent)
    varsel$y_wobs_test$y_prop <- varsel$y_wobs_test$y / varsel$y_wobs_test$wobs
  }

  ## fetch the mu and lppd for the baseline model
  if (is.null(nfeat_baseline)) {
    ## no baseline model, i.e, compute the statistics on the actual
    ## (non-relative) scale
    summaries_baseline <- NULL
    delta <- FALSE
  } else {
    if (nfeat_baseline == Inf) {
      summaries_baseline <- summaries_ref
    } else {
      summaries_baseline <- summaries_sub[[nfeat_baseline + 1]]
    }
    delta <- TRUE
  }

  for (s in seq_along(stats)) {
    stat <- stats[s]

    ## reference model statistics
    res <- get_stat(summaries = summaries_ref,
                    summaries_baseline = summaries_baseline,
                    summaries_fast = NULL,
                    loo_inds = NULL,
                    varsel$y_wobs_test, stat, alpha = alpha, ...)
    row <- data.frame(
      data = varsel$type_test, size = Inf, delta = delta, statistic = stat,
      value = res$value, lq = res$lq, uq = res$uq, se = res$se, diff = NA,
      diff.se = NA
    )
    stat_tab <- rbind(stat_tab, row)

    ## submodel statistics
    for (k in seq_along(summaries_sub)) {
      res <- get_stat(summaries = summaries_sub[[k]],
                      summaries_baseline = summaries_baseline,
                      summaries_fast = summaries_fast_sub[[k]],
                      loo_inds = varsel$loo_inds,
                      varsel$y_wobs_test, stat, alpha = alpha, ...)
      diff <- get_stat(summaries = summaries_sub[[k]],
                       summaries_baseline = summaries_ref,
                       summaries_fast = summaries_fast_sub[[k]],
                       loo_inds = varsel$loo_inds,
                       varsel$y_wobs_test, stat, alpha = alpha, ...)
      row <- data.frame(
        data = varsel$type_test, size = k - 1, delta = delta, statistic = stat,
        value = res$value, lq = res$lq, uq = res$uq, se = res$se,
        diff = diff$value, diff.se = diff$se
      )
      stat_tab <- rbind(stat_tab, row)
    }
  }

  return(stat_tab)
}

# Helper function checking whether all entries of a summaries vector are `NA`.
#
# @param summaries_sub_k Typically `<vsel_object>$summaries$sub[[k]]`.
# @param el_nm A single character string, giving the name of the subelement of
#   `summaries_sub_k` to check for `NA`s.
#
# @return A single logical value, indicating whether all entries of
#   `summaries_sub_k[[el_nm]]` are `NA`.
check_sub_NA <- function(summaries_sub_k, el_nm) {
  all(is.na(summaries_sub_k[[el_nm]]))
}

## Calculates given statistic stat with standard error and confidence bounds.
## `summaries_baseline` contains the pointwise mu and lppd for another model
## that is used as a baseline for computing the difference (ratio in case of the
## GMPD) in the given statistic. If these arguments are not given (NULL) then
## the actual (non-relative) value is computed. The actual observation weights
## (specified by the user) are contained in `y_wobs_test$wobs`. These are
## already taken into account by `<refmodel_object>$family$ll_fun()` (or
## `<refmodel_object>$family$latent_ll_oscale()`) and are thus already taken
## into account in `lppd`. However, `mu` does not take them into account, so
## some further adjustments are necessary below.
get_stat <- function(summaries, summaries_baseline = NULL,
                     summaries_fast = NULL, loo_inds = NULL,
                     y_wobs_test, stat, alpha = 0.1, ...) {
  mu <- summaries$mu
  lppd <- summaries$lppd
  n_full <- length(lppd)
  n_loo <- if (is.null(loo_inds)) n_full else length(loo_inds)
  alpha_half <- alpha / 2
  one_minus_alpha_half <- 1 - alpha_half

  if (stat %in% c("elpd", "mlpd", "gmpd")) {
    if (is.null(summaries_baseline)) {
      lppd_baseline <- 0
    } else {
      lppd_baseline <- summaries_baseline$lppd
    }
    if (n_loo < n_full) {
      # subsampling difference estimator (Magnusson et al., 2020)
      srs_diffe <- .srs_diff_est_w(y_approx = summaries_fast$lppd - lppd_baseline,
                                   y = (lppd - lppd_baseline)[loo_inds],
                                   y_idx = loo_inds)
      value <- srs_diffe$y_hat
      # combine estimates of var(y_hat) and var(y)
      value_se <- sqrt(srs_diffe$v_y_hat + srs_diffe$hat_v_y)
    } else {
      # full LOO estimator
      value <- sum(lppd - lppd_baseline)
      value_se <- sd(lppd - lppd_baseline) * sqrt(n_full)
    }
    if (stat %in% c("mlpd", "gmpd")) {
      value <- value / n_full
      value_se <- value_se / n_full
      if (stat == "gmpd") {
        value_gmpd <- exp(value)
        # delta method
        value_gmpd_se <- value_se * value_gmpd
      }
    }
  } else if (stat %in% c("mse", "rmse", "R2")) {
    y <- y_wobs_test$y_prop %||% y_wobs_test$y
    wobs <- y_wobs_test$wobs
    wobs <- n_full * wobs / sum(wobs)
    if (!is.null(summaries_baseline)) {
      mu_baseline <- summaries_baseline$mu
    }
    # Use exact standard error for mse and delta method for rmse and R2
    if (n_loo < n_full) {
      # subsampling difference estimator (Magnusson et al., 2020)
      srs_diffe <- .srs_diff_est_w(y_approx = (summaries_fast$mu - y)^2,
                                   y = ((mu - y)^2)[loo_inds],
                                   y_idx = loo_inds,
                                   wobs = wobs)
      value <- srs_diffe$y_hat / n_full
      # combine estimates of var(y_hat) and var(y)
      value_se <- sqrt(srs_diffe$v_y_hat + srs_diffe$hat_v_y) / n_full
    } else {
      # full LOO estimator
      value <- mean(wobs * (mu - y)^2)
      value_se <- .weighted_sd((mu - y)^2, wobs) / sqrt(n_full)
    }
    # store for later calculations
    mse_e <- value
    if (!is.null(summaries_baseline)) {
      # delta=TRUE, variance of difference of two normally distributed
      # quantities (log-normally in case of MSE and RMSE, although the central
      # limit theorem would ensure convergence -- probably slower, though -- to
      # a normal distribution even for MSE and RMSE)
      mse_b <- mean(wobs * (mu_baseline - y)^2)
      var_mse_b <- .weighted_sd((mu_baseline - y)^2, wobs)^2 / n_full
      if (n_loo < n_full) {
        mse_e_fast <- mean(wobs * (summaries_fast$mu - y)^2)
        srs_diffe <-
          .srs_diff_est_w(y_approx = ((summaries_fast$mu - y)^2 - mse_e_fast) *
                            ((mu_baseline - y)^2 - mse_b),
                          y = (((mu - y)^2 - mse_e) *
                                 ((mu_baseline - y)^2 - mse_b))[loo_inds],
                          y_idx = loo_inds,
                          wobs = wobs)
        cov_mse_e_b <- srs_diffe$y_hat / (n_full * (n_full - 1))
      } else {
        cov_mse_e_b <- mean(wobs * ((mu - y)^2 - mse_e) *
                              ((mu_baseline - y)^2 - mse_b)) / (n_full - 1)
      }
      if (stat != "rmse") {
        value_se <- sqrt(value_se^2 - 2 * cov_mse_e_b + var_mse_b)
      }
    }
    if (stat == "mse") {
      value <- mse_e - ifelse(is.null(summaries_baseline), 0, mse_b)
    } else if (stat == "rmse") {
      # simple transformation of mse
      value <- sqrt(mse_e) - ifelse(is.null(summaries_baseline), 0, sqrt(mse_b))
      # the first-order Taylor approximation of the variance
      if (is.null(summaries_baseline)) {
        value_se <- sqrt(value_se^2 / mse_e / 4)
      } else {
        value_se <- sqrt((value_se^2 / mse_e -
                            2 * cov_mse_e_b / sqrt(mse_e * mse_b) +
                            var_mse_b / mse_b) / 4)
      }
    } else if (stat == "R2") {
      y_mean_w <- mean(wobs * y)
      # simple transformation of mse
      mse_y <- mean(wobs * (y_mean_w - y)^2)
      value <- 1 - mse_e / mse_y - ifelse(is.null(summaries_baseline), 0, 1 - mse_b / mse_y)
      # the first-order Taylor approximation of the variance
      var_mse_y <- .weighted_sd((y_mean_w - y)^2, wobs)^2 / n_full
      if (n_loo < n_full) {
        mse_e_fast <- mean(wobs * (summaries_fast$mu - y)^2)
        if (is.null(summaries_baseline)) {
          srs_diffe <-
            .srs_diff_est_w(y_approx = ((summaries_fast$mu - y)^2 - mse_e_fast) *
                              ((y_mean_w - y)^2 - mse_y),
                            y = (((mu - y)^2 - mse_e) *
                                   ((y_mean_w - y)^2 - mse_y))[loo_inds],
                            y_idx = loo_inds,
                            wobs = wobs)
        } else {
          srs_diffe <-
            .srs_diff_est_w(y_approx = ((summaries_fast$mu - y)^2 - mse_e_fast -
                                          ((mu_baseline - y)^2 - mse_b)) *
                              ((y_mean_w - y)^2 - mse_y),
                            y = (((mu - y)^2 - mse_e -
                                    ((mu_baseline - y)^2 - mse_b)) *
                                   ((y_mean_w - y)^2 - mse_y))[loo_inds],
                            y_idx = loo_inds,
                            wobs = wobs)
        }
        cov_mse_e_y <- srs_diffe$y_hat / (n_full * (n_full - 1))
      } else {
        if (is.null(summaries_baseline)) {
          cov_mse_e_y <- mean(wobs * ((mu - y)^2 - mse_e) *
                                ((y_mean_w - y)^2 - mse_y)) / (n_full - 1)
        } else {
          cov_mse_e_y <- mean(wobs * ((mu - y)^2 - mse_e -
                                        ((mu_baseline - y)^2 - mse_b)) *
                                ((y_mean_w - y)^2 - mse_y)) / (n_full - 1)
        }
      }
      # part of delta se comes automatically via mse
      if (!is.null(summaries_baseline)) {
        # delta=TRUE
        mse_e <- mse_e - mse_b
      }
      value_se_sq <- (value_se^2 -
                        2 * mse_e / mse_y * cov_mse_e_y +
                        (mse_e / mse_y)^2 * var_mse_y) / mse_y^2
      if (!is.na(value_se_sq) && sign(value_se_sq) == -1) {
        if (abs(value_se_sq) < sqrt(.Machine$double.eps)) {
          value_se_sq <- 0
        } else {
          stop("Negative (and numerically non-zero) `value_se_sq`.")
        }
      }
      value_se <- sqrt(value_se_sq)
    }
  } else if (stat %in% c("acc", "pctcorr", "auc")) {
    y <- y_wobs_test$y
    # In this case (`stat %in% c("acc", "pctcorr", "auc")`), we hard-code `wobs`
    # to be full of ones because currently, the user-supplied observation
    # weights are required to be (positive) whole numbers, so these observation
    # weights are incorporated by "de-aggregating" the aggregated dataset that
    # was supplied by the user (the term "de-aggregation" refers to the
    # de-aggregation of the multiple Bernoulli trials belonging to one row in
    # the aggregated dataset). Currently, `wobs` is not really useful and could
    # be removed, but we leave it here for the future (perhaps one day, we will
    # not require the user-supplied observation weights to be whole numbers
    # anymore).
    wobs <- rep(1, n_full)
    if (!is.null(y_wobs_test$y_prop)) {
      # CAUTION: The following checks also ensure that `y` does not have `NA`s
      # (see the other "CAUTION" comments below for changes that are needed if
      # `y` is allowed to have `NA`s here):
      stopifnot(all(is_wholenumber(y_wobs_test$wobs)))
      stopifnot(all(is_wholenumber(y)))
      stopifnot(all(0 <= y & y <= y_wobs_test$wobs))
      y <- unlist(lapply(seq_along(y), function(i_short) {
        c(rep(0L, y_wobs_test$wobs[i_short] - y[i_short]),
          rep(1L, y[i_short]))
      }))
      mu <- rep(mu, y_wobs_test$wobs)
      mu_baseline <- rep(summaries_baseline$mu, y_wobs_test$wobs)
      mu_fast <- rep(summaries_fast$mu, y_wobs_test$wobs)
      # CAUTION: If `y` is allowed to have `NA`s here, then the following
      # definition of `n_full` needs to be adapted:
      n_full <- length(mu)
      if (!is.null(loo_inds)) {
        stopifnot(all(y_wobs_test$wobs > 0))
        loo_inds <- unlist(lapply(loo_inds, function(loo_idx) {
          cumsum_wobs <- cumsum(y_wobs_test$wobs)
          if (loo_idx == 1) {
            lower_idx_new <- 1L
          } else {
            lower_idx_new <- cumsum_wobs[loo_idx - 1L] + 1L
          }
          upper_idx_new <- cumsum_wobs[loo_idx]
          return(lower_idx_new:upper_idx_new)
        }))
      }
      n_loo <- if (is.null(loo_inds)) n_full else length(loo_inds)
      wobs <- rep(wobs, y_wobs_test$wobs)
      wobs <- n_full * wobs / sum(wobs)
    } else {
      stopifnot(all(y_wobs_test$wobs == 1))
      mu_baseline <- summaries_baseline$mu
      mu_fast <- summaries_fast$mu
    }
    if (stat %in% c("acc", "pctcorr")) {
      # Find out whether each observation was classified correctly or not:
      if (!is.factor(mu)) {
        mu <- round(mu)
      }
      correct <- mu == y

      if (!is.null(mu_baseline)) {
        if (!is.factor(mu_baseline)) {
          mu_baseline <- round(mu_baseline)
        }
        correct_baseline <- mu_baseline == y
      } else {
        correct_baseline <- 0
      }

      if (n_loo < n_full) {
        # subsampling difference estimator (Magnusson et al., 2020)
        if (!is.factor(mu_fast)) {
          mu_fast <- round(mu_fast)
        }
        correct_fast <- mu_fast == y
        srs_diffe <- .srs_diff_est_w(y_approx = correct_fast - correct_baseline,
                                     y = (correct - correct_baseline)[loo_inds],
                                     y_idx = loo_inds,
                                     wobs = wobs)
        value <- srs_diffe$y_hat / n_full
        # combine estimates of var(y_hat) and var(y)
        value_se <- sqrt(srs_diffe$v_y_hat + srs_diffe$hat_v_y) / n_full
      } else {
        # full LOO estimator
        value <- mean(wobs * correct) - mean(wobs * correct_baseline)
        value_se <- .weighted_sd(correct - correct_baseline, wobs) / sqrt(n_full)
      }
    } else if (stat == "auc") {
      if (n_loo < n_full) {
        # Note: Previously, subsampled LOO with AUC caused the fast LOO results
        # to be used automatically (via `mu <- mu_fast`), see PR #496.
        stop("Subsampled LOO-CV with AUC not implemented.")
      }
      if (!is.null(mu_baseline)) {
        auc_data <- cbind(y, mu, wobs)
        auc_data_baseline <- cbind(y, mu_baseline, wobs)
        value <- .auc(auc_data) - .auc(auc_data_baseline)
        idxs_cols <- seq_len(ncol(auc_data))
        idxs_cols_bs <- setdiff(seq_len(ncol(auc_data) + ncol(auc_data_baseline)),
                                idxs_cols)
        diffvalue.bootstrap <- bootstrap(
          cbind(auc_data, auc_data_baseline),
          function(x) {
            .auc(x[, idxs_cols, drop = FALSE]) -
              .auc(x[, idxs_cols_bs, drop = FALSE])
          },
          ...
        )
        value_se <- sd(diffvalue.bootstrap)
        if (any(is.na(diffvalue.bootstrap))) {
          # quantile() is not able to deal with `NA`s
          lq_uq <- rep(NA_real_, 2)
        } else {
          lq_uq <- quantile(diffvalue.bootstrap,
                            probs = c(alpha_half, one_minus_alpha_half),
                            names = FALSE)
        }
      } else {
        auc_data <- cbind(y, mu, wobs)
        value <- .auc(auc_data)
        value.bootstrap <- bootstrap(auc_data, .auc, ...)
        value_se <- sd(value.bootstrap)
        if (any(is.na(value.bootstrap))) {
          # quantile() is not able to deal with `NA`s
          lq_uq <- rep(NA_real_, 2)
        } else {
          lq_uq <- quantile(value.bootstrap,
                            probs = c(alpha_half, one_minus_alpha_half),
                            names = FALSE)
        }
      }
    }
  }

  if (stat %in% c("mse", "rmse") && is.null(summaries_baseline)) {
    # Compute mean and variance in log scale by matching the variance of a
    # log-normal approximation
    # https://en.wikipedia.org/wiki/Log-normal_distribution#Arithmetic_moments
    mul <- log(value^2 / sqrt(value_se^2 + value^2))
    varl <- log1p(value_se^2 / value^2)
    lq <- qnorm(alpha_half, mean = mul, sd = sqrt(varl))
    uq <- qnorm(one_minus_alpha_half, mean = mul, sd = sqrt(varl))
    # Go back to linear scale
    lq <- exp(lq)
    uq <- exp(uq)
  } else if (stat %in% c("auc")) {
    lq <- lq_uq[1]
    uq <- lq_uq[2]
  } else {
    lq <- qnorm(alpha_half, mean = value, sd = value_se)
    uq <- qnorm(one_minus_alpha_half, mean = value, sd = value_se)
  }

  if (stat == "gmpd") {
    lq <- exp(lq)
    uq <- exp(uq)
    value <- value_gmpd
    value_se <- value_gmpd_se
  }

  return(list(value = value, se = value_se, lq = lq, uq = uq))
}

is_util <- function(stat) {
  ## a simple function to determine whether a given statistic (string) is
  ## a utility (we want to maximize) or loss (we want to minimize)
  ## by the time we get here, stat should have already been validated
  return(!stat %in% c("rmse", "mse"))
}

get_nfeat_baseline <- function(object, baseline, stat, ...) {
  ## get model size that is used as a baseline in comparisons. baseline is one
  ## of 'best' or 'ref', stat is the statistic according to which the selection
  ## is done
  if (baseline == "best") {
    ## find number of features that maximizes the utility (or minimizes the
    ## loss)
    tab <- .tabulate_stats(object, stat, B = 2, ...)
    stats_table <- subset(tab, tab$size != Inf)
    ## tab <- .tabulate_stats(object, B = 2, ...)
    ## stats_table <- subset(tab, tab$delta == FALSE &
    ##   tab$statistic == stat & tab$size != Inf)
    optfun <- ifelse(is_util(stat), which.max, which.min)
    nfeat_baseline <- stats_table$size[optfun(stats_table$value)]
  } else {
    ## use reference model
    nfeat_baseline <- Inf
  }
  return(nfeat_baseline)
}

## Difference estimation using SRS-WOR sampling (Magnusson et al., 2020)
## copied from loo:::srs_diff_est and added weights
.srs_diff_est_w <- function(y_approx, y, y_idx, wobs = 1) {

  N <- length(y_approx)
  m <- length(y)
  y_approx_m <- y_approx[y_idx]
  wobs_m <- 1
  if (length(wobs) > 1) {
    wobs_m <- wobs[y_idx]
    wobs_m <- length(wobs_m) * wobs_m / sum(wobs_m)
    wobs <- length(wobs) * wobs / sum(wobs)
  }

  e_i <- y - y_approx_m
  t_pi_tilde <- sum(wobs * y_approx)
  t_pi2_tilde <- sum(wobs * y_approx^2)
  t_e <- N * mean(wobs_m * e_i)
  t_hat_epsilon <- N * mean(wobs_m * (y^2 - y_approx_m^2))

  est_list <- nlist(m, N)
  # eq (7)
  est_list$y_hat <- t_pi_tilde + t_e
  # eq (8)
  var_e_i <- .weighted_sd(e_i, w = wobs_m)^2
  est_list$v_y_hat <- N^2 * (1 - m / N) * var_e_i / m
  # eq (9) first row second `+` should be `-`
  # Supplementary material eq (6) has this correct
  # Here the variance is for sum, while in the paper the variance is for mean
  # which explains the proportional difference of 1/N
  est_list$hat_v_y <- (t_pi2_tilde + t_hat_epsilon) - # a (has been checked)
    (1 / N) * (t_e^2 - est_list$v_y_hat + 2 * t_pi_tilde * est_list$y_hat - t_pi_tilde^2) # b
  return(est_list)
}
