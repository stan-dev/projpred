.get_sub_summaries <- function(submodels, refmodel, test_points, newdata = NULL,
                               offset = refmodel$offset[test_points],
                               wobs = refmodel$wobs[test_points],
                               y = refmodel$y[test_points],
                               yOrig = refmodel$yOrig[test_points]) {
  lapply(submodels, function(initsubmodl) {
    .weighted_summary_means(
      y_test = list(y = y, yOrig = yOrig, weights = wobs),
      family = refmodel$family,
      wsample = initsubmodl$weights,
      mu = refmodel$family$mu_fun(initsubmodl$submodl, obs = test_points,
                                  newdata = newdata, offset = offset),
      dis = initsubmodl$dis,
      cl_ref = initsubmodl$cl_ref,
      wdraws_ref = initsubmodl$wdraws_ref
    )
  })
}

# Calculate log predictive density values and average them across parameter
# draws (together with the corresponding expected response values).
#
# @param y_test A `list`, at least with elements `y` (response values) and
#   `weights` (observation weights). In case of the latent projection, this
#   `list` also needs to contain `yOrig` (response values on the original
#   response scale, i.e., the non-latent response values).
# @param family A `family` object.
# @param wsample A vector of weights for the parameter draws.
# @param mu A matrix of expected values for `y`.
# @param dis A vector of dispersion parameter draws.
# @param cl_ref A numeric vector of length \eqn{S} (with \eqn{S} denoting the
#   number of parameter draws in the reference model), giving the cluster
#   indices for the parameter draws in the reference model. Draws that should be
#   dropped (e.g., because of thinning by `ndraws` or `ndraws_pred`) need to
#   have an `NA` in `cl_ref`. Caution: This always refers to the reference
#   model's parameter draws, not necessarily to the columns of `mu`, the entries
#   of `wsample`, or the entries of `dis`!
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
.weighted_summary_means <- function(y_test, family, wsample, mu, dis, cl_ref,
                                    wdraws_ref = rep(1, length(cl_ref))) {
  if (!is.matrix(mu)) {
    stop("Unexpected structure for `mu`. Do the return values of ",
         "`proj_predfun` and `ref_predfun` have the correct structure?")
  }
  loglik <- family$ll_fun(mu, dis, y_test$y, y_test$weights)
  if (!is.matrix(loglik)) {
    stop("Unexpected structure for `loglik`. Please notify the package ",
         "maintainer.")
  }
  # Average over the draws, taking their weights into account:
  avg <- list(
    mu = structure(c(mu %*% wsample),
                   nobs_orig = attr(mu, "nobs_orig"),
                   class = sub("augmat", "augvec", oldClass(mu), fixed = TRUE)),
    lppd = apply(loglik, 1, log_weighted_mean_exp, wsample)
  )
  if (family$for_latent && family$lat2resp_possible) {
    mu_Orig <- family$latent_ilink(t(mu), cl_ref = cl_ref,
                                   wdraws_ref = wdraws_ref)
    if (length(dim(mu_Orig)) < 2) {
      stop("Unexpected structure for `mu_Orig`. Does the return value of ",
           "`latent_ilink` have the correct structure?")
    }
    loglik_Orig <- family$latent_llOrig(mu_Orig, yOrig = y_test$yOrig,
                                        wobs = y_test$weights)
    if (!is.matrix(loglik_Orig)) {
      stop("Unexpected structure for `loglik_Orig`. Does the return value of ",
           "`latent_llOrig` have the correct structure?")
    }
    if (length(dim(mu_Orig)) == 3) {
      # In this case, `mu_Orig` is a 3-dimensional array (S x N x C), so coerce
      # it to an augmented-rows matrix:
      mu_Orig <- arr2augmat(mu_Orig, margin_draws = 1)
      mu_Orig_avg <- structure(
        c(mu_Orig %*% wsample),
        nobs_orig = attr(mu_Orig, "nobs_orig"),
        class = sub("augmat", "augvec", oldClass(mu_Orig), fixed = TRUE)
      )
    } else {
      # In principle, we could use the same code for `mu_Orig_avg` as above.
      # However, that would require `mu_Orig <- t(mu_Orig)` beforehand, so the
      # following should be more efficient:
      mu_Orig_avg <- c(wsample %*% mu_Orig)
    }
    avg$resp <- list(
      mu = mu_Orig_avg,
      lppd = apply(loglik_Orig, 2, log_weighted_mean_exp, wsample)
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
                            nfeat_baseline = NULL, lat2resp = FALSE, ...) {
  stat_tab <- data.frame()
  summ_ref <- varsel$summaries$ref
  summ_sub <- varsel$summaries$sub
  if (!varsel$refmodel$family$for_latent && lat2resp) {
    stop("`lat2resp = TRUE` can only be used in case of the latent projection.")
  }
  if (lat2resp) {
    summ_sub_Orig <- lapply(summ_sub, "[[", "Orig")
    # `lat2resp = TRUE` only makes sense if element `"Orig"` is available:
    if (is.null(summ_ref$resp) || any(sapply(summ_sub_Orig, is.null))) {
      stop("Cannot calculate the performance statistics on response scale if ",
           "`latent_ilink` or `latent_llOrig` are missing. Use ",
           "`lat2resp = FALSE` or provide the missing functions when creating ",
           "the reference model (see the documentation of extend_family() ",
           "which is called by init_refmodel()).")
    }
    summ_ref <- summ_ref$resp
    summ_sub <- summ_sub_Orig
  }
  if ((!varsel$refmodel$family$for_latent || lat2resp) &&
      !is.null(varsel$refmodel$family$cats) &&
      any(stats %in% c("acc", "pctcorr"))) {
    summ_ref$mu <- catmaxprb(summ_ref$mu, lvls = varsel$refmodel$family$cats)
    summ_sub <- lapply(summ_sub, function(summ_sub_k) {
      summ_sub_k$mu <- catmaxprb(summ_sub_k$mu,
                                 lvls = varsel$refmodel$family$cats)
      return(summ_sub_k)
    })
    # Since `mu` is an unordered factor, `y` needs to be unordered, too (or both
    # would need to be ordered; however, unordered is the simpler type):
    varsel$d_test$y <- factor(varsel$d_test$y, ordered = FALSE)
    varsel$d_test$yOrig <- factor(varsel$d_test$yOrig, ordered = FALSE)
  }
  if (lat2resp) {
    varsel$d_test$y <- varsel$d_test$yOrig
  }
  # Just to avoid that `$y` gets expanded to `$yOrig` if element `"y"` does not
  # exist (for whatever reason; actually, it should always exist):
  varsel$d_test$yOrig <- NULL

  if (varsel$refmodel$family$family == "binomial" &&
      !all(varsel$d_test$weights == 1)) {
    # This case should not occur (yet) for the augmented-data or the latent
    # projection:
    stopifnot(!varsel$refmodel$family$for_augdat)
    stopifnot(!varsel$refmodel$family$for_latent)
    varsel$d_test$y_prop <- varsel$d_test$y / varsel$d_test$weights
  }

  ## fetch the mu and lppd for the baseline model
  if (is.null(nfeat_baseline)) {
    ## no baseline model, i.e, compute the statistics on the actual
    ## (non-relative) scale
    mu.bs <- NULL
    lppd.bs <- NULL
    delta <- FALSE
  } else {
    if (nfeat_baseline == Inf) {
      summ.bs <- summ_ref
    } else {
      summ.bs <- summ_sub[[nfeat_baseline + 1]]
    }
    mu.bs <- summ.bs$mu
    lppd.bs <- summ.bs$lppd
    delta <- TRUE
  }

  for (s in seq_along(stats)) {
    stat <- stats[s]

    ## reference model statistics
    summ <- summ_ref
    res <- get_stat(summ$mu, summ$lppd, varsel$d_test, stat, mu.bs = mu.bs,
                    lppd.bs = lppd.bs, wcv = summ$wcv, alpha = alpha, ...)
    row <- data.frame(
      data = varsel$d_test$type, size = Inf, delta = delta, statistic = stat,
      value = res$value, lq = res$lq, uq = res$uq, se = res$se, diff = NA,
      diff.se = NA
    )
    stat_tab <- rbind(stat_tab, row)

    ## submodel statistics
    for (k in seq_along(summ_sub)) {
      summ <- summ_sub[[k]]
      if (delta == FALSE && sum(!is.na(summ_ref$mu)) > sum(!is.na(summ$mu))) {
        ## special case (subsampling loo): reference model summaries computed
        ## for more points than for the submodel, so utilize the reference model
        ## results to get more accurate statistic fot the submodel on the actual
        ## scale
        res_ref <- get_stat(summ_ref$mu, summ_ref$lppd, varsel$d_test,
                            stat, mu.bs = NULL, lppd.bs = NULL,
                            wcv = summ_ref$wcv, alpha = alpha, ...)
        res_diff <- get_stat(summ$mu, summ$lppd, varsel$d_test, stat,
                             mu.bs = summ_ref$mu, lppd.bs = summ_ref$lppd,
                             wcv = summ$wcv, alpha = alpha, ...)
        val <- res_ref$value + res_diff$value
        val.se <- sqrt(res_ref$se^2 + res_diff$se^2)
        lq <- qnorm(alpha / 2, mean = val, sd = val.se)
        uq <- qnorm(1 - alpha / 2, mean = val, sd = val.se)
        row <- data.frame(
          data = varsel$d_test$type, size = k - 1, delta = delta,
          statistic = stat, value = val, lq = lq, uq = uq, se = val.se,
          diff = res_diff$value, diff.se = res_diff$se
        )
      } else {
        ## normal case
        res <- get_stat(summ$mu, summ$lppd, varsel$d_test, stat, mu.bs = mu.bs,
                        lppd.bs = lppd.bs, wcv = summ$wcv, alpha = alpha, ...)
        diff <- get_stat(summ$mu, summ$lppd, varsel$d_test, stat,
                         mu.bs = summ_ref$mu, lppd.bs = summ_ref$lppd,
                         wcv = summ$wcv, alpha = alpha, ...)
        row <- data.frame(
          data = varsel$d_test$type, size = k - 1, delta = delta,
          statistic = stat, value = res$value, lq = res$lq, uq = res$uq,
          se = res$se, diff = diff$value, diff.se = diff$se
        )
      }
      stat_tab <- rbind(stat_tab, row)
    }
  }

  return(stat_tab)
}

## Calculates given statistic stat with standard error and confidence bounds.
## mu.bs and lppd.bs are the pointwise mu and lppd for another model that is
## used as a baseline for computing the difference in the given statistic,
## for example the relative elpd. If these arguments are not given (NULL) then
## the actual (non-relative) value is computed.
## NOTE: Element `wcv[i]` (with i = 1, ..., N and N denoting the number of
## observations) contains the weight of the CV fold that observation i is in. In
## case of varsel() output, this is `NULL`. Currently, these `wcv` are
## nonconstant (and not `NULL`) only in case of subsampled LOO CV. The actual
## observation weights (specified by the user) are contained in
## `d_test$weights`. These are already taken into account by
## `<refmodel_object>$family$ll_fun()` (or
## `<refmodel_object>$family$latent_llOrig()`) and are thus already taken into
## account in `lppd`. However, `mu` does not take them into account, so some
## further adjustments are necessary below.
get_stat <- function(mu, lppd, d_test, stat, mu.bs = NULL, lppd.bs = NULL,
                     wcv = NULL, alpha = 0.1, ...) {
  if (stat %in% c("mlpd", "elpd")) {
    n <- length(lppd)
    n_notna <- sum(!is.na(lppd))
  } else {
    n <- length(mu)
    n_notna <- sum(!is.na(mu))
  }

  if (is.null(wcv)) {
    ## set default CV fold weights if not given
    wcv <- rep(1, n)
  }
  ## ensure the CV fold weights sum to n_notna
  wcv <- n_notna * wcv / sum(wcv)

  if (stat %in% c("mlpd", "elpd")) {
    if (!is.null(lppd.bs)) {
      value <- sum((lppd - lppd.bs) * wcv, na.rm = TRUE)
      value.se <- weighted.sd(lppd - lppd.bs, wcv, na.rm = TRUE) *
        sqrt(n_notna)
    } else {
      value <- sum(lppd * wcv, na.rm = TRUE)
      value.se <- weighted.sd(lppd, wcv, na.rm = TRUE) *
        sqrt(n_notna)
    }
    if (stat == "mlpd") {
      value <- value / n_notna
      value.se <- value.se / n_notna
    }
  } else if (stat %in% c("mse", "rmse")) {
    if (is.null(d_test$y_prop)) {
      y <- d_test$y
    } else {
      y <- d_test$y_prop
    }
    if (!all(d_test$weights == 1)) {
      wcv <- wcv * d_test$weights
      wcv <- n_notna * wcv / sum(wcv)
    }
    if (stat == "mse") {
      if (!is.null(mu.bs)) {
        value <- mean(wcv * ((mu - y)^2 - (mu.bs - y)^2), na.rm = TRUE)
        value.se <- weighted.sd((mu - y)^2 - (mu.bs - y)^2, wcv,
                                na.rm = TRUE) /
          sqrt(n_notna)
      } else {
        value <- mean(wcv * (mu - y)^2, na.rm = TRUE)
        value.se <- weighted.sd((mu - y)^2, wcv, na.rm = TRUE) /
          sqrt(n_notna)
      }
    } else if (stat == "rmse") {
      if (!is.null(mu.bs)) {
        ## make sure the relative rmse is computed using only those points for
        ## which
        mu.bs[is.na(mu)] <- NA
        mu[is.na(mu.bs)] <- NA # both mu and mu.bs are non-NA
        value <- sqrt(mean(wcv * (mu - y)^2, na.rm = TRUE)) -
          sqrt(mean(wcv * (mu.bs - y)^2, na.rm = TRUE))
        value.bootstrap1 <- bootstrap(
          (mu - y)^2,
          function(resid2) {
            sqrt(mean(wcv * resid2, na.rm = TRUE))
          },
          ...
        )
        value.bootstrap2 <- bootstrap(
          (mu.bs - y)^2,
          function(resid2) {
            sqrt(mean(wcv * resid2, na.rm = TRUE))
          },
          ...
        )
        value.se <- sd(value.bootstrap1 - value.bootstrap2)
      } else {
        value <- sqrt(mean(wcv * (mu - y)^2, na.rm = TRUE))
        value.bootstrap <- bootstrap(
          (mu - y)^2,
          function(resid2) {
            sqrt(mean(wcv * resid2, na.rm = TRUE))
          },
          ...
        )
        value.se <- sd(value.bootstrap)
      }
    }
  } else if (stat %in% c("acc", "pctcorr", "auc")) {
    y <- d_test$y
    if (!is.null(d_test$y_prop)) {
      # In fact, the following stopifnot() checks should not be necessary
      # because this case should only occur for the binomial family (where
      # `d_test$weights` contains the numbers of trials) with more than 1 trial
      # for at least one observation:
      stopifnot(all(.is.wholenumber(d_test$weights)))
      stopifnot(all(.is.wholenumber(y)))
      stopifnot(all(0 <= y & y <= d_test$weights))
      y <- unlist(lapply(seq_along(y), function(i_short) {
        c(rep(0L, d_test$weights[i_short] - y[i_short]),
          rep(1L, y[i_short]))
      }))
      mu <- rep(mu, d_test$weights)
      if (!is.null(mu.bs)) {
        mu.bs <- rep(mu.bs, d_test$weights)
      }
      n_notna <- sum(d_test$weights)
      wcv <- rep(wcv, d_test$weights)
      wcv <- n_notna * wcv / sum(wcv)
    } else {
      stopifnot(all(d_test$weights == 1))
    }
    if (stat %in% c("acc", "pctcorr")) {
      # Find out whether each observation was classified correctly or not:
      if (!is.factor(mu)) {
        mu <- round(mu)
      }
      crrct <- mu == y

      if (!is.null(mu.bs)) {
        if (!is.factor(mu.bs)) {
          mu.bs <- round(mu.bs)
        }
        crrct.bs <- mu.bs == y

        value <- mean(wcv * (crrct - crrct.bs), na.rm = TRUE)
        value.se <- weighted.sd(crrct - crrct.bs, wcv, na.rm = TRUE) /
          sqrt(n_notna)
      } else {
        value <- mean(wcv * crrct, na.rm = TRUE)
        value.se <- weighted.sd(crrct, wcv, na.rm = TRUE) / sqrt(n_notna)
      }
    } else if (stat == "auc") {
      auc.data <- cbind(y, mu, wcv)
      if (!is.null(mu.bs)) {
        mu.bs[is.na(mu)] <- NA # compute the AUCs using only those points
        mu[is.na(mu.bs)] <- NA # for which both mu and mu.bs are non-NA
        auc.data.bs <- cbind(y, mu.bs, wcv)
        value <- auc(auc.data) - auc(auc.data.bs)
        value.bootstrap1 <- bootstrap(auc.data, auc, ...)
        value.bootstrap2 <- bootstrap(auc.data.bs, auc, ...)
        value.se <- sd(value.bootstrap1 - value.bootstrap2, na.rm = TRUE)
      } else {
        value <- auc(auc.data)
        value.bootstrap <- bootstrap(auc.data, auc, ...)
        value.se <- sd(value.bootstrap, na.rm = TRUE)
      }
    }
  }

  lq <- qnorm(alpha / 2, mean = value, sd = value.se)
  uq <- qnorm(1 - alpha / 2, mean = value, sd = value.se)

  return(list(value = value, se = value.se, lq = lq, uq = uq))
}

.is_util <- function(stat) {
  ## a simple function to determine whether a given statistic (string) is
  ## a utility (we want to maximize) or loss (we want to minimize)
  ## by the time we get here, stat should have already been validated
  return(!stat %in% c("rmse", "mse"))
}

.get_nfeat_baseline <- function(object, baseline, stat, ...) {
  ## get model size that is used as a baseline in comparisons. baseline is one
  ## of 'best' or 'ref', stat is the statistic according to which the selection
  ## is done
  if (baseline == "best") {
    ## find number of features that maximizes the utility (or minimizes the
    ## loss)
    tab <- .tabulate_stats(object, stat, ...)
    stats_table <- subset(tab, tab$size != Inf)
    ## tab <- .tabulate_stats(object, ...)
    ## stats_table <- subset(tab, tab$delta == FALSE &
    ##   tab$statistic == stat & tab$size != Inf)
    optfun <- ifelse(.is_util(stat), which.max, which.min)
    nfeat_baseline <- stats_table$size[optfun(stats_table$value)]
  } else {
    ## use reference model
    nfeat_baseline <- Inf
  }
  return(nfeat_baseline)
}
