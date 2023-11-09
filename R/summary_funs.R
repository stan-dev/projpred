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
  if (!is.matrix(mu)) {
    stop("Unexpected structure for `mu`. Do the return values of ",
         "`proj_predfun` and `ref_predfun` have the correct structure?")
  }
  loglik <- family$ll_fun(mu, dis, y_wobs_test$y, y_wobs_test$wobs)
  if (!is.matrix(loglik)) {
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
  summ_ref <- varsel$summaries$ref
  summ_sub <- varsel$summaries$sub

  if (!varsel$refmodel$family$for_latent && !resp_oscale) {
    stop("`resp_oscale = FALSE` can only be used in case of the latent ",
         "projection.")
  }
  if (varsel$refmodel$family$for_latent) {
    if (resp_oscale) {
      summ_ref <- summ_ref$oscale
      summ_sub <- lapply(summ_sub, "[[", "oscale")
      ref_lppd_NA <- all(is.na(summ_ref$lppd))
      sub_lppd_NA <- any(sapply(summ_sub, check_sub_NA, el_nm = "lppd"))
      ref_mu_NA <- all(is.na(summ_ref$mu))
      sub_mu_NA <- any(sapply(summ_sub, check_sub_NA, el_nm = "mu"))
      if (ref_mu_NA || sub_mu_NA) {
        message(
          "`latent_ilink` returned only `NA`s, so all performance statistics ",
          "will also be `NA` as long as `resp_oscale = TRUE`."
        )
      } else if (any(stats %in% c("elpd", "mlpd")) &&
                 (ref_lppd_NA || sub_lppd_NA)) {
        message(
          "`latent_ll_oscale` returned only `NA`s, so ELPD and MLPD will also ",
          "be `NA` as long as `resp_oscale = TRUE`."
        )
      }
      varsel$y_wobs_test$y <- varsel$y_wobs_test$y_oscale
    } else {
      if (all(is.na(varsel$refmodel$dis)) &&
          any(stats %in% c("elpd", "mlpd"))) {
        message(
          "Cannot calculate ELPD or MLPD if `resp_oscale = FALSE` and ",
          "`<refmodel>$dis` consists of only `NA`s. If it's not possible to ",
          "supply a suitable argument `dis` to init_refmodel(), consider (i) ",
          "switching to `resp_oscale = TRUE` (which might require the ",
          "specification of functions needed by extend_family()) or (ii) ",
          "using a performance statistic other than ELPD or MLPD."
        )
      }
      if (all(is.na(varsel$y_wobs_test$y))) {
        message(
          "Cannot calculate performance statistics if `resp_oscale = FALSE` ",
          "and `<vsel>$y_wobs_test$y` consists of only `NA`s. The reason for ",
          "these `NA`s is probably that `<vsel>` was created by cv_varsel() ",
          "with `cv_method = \"kfold\"`. (In case of K-fold cross-validation, ",
          "the latent response values for the test datasets cannot be defined ",
          "in a straightforward manner without inducing dependencies between ",
          "training and test datasets.)"
        )
      }
    }
  }
  # Just to avoid that `$y` gets expanded to `$y_oscale` if element `y` does not
  # exist (for whatever reason; actually, it should always exist):
  varsel$y_wobs_test$y_oscale <- NULL

  if (resp_oscale && !is.null(varsel$refmodel$family$cats) &&
      any(stats %in% c("acc", "pctcorr"))) {
    summ_ref$mu <- catmaxprb(summ_ref$mu, lvls = varsel$refmodel$family$cats)
    summ_sub <- lapply(summ_sub, function(summ_sub_k) {
      summ_sub_k$mu <- catmaxprb(summ_sub_k$mu,
                                 lvls = varsel$refmodel$family$cats)
      return(summ_sub_k)
    })
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
    res <- get_stat(summ$mu, summ$lppd, varsel$y_wobs_test, stat, mu.bs = mu.bs,
                    lppd.bs = lppd.bs, wcv = summ$wcv, alpha = alpha, ...)
    row <- data.frame(
      data = varsel$type_test, size = Inf, delta = delta, statistic = stat,
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
        res_ref <- get_stat(summ_ref$mu, summ_ref$lppd, varsel$y_wobs_test,
                            stat, mu.bs = NULL, lppd.bs = NULL,
                            wcv = summ_ref$wcv, alpha = alpha, ...)
        res_diff <- get_stat(summ$mu, summ$lppd, varsel$y_wobs_test, stat,
                             mu.bs = summ_ref$mu, lppd.bs = summ_ref$lppd,
                             wcv = summ$wcv, alpha = alpha, ...)
        val <- res_ref$value + res_diff$value
        # TODO (subsampled PSIS-LOO CV): Is `val.se` really computed correctly
        # or do we need to take into account that `res_ref$se` and `res_diff$se`
        # might be stochastically dependent?
        val.se <- sqrt(res_ref$se^2 + res_diff$se^2)
        if (stat %in% c("rmse", "auc")) {
          # TODO (subsampled PSIS-LOO CV): Use bootstrap for lower and upper
          # confidence interval bounds as well as for the standard error.
          warning("Lower and upper confidence interval bounds of performance ",
                  "statistic `", stat, "` are based on a normal ",
                  "approximation, not the bootstrap. The standard error of ",
                  "performance statistic `", stat, "` is also not based on a ",
                  "bootstrap.")
        }
        lq <- qnorm(alpha / 2, mean = val, sd = val.se)
        uq <- qnorm(1 - alpha / 2, mean = val, sd = val.se)
        row <- data.frame(
          data = varsel$type_test, size = k - 1, delta = delta,
          statistic = stat, value = val, lq = lq, uq = uq, se = val.se,
          diff = res_diff$value, diff.se = res_diff$se
        )
      } else {
        ## normal case
        res <- get_stat(summ$mu, summ$lppd, varsel$y_wobs_test, stat,
                        mu.bs = mu.bs, lppd.bs = lppd.bs, wcv = summ$wcv,
                        alpha = alpha, ...)
        diff <- get_stat(summ$mu, summ$lppd, varsel$y_wobs_test, stat,
                         mu.bs = summ_ref$mu, lppd.bs = summ_ref$lppd,
                         wcv = summ$wcv, alpha = alpha, ...)
        row <- data.frame(
          data = varsel$type_test, size = k - 1, delta = delta,
          statistic = stat, value = res$value, lq = res$lq, uq = res$uq,
          se = res$se, diff = diff$value, diff.se = diff$se
        )
      }
      stat_tab <- rbind(stat_tab, row)
    }
  }

  return(stat_tab)
}

# Helper function checking whether all entries of a summaries vector are `NA`.
#
# @param summ_sub_k Typically `<vsel_object>$summaries$sub[[k]]`.
# @param el_nm A single character string, giving the name of the subelement of
#   `summ_sub_k` to check for `NA`s.
#
# @return A single logical value, indicating whether all entries of
#   `summ_sub_k[[el_nm]]` are `NA`.
check_sub_NA <- function(summ_sub_k, el_nm) {
  all(is.na(summ_sub_k[[el_nm]]))
}

## Calculates given statistic stat with standard error and confidence bounds.
## mu.bs and lppd.bs are the pointwise mu and lppd for another model that is
## used as a baseline for computing the difference in the given statistic, for
## example the relative elpd. If these arguments are not given (NULL) then the
## actual (non-relative) value is computed. NOTE: Element `wcv[i]` (with i = 1,
## ..., N and N denoting the number of observations) contains the weight of the
## CV fold that observation i is in. In case of varsel() output, this is `NULL`.
## Currently, these `wcv` are nonconstant (and not `NULL`) only in case of
## subsampled PSIS-LOO CV. The actual observation weights (specified by the
## user) are contained in `y_wobs_test$wobs`. These are already taken into
## account by `<refmodel_object>$family$ll_fun()` (or
## `<refmodel_object>$family$latent_ll_oscale()`) and are thus already taken
## into account in `lppd`. However, `mu` does not take them into account, so
## some further adjustments are necessary below.
get_stat <- function(mu, lppd, y_wobs_test, stat, mu.bs = NULL, lppd.bs = NULL,
                     wcv = NULL, alpha = 0.1, ...) {
  n_notna.bs <- NULL
  if (stat %in% c("mlpd", "elpd")) {
    if (!is.null(lppd.bs)) {
      # Compute the performance statistics using only those observations for
      # which both `lppd` and `lppd.bs` are not `NA`:
      lppd[is.na(lppd.bs)] <- NA
      lppd.bs[is.na(lppd)] <- NA
      n_notna.bs <- sum(!is.na(lppd.bs))
    }
    n_notna <- sum(!is.na(lppd))
    n <- length(lppd)
  } else {
    hasNA_y <- is.na(y_wobs_test$y_prop %||% y_wobs_test$y)
    if (!is.null(mu.bs)) {
      # Compute the performance statistics using only those observations for
      # which both `mu` and `mu.bs` are not `NA`:
      mu[is.na(mu.bs)] <- NA
      mu.bs[is.na(mu)] <- NA
      n_notna.bs <- sum(!is.na(mu.bs) & !hasNA_y)
    }
    n_notna <- sum(!is.na(mu) & !hasNA_y)
    n <- length(mu)
  }
  if (!is.null(n_notna.bs) && getOption("projpred.additional_checks", FALSE)) {
    stopifnot(n_notna == n_notna.bs)
  }
  if (n_notna == 0) {
    return(list(value = NA, se = NA, lq = NA, uq = NA))
  }

  if (is.null(wcv)) {
    ## set default CV fold weights if not given
    wcv <- rep(1, n)
  }
  ## ensure the CV fold weights sum to n_notna
  wcv <- n_notna * wcv / sum(wcv)

  alpha_half <- alpha / 2
  one_minus_alpha_half <- 1 - alpha_half

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
    y <- y_wobs_test$y_prop %||% y_wobs_test$y
    if (!all(y_wobs_test$wobs == 1)) {
      wcv <- wcv * y_wobs_test$wobs
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
        value <- sqrt(mean(wcv * (mu - y)^2, na.rm = TRUE)) -
          sqrt(mean(wcv * (mu.bs - y)^2, na.rm = TRUE))
        diffvalue.bootstrap <- bootstrap(
          cbind((mu - y)^2, (mu.bs - y)^2),
          function(resid2) {
            sqrt(mean(wcv * resid2[, 1], na.rm = TRUE)) -
              sqrt(mean(wcv * resid2[, 2], na.rm = TRUE))
          },
          ...
        )
        value.se <- sd(diffvalue.bootstrap)
        lq_uq <- quantile(diffvalue.bootstrap,
                          probs = c(alpha_half, one_minus_alpha_half),
                          names = FALSE, na.rm = TRUE)
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
        lq_uq <- quantile(value.bootstrap,
                          probs = c(alpha_half, one_minus_alpha_half),
                          names = FALSE, na.rm = TRUE)
      }
    }
  } else if (stat %in% c("acc", "pctcorr", "auc")) {
    y <- y_wobs_test$y
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
      if (!is.null(mu.bs)) {
        mu.bs <- rep(mu.bs, y_wobs_test$wobs)
        # CAUTION: If `y` is allowed to have `NA`s here, then `n_notna.bs` needs
        # to be adapted:
        n_notna.bs <- sum(!is.na(mu.bs))
      }
      # CAUTION: If `y` is allowed to have `NA`s here, then `n_notna` needs to
      # be adapted:
      n_notna <- sum(!is.na(mu))
      if (!is.null(n_notna.bs) &&
          getOption("projpred.additional_checks", FALSE)) {
        stopifnot(n_notna == n_notna.bs)
      }
      wcv <- rep(wcv, y_wobs_test$wobs)
      wcv <- n_notna * wcv / sum(wcv)
    } else {
      stopifnot(all(y_wobs_test$wobs == 1))
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
      if (!is.null(mu.bs)) {
        auc.data <- cbind(y, mu, wcv)
        auc.data.bs <- cbind(y, mu.bs, wcv)
        value <- auc(auc.data) - auc(auc.data.bs)
        idxs_cols <- seq_len(ncol(auc.data))
        idxs_cols_bs <- setdiff(seq_len(ncol(auc.data) + ncol(auc.data.bs)),
                                idxs_cols)
        diffvalue.bootstrap <- bootstrap(
          cbind(auc.data, auc.data.bs),
          function(x) {
            auc(x[, idxs_cols, drop = FALSE]) -
              auc(x[, idxs_cols_bs, drop = FALSE])
          },
          ...
        )
        value.se <- sd(diffvalue.bootstrap, na.rm = TRUE)
        lq_uq <- quantile(diffvalue.bootstrap,
                          probs = c(alpha_half, one_minus_alpha_half),
                          names = FALSE, na.rm = TRUE)
      } else {
        auc.data <- cbind(y, mu, wcv)
        value <- auc(auc.data)
        value.bootstrap <- bootstrap(auc.data, auc, ...)
        value.se <- sd(value.bootstrap, na.rm = TRUE)
        lq_uq <- quantile(value.bootstrap,
                          probs = c(alpha_half, one_minus_alpha_half),
                          names = FALSE, na.rm = TRUE)
      }
    }
  }

  if (!stat %in% c("rmse", "auc")) {
    lq <- qnorm(alpha_half, mean = value, sd = value.se)
    uq <- qnorm(one_minus_alpha_half, mean = value, sd = value.se)
  } else {
    lq <- lq_uq[1]
    uq <- lq_uq[2]
  }

  return(list(value = value, se = value.se, lq = lq, uq = uq))
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
