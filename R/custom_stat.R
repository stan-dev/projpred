#' Custom predictive performance statistic
#'
#' Defines a predictive performance statistic (a utility or a loss) from a user
#' function. Use the result in argument `stats` of [summary.vsel()],
#' [plot.vsel()] and [performances()], or in argument `stat` of
#' [suggest_size()]. The statistic affects the evaluation of the submodels and
#' the suggested size, but not the search.
#'
#' @param fun A function with arguments `y`, `mu` and `wobs` (see section
#'   "Details").
#' @param name A single character string. The output of [summary.vsel()] and
#'   [performances()] uses it in the column names and adds suffixes such as
#'   `.se`, so it must start with a letter and contain only letters, digits and
#'   underscores. It must not be the name of a built-in statistic or of another
#'   column of these tables (`"size"`, `"ranking_fulldata"`,
#'   `"cv_proportions_diag"`).
#' @param type Either `"pointwise"` or `"aggregate"` (see section "Details").
#' @param higher_better `TRUE` for a utility (higher values are better), `FALSE`
#'   for a loss. `baseline = "best"` and [suggest_size()] use it.
#' @param label A single character string for the facet label in
#'   [plot.vsel()].
#'
#' @details Function `fun` gets these arguments for each model (the reference
#'   model and each submodel size):
#'   * `y`: the observed response values of the evaluation (test) set. For the
#'   [binomial()] family, `y` holds the proportions of successes (the number of
#'   successes divided by the number of trials). For the latent projection,
#'   `y` holds the response values on the original response scale if
#'   `resp_oscale = TRUE` and the latent response values if
#'   `resp_oscale = FALSE`.
#'   * `mu`: the predictive means for the observations in `y`, averaged over
#'   the parameter draws of the reference model or over the projected draws of
#'   the submodel. For the [binomial()] family, these are the predicted
#'   probabilities of success. With cross-validation, `mu` holds the
#'   cross-validated predictions. For the latent projection, `mu` is on the
#'   same scale as `y`. `mu` can contain `NA`s. With subsampled LOO CV (see
#'   argument `nloo` of [cv_varsel()]), the submodel predictions are `NA` for
#'   the observations outside the subsample. For a reference model of class
#'   `datafit`, all reference model predictions are `NA`.
#'   * `wobs`: the observation weights (for the [binomial()] family: the
#'   number of trials).
#'
#'   For `type = "pointwise"`, `fun` returns a numeric (or logical) vector of
#'   the same length as `y` with the utility (or loss) of each observation.
#'   Element `i` must depend only on `y[i]`, `mu[i]` and `wobs[i]`. Do not apply
#'   `wobs` in `fun`: the statistic is the mean of these values, weighted by
#'   `wobs`. As for `"mse"`, the standard error treats each row of the
#'   evaluation set as one observation. For the [binomial()] family, a row with
#'   several trials thus counts once, whereas `"acc"` splits it into single
#'   trials. The difference from the baseline model pairs the two models
#'   observation by observation. The uncertainty intervals use a normal
#'   approximation. Subsampled LOO CV works with this type.
#'
#'   For `type = "aggregate"`, `fun` returns a single number for the whole
#'   evaluation set. As for `"auc"`, the standard error and the uncertainty
#'   intervals come from a bootstrap over the rows of the evaluation set (see
#'   argument `...` of [summary.vsel()] for the number of bootstrap samples and
#'   the seed). The difference from the baseline model uses the same bootstrap
#'   samples for both models. For the [binomial()] family, a row with several
#'   trials counts once, whereas `"auc"` splits it into single trials. This type
#'   does not support subsampled LOO CV.
#'
#'   Custom statistics are not available for the augmented-data projection, or
#'   for the latent projection with `resp_oscale = TRUE` and a non-`NULL`
#'   element `family$cats`.
#'
#' @return An object of class `custom_stat`.
#'
#' @examplesIf requireNamespace("rstanarm", quietly = TRUE)
#' # Net benefit at the decision threshold `thres` (a utility). Each observation
#' # contributes one pointwise term:
#' net_benefit <- function(thres) {
#'   custom_stat(
#'     function(y, mu, wobs) {
#'       treat <- mu >= thres
#'       treat * (y - (1 - y) * thres / (1 - thres))
#'     },
#'     name = paste0("nb", 100 * thres),
#'     label = paste0("Net benefit \u2265 ", thres)
#'   )
#' }
#'
#' # Data:
#' dat_binom <- data.frame(y = df_binom$y, df_binom$x)
#'
#' # The `stanreg` fit which will be used as the reference model (with small
#' # values for `chains` and `iter`, but only for technical reasons in this
#' # example; this is not recommended in general):
#' fit <- rstanarm::stan_glm(
#'   y ~ X1 + X2 + X3 + X4 + X5, family = binomial(), data = dat_binom,
#'   chains = 2, iter = 500, refresh = 0, seed = 9876
#' )
#'
#' # Run varsel() (here without cross-validation and with small values for
#' # `nterms_max` and `nclusters_pred`, but only for the sake of speed in this
#' # example; this is not recommended in general):
#' vs <- varsel(fit, nterms_max = 3, nclusters_pred = 10, seed = 5555)
#' print(summary(vs, stats = list("elpd", net_benefit(0.2))), digits = 2)
#'
#' @export
custom_stat <- function(fun, name, type = c("pointwise", "aggregate"),
                        higher_better = TRUE, label = name) {
  if (!is.function(fun)) {
    stop("Argument `fun` must be a function.")
  }
  if (!all(c("y", "mu", "wobs") %in% names(formals(fun))) &&
      !"..." %in% names(formals(fun))) {
    stop("Function `fun` must have arguments `y`, `mu` and `wobs`.")
  }
  if (!is.character(name) || length(name) != 1 || is.na(name) ||
      !grepl("^[A-Za-z][A-Za-z0-9_]*$", name)) {
    stop("Argument `name` must be a single character string that starts ",
         "with a letter and contains only letters, digits and underscores.")
  }
  reserved <- c(builtin_stats(), "size", "ranking_fulldata",
                "cv_proportions_diag")
  if (name %in% reserved) {
    stop("Argument `name` must not be one of ",
         paste0("`\"", reserved, "\"`", collapse = ", "), ".")
  }
  type <- match.arg(type)
  if (!isTRUE(higher_better) && !isFALSE(higher_better)) {
    stop("Argument `higher_better` must be `TRUE` or `FALSE`.")
  }
  if (!is.character(label) || length(label) != 1 || is.na(label)) {
    stop("Argument `label` must be a single character string.")
  }
  return(structure(nlist(fun, name, type, higher_better, label),
                   class = "custom_stat"))
}

#' @noRd
#' @export
print.custom_stat <- function(x, ...) {
  cat("Custom performance statistic `", x$name, "` (", x$type, ", ",
      if (x$higher_better) "higher" else "lower", " is better)\n", sep = "")
  return(invisible(x))
}

builtin_stats <- function() {
  c("elpd", "mlpd", "gmpd", "mse", "rmse", "R2", "acc", "pctcorr", "auc")
}

# Split argument `stats` into the names of the statistics and a named list of
# the `custom_stat` objects. `stats` is an atomic vector (normally character),
# one `custom_stat` object, or a list of single strings and `custom_stat`
# objects.
parse_stats <- function(stats) {
  if (is.null(stats)) {
    stop("Statistic specified as NULL.")
  }
  if (inherits(stats, "custom_stat")) {
    stats <- list(stats)
  }
  if (is.atomic(stats)) {
    # validate_vsel_object_stats() reports invalid strings and non-character
    # input such as `NA`:
    return(list(nms = as.character(stats), custom = list()))
  }
  if (!is.list(stats) || length(stats) == 0) {
    stop("Argument `stats` must be a character vector, an object created by ",
         "custom_stat() or a list of these.")
  }
  is_custom <- vapply(stats, inherits, logical(1), what = "custom_stat")
  is_chr <- vapply(stats, function(s) {
    is.character(s) && length(s) == 1
  }, logical(1))
  if (!all(is_custom | is_chr)) {
    stop("Each element of a list `stats` must be a single character string ",
         "or an object created by custom_stat().")
  }
  nms <- vapply(stats, function(s) {
    if (inherits(s, "custom_stat")) s$name else s
  }, character(1))
  custom <- stats[is_custom]
  names(custom) <- nms[is_custom]
  if (anyDuplicated(names(custom))) {
    stop("The custom statistics in `stats` must have different names.")
  }
  return(list(nms = unname(nms), custom = custom))
}

# Compute a custom statistic (arguments as for get_stat()). The pointwise type
# follows the `acc` code of get_stat() but does not split binomial rows into
# single trials. The aggregate type follows the `auc` code.
get_stat_custom <- function(cstat, mu, mu_baseline, mu_fast, loo_inds,
                            y_wobs_test, alpha, ...) {
  y <- y_wobs_test$y_prop %||% y_wobs_test$y
  wobs_orig <- y_wobs_test$wobs
  n_full <- length(mu)
  n_loo <- if (is.null(loo_inds)) n_full else length(loo_inds)
  call_fun <- function(mu_k, y_k = y, wobs_k = wobs_orig) {
    out <- cstat$fun(y = y_k, mu = mu_k, wobs = wobs_k)
    if (is.logical(out)) {
      out <- as.numeric(out)
    }
    return(out)
  }
  if (cstat$type == "pointwise") {
    wobs <- n_full * wobs_orig / sum(wobs_orig)
    check_pw <- function(u) {
      if (!is.numeric(u) || length(u) != n_full) {
        stop("The function of custom statistic `\"", cstat$name, "\"` must ",
             "return a numeric vector of length ", n_full, " (the number of ",
             "observations).")
      }
      return(u)
    }
    u <- check_pw(call_fun(mu))
    if (!is.null(mu_baseline)) {
      u_baseline <- check_pw(call_fun(mu_baseline))
    } else {
      u_baseline <- 0
    }
    if (n_loo < n_full) {
      # subsampling difference estimator (Magnusson et al., 2020)
      u_fast <- check_pw(call_fun(mu_fast))
      srs_diffe <- .srs_diff_est_w(y_approx = u_fast - u_baseline,
                                   y = (u - u_baseline)[loo_inds],
                                   y_idx = loo_inds,
                                   wobs = wobs)
      value <- srs_diffe$y_hat / n_full
      # combine estimates of var(y_hat) and var(y)
      value_se <- sqrt(srs_diffe$v_y_hat + srs_diffe$hat_v_y) / n_full
    } else {
      value <- mean(wobs * u) - mean(wobs * u_baseline)
      value_se <- .weighted_sd(u - u_baseline, wobs) / sqrt(n_full)
    }
    lq <- qnorm(alpha / 2, mean = value, sd = value_se)
    uq <- qnorm(1 - alpha / 2, mean = value, sd = value_se)
  } else {
    if (n_loo < n_full) {
      stop("Subsampled LOO-CV with an aggregate custom statistic is not ",
           "implemented.")
    }
    check_agg <- function(v) {
      if (!is.numeric(v) || length(v) != 1) {
        stop("The function of custom statistic `\"", cstat$name, "\"` must ",
             "return a single number.")
      }
      return(v)
    }
    stat_data <- cbind(y, mu, wobs_orig)
    if (!is.null(mu_baseline)) {
      stat_data <- cbind(stat_data, mu_baseline)
      fun_bs <- function(x) {
        call_fun(x[, 2], y_k = x[, 1], wobs_k = x[, 3]) -
          call_fun(x[, 4], y_k = x[, 1], wobs_k = x[, 3])
      }
      value <- check_agg(call_fun(mu)) - check_agg(call_fun(mu_baseline))
    } else {
      fun_bs <- function(x) {
        call_fun(x[, 2], y_k = x[, 1], wobs_k = x[, 3])
      }
      value <- check_agg(call_fun(mu))
    }
    value_bs <- bootstrap(stat_data, fun_bs, ...)
    value_se <- sd(value_bs)
    if (any(is.na(value_bs))) {
      # quantile() is not able to deal with `NA`s
      lq_uq <- rep(NA_real_, 2)
    } else {
      lq_uq <- quantile(value_bs, probs = c(alpha / 2, 1 - alpha / 2),
                        names = FALSE)
    }
    lq <- lq_uq[1]
    uq <- lq_uq[2]
  }
  return(list(value = value, se = value_se, lq = lq, uq = uq))
}
