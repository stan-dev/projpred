.onAttach <- function(...) {
  ver <- utils::packageVersion("projpred")
  msg <- paste0("This is projpred version ", ver, ".")
  packageStartupMessage(msg)
}

nms_d_test <- function() {
  c("type", "data", "offset", "weights", "y", "y_oscale")
}

nms_y_wobs_test <- function(wobs_nm = "wobs") {
  c("y", "y_oscale", wobs_nm)
}

.weighted_sd <- function(x, w, na.rm = FALSE) {
  if (length(w) == 1 && length(x) > 1) {
    w <- rep_len(w, length.out = length(x))
  }
  if (na.rm) {
    ind <- !is.na(w) & !is.na(x)
    n <- sum(ind)
  } else {
    n <- length(x)
    ind <- rep(TRUE, n)
  }
  if (n %in% c(0, 1)) {
    return(NA_real_)
  }
  w <- w / sum(w[ind])
  m <- sum(x[ind] * w[ind])
  sqrt(n / (n - 1) * sum(w[ind] * (x[ind] - m)^2))
}

log_weighted_mean_exp <- function(x, w) {
  log_sum_exp(x + log(w))
}

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# Evaluate `expr` and then restore the RNG state (and the RNG kind if `expr`
# changed it).
with_preserved_seed <- function(expr) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit({
    if (!identical(RNGkind(), old_kind)) {
      suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    }
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  expr
}

# Base seed for the weights of the weight perturbation, drawn from the current
# RNG stream without advancing it.
perturb_base_seed <- function() {
  with_preserved_seed(sample.int(.Machine$integer.max, 1))
}

# Observation weights n * Dirichlet(alpha, ..., alpha). They use their own seed
# and generator, so sequential and parallel runs give the same weights.
perturb_weights <- function(n, alpha, seed) {
  g <- with_preserved_seed({
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    stats::rgamma(n, alpha)
  })
  w <- n * g / sum(g)
  if (!all(is.finite(w))) {
    stop("Argument `loo_perturb_alpha` is too small: the perturbation ",
         "weights are not finite.")
  }
  w
}

# `TRUE` if all values of the numeric array `x` are finite. min() and max() do
# not copy `x`. They return NA, NaN or +-Inf if a value of `x` is not finite:
all_finite <- function(x) {
  length(x) == 0 || (is.finite(min(x)) && is.finite(max(x)))
}

# Reference model for the search of LOO fold i with weight perturbation. It
# multiplies `wobs` by w = n * Dirichlet(alpha, ..., alpha), with w_i set to 0.
# Each draw of the linear predictor of observation k moves by the first-order
# change of its posterior mean under these weights:
#   c_s = sum_{j != i} (w_j - 1) log p(y_j | theta_s),
#   cbar = sum_s omega_s c_s,
#   delta_k = sum_s omega_s (c_s - cbar) eta_ks,
# where omega are the PSIS-LOO weights of fold i (`lw_i`, log scale). The log
# of the dispersion moves in the same way. `ll` is the S x n log-likelihood
# matrix of the reference draws. Returns `NULL` if the shifted linear
# predictors, means (with and without offsets) or dispersions are not finite.
perturb_refmodel <- function(refmodel, i, ll, lw_i, alpha, seed) {
  n <- ncol(ll)
  w <- perturb_weights(n, alpha, seed)
  w[i] <- 0
  dw <- w - 1
  dw[i] <- 0
  om <- exp(lw_i - log_sum_exp(lw_i))
  cc <- drop(ll %*% dw)
  cd <- om * (cc - sum(om * cc))
  if (!all(is.finite(cd))) {
    return(NULL)
  }
  refmodel$eta <- refmodel$eta + drop(refmodel$eta %*% cd)
  # With a bounded inverse link (logit, probit), an infinite `eta` gives a
  # finite `mu`, so check `eta` itself:
  if (!all_finite(refmodel$eta)) {
    return(NULL)
  }
  refmodel$mu <- refmodel$family$linkinv(refmodel$eta)
  if (!all_finite(refmodel$mu)) {
    return(NULL)
  }
  if (!all(refmodel$offset == 0)) {
    if (refmodel$family$family %in% fams_neg_linpred()) {
      eta_offs <- refmodel$eta - refmodel$offset
    } else {
      eta_offs <- refmodel$eta + refmodel$offset
    }
    refmodel$mu_offs <- refmodel$family$linkinv(eta_offs)
    if (!all_finite(refmodel$mu_offs)) {
      return(NULL)
    }
  } else {
    refmodel$mu_offs <- refmodel$mu
  }
  if (!any(is.na(refmodel$dis)) && all(refmodel$dis > 0)) {
    ldis <- log(refmodel$dis)
    refmodel$dis <- exp(ldis + sum(ldis * cd))
  }
  if (!any(is.na(refmodel$dis)) && !all(is.finite(refmodel$dis))) {
    return(NULL)
  }
  refmodel$wobs <- refmodel$wobs * w
  refmodel
}

linkfun_raw <- function(x, link_nm) {
  if (link_nm %in% c("logistic")) {
    link_nm <- "logit"
  } else if (link_nm %in% c("probit_approx")) {
    link_nm <- "probit"
  } else if (link_nm == "cloglog") {
    # The `"cloglog"` link is also supported by binomial(), but the following
    # should be numerically more stable:
    return(log(-log1p(-x)))
  }
  basic_link <- binomial(link = link_nm)$linkfun
  return(basic_link(x))
}

ilinkfun_raw <- function(x, link_nm) {
  if (link_nm %in% c("logistic")) {
    link_nm <- "logit"
  } else if (link_nm %in% c("probit_approx")) {
    link_nm <- "probit"
  }
  basic_ilink <- binomial(link = link_nm)$linkinv
  return(basic_ilink(x))
}

.auc <- function(x) {
  resp <- x[, 1]
  pred <- x[, 2]
  wobs <- x[, 3]

  # Make it explicit that `x` should not be used anymore:
  rm(x)

  ord <- order(pred, decreasing = TRUE, na.last = FALSE)
  n <- length(ord)

  resp <- resp[ord]
  pred <- pred[ord]
  wobs <- wobs[ord]

  w0 <- w1 <- wobs
  stopifnot(all(na.omit(resp) %in% c(0, 1)))
  w0[is.na(resp)] <- NA # ensure that `NA`s in `resp` propagate to the output
  w1[is.na(resp)] <- NA # ensure that `NA`s in `resp` propagate to the output
  w0[resp == 1] <- 0 # for calculating the false positive rate (fpr)
  w1[resp == 0] <- 0 # for calculating the true positive rate (tpr)
  cum_w0 <- cumsum(w0)
  cum_w1 <- cumsum(w1)

  ## ignore tied predicted probabilities, keeping only the rightmost one
  rightmost.prob <- c(diff(pred) != 0, TRUE)
  fpr <- c(0, cum_w0[rightmost.prob]) / cum_w0[n]
  tpr <- c(0, cum_w1[rightmost.prob]) / cum_w1[n]
  delta_fpr <- c(diff(fpr), 0)
  delta_tpr <- c(diff(tpr), 0)

  ## sum the area of the rectangles that fall completely below the ROC curve
  ## plus half the area of the rectangles that are cut in two by the curve
  return(sum(delta_fpr * tpr) + sum(delta_fpr * delta_tpr) / 2)
}

# Bootstrap an arbitrary quantity `fun` that takes the sample `x` as the first
# input. Other arguments of `fun` can be passed by `...`. Example:
# `boostrap(x, mean)`.
bootstrap <- function(x, fun = mean, B = 2000, seed = NA, ...) {
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

  seq_x <- seq_len(NROW(x))
  is_vector <- NCOL(x) == 1
  bsstat <- rep(NA, B)
  for (i in 1:B) {
    bsind <- sample(seq_x, replace = TRUE)
    bsstat[i] <- fun(if (is_vector) x[bsind] else x[bsind, , drop = FALSE], ...)
  }
  return(bsstat)
}

# From `?is.integer` (slightly modified):
is_wholenumber <- function(x) {
  abs(x - round(x)) < .Machine$double.eps^0.5
}

validate_num_folds <- function(k, n) {
  if (!is.numeric(k) || length(k) != 1 || !is_wholenumber(k)) {
    stop("Number of folds must be a single integer value.")
  }
  if (k < 2) {
    stop("Number of folds must be at least 2.")
  }
  if (k > n) {
    stop("Number of folds cannot exceed n.")
  }
}

validate_vsel_object_stats <- function(object, stats, resp_oscale = TRUE,
                                       custom_stats = list()) {
  if (!inherits(object, c("vsel"))) {
    stop("The object is not a variable selection object. Run variable ",
         "selection first")
  }
  if (!object$refmodel$family$for_latent && !resp_oscale) {
    stop("`resp_oscale = FALSE` can only be used in case of the latent ",
         "projection.")
  }
  resp_oscale <- object$refmodel$family$for_latent && resp_oscale

  trad_stats <- builtin_stats()
  trad_stats_binom_only <- c("acc", "pctcorr", "auc")
  augdat_stats <- c("elpd", "mlpd", "gmpd", "acc", "pctcorr")
  resp_oscale_stats_fac <- augdat_stats

  if (is.null(stats)) {
    stop("Statistic specified as NULL.")
  }
  if (resp_oscale) {
    fam_ch <- object$refmodel$family$family_oscale
  } else {
    fam_ch <- object$refmodel$family$family
  }
  for (stat in stats) {
    if (stat %in% names(custom_stats)) {
      if (object$refmodel$family$for_augdat) {
        stop("Currently, the augmented-data projection may not be combined ",
             "with a custom performance statistic.")
      }
      if (resp_oscale && !is.null(object$refmodel$family$cats)) {
        stop("Currently, the latent projection with `resp_oscale = TRUE` and ",
             "a non-`NULL` element `family$cats` may not be combined with a ",
             "custom performance statistic.")
      }
    } else if (object$refmodel$family$for_augdat) {
      if (!stat %in% augdat_stats) {
        stop("Currently, the augmented-data projection may not be combined ",
             "with performance statistic `\"", stat, "\"`.")
      }
    } else if (resp_oscale && !is.null(object$refmodel$family$cats)) {
      if (!stat %in% resp_oscale_stats_fac) {
        stop("Currently, the latent projection with `resp_oscale = TRUE` and ",
             "a non-`NULL` element `family$cats` may not be combined with ",
             "performance statistic `\"", stat, "\"`.")
      }
    } else {
      if (!stat %in% trad_stats) {
        stop(sprintf("Statistic '%s' not recognized.", stat))
      }
      if (stat %in% trad_stats_binom_only && fam_ch != "binomial") {
        stop("In case of (i) the traditional projection or (ii) the latent ",
             "projection with `resp_oscale = TRUE` and a `NULL` element ",
             "`family$cats`, the performance statistic `\"", stat, "\"` is ",
             "available only for the binomial family. This also explains why ",
             "performance statistic `\"", stat, "\"` is not available in case ",
             "of the latent projection with `resp_oscale = FALSE` (because a ",
             "latent Gaussian distribution is used there).")
      }
    }
    if (object$refmodel$family$for_latent &&
        !is.null(attr(object$refmodel$family$latent_ll_oscale, "cens_var")) &&
        resp_oscale &&
        !stat %in% c("elpd", "mlpd", "gmpd")) {
      warning("Performance statistic `\"", stat, "\"` does not take ",
              "censoring into account.")
    }
  }
  return(invisible(TRUE))
}

validate_baseline <- function(vsel_obj, baseline, deltas) {
  stopifnot(!is.null(baseline))
  if (!(baseline %in% c("ref", "best"))) {
    stop("Argument 'baseline' must be either 'ref' or 'best'.")
  }
  if (baseline == "ref" && (identical(deltas, "mixed") || isTRUE(deltas)) &&
      inherits(vsel_obj$refmodel, "datafit")) {
    # no reference model (or the results missing for some other reason),
    # so cannot compute differences (or ratios) vs. the reference model
    stop("Cannot use `deltas = TRUE` (or `deltas = \"mixed\"`) and ",
         "`baseline = \"ref\"` when there is no reference model.")
  }
  if (baseline == "best" && vsel_obj$cv_method == "LOO" &&
      isTRUE(vsel_obj$nloo < vsel_obj$refmodel$nobs)) {
    stop("Cannot use `baseline = \"best\"` in case of subsampled LOO-CV.")
  }
  return(baseline)
}

# A function for retrieving `y` and the corresponding observation weights
# `weights` in their "standard" forms:
#   * If `NCOL(y) == 2`: `y` is the first column and `weights` the second.
#   * If `NCOL(y) == 1`: `weights` is basically unchanged (unless of length zero
#     in which case it is replaced by a vector of ones). For a binomial family,
#     if `is.factor(y)`, `y` is transformed into a zero-one vector (i.e., with
#     values in the set {0, 1}).
get_standard_y <- function(y, weights, fam) {
  if (NCOL(y) == 1) {
    if (length(weights) > 0) {
      weights <- unname(weights)
    } else {
      weights <- rep(1, length(y))
    }
    if (fam$family == "binomial") {
      if (is.factor(y) && !fam$for_augdat) {
        if (nlevels(y) > 2) {
          stop("y cannot contain more than two classes if specified as factor.")
        }
        y <- as.vector(y, mode = "integer") - 1L # zero-one vector
      }
    } else {
      if (is.factor(y) && !fam$for_augdat && !fam$for_latent) {
        stop("y cannot be a factor for models other than the binomial model.")
      }
    }
  } else if (NCOL(y) == 2) {
    if (fam$family != "binomial") {
      stop("For non-binomial families, a two-column response is not allowed.")
    }
    weights <- unname(y[, 1] + y[, 2])
    y <- unname(y[, 1])
  } else {
    stop("The response is not allowed to have more than two columns.")
  }
  return(nlist(y, weights))
}

# Create the "reference distribution", i.e., reduce the number of posterior
# draws from the reference model by clustering, thinning, or subsampling them
#
# @param refmodel An object of class `refmodel`.
# @param nclusters The desired number of clusters of draws. If
#   `!is.null(nclusters)`, then clustering is used and `ndraws` is ignored.
# @param ndraws The desired number of draws. If `!is.null(nclusters)`, then
#   clustering is used and `ndraws` is ignored.
# @param thinning A single logical value indicating whether in the case where
#   `ndraws` is used, the reference model's draws should be thinned or
#   subsampled (without replacement).
#
# @return Let \eqn{y} denote the response (vector), \eqn{N} the number of
#   observations (for the traditional or the latent projection) or the number of
#   augmented observations (for augmented-data projection), and
#   \eqn{S_{\mathrm{prj}}}{S_prj} the number of projected draws (= either
#   `nclusters` or `ndraws`, depending on which one is used). Then the return
#   value is a list with elements:
#
#   * `mu`: An \eqn{N \times S_{\mathrm{prj}}}{N x S_prj} matrix of expected
#   values for \eqn{y} (probabilities for the response categories in case of the
#   augmented-data projection) for each draw/cluster.
#   * `mu_offs`: Same as `mu`, but taking offsets into account.
#   * `var`: An \eqn{N \times S_{\mathrm{prj}}}{N x S_prj} matrix of predictive
#   variances for \eqn{y} for each draw/cluster which are needed for projecting
#   the dispersion parameter (the predictive variances are NA for those families
#   that do not have a dispersion parameter).
#   * `dis`: A vector of length \eqn{S_{\mathrm{prj}}}{S_prj} containing the
#   reference model's dispersion parameter value for each draw/cluster (NA for
#   those families that do not have a dispersion parameter).
#   * `wdraws_prj`: A vector of length \eqn{S_{\mathrm{prj}}}{S_prj} containing
#   the weights for the projected draws/clusters.
#   * `const_wdraws_prj`: A single logical value indicating whether all
#   projected draws have the same weight.
#   * `nprjdraws`: A single value: \eqn{S_{\mathrm{prj}}}{S_prj}.
#   * `cl`: Cluster assignment for each posterior draw, that is, a vector that
#   has length equal to the number of posterior draws and each value is an
#   integer between 1 and \eqn{S_{\mathrm{prj}}}{S_prj}.
#   * `wdraws_orig`: A numeric vector of length equal to the number of
#   posterior draws, giving the weights of these draws. These weights should be
#   treated as not being normalized (i.e., they don't necessarily sum to `1`).
#   Currently, this element could be named `wdraws_ref` instead because
#   get_p_clust() is always applied to inputs that are specific to a `refmodel`
#   object (either the initial reference model or a K-fold-specific `refmodel`
#   object) (and get_refdist() is applied to inputs that are specific to a
#   `refmodel` object anyway). However, get_p_clust() intentionally seems to
#   have been kept as general as possible and `wdraws_orig` is more general than
#   `wdraws_ref`.
#   * `clust_used`: A single logical value indicating whether clustering (i.e.,
#   get_p_clust(), with the possibility that its output element `clust_used` is
#   overwritten by its argument `clust_used_forced`) has been used.
get_refdist <- function(refmodel, ndraws = NULL, nclusters = NULL,
                        thinning = TRUE,
                        throw_mssg_ndraws = getOption("projpred.mssg_ndraws",
                                                      TRUE)) {
  # Number of draws in the reference model:
  S <- NCOL(refmodel$mu)

  if (!is.null(nclusters)) {
    # use clustering (ignore ndraws argument)
    nclusters <- min(S, nclusters)
    if (nclusters == S) {
      # number of clusters equal to the number of draws, so return the draws
      return(get_refdist(refmodel, ndraws = nclusters,
                         throw_mssg_ndraws = FALSE))
    } else if (nclusters == 1) {
      # special case, only one cluster
      p_ref <- get_p_clust(family = refmodel$family, eta = refmodel$eta,
                           mu = refmodel$mu, mu_offs = refmodel$mu_offs,
                           dis = refmodel$dis, wobs = refmodel$wobs,
                           cl = rep(1, S))
    } else {
      # several clusters
      p_ref <- get_p_clust(family = refmodel$family, eta = refmodel$eta,
                           mu = refmodel$mu, mu_offs = refmodel$mu_offs,
                           dis = refmodel$dis, wobs = refmodel$wobs,
                           nclusters = nclusters)
    }
  } else {
    if (length(unique(refmodel$wdraws_ref)) != 1) {
      stop("Currently, projpred requires the reference model's posterior ",
           "draws to have constant weights.")
    }
    ndraws <- min(S, ndraws)
    if (ndraws <= 20 && throw_mssg_ndraws) {
      message("The number of draws to project is quite small (<= 20). In such ",
              "cases, it is usually better to use clustering.")
    }
    if (thinning) {
      s_ind <- round(seq(from = 1, to = S, length.out = ndraws))
    } else {
      s_ind <- draws_subsample(S = S, ndraws = ndraws)
    }
    cl <- rep(NA, S)
    cl[s_ind] <- 1:ndraws
    predvar <- do.call(cbind, lapply(s_ind, function(j) {
      refmodel$family$predvar(refmodel$mu_offs[, j, drop = FALSE],
                              refmodel$dis[j], refmodel$wdraws_ref[j])
    }))
    p_ref <- list(
      mu = refmodel$mu[, s_ind, drop = FALSE],
      mu_offs = refmodel$mu_offs[, s_ind, drop = FALSE],
      var = structure(predvar,
                      ndiscrete = attr(refmodel$mu, "ndiscrete"),
                      class = oldClass(refmodel$mu)),
      dis = refmodel$dis[s_ind],
      wdraws_prj = rep(1 / ndraws, ndraws),
      const_wdraws_prj = TRUE,
      nprjdraws = ndraws,
      cl = cl,
      wdraws_orig = rep(1, S),
      clust_used = FALSE
    )
  }

  return(p_ref)
}

# Returns the draws (columns) of `eta` as an S x r matrix of coordinates in an
# orthonormal basis of their centered span. These coordinates give the same
# Euclidean distances as `t(eta)`, so kmeans() finds the same clustering (up to
# rounding). The cost of kmeans() then does not grow with the number of
# observations. The draws of a GLM-type reference model span at most as many
# dimensions as the model has coefficients. The basis comes from a
# deterministic subset of the draws, so the state of the random number
# generator does not change. The function uses the basis only if it reproduces
# all draws (relative tolerance `100 * tol`; `tol` is the rank tolerance of
# qr()).
# Otherwise, the function returns `t(eta)`.
eta_for_kmeans <- function(eta, rank_max = floor(min(dim(eta)) / 4),
                           tol = 1e-10) {
  S <- ncol(eta)
  if (rank_max < 2 || anyNA(eta)) {
    return(t(eta))
  }
  eta_c <- eta - rowMeans(eta)
  scale_c <- max(abs(eta_c))
  if (scale_c == 0) {
    return(t(eta))
  }
  k <- min(64L, rank_max)
  repeat {
    idx <- unique(round(seq(1, S, length.out = k)))
    qr_sub <- qr(eta_c[, idx, drop = FALSE], tol = tol)
    if (qr_sub$rank < length(idx)) {
      # The subset has a lower rank than its number of draws, so its span can
      # be the span of all draws (checked next):
      q <- qr.Q(qr_sub)[, seq_len(qr_sub$rank), drop = FALSE]
      coords <- crossprod(q, eta_c)
      if (max(abs(eta_c - q %*% coords)) <= tol * scale_c * 1e2) {
        return(t(coords))
      }
      return(t(eta))
    }
    if (k >= rank_max) {
      return(t(eta))
    }
    k <- min(2L * k, rank_max)
  }
}

# Function for clustering the parameter draws:
get_p_clust <- function(family, eta, mu, mu_offs, dis, nclusters = 10,
                        wobs = rep(1, dim(mu)[1]),
                        wdraws = rep(1, dim(mu)[2]), cl = NULL,
                        clust_used_forced = NULL, as_moments = FALSE) {
  # With `as_moments = TRUE`, the function returns `NULL` or, for a given `cl`
  # with one draw per cluster, finite nonzero weights and a
  # `family$predvar_mix`, a list without the observations x clusters matrices
  # `mu`, `mu_offs` and `var`. Instead, the list has `draws_offs` (the input
  # `mu_offs`), the draw `idx` and the factor `ws` of each cluster (cluster
  # center `ws[j] * mu_offs[, idx[j]]`), and `draws_predvar` (the predictive
  # variance of cluster j for all observations). The function also returns
  # the list for non-finite projected draws, where `draws_predvar` can be
  # wrong (see `mu_var` below). perf_eval_gauss() detects them with the
  # weighted sums of squares and then uses the full result.
  # (The names do not start with "mu" or "var", because `$` matches partial
  # names.)
  if (as_moments && is.null(cl)) {
    # Do not run kmeans() twice (here and in the call without `as_moments`):
    return(NULL)
  }
  # cluster the draws in the latent space if no clustering provided
  if (is.null(cl)) {
    # Note: A seed is not set here because this function is not exported and has
    # a calling stack at the beginning of which a seed is set.

    out <- kmeans(eta_for_kmeans(eta), nclusters, iter.max = 50)
    cl <- out$cluster # cluster indices for each draw
  } else if (typeof(cl) == "list") {
    # old clustering solution provided, so fetch the cluster indices
    if (is.null(cl$cluster)) {
      stop("argument cl must be a vector of cluster indices or a clustering ",
           "object returned by k-means.")
    }
    cl <- cl$cluster
  }

  # (Re)compute the cluster centers, because they may be different from the ones
  # returned by kmeans() if the draws have differing weights. `cl` has one
  # element per draw (column of `mu`).
  # Number of clusters (assumes labeling "1, ..., nclusters"):
  nclusters <- max(cl, na.rm = TRUE)
  # Draw indices of each cluster (one pass over `cl`; draws with `cl` `NA` are
  # in no cluster). match() compares numbers, as `cl == j` does:
  cl_ind <- match(cl, seq_len(nclusters))
  ind_by_cl <- split(seq_along(cl),
                     factor(cl_ind, levels = seq_len(nclusters)))
  eps <- 1e-10
  nobs <- dim(mu)[1]
  # Use `family$predvar_mix` only if `family$predvar` is the function that it
  # belongs to (a user can replace `family$predvar`):
  use_mix <- !is.null(family$predvar_mix) &&
    identical(family$predvar, attr(family$predvar_mix, "predvar"))
  # The weighted mean of `dis^2` and the weighted variance of `mu_offs` within
  # the clusters are needed only if `family$predvar_mix` uses them:
  need_var <- use_mix && !isFALSE(attr(family$predvar_mix, "uses_moments"))
  if (as_moments && !all(lengths(ind_by_cl) == 1L)) {
    return(NULL)
  }
  if (all(lengths(ind_by_cl) == 1L)) {
    # One draw per cluster (thinned draws, or as many clusters as draws). This
    # gives the same results as the cluster loop at the end, for all clusters
    # at once.
    idx <- unlist(ind_by_cl, use.names = FALSE)
    wcluster <- unname(as.numeric(wdraws[idx]))
    ws <- wcluster / wcluster * (1 - eps)
    # Usually, all draws are used in their order, and `w / w` is exactly 1
    # (for a finite nonzero weight), so that `ws` is `1 - eps`. Then one
    # product with a scalar gives the centers, without copies of the draws:
    all_draws <- identical(idx, seq_len(ncol(mu)))
    if (as_moments) {
      if (!use_mix || anyNA(ws)) {
        return(NULL)
      }
      return(list(
        moments_only = TRUE, draws_offs = mu_offs, idx = idx, ws = ws,
        draws_predvar = drop(family$predvar_mix(as.numeric(dis[idx])^2, 0,
                                                1L)),
        dis = unname(as.numeric(dis[idx]) * ws),
        wdraws_prj = wcluster / sum(wcluster),
        const_wdraws_prj = length(unique(wcluster / sum(wcluster))) == 1,
        nprjdraws = nclusters,
        cl = cl, wdraws_orig = wdraws,
        clust_used = clust_used_forced %||% TRUE
      ))
    }
    ws_mult <- if (!anyNA(ws)) 1 - eps else rep(ws, each = nobs)
    scale_draws <- function(m) {
      if (!all_draws) {
        m <- unclass(m)[, idx, drop = FALSE]
      }
      out <- m * ws_mult
      attributes(out) <- list(dim = dim(out))
      out
    }
    centers <- scale_draws(mu)
    centers_offs <- scale_draws(mu_offs)
    dis_agg <- unname(as.numeric(dis[idx]) * ws)
    if (use_mix) {
      # The variance across one draw is 0 (NaN for a non-finite `mu_offs`, as
      # in rowvar_wdraws()):
      mu_var <- 0
      if (need_var) {
        mu_offs_idx <- unclass(mu_offs)
        if (!all_draws) {
          mu_offs_idx <- mu_offs_idx[, idx, drop = FALSE]
        }
        if (length(mu_offs_idx) > 0 && !all(is.finite(range(mu_offs_idx)))) {
          mu_var <- (mu_offs_idx - mu_offs_idx)^2
          attributes(mu_var) <- list(dim = dim(mu_var))
        }
      }
      predvar <- family$predvar_mix(as.numeric(dis[idx])^2, mu_var, nobs)
    }
  } else {
    # Several draws in a cluster: weighted moments of the draws in compiled
    # code. The results differ from the cluster loop only by rounding.
    wcluster <- vapply(ind_by_cl, function(ind) sum(wdraws[ind]), numeric(1),
                       USE.NAMES = FALSE)
    wn <- as.numeric(wdraws) / wcluster[cl_ind]
    wn[is.na(cl_ind)] <- 0
    mom <- clust_moments_c(
      mu, mu_offs, as.numeric(dis)[seq_along(cl)], wn,
      ifelse(is.na(cl_ind), -1L, cl_ind - 1L), nclusters, need_var
    )
    centers <- mom$mu * (1 - eps)
    centers_offs <- mom$mu_offs * (1 - eps)
    dis_agg <- mom$dis * (1 - eps)
    if (use_mix) {
      predvar <- family$predvar_mix(mom$dis2, mom$var_offs, nobs)
    }
  }
  # The cluster loop: for the predictive variances without
  # `family$predvar_mix`, and for the clusters with a weight sum of 0, Inf or
  # NaN, or without draws (so that NA and NaN propagate as in the loop):
  ind_loop <- which(!is.finite(wcluster) | wcluster == 0)
  if (!use_mix) {
    predvar <- matrix(nrow = nobs, ncol = nclusters)
    ind_loop <- seq_len(nclusters)
  }
  for (j in ind_loop) {
    ind <- ind_by_cl[[j]]
    # Compute normalized weights within the j-th cluster; `1 - eps` is for
    # numerical stability:
    ws_j <- wdraws[ind] / sum(wdraws[ind]) * (1 - eps)
    # Predictive variance for the j-th cluster:
    predvar[, j] <- family$predvar(mu_offs[, ind, drop = FALSE], dis[ind],
                                   ws_j)
    if (is.finite(wcluster[j]) && wcluster[j] != 0) {
      next
    }
    # Center of the j-th cluster:
    centers[, j] <- mu[, ind, drop = FALSE] %*% ws_j
    # The same centers, but taking offsets into account:
    centers_offs[, j] <- mu_offs[, ind, drop = FALSE] %*% ws_j
    # Aggregated dispersion parameter for the j-th cluster:
    dis_agg[j] <- crossprod(dis[ind], ws_j)
  }
  predvar <- unname(predvar)
  wcluster <- wcluster / sum(wcluster)

  # combine the results
  return(list(
    mu = structure(unname(centers),
                   ndiscrete = attr(mu, "ndiscrete"),
                   class = oldClass(mu)),
    mu_offs = structure(unname(centers_offs),
                        ndiscrete = attr(mu_offs, "ndiscrete"),
                        class = oldClass(mu_offs)),
    var = structure(predvar,
                    ndiscrete = attr(mu, "ndiscrete"),
                    class = oldClass(mu)),
    dis = dis_agg,
    wdraws_prj = wcluster,
    const_wdraws_prj = length(unique(wcluster)) == 1,
    nprjdraws = nclusters,
    cl = cl,
    wdraws_orig = wdraws,
    clust_used = clust_used_forced %||% TRUE
  ))
}

draws_subsample <- function(S, ndraws) {
  # Note: A seed is not set here because this function is not exported and has a
  # calling stack at the beginning of which a seed is set.

  return(sample.int(S, size = ndraws))
}

# If no get_refdist() (or get_p_clust()) output is available to infer elements
# `clust_used` (a single logical value indicating whether the given combination
# of `ndraws` and `nclusters` will lead to clustered projected draws or not) and
# `nprjdraws` (a single numeric value giving the number of (possibly clustered)
# projected draws resulting from the given combination of `ndraws` and
# `nclusters`), the following function can be used as a workaround (argument `S`
# denotes the number of posterior draws in the reference model object):
clust_info <- function(ndraws, nclusters, S) {
  clust_used_out <- !is.null(nclusters) && nclusters < S
  nprjdraws_out <- if (clust_used_out) nclusters else ndraws
  return(list(
    clust_used = clust_used_out,
    nprjdraws = nprjdraws_out
  ))
}

is_proj_list <- function(proj) {
  # Better use a formal class `proj_list`, but for now, use this workaround:
  is.list(proj) && length(proj) && all(sapply(proj, inherits, "projection"))
}

unlist_proj <- function(p) {
  if (length(p) == 1) p[[1]] else p
}

## create a named list using object names
nlist <- function(...) {
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE else nzchar(names(dots))
  if (all(has_name)) {
    return(dots)
  }
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

# The `%||%` special binary (infix) operator from brms (equivalent to the
# `%ORifNULL%` operator from rstanarm):
`%||%` <- function(x, y) {
  if (is.null(x)) x <- y
  x
}

#' Execute a function call
#'
#' Execute a function call similar to [do.call()], but without deparsing
#' function arguments.
#'
#' @param what Either a function or a non-empty character string naming the
#'   function to be called.
#' @param args A `list` of arguments to the function call. The [`names`]
#'   attribute of `args` gives the argument names.
#' @param pkg Optional name of the package in which to search for the function
#'   if `what` is a character string.
#'
#' @return The result of the (evaluated) function call.
#'
#' @keywords internal
#' @export
do_call <- function(what, args, pkg = NULL) {
  call <- ""
  if (length(args)) {
    if (!is.list(args)) {
      stop2("'args' must be a list.")
    }
    fun_args <- names(args)
    if (is.null(fun_args)) {
      fun_args <- rep("", length(args))
    } else {
      nzc <- nzchar(fun_args)
      fun_args[nzc] <- paste0("`", fun_args[nzc], "` = ")
    }
    names(args) <- paste0(".x", seq_along(args))
    call <- paste0(fun_args, names(args), collapse = ",")
  } else {
    args <- list()
  }
  if (is.function(what)) {
    args$.fun <- what
    what <- ".fun"
  } else {
    what <- paste0("`", as_one_character(what), "`")
    if (!is.null(pkg)) {
      what <- paste0(as_one_character(pkg), "::", what)
    }
  }
  call <- paste0(what, "(", call, ")")
  eval2(call, envir = args, enclos = parent.frame())
}

# like 'eval' but parses characters before evaluation
eval2 <- function(expr, envir = parent.frame(), ...) {
  if (is.character(expr)) {
    expr <- parse(text = expr)
  }
  eval(expr, envir, ...)
}

# coerce `x` to a single character string
as_one_character <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.character(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- deparse_combine(s, max_char = 100L)
    stop2("Cannot coerce '", s, "' to a single character value.")
  }
  x
}

stop2 <- function(...) {
  stop(..., call. = FALSE)
}

# combine deparse lines into one string
deparse_combine <- function(x, max_char = NULL) {
  out <- paste0(deparse(x), collapse = "")
  if (isTRUE(max_char > 0)) {
    out <- substr(out, 1L, max_char)
  }
  out
}

# `R CMD check` throws a note when using <package>:::<function>() (for accessing
# <function> which is not exported by its <package>). Of course, usage of
# non-exported functions should be avoided, but sometimes there's no way around
# that. Thus, with the following helper operator, it is possible to redefine
# such functions here in projpred:
`%:::%` <- function(pkg, fun) {
  # Note: `utils::getFromNamespace(fun, pkg)` could probably be used, too (but
  # its documentation is unclear about the inheritance from parent
  # environments).
  get(fun, envir = asNamespace(pkg), inherits = FALSE)
}

# Helper function to combine separate `list`s into a single `list`:
rbind2list <- function(x) {
  is_augeval <- any(sapply(x, function(x_i) {
    is.list(x_i) &&
      identical(names(x_i), c("mu", "lppd")) &&
      inherits(x_i$mu, "augvec")
  }))
  if (is_augeval) {
    mu_arr <- abind::abind(
      lapply(x, function(x_i) {
        augmat2arr(augvec2augmat(x_i$mu))
      }),
      along = 1
    )
    return(c(
      list(mu = augmat2augvec(arr2augmat(mu_arr))),
      rbind2list(lapply(x, "[", "lppd"))
    ))
  }
  binded_list <- as.list(do.call(rbind, lapply(x, function(x_i) {
    as.data.frame(x_i[setdiff(names(x_i), "oscale")])
  })))
  is_lateval_oscale <- any(sapply(x, function(x_i) {
    is.list(x_i) &&
      identical(names(x_i), c("mu", "lppd", "oscale"))
  }))
  if (is_lateval_oscale) {
    binded_list$oscale <- rbind2list(lapply(x, "[[", "oscale"))
  }
  return(binded_list)
}

cat_cls <- function(x) {
  cls <- paste0("`", class(x), "`")
  cat("Object of class", if (length(cls) > 1) "es" else "", " ",
      paste(cls, collapse = ", "), "\n\n", sep = "")
}

# Print out text via message() if `verbose = TRUE`:
verb_out <- function(..., verbose = TRUE) {
  if (verbose) {
    message(...)
  }
}

# Helper function intended for use in verbose messages:
txt_clust_draws <- function(clust_used, nprjdraws) {
  out <- paste0(nprjdraws, " ")
  if (clust_used) {
    out <- paste0(out, "cluster")
  } else {
    out <- paste0(out, "draw")
  }
  if (nprjdraws > 1) {
    out <- paste0(out, "s")
  }
  return(out)
}

# Ensure that stderr() is used for throwing an error, even while sink()-ing or
# capture.output()-ing with `type = "message"`:
throw_err <- function(e) {
  sink(type = "message")
  stop(e)
}

# A wrapper for capture.output() with `type = "message"`, but throwing error
# messages appropriately (i.e., only messages and warnings are captured).
# Note: This function should only be used to filter out messages or warnings
# (not to make downstream code dependent on catched messages or warnings), see
# <https://github.com/stan-dev/loo/issues/227#issuecomment-1663499985>.
capt_mssgs_warns <- function(expr) {
  if (getOption("warn") == 0) {
    warn_orig <- options(warn = 1)
    on.exit(options(warn_orig))
  }
  mssgs_warns_capt <- utils::capture.output(tryCatch(expr, error = throw_err),
                                            type = "message")
  if (identical(Sys.getenv("RSTUDIO"), "1")) {
    # Some RStudio versions (observed in RStudio 2025.05.0+496) add special
    # characters/symbols, so remove them:
    mssgs_warns_capt <- gsub("\033[GHgh]{1}2?;?", "", mssgs_warns_capt)
  }
  mssgs_warns_capt <- setdiff(mssgs_warns_capt, "")
  return(mssgs_warns_capt)
}

# Cache of the parts of a projection that depend only on the reference model and
# on the predictor terms of the submodel (see proj_to_submodl()). varsel(),
# cv_varsel() and project() activate it for the duration of their call; nested
# calls share it (counted by `depth`), and the outermost call empties it at its
# exit. Each reference model gets its own store, identified by the environment
# of its `fetch_data()` function. Global option `projpred.proj_cache` (default
# `TRUE`) switches the cache off.
#
# The size of the cache is limited: a forward search with a validated search
# adds an entry for each candidate submodel of each fold, and most of these
# entries are not used again. Each store has two generations of entries
# (`cur` and `old`). When the entries added to all `cur` generations exceed
# global option `projpred.proj_cache_mb` (default 512, in units of 2^20 bytes
# as estimated by utils::object.size(); read at the outermost call), every
# `cur` generation becomes the `old` one and the previous `old` one is dropped.
# A hit in `old` moves the entry to `cur`, so the entries that each fold uses
# (the small submodels, for example) stay in the cache. Thus, the cache holds
# at most about twice the larger of this limit and the largest entry. The limit
# applies to each R process; forked parallel workers have their own copy of the
# cache.
proj_cache <- new.env(parent = emptyenv())
proj_cache$depth <- 0L
proj_cache$anchors <- list()
proj_cache$stores <- list()
proj_cache$bytes <- 0
proj_cache$limit <- 512 * 2^20

proj_cache_enter <- function() {
  if (proj_cache$depth == 0L && getOption("projpred.proj_cache", TRUE)) {
    limit_mb <- getOption("projpred.proj_cache_mb", 512)
    if (!is.numeric(limit_mb) || length(limit_mb) != 1 || is.na(limit_mb) ||
        limit_mb < 0) {
      stop("Global option `projpred.proj_cache_mb` has to be a single ",
           "non-negative number.")
    }
    proj_cache$limit <- limit_mb * 2^20
  }
  proj_cache$depth <- proj_cache$depth + 1L
  return(invisible())
}

proj_cache_exit <- function() {
  proj_cache$depth <- max(proj_cache$depth - 1L, 0L)
  if (proj_cache$depth == 0L) {
    proj_cache$anchors <- list()
    proj_cache$stores <- list()
    proj_cache$bytes <- 0
  }
  return(invisible())
}

# Whether `key` can be the name of an entry. `[[` and assign() make a symbol
# from the key, and a symbol has a limit of 10000 bytes, so long keys are not
# cached:
proj_cache_key_ok <- function(key) {
  nchar(key, type = "bytes") < 10000L
}

# The entry `key` of store `store`, or `NULL` if there is none:
proj_cache_get <- function(store, key) {
  if (!proj_cache_key_ok(key)) {
    return(NULL)
  }
  ent <- store$cur[[key]]
  if (!is.null(ent)) {
    return(ent$value)
  }
  ent <- store$old[[key]]
  if (is.null(ent)) {
    return(NULL)
  }
  rm(list = key, envir = store$old)
  proj_cache_add(store, key, ent)
  return(ent$value)
}

# Stores `value` as entry `key` of store `store`:
proj_cache_put <- function(store, key, value) {
  if (!proj_cache_key_ok(key)) {
    return(invisible())
  }
  proj_cache_add(store, key, list(value = value,
                                  bytes = as.numeric(utils::object.size(value))))
  return(invisible())
}

# Adds the entry `ent` (with elements `value` and `bytes`) to the `cur`
# generation of store `store`, after the change of generations if the limit
# would be exceeded:
proj_cache_add <- function(store, key, ent) {
  if (proj_cache$bytes + ent$bytes > proj_cache$limit) {
    for (st in proj_cache$stores) {
      st$old <- st$cur
      st$cur <- new.env(parent = emptyenv())
    }
    proj_cache$bytes <- 0
  }
  assign(key, ent, envir = store$cur)
  proj_cache$bytes <- proj_cache$bytes + ent$bytes
  return(invisible())
}

# The store of reference model `refmodel`, or `NULL` if the cache is not active.
# A store is an environment with the generations `cur` and `old` (see above):
proj_cache_store <- function(refmodel) {
  if (proj_cache$depth == 0L || !getOption("projpred.proj_cache", TRUE)) {
    return(NULL)
  }
  anchor <- environment(refmodel$fetch_data)
  if (!is.environment(anchor)) {
    return(NULL)
  }
  for (i in seq_along(proj_cache$anchors)) {
    if (identical(proj_cache$anchors[[i]], anchor)) {
      return(proj_cache$stores[[i]])
    }
  }
  store <- new.env(parent = emptyenv())
  store$cur <- new.env(parent = emptyenv())
  store$old <- new.env(parent = emptyenv())
  proj_cache$anchors <- c(proj_cache$anchors, list(anchor))
  proj_cache$stores <- c(proj_cache$stores, list(store))
  return(store)
}

# Like capt_mssgs_warns(), but collects the messages and warnings by condition
# handlers. Faster, but does not capture text that is written directly to the
# message stream, so use it only for fitters that throw R conditions.
capt_conds <- function(expr) {
  capt <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      capt <<- c(capt, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      capt <<- c(capt, sub("\n$", "", conditionMessage(m)))
      invokeRestart("muffleMessage")
    }
  )
  capt <- as.character(unlist(strsplit(capt, "\n", fixed = TRUE)))
  return(setdiff(capt, ""))
}

# Parse the argument containing the observation weights (`wobs` or `weights`)
# for the <family_object>$ppd() functions used by proj_predict():
parse_wobs_ppd <- function(wobs, n_obs) {
  if (length(wobs) == 0) {
    wobs <- rep(1, n_obs)
  } else if (length(wobs) == 1) {
    wobs <- rep(wobs, n_obs)
  } else if (length(wobs) != n_obs) {
    stop("Argument `wobs` needs to be of length 0, 1, or the number of ",
         "observations.")
  }
  if (!all(wobs == 1) && getOption("projpred.warn_wobs_ppd", TRUE)) {
    warning("Currently, proj_predict() ignores observation weights not equal ",
            "to `1`.")
  }
  return(wobs)
}

# Constructs all possible permutations of a single interaction term `ia` (given
# as a single character string):
all_ia_perms <- function(ia, is_split = FALSE) {
  if (is_split) {
    ia_split <- ia
  } else {
    ia_split <- strsplit(ia, ":")[[1]]
  }
  ia_perms_split <- gtools::permutations(n = length(ia_split),
                                         r = length(ia_split),
                                         v = ia_split, set = FALSE)
  return(unlist(
    apply(ia_perms_split, 1, paste, collapse = ":", simplify = FALSE)
  ))
}

# Reorders the main-effect terms involved in a single interaction term `ia`
# (given as a single character string) such that it matches a corresponding
# interaction term in a character vector `y` (e.g., `a:b` is considered to be
# the same as `b:a`):
reorder_ia <- function(ia, y) {
  ia_perms <- all_ia_perms(ia)
  ia_reordered <- intersect(ia_perms, y)
  if (length(ia_reordered) == 0) {
    ia_reordered <- NA_character_
  } else if (length(ia_reordered) > 1) {
    stop("Unexpected length of `ia_reordered`. Please notify the package ",
         "maintainer.")
  }
  return(ia_reordered)
}

# For each interaction term in a character vector `x`, this function reorders
# the main-effect terms involved in it such that the reordered interaction term
# matches a corresponding interaction term in a character vector `y` (e.g.,
# `a:b` is considered to be the same as `b:a`):
reorder_ias <- function(x, y) {
  ia_idxs <- grep(":", x)
  ia_idxs <- ia_idxs[!x[ia_idxs] %in% y]
  for (ia_idx in ia_idxs) {
    x[ia_idx] <- reorder_ia(x[ia_idx], y)
  }
  return(x)
}

# Retrieves an element (with name given in argument `nm`) that is duplicated
# across the elements of a `list`, typically a `list` of `submodl`s:
element_unq <- function(list_obj, nm) {
  if (getOption("projpred.additional_checks", FALSE)) {
    el_unq <- unique(unlist(lapply(list_obj, "[[", nm)))
    stopifnot(length(el_unq) == 1)
  } else {
    el_unq <- list_obj[[1]][[nm]]
  }
  return(el_unq)
}

# Returns `fun` with a new enclosing environment that holds only the objects
# `nms` of `env` (its parent is `topenv(env)`, here the projpred namespace). A
# parallel backend then serializes these objects together with `fun`, not the
# whole `env`. `nms` must list every free variable of `fun` that is a local
# of `env`; a missing name is looked up in the namespace instead (test
# "the object lists of the loop functions ..." in test_parallel.R).
fun_with_objs <- function(fun, nms, env = parent.frame()) {
  environment(fun) <- list2env(mget(nms, envir = env, inherits = FALSE),
                               parent = topenv(env))
  return(fun)
}

# The objects of the init_refmodel() environment that mu_fun() and
# fetch_data() of a reference model use (see refmodel_for_folds()):
refmodel_fun_objs <- c("data", "family", "proj_predfun",
                       "proj_predfun_is_subprd")

# A copy of `refmodel` for the parallel CV folds (one_obs() in loo_varsel();
# one_fold() in kfold_varsel(), for the full-data reference model and for the
# reference models of the folds). The folds use the data, the draws, `family`,
# `div_minimizer`, mu_fun() and fetch_data() of the reference model, but not
# the reference model fit, its fold fits (`cvfits`), and the functions that
# need the fit (ref_predfun(), extract_model_data(), cvfun(), cvrefbuilder();
# kfold_varsel() calls ref_predfun() before the parallel loop). Without the
# copy, the workers get the fit and load the packages that it needs (e.g.,
# brms or rstanarm), which takes seconds and memory. The copy has `NULL` for the fits and
# not_in_folds() for these functions. mu_fun() and fetch_data() of
# init_refmodel() share an environment that also holds the fit, so they get
# an environment with only the objects that they use. Other mu_fun() and
# fetch_data() functions stay as they are. With `formula_env`, the formula of
# the copy gets a small environment (see formula_for_folds(); for the
# reference model of a K-fold fold).
refmodel_for_folds <- function(refmodel, formula_env = NULL) {
  # The family functions of the latent and the augmented-data projection, and
  # family functions from other packages, can hold the fit in their
  # environments. The copy cannot remove the fit from there, so the workers get
  # the reference model as it is, also its formula. (Without `formula_env`,
  # the environment of the formula stays as it is. Functions of the user, e.g.,
  # a custom `proj_predfun`, also stay as they are. They can also hold the
  # fit.)
  fam_pkgs <- vapply(Filter(is.function, refmodel$family), function(f) {
    environmentName(topenv(environment(f) %||% baseenv()))
  }, character(1))
  if (isTRUE(refmodel$family$for_latent) ||
      isTRUE(refmodel$family$for_augdat) ||
      !all(fam_pkgs %in% c("base", "stats", "projpred"))) {
    return(refmodel)
  }
  out <- refmodel
  out["fit"] <- list(NULL)
  out["cvfits"] <- list(NULL)
  for (nm in intersect(c("ref_predfun", "extract_model_data", "cvfun",
                         "cvrefbuilder"), names(out))) {
    out[[nm]] <- not_in_folds
  }
  env <- environment(refmodel$mu_fun)
  if (is.environment(env) &&
      identical(environment(refmodel$fetch_data), env) &&
      all(vapply(refmodel_fun_objs, exists, logical(1), envir = env,
                 inherits = FALSE))) {
    env_new <- list2env(mget(refmodel_fun_objs, envir = env, inherits = FALSE),
                        parent = topenv(env))
    # The data can have formula attributes (e.g., "terms" and "formula" for
    # rstanarm::stan_gamm4()) with the environment of the model fitting, which
    # can hold the fit. projpred does not use these attributes, so give them
    # the environment of the reference model formula:
    for (anm in names(attributes(env_new$data))) {
      if (inherits(attr(env_new$data, anm), "formula")) {
        environment(attr(env_new$data, anm)) <- environment(refmodel$formula)
      }
    }
    environment(out$mu_fun) <- env_new
    environment(out$fetch_data) <- env_new
  }
  if (!is.null(formula_env)) {
    out$formula <- formula_for_folds(out$formula,
                                     data_nms = names(out$fetch_data(obs = 1L)),
                                     env_parent = formula_env)
  }
  return(out)
}

# The formula of the reference model of a K-fold fold can have the environment
# of the refit call (e.g., of rstanarm::kfold(), if the call of the fit has
# the formula as an expression). That environment can hold the fits, and its
# ancestors can be the namespace of the package of the fit. This function
# returns `formula` with a new environment that holds only the objects that
# the formula uses from the old environment and its ancestors, up to the first
# namespace or the global environment (without the names of the data columns
# `data_nms`, which the data provides). The parent of the new environment is
# `env_parent` (the environment of the full-data formula). So a name that the
# old environment found only in a namespace among its ancestors (e.g., a
# function of the package of the fit that is not exported) is looked up from
# `env_parent` instead. A formula whose environment is a namespace, the
# global, the base or the empty environment, or `env_parent` stays as it is.
formula_for_folds <- function(formula, data_nms, env_parent) {
  is_top <- function(env) {
    isNamespace(env) || identical(env, globalenv()) ||
      identical(env, baseenv()) || identical(env, emptyenv())
  }
  env <- environment(formula)
  if (!is.environment(env) || is_top(env) || identical(env, env_parent)) {
    return(formula)
  }
  nms <- setdiff(all.names(formula, unique = TRUE), data_nms)
  objs <- list()
  while (!is_top(env) && !identical(env, env_parent)) {
    for (nm in setdiff(nms, names(objs))) {
      if (exists(nm, envir = env, inherits = FALSE)) {
        # A function argument without a value gives an error here; the
        # formula cannot use it:
        obj <- tryCatch(list(get(nm, envir = env, inherits = FALSE)),
                        error = function(e) NULL)
        if (!is.null(obj)) {
          objs[nm] <- obj
        }
      }
    }
    env <- parent.env(env)
  }
  environment(formula) <- list2env(objs, parent = env_parent)
  return(formula)
}

not_in_folds <- function(...) {
  stop("This function of the reference model is not available in the ",
       "parallel CV folds.")
}

use_progressr <- function() {
  getOption("projpred.use_progressr",
            requireNamespace("progressr", quietly = TRUE) &&
              interactive() &&
              identical(foreach::getDoParName(), "doFuture"))
}

sqrt_cut0 <- function(x) {
  if (!is.na(x) && sign(x) == -1) {
    if (abs(x) < 1e3 * sqrt(.Machine$double.eps)) {
      x <- 0
    } else {
      stop("Negative (and numerically non-zero) value used as input to ",
           "sqrt_cut0().")
    }
  }
  return(sqrt(x))
}

verbose_from_deprecated_options <- function(verbose, with_cv = FALSE,
                                            proj_only = FALSE) {
  if (!is.logical(verbose) || !is.null(getOption("projpred.verbose"))) {
    return(verbose)
  }
  if (!is.null(getOption("projpred.extra_verbose")) &&
      !is.null(getOption("projpred.verbose_project"))) {
    warning(
      "Global options `projpred.extra_verbose` and `projpred.verbose_project` ",
      "are deprecated. Please use argument ",
      "`verbose` instead. This argument can also be controlled via a ",
      "corresponding global option. Now trying to find an appropriate value ",
      "for argument `verbose`."
    )
    verbose <- as.logical(verbose) * (
      proj_only * 2L +
        (!proj_only) * (
          (getOption("projpred.extra_verbose") ||
             getOption("projpred.verbose_project")) *
            (1L + with_cv + getOption("projpred.verbose_project")) +
            1L
        )
    )
  } else {
    verbose <- verbose_from_deprecated_option(
      verbose,
      with_cv = with_cv,
      proj_only = proj_only,
      nm_option = "projpred.extra_verbose"
    )
    verbose <- verbose_from_deprecated_option(
      verbose,
      with_cv = with_cv,
      proj_only = proj_only,
      nm_option = "projpred.verbose_project"
    )
  }
  return(verbose)
}

verbose_from_deprecated_option <- function(verbose, with_cv = FALSE,
                                           proj_only = FALSE, nm_option) {
  if (!is.null(getOption(nm_option))) {
    warning(
      "Global option `", nm_option, "` is deprecated. Please use argument ",
      "`verbose` instead. This argument can also be controlled via a ",
      "corresponding global option. Now trying to find an appropriate value ",
      "for argument `verbose`."
    )
    verbose <- as.logical(verbose) * (
      proj_only * 2L +
        (!proj_only) * (
          getOption(nm_option) *
            (1L + with_cv + identical(nm_option, "projpred.verbose_project")) +
            1L
        )
    )
  }
  return(verbose)
}
