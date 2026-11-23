# Projection onto a negative binomial (NB) submodel with mean `mu` and shape
# `phi` (variance mu + mu^2 / phi).
#
# Notation: for one projected draw (or cluster), the reference model gives the
# NB distribution of observation i with mean `mu_r[i]` and shape `shape_r[i]`.
# The submodel gives the mean `mu[i]` and one shape `phi`. With observation
# weights `w`, the projection maximizes
#   L(mu, phi) = sum_i w[i] E_r[log NB(y_i | mu[i], phi)].
# Per observation, the expectation is
#   -phi log1p(mu / phi) - mu_r log1p(phi / mu)
#   + E_r[lgamma(y + phi) - lgamma(phi) - lgamma(y + 1)].
# For a fixed `phi`, only `mu_r` enters (so the coefficients are the NB fit to
# the response `mu_r`). The last term enters through the weighted sum over the
# observations, which depends on the reference only through the aggregated
# probability mass function P(y) = sum_i w[i] p_i(y). nb_gsum() computes P
# exactly for y < `y0` and, for larger y, as bins with their exact mass, mean
# and variance; nb_geval() evaluates the sum and its derivatives in `phi`.

# The probability of the upper tail beyond which nb_gsum() ignores the mass:
nb_tail_prob <- function() 1e-17

# The method of nb_gsum() from global option `projpred.nb_gsum`: "quadrature"
# (the default; exact for y < y0, a 4-point quadrature for the bins) or
# "exact" (the probabilities of all y in the support).
nb_gsum_method <- function() {
  method <- getOption("projpred.nb_gsum", "quadrature")
  if (!isTRUE(method %in% c("exact", "quadrature"))) {
    stop("Global option `projpred.nb_gsum` must be \"exact\" or ",
         "\"quadrature\".")
  }
  method
}

# Upper bound of the projected shape (the Poisson limit):
nb_phi_max <- function() 1e8

# The aggregated reference pmf for the observation weights `w` and the draw
# weights `v` (summing to 1): `p0` (P(y) for y = 0, ..., y0 - 1) and the bins
# for y >= y0 (columns `m` (mass), `mean`, `var`). `mu` and `shape` are
# observations x draws matrices (a vector is one draw; `shape` is recycled).
# The bins have the edges a_1 = y0, a_{j+1} = a_j + max(1, floor(c_bin *
# a_j^(3/4))); nb_gsum_c() computes their mass, mean and variance. The
# mass beyond the upper nb_tail_prob() quantile of each observation and draw
# is ignored. Only observations with a positive weight enter.
nb_gsum <- function(mu, shape, w, v = 1, y0 = 100L, c_bin = 0.25,
                    method = nb_gsum_method()) {
  mu <- as.matrix(mu)
  shape <- matrix(shape, nrow = nrow(mu), ncol = ncol(mu))
  if (anyNA(mu) || any(!is.finite(mu)) || any(mu < 0) || anyNA(shape) ||
      any(shape <= 0)) {
    stop("The negative binomial projection needs finite reference means ",
         ">= 0 and positive shapes (argument `dis`).", call. = FALSE)
  }
  w <- rep_len(as.numeric(w), nrow(mu))
  v <- rep_len(as.numeric(v), ncol(mu))
  out <- nb_gsum_c(mu, shape, w, v, as.integer(y0), c_bin, nb_tail_prob(),
                   match(method, c("exact", "quadrature")) - 1L)
  sums <- out$sums
  a_b <- out$edges[-length(out$edges)]
  b_b <- out$edges[-1] - 1
  # nb_gsum_c() gives the sums of the mass, of (y - a) and of (y - a)^2 for
  # bin [a, b]:
  pos <- sums[, 1] > 0
  d1 <- ifelse(pos, sums[, 2] / sums[, 1], 0)
  bmean <- pmin(pmax(a_b + d1, a_b), b_b)
  bvar <- ifelse(pos, pmin(pmax(sums[, 3] / sums[, 1] - d1^2, 0),
                           (b_b - a_b)^2 / 4), 0)
  bins <- cbind(m = sums[, 1], mean = bmean, var = bvar)[pos, , drop = FALSE]
  rownames(bins) <- NULL
  list(p0 = out$p0, bins = bins)
}

# For the aggregated pmf `gs` (from nb_gsum()), the sum over y of
# P(y) (lgamma(y + phi) - lgamma(phi) - lgamma(y + 1)) and its first and
# second derivatives in `phi` (scalar). A bin uses the second-order expansion
# of the summand around its mean.
nb_geval <- function(gs, phi, nderiv = 2L) {
  y <- seq_along(gs$p0) - 1
  m <- gs$bins[, "m"]
  x <- gs$bins[, "mean"]
  v <- gs$bins[, "var"]
  out <- sum(gs$p0 * (lgamma(y + phi) - lgamma(phi) - lgamma(y + 1))) +
    sum(m * (lgamma(x + phi) - lgamma(phi) - lgamma(x + 1) +
               (trigamma(x + phi) - trigamma(x + 1)) * v / 2))
  if (nderiv >= 1) {
    out <- c(out,
             sum(gs$p0 * (digamma(y + phi) - digamma(phi))) +
               sum(m * (digamma(x + phi) - digamma(phi) +
                          psigamma(x + phi, 2) * v / 2)))
  }
  if (nderiv >= 2) {
    out <- c(out,
             sum(gs$p0 * (trigamma(y + phi) - trigamma(phi))) +
               sum(m * (trigamma(x + phi) - trigamma(phi) +
                          psigamma(x + phi, 3) * v / 2)))
  }
  out
}

# The objective L (see the top of the file), the expected log density
# summed over the observations with weights `w`, and its first and second
# derivatives in `phi`:
nb_obj <- function(mu, mu_r, w, gs, phi, nderiv = 2L) {
  ge <- nb_geval(gs, phi, nderiv = nderiv)
  out <- sum(w * (-phi * log1p(mu / phi) - mu_r * log1p(phi / mu))) + ge[1]
  if (nderiv >= 1) {
    out <- c(out, sum(w * (-log1p(mu / phi) + (mu - mu_r) / (mu + phi))) +
               ge[2])
  }
  if (nderiv >= 2) {
    out <- c(out, sum(w * (1 / phi - 1 / (mu + phi) -
                             (mu - mu_r) / (mu + phi)^2)) + ge[3])
  }
  out
}

# The projected shape for the submodel means `mu`: maximizes nb_obj() in
# u = log(phi) with Newton steps (at most 2 in u), a gradient step where the
# objective is not concave, and step halving if the objective does not
# increase. The result is at most nb_phi_max().
nb_phi_step <- function(mu, mu_r, w, gs, phi0 = 1, tol = 1e-10,
                        maxit = 100L) {
  w <- rep_len(w, length(mu))
  keep <- w > 0
  mu <- mu[keep]
  mu_r <- mu_r[keep]
  w <- w[keep]
  u_max <- log(nb_phi_max())
  u <- min(log(if (is.finite(phi0) && phi0 > 0) phi0 else 1), u_max)
  cur <- nb_obj(mu, mu_r, w, gs, exp(u))
  for (it in seq_len(maxit)) {
    phi <- exp(u)
    d1 <- phi * cur[2]
    d2 <- phi * cur[2] + phi^2 * cur[3]
    step <- if (is.finite(d2) && d2 < 0) -d1 / d2 else sign(d1)
    step <- max(min(step, 2), -2)
    if (u >= u_max && step > 0) {
      break
    }
    for (k in 1:40) {
      u_new <- min(u + step, u_max)
      new <- nb_obj(mu, mu_r, w, gs, exp(u_new))
      if (is.finite(new[1]) && new[1] >= cur[1]) {
        break
      }
      step <- step / 2
    }
    if (!(is.finite(new[1]) && new[1] >= cur[1])) {
      break
    }
    done <- abs(u_new - u) < tol
    u <- u_new
    cur <- new
    if (done) {
      break
    }
  }
  exp(u)
}

# The names of the NB families of rstanarm (`neg_binomial_2`) and MASS
# (`negative.binomial(theta)`, family name "Negative Binomial(<theta>)") become
# "negbinomial" (the brms name):
nb_canonical_name <- function(family) {
  if (isTRUE(family$family == "neg_binomial_2") ||
      isTRUE(grepl("^Negative Binomial\\(", family$family))) {
    family$family <- "negbinomial"
  }
  family
}

# A family object for the fit of the coefficients with the shape `phi` fixed
# (from MASS::negative.binomial()). Its `aic` uses lgamma() instead of the
# factorial, so it is finite for a non-integer response.
nb_fixed_family <- function(phi, link) {
  fam <- MASS::negative.binomial(theta = phi, link = link)
  fam$aic <- function(y, n, mu, wt, dev) {
    term <- lgamma(phi + y) - lgamma(phi) - lgamma(y + 1) +
      phi * log(phi / (phi + mu)) + ifelse(y == 0, 0, y * log(mu / (phi + mu)))
    -2 * sum(term * wt)
  }
  environment(fam$aic) <- list2env(list(phi = phi), parent = baseenv())
  # pseudo_data() uses this flag. glm_ridge() calls extend_family(); the family
  # has all elements that glm_ridge() and pseudo_data() use, so it is marked as
  # extended:
  fam$nb_fixed <- TRUE
  fam$is_extended <- TRUE
  fam
}

extend_family_negbinomial <- function(family) {
  if (!isTRUE(family$link %in% c("log", "sqrt", "identity"))) {
    stop("For the negative binomial family, projpred supports only the links ",
         "\"log\", \"sqrt\" and \"identity\".", call. = FALSE)
  }
  link <- family$link
  # MASS::negative.binomial() gives the link functions for all three links
  # (a brms family object does not have them):
  fam_mass <- MASS::negative.binomial(theta = 1, link = link)
  family$linkfun <- fam_mass$linkfun
  family$linkinv <- fam_mass$linkinv
  family$mu.eta <- fam_mass$mu.eta
  family$valideta <- fam_mass$valideta

  ce_nb <- function(pref, data, psub) {
    wobs <- rep_len(data$weights, NROW(pref$mu))
    ce_sums <- vapply(seq_len(NCOL(pref$mu)), function(j) {
      gs <- nb_pref_gsum(pref, j, wobs)
      keep <- wobs > 0
      -nb_obj(psub$mu[keep, j], pref$mu[keep, j], wobs[keep], gs,
              psub$dis[j], nderiv = 0L)
    }, numeric(1))
    return(ce_sums / sum(wobs))
  }
  # The start of the shape step is `psub$dis` if present (the shapes of the
  # last iteration in nb_alternate()), else the reference shape:
  dis_nb <- function(pref, psub, wobs = 1) {
    wobs <- rep_len(wobs, NROW(pref$mu))
    vapply(seq_len(NCOL(pref$mu)), function(j) {
      gs <- nb_pref_gsum(pref, j, wobs)
      nb_phi_step(psub$mu[, j], pref$mu[, j], wobs, gs,
                  phi0 = psub$dis[j] %||% pref$dis[j] %||% 1)
    }, numeric(1))
  }
  predvar_na <- function(mu, dis, wdraws = 1) {
    rep(NA, NROW(mu))
  }
  ll_nb <- function(mu, dis, y, weights = 1, eta = NULL) {
    y <- as.matrix(y)
    dis <- matrix(rep(dis, each = NROW(mu)), ncol = NCOL(mu))
    weights * dnbinom(y, size = dis, mu = mu, log = TRUE)
  }
  # The deviance without the terms that do not depend on `mu` (so it is valid
  # for a non-integer `y`; it is for fits with a fixed `dis`). `dis` is
  # required:
  dev_nb <- function(mu, y, weights = 1, dis = NULL) {
    if (is.null(dis)) {
      stop("The deviance of the negative binomial family needs `dis`.")
    }
    dis <- matrix(rep(dis, each = NROW(mu)), ncol = NCOL(mu))
    if (NCOL(y) < NCOL(mu)) {
      y <- matrix(y, nrow = length(y), ncol = NCOL(mu))
    }
    2 * weights * (dis * log1p(mu / dis) + y * log1p(dis / mu))
  }
  ppd_nb <- function(mu, dis, weights = 1) {
    weights <- parse_wobs_ppd(weights, n_obs = length(mu))
    rnbinom(length(mu), size = dis, mu = mu)
  }

  family$ce <- ce_nb
  family$dis_fun <- dis_nb
  family$predvar <- predvar_na
  family <- add_predvar_mix(family, predvar_mix_na, uses_moments = FALSE)
  family$ll_fun <- ll_nb
  family$deviance <- dev_nb
  family$ppd <- ppd_nb
  family$fixed_dis_family <- function(phi) nb_fixed_family(phi, link)

  return(family)
}

# The aggregated reference pmf (nb_gsum()) of projected draw (or cluster) `j`
# of `pref` for the observation weights `wobs`. A cluster with several draws
# uses the mixture of the NB distributions of its draws: `pref$nb_mix` (see
# get_p_clust()) holds the draws (`mu`, observations x draws, with offsets),
# their shapes `dis`, the draw indices `ind` of each cluster and the
# normalized draw weights `v` of each cluster. Otherwise, column j of
# `pref$mu` is one draw with shape `pref$dis[j]`. `pref$nb_cache` (an
# environment, if present) keeps the last result of each column for its
# weights and means.
nb_pref_gsum <- function(pref, j, wobs) {
  wobs <- rep_len(as.numeric(wobs), NROW(pref$mu))
  mu_j <- pref$mu[, j]
  cache <- pref$nb_cache
  key <- as.character(j)
  method <- nb_gsum_method()
  if (is.environment(cache)) {
    hit <- cache[[key]]
    if (!is.null(hit) && identical(hit$wobs, wobs) &&
        identical(hit$mu_j, mu_j) && identical(hit$dis_j, pref$dis[j]) &&
        identical(hit$method, method)) {
      return(hit$gs)
    }
  }
  mix <- pref$nb_mix
  if (!is.null(mix) && length(mix$ind[[j]]) > 1L) {
    ind <- mix$ind[[j]]
    gs <- nb_gsum(mix$mu[, ind, drop = FALSE],
                  rep(mix$dis[ind], each = nrow(mix$mu)), wobs,
                  v = mix$v[[j]], method = method)
  } else {
    gs <- nb_gsum(mu_j, pref$dis[j], wobs, method = method)
  }
  if (is.environment(cache)) {
    assign(key, list(wobs = wobs, mu_j = mu_j, dis_j = pref$dis[j],
                     method = method, gs = gs),
           envir = cache)
  }
  gs
}
