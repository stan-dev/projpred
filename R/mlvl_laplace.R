# Laplace fit of multilevel submodels with random intercepts
#
# mlvl_laplace() fits binomial (logit) and Poisson (log) submodels whose
# group-level terms are all `(1 | <group>)`, with the objective of
# lme4::glmer() with `nAGQ = 1` (the Laplace deviance, see `src/mlvlfit.cpp`).
# It minimizes this objective with nlminb() and the exact gradient.  Unlike
# glmer() with lme4 (<= 2.1-0), it uses the likelihood of non-integer
# responses (projected means) without rounding.  The returned object has class
# `projpred_mlvl`.  divmin() calls it through fit_glmer_callback(); other
# multilevel submodels keep lme4::glmer().
#
# The stored `deviance` uses the deviance residuals.  For the Poisson family,
# lme4 uses the log-likelihood instead, so its value differs by a constant
# that depends only on the data.

# Family code of mlvl_laplace_c(), or `NULL` if mlvl_laplace() does not handle
# `family`:
mlvl_laplace_fam <- function(family) {
  if (identical(family$family, "binomial") && identical(family$link, "logit")) {
    return(0L)
  }
  if (identical(family$family, "poisson") && identical(family$link, "log")) {
    return(1L)
  }
  return(NULL)
}

# Names of the grouping variables of `formula`, or `NULL` if a group-level term
# is not of the form `(1 | <variable>)` or a grouping variable occurs twice:
mlvl_laplace_grps <- function(formula) {
  bars <- reformulas::findbars(formula)
  if (!length(bars)) return(NULL)
  grps <- vapply(bars, function(b) {
    if (!identical(b[[2]], 1) || !is.name(b[[3]])) return(NA_character_)
    as.character(b[[3]])
  }, character(1))
  if (anyNA(grps) || anyDuplicated(grps)) return(NULL)
  return(grps)
}

# Largest number of group levels (summed over the grouping variables) for the
# dense matrices of mlvl_laplace_c():
mlvl_laplace_max_levels <- 500L

# Design of mlvl_laplace(): the fixed-effects model matrix, the response, the
# weights and the level codes of the grouping variables, in the term order of
# lme4 (reformulas::mkReTrms() sorts the terms by decreasing number of levels if
# they are not in this order).  Returns `NULL` for cases that the fit does not
# handle (missing values, offsets, a rank-deficient model matrix, a response
# outside the range of the family, more than `mlvl_laplace_max_levels` group
# levels).
mlvl_laplace_design <- function(formula, data, weights, grps, fam) {
  if (!is.data.frame(data) || !all(grps %in% names(data))) return(NULL)
  n <- nrow(data)
  fixed <- reformulas::nobars(formula)
  mf <- tryCatch(
    model.frame(fixed, data = data, drop.unused.levels = TRUE,
                na.action = na.pass),
    error = function(e) NULL
  )
  if (is.null(mf) || nrow(mf) != n || !is.null(model.offset(mf))) return(NULL)
  trms <- attr(mf, "terms")
  if (attr(trms, "intercept") != 1) return(NULL)
  x <- model.matrix(trms, mf)
  y <- model.response(mf)
  # A factor response (the augmented-data projection) as in glm(): the first
  # level is failure, all others are success:
  if (fam == 0L && (is.factor(y) || is.logical(y))) {
    y <- if (is.factor(y)) as.numeric(y != levels(y)[1]) else as.numeric(y)
  }
  if (anyNA(x) || !is.numeric(y) || !is.null(dim(y)) || anyNA(y) ||
      qr(x)$rank < ncol(x)) {
    return(NULL)
  }
  if (any(y < 0) || (fam == 0L && any(y > 1))) return(NULL)
  if (is.null(weights)) weights <- rep(1, n)
  if (length(weights) != n || anyNA(weights) || any(weights < 0)) return(NULL)
  facs <- lapply(setNames(nm = grps), function(g) {
    droplevels(as.factor(data[[g]]))
  })
  if (any(vapply(facs, anyNA, logical(1)))) return(NULL)
  nlev <- vapply(facs, nlevels, integer(1))
  if (sum(nlev) > mlvl_laplace_max_levels) return(NULL)
  ord <- seq_along(nlev)
  if (any(diff(nlev) > 0)) ord <- rev(order(nlev))
  facs <- facs[ord]
  return(list(
    x = x, y = as.vector(y), w = as.vector(weights),
    grp = matrix(vapply(facs, as.integer, integer(n)) - 1L, nrow = n),
    nlev = nlev[ord], lvls = lapply(facs, levels),
    terms = delete.response(trms), xlevels = .getXlevels(trms, mf),
    contrasts = attr(x, "contrasts")
  ))
}

# Objective of the design `des` at (theta, beta), with the conditional modes
# found from u = 0:
mlvl_laplace_dev <- function(des, fam, theta, beta) {
  return(mlvl_laplace_c(theta, beta, des$x, des$grp, des$nlev, des$y, des$w,
                        rep(0, nrow(des$x)), rep(0, sum(des$nlev)), fam,
                        FALSE)$D)
}

# Minimize the Laplace deviance over (theta, beta) with nlminb().  Start: theta
# 1, beta from glm.fit().  The deviance is even in each theta_j, so theta_j = 0
# is a stationary point where the gradient-based search can stop: if a theta_j
# ends below `small_theta`, restart from `restart_theta` and keep the lower
# deviance.
mlvl_laplace_optim <- function(des, fam, rel_tol = 1e-10, small_theta = 1e-3,
                               restart_theta = 0.3) {
  k <- length(des$nlev)
  p <- ncol(des$x)
  n <- nrow(des$x)
  offs <- rep(0, n)
  # Start values from glm.fit(), without its aic(): for non-integer responses,
  # dpois() warns once per observation, which takes most of the time of a fit:
  fam_start <- if (fam == 1L) poisson() else binomial()
  fam_start$aic <- function(...) 0
  g0 <- suppressWarnings(glm.fit(des$x, des$y, weights = des$w,
                                 family = fam_start))
  beta0 <- coef(g0)
  if (anyNA(beta0) || !all(is.finite(beta0))) beta0 <- rep(0, p)
  cache <- new.env(parent = emptyenv())
  cache$u <- rep(0, sum(des$nlev))
  run <- function(par) {
    ev <- mlvl_laplace_c(par[seq_len(k)], par[k + seq_len(p)], des$x,
                         des$grp, des$nlev, des$y, des$w, offs, cache$u,
                         fam, TRUE)
    if (is.finite(ev$D)) cache$u <- ev$u
    cache$par <- par
    cache$ev <- ev
    ev
  }
  fn <- function(par) run(par)$D
  # nlminb() asks for the gradient only at points with a finite objective; at
  # other points, the gradient is `NaN`, and nlminb() stops with an error:
  gr <- function(par) {
    if (!identical(cache$par, par)) run(par)
    as.vector(cache$ev$grad %||% rep(NaN, length(par)))
  }
  lower <- c(rep(0, k), rep(-Inf, p))
  ctrl <- list(rel.tol = rel_tol, iter.max = 300, eval.max = 600)
  opt <- nlminb(c(rep(1, k), beta0), fn, gr, lower = lower, control = ctrl)
  n_restart <- 0L
  small <- opt$par[seq_len(k)] < small_theta
  if (any(small) && restart_theta > 0) {
    n_restart <- 1L
    start2 <- opt$par
    start2[which(small)] <- restart_theta
    cache$u <- rep(0, sum(des$nlev))
    opt2 <- nlminb(start2, fn, gr, lower = lower, control = ctrl)
    if (is.finite(opt2$objective) && opt2$objective < opt$objective) {
      opt <- opt2
    }
  }
  # The conditional modes at the optimum:
  theta <- opt$par[seq_len(k)]
  beta <- opt$par[k + seq_len(p)]
  ev <- mlvl_laplace_c(theta, beta, des$x, des$grp, des$nlev, des$y, des$w,
                       offs, cache$u, fam, FALSE)
  return(list(theta = theta, beta = beta, u = ev$u, deviance = ev$D,
              convergence = opt$convergence, message = opt$message,
              iterations = opt$iterations, n_restart = n_restart))
}

# The fit object.  `sd`: standard deviation of the random intercepts per
# grouping variable; `ranef`: the random intercepts per grouping variable,
# named by level (both in the term order of lme4).
new_projpred_mlvl <- function(formula, family, link, terms, xlevels,
                              contrasts, coefficients, sd, ranef, deviance,
                              convergence, optinfo = list()) {
  return(structure(
    list(formula = formula, family = family, link = link, terms = terms,
         xlevels = xlevels, contrasts = contrasts, coefficients = coefficients,
         sd = sd, ranef = ranef, deviance = deviance,
         convergence = convergence, optinfo = optinfo),
    class = "projpred_mlvl"
  ))
}

# Whether fit_glmer_callback() uses mlvl_laplace(): `use` (global option
# `projpred.mlvl_laplace`, read by divmin() in the main R process) is `TRUE`,
# global option `projpred.PQL` is not set, the family and the group-level
# terms are supported, and `dot_args` has no lme4::glmer() arguments other
# than `data` and `weights`.
mlvl_laplace_use <- function(formula, family, dot_args, use = TRUE) {
  if (!isTRUE(use) || getOption("projpred.PQL", FALSE) ||
      is.null(mlvl_laplace_fam(family)) ||
      is.null(mlvl_laplace_grps(formula))) {
    return(FALSE)
  }
  glmer_args <- intersect(names(dot_args), methods::formalArgs(lme4::glmer))
  return(all(glmer_args %in% c("data", "weights")))
}

# The Laplace fit.  Returns `NULL` if the design does not apply (see
# mlvl_laplace_design()), and an object with `convergence != 0` if nlminb()
# failed or did not converge (then `optinfo$message` says why).  The design is
# in attribute "des" of an object with `convergence != 0`.
mlvl_laplace <- function(formula, data, family, weights = NULL) {
  fam <- mlvl_laplace_fam(family)
  grps <- mlvl_laplace_grps(formula)
  if (is.null(fam) || is.null(grps)) return(NULL)
  des <- mlvl_laplace_design(formula, data, weights, grps, fam)
  if (is.null(des)) return(NULL)
  opt <- tryCatch(mlvl_laplace_optim(des, fam), error = function(e) e)
  if (inherits(opt, "error")) {
    out <- new_projpred_mlvl(formula, family$family, family$link, des$terms,
                             des$xlevels, des$contrasts, NULL, NULL, NULL,
                             NA_real_, convergence = -1L,
                             optinfo = list(message = conditionMessage(opt)))
    return(structure(out, des = des))
  }
  ok <- all(is.finite(c(opt$theta, opt$beta, opt$deviance, opt$u)))
  k <- length(des$nlev)
  offs_u <- c(0L, cumsum(des$nlev))
  ranef <- lapply(seq_len(k), function(j) {
    setNames(opt$theta[j] * opt$u[offs_u[j] + seq_len(des$nlev[j])],
             des$lvls[[j]])
  })
  out <- new_projpred_mlvl(
    formula = formula, family = family$family, link = family$link,
    terms = des$terms, xlevels = des$xlevels, contrasts = des$contrasts,
    coefficients = setNames(opt$beta, colnames(des$x)),
    sd = setNames(opt$theta, names(des$lvls)),
    ranef = setNames(ranef, names(des$lvls)), deviance = opt$deviance,
    convergence = if (ok) opt$convergence else -1L,
    optinfo = opt[c("message", "iterations", "n_restart")]
  )
  if (out$convergence != 0) attr(out, "des") <- des
  return(out)
}

# Convert an lme4::glmer() fit with random intercepts only to class
# `projpred_mlvl` (for the fallback of fit_glmer_callback(), so that all draws
# of a submodel have the same class).  `failed`: the mlvl_laplace() object that
# did not converge; its design gives the deviance of the glmer() fit on the
# scale of mlvl_laplace().  Returns `fit` unchanged if it has other group-level
# terms.
mlvl_from_merMod <- function(fit, formula, failed) {
  if (!inherits(fit, "glmerMod")) return(fit)
  re <- lme4::ranef(fit, condVar = FALSE)
  des <- attr(failed, "des")
  if (!all(vapply(re, function(r) identical(names(r), "(Intercept)"),
                  logical(1))) ||
      !setequal(names(re), names(des$lvls))) {
    return(fit)
  }
  re <- re[names(des$lvls)]
  vc <- lme4::VarCorr(fit)
  sd <- vapply(vc[names(re)], function(v) unname(attr(v, "stddev")),
               numeric(1))
  fam <- stats::family(fit)
  # The family of glmer_family() has the suffix "_projpred":
  fam$family <- sub("_projpred$", "", fam$family)
  dev <- tryCatch(
    mlvl_laplace_dev(des, mlvl_laplace_fam(fam), sd,
                     unname(lme4::fixef(fit))),
    error = function(e) NA_real_
  )
  trms <- delete.response(terms(fit, fixed.only = TRUE))
  attr(trms, ".Environment") <- environment(formula)
  return(new_projpred_mlvl(
    formula = formula, family = fam$family, link = fam$link, terms = trms,
    xlevels = .getXlevels(terms(fit, fixed.only = TRUE), fit@frame),
    contrasts = attr(lme4::getME(fit, "X"), "contrasts"),
    coefficients = lme4::fixef(fit), sd = sd,
    ranef = lapply(re, function(r) setNames(r[["(Intercept)"]], rownames(r))),
    deviance = dev, convergence = if (check_conv_s(fit)) 0L else 1L,
    optinfo = list(message = paste0("lme4::glmer() fit (Laplace fit: ",
                                    failed$optinfo$message, ")"))
  ))
}

# Linear predictors of `projpred_mlvl` fits (the draws of one submodel) for
# `newdata`: X beta_s plus the random intercepts of the existing group levels;
# for new group levels (and for all levels if global option
# `projpred.mlvl_pred_new` is `TRUE`), random intercepts drawn from their
# Gaussian distribution (per draw, per grouping variable, per level in the
# order of the factor levels), as repair_re.merMod() does for lme4 fits (which
# uses another order of the draws).
subprd_mlvl <- function(fits, newdata) {
  stopifnot(!is.null(newdata))
  # Cumulative fits have their own function (subprd_mlvl_cumul()):
  stopifnot(!inherits(fits[[1]], "projpred_mlvl_cumul"))
  pred_new <- getOption("projpred.mlvl_pred_new", FALSE)
  # Group levels, as in repair_re.merMod():
  grps <- names(fits[[1]]$sd)
  for (g in grps) {
    if (!g %in% names(newdata)) {
      stop("Could not find column `", g, "` in `newdata`.")
    }
  }
  lvls_new <- lapply(setNames(nm = grps), function(g) {
    levels(factor(newdata[[g]]))
  })
  lvls_chk <- if (pred_new) {
    lvls_new
  } else {
    lapply(setNames(nm = grps), function(g) {
      union(names(fits[[1]]$ranef[[g]]), lvls_new[[g]])
    })
  }
  if (length(lvls_chk) >= 2 &&
      anyDuplicated(unlist(lvls_chk, use.names = FALSE))) {
    stop("Currently, projpred requires all variables with group-level effects ",
         "to have disjoint level sets.")
  }
  prd_list <- lapply(fits, function(fit) {
    mfr <- model.frame(fit$terms, newdata, na.action = na.pass,
                       xlev = fit$xlevels)
    x <- model.matrix(fit$terms, mfr, contrasts.arg = fit$contrasts)
    eta <- drop(x[, names(fit$coefficients), drop = FALSE] %*%
                  fit$coefficients)
    for (g in names(fit$sd)) {
      lv <- as.character(newdata[[g]])
      b <- fit$ranef[[g]]
      lv_new <- if (pred_new) lvls_new[[g]] else setdiff(lvls_new[[g]], names(b))
      b_lv <- if (pred_new) rep(0, length(lv)) else unname(b[lv])
      b_lv[is.na(b_lv)] <- 0
      if (length(lv_new)) {
        b_new <- setNames(rnorm(length(lv_new), sd = fit$sd[[g]]), lv_new)
        is_new <- lv %in% lv_new
        b_lv[is_new] <- b_new[lv[is_new]]
      }
      eta <- eta + b_lv
    }
    eta
  })
  prd <- do.call(cbind, prd_list)
  dimnames(prd) <- list(rownames(newdata), NULL)
  return(prd)
}

#' Convert a submodel fit of class `projpred_mlvl` to class `glmerMod`
#'
#' For binomial and Poisson multilevel submodels with random intercepts,
#' \pkg{projpred} uses its own fit, which has class `projpred_mlvl` (see section
#' "Draw-wise divergence minimizers" of [projpred-package]). This function
#' builds the corresponding \pkg{lme4} object, for use with the methods of
#' \pkg{lme4}. It does not fit the submodel again: it evaluates the objective of
#' [lme4::glmer()] (with `nAGQ = 1`) once, at the parameter estimates of
#' `object`.
#'
#' @param object An object of class `projpred_mlvl`, for example an element of
#'   element `outdmin` of an object returned by [project()].
#' @param data A `data.frame` with the data that `object` was fitted to: the
#'   predictors, the grouping variables and the response, with the names used
#'   in `formula(object)`. The object does not store these data. In the
#'   traditional projection, the response of a submodel fit is the vector of
#'   the reference model's predicted means (success probabilities for the
#'   binomial family) for the projected draw or cluster. In the augmented-data
#'   projection, the data are the augmented rows (one row per observation and
#'   response category), and the predicted probabilities are in the weights.
#' @param weights A numeric vector of the observation weights that `object` was
#'   fitted with, or `NULL` for weights of 1. Rows of `data` with weight 0 are
#'   removed, as \pkg{projpred} does before it fits a multilevel submodel.
#'
#' @details The group-level intercepts of the returned object are computed by
#'   \pkg{lme4} from `data` and `weights`. If they or the value of the
#'   objective function at the parameter estimates differ from those of
#'   `object`, then `data` or `weights` are not those of the fit, and this
#'   function throws a warning.
#'
#'   For non-integer responses, \pkg{lme4} versions up to 2.1-0 round binomial
#'   responses in the log-likelihood and give a log-likelihood of `-Inf` for
#'   the Poisson family. With these versions, the log-likelihood of the
#'   returned object (and what depends on it, for example [AIC()]) is not that
#'   of the fit; the parameter estimates and the predictions are not affected.
#'
#'   The call stored in the returned object is the [lme4::glmer()] call for
#'   `data` and `weights`, with the default settings of [lme4::glmer()].
#'
#' @return An object of class `glmerMod` (see [lme4::glmer()]) with the
#'   fixed effects, the group-level standard deviations and the group-level
#'   intercepts of `object`.
#'
#' @examplesIf requireNamespace("rstanarm", quietly = TRUE)
#' # Data:
#' set.seed(1)
#' dat <- data.frame(x1 = rnorm(120), x2 = rnorm(120),
#'                   g = factor(rep(1:8, each = 15)))
#' dat$y <- rbinom(120, 1, plogis(0.8 * dat$x1 + rnorm(8)[dat$g]))
#'
#' # The `stanreg` fit which will be used as the reference model (with small
#' # values for `chains` and `iter`, but only for technical reasons in this
#' # example; this is not recommended in general):
#' fit <- rstanarm::stan_glmer(
#'   y ~ x1 + x2 + (1 | g), family = binomial(), data = dat,
#'   chains = 2, iter = 500, refresh = 0, seed = 9876
#' )
#'
#' # Projection onto a multilevel submodel, with a single cluster (so the
#' # response of the submodel fit is the mean of the reference model's
#' # predicted success probabilities over the posterior draws):
#' refm <- suppressWarnings(get_refmodel(fit))
#' prj <- project(refm, predictor_terms = c("x1", "(1 | g)"), nclusters = 1)
#' fit_sub <- prj$outdmin[[1]]
#'
#' # Conversion (the name of the response is taken from the formula of the
#' # submodel fit):
#' dat_sub <- dat
#' dat_sub[[all.vars(formula(fit_sub))[1]]] <- rowMeans(refm$mu)
#' fit_lme4 <- as_glmerMod(fit_sub, data = dat_sub, weights = refm$wobs)
#' lme4::ranef(fit_lme4)
#'
#' @export
as_glmerMod <- function(object, data, weights = NULL) {
  if (!inherits(object, "projpred_mlvl")) {
    stop("`object` must be of class `projpred_mlvl`.")
  }
  if (inherits(object, "projpred_mlvl_cumul")) {
    stop("as_glmerMod() does not support fits of the cumulative family.")
  }
  if (is.null(object$coefficients) || is.null(object$sd)) {
    stop("`object` has no parameter estimates.")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a `data.frame`.")
  }
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (!is.numeric(weights) || length(weights) != nrow(data) ||
      anyNA(weights) || any(weights < 0)) {
    stop("`weights` must be a numeric vector of non-negative values with one ",
         "element per row of `data`.")
  }
  if (".projpred_wobs" %in% names(data)) {
    stop("`data` must not have a column named `.projpred_wobs`.")
  }
  fam <- switch(object$family,
                binomial = binomial(link = object$link),
                poisson = poisson(link = object$link),
                stop("Unexpected family of `object`."))
  # The call of the returned object: the lme4::glmer() call for these data.
  mc_in <- match.call()
  mc <- as.call(c(
    list(quote(lme4::glmer), formula = object$formula, data = mc_in$data,
         family = call(object$family, link = object$link)),
    if (!is.null(mc_in$weights)) list(weights = mc_in$weights)
  ))
  # As divmin() does before the fit:
  if (any(weights == 0)) {
    data <- data[weights != 0, , drop = FALSE]
    weights <- weights[weights != 0]
  }
  # The objective of mlvl_laplace() for these data, at the estimates of
  # `object`:
  des <- mlvl_laplace_design(object$formula, data, weights,
                             mlvl_laplace_grps(object$formula),
                             mlvl_laplace_fam(fam))
  dev <- if (is.null(des)) {
    NA_real_
  } else {
    tryCatch(
      mlvl_laplace_dev(des, mlvl_laplace_fam(fam),
                       unname(object$sd[names(des$lvls)]),
                       unname(object$coefficients[colnames(des$x)])),
      error = function(e) NA_real_
    )
  }
  data$.projpred_wobs <- weights
  # The responses of a projection are not integers, so the warnings of the
  # family functions about non-integer responses are expected:
  out <- withCallingHandlers({
    glmod <- do.call(lme4::glFormula, list(
      formula = object$formula, data = data, family = fam,
      weights = quote(.projpred_wobs)
    ))
    nms_x <- colnames(glmod$X)
    nms_g <- names(glmod$reTrms$cnms)
    if (!setequal(nms_x, names(object$coefficients)) ||
        !setequal(nms_g, names(object$sd))) {
      stop("The fixed effects or the grouping variables that lme4 builds ",
           "from `data` are not those of `object`.")
    }
    devfun <- do.call(lme4::mkGlmerDevfun, c(glmod, list(nAGQ = 0L)))
    devfun <- lme4::updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = 1L)
    par <- unname(c(object$sd[nms_g], object$coefficients[nms_x]))
    opt <- list(par = par, fval = devfun(par), conv = 0L, feval = 1L,
                message = "parameter estimates of projpred's Laplace fit")
    lme4::mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr,
                   mc = mc)
  }, warning = function(w) {
    if (grepl("^non-integer ", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
  # `data` and `weights` are those of the fit if they give the objective value
  # and the group-level intercepts of `object`:
  re <- lme4::ranef(out, condVar = FALSE)
  differs <- vapply(nms_g, function(g) {
    b <- object$ranef[[g]]
    b_out <- setNames(re[[g]][["(Intercept)"]], rownames(re[[g]]))
    !setequal(names(b), names(b_out)) ||
      max(abs(b - b_out[names(b)])) > 1e-3 * (1 + max(abs(b)))
  }, logical(1))
  dev_differs <- is.finite(object$deviance) &&
    (!is.finite(dev) ||
       abs(dev - object$deviance) > 1e-6 * (1 + abs(object$deviance)))
  if (any(differs) || dev_differs) {
    warning("The converted fit differs from `object` in the group-level ",
            "intercepts or in the value of the objective function. Check that ",
            "`data` (including the response) and `weights` are those that ",
            "`object` was fitted to.")
  }
  return(out)
}

#' @noRd
#' @export
formula.projpred_mlvl <- function(x, ...) {
  return(x$formula)
}

#' @noRd
#' @export
coef.projpred_mlvl <- function(object, ...) {
  return(object$coefficients)
}

# Linear predictor for `newdata` with the random intercepts of the fitted group
# levels and 0 for new levels:
#' @noRd
#' @export
predict.projpred_mlvl <- function(object, newdata, ...) {
  mfr <- model.frame(object$terms, newdata, na.action = na.pass,
                     xlev = object$xlevels)
  x <- model.matrix(object$terms, mfr, contrasts.arg = object$contrasts)
  eta <- drop(x[, names(object$coefficients), drop = FALSE] %*%
                object$coefficients)
  for (g in names(object$sd)) {
    b <- unname(object$ranef[[g]][as.character(newdata[[g]])])
    b[is.na(b)] <- 0
    eta <- eta + b
  }
  return(setNames(eta, rownames(newdata)))
}

#' @noRd
#' @export
print.projpred_mlvl <- function(x, ...) {
  cat("Multilevel submodel fit (projpred's Laplace fit)\n")
  cat("Family: ", x$family, " (link: ", x$link, ")\n", sep = "")
  cat("Formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat("Fixed effects:\n")
  print(x$coefficients)
  cat("Standard deviations of the group-level intercepts:\n")
  print(x$sd)
  return(invisible(x))
}
