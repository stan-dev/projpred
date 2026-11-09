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

# Design of mlvl_laplace(): the fixed-effects model matrix, the response, the
# weights and the level codes of the grouping variables, in the term order of
# lme4 (reformulas::mkReTrms() sorts the terms by decreasing number of levels if
# they are not in this order).  Returns `NULL` for cases that the fit does not
# handle (missing values, offsets, a rank-deficient model matrix, a response
# outside the range of the family).
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
  gr <- function(par) {
    if (!identical(cache$par, par)) run(par)
    cache$ev$grad %||% rep(NaN, length(par))
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

# Whether fit_glmer_callback() uses mlvl_laplace(): the global option
# `projpred.mlvl_laplace` is `TRUE`, `projpred.PQL` is not set, the family and
# the group-level terms are supported, and `dot_args` has no lme4::glmer()
# arguments other than `data` and `weights`.
mlvl_laplace_use <- function(formula, family, dot_args) {
  if (!getOption("projpred.mlvl_laplace", TRUE) ||
      getOption("projpred.PQL", FALSE) || is.null(mlvl_laplace_fam(family)) ||
      is.null(mlvl_laplace_grps(formula))) {
    return(FALSE)
  }
  glmer_args <- intersect(names(dot_args), methods::formalArgs(lme4::glmer))
  return(all(glmer_args %in% c("data", "weights")))
}

# The Laplace fit.  Returns `NULL` if the design does not apply (see
# mlvl_laplace_design()) and the fit with `convergence != 0` if nlminb() did not
# converge or gave non-finite values.
mlvl_laplace <- function(formula, data, family, weights = NULL) {
  fam <- mlvl_laplace_fam(family)
  grps <- mlvl_laplace_grps(formula)
  if (is.null(fam) || is.null(grps)) return(NULL)
  des <- mlvl_laplace_design(formula, data, weights, grps, fam)
  if (is.null(des)) return(NULL)
  opt <- tryCatch(mlvl_laplace_optim(des, fam), error = function(e) NULL)
  if (is.null(opt)) {
    return(new_projpred_mlvl(formula, family$family, family$link, des$terms,
                             des$xlevels, des$contrasts, NULL, NULL, NULL,
                             NA_real_, convergence = -1L))
  }
  ok <- all(is.finite(c(opt$theta, opt$beta, opt$deviance, opt$u)))
  k <- length(des$nlev)
  offs_u <- c(0L, cumsum(des$nlev))
  ranef <- lapply(seq_len(k), function(j) {
    setNames(opt$theta[j] * opt$u[offs_u[j] + seq_len(des$nlev[j])],
             des$lvls[[j]])
  })
  return(new_projpred_mlvl(
    formula = formula, family = family$family, link = family$link,
    terms = des$terms, xlevels = des$xlevels, contrasts = des$contrasts,
    coefficients = setNames(opt$beta, colnames(des$x)),
    sd = setNames(opt$theta, names(des$lvls)),
    ranef = setNames(ranef, names(des$lvls)), deviance = opt$deviance,
    convergence = if (ok) opt$convergence else -1L,
    optinfo = opt[c("message", "iterations", "n_restart")]
  ))
}

# Convert an lme4::glmer() fit with random intercepts only to class
# `projpred_mlvl` (for the fallback of fit_glmer_callback(), so that all draws
# of a submodel have the same class).  Returns `fit` unchanged if it has other
# group-level terms.
mlvl_from_merMod <- function(fit, formula) {
  if (!inherits(fit, "glmerMod")) return(fit)
  re <- lme4::ranef(fit, condVar = FALSE)
  if (!all(vapply(re, function(r) identical(names(r), "(Intercept)"),
                  logical(1)))) {
    return(fit)
  }
  vc <- lme4::VarCorr(fit)
  trms <- delete.response(terms(fit, fixed.only = TRUE))
  attr(trms, ".Environment") <- environment(formula)
  fam <- stats::family(fit)
  return(new_projpred_mlvl(
    formula = formula, family = fam$family, link = fam$link, terms = trms,
    xlevels = .getXlevels(terms(fit, fixed.only = TRUE), fit@frame),
    contrasts = attr(lme4::getME(fit, "X"), "contrasts"),
    coefficients = lme4::fixef(fit),
    sd = vapply(vc[names(re)], function(v) unname(attr(v, "stddev")),
                numeric(1)),
    ranef = lapply(re, function(r) setNames(r[["(Intercept)"]], rownames(r))),
    deviance = stats::deviance(fit),
    convergence = if (check_conv_s(fit)) 0L else 1L,
    optinfo = list(message = "lme4::glmer() fit")
  ))
}

# Linear predictors of `projpred_mlvl` fits (the draws of one submodel) for
# `newdata`: X beta_s plus the random intercepts of the existing group levels;
# for new group levels (and for all levels if global option
# `projpred.mlvl_pred_new` is `TRUE`), random intercepts drawn from their
# Gaussian distribution, as repair_re.merMod() does for lme4 fits.
subprd_mlvl <- function(fits, newdata) {
  stopifnot(!is.null(newdata))
  pred_new <- getOption("projpred.mlvl_pred_new", FALSE)
  # Group levels, as in repair_re.merMod():
  grps <- names(fits[[1]]$sd)
  for (g in grps) {
    if (!g %in% names(newdata)) {
      stop("Could not find column `", g, "` in `newdata`.")
    }
  }
  lvls_new <- lapply(setNames(nm = grps), function(g) {
    as.character(unique(newdata[[g]]))
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
