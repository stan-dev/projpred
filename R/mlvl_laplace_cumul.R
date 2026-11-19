# Laplace fit of cumulative-link multilevel submodels with one group-level
# intercept term
#
# mlvl_laplace_cumul() fits the submodels of the augmented-data projection for
# the brms::cumulative() family whose only group-level term is `(1 | <group>)`.
# The model and the objective are those of ordinal::clmm() (the Laplace
# approximation, see `src/mlvlcumul.cpp`): with thresholds alpha, coefficients
# beta and group-level intercepts b,
#   P(y <= j) = F(alpha_j - x' beta - b_g).
# For the links cloglog and loglog, they are those of the code of clmm() for
# `useMatrix = TRUE`; its default code for a single group-level intercept term
# swaps these two links.
# It minimizes the objective with nlminb() and the exact gradient. The returned
# object has class `c("projpred_mlvl_cumul", "projpred_mlvl")`.
# fit_cumul_mlvl() calls it; other multilevel submodels of this family keep
# ordinal::clmm().

# Link code of mlvl_cumul_c() and link name in the names of ordinal::clmm(), or
# `NULL` if mlvl_laplace_cumul() does not handle `link`:
mlvl_cumul_link <- function(link) {
  if (identical(link, "probit_approx")) link <- "probit"
  links <- c(logit = 1L, probit = 2L, loglog = 3L, cloglog = 4L, cauchit = 5L)
  if (!is.character(link) || length(link) != 1 || !link %in% names(links)) {
    return(NULL)
  }
  return(list(code = links[[link]], name = link))
}

# Quantile function of the latent distribution of link `name`:
mlvl_cumul_quantile <- function(p, name) {
  switch(name,
         logit = qlogis(p),
         probit = qnorm(p),
         loglog = -log(-log(p)),
         cloglog = log(-log1p(-p)),
         cauchit = qcauchy(p))
}

# Design of mlvl_laplace_cumul(): the model matrix without the intercept, the
# category codes of the response, the weights and the level codes of the
# grouping variable. The rows are those of `data` (in the augmented-data
# projection: one row per observation and response category, with the
# probability of the category as the weight). Returns `NULL` for cases that the
# fit does not handle: group-level terms other than one `(1 | <group>)`, a
# response that is not a factor, missing values, offsets, a rank-deficient
# model matrix, a response category without weight.
mlvl_cumul_design <- function(formula, data, weights) {
  bars <- reformulas::findbars(formula)
  if (length(bars) != 1 || !identical(bars[[1]][[2]], 1) ||
      !is.name(bars[[1]][[3]])) {
    return(NULL)
  }
  grp_nm <- as.character(bars[[1]][[3]])
  if (!is.data.frame(data) || !grp_nm %in% names(data)) return(NULL)
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
  cons <- attr(x, "contrasts")
  y <- model.response(mf)
  if (!is.factor(y) || nlevels(y) < 2 || anyNA(y) || anyNA(x) ||
      !all(is.finite(x)) || qr(x)$rank < ncol(x)) {
    return(NULL)
  }
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  if (is.null(weights)) weights <- rep(1, n)
  if (length(weights) != n || !all(is.finite(weights)) || any(weights < 0)) {
    return(NULL)
  }
  w_cat <- vapply(seq_len(nlevels(y)), function(j) {
    sum(weights[as.integer(y) == j])
  }, numeric(1))
  if (!all(w_cat > 0)) return(NULL)
  fac <- droplevels(as.factor(data[[grp_nm]]))
  if (anyNA(fac)) return(NULL)
  return(list(
    x = unname(x), x_nms = colnames(x), y = as.integer(y), lev = levels(y),
    w = as.vector(weights), w_cat = w_cat, grp = as.integer(fac) - 1L,
    nlev = nlevels(fac), lvls = levels(fac), grp_nm = grp_nm,
    terms = delete.response(trms), xlevels = .getXlevels(trms, mf),
    contrasts = cons
  ))
}

# Objective of the design `des` at (alpha, beta, tau), with the conditional
# modes found from u = 0:
mlvl_cumul_obj <- function(des, link_code, alpha, beta, tau) {
  return(mlvl_cumul_c(alpha, beta, tau, des$x, des$grp, des$nlev, des$y,
                      des$w, rep(0, des$nlev), link_code, FALSE)$L)
}

# Minimize the objective over (alpha, beta, tau) with nlminb(). Start: the
# thresholds from the weighted proportions of the categories, beta 0, tau 1.
# The objective is even in tau, so tau = 0 is a stationary point where the
# gradient-based search can stop: if tau ends below `small_tau`, restart from
# `restart_tau` and keep the lower objective. Each evaluation starts the search
# for the conditional modes at those of the last evaluation. The fit has not
# converged (code 2) if the gradient with respect to the conditional modes is
# above `tol_modes` at the optimum.
mlvl_cumul_optim <- function(des, link, rel_tol = 1e-10, small_tau = 1e-3,
                             restart_tau = 0.3, tol_modes = 1e-6) {
  n_thres <- length(des$lev) - 1L
  p <- ncol(des$x)
  idx_alpha <- seq_len(n_thres)
  idx_beta <- n_thres + seq_len(p)
  idx_tau <- n_thres + p + 1L
  cum_prop <- cumsum(des$w_cat) / sum(des$w_cat)
  alpha0 <- mlvl_cumul_quantile(cum_prop[idx_alpha], link$name)
  cache <- new.env(parent = emptyenv())
  cache$u <- rep(0, des$nlev)
  run <- function(par) {
    ev <- mlvl_cumul_c(par[idx_alpha], par[idx_beta], par[idx_tau], des$x,
                       des$grp, des$nlev, des$y, des$w, cache$u, link$code,
                       TRUE)
    if (is.finite(ev$L)) cache$u <- ev$u
    cache$par <- par
    cache$ev <- ev
    ev
  }
  fn <- function(par) run(par)$L
  # nlminb() asks for the gradient only at points with a finite objective; at
  # other points, the gradient is `NaN`, and nlminb() stops with an error:
  gr <- function(par) {
    if (!identical(cache$par, par)) run(par)
    as.vector(cache$ev$grad %||% rep(NaN, length(par)))
  }
  lower <- c(rep(-Inf, n_thres + p), 0)
  ctrl <- list(rel.tol = rel_tol, iter.max = 300, eval.max = 600)
  opt <- nlminb(c(alpha0, rep(0, p), 1), fn, gr, lower = lower,
                control = ctrl)
  n_restart <- 0L
  if (opt$par[idx_tau] < small_tau && restart_tau > 0) {
    n_restart <- 1L
    start2 <- opt$par
    start2[idx_tau] <- restart_tau
    cache$u <- rep(0, des$nlev)
    opt2 <- nlminb(start2, fn, gr, lower = lower, control = ctrl)
    if (is.finite(opt2$objective) && opt2$objective < opt$objective) {
      opt <- opt2
    }
  }
  # The conditional modes at the optimum:
  alpha <- opt$par[idx_alpha]
  beta <- opt$par[idx_beta]
  tau <- opt$par[idx_tau]
  ev <- mlvl_cumul_c(alpha, beta, tau, des$x, des$grp, des$nlev, des$y, des$w,
                     cache$u, link$code, FALSE)
  convergence <- opt$convergence
  message <- opt$message
  if (convergence == 0 && !is.finite(ev$L)) {
    convergence <- 2L
    message <- "the objective is not finite at the optimum"
  } else if (convergence == 0 && !isTRUE(ev$gmax < tol_modes)) {
    convergence <- 2L
    message <- "the search for the conditional modes did not converge"
  }
  return(list(alpha = alpha, beta = beta, tau = tau, u = ev$u, objective = ev$L,
              convergence = convergence, message = message,
              iterations = opt$iterations, n_restart = n_restart))
}

# The fit object. `alpha`: thresholds, named as by ordinal::clmm();
# `coefficients`: regression coefficients (no intercept); `sd`: standard
# deviation of the group-level intercepts; `ranef`: the group-level intercepts,
# named by level; `deviance`: twice the objective of `src/mlvlcumul.cpp`.
new_projpred_mlvl_cumul <- function(formula, link, des, alpha, coefficients,
                                    sd, ranef, deviance, convergence,
                                    optinfo = list()) {
  n_cat <- length(des$lev)
  nms <- function(v, nm) if (is.null(v)) NULL else setNames(v, nm)
  return(structure(
    list(formula = formula, family = "cumulative", link = link,
         terms = des$terms, xlevels = des$xlevels, contrasts = des$contrasts,
         lev = des$lev,
         alpha = nms(alpha, paste(des$lev[-n_cat], des$lev[-1], sep = "|")),
         coefficients = nms(coefficients, des$x_nms),
         sd = nms(sd, des$grp_nm),
         ranef = if (is.null(ranef)) {
           NULL
         } else {
           setNames(list(setNames(ranef, des$lvls)), des$grp_nm)
         },
         deviance = deviance, convergence = convergence, optinfo = optinfo),
    class = c("projpred_mlvl_cumul", "projpred_mlvl")
  ))
}

# The Laplace fit. Returns `NULL` if the link or the design does not apply (see
# mlvl_cumul_design()), and an object with `convergence != 0` if nlminb()
# failed or did not converge (then `optinfo$message` says why). The design is
# in attribute "des" of an object with `convergence != 0`.
mlvl_laplace_cumul <- function(formula, data, family, weights = NULL) {
  link <- mlvl_cumul_link(family$link)
  if (is.null(link)) return(NULL)
  des <- mlvl_cumul_design(formula, data, weights)
  if (is.null(des)) return(NULL)
  opt <- tryCatch(mlvl_cumul_optim(des, link), error = function(e) e)
  if (inherits(opt, "error")) {
    out <- new_projpred_mlvl_cumul(
      formula, link$name, des, NULL, NULL, NULL, NULL, NA_real_,
      convergence = -1L, optinfo = list(message = conditionMessage(opt))
    )
    return(structure(out, des = des))
  }
  ok <- all(is.finite(c(opt$alpha, opt$beta, opt$tau, opt$u, opt$objective)))
  out <- new_projpred_mlvl_cumul(
    formula, link$name, des, alpha = opt$alpha, coefficients = opt$beta,
    sd = opt$tau, ranef = opt$tau * opt$u, deviance = 2 * opt$objective,
    convergence = if (ok) opt$convergence else -1L,
    optinfo = opt[c("message", "iterations", "n_restart")]
  )
  if (out$convergence != 0) attr(out, "des") <- des
  return(out)
}

# Convert an ordinal::clmm() fit with one group-level intercept term to class
# `projpred_mlvl_cumul` (for the fallback of fit_cumul_mlvl(), so that all
# draws of a submodel have the same class). `failed`: the mlvl_laplace_cumul()
# object that did not converge; its design gives the names and the objective
# value of the clmm() fit on the scale of mlvl_laplace_cumul(). Returns `fit`
# unchanged if its parameters are not those of the design.
mlvl_cumul_from_clmm <- function(fit, formula, failed) {
  if (!inherits(fit, "clmm")) return(fit)
  des <- attr(failed, "des")
  re <- ordinal::ranef(fit)
  vc <- ordinal::VarCorr(fit)
  beta <- fit$beta %||% numeric()
  if (!identical(names(re), des$grp_nm) ||
      !identical(names(re[[1]]), "(Intercept)") ||
      !setequal(rownames(re[[1]]), des$lvls) ||
      length(fit$alpha) != length(des$lev) - 1L ||
      !setequal(names(beta), des$x_nms)) {
    return(fit)
  }
  beta <- unname(beta[des$x_nms])
  alpha <- unname(fit$alpha)
  tau <- unname(attr(vc[[des$grp_nm]], "stddev"))
  link <- mlvl_cumul_link(fit$link)
  obj <- mlvl_cumul_obj(des, link$code, alpha, beta, tau)
  return(new_projpred_mlvl_cumul(
    formula, link$name, des, alpha = alpha, coefficients = beta, sd = tau,
    ranef = re[[1]][des$lvls, "(Intercept)"], deviance = 2 * obj,
    convergence = if (check_conv_s(fit)) 0L else 1L,
    optinfo = list(message = paste0("ordinal::clmm() fit (Laplace fit: ",
                                    failed$optinfo$message, ")"))
  ))
}

# Linear predictors of `projpred_mlvl_cumul` fits (the draws of one submodel)
# for `newdata`, on the link scale that subprd_augdat() uses: for response
# category j (without the last one), alpha_j - x' beta - b_g. For new group
# levels (and for all levels if global option `projpred.mlvl_pred_new` is
# `TRUE`), the group-level intercepts are drawn from their Gaussian
# distribution, with the random numbers that repair_re.clmm() uses for
# ordinal::clmm() fits (per draw, one number per level: first the fitted
# levels in their order, if they are drawn, then the new levels in the order
# in which they occur in `newdata`). Returns an array with the dimensions
# observations, categories (without the last one) and draws, or (for `as_array
# = FALSE`) the list of the matrices of the draws.
subprd_mlvl_cumul <- function(fits, newdata, as_array = TRUE) {
  stopifnot(!is.null(newdata))
  pred_new <- getOption("projpred.mlvl_pred_new", FALSE)
  fit_1 <- fits[[1]]
  grp_nm <- names(fit_1$sd)
  if (!grp_nm %in% names(newdata)) {
    stop("Could not find column `", grp_nm, "` in `newdata`.")
  }
  mfr <- model.frame(fit_1$terms, newdata, na.action = na.pass,
                     xlev = fit_1$xlevels)
  x <- model.matrix(fit_1$terms, mfr, contrasts.arg = fit_1$contrasts)
  lv <- as.character(newdata[[grp_nm]])
  lvls_new <- unique(lv)
  n <- nrow(x)
  prd_list <- lapply(fits, function(fit) {
    x_s <- x[, names(fit$coefficients), drop = FALSE]
    eta <- drop(x_s %*% fit$coefficients)
    b <- fit$ranef[[grp_nm]]
    lvls_draw <- setdiff(lvls_new, names(b))
    if (pred_new) lvls_draw <- c(intersect(names(b), lvls_new), lvls_draw)
    b_lv <- if (pred_new) rep(0, n) else unname(b[lv])
    b_lv[is.na(b_lv)] <- 0
    if (length(lvls_draw)) {
      b_draw <- setNames(rnorm(length(lvls_draw)) * fit$sd[[grp_nm]],
                         lvls_draw)
      is_draw <- lv %in% lvls_draw
      b_lv[is_draw] <- b_draw[lv[is_draw]]
    }
    prd <- matrix(fit$alpha, nrow = n, ncol = length(fit$alpha),
                  byrow = TRUE) - (eta + b_lv)
    # The dimnames of the ordinal::clmm() path:
    dimnames(prd) <- list(as.character(seq_len(n)), fit$lev[-length(fit$lev)])
    prd
  })
  if (!as_array) return(prd_list)
  return(do.call(abind::abind, list(prd_list, rev.along = 0)))
}

# Linear predictors alpha_j - x' beta - b_g for `newdata` (a matrix with one
# column per response category, without the last one), with the group-level
# intercepts of the fitted group levels and 0 for new levels:
#' @noRd
#' @export
predict.projpred_mlvl_cumul <- function(object, newdata, ...) {
  grp_nm <- names(object$sd)
  mfr <- model.frame(object$terms, newdata, na.action = na.pass,
                     xlev = object$xlevels)
  x <- model.matrix(object$terms, mfr, contrasts.arg = object$contrasts)
  eta <- drop(x[, names(object$coefficients), drop = FALSE] %*%
                object$coefficients)
  b <- unname(object$ranef[[grp_nm]][as.character(newdata[[grp_nm]])])
  b[is.na(b)] <- 0
  prd <- matrix(object$alpha, nrow = nrow(x), ncol = length(object$alpha),
                byrow = TRUE) - (eta + b)
  dimnames(prd) <- list(rownames(newdata), names(object$alpha))
  return(prd)
}

# The thresholds and the regression coefficients:
#' @noRd
#' @export
coef.projpred_mlvl_cumul <- function(object, ...) {
  return(c(object$alpha, object$coefficients))
}

#' @noRd
#' @export
print.projpred_mlvl_cumul <- function(x, ...) {
  cat("Multilevel submodel fit (projpred's Laplace fit)\n")
  cat("Family: ", x$family, " (link: ", x$link, ")\n", sep = "")
  cat("Formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat("Thresholds:\n")
  print(x$alpha)
  cat("Coefficients:\n")
  print(x$coefficients)
  cat("Standard deviation of the group-level intercepts:\n")
  print(x$sd)
  return(invisible(x))
}
