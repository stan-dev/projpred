# Divergence minimizers ---------------------------------------------------

## Traditional (and latent) projection ------------------------------------

# Needed to avoid a NOTE in `R CMD check`:
if (getRversion() >= package_version("2.15.1")) {
  utils::globalVariables("formula_s")
  utils::globalVariables("projpred_var_s")
  utils::globalVariables("projpred_formula_no_random_s")
}

divmin <- function(
    formula,
    projpred_var,
    projpred_ws_aug,
    verbose_divmin,
    throw_warn_sdivmin = getOption("projpred.warn_proj_drawwise", TRUE),
    do_check_conv = getOption("projpred.check_convergence", TRUE),
    ...
) {
  if (missing(throw_warn_sdivmin) &&
      is.null(getOption("projpred.warn_proj_drawwise")) &&
      !is.null(getOption("projpred.warn_prj_drawwise"))) {
    warning(
      "Global option `projpred.warn_prj_drawwise` is deprecated. Please use ",
      "global option `projpred.warn_proj_drawwise` instead. Now using the ",
      "value from global option `projpred.warn_prj_drawwise`."
    )
    throw_warn_sdivmin <- getOption("projpred.warn_prj_drawwise")
  }
  if (missing(do_check_conv) &&
      is.null(getOption("projpred.check_convergence")) &&
      !is.null(getOption("projpred.check_conv"))) {
    warning(
      "Global option `projpred.check_conv` is deprecated. Please use ",
      "global option `projpred.check_convergence` instead. Now using the ",
      "value from global option `projpred.check_conv`."
    )
    do_check_conv <- getOption("projpred.check_conv")
  }
  trms_all <- extract_terms_response(formula)
  has_grp <- length(trms_all$group_terms) > 0
  has_add <- length(trms_all$additive_terms) > 0
  projpred_formulas_no_random <- NA
  projpred_random <- NA
  # Define sdivmin(), the divergence minimizer for each draw s = 1, ..., S (and
  # perform other actions, if necessary):
  if (!has_grp && !has_add) {
    sdivmin <- get(getOption("projpred.glm_fitter", "fit_glm_ridge_callback"),
                   mode = "function")
  } else if (has_grp && !has_add) {
    sdivmin <- fit_glmer_callback
    if (getOption("projpred.PQL", FALSE)) {
      # Split up the formula into a fixed and a random part (note: we could also
      # use reformulas::nobars() and reformulas::findbars() here):
      formula_random <- split_formula_random_gamm4(formula)
      projpred_formulas_no_random <- validate_response_formula(
        formula_random$formula
      )
      projpred_random <- formula_random$random
    }
  } else if (!has_grp && has_add) {
    sdivmin <- fit_gam_callback
  } else if (has_grp && has_add) {
    sdivmin <- fit_gamm_callback
    formula_random <- split_formula_random_gamm4(formula)
    projpred_formulas_no_random <- validate_response_formula(
      formula_random$formula
    )
    projpred_random <- formula_random$random
  }
  formulas <- validate_response_formula(formula)
  # Gaussian multilevel submodels: one lFormula() for all draws (see
  # lmer_shared_callback()). The options and the lme4 check run here, because
  # the draw-wise calls may run in parallel workers:
  if (identical(sdivmin, fit_glmer_callback) && length(formulas) > 1 &&
      !getOption("projpred.PQL", FALSE) &&
      getOption("projpred.gaussian_not_as_generalized", TRUE) &&
      identical(list(...)$family$family, "gaussian") &&
      identical(list(...)$family$link, "identity") &&
      getOption("projpred.lmer_shared", TRUE) && lmer_shared_ok()) {
    sdivmin <- shared_setup_callback(lmer_shared_callback, formulas, sdivmin)
  }
  # Additive submodels without group-level terms: one mgcv setup for all draws
  # (see gam_shared_callback()):
  if (identical(sdivmin, fit_gam_callback) && length(formulas) > 1 &&
      getOption("projpred.gam_shared", TRUE)) {
    sdivmin <- shared_setup_callback(gam_shared_callback, formulas, sdivmin)
  }
  if (is.list(projpred_formulas_no_random)) {
    stopifnot(length(projpred_formulas_no_random) == length(formulas))
  } else {
    projpred_formulas_no_random <- as.list(
      rep(projpred_formulas_no_random, length(formulas))
    )
  }

  if (is.null(getOption("projpred.parallel_proj_trigger")) &&
      !is.null(getOption("projpred.prll_prj_trigger"))) {
    warning(
      "Global option `projpred.prll_prj_trigger` is deprecated. Please use ",
      "global option `projpred.parallel_proj_trigger` instead. Now using the ",
      "value from global option `projpred.prll_prj_trigger`."
    )
    prll_prj_trigger <- getOption("projpred.prll_prj_trigger")
  } else {
    prll_prj_trigger <- getOption("projpred.parallel_proj_trigger", Inf)
  }
  if (length(formulas) < prll_prj_trigger) {
    # Sequential case. Actually, we could simply use ``%do_projpred%` <-
    # foreach::`%do%`` here and then proceed as in the parallel case, but that
    # would require adding more "hard" dependencies (because packages 'foreach'
    # and 'iterators' would have to be moved from `Suggests:` to `Imports:`).
    if (verbose_divmin) {
      pb <- utils::txtProgressBar(min = 0, max = length(formulas), style = 3,
                                  initial = 0, file = stderr())
      on.exit(close(pb))
    }
    outdmin <- lapply(seq_along(formulas), function(s) {
      if (verbose_divmin) {
        on.exit(utils::setTxtProgressBar(pb, s))
      }
      mssgs_warns_capt <- capt_mssgs_warns(
        soutdmin <- sdivmin(
          formula = formulas[[s]],
          projpred_var = projpred_var[, s, drop = FALSE],
          projpred_formula_no_random = projpred_formulas_no_random[[s]],
          projpred_random = projpred_random,
          ...
        )
      )
      return(nlist(soutdmin, mssgs_warns_capt))
    })
  } else {
    # Parallel case.
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Please install the 'foreach' package.")
    }
    if (!requireNamespace("iterators", quietly = TRUE)) {
      stop("Please install the 'iterators' package.")
    }
    if (verbose_divmin && use_progressr()) {
      progressor_obj <- progressr::progressor(length(formulas))
    } else {
      progressor_obj <- NULL
    }
    dot_args <- list(...)
    `%do_projpred%` <- foreach::`%dopar%`
    outdmin <- foreach::foreach(
      formula_s = formulas,
      projpred_var_s = iterators::iter(projpred_var, by = "column"),
      projpred_formula_no_random_s = projpred_formulas_no_random,
      .packages = c("projpred"),
      .export = c("sdivmin", "projpred_random", "dot_args", "progressor_obj"),
      .noexport = c(
        "object", "p_sel", "search_path", "p_ref", "refmodel", "formulas",
        "projpred_var", "projpred_ws_aug", "projpred_formulas_no_random"
      ),
      .errorhandling = getOption("projpred.foreach_errorhandling", "stop"),
      .verbose = getOption("projpred.foreach_verbose", FALSE)
    ) %do_projpred% {
      mssgs_warns_capt <- capt_mssgs_warns(
        soutdmin <- do.call(
          sdivmin,
          c(list(formula = formula_s,
                 projpred_var = projpred_var_s,
                 projpred_formula_no_random = projpred_formula_no_random_s,
                 projpred_random = projpred_random),
            dot_args)
        )
      )
      if (!is.null(progressor_obj)) progressor_obj()
      return(nlist(soutdmin, mssgs_warns_capt))
    }
  }
  mssgs_warns_capts <- lapply(outdmin, "[[", "mssgs_warns_capt")
  outdmin <- lapply(outdmin, "[[", "soutdmin")
  mssgs_warns_capts <- lapply(mssgs_warns_capts, function(mssgs_warns_capt) {
    # Filter out some warnings.
    mssgs_warns_capt <- grep("Warning in [^:]*:$",
                             mssgs_warns_capt, value = TRUE, invert = TRUE)
    return(mssgs_warns_capt)
  })
  warn_prj_drawwise(mssgs_warns_capts, throw_warn = throw_warn_sdivmin)
  check_conv(outdmin, lengths(mssgs_warns_capts), do_check = do_check_conv)
  return(outdmin)
}

# Draw-wise divergence minimizer for divmin() that builds the shared setup
# (`make_shared`: lmer_shared_callback() or gam_shared_callback())
# from the arguments of its first call. The setup then uses the data of the
# draw-wise fits, also if divmin() changes the data before these fits. A call
# with other data, or a setup that does not apply, uses `sdivmin_plain`. The
# setup is built once, or once per task in parallel runs.
shared_setup_callback <- function(make_shared, formulas, sdivmin_plain) {
  # Force the arguments: divmin() assigns the result to `sdivmin`, so a
  # promise for `sdivmin_plain` would return this function (infinite
  # recursion):
  force(sdivmin_plain)
  force(make_shared)
  force(formulas)
  sdivmin_shared <- NULL
  data_shared <- NULL
  done <- FALSE
  return(function(formula, ...) {
    if (!done) {
      done <<- TRUE
      dot_args <- list(...)
      sdivmin_shared <<- make_shared(formulas, dot_args = dot_args)
      data_shared <<- dot_args[["data"]]
    }
    # identical() returns fast if both are the same object (the usual case):
    if (is.null(sdivmin_shared) ||
        !identical(list(...)[["data"]], data_shared)) {
      return(sdivmin_plain(formula = formula, ...))
    }
    return(sdivmin_shared(formula = formula, ...))
  })
}

# Additive submodels (mgcv::gam()) of all draws of one projection
# ---------------------------------------------------------------
#
# The draw-wise formulas differ only in the response, so all draws can share the
# setup of mgcv::gam() (model frame, smooth bases and penalties;
# mgcv::gam(fit = FALSE)). gam_shared_callback() fits each draw with
# mgcv::gam(G = ), the public interface for a prepared setup, with the response
# of the draw. It sets the elements of the fit that contain the response name
# (`formula`, `terms`, `pterms`, `model`, `pred.formula`) and `call` as in a fit
# by fit_gam_callback(), so the fits are equal (all.equal()) to those of
# fit_gam_callback().

# Draw-wise divergence minimizer for divmin() that uses one mgcv::gam() setup
# for all draws, or `NULL` if this does not apply.
gam_shared_callback <- function(formulas, dot_args) {
  if (!getOption("projpred.gam_shared", TRUE)) {
    return(NULL)
  }
  data <- dot_args$data
  resp_nms <- vapply(formulas, function(fml) {
    if (is.name(fml[[2]])) as.character(fml[[2]]) else NA_character_
  }, character(1))
  if (!is.data.frame(data) || anyNA(resp_nms) ||
      !all(resp_nms %in% names(data)) || anyNA(data[resp_nms])) {
    return(NULL)
  }
  # As in fit_gam_callback():
  gam_dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(mgcv::gam), methods::formalArgs(mgcv::gam.fit))
  )]
  if (any(c("G", "fit") %in% names(gam_dot_args))) {
    return(NULL)
  }
  g_setup <- tryCatch(
    do.call(mgcv::gam, c(list(formula = formulas[[1]]), gam_dot_args,
                         list(fit = FALSE))),
    error = function(e) NULL
  )
  if (is.null(g_setup) || length(g_setup$y) != nrow(data)) {
    return(NULL)
  }
  resp_old <- resp_nms[1]
  return(function(formula, ...) {
    resp_new <- as.character(formula[[2]])
    tryCatch({
      g_s <- g_setup
      g_s$y <- data[[resp_new]]
      g_s$mf[[resp_old]] <- data[[resp_new]]
      fit <- mgcv::gam(G = g_s)
      if (resp_new != resp_old) {
        for (el in c("formula", "terms", "pterms", "model", "pred.formula")) {
          if (!is.null(fit[[el]])) {
            fit[[el]] <- rename_sym(fit[[el]], resp_old, resp_new)
          }
        }
      }
      fit$call <- match.call(mgcv::gam, as.call(c(
        list(mgcv::gam), list(formula = formula), gam_dot_args
      )))
      fit
    }, error = function(e) fit_gam_callback(formula = formula, ...))
  })
}

# Linear predictors of `gam` fits (without group-level terms) of the same
# submodel for `newdata`, from one linear-predictor matrix (predict.gam() with
# `type = "lpmatrix"`) for all draws, instead of one predict.gam() call per
# draw. Returns `NULL` for other cases.
subprd_gam_direct <- function(fits, newdata) {
  if (is.null(newdata) || length(fits) < 2 ||
      !all(vapply(fits, function(fit) {
        inherits(fit, "gam") && !inherits(fit, "gamm4")
      }, logical(1)))) {
    return(NULL)
  }
  fit_1 <- fits[[1]]
  terms_noenv <- function(trms) {
    trms <- delete.response(trms)
    attr(trms, ".Environment") <- NULL
    trms
  }
  trms_1 <- terms_noenv(fit_1$pterms)
  smooth_1 <- fit_1$smooth
  if (!all(vapply(fits, function(fit) {
    identical(terms_noenv(fit$pterms), trms_1) &&
      identical(fit$smooth, smooth_1) &&
      identical(names(stats::coef(fit)), names(stats::coef(fit_1)))
  }, logical(1)))) {
    return(NULL)
  }
  if (length(attr(fit_1$terms, "offset")) ||
      length(attr(fit_1$pterms, "offset"))) {
    return(NULL)
  }
  x_lp <- stats::predict(fit_1, newdata = newdata, type = "lpmatrix")
  prd <- x_lp %*% do.call(cbind, lapply(fits, stats::coef))
  dimnames(prd) <- list(rownames(x_lp), NULL)
  return(prd)
}

# Gaussian multilevel submodels (lme4::lmer()) of all draws of one projection
# --------------------------------------------------------------------------
#
# lme4::lmer() is lme4::lFormula() (model frame, fixed-effects model matrix,
# random-effects structure), followed by mkLmerDevfun(), optimizeLmer(),
# checkConv() and mkMerMod(). The draw-wise formulas of one submodel differ only
# in the response, so all draws can share one lFormula() result.
# lmer_shared_fit() repeats the other steps of lmer() for one draw: it
# transcribes the body of lmer() in lme4 2.0-6 to 2.1-1 (R package lme4,
# GPL-2 | GPL-3). lme4's C++ objects change some vectors of the lFormula()
# result in place, so each draw gets fresh copies of them. The fits are equal
# (all.equal()) to those of lmer(). lmer_shared_ok() checks this once per R
# session on a small fixed problem, and the test suite checks it on CRAN. If
# the check fails (for example after a change in lme4), divmin() fits each draw
# with lmer().

# Replaces the name `old` by `new` in all language objects of `x`, including
# its attributes and names (used for the response of a shared model frame).
rename_sym <- function(x, old, new) {
  if (is.name(x)) {
    if (identical(x, as.name(old))) {
      return(as.name(new))
    }
    return(x)
  }
  if (is.call(x) || is.pairlist(x) || (is.list(x) && !is.data.frame(x))) {
    at <- attributes(x)
    for (i in seq_along(x)) {
      if (!is.null(x[[i]])) {
        x[[i]] <- rename_sym(x[[i]], old, new)
      }
    }
    attributes(x) <- at
  }
  at <- attributes(x)
  if (!is.null(at)) {
    for (nm in setdiff(names(at), c("names", "dim", "dimnames", "row.names",
                                    "class"))) {
      attr(x, nm) <- rename_sym(at[[nm]], old, new)
    }
    if (!is.null(at$names)) {
      names(x)[names(x) == old] <- new
    }
    if (!is.null(at$dimnames)) {
      attr(x, "dimnames") <- lapply(at$dimnames, function(dnms) {
        if (!is.null(dnms)) {
          dnms[dnms == old] <- new
        }
        dnms
      })
    }
  }
  return(x)
}

# Fits the lmer() submodel `formula` (response: a column of `data`) from the
# shared lFormula() result `lmod0` (of the same submodel for another draw).
# `lmer_args` are the arguments of the lmer() call that this fit replaces.
lmer_shared_fit <- function(lmod0, formula, data, control, lmer_args) {
  resp_old <- as.character(lmod0$formula[[2]])
  resp_new <- as.character(formula[[2]])
  lmod <- lmod0
  # Fresh copies of the vectors that lme4's C++ objects change in place:
  lmod$reTrms$Lambdat@x <- lmod0$reTrms$Lambdat@x + 0
  lmod$reTrms$theta <- lmod0$reTrms$theta + 0
  lmod$reTrms$Zt@x <- lmod0$reTrms$Zt@x + 0
  lmod$X <- lmod0$X + 0
  fr <- lmod$fr
  fr[[resp_old]] <- data[[resp_new]]
  if (resp_new != resp_old) {
    fr <- rename_sym(fr, resp_old, resp_new)
  }
  lmod$fr <- fr
  # The following steps are those of lme4::lmer():
  mcout <- match.call(lme4::lmer, as.call(c(list(lme4::lmer), lmer_args)))
  mcout$formula <- formula
  lmod$formula <- NULL
  devfun <- do.call(lme4::mkLmerDevfun,
                    c(lmod, list(start = NULL, verbose = 0L, control = control)))
  rho <- environment(devfun)
  nobs <- nrow(lmod$fr)
  npar <- length(rho$lower)
  calc.derivs <- control$calc.derivs %||%
    (nobs < control$checkConv$check.conv.nobsmax &&
       npar < control$checkConv$check.conv.nparmax)
  opt <- lme4::optimizeLmer(
    devfun, optimizer = control$optimizer,
    restart_edge = control$restart_edge, boundary.tol = control$boundary.tol,
    control = control$optCtrl, verbose = 0L, start = NULL,
    calc.derivs = calc.derivs,
    force.calc.derivs = control$calc.derivs %||% FALSE,
    use.last.params = control$use.last.params
  )
  # Same name as in lmer(), so that its warnings are the same:
  checkConv <- "lme4" %:::% "checkConv"
  cc <- checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
                  lbound = rho$lower, ubound = rho$upper, nobs = nobs,
                  ndim = npar)
  return(lme4::mkMerMod(rho, opt, lmod$reTrms, fr = lmod$fr, mc = mcout,
                        lme4conv = cc))
}

# Draw-wise divergence minimizer for divmin() that uses lmer_shared_fit() with
# the lFormula() result of the first draw, or `NULL` if this does not apply.
lmer_shared_callback <- function(formulas, dot_args) {
  if (!getOption("projpred.lmer_shared", TRUE) || !lmer_shared_ok()) {
    return(NULL)
  }
  data <- dot_args$data
  resp_nms <- vapply(formulas, function(fml) {
    if (is.name(fml[[2]])) as.character(fml[[2]]) else NA_character_
  }, character(1))
  if (!is.data.frame(data) || anyNA(resp_nms) ||
      !all(resp_nms %in% names(data)) || anyNA(data[resp_nms])) {
    return(NULL)
  }
  control <- dot_args$control %||% control_callback(dot_args$family)
  lmer_dot_args <- dot_args[intersect(names(dot_args),
                                      methods::formalArgs(lme4::lmer))]
  lmer_dot_args$control <- NULL
  lmod0 <- tryCatch(
    do.call(lme4::lFormula, c(list(formula = formulas[[1]], control = control),
                              lmer_dot_args)),
    error = function(e) NULL
  )
  if (is.null(lmod0) || nrow(lmod0$fr) != nrow(data)) {
    return(NULL)
  }
  force(lmod0)
  return(function(formula, ...) {
    lmer_args <- c(list(formula = formula, control = control), lmer_dot_args)
    tryCatch(
      lmer_shared_fit(lmod0, formula = formula, data = data, control = control,
                      lmer_args = lmer_args),
      error = function(e) fit_glmer_callback(formula = formula, ...)
    )
  })
}

# Whether lmer_shared_fit() gives the fits of lme4::lmer() with the installed
# lme4 (checked once per R session on a small fixed problem, without random
# numbers).
lmer_shared_state <- new.env(parent = emptyenv())
lmer_shared_ok <- function() {
  if (is.null(lmer_shared_state$ok)) {
    lmer_shared_state$ok <- isTRUE(tryCatch(lmer_shared_check(),
                                            error = function(e) FALSE))
  }
  return(lmer_shared_state$ok)
}
lmer_shared_check <- function() {
  n <- 60L
  idx <- seq_len(n)
  dat <- data.frame(x1 = sin(idx), x2 = cos(1.7 * idx),
                    g = factor(paste0("g", idx %% 6L)),
                    h = factor(paste0("h", idx %% 4L)))
  u <- c(0.8, -0.5, 0.3, -0.9, 0.6, -0.2)
  for (s in 1:3) {
    dat[[paste0(".y.", s)]] <- dat$x1 + 0.5 * dat$x2 +
      u[as.integer(dat$g)] * (1 + 0.3 * s) + 0.4 * sin(s * idx) +
      0.2 * dat$x1 * u[as.integer(dat$g)]
  }
  wts <- 1 + 0.5 * cos(idx)
  control <- lme4::lmerControl()
  fmls <- lapply(1:3, function(s) {
    as.formula(paste0(".y.", s, " ~ x1 + x2 + (1 + x1 | g) + (1 | h)"))
  })
  lmod0 <- lme4::lFormula(fmls[[1]], data = dat, weights = wts,
                          control = control)
  same <- logical(3)
  for (s in 1:3) {
    lmer_args <- list(formula = fmls[[s]], control = control, data = dat,
                      weights = wts)
    suppressMessages(suppressWarnings({
      fit_ref <- do.call(lme4::lmer, lmer_args)
      fit_shared <- lmer_shared_fit(lmod0, formula = fmls[[s]], data = dat,
                                    control = control, lmer_args = lmer_args)
    }))
    same[s] <- isTRUE(all.equal(fit_ref, fit_shared, tolerance = 0)) &&
      identical(lme4::getME(fit_ref, "theta"),
                lme4::getME(fit_shared, "theta"))
  }
  return(all(same))
}

# Use projpred's own implementation to fit non-multilevel non-additive
# submodels:
fit_glm_ridge_callback <- function(formula, data,
                                   projpred_var = matrix(nrow = nrow(data)),
                                   regul = 1e-4, thresh_conv = 1e-7, ...) {
  # Preparations:
  fr <- model.frame(formula, data = data, drop.unused.levels = TRUE)
  da_classes <- attr(attr(fr, "terms"), "dataClasses")
  nms_chr_fac <- names(da_classes)[da_classes %in% c("character", "factor")]
  resp_nm <- all.vars(attr(fr, "terms"))[attr(attr(fr, "terms"), "response")]
  nms_chr_fac <- setdiff(nms_chr_fac, resp_nm)
  if (length(nms_chr_fac) > 0) {
    xlvls <- lapply(setNames(nm = nms_chr_fac), function(nm_chr_fac) {
      levels(as.factor(fr[[nm_chr_fac]]))
    })
  } else {
    xlvls <- NULL
  }
  # TODO: In the following model.matrix() call, allow user-specified contrasts
  # to be passed to argument `contrasts.arg`. The `contrasts.arg` default
  # (`NULL`) uses `options("contrasts")` internally, but it might be more
  # convenient to let users specify contrasts directly. At that occasion,
  # contrasts should also be tested thoroughly (not done until now).
  x <- model.matrix(formula, data = fr)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  y <- model.response(fr)
  # Exclude arguments from `...` which cannot be passed to glm_ridge():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    methods::formalArgs(glm_ridge)
  )]
  # Call the submodel fitter:
  out_capt <- utils::capture.output(
    fit <- do.call(glm_ridge, c(
      list(x = x, y = y, obsvar = projpred_var, lambda = regul,
           thresh = thresh_conv),
      dot_args
    ))
  )
  out_capt <- unique(grep("[Ww]arning|bug", out_capt, value = TRUE))
  if (length(out_capt) > 0) {
    warning(paste(out_capt, collapse = "\n"))
  }
  # Post-processing:
  rownames(fit$beta) <- colnames(x)
  sub <- nlist(alpha = fit$beta0, beta = fit$beta, w = fit$w, formula, x, y,
               xlvls)
  class(sub) <- "subfit"
  return(sub)
}

# Alternative to fit_glm_ridge_callback() (may be used via global option
# `projpred.glm_fitter`):
fit_glm_callback <- function(formula, family, ...) {
  if (family$family == "gaussian" && family$link == "identity" &&
      getOption("projpred.gaussian_not_as_generalized", TRUE)) {
    # Exclude arguments from `...` which cannot be passed to stats::lm():
    dot_args <- list(...)
    dot_args <- dot_args[intersect(
      names(dot_args),
      c(methods::formalArgs(stats::lm),
        methods::formalArgs(stats::lm.fit),
        methods::formalArgs(stats::lm.wfit))
    )]
    # Call the submodel fitter:
    return(do.call(stats::lm, c(list(formula = formula), dot_args)))
  } else {
    # Exclude arguments from `...` which cannot be passed to stats::glm():
    dot_args <- list(...)
    dot_args <- dot_args[intersect(
      names(dot_args),
      c(methods::formalArgs(stats::glm),
        methods::formalArgs(stats::glm.control))
    )]
    # Call the submodel fitter:
    return(do.call(stats::glm, c(list(formula = formula, family = family),
                                 dot_args)))
  }
}

# Use package "mgcv" to fit additive non-multilevel submodels:
fit_gam_callback <- function(formula, ...) {
  # Exclude arguments from `...` which cannot be passed to mgcv::gam():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(mgcv::gam),
      methods::formalArgs(mgcv::gam.fit))
  )]
  # Call the submodel fitter:
  return(do.call(mgcv::gam, c(list(formula = formula), dot_args)))
}

# Use package "gamm4" to fit additive multilevel submodels:
fit_gamm_callback <- function(formula, projpred_formula_no_random,
                              projpred_random, data, family,
                              control = control_callback(family), ...) {
  # Exclude arguments from `...` which cannot be passed to gamm4::gamm4():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(gamm4::gamm4),
      methods::formalArgs(lme4::lFormula),
      methods::formalArgs(lme4::glFormula))
  )]
  # Call the submodel fitter:
  fit <- tryCatch({
    do.call(gamm4::gamm4, c(
      list(formula = projpred_formula_no_random, random = projpred_random,
           data = data, family = family, control = control),
      dot_args
    ))
  }, error = function(e) {
    if (grepl("not positive definite", as.character(e))) {
      if ("optimx" %in% control$optimizer &&
          length(control$optCtrl$method) > 0 &&
          control$optCtrl$method == "nlminb") {
        stop("Encountering the `not positive definite` error while running ",
             "the lme4 fitting procedure, but cannot fix this automatically ",
             "anymore. You will probably have to tweak gamm4 tuning ",
             "parameters manually (via `...`).")
      }
      if (!requireNamespace("optimx", quietly = TRUE)) {
        stop("Trying to fix an lme4 error with the help of an optimx ",
             "optimizer, but cannot find the optimx package. Please install ",
             "it to be able to try this fix.")
      }
      return(fit_gamm_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        data = data,
        family = family,
        control = control_callback(family,
                                   optimizer = "optimx",
                                   optCtrl = list(method = "nlminb")),
        ...
      ))
    } else if (grepl("pwrssUpdate", as.character(e))) {
      # lme4's PIRLS could not reduce the deviance any further, so relax the
      # convergence tolerance for the PIRLS step by a factor of 10 and retry
      # (analogous to the automatic fix in fit_glmer_callback()):
      tolPwrss_new <- 10 * (control$tolPwrss %||%
                            control_callback(family)$tolPwrss)
      if (tolPwrss_new > 1e-4) {
        stop("Encountering a pwrssUpdate error while running the gamm4 ",
             "fitting procedure, but cannot fix this automatically anymore. ",
             "You will probably have to tweak lme4 tuning parameters manually ",
             "(via `...`).")
      }
      return(fit_gamm_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        data = data,
        family = family,
        control = control_callback(family, tolPwrss = tolPwrss_new),
        ...
      ))
    } else {
      stop(e)
    }
  })

  fit$random <- projpred_random
  fit$formula <- projpred_formula_no_random
  class(fit) <- c("gamm4")
  return(fit)
}

# Use package "lme4" to fit multilevel submodels:
fit_glmer_callback <- function(formula, projpred_formula_no_random,
                               projpred_random, family,
                               control = control_callback(family), ...) {
  tryCatch({
    if (getOption("projpred.PQL", FALSE)) {
      # Exclude arguments from `...` which cannot be passed to MASS::glmmPQL():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(MASS::glmmPQL)
      )]
      # Strip parentheses from group-level terms:
      random_trms <- labels(terms(projpred_random))
      random_fmls <- lapply(random_trms, function(random_trm) {
        as.formula(paste("~", random_trm))
      })
      if (length(random_fmls) == 1) {
        random_fmls <- random_fmls[[1]]
      } else if (length(random_fmls) > 1) {
        # Use the workaround(s) from <https://stackoverflow.com/questions/
        # 36643713/how-to-specify-different-random-effects-in-nlme-vs-lme4/
        # 38805602#38805602>:
        random_trms_nogrp <- lapply(random_trms, function(random_trm) {
          random_fml_nogrp <- as.formula(
            paste("~", sub("[[:blank:]]*\\|.*$", "", random_trm))
          )
          return(labels(terms(random_fml_nogrp)))
        })
        ### TODO (glmmPQL): Do this via adding argument `data` explicitly:
        stopifnot(!is.null(dot_args$data))
        if ("projpred_internal_dummy1s_PQL" %in% names(dot_args$data)) {
          stop("Need to write to column `projpred_internal_dummy1s_PQL` of ",
               "`data`, but that column already exists. Please rename this ",
               "column in `data` and try again.")
        }
        dot_args$data$projpred_internal_dummy1s_PQL <- factor(1)
        ###
        list_pdIdent <- lapply(random_trms, function(random_trm) {
          return(nlme::pdIdent(as.formula(
            paste("~", sub("^.*\\|[[:blank:]]*", "", random_trm), "- 1")
          )))
        })
        idxs_multi_trms <- which(lengths(random_trms_nogrp) > 0)
        if (length(idxs_multi_trms) > 0) {
          warning("In order to be able to use MASS::glmmPQL(), we have to use ",
                  "a workaround by which the random intercepts and random ",
                  "slopes will be assumed to have a correlation of zero.")
          list_pdIdent_add <- lapply(idxs_multi_trms, function(idx_multi_trms) {
            random_grp <- sub("^.*\\|[[:blank:]]*", "",
                              random_trms[[idx_multi_trms]])
            random_IAs <- paste0(random_trms_nogrp[[idx_multi_trms]], ":",
                                 random_grp)
            return(lapply(random_IAs, function(random_IA) {
              nlme::pdIdent(as.formula(paste("~", random_IA)))
            }))
          })
          list_pdIdent_add <- unlist(list_pdIdent_add, recursive = FALSE)
          list_pdIdent <- c(list_pdIdent, list_pdIdent_add)
        }
        random_fmls <- list(
          projpred_internal_dummy1s_PQL = nlme::pdBlocked(list_pdIdent)
        )
      } else if (length(random_fmls) == 0) {
        stop("Unexpected length of `random_fmls`.")
      }
      # Call the submodel fitter:
      return(do.call(MASS::glmmPQL, c(
        list(fixed = projpred_formula_no_random, random = random_fmls,
             family = family, control = control),
        dot_args
      )))
    } else if (family$family == "gaussian" && family$link == "identity" &&
               getOption("projpred.gaussian_not_as_generalized", TRUE)) {
      # Exclude arguments from `...` which cannot be passed to lme4::lmer():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(lme4::lmer)
      )]
      # Call the submodel fitter:
      return(do.call(lme4::lmer, c(list(formula = formula, control = control),
                                   dot_args)))
    } else {
      # Exclude arguments from `...` which cannot be passed to lme4::glmer():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(lme4::glmer)
      )]
      # Call the submodel fitter:
      return(do.call(lme4::glmer, c(list(formula = formula, family = family,
                                         control = control),
                                    dot_args)))
    }
  }, error = function(e) {
    if (grepl("No random effects", as.character(e))) {
      # This case should not occur anymore (because divmin() should pick the
      # correct submodel fitter based on the submodel's formula), but leave it
      # here in case user-specified divergence minimizers make use of it.
      glm_fitter <- get(getOption("projpred.glm_fitter",
                                  "fit_glm_ridge_callback"),
                        mode = "function")
      return(glm_fitter(
        formula = formula, family = family, ...
      ))
    } else if (grepl("not positive definite", as.character(e))) {
      if ("optimx" %in% control$optimizer &&
          length(control$optCtrl$method) > 0 &&
          control$optCtrl$method == "nlminb") {
        stop("Encountering the `not positive definite` error while running ",
             "the lme4 fitting procedure, but cannot fix this automatically ",
             "anymore. You will probably have to tweak lme4 tuning parameters ",
             "manually (via `...`).")
      }
      if (!requireNamespace("optimx", quietly = TRUE)) {
        stop("Trying to fix an lme4 error with the help of an optimx ",
             "optimizer, but cannot find the optimx package. Please install ",
             "it to be able to try this fix.")
      }
      return(fit_glmer_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        family = family,
        control = control_callback(family,
                                   optimizer = "optimx",
                                   optCtrl = list(method = "nlminb")),
        ...
      ))
    } else if (grepl("PIRLS", as.character(e))) {
      if (length(dot_args$nAGQ) > 0) {
        nAGQ_new <- dot_args$nAGQ + 1L
      } else {
        nAGQ_new <- 20L
      }
      if (nAGQ_new > 30L) {
        stop("Encountering a PIRLS error while running the lme4 fitting ",
             "procedure, but cannot fix this automatically anymore. You will ",
             "probably have to tweak lme4 tuning parameters manually (via ",
             "`...`).")
      }
      return(fit_glmer_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        family = family,
        control = control,
        nAGQ = nAGQ_new,
        ...
      ))
    } else if (grepl("pwrssUpdate did not converge in \\(maxit\\) iterations",
                     as.character(e))) {
      tolPwrss_new <- 10 * control$tolPwrss
      optCtrl_new <- list()
      if (length(control$optimizer) == 0) {
        stop("Unexpected length of `control$optimizer`. Please notify the ",
             "package maintainer.")
      }
      if (any(control$optimizer %in% c("Nelder_Mead", "bobyqa"))) {
        if (length(control$optCtrl$maxfun) > 0) {
          maxfun_new <- 10 * control$optCtrl$maxfun
        } else {
          maxfun_new <- 1e4
        }
        optCtrl_new <- c(optCtrl_new, list(maxfun = maxfun_new))
      }
      if (any(!control$optimizer %in% c("Nelder_Mead", "bobyqa"))) {
        if (length(control$optCtrl$maxit) > 0) {
          maxit_new <- 10 * control$optCtrl$maxit
        } else {
          maxit_new <- 1e4
        }
        optCtrl_new <- c(optCtrl_new, list(maxit = maxit_new))
      }
      if (tolPwrss_new > 1e-4 &&
          (optCtrl_new$maxfun %||% -Inf > 1e7 ||
           optCtrl_new$maxit %||% -Inf > 1e7)) {
        stop("Encountering the ",
             "`pwrssUpdate did not converge in (maxit) iterations` error ",
             "while running the lme4 fitting procedure, but cannot fix this ",
             "automatically anymore. You will probably have to tweak lme4 ",
             "tuning parameters manually (via `...`).")
      }
      return(fit_glmer_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        family = family,
        control = control_callback(family, tolPwrss = tolPwrss_new,
                                   optCtrl = optCtrl_new),
        ...
      ))
    } else if (getOption("projpred.PQL", FALSE) &&
               grepl("iteration limit reached without convergence",
                     as.character(e))) {
      if (length(control$msMaxIter) > 0 && control$msMaxIter >= 100) {
        stop("Encountering the `iteration limit reached without convergence` ",
             "error while running the MASS::glmmPQL() fitting procedure, but ",
             "cannot fix this automatically anymore. You will probably have ",
             "to tweak MASS::glmmPQL() tuning parameters manually (via `...`).")
      }
      return(fit_glmer_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        family = family,
        control = control_callback(msMaxIter = 100),
        ...
      ))
    } else if (getOption("projpred.PQL", FALSE) &&
               grepl("false convergence", as.character(e))) {
      if (length(control$niterEM) > 0 && control$niterEM >= 50) {
        stop("Encountering the `false convergence` ",
             "error while running the MASS::glmmPQL() fitting procedure, but ",
             "cannot fix this automatically anymore. You will probably have ",
             "to tweak MASS::glmmPQL() tuning parameters manually (via `...`).")
      }
      return(fit_glmer_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        family = family,
        control = control_callback(niterEM = 50),
        ...
      ))
    } else if (getOption("projpred.PQL", FALSE) &&
               grepl("fewer observations than random effects",
                     as.character(e))) {
      if (length(control$allow.n.lt.q) > 0 && isTRUE(control$allow.n.lt.q)) {
        stop("Encountering the `fewer observations than random effects` ",
             "error while running the MASS::glmmPQL() fitting procedure, but ",
             "cannot fix this automatically anymore. You will probably have ",
             "to tweak MASS::glmmPQL() tuning parameters manually (via `...`).")
      }
      return(fit_glmer_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        family = family,
        control = control_callback(allow.n.lt.q = TRUE),
        ...
      ))
    } else {
      stop(e)
    }
  })
}

# Helper function for fit_glmer_callback() and fit_gamm_callback() to get the
# appropriate control options depending on the family:
control_callback <- function(family, ...) {
  if (getOption("projpred.PQL", FALSE)) {
    return(nlme::lmeControl(...))
  } else if (family$family == "gaussian" && family$link == "identity" &&
             getOption("projpred.gaussian_not_as_generalized", TRUE)) {
    return(lme4::lmerControl(...))
  } else {
    return(lme4::glmerControl(...))
  }
}

## Augmented-data projection ----------------------------------------------

# Needed to avoid a NOTE in `R CMD check`:
if (getRversion() >= package_version("2.15.1")) {
  utils::globalVariables("projpred_w_aug_s")
  utils::globalVariables("projpred_internal_w_aug")
}

divmin_augdat <- function(
    formula,
    data,
    family,
    weights,
    projpred_var,
    projpred_ws_aug,
    verbose_divmin,
    throw_warn_sdivmin = getOption("projpred.warn_proj_drawwise", TRUE),
    do_check_conv = getOption("projpred.check_convergence", TRUE),
    ...
) {
  if (missing(throw_warn_sdivmin) &&
      is.null(getOption("projpred.warn_proj_drawwise")) &&
      !is.null(getOption("projpred.warn_prj_drawwise"))) {
    warning(
      "Global option `projpred.warn_prj_drawwise` is deprecated. Please use ",
      "global option `projpred.warn_proj_drawwise` instead. Now using the ",
      "value from global option `projpred.warn_prj_drawwise`."
    )
    throw_warn_sdivmin <- getOption("projpred.warn_prj_drawwise")
  }
  if (missing(do_check_conv) &&
      is.null(getOption("projpred.check_convergence")) &&
      !is.null(getOption("projpred.check_conv"))) {
    warning(
      "Global option `projpred.check_conv` is deprecated. Please use ",
      "global option `projpred.check_convergence` instead. Now using the ",
      "value from global option `projpred.check_conv`."
    )
    do_check_conv <- getOption("projpred.check_conv")
  }
  trms_all <- extract_terms_response(formula)
  has_grp <- length(trms_all$group_terms) > 0
  projpred_formula_no_random <- NA
  projpred_random <- NA
  if (family$family == "binomial") {
    # Switch back to the traditionally extended binomial family (so that the
    # original link and inverse-link functions are used, for example):
    family <- extend_family(
      get(family$family, mode = "function")(link = family$link)
    )
    if (!has_grp) {
      sdivmin <- get(getOption("projpred.glm_fitter", "fit_glm_ridge_callback"),
                     mode = "function")
      if (getOption("projpred.glm_fitter", "fit_glm_ridge_callback") ==
          "fit_glm_ridge_callback") {
        # glm_ridge() cannot handle `factor` responses:
        response_name <- extract_terms_response(formula)$response
        stopifnot(is.factor(data[[response_name]]))
        data[[response_name]] <- as.integer(data[[response_name]]) - 1L
      }
    } else {
      sdivmin <- fit_glmer_callback
      if (getOption("projpred.PQL", FALSE)) {
        formula_random <- split_formula_random_gamm4(formula)
        projpred_formula_no_random <- formula_random$formula
        projpred_random <- formula_random$random
      }
    }
  } else if (family$family %in% c("cumulative", "cumulative_rstanarm")) {
    if (!has_grp) {
      sdivmin <- fit_cumul
    } else {
      sdivmin <- fit_cumul_mlvl
    }
  } else if (family$family == "categorical") {
    if (!has_grp) {
      sdivmin <- fit_categ
    } else {
      sdivmin <- fit_categ_mlvl
      formula_random <- split_formula_random_gamm4(formula)
      projpred_formula_no_random <- formula_random$formula
      projpred_random <- formula_random$random
    }
  } else {
    stop("Family `", family$family, "` is not supported by divmin_augdat().")
  }

  # Coerce augmented-data matrices to ordinary matrices and augmented-length
  # vectors to ordinary vectors so that external packages may subset as usual:
  projpred_ws_aug <- unclass(projpred_ws_aug)
  attr(projpred_ws_aug, "ndiscrete") <- NULL

  if (is.null(getOption("projpred.parallel_proj_trigger")) &&
      !is.null(getOption("projpred.prll_prj_trigger"))) {
    warning(
      "Global option `projpred.prll_prj_trigger` is deprecated. Please use ",
      "global option `projpred.parallel_proj_trigger` instead. Now using the ",
      "value from global option `projpred.prll_prj_trigger`."
    )
    prll_prj_trigger <- getOption("projpred.prll_prj_trigger")
  } else {
    prll_prj_trigger <- getOption("projpred.parallel_proj_trigger", Inf)
  }
  if (ncol(projpred_ws_aug) < prll_prj_trigger) {
    # Sequential case. Actually, we could simply use ``%do_projpred%` <-
    # foreach::`%do%`` here and then proceed as in the parallel case, but that
    # would require adding more "hard" dependencies (because packages 'foreach'
    # and 'iterators' would have to be moved from `Suggests:` to `Imports:`).
    if (verbose_divmin) {
      pb <- utils::txtProgressBar(min = 0, max = ncol(projpred_ws_aug),
                                  style = 3, initial = 0, file = stderr())
      on.exit(close(pb))
    }
    outdmin <- lapply(seq_len(ncol(projpred_ws_aug)), function(s) {
      if (verbose_divmin) {
        on.exit(utils::setTxtProgressBar(pb, s))
      }
      mssgs_warns_capt <- capt_mssgs_warns(
        soutdmin <- sdivmin(
          formula = formula,
          data = data,
          family = family,
          weights = projpred_ws_aug[, s],
          projpred_formula_no_random = projpred_formula_no_random,
          projpred_random = projpred_random,
          ...
        )
      )
      return(nlist(soutdmin, mssgs_warns_capt))
    })
  } else {
    # Parallel case.
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Please install the 'foreach' package.")
    }
    if (!requireNamespace("iterators", quietly = TRUE)) {
      stop("Please install the 'iterators' package.")
    }
    if (verbose_divmin && use_progressr()) {
      progressor_obj <- progressr::progressor(ncol(projpred_ws_aug))
    } else {
      progressor_obj <- NULL
    }
    dot_args <- list(...)
    `%do_projpred%` <- foreach::`%dopar%`
    outdmin <- foreach::foreach(
      projpred_w_aug_s = iterators::iter(projpred_ws_aug, by = "column"),
      .packages = c("projpred"),
      .export = c(
        "sdivmin", "formula", "data", "family", "projpred_formula_no_random",
        "projpred_random", "dot_args", "progressor_obj"
      ),
      .noexport = c(
        "object", "p_sel", "search_path", "p_ref", "refmodel", "projpred_var",
        "projpred_ws_aug", "linkobjs"
      ),
      .errorhandling = getOption("projpred.foreach_errorhandling", "stop"),
      .verbose = getOption("projpred.foreach_verbose", FALSE)
    ) %do_projpred% {
      mssgs_warns_capt <- capt_mssgs_warns(
        soutdmin <- do.call(
          sdivmin,
          c(list(formula = formula,
                 data = data,
                 family = family,
                 weights = as.vector(projpred_w_aug_s),
                 projpred_formula_no_random = projpred_formula_no_random,
                 projpred_random = projpred_random),
            dot_args)
        )
      )
      if (!is.null(progressor_obj)) progressor_obj()
      return(nlist(soutdmin, mssgs_warns_capt))
    }
  }
  mssgs_warns_capts <- lapply(outdmin, "[[", "mssgs_warns_capt")
  outdmin <- lapply(outdmin, "[[", "soutdmin")
  mssgs_warns_capts <- lapply(mssgs_warns_capts, function(mssgs_warns_capt) {
    # Filter out some warnings.
    mssgs_warns_capt <- grep("Warning in [^:]*:$",
                             mssgs_warns_capt, value = TRUE, invert = TRUE)
    # For MASS::polr():
    mssgs_warns_capt <- grep("non-integer #successes in a binomial glm!$",
                             mssgs_warns_capt, value = TRUE, invert = TRUE)
    # For ordinal::clmm():
    mssgs_warns_capt <- grep(paste("Using formula\\(x\\) is deprecated when x",
                                   "is a character vector of length > 1\\.$"),
                             mssgs_warns_capt, value = TRUE, invert = TRUE)
    # For ordinal::clmm():
    mssgs_warns_capt <- grep(
      "Consider formula\\(paste\\(x, collapse = .*\\)\\) instead\\.$",
      mssgs_warns_capt, value = TRUE, invert = TRUE
    )
    return(mssgs_warns_capt)
  })
  warn_prj_drawwise(mssgs_warns_capts, throw_warn = throw_warn_sdivmin)
  check_conv(outdmin, lengths(mssgs_warns_capts), do_check = do_check_conv)
  return(outdmin)
}

# Use MASS::polr() to fit submodels for the brms::cumulative() family:
fit_cumul <- function(formula, data, family, weights, ...) {
  # Argument `weights` of MASS::polr() expects a column name (as a symbol):
  if ("projpred_internal_w_aug" %in% names(data)) {
    stop("Need to write to column `projpred_internal_w_aug` of `data`, but ",
         "that column already exists. Please rename this column in `data` and ",
         "try again.")
  }
  stopifnot(length(weights) == nrow(data))
  data$projpred_internal_w_aug <- weights
  # Exclude arguments from `...` which cannot be passed to MASS::polr():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(MASS::polr),
      methods::formalArgs(stats::optim))
  )]
  # Adapt `link`:
  link_nm <- family$link
  if (link_nm == "logit") {
    link_nm <- "logistic"
  } else if (link_nm == "probit_approx") {
    link_nm <- "probit"
  }
  # Call the submodel fitter:
  fitobj <- try(do.call(MASS::polr, c(
    list(formula = formula,
         data = data,
         weights = quote(projpred_internal_w_aug),
         model = FALSE,
         method = link_nm),
    dot_args
  )), silent = TRUE)
  if (inherits(fitobj, "try-error") &&
      grepl(paste("initial value in 'vmmin' is not finite",
                  "attempt to find suitable starting values failed",
                  sep = "|"),
            attr(fitobj, "condition")$message)) {
    # Try to fix this automatically by specifying `start` values (the check for
    # the intercept being present is just performed to show that *by
    # construction*, we expect an intercept here; later code might be general
    # enough to deal with a missing intercept, but that's not the point).
    stopifnot(attr(terms(formula), "intercept") == 1)
    ncoefs <- ncol(model.matrix(formula, data = data)) -
      attr(terms(formula), "intercept")
    start_coefs <- rep(0, ncoefs)
    ncats <- length(unique(eval_lhs(formula = formula, data = data)))
    nthres <- ncats - 1L
    # Start with thresholds which imply equal probabilities for the response
    # categories:
    start_thres <- linkfun_raw(seq_len(nthres) / ncats, link_nm = link_nm)
    fitobj <- do.call(MASS::polr, c(
      list(formula = formula,
           data = data,
           weights = quote(projpred_internal_w_aug),
           model = FALSE,
           method = link_nm,
           start = c(start_coefs, start_thres)),
      dot_args
    ))
  } else if (inherits(fitobj, "try-error")) {
    stop(attr(fitobj, "condition")$message)
  }
  return(fitobj)
}

# Use ordinal::clmm() to fit multilevel submodels for the brms::cumulative()
# family:
fit_cumul_mlvl <- function(formula, data, family, weights, ...) {
  # Argument `weights` of ordinal::clmm() expects a column name (as a symbol):
  if ("projpred_internal_w_aug" %in% names(data)) {
    stop("Need to write to column `projpred_internal_w_aug` of `data`, but ",
         "that column already exists. Please rename this column in `data` and ",
         "try again.")
  }
  stopifnot(length(weights) == nrow(data))
  data$projpred_internal_w_aug <- weights
  # Exclude arguments from `...` which cannot be passed to ordinal::clmm():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(ordinal::clmm),
      methods::formalArgs(ordinal::clm.control),
      if (requireNamespace("ucminf", quietly = TRUE)) {
        methods::formalArgs(ucminf::ucminf)
      } else {
        character()
      },
      methods::formalArgs(stats::nlminb),
      methods::formalArgs(stats::optim))
  )]
  # Since ordinal::clmm() sometimes seems to have issues with formulas lacking
  # an intercept, add it here (note that this is a quite "hacky" solution;
  # perhaps reformulate() could be used instead):
  stopifnot(attr(terms(formula), "intercept") == 1)
  formula[[3]] <- str2lang(paste0("1 + ", as.character(formula)[[3]]))
  # Adapt `link`:
  link_nm <- family$link
  if (link_nm == "probit_approx") {
    link_nm <- "probit"
  }
  # Call the submodel fitter:
  fitobj <- do.call(ordinal::clmm, c(
    list(formula = formula,
         data = data,
         weights = quote(projpred_internal_w_aug),
         contrasts = NULL,
         Hess = FALSE,
         model = FALSE,
         link = link_nm),
    dot_args
  ))
  # Needed for the ordinal:::predict.clm() workaround (the value `"negative"` is
  # the default, see `?ordinal::clm.control`):
  fitobj$control$sign.location <- "negative"
  return(fitobj)
}

# Use nnet::multinom() to fit submodels for the brms::categorical() family:
fit_categ <- function(formula, data, family, weights, ...) {
  # Exclude arguments from `...` which cannot be passed to nnet::multinom():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(nnet::multinom),
      methods::formalArgs(nnet::nnet))
  )]
  # Adapt `link`:
  link_nm <- family$link
  if (link_nm != "logit") {
    stop("For the brms::categorical() family, projpred only supports the ",
         "logit link.")
  }
  # Handle offsets (but note that currently, offsets in the submodel fitter are
  # unexpected; this case could only occur in `test_div_minimizer.R`):
  trms <- terms(formula, data = data)
  offs_attr <- attr(trms, "offset")
  if (length(offs_attr)) {
    mfr <- model.frame(formula = formula, data = data)
    offs_obj <- model.offset(mfr)
    if (is.null(dim(offs_obj))) {
      # Recycle the vector of offsets to a matrix with number of columns equal
      # to the number of response categories:
      resp_vec <- model.response(mfr)
      ncats <- length(unique(resp_vec))
      offs_mat <- matrix(offs_obj, nrow = nrow(data), ncol = ncats)
      offs_nm <- as.character(attr(trms, "variables")[[offs_attr + 1L]][[2]])
      data[[offs_nm]] <- offs_mat
    }
  }
  # Call the submodel fitter:
  out_capt <- utils::capture.output(
    fitobj <- do.call(nnet::multinom, c(
      list(formula = formula,
           data = data,
           weights = weights),
      dot_args
    ))
  )
  if (utils::tail(out_capt, 1) != "converged") {
    warning(
      "Could not find the string \"converged\" at the end of the `stdout` ",
      "output from the nnet::multinom() submodel fit, so perhaps this fitting ",
      "procedure did not converge."
    )
  }
  return(fitobj)
}

# Use mclogit::mblogit() to fit multilevel submodels for the brms::categorical()
# family:
fit_categ_mlvl <- function(formula, projpred_formula_no_random,
                           projpred_random, data, family, weights, ...) {
  # Exclude arguments from `...` which cannot be passed to mclogit::mblogit():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(mclogit::mblogit),
      methods::formalArgs(mclogit::mmclogit.control))
  )]
  # Adapt `link`:
  link_nm <- family$link
  if (link_nm != "logit") {
    stop("For the brms::categorical() family, projpred only supports the ",
         "logit link.")
  }
  # Handle offsets (but note that currently, offsets in the submodel fitter are
  # unexpected; this case could only occur in `test_div_minimizer.R`):
  trms <- terms(projpred_formula_no_random, data = data)
  offs_attr <- attr(trms, "offset")
  if (length(offs_attr)) {
    # It is not clear whether mclogit::mblogit() supports offsets correctly (but
    # currently, offsets are unexpected here anyway):
    warning("Offsets for a multilevel submodel of a brms::categorical() ",
            "reference model are currently experimental.")
    mfr <- model.frame(formula = projpred_formula_no_random, data = data)
    offs_obj <- model.offset(mfr)
    if (is.null(dim(offs_obj))) {
      # Recycle the vector of offsets to a matrix with number of columns equal
      # to the number of response categories:
      resp_vec <- model.response(mfr)
      ncats <- length(unique(resp_vec))
      offs_mat <- matrix(offs_obj, nrow = nrow(data), ncol = ncats)
      offs_nm <- as.character(attr(trms, "variables")[[offs_attr + 1L]][[2]])
      data[[offs_nm]] <- offs_mat
    }
  }
  # Strip parentheses from group-level terms:
  random_trms <- labels(terms(projpred_random))
  if (utils::packageVersion("mclogit") < "0.9" && length(random_trms) > 1) {
    stop("'mclogit' versions < 0.9 can only handle a single group-level term.")
  }
  random_fmls <- lapply(random_trms, function(random_trm) {
    as.formula(paste("~", random_trm))
  })
  if (length(random_fmls) == 1) {
    random_fmls <- random_fmls[[1]]
  }
  # Call the submodel fitter:
  out_capt <- utils::capture.output(
    fitobj <- do.call(mclogit::mblogit, c(
      list(formula = projpred_formula_no_random,
           data = data,
           random = random_fmls,
           weights = weights,
           model = FALSE,
           x = FALSE,
           y = FALSE,
           method = "PQL"),
      dot_args
    ))
  )
  if (utils::tail(out_capt, 1) != "converged") {
    warning(
      "Could not find the string \"converged\" at the end of the `stdout` ",
      "output from the mclogit::mblogit() submodel fit, so perhaps this ",
      "fitting procedure did not converge."
    )
  }
  return(fitobj)
}

# Convergence issues ------------------------------------------------------

# Throw unique messages and warnings from a list of messages and warnings
# retrieved during draw-wise projections:
warn_prj_drawwise <- function(mssgs_warns_capts, throw_warn = TRUE) {
  if (!throw_warn) return()
  mssgs_warns_capts_unq <- unique(unlist(mssgs_warns_capts))
  if (length(mssgs_warns_capts_unq) > 0) {
    warning(paste(
      c(paste0("The following messages and/or warnings have been thrown by ",
               "the current submodel fitter (i.e., the current draw-wise ",
               "divergence minimizer):"),
        "---", mssgs_warns_capts_unq, "---"),
      collapse = "\n"
    ))
  }
  return()
}

# Check the convergence of the submodel fits from a whole `outdmin` object, also
# taking into account whether messages and warnings were thrown (indicated by
# argument `lengths_mssgs_warns` which must be of the same length as `outdmin`):
check_conv <- function(outdmin, lengths_mssgs_warns, do_check = TRUE) {
  if (!do_check) return()
  is_conv <- tryCatch(
    unlist(lapply(outdmin, check_conv_s)),
    error = function(e) {
      warning("The draw-wise convergence checker errored with message \n```\n",
              e, "```\nso the convergence of the current submodel fits cannot ",
              "be checked. Please notify the package maintainer. Current ",
              "submodel formula (right-hand side): ",
              update(formula(outdmin[[1]]), NULL ~ .))
      return(rep(TRUE, length(outdmin)))
    }
  )
  is_conv <- is_conv & (lengths_mssgs_warns == 0)
  not_conv <- sum(!is_conv)
  if (not_conv > 0) {
    if (getOption("projpred.additional_checks", FALSE)) {
      cls <- unique(lapply(outdmin, class))
      stopifnot(length(cls) == 1)
      cls <- cls[[1]]
    } else {
      cls <- class(outdmin[[1]])
    }
    cls <- paste0("c(", paste(paste0("\"", cls, "\""), collapse = ", "), ")")
    warning(
      not_conv, " out of ", length(is_conv), " submodel fits (there is one ",
      "submodel fit per projected draw) might not have converged ",
      "(appropriately). It is recommended to inspect this in detail and (if ",
      "necessary) to adjust tuning parameters via `...` (the ellipsis of the ",
      "employed top-level function such as project(), varsel(), or ",
      "cv_varsel()) or via a custom `div_minimizer` function (an argument of ",
      "init_refmodel()). In the present case, the submodel fits are of ",
      "class(es) `", cls, "`. Documentation for corresponding tuning ",
      "parameters is linked in section \"Draw-wise divergence minimizers\" of ",
      "`` ?`projpred-package` ``. Current submodel formula (right-hand side): ",
      update(formula(outdmin[[1]]), NULL ~ .)
    )
  }
  return()
}

# Helper function for checking the convergence of a single submodel fit (not of
# a whole `outdmin` object):
check_conv_s <- function(fit_s) {
  if (inherits(fit_s, "mmblogit")) {
    return(fit_s$converged)
  } else if (inherits(fit_s, "multinom")) {
    return(fit_s$convergence == 0)
  } else if (inherits(fit_s, "clmm")) {
    return(fit_s$optRes$convergence == 0)
  } else if (inherits(fit_s, "polr")) {
    return(fit_s$convergence == 0)
  } else if (inherits(fit_s, "gam")) {
    # TODO (GAMs): Is this correct?:
    return(fit_s$converged && fit_s$mgcv.conv$fully.converged %||% TRUE)
  } else if (inherits(fit_s, "gamm4")) {
    # TODO (GAMMs): I couldn't find any convergence-related information in
    # element `fit_s$gam`, so the GAM part is currently not checked for
    # convergence. For now, all we can check is the GLMM part from element
    # `fit_s$mer`:
    return(check_conv_s(fit_s$mer))
  } else if (inherits(fit_s, c("lmerMod", "glmerMod"))) {
    # The following was inferred from the source code of lme4::checkConv() and
    # lme4::.prt.warn() (see also `?lme4::mkMerMod`).
    return(fit_s@optinfo$conv$opt == 0 && (
      # Since lme4::.prt.warn() does not refer to `optinfo$conv$lme4$code`,
      # that element might not always exist:
      (!is.null(fit_s@optinfo$conv$lme4$code) &&
         all(fit_s@optinfo$conv$lme4$code == 0)) ||
        is.null(fit_s@optinfo$conv$lme4$code)
    ) && length(unlist(fit_s@optinfo$conv$lme4$messages)) == 0 &&
      length(fit_s@optinfo$warnings) == 0)
  } else if (inherits(fit_s, "glm")) {
    return(fit_s$converged)
  } else if (inherits(fit_s, "lm")) {
    # There doesn't seem to be a better way to check for convergence other than
    # checking `NA` coefficients:
    return(all(!is.na(coef(fit_s))))
  } else if (inherits(fit_s, "subfit")) {
    # For a submodel of class `subfit`, non-convergence is only indicated by
    # output written to the console (which is converted to a warning in
    # fit_glm_ridge_callback() and then checked in check_conv()), so we need to
    # return `TRUE` here:
    return(TRUE)
  } else {
    warning("Unrecognized submodel fit. Please notify the package maintainer.")
    return(TRUE)
  }
}

# Prediction functions for submodels --------------------------------------

subprd <- function(fits, newdata) {
  # lme4 fits of one submodel and `newdata` without new group levels: all draws
  # at once (see subprd_mer_direct()):
  prd_direct <- subprd_mer_direct(fits, newdata)
  if (!is.null(prd_direct)) {
    return(prd_direct)
  }
  # mgcv fits of one submodel: all draws at once (see subprd_gam_direct()):
  if (!is.null(newdata) && length(fits) > 1 && inherits(fits[[1]], "gam")) {
    prd_direct <- subprd_gam_direct(
      fits, cbind(`(Intercept)` = rep(1, NROW(newdata)), newdata)
    )
    if (!is.null(prd_direct)) {
      return(prd_direct)
    }
  }
  prd_list <- lapply(fits, function(fit) {
    is_glmmPQL <- inherits(fit, "glmmPQL")
    is_glmm <- inherits(fit, c("lmerMod", "glmerMod"))
    is_gam_gamm <- inherits(fit, c("gam", "gamm4"))
    if (is_gam_gamm && !is.null(newdata)) {
      newdata <- cbind(`(Intercept)` = rep(1, NROW(newdata)), newdata)
    }
    if (is_glmmPQL) {
      ### TODO (glmmPQL): Remove this as soon as a repair_re.glmmPQL() method
      ### has been added:
      if (!is.null(newdata)) {
        has_new_grps <- sapply(names(fit$groups), function(grp_nm) {
          any(!unique(newdata[[grp_nm]]) %in% unique(fit$groups[[grp_nm]]))
        })
        if (any(has_new_grps)) {
          stop("Under construction (a repair_re.glmmPQL() method needs to be ",
               "added to projpred.")
        }
      }
      ###
      if ("projpred_internal_dummy1s_PQL" %in% names(fit$data) &&
          !is.null(newdata)) {
        if ("projpred_internal_dummy1s_PQL" %in% names(newdata)) {
          stop("Need to write to column `projpred_internal_dummy1s_PQL` of ",
               "`newdata`, but that column already exists. Please rename this ",
               "column in `newdata` and try again.")
        }
        newdata$projpred_internal_dummy1s_PQL <- factor(1)
      }
      return(
        predict(fit, newdata = newdata)
        ### TODO (glmmPQL): Add a repair_re.glmmPQL() method for this:
        # predict(fit, newdata = newdata, level = 0) +
        #   repair_re(fit, newdata = newdata)
        ###
      )
    } else if (is_glmm) {
      if (getOption("projpred.mlvl_pred_new", FALSE)) {
        re_fml_predict <- ~0
      } else {
        re_fml_predict <- NULL
      }
      return(
        predict(fit, newdata = newdata, re.form = re_fml_predict,
                allow.new.levels = TRUE) +
          repair_re(fit, newdata = newdata)
      )
    } else {
      return(predict(fit, newdata = newdata))
    }
  })
  return(do.call(cbind, prd_list))
}

# Linear predictors of lme4 fits (`lmerMod` or `glmerMod`) of the same submodel
# for `newdata` whose group levels all occur in the fits. In this case,
# repair_re() adds zero and draws no random numbers, and predict.merMod() with
# `re.form = NULL` computes X beta_s + Z b_s. This function builds the
# fixed-effects model matrix X and the random-effects design Z for `newdata`
# once, with the steps of predict.merMod() (R package lme4, GPL-2 | GPL-3). It
# maps the conditional modes b_s of each draw (in the row order of the fitted
# Zt) onto the rows of Z by the names of the group levels. These names are
# unique across the grouping variables (see repair_re.merMod()), and within a
# grouping variable they have the same order of terms and coefficients in both.
# The results differ from those of subprd() only by rounding. Returns `NULL`
# for other cases.
subprd_mer_direct <- function(fits, newdata) {
  if (is.null(newdata) || length(fits) == 0 ||
      getOption("projpred.mlvl_pred_new", FALSE) ||
      !all(vapply(fits, inherits, logical(1), c("lmerMod", "glmerMod")))) {
    return(NULL)
  }
  fit_1 <- fits[[1]]
  flist_1 <- lapply(fit_1@flist, levels)
  cnms_1 <- lme4::getME(fit_1, "cnms")
  trms_fix <- delete.response(terms(fit_1, fixed.only = TRUE))
  # The terms are compared without their environment:
  terms_noenv <- function(trms) {
    attr(trms, ".Environment") <- NULL
    trms
  }
  trms_fix_noenv <- terms_noenv(trms_fix)
  same_structure <- vapply(fits, function(fit) {
    trms_s <- delete.response(terms(fit, fixed.only = TRUE))
    identical(lapply(fit@flist, levels), flist_1) &&
      identical(lme4::getME(fit, "cnms"), cnms_1) &&
      identical(terms_noenv(trms_s), trms_fix_noenv) &&
      identical(names(lme4::fixef(fit)), names(lme4::fixef(fit_1)))
  }, logical(1))
  if (!all(same_structure) || !is.null(attr(trms_fix, "offset"))) {
    return(NULL)
  }
  # Group levels: all have to occur in the fits, and the level sets have to be
  # disjoint (otherwise repair_re.merMod() handles the case):
  vnms <- names(flist_1)
  if (!all(vnms %in% names(newdata))) {
    return(NULL)
  }
  for (vnm in vnms) {
    if (!all(as.character(unique(newdata[[vnm]])) %in% flist_1[[vnm]])) {
      return(NULL)
    }
  }
  lvls_all <- unlist(flist_1, use.names = FALSE)
  if (anyDuplicated(lvls_all)) {
    return(NULL)
  }
  # Fixed effects, as in predict.merMod():
  x_contr <- attr(lme4::getME(fit_1, "X"), "contrasts")
  mfr <- model.frame(trms_fix, newdata, na.action = na.pass,
                     xlev = .getXlevels(trms_fix, fit_1@frame))
  x_new <- model.matrix(trms_fix, mfr, contrasts.arg = x_contr)
  fix_nms <- names(lme4::fixef(fit_1))
  if (!all(fix_nms %in% colnames(x_new)) || anyNA(x_new) ||
      nrow(x_new) != nrow(newdata)) {
    return(NULL)
  }
  x_new <- x_new[, fix_nms, drop = FALSE]
  # Random effects:
  re_new <- lme4::mkNewReTrms(fit_1, newdata = newdata,
                              re.form = ("lme4" %:::% "reOnly")(formula(fit_1)),
                              allow.new.levels = TRUE)
  zt_new <- re_new$Zt
  occ_key <- function(nms) paste(nms, stats::ave(seq_along(nms), nms,
                                                 FUN = seq_along), sep = "\r")
  idx_b <- match(occ_key(rownames(zt_new)),
                 occ_key(rownames(lme4::getME(fit_1, "Zt"))))
  if (anyNA(idx_b)) {
    return(NULL)
  }
  b_1 <- as.vector(lme4::getME(fit_1, "b"))
  if (!isTRUE(all.equal(unname(b_1[idx_b]), unname(as.vector(re_new$b)),
                        tolerance = 1e-12))) {
    return(NULL)
  }
  betas <- do.call(cbind, lapply(fits, lme4::fixef))
  bs <- do.call(cbind, lapply(fits, function(fit) {
    as.vector(lme4::getME(fit, "b"))[idx_b]
  }))
  prd <- x_new %*% betas + t(as(t(bs) %*% zt_new, "matrix"))
  # As cbind() of the named vectors from predict.merMod():
  dimnames(prd) <- list(rownames(newdata), NULL)
  return(prd)
}

subprd_augdat <- function(fits, newdata) {
  if (inherits(fits[[1]], "clmm")) {
    y_nm <- extract_terms_response(formula(fits[[1]]))$response
  }
  prd_list <- lapply(fits, function(fit) {
    # Basic predictions (if possible, these should be already on link scale; if
    # this is not possible, the transformation to link scale has to follow):
    if (inherits(fit, "polr")) {
      prbs <- predict(fit, newdata = newdata, type = "probs")
      if (nrow(newdata) == 1) {
        prbs <- t(prbs)
      }
      link_nm <- fit$method
    } else if (inherits(fit, "clmm")) {
      # A predict() method for objects of class `clmm` does not seem to exist,
      # so use a workaround based on ordinal:::predict.clm().
      # Note: There is also `type = "linear.predictor"`, but that returns the
      # predictions for each j and j - 1 (with j = 1, ..., C_cat indexing the
      # response categories), so performs redundant calculations.
      prbs <- predict(structure(fit, class = c(oldClass(fit), "clm")),
                      # Remove the response so ordinal:::predict.clm() returns
                      # predictions for *all* response categories:
                      newdata = within(newdata, assign(y_nm, NULL)),
                      type = "prob")$fit
      link_nm <- fit$link
    } else if (inherits(fit, "multinom")) {
      prbs <- predict(fit, newdata = newdata, type = "probs")
      if (nrow(newdata) == 1) {
        prbs <- t(prbs)
      }
    } else if (inherits(fit, "mmblogit")) {
      # Note: `conditional = FALSE` is used here because
      # mclogit:::predict.mmblogit() doesn't work for new group levels. Instead,
      # repair_re() handles both, existing and new group levels.
      lpreds <- predict(fit, newdata = newdata, conditional = FALSE) +
        repair_re(fit, newdata = newdata)
    } else {
      stop("Unknown `fit`.")
    }
    # Where necessary, transform probabilities to link scale (i.e., to linear
    # predictors):
    if (inherits(fit, c("polr", "clmm"))) {
      prbs <- do.call(rbind, apply(prbs, 1, cumsum, simplify = FALSE))
      prbs <- prbs[, -ncol(prbs), drop = FALSE]
      lpreds <- linkfun_raw(prbs, link_nm = link_nm)
      if (inherits(fit, "clmm")) {
        lpreds <- lpreds + repair_re(fit, newdata = newdata)
      }
    } else if (inherits(fit, "multinom")) {
      lpreds <- log(sweep(prbs[, -1, drop = FALSE], 1, prbs[, 1], "/"))
    }
    return(lpreds)
  })
  return(do.call(abind::abind, list(prd_list, rev.along = 0)))
}

subprd_augdat_binom <- function(fits, newdata) {
  augprd <- subprd(fits, newdata = newdata)
  return(array(augprd, dim = c(nrow(augprd), 1L, ncol(augprd))))
}

## FIXME: find a way that allows us to remove this
#' @noRd
#' @export
predict.subfit <- function(object, newdata = NULL, ...) {
  beta <- object$beta
  alpha <- object$alpha
  if (is.null(newdata)) {
    if (is.null(beta)) {
      return(as.matrix(rep(alpha, nrow(object$x))))
    } else {
      return(cbind(1, object$x) %*% rbind(alpha, beta))
    }
  } else {
    if (is.null(beta)) {
      return(as.matrix(rep(alpha, nrow(newdata))))
    } else {
      # TODO: In the following model.matrix() call, allow user-specified
      # contrasts to be passed to argument `contrasts.arg`. The `contrasts.arg`
      # default (`NULL`) uses `options("contrasts")` internally, but it might be
      # more convenient to let users specify contrasts directly. At that
      # occasion, contrasts should also be tested thoroughly (not done until
      # now).
      x <- model.matrix(delete.response(terms(object$formula)), data = newdata,
                        xlev = object$xlvls)
      # Need to ensure that the columns of `x` match the coefficients (note that
      # `rownames(beta)` is the same as `colnames(object$x)`, see the tests):
      nms_to_use <- c("(Intercept)", colnames(object$x))
      nms_to_use <- reorder_ias(nms_to_use, colnames(x))
      if (all(!is.na(nms_to_use)) && all(nms_to_use %in% colnames(x))) {
        x <- x[, nms_to_use, drop = FALSE]
      } else {
        stop("The column names of the new model matrix don't match the ",
             "column names of the original model matrix. Please notify the ",
             "package maintainer.")
      }
      return(x %*% rbind(alpha, beta))
    }
  }
}

#' @noRd
#' @export
predict.gamm4 <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    newdata <- model.frame(object$mer)
  }
  formula <- object$formula
  random <- object$random
  gamm_struct <- model.matrix_gamm4(delete.response(terms(formula)),
                                    random = random, data = newdata)
  ranef <- lme4::ranef(object$mer) # TODO (GAMMs): Add `, condVar = FALSE` here?
  b <- gamm_struct$b
  mf <- gamm_struct$mf

  ## base pred only smooth and fixed effects
  gamm_pred <- predict(object$mer, newdata = mf, re.form = NA)

  ## gamm4 trick to replace dummy smooth variables with actual smooth terms
  sn <- names(ranef)
  tn <- names(b$reTrms$cnms)
  ind <- seq_along(tn)
  for (i in seq_along(tn)) { ## loop through random effect smooths
    k <- ind[sn[i] == tn] ## which term should contain G$random[[i]]
    ii <- (b$reTrms$Gp[k] + 1):b$reTrms$Gp[k + 1]
    r_pred <- t(as.matrix(b$reTrms$Zt[ii, ])) %*%
      as.matrix(c(as.matrix(ranef[[i]])))
    gamm_pred <- gamm_pred + r_pred
  }
  return(as.matrix(unname(gamm_pred)))
}

## Random-effects adjustments ---------------------------------------------

empty_intersection <- function(x, el_nm = "new") {
  length(intersect(x[[1]][[el_nm]], x[[2]][[el_nm]])) == 0
}

# License/copyright notice: mkNewReTrms_man() is strongly based on (i.e., it was
# only slightly modified from) lme4::mkNewReTrms() from lme4 version 1.1-28 (see
# <https://CRAN.R-project.org/package=lme4>).
#
# The copyright statement for lme4 version 1.1-28 is:
# Copyright (C) 2003-2022 The LME4 Authors (see
# <https://CRAN.R-project.org/package=lme4>).
#
# The license of lme4 version 1.1-28 is:
# "GPL (>=2)" (see <https://CRAN.R-project.org/package=lme4>).
mkNewReTrms_man <- function(re.form, newdata, xlevels, re, D = NULL) {
  stopifnot(!is.null(newdata))
  tt <- terms(reformulas::subbars(re.form))
  rfd <- suppressWarnings(
    model.frame(tt, newdata, na.action = na.pass, xlev = xlevels)
  )
  ReTrms <- reformulas::mkReTrms(reformulas::findbars(re.form[[2]]), rfd)
  ns.re <- names(re)
  nRnms <- names(Rcnms <- ReTrms$cnms)
  if (!all(nRnms %in% ns.re)) {
    stop("grouping factors specified in re.form that were not present in ",
         "original model")
  }
  if (!is.null(D)) {
    # categorical() case
    Zt_old_rnms <- rownames(ReTrms$Zt)
    ReTrms$Zt <- as(as.matrix(ReTrms$Zt %x% t(D[-1, , drop = FALSE])),
                    class(ReTrms$Zt))
    rownames(ReTrms$Zt) <- rep(Zt_old_rnms, each = ncol(D))
    for (vnm in nRnms) {
      ReTrms$cnms[[vnm]] <- names(re[[vnm]])
    }
    Rcnms <- ReTrms$cnms
  }
  new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
  re_x <- Map(
    function(r, n) {
      ("lme4" %:::% "levelfun")(r, n, allow.new.levels = TRUE)
    },
    re[names(new_levels)],
    new_levels
  )
  get_re <- function(rname, cnms) {
    nms <- names(re[[rname]])
    miss_names <- setdiff(cnms, nms)
    if (length(miss_names) > 0) {
      stop("random effects specified in re.form that were not present in ",
           "original model ", paste(miss_names, collapse = ", "))
    }
    t(re_x[[rname]][, cnms])
  }
  re_new <- unlist(Map(get_re, nRnms, Rcnms))
  Zt <- ReTrms$Zt
  attr(Zt, "na.action") <- attr(re_new, "na.action") <- NULL
  list(Zt = Zt, b = re_new)
}

repair_re <- function(object, newdata) {
  UseMethod("repair_re")
}

# For objects of class `merMod`, the following repair_re() method will draw the
# random effects for new group levels from a (multivariate) Gaussian
# distribution, but option `projpred.mlvl_pred_new` determines whether existing
# group levels are also considered as new group levels.
#
# License/copyright notice: repair_re.merMod() is inspired by and uses code
# snippets from lme4:::predict.merMod() from lme4 version 1.1-28 (see
# <https://CRAN.R-project.org/package=lme4>). See the `LICENSE` file in
# projpred's root directory for details.
#
# The copyright statement for lme4 version 1.1-28 is:
# Copyright (C) 2003-2022 The LME4 Authors (see
# <https://CRAN.R-project.org/package=lme4>).
#
# The license of lme4 version 1.1-28 is:
# "GPL (>=2)" (see <https://CRAN.R-project.org/package=lme4>).
#' @noRd
#' @export
repair_re.merMod <- function(object, newdata) {
  stopifnot(!is.null(newdata))
  ranef_tmp <- lme4::ranef(object, condVar = FALSE)
  vnms <- names(ranef_tmp)
  lvls_list <- lapply(setNames(nm = vnms), function(vnm) {
    from_fit <- rownames(ranef_tmp[[vnm]])
    if (!vnm %in% names(newdata)) {
      for (special_char in c("/", ":")) {
        if (any(grepl(paste0("\\|.+", special_char),
                      labels(terms(formula(object)))))) {
          if (special_char == "/") {
            add_hint <- paste0(
              " implied by the `", special_char, "` syntax (see lme4 vignette ",
              "\"Fitting Linear Mixed-Effects Models Using lme4\", Table 2)"
            )
          } else {
            add_hint <- ""
          }
          stop("The `", special_char, "` syntax for grouping variables is ",
               "currently not supported. Please write out (i.e., create ",
               "corresponding columns in the dataset manually and then adapt ",
               "the model formula) interaction terms", add_hint, ".")
        }
      }
      stop("Could not find column `", vnm, "` in `newdata`.")
    }
    from_new <- unique(newdata[, vnm])
    if (is.factor(from_new)) {
      # Strictly speaking, this is not necessary (currently), but include it for
      # safety reasons, in case downstream code is changed in the future (or in
      # case the behavior of `factor`s in R is changed in general):
      from_new <- as.character(from_new)
    }
    list(comb = union(from_fit, from_new),
         exist = intersect(from_new, from_fit),
         new = from_new)
  })
  # In case of duplicated levels across group variables, later code would have
  # to be adapted:
  if (getOption("projpred.mlvl_pred_new", FALSE)) {
    el_nm_mer <- "new"
  } else {
    el_nm_mer <- "comb"
  }
  if (length(lvls_list) >= 2 &&
      !all(utils::combn(lvls_list, 2, empty_intersection, el_nm = el_nm_mer))) {
    stop("Currently, projpred requires all variables with group-level effects ",
         "to have disjoint level sets.")
  }
  re_fml <- ("lme4" %:::% "reOnly")(formula(object))
  ranefs_prep <- lme4::mkNewReTrms(object,
                                   newdata = newdata,
                                   re.form = re_fml,
                                   allow.new.levels = TRUE)
  names(ranefs_prep$b) <- rownames(ranefs_prep$Zt)

  VarCorr_tmp <- lme4::VarCorr(object)
  for (vnm in vnms) {
    lvls_new <- lvls_list[[vnm]]$new
    if (!getOption("projpred.mlvl_pred_new", FALSE)) {
      lvls_exist <- lvls_list[[vnm]]$exist
      lvls_new <- setdiff(lvls_new, lvls_exist)
      ranefs_prep$b[names(ranefs_prep$b) %in% lvls_exist] <- 0
    }
    if (length(lvls_new) > 0) {
      ranefs_prep$b[names(ranefs_prep$b) %in% lvls_new] <- t(mvtnorm::rmvnorm(
        n = length(lvls_new),
        # Add `[, , drop = FALSE]` to drop attributes:
        sigma = VarCorr_tmp[[vnm]][, , drop = FALSE],
        checkSymmetry = FALSE
      ))
    }
  }
  return(drop(as(ranefs_prep$b %*% ranefs_prep$Zt, "matrix")))
}

# For objects of class `clmm`, the following repair_re() method will re-use the
# estimated random effects for existing group levels and will draw the random
# effects for new group levels from a (multivariate) Gaussian distribution, but
# option `projpred.mlvl_pred_new` determines whether existing group levels are
# also considered as new group levels.
#
# License/copyright notice: repair_re.clmm() is inspired by and uses code
# snippets from lme4:::predict.merMod() from lme4 version 1.1-28 (see
# <https://CRAN.R-project.org/package=lme4>). See the `LICENSE` file in
# projpred's root directory for details.
#
# The copyright statement for lme4 version 1.1-28 is:
# Copyright (C) 2003-2022 The LME4 Authors (see
# <https://CRAN.R-project.org/package=lme4>).
#
# The license of lme4 version 1.1-28 is:
# "GPL (>=2)" (see <https://CRAN.R-project.org/package=lme4>).
#' @noRd
#' @export
repair_re.clmm <- function(object, newdata) {
  stopifnot(!is.null(newdata))
  ranef_tmp <- ordinal::ranef(object)
  vnms <- names(ranef_tmp)
  lvls_list <- lapply(setNames(nm = vnms), function(vnm) {
    from_fit <- rownames(ranef_tmp[[vnm]])
    if (!vnm %in% names(newdata)) {
      for (special_char in c("/", ":")) {
        if (any(grepl(paste0("\\|.+", special_char),
                      labels(terms(formula(object)))))) {
          if (special_char == "/") {
            add_hint <- paste0(
              " implied by the `", special_char, "` syntax (see lme4 vignette ",
              "\"Fitting Linear Mixed-Effects Models Using lme4\", Table 2)"
            )
          } else {
            add_hint <- ""
          }
          stop("The `", special_char, "` syntax for grouping variables is ",
               "currently not supported. Please write out (i.e., create ",
               "corresponding columns in the dataset manually and then adapt ",
               "the model formula) interaction terms", add_hint, ".")
        }
      }
      stop("Could not find column `", vnm, "` in `newdata`.")
    }
    from_new <- unique(newdata[, vnm])
    if (is.factor(from_new)) {
      # Strictly speaking, this is not necessary (currently), but include it for
      # safety reasons, in case downstream code is changed in the future (or in
      # case the behavior of `factor`s in R is changed in general):
      from_new <- as.character(from_new)
    }
    list(comb = union(from_fit, from_new),
         exist = intersect(from_new, from_fit),
         new = from_new)
  })
  # In case of duplicated levels across group variables, later code would have
  # to be adapted:
  if (length(lvls_list) >= 2 &&
      !all(utils::combn(lvls_list, 2, empty_intersection))) {
    stop("Currently, projpred requires all variables with group-level effects ",
         "to have disjoint level sets.")
  }
  re_fml <- ("lme4" %:::% "reOnly")(formula(object))
  ranefs_prep <- mkNewReTrms_man(re.form = re_fml,
                                 newdata = newdata,
                                 xlevels = c(lapply(lvls_list, "[[", "comb"),
                                             object$xlevels),
                                 re = ranef_tmp)
  names(ranefs_prep$b) <- rownames(ranefs_prep$Zt)

  VarCorr_tmp <- ordinal::VarCorr(object)
  for (vnm in vnms) {
    lvls_new <- lvls_list[[vnm]]$new
    if (!getOption("projpred.mlvl_pred_new", FALSE)) {
      lvls_exist <- lvls_list[[vnm]]$exist
      lvls_new <- setdiff(lvls_new, lvls_exist)
    }
    if (length(lvls_new) > 0) {
      ranefs_prep$b[names(ranefs_prep$b) %in% lvls_new] <- t(mvtnorm::rmvnorm(
        n = length(lvls_new),
        # Add `[, , drop = FALSE]` to drop attributes:
        sigma = VarCorr_tmp[[vnm]][, , drop = FALSE],
        checkSymmetry = FALSE
      ))
    }
  }
  return(-drop(as(ranefs_prep$b %*% ranefs_prep$Zt, "matrix")))
}

# For objects of class `mmblogit`, the following repair_re() method will re-use
# the estimated random effects for existing group levels and will draw the
# random effects for new group levels from a (multivariate) Gaussian
# distribution, but option `projpred.mlvl_pred_new` determines whether existing
# group levels are also considered as new group levels.
#
# License/copyright notice: repair_re.mmblogit() is inspired by and uses code
# snippets from lme4:::predict.merMod() from lme4 version 1.1-28 (see
# <https://CRAN.R-project.org/package=lme4>). See the `LICENSE` file in
# projpred's root directory for details.
#
# The copyright statement for lme4 version 1.1-28 is:
# Copyright (C) 2003-2022 The LME4 Authors (see
# <https://CRAN.R-project.org/package=lme4>).
#
# The license of lme4 version 1.1-28 is:
# "GPL (>=2)" (see <https://CRAN.R-project.org/package=lme4>).
#' @noRd
#' @export
repair_re.mmblogit <- function(object, newdata) {
  stopifnot(!is.null(newdata))
  vnms <- names(object$groups)
  stopifnot(length(vnms) == length(object$random.effects))
  if (utils::packageVersion("mclogit") < "0.9") {
    stopifnot(length(vnms) == 1)
  } else if (length(vnms) < length(object$random)) {
    stop("The length of `<mmblogit_object>$random` is greater than that of ",
         "`vnms = c(\"", paste(vnms, collapse = "\", \""), "\")`.")
  }
  # The number of latent response categories:
  nlats <- ncol(object$D)
  # Coerce the random effects into the same format as the output of ranef() from
  # packages 'lme4' and 'ordinal':
  ranef_tmp <- lapply(
    setNames(seq_along(object$random.effects), vnms),
    function(grp_idx) {
      if (utils::packageVersion("mclogit") < "0.9") {
        obj_rand <- object$random
        VarCov_idx <- 1
      } else {
        if (length(object$random) < grp_idx) {
          stop("Unexpected length of `<mmblogit_object>$random`.")
        }
        obj_rand <- object$random[[grp_idx]]
        VarCov_idx <- vnms[grp_idx]
      }
      if (length(obj_rand$groups) > 1) {
        stop("It seems like you fitted a model with a nested group-level ",
             "term. Currently, these are not supported.")
      }
      # The number of "random slopes", including the intercept (the intercept is
      # the reason for `+ 1L`):
      ncoefs <- length(all.vars(obj_rand$formula)) + 1L
      # The coercion itself:
      ranef_i <- matrix(
        object$random.effects[[grp_idx]],
        nrow = nlats * ncoefs,
        ncol = nlevels(object$groups[[vnms[grp_idx]]]),
        dimnames = list(colnames(object$VarCov[[VarCov_idx]]),
                        levels(object$groups[[vnms[grp_idx]]]))
      )
      return(as.data.frame(t(ranef_i)))
    }
  )
  lvls_list <- lapply(setNames(nm = vnms), function(vnm) {
    from_fit <- levels(object$groups[[vnm]])
    if (!vnm %in% names(newdata)) {
      stop("Could not find column `", vnm, "` in `newdata`.")
    }
    from_new <- unique(newdata[, vnm])
    if (is.factor(from_new)) {
      # Strictly speaking, this is not necessary (currently), but include it for
      # safety reasons, in case downstream code is changed in the future (or in
      # case the behavior of `factor`s in R is changed in general):
      from_new <- as.character(from_new)
    }
    list(comb = union(from_fit, from_new),
         exist = intersect(from_new, from_fit),
         new = from_new)
  })
  # In case of duplicated levels across group variables, later code would have
  # to be adapted:
  if (length(lvls_list) >= 2 &&
      !all(utils::combn(lvls_list, 2, empty_intersection))) {
    stop("Currently, projpred requires all variables with group-level effects ",
         "to have disjoint level sets.")
  }
  # Create the lme4-type random-effects formula needed by mkNewReTrms_man():
  if (utils::packageVersion("mclogit") < "0.9") {
    re_fml <- update(object$random$formula,
                     paste("~ . |", object$random$groups))
  } else {
    for (grp_idx in seq_along(object$random)) {
      re_fml_i <- update(object$random[[grp_idx]]$formula,
                         paste("~ . |", object$random[[grp_idx]]$groups))
      if (grp_idx == 1) {
        re_fml <- re_fml_i
      } else {
        re_fml <- update(re_fml, paste("~ . +", as.character(re_fml_i[2])))
      }
    }
  }
  # Get random effects for existing group levels (and zero for new group levels)
  # as well as the random-effects model matrix:
  ranefs_prep <- mkNewReTrms_man(re.form = re_fml,
                                 newdata = newdata,
                                 xlevels = c(lapply(lvls_list, "[[", "comb"),
                                             object$xlevels),
                                 re = ranef_tmp,
                                 D = object$D)
  names(ranefs_prep$b) <- rownames(ranefs_prep$Zt)

  VarCorr_tmp <- object$VarCov
  if (utils::packageVersion("mclogit") < "0.9") {
    VarCorr_tmp <- setNames(VarCorr_tmp, vnms)
  }
  for (vnm in vnms) {
    lvls_new <- lvls_list[[vnm]]$new
    if (!getOption("projpred.mlvl_pred_new", FALSE)) {
      lvls_exist <- lvls_list[[vnm]]$exist
      lvls_new <- setdiff(lvls_new, lvls_exist)
    }
    if (length(lvls_new) > 0) {
      ranefs_prep$b[names(ranefs_prep$b) %in% lvls_new] <- t(mvtnorm::rmvnorm(
        n = length(lvls_new),
        sigma = VarCorr_tmp[[vnm]],
        checkSymmetry = FALSE
      ))
    }
  }
  re_vec <- drop(as(ranefs_prep$b %*% ranefs_prep$Zt, "matrix"))
  re_mat <- matrix(re_vec, nrow = nlats, ncol = nrow(newdata),
                   dimnames = list(colnames(object$D), NULL))
  return(t(re_mat))
}
