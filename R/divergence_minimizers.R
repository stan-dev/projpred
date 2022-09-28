# Divergence minimizers ---------------------------------------------------

## Traditional (and latent) projection ------------------------------------

# Needed to avoid a NOTE in `R CMD check`:
if (getRversion() >= package_version("2.15.1")) {
  utils::globalVariables("formula_s")
  utils::globalVariables("projpred_var_s")
  utils::globalVariables("projpred_formula_no_random_s")
}

divmin <- function(formula, projpred_var, ...) {
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
      # use lme4::nobars() and lme4::findbars() here):
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
  if (is.list(projpred_formulas_no_random)) {
    stopifnot(length(projpred_formulas_no_random) == length(formulas))
  } else {
    projpred_formulas_no_random <- as.list(
      rep(projpred_formulas_no_random, length(formulas))
    )
  }

  if (length(formulas) < getOption("projpred.prll_prj_trigger", Inf)) {
    # Sequential case. Actually, we could simply use ``%do_projpred%` <-
    # foreach::`%do%`` here and then proceed as in the parallel case, but that
    # would require adding more "hard" dependencies (because packages 'foreach'
    # and 'iterators' would have to be moved from `Suggests:` to `Imports:`).
    return(lapply(seq_along(formulas), function(s) {
      sdivmin(
        formula = formulas[[s]],
        projpred_var = projpred_var[, s, drop = FALSE],
        projpred_formula_no_random = projpred_formulas_no_random[[s]],
        projpred_random = projpred_random,
        ...
      )
    }))
  } else {
    # Parallel case.
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Please install the 'foreach' package.")
    }
    if (!requireNamespace("iterators", quietly = TRUE)) {
      stop("Please install the 'iterators' package.")
    }
    dot_args <- list(...)
    `%do_projpred%` <- foreach::`%dopar%`
    return(foreach::foreach(
      formula_s = formulas,
      projpred_var_s = iterators::iter(projpred_var, by = "column"),
      projpred_formula_no_random_s = projpred_formulas_no_random,
      .export = c("sdivmin", "projpred_random", "dot_args"),
      .noexport = c(
        "object", "p_sel", "p_pred", "search_path", "p_ref", "refmodel",
        "formulas", "projpred_var", "projpred_formulas_no_random"
      )
    ) %do_projpred% {
      do.call(
        sdivmin,
        c(list(formula = formula_s,
               projpred_var = projpred_var_s,
               projpred_formula_no_random = projpred_formula_no_random_s,
               projpred_random = projpred_random),
          dot_args)
      )
    })
  }
}

# Use projpred's own implementation to fit non-multilevel non-additive
# submodels:
fit_glm_ridge_callback <- function(formula, data,
                                   projpred_var = matrix(nrow = nrow(data)),
                                   projpred_regul = 1e-4, ...) {
  # Preparations:
  fr <- model.frame(formula, data = data)
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
  fit <- do.call(glm_ridge, c(
    list(x = x, y = y, lambda = projpred_regul, obsvar = projpred_var),
    dot_args
  ))
  # Post-processing:
  rownames(fit$beta) <- colnames(x)
  sub <- nlist(
    alpha = fit$beta0,
    beta = fit$beta,
    w = fit$w,
    formula,
    x, y
  )
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
    return(suppressMessages(suppressWarnings(do.call(stats::lm, c(
      list(formula = formula),
      dot_args
    )))))
  } else {
    # Exclude arguments from `...` which cannot be passed to stats::glm():
    dot_args <- list(...)
    dot_args <- dot_args[intersect(
      names(dot_args),
      c(methods::formalArgs(stats::glm),
        methods::formalArgs(stats::glm.control))
    )]
    # Call the submodel fitter:
    return(suppressMessages(suppressWarnings(do.call(stats::glm, c(
      list(formula = formula, family = family),
      dot_args
    )))))
  }
}

# Use package "mgcv" to fit additive non-multilevel submodels:
#' @importFrom mgcv gam
fit_gam_callback <- function(formula, ...) {
  # Exclude arguments from `...` which cannot be passed to mgcv::gam():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(gam),
      methods::formalArgs(mgcv::gam.fit))
  )]
  # Call the submodel fitter:
  return(suppressMessages(suppressWarnings(do.call(gam, c(
    list(formula = formula),
    dot_args
  )))))
}

# Use package "gamm4" to fit additive multilevel submodels:
#' @importFrom gamm4 gamm4
fit_gamm_callback <- function(formula, projpred_formula_no_random,
                              projpred_random, data, family,
                              control = control_callback(family), ...) {
  # Exclude arguments from `...` which cannot be passed to gamm4::gamm4():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    c(methods::formalArgs(gamm4),
      methods::formalArgs(lme4::lFormula),
      methods::formalArgs(lme4::glFormula))
  )]
  # Call the submodel fitter:
  fit <- tryCatch({
    suppressMessages(suppressWarnings(do.call(gamm4, c(
      list(formula = projpred_formula_no_random, random = projpred_random,
           data = data, family = family, control = control),
      dot_args
    ))))
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
      return(suppressMessages(suppressWarnings(do.call(MASS::glmmPQL, c(
        list(fixed = projpred_formula_no_random, random = random_fmls,
             family = family, control = control),
        dot_args
      )))))
    } else if (family$family == "gaussian" && family$link == "identity" &&
               getOption("projpred.gaussian_not_as_generalized", TRUE)) {
      # Exclude arguments from `...` which cannot be passed to lme4::lmer():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(lme4::lmer)
      )]
      # Call the submodel fitter:
      return(suppressMessages(suppressWarnings(do.call(lme4::lmer, c(
        list(formula = formula, control = control),
        dot_args
      )))))
    } else {
      # Exclude arguments from `...` which cannot be passed to lme4::glmer():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(lme4::glmer)
      )]
      # Call the submodel fitter:
      return(suppressMessages(suppressWarnings(do.call(lme4::glmer, c(
        list(formula = formula, family = family, control = control),
        dot_args
      )))))
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
      if (length(control$optCtrl$maxfun) > 0) {
        maxfun_new <- 10 * control$optCtrl$maxfun
      } else {
        maxfun_new <- 1e4
      }
      if (length(control$optCtrl$maxit) > 0) {
        maxit_new <- 10 * control$optCtrl$maxit
      } else {
        maxit_new <- 1e4
      }
      if (tolPwrss_new > 1e-4 && maxfun_new > 1e7 && maxit_new > 1e7) {
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
                                   optCtrl = list(maxfun = maxfun_new,
                                                  maxit = maxit_new)),
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
  utils::globalVariables("s")
  utils::globalVariables("projpred_internal_w_aug")
}

divmin_augdat <- function(formula, data, family, weights, projpred_var,
                          projpred_ws_aug, ...) {
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

  if (ncol(projpred_ws_aug) < getOption("projpred.prll_prj_trigger", Inf)) {
    # Sequential case. Actually, we could simply use ``%do_projpred%` <-
    # foreach::`%do%`` here and then proceed as in the parallel case, but that
    # would require adding more "hard" dependencies (because packages 'foreach'
    # and 'iterators' would have to be moved from `Suggests:` to `Imports:`).
    return(lapply(seq_len(ncol(projpred_ws_aug)), function(s) {
      sdivmin(
        formula = formula,
        data = data,
        family = family,
        weights = projpred_ws_aug[, s],
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        ...
      )
    }))
  } else {
    # Parallel case.
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Please install the 'foreach' package.")
    }
    # Unfortunately, iterators::iter() seems to conflict with augmented-data
    # matrices (see note below). Thus, `requireNamespace("iterators")` is not
    # necessary here.
    dot_args <- list(...)
    `%do_projpred%` <- foreach::`%dopar%`
    return(foreach::foreach(
      # Unfortunately, iterators::iter() seems to conflict with augmented-data
      # matrices, even if using unclass(). Thus, iterate over the indices:
      s = seq_len(ncol(projpred_ws_aug)),
      .export = c(
        "sdivmin", "formula", "data", "family", "projpred_ws_aug",
        "projpred_formula_no_random", "projpred_random", "dot_args"
      ),
      .noexport = c(
        "object", "p_sel", "p_pred", "search_path", "p_ref", "refmodel",
        "linkobjs"
      )
    ) %do_projpred% {
      do.call(
        sdivmin,
        c(list(formula = formula,
               data = data,
               family = family,
               weights = projpred_ws_aug[, s],
               projpred_formula_no_random = projpred_formula_no_random,
               projpred_random = projpred_random),
          dot_args)
      )
    })
  }
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
  # For catching warnings via capture.output() (which is necessary to filter out
  # the warning "non-integer #successes in a binomial glm!"):
  warn_orig <- options(warn = 1)
  on.exit(options(warn_orig))
  # Call the submodel fitter:
  warn_capt <- utils::capture.output({
    fitobj <- try(do.call(MASS::polr, c(
      list(formula = formula,
           data = data,
           weights = quote(projpred_internal_w_aug),
           model = FALSE,
           method = link_nm),
      dot_args
    )), silent = TRUE)
  }, type = "message")
  if (inherits(fitobj, "try-error") &&
      grepl(paste("initial value in 'vmmin' is not finite",
                  "attempt to find suitable starting values failed",
                  sep = "|"),
            attr(fitobj, "condition")$message)) {
    # Try to fix this automatically by specifying `start` values.
    ncoefs <- count_terms_in_formula(formula) -
      attr(terms(formula), "intercept")
    start_coefs <- rep(0, ncoefs)
    ncats <- length(unique(eval(formula[[2]], data, environment(formula))))
    nthres <- ncats - 1L
    # Start with thresholds which imply equal probabilities for the response
    # categories:
    start_thres <- linkfun_raw(seq_len(nthres) / ncats, link_nm = link_nm)
    warn_capt <- utils::capture.output({
      fitobj <- try(do.call(MASS::polr, c(
        list(formula = formula,
             data = data,
             weights = quote(projpred_internal_w_aug),
             model = FALSE,
             method = link_nm,
             start = c(start_coefs, start_thres)),
        dot_args
      )), silent = TRUE)
    }, type = "message")
  }
  if (inherits(fitobj, "try-error")) {
    stop(attr(fitobj, "condition")$message)
  }
  warn_capt <- setdiff(
    warn_capt,
    c("Warning in eval(family$initialize) :",
      "  non-integer #successes in a binomial glm!")
  )
  if (length(warn_capt) > 0) {
    warning(warn_capt)
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
      methods::formalArgs(ucminf::ucminf),
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
  # For catching warnings via capture.output() (which is necessary to filter out
  # the warning "Using formula(x) is deprecated when x is a character vector of
  # length > 1. [...]"):
  warn_orig <- options(warn = 1)
  on.exit(options(warn_orig))
  # Call the submodel fitter:
  warn_capt <- utils::capture.output({
    fitobj <- try(do.call(ordinal::clmm, c(
      list(formula = formula,
           data = data,
           weights = quote(projpred_internal_w_aug),
           contrasts = NULL,
           Hess = FALSE,
           model = FALSE,
           link = link_nm),
      dot_args
    )), silent = TRUE)
  }, type = "message")
  if (inherits(fitobj, "try-error")) {
    stop(attr(fitobj, "condition")$message)
  }
  warn_capt <- setdiff(
    warn_capt,
    c(paste("Warning: Using formula(x) is deprecated when x is a character",
            "vector of length > 1."),
      "  Consider formula(paste(x, collapse = \" \")) instead.")
  )
  if (length(warn_capt) > 0) {
    warning(warn_capt)
  }
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
  out_capt <- utils::capture.output({
    fitobj <- do.call(nnet::multinom, c(
      list(formula = formula,
           data = data,
           weights = weights),
      dot_args
    ))
  })
  if (utils::tail(out_capt, 1) != "converged") {
    warning("The nnet::multinom() submodel fit did not converge.")
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
  out_capt <- utils::capture.output({
    fitobj <- do.call(mclogit::mblogit, c(
      list(formula = projpred_formula_no_random,
           data = data,
           random = random_fmls,
           weights = weights,
           model = FALSE,
           x = FALSE,
           y = FALSE,
           estimator = "ML",
           dispersion = FALSE,
           from.table = FALSE,
           groups = NULL),
      dot_args
    ))
  })
  if (utils::tail(out_capt, 1) != "converged") {
    warning("The mclogit::mblogit() submodel fit did not converge.")
  }
  return(fitobj)
}

# Convergence checker -----------------------------------------------------

check_conv <- function(fit) {
  conv_info <- do.call(cbind, lapply(fit, function(fit_s) {
    if (inherits(fit_s, "gam")) {
      # TODO (GAMs):
      #   1. For GAMs, there is also `fit_s$mgcv.conv` (see
      #   `?mgcv::gamObject`). Do we need to take this into account?
      #   2. If there is a (convenient) way to retrieve warnings, then this
      #   should be done to get a sensible value for `no_warnings` below.
      return(c(no_gross_fail = fit_s$converged, no_warnings = TRUE))
    } else if (inherits(fit_s, "gamm4")) {
      # TODO (GAMMs): Both, `no_gross_fail` and `no_warnings` need to be
      # implemented. Return `TRUE` for now.
      return(c(no_gross_fail = TRUE, no_warnings = TRUE))
    } else if (inherits(fit_s, c("lmerMod", "glmerMod"))) {
      # The following was inferred from the source code of lme4::checkConv() and
      # lme4::.prt.warn() (see also `?lme4::mkMerMod`).
      return(c(
        no_gross_fail = fit_s@optinfo$conv$opt == 0 && (
          # Since lme4::.prt.warn() does not refer to `optinfo$conv$lme4$code`,
          # that element might not always exist:
          (!is.null(fit_s@optinfo$conv$lme4$code) &&
             fit_s@optinfo$conv$lme4$code >= 0) ||
            is.null(fit_s@optinfo$conv$lme4$code)
        ),
        no_warnings = length(fit_s@optinfo$warnings) &&
          length(unlist(fit_s@optinfo$conv$lme4$messages)) == 0 && (
            # Since lme4::.prt.warn() does not refer to `optinfo$conv$lme4$code`,
            # that element might not always exist:
            (!is.null(fit_s@optinfo$conv$lme4$code) &&
               fit_s@optinfo$conv$lme4$code == 0) ||
              is.null(fit_s@optinfo$conv$lme4$code)
          )
      ))
    } else if (inherits(fit_s, "glm")) {
      # TODO (GLMs): If there is a (convenient) way to retrieve warnings, then
      # this should be done to get a sensible value for `no_warnings` below.
      return(c(no_gross_fail = fit_s$converged, no_warnings = TRUE))
    } else if (inherits(fit_s, "lm")) {
      # Note: There doesn't seem to be a better way to check for convergence
      # other than checking `NA` coefficients (see below).
      # TODO (LMs): If there is a (convenient) way to retrieve warnings, then
      # this should be done to get a sensible value for `no_warnings` below.
      return(c(no_gross_fail = all(!is.na(coef(fit_s))), no_warnings = TRUE))
    } else if (inherits(fit_s, "subfit")) {
      # Note: There doesn't seem to be any way to check for convergence, so
      # return `TRUE` for now.
      # TODO (GLMs with ridge regularization): Add a logical indicating
      # convergence to objects of class `subfit` (i.e., from glm_ridge())?
      return(c(no_gross_fail = TRUE, no_warnings = TRUE))
    } else {
      stop("Unrecognized submodel fit. Please notify the package maintainer.")
    }
  }))
  is_conv <- conv_info["no_gross_fail", ]
  if (any(!is_conv)) {
    warning(sum(!is_conv), " out of ", length(is_conv), " submodel fits ",
            "(there is one submodel fit per projected draw) did not converge. ",
            "Formula (right-hand side): ", update(formula(fit[[1]]), NULL ~ .))
  }
  no_warns <- conv_info["no_warnings", ]
  if (any(!no_warns)) {
    warning(sum(!no_warns), " out of ", length(no_warns), " submodel fits ",
            "(there is one submodel fit per projected draw) threw a warning ",
            "which might be relevant for convergence. ",
            "Formula (right-hand side): ", update(formula(fit[[1]]), NULL ~ .))
  }
  return(invisible(TRUE))
}

# Prediction functions for submodels --------------------------------------

subprd <- function(fits, newdata) {
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
      return(predict(fit, newdata = newdata, allow.new.levels = TRUE) +
               repair_re(fit, newdata = newdata))
    } else {
      return(predict(fit, newdata = newdata))
    }
  })
  return(do.call(cbind, prd_list))
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
predict.subfit <- function(subfit, newdata = NULL) {
  beta <- subfit$beta
  alpha <- subfit$alpha
  if (is.null(newdata)) {
    if (is.null(beta)) {
      return(as.matrix(rep(alpha, NROW(subfit$x))))
    } else {
      return(cbind(1, subfit$x) %*% rbind(alpha, beta))
    }
  } else {
    # TODO: In the following model.matrix() call, allow user-specified contrasts
    # to be passed to argument `contrasts.arg`. The `contrasts.arg` default
    # (`NULL`) uses `options("contrasts")` internally, but it might be more
    # convenient to let users specify contrasts directly. At that occasion,
    # contrasts should also be tested thoroughly (not done until now).
    x <- model.matrix(delete.response(terms(subfit$formula)), data = newdata)
    if (is.null(beta)) {
      return(as.matrix(rep(alpha, NROW(x))))
    } else {
      if (ncol(x) != length(beta) + 1L) {
        stop("The number of columns in the model matrix (\"X\") doesn't match ",
             "the number of coefficients.")
      }
      return(x %*% rbind(alpha, beta))
    }
  }
}

predict.gamm4 <- function(fit, newdata = NULL) {
  if (is.null(newdata)) {
    newdata <- model.frame(fit$mer)
  }
  formula <- fit$formula
  random <- fit$random
  gamm_struct <- model.matrix_gamm4(delete.response(terms(formula)),
                                    random = random, data = newdata)
  ranef <- lme4::ranef(fit$mer) # TODO (GAMMs): Add `, condVar = FALSE` here?
  b <- gamm_struct$b
  mf <- gamm_struct$mf

  ## base pred only smooth and fixed effects
  gamm_pred <- predict(fit$mer, newdata = mf, re.form = NA)

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

empty_intersection_comb <- function(x) {
  length(intersect(x[[1]]$comb, x[[2]]$comb)) == 0
}

empty_intersection_new <- function(x) {
  length(intersect(x[[1]]$new, x[[2]]$new)) == 0
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
  tt <- terms(lme4::subbars(re.form))
  rfd <- suppressWarnings(
    model.frame(tt, newdata, na.action = na.pass, xlev = xlevels)
  )
  ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
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
# distribution.
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
      if (any(grepl("\\|.+/", labels(terms(formula(object)))))) {
        stop("The `/` syntax for nested group-level terms is currently not ",
             "supported. Please try to write out the interaction term implied ",
             "by the `/` syntax (see Table 2 in lme4's vignette called ",
             "\"Fitting Linear Mixed-Effects Models Using lme4\").")
      } else {
        stop("Could not find column `", vnm, "` in `newdata`.")
      }
    }
    from_new <- levels(as.factor(newdata[, vnm]))
    list(comb = union(from_fit, from_new),
         exist = intersect(from_new, from_fit),
         new = setdiff(from_new, from_fit))
  })
  # In case of duplicated levels across group variables, later code would have
  # to be adapted:
  if (length(lvls_list) >= 2 &&
      !all(utils::combn(lvls_list, 2, empty_intersection_comb))) {
    stop("Currently, projpred requires all variables with group-level effects ",
         "to have disjoint level sets.")
  }
  re_fml <- ("lme4" %:::% "reOnly")(formula(object))
  # Note: Calling lme4::mkNewReTrms() with `re.form = NULL` fails.
  ranefs_prep <- lme4::mkNewReTrms(object,
                                   newdata = newdata,
                                   re.form = re_fml,
                                   allow.new.levels = TRUE)
  names(ranefs_prep$b) <- rownames(ranefs_prep$Zt)

  VarCorr_tmp <- lme4::VarCorr(object)
  for (vnm in vnms) {
    lvls_exist <- lvls_list[[vnm]]$exist
    lvls_new <- lvls_list[[vnm]]$new
    ranefs_prep$b[names(ranefs_prep$b) %in% lvls_exist] <- 0
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
# effects for new group levels from a (multivariate) Gaussian distribution.
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
      if (any(grepl("\\|.+/", labels(terms(formula(object)))))) {
        stop("The `/` syntax for nested group-level terms is currently not ",
             "supported. Please try to write out the interaction term implied ",
             "by the `/` syntax (see Table 2 in lme4's vignette called ",
             "\"Fitting Linear Mixed-Effects Models Using lme4\").")
      } else {
        stop("Could not find column `", vnm, "` in `newdata`.")
      }
    }
    from_new <- levels(as.factor(newdata[, vnm]))
    list(comb = union(from_fit, from_new),
         new = setdiff(from_new, from_fit))
  })
  # In case of duplicated levels across group variables, later code would have
  # to be adapted:
  if (length(lvls_list) >= 2 &&
      !all(utils::combn(lvls_list, 2, empty_intersection_new))) {
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
# distribution.
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
    from_new <- levels(as.factor(newdata[, vnm]))
    list(comb = union(from_fit, from_new),
         new = setdiff(from_new, from_fit))
  })
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
