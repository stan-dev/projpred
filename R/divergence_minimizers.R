# Divergence minimizers ---------------------------------------------------

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
  if (!has_grp && !has_add) {
    sdivmin <- get(getOption("projpred.glm_fitter", "fit_glm_ridge_callback"),
                   mode = "function")
  } else if (has_grp && !has_add) {
    sdivmin <- fit_glmer_callback
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

fit_glm_ridge_callback <- function(formula, data,
                                   projpred_var = matrix(nrow = nrow(data)),
                                   projpred_regul = 1e-4, ...) {
  fr <- model.frame(delete.intercept(formula), data = data)
  contrasts_arg <- get_contrasts_arg_list(formula, data = data)
  x <- model.matrix(fr, data = data, contrasts.arg = contrasts_arg)
  y <- model.response(fr)
  # Exclude arguments from `...` which cannot be passed to glm_ridge():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    methods::formalArgs(glm_ridge)
  )]
  fit <- do.call(glm_ridge, c(
    list(x = x, y = y, lambda = projpred_regul, obsvar = projpred_var),
    dot_args
  ))
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
fit_glm_callback <- function(formula, family, projpred_var, projpred_regul,
                             ...) {
  ## make sure correct 'weights' can be found
  environment(formula) <- environment()
  tryCatch({
    if (family$family == "gaussian" && family$link == "identity") {
      # Exclude arguments from `...` which cannot be passed to stats::lm():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(stats::lm)
      )]
      # Use suppressMessages(suppressWarnings()) here?:
      return(do.call(stats::lm, c(
        list(formula = formula),
        dot_args
      )))
    } else {
      # Exclude arguments from `...` which cannot be passed to stats::glm():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(stats::glm)
      )]
      # Use suppressMessages(suppressWarnings()) here?:
      return(do.call(stats::glm, c(
        list(formula = formula, family = family),
        dot_args
      )))
    }
  }, error = function(e) {
    # May be used to handle errors.
    stop(e)
  })
}

# Use package "mgcv" to fit additive non-multilevel submodels:
#' @importFrom mgcv gam
fit_gam_callback <- function(formula, ...) {
  # make sure correct 'weights' can be found
  environment(formula) <- environment()
  # Exclude arguments from `...` which cannot be passed to mgcv::gam():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    union(methods::formalArgs(gam),
          methods::formalArgs(mgcv::gam.fit))
  )]
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
  # make sure correct 'weights' can be found
  environment(projpred_formula_no_random) <- environment()
  # Exclude arguments from `...` which cannot be passed to gamm4::gamm4():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    union(union(methods::formalArgs(gamm4),
                methods::formalArgs(lme4::lFormula)),
          methods::formalArgs(lme4::glFormula))
  )]
  fit <- tryCatch({
    suppressMessages(suppressWarnings(do.call(gamm4, c(
      list(formula = projpred_formula_no_random, random = projpred_random,
           data = data, family = family, control = control),
      dot_args
    ))))
  }, error = function(e) {
    if (grepl("not positive definite", as.character(e))) {
      scaled_data <- preprocess_data(data, projpred_formula_no_random)
      fit_gamm_callback(
        formula, projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random, data = scaled_data, family = family,
        control = control_callback(family,
                                   optimizer = "optimx",
                                   optCtrl = list(method = "nlminb")),
        ...
      )
    } else {
      stop(e)
    }
  })

  fit$random <- projpred_random
  fit$formula <- projpred_formula_no_random
  class(fit) <- c("gamm4")
  return(fit)
}

# Use package "lme4" to fit submodels for multilevel reference models (with a
# fallback to "projpred"'s own implementation for fitting non-multilevel (and
# non-additive) submodels):
fit_glmer_callback <- function(formula, family,
                               control = control_callback(family), ...) {
  ## make sure correct 'weights' can be found
  environment(formula) <- environment()
  tryCatch({
    if (family$family == "gaussian" && family$link == "identity") {
      # Exclude arguments from `...` which cannot be passed to lme4::lmer():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(lme4::lmer)
      )]
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
      return(suppressMessages(suppressWarnings(do.call(lme4::glmer, c(
        list(formula = formula, family = family,
             control = control),
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
        formula, family = family, ...
      ))
    } else if (grepl("not positive definite", as.character(e))) {
      if ("optimx" %in% control$optimizer &&
          length(control$optCtrl$method) > 0 &&
          control$optCtrl$method == "nlminb") {
        stop("Encountering the `not positive definite` error while running ",
             "the lme4 fitting procedure, but cannot fix this automatically ",
             "anymore.")
      }
      return(fit_glmer_callback(
        formula, family = family,
        control = control_callback(family,
                                   optimizer = "optimx",
                                   optCtrl = list(method = "nlminb")),
        ...
      ))
    } else if (grepl("PIRLS step-halvings", as.character(e))) {
      if (length(dot_args$nAGQ) > 0) {
        nAGQ_new <- dot_args$nAGQ + 1L
      } else {
        nAGQ_new <- 20L
      }
      if (nAGQ_new > 30L) {
        stop("Encountering the `PIRLS step-halvings` error while running the ",
             "lme4 fitting procedure, but cannot fix this automatically ",
             "anymore.")
      }
      return(fit_glmer_callback(
        formula, family = family, control = control, nAGQ = nAGQ_new, ...
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
             "automatically anymore.")
      }
      return(fit_glmer_callback(
        formula, family = family,
        control = control_callback(family, tolPwrss = tolPwrss_new,
                                   optCtrl = list(maxfun = maxfun_new,
                                                  maxit = maxit_new)),
        ...
      ))
    } else {
      stop(e)
    }
  })
}

preprocess_data <- function(data, formula) {
  tt <- extract_terms_response(formula)
  non_group_terms <- c(tt$individual_terms, tt$interaction_terms)
  X <- data %>%
    dplyr::select(non_group_terms) %>%
    scale()
  data[, non_group_terms] <- X
  return(data)
}

# Helper function for fit_glmer_callback() and fit_gamm_callback() to get the
# appropriate control options depending on the family:
control_callback <- function(family, ...) {
  if (family$family == "gaussian" && family$link == "identity") {
    return(lme4::lmerControl(...))
  } else {
    return(lme4::glmerControl(...))
  }
}

# Prediction functions for submodels --------------------------------------

subprd <- function(fit, newdata) {
  return(do.call(cbind, lapply(fit, function(fit) {
    # Only pass argument `allow.new.levels` to the predict() generic if the fit
    # is multilevel:
    has_grp <- inherits(fit, c("lmerMod", "glmerMod"))
    has_add <- inherits(fit, c("gam", "gamm4"))
    if (has_add && !is.null(newdata)) {
      newdata <- cbind(`(Intercept)` = rep(1, NROW(newdata)), newdata)
    }
    if (!has_grp) {
      return(predict(fit, newdata = newdata))
    } else {
      return(predict(fit, newdata = newdata, allow.new.levels = TRUE))
    }
  })))
}

## FIXME: find a way that allows us to remove this
predict.subfit <- function(subfit, newdata = NULL) {
  beta <- subfit$beta
  alpha <- subfit$alpha
  if (is.null(newdata)) {
    if (is.null(beta)) {
      return(as.matrix(rep(alpha, NROW(subfit$x))))
    } else {
      return(subfit$x %*% rbind(alpha, beta))
    }
  } else {
    contrasts_arg <- get_contrasts_arg_list(subfit$formula, newdata)
    x <- model.matrix(delete.response(terms(subfit$formula)), newdata,
                      contrasts.arg = contrasts_arg)
    if (is.null(beta)) {
      return(as.matrix(rep(alpha, NROW(x))))
    } else {
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
  gamm_struct <- model.matrix.gamm4(delete.response(terms(formula)),
                                    random = random, data = newdata)
  ranef <- lme4::ranef(fit$mer)
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
