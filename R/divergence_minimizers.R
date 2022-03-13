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
  # Define sdivmin(), the divergence minimizer for each draw s = 1, ..., S (and
  # perform other actions, if necessary):
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

# Use projpred's own implementation to fit non-multilevel non-additive
# submodels:
fit_glm_ridge_callback <- function(formula, data,
                                   projpred_var = matrix(nrow = nrow(data)),
                                   projpred_regul = 1e-4, contrasts = NULL,
                                   ...) {
  # Preparations:
  fr <- model.frame(formula, data = data)
  x <- model.matrix(fr, data = data, contrasts.arg = contrasts)
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
fit_glm_callback <- function(formula, family, projpred_var, projpred_regul,
                             ...) {
  if (family$family == "gaussian" && family$link == "identity") {
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
      scaled_data <- preprocess_data(data, projpred_formula_no_random)
      fit_gamm_callback(
        formula = formula,
        projpred_formula_no_random = projpred_formula_no_random,
        projpred_random = projpred_random,
        data = scaled_data,
        family = family,
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

# Use package "lme4" to fit multilevel submodels:
fit_glmer_callback <- function(formula, family,
                               control = control_callback(family), ...) {
  tryCatch({
    if (family$family == "gaussian" && family$link == "identity") {
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
        formula = formula, family = family, ...
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
        formula = formula,
        family = family,
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
        formula = formula,
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
             "automatically anymore.")
      }
      return(fit_glmer_callback(
        formula = formula,
        family = family,
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
    is_glmm <- inherits(fit, c("lmerMod", "glmerMod"))
    is_gam_gamm <- inherits(fit, c("gam", "gamm4"))
    if (is_gam_gamm && !is.null(newdata)) {
      newdata <- cbind(`(Intercept)` = rep(1, NROW(newdata)), newdata)
    }
    if (is_glmm) {
      return(predict(fit, newdata = newdata, allow.new.levels = TRUE))
    } else {
      return(predict(fit, newdata = newdata))
    }
  })
  return(do.call(cbind, prd_list))
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
    x <- model.matrix(delete.response(terms(subfit$formula)), newdata,
                      # TODO: Allow user-specified contrasts here:
                      contrasts.arg = NULL)
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
  gamm_struct <- model.matrix.gamm4(delete.response(terms(formula)),
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
