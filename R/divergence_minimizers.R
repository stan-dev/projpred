fit_glm_ridge_callback <- function(formula, data, family, weights = NULL,
                                   projpred_var = 0, projpred_regul = 1e-4,
                                   ...) {
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
    list(x = x, y = y, family = family, lambda = projpred_regul,
         weights = weights, obsvar = projpred_var),
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

# Use package "mgcv" to fit submodels for additive reference models. Use package
# "gamm4" to fit submodels for additive multilevel reference models:
fit_gam_gamm_callback <- function(formula, data, family, weights = NULL,
                                  projpred_random, ...) {
  stopifnot(is.null(projpred_random) || inherits(projpred_random, "formula"))
  if (is.null(projpred_random)) {
    return(fit_gam_callback(
      formula, data = data, family = family, weights = weights, ...
    ))
  } else {
    return(fit_gamm_callback(
      formula, random = projpred_random, data = data, family = family,
      weights = weights, ...
    ))
  }
}

# Helper function for fit_gam_gamm_callback():
#' @importFrom mgcv gam
fit_gam_callback <- function(formula, data, family, weights = NULL, ...) {
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
    list(formula = formula, data = data, family = family, weights = weights),
    dot_args
  )))))
}

# Helper function for fit_gam_gamm_callback():
#' @importFrom gamm4 gamm4
fit_gamm_callback <- function(formula, random, data, family, weights = NULL,
                              control = control_callback(family), ...) {
  # make sure correct 'weights' can be found
  environment(formula) <- environment()
  # Exclude arguments from `...` which cannot be passed to gamm4::gamm4():
  dot_args <- list(...)
  dot_args <- dot_args[intersect(
    names(dot_args),
    union(union(methods::formalArgs(gamm4),
                methods::formalArgs(lme4::lFormula)),
          methods::formalArgs(lme4::glFormula))
  )]
  fit <- suppressMessages(suppressWarnings(tryCatch({
    do.call(gamm4, c(
      list(formula = formula, random = random, data = data, family = family,
           weights = weights, control = control),
      dot_args
    ))
  }, error = function(e) {
    if (grepl("not positive definite", as.character(e))) {
      scaled_data <- preprocess_data(data, formula)
      fit_gamm_callback(
        formula, random = random, data = scaled_data, weights = weights,
        family = family,
        control = control_callback(family,
                                   optimizer = "optimx",
                                   optCtrl = list(method = "nlminb")),
        ...
      )
    } else {
      stop(e)
    }
  })))

  fit$random <- random
  fit$formula <- formula
  class(fit) <- c("gamm4")
  return(fit)
}

# Use package "lme4" to fit submodels for multilevel reference models (with a
# fallback to "projpred"'s own implementation for fitting non-multilevel (and
# non-additive) submodels):
fit_glmer_callback <- function(formula, data, family, weights = NULL,
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
        list(formula = formula, data = data, weights = weights,
             control = control),
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
        list(formula = formula, data = data, family = family, weights = weights,
             control = control),
        dot_args
      )))))
    }
  }, error = function(e) {
    if (grepl("No random effects", as.character(e))) {
      return(fit_glm_ridge_callback(
        formula, data = data, family = family, weights = weights, ...
      ))
    } else if (grepl("not positive definite", as.character(e))) {
      return(fit_glmer_callback(
        formula, data = data, family = family, weights = weights,
        control = control_callback(family,
                                   optimizer = "optimx",
                                   optCtrl = list(method = "nlminb")),
        ...
      ))
    } else if (grepl("PIRLS step-halvings", as.character(e))) {
      return(fit_glmer_callback(
        formula, data = data, family = family, weights = weights,
        control = control, nAGQ = 20L, ...
      ))
    } else if (grepl("pwrssUpdate did not converge in \\(maxit\\) iterations",
                     as.character(e))) {
      return(fit_glmer_callback(
        formula, data = data, family = family, weights = weights,
        control = control_callback(family, tolPwrss = 1e-6), ...
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

linear_multilevel_proj_predfun <- function(fit, newdata = NULL) {
  return(do.call(cbind, lapply(fit, function(fit) {
    # Only pass argument `allow.new.levels` to the predict() generic if the fit
    # is multilevel:
    if (inherits(fit, c("lmerMod", "glmerMod"))) {
      return(predict(fit, newdata = newdata, allow.new.levels = TRUE))
    } else {
      return(predict(fit, newdata = newdata))
    }
  })))
}

linear_proj_predfun <- function(fit, newdata = NULL) {
  return(do.call(cbind, lapply(fit, function(fit) {
    predict(fit, newdata = newdata)
  })))
}

additive_proj_predfun <- function(fit, newdata = NULL) {
  if (!is.null(newdata)) {
    newdata <- cbind(`(Intercept)` = rep(1, NROW(newdata)), newdata)
  }
  return(linear_multilevel_proj_predfun(fit, newdata))
}

## FIXME: find a way that allows us to remove this
predict.subfit <- function(subfit, newdata = NULL) {
  beta <- subfit$beta
  alpha <- subfit$alpha
  x <- subfit$x
  if (is.null(newdata)) {
    if (is.null(beta)) {
      return(as.matrix(rep(alpha, NROW(subfit$x))))
    } else {
      return(x %*% rbind(alpha, beta))
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
