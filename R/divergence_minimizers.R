linear_mle <- function(formula, data, family, weights = NULL, var = 0, ...) {
  formula <- validate_response_formula(formula)
  if (inherits(formula, "formula")) {
    return(fit_glm_ridge_callback(
      formula, data = data, family = family, weights = weights, var = var,
      ...
    ))
  } else if (inherits(formula, "list")) {
    return(lapply(seq_along(formula), function(s) {
      fit_glm_ridge_callback(
        formula[[s]], data = data, family = family, weights = weights,
        var = var[, s, drop = FALSE], ...
      )
    }))
  } else {
    stop("The provided formula is neither a formula object nor a list")
  }
}

fit_glm_ridge_callback <- function(formula, data, family, weights, var = 0,
                                   regul = 1e-4, ...) {
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
    list(x = x, y = y, family = family, lambda = regul, weights = weights,
         obsvar = var),
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

# Use mgcv to fit the projection to the posterior draws for additive multilevel
# models.
additive_mle <- function(formula, data, family, weights = NULL, ...) {
  f <- split_formula_random_gamm4(formula)
  formula <- f$formula
  random <- f$random
  formula <- validate_response_formula(formula)
  if (inherits(formula, "formula")) {
    if (is.null(random)) {
      return(fit_gam_callback(
        formula, data = data, family = family, weights = weights, ...
      ))
    } else {
      return(fit_gamm_callback(
        formula, random = random, data = data, family = family,
        weights = weights, ...
      ))
    }
  } else if (inherits(formula, "list")) {
    if (is.null(random)) {
      return(lapply(
        formula, fit_gam_callback, data = data, family = family,
        weights = weights, ...
      ))
    } else {
      return(lapply(
        formula, fit_gamm_callback, random = random, data = data,
        family = family, weights = weights, ...
      ))
    }
  } else {
    stop("The provided formula is neither a formula object nor a list")
  }
}

# helper function of 'additive_mle'
#' @importFrom mgcv gam
fit_gam_callback <- function(formula, data, family, weights, ...) {
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

# helper function of 'additive_mle'
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

# Use lmer to fit the projection to the posterior draws for multilevel models.
linear_multilevel_mle <- function(formula, data, family, weights = NULL,
                                  var = 0, ...) {
  formula <- validate_response_formula(formula)
  if (inherits(formula, "formula")) {
    return(fit_glmer_callback(
      formula, data = data, family = family, weights = weights, var = var, ...
    ))
  } else if (inherits(formula, "list")) {
    return(lapply(seq_along(formula), function(s) {
      fit_glmer_callback(
        formula[[s]], data = data, family = family, weights = weights,
        var = var[, s, drop = FALSE], ...
      )
    }))
  } else {
    stop("The provided formula is neither a formula object nor a list")
  }
}

# helper function of 'linear_multilevel_mle'
fit_glmer_callback <- function(formula, data, family, weights,
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
      return(do.call(lme4::lmer, c(
        list(formula = formula, data = data, weights = weights,
             control = control),
        dot_args
      )))
    } else {
      # Exclude arguments from `...` which cannot be passed to lme4::glmer():
      dot_args <- list(...)
      dot_args <- dot_args[intersect(
        names(dot_args),
        methods::formalArgs(lme4::glmer)
      )]
      return(do.call(lme4::glmer, c(
        list(formula = formula, data = data, family = family, weights = weights,
             control = control),
        dot_args
      )))
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

# helper function of fit_glmer_callback to pass the proper kind of control
# options depending on the family
control_callback <- function(family, ...) {
  if (family$family == "gaussian" && family$link == "identity") {
    return(lme4::lmerControl(...))
  } else {
    return(lme4::glmerControl(...))
  }
}

# helper function for linear_multilevel_proj_predfun to only pass
# allow.new.levels if the fit is multilevel
predict_multilevel_callback <- function(fit, newdata = NULL) {
  if (inherits(fit, c("lmerMod", "glmerMod"))) {
    return(predict(fit, newdata = newdata, allow.new.levels = TRUE))
  } else {
    return(predict(fit, newdata = newdata))
  }
}

linear_multilevel_proj_predfun <- function(fit, newdata = NULL) {
  if (inherits(fit, "list")) {
    return(do.call(cbind, lapply(fit, function(fit) {
      predict_multilevel_callback(fit, newdata)
    })))
  } else {
    return(as.matrix(predict_multilevel_callback(fit, newdata)))
  }
}

linear_proj_predfun <- function(fit, newdata = NULL) {
  if (inherits(fit, "list")) {
    if (!is.null(newdata)) {
      return(do.call(cbind, lapply(fit, function(fit) {
        predict(fit, newdata = newdata)
      })))
    } else {
      return(do.call(cbind, lapply(fit, function(fit) {
        predict(fit)
      })))
    }
  }
  else {
    if (!is.null(newdata)) {
      return(predict(fit, newdata = newdata))
    } else {
      return(predict(fit))
    }
  }
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
