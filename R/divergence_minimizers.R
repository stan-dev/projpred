fetch_data <- function(data, obs = NULL, newdata = NULL) {
  if (is.null(obs)) {
    if (is.null(newdata)) {
      return(data)
    } else {
      return(newdata)
    }
  } else if (is.null(newdata)) {
    return(data[obs, , drop = FALSE])
  } else {
    return(newdata[obs, , drop = FALSE])
  }
}

linear_mle <- function(formula, data, family, weights = NULL, regul = NULL) {
  formula <- validate_response_formula(formula)
  if (inherits(formula, "formula")) {
    return(fit_glm_callback(formula, data, family, weights))
  } else if (inherits(formula, "list")) {
    return(lapply(formula, fit_glm_callback, data, family, weights))
  } else {
    stop("The provided formula is neither a formula object nor a list")
  }
}

# helper function of 'linear_mle'
fit_glm_callback <- function(formula, data, family, weights) {
  # make sure correct 'weights' can be found
  environment(formula) <- environment()
  if (family$family == "gaussian" && family$link == "identity") {
    return(lm(formula, data = data, weights = weights))
  } else {
    return(glm(formula, data = data, family = family, weights = weights))
  }
}

#' Use lmer to fit the projection to the posterior draws for multilevel models.
#' Note that we don't use glmer because the target is a pseudo-Gaussian
#' transformation.
linear_multilevel_mle <- function(formula, data, family, weights = NULL,
                                  regul = NULL) {
  formula <- validate_response_formula(formula)
  if (inherits(formula, "formula")) {
    return(fit_glmer_callback(formula, data, family, weights))
  } else if (inherits(formula, "list")) {
    return(lapply(formula, fit_glmer_callback, data, family, weights))
  } else {
    stop("The provided formula is neither a formula object nor a list")
  }
}

# helper function of 'linear_multilevel_mle'
fit_glmer_callback <- function(formula, data, family, weights,
                               control = control_callback(family)) {
  ## make sure correct 'weights' can be found
  environment(formula) <- environment()
  tryCatch({
      if (family$family == "gaussian" && family$link == "identity") {
        return(lme4::lmer(formula,
          data = data, weights = weights,
          control = control
        ))
      } else {
        return(lme4::glmer(formula,
          data = data, family = family, weights = weights,
          control = control
        ))
      }
    },
    error = function(e) {
      if (grepl("No random effects", as.character(e))) {
        return(fit_glm_callback(formula,
          data = data, family = family,
          weights = weights
        ))
      } else if (grepl("not positive definite", as.character(e))) {
        return(fit_glmer_callback(formula,
          data = data, weights = weights, family = family,
          control = control_callback(family,
            optimizer = "optimx",
            optCtrl = list(method = "nlminb")
          )
        ))
      } else {
        return(e)
      }
    }
  )
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
predict_multilevel_callback <- function(fit, newdata = NULL, weights = NULL) {
  if (inherits(fit, "lmerMod")) {
    return(predict(fit,
      newdata = newdata, allow.new.levels = TRUE,
      weights = weights
    ))
  } else {
    return(predict(fit, newdata = newdata, weights = weights))
  }
}

linear_multilevel_proj_predfun <- function(fit, newdata = NULL,
                                           weights = NULL) {
  if (is.null(weights)) {
    weights <- 1
  }
  if (inherits(fit, "list")) {
    return(do.call(cbind, lapply(fit, function(fit) {
      predict_multilevel_callback(fit, newdata, weights)
    })))
  } else {
    return(predict_multilevel_callback(fit, newdata, weights))
  }
}

linear_proj_predfun <- function(fit, newdata = NULL, weights = NULL) {
  if (is.null(weights)) {
    weights <- 1
  }
  if (inherits(fit, "list")) {
    if (!is.null(newdata)) {
      return(do.call(cbind, lapply(fit, function(fit) {
        predict(fit, newdata = newdata, weights = weights)
      })))
    } else {
      return(do.call(cbind, lapply(fit, function(fit) {
        predict(fit)
      })))
    }
  }
  else {
    if (!is.null(newdata)) {
      return(predict(fit, newdata = newdata, weights = weights))
    } else {
      return(predict(fit))
    }
  }
}
