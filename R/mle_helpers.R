fetch_data <- function(data, obs=NULL, newdata=NULL) {
  if (is.null(obs))
    if (is.null(newdata))
      return(data)
    else
      return(newdata)
  else if (is.null(newdata))
    return(data[obs,, drop = FALSE])
  else
    return(newdata[obs,, drop = FALSE])
}

linear_mle <- function(formula, data, weights=NULL, regul=NULL)
  lm(formula, data = data, weights = weights)

#' Use lmer to fit the projection to the posterior draws for multilevel models.
#' Note that we don't use glmer because the target is a pseudo-Gaussian
#' transformation.
linear_multilevel_mle <- function(formula, data, weights = NULL, regul=NULL) {
  formula <- validate_response_formula(formula)
  fit_lmer_callback <- function(f) {
    tryCatch(lme4::lmer(f, data = data, weigts = weights),
             error=function(e) {
               if (grepl("No random effects", as.character(e)))
                 lm(f, data = data, weights = weights)
               else if (grepl("not positive definite", as.character(e)))
                 lme4::lmer(f, data = data, weights = weights,
                            control = lmerControl(optimizer="optimx",
                                                  optCtrl = list(method="nlminb")))
               else
                 e
             })
  }
  if (inherits(formula, "formula"))
    return(fit_lmer_callback(formula))
  else if (inherits(formula, "list"))
    return(lapply(formula, function(f) (fit_lmer_callback(f))))
  else
    stop("The provided formula is neither a formula object nor a list")
}

linear_multilevel_proj_predfun <- function(fit, newdata=NULL) {
  if (inherits(fit, "list")) {
    return(do.call(cbind, lapply(fit, function(fit)
      predict(fit, newdata=newdata, allow.new.levels=TRUE))))
  } else
    return(predict(fit, newdata=newdata, allow.new.levels=TRUE))
}

linear_proj_predfun <- function(fit, newdata=NULL) {
  if (!is.null(newdata))
    return(predict(fit, newdata=newdata))
  else
    return(predict(fit))
}
