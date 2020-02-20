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

linear_mle <- function(formula, data, regul=NULL)
  lm(formula, data = data)

penalized_linear_mle <- function(formula, data, regul=1e-2) {
  formula <- validate_response_formula(formula)
  fit_penalized <- function(f) {
    all <- extract_terms_response(f)
    response <- data[, as.character(all$response)]
    ## no group terms allowed for penalized MLE
    covariates <- c(all$individual_terms, all$interaction_terms)
    if (length(covariates) == 0)
      linear_mle(f, data)
    else {
      covariates_formula <- as.formula(paste0("~ ", paste(covariates, collapse=" + ")))
      penalized::penalized(response=response, penalized=covariates_formula,
                           lambda2=regul, data=data, model="linear",
                           standardize=TRUE) # standard l2 regularisation
    }
  }
  if (inherits(formula, "formula"))
    return(fit_penalized(formula))
  else if (inherits(formula, "list"))
    return(lapply(formula, function(f) fit_penalized(f)))
  else
    stop("The provided formula is neither a formula nor a list")
}

#' Use lmer to fit the projection to the posterior draws for multilevel models.
#' Note that we don't use #' glmer because the target is a pseudo-Gaussian
#' transformation.
linear_multilevel_mle <- function(formula, data, regul=NULL) {
  formula <- validate_response_formula(formula)
  fit_lmer_callback <- function(f) {
    tryCatch(lme4::lmer(f, data = data),
             error=function(e) {
               if (grepl("No random effects", as.character(e)))
                 lm(f, data = data)
               else if (grepl("not positive definite", as.character(e)))
                 lme4::lmer(f, data = data,
                            control=lmerControl(optimizer="optimx",
                                                optCtrl=list( method="nlminb" )))
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

penalized_linear_proj_predfun <- function(fit, newdata=NULL) {
  if (inherits(fit, "list")) {
    return(do.call(cbind, lapply(fit, function(fit)
      penalized_linear_proj_predfun(fit, newdata=newdata))))
  } else if (inherits(fit, "penfit"))
    if (is.null(newdata))
      return(fit@fitted)
    else {
      covariates <- as.formula(paste0("~ ", paste(names(fit@penalized), collapse=" + ")))
      return(matrix(penalized::predict(fit, penalized=covariates, data=newdata), nrow(newdata), 2)[, 1])
    }
  else
    return(linear_proj_predfun(fit, newdata))
}
