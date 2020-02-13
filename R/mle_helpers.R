library(lme4)
library(penalized)
library(optimx)

fetch_data <- function(data, data_points=NULL, newdata=NULL) {
  if (is.null(data_points))
    if (is.null(newdata))
      return(data)
    else
      return(newdata)
  else if (is.null(newdata))
    return(data[data_points,, drop = FALSE])
  else
    return(newdata[data_points,, drop = FALSE])
}

linear_mle <- function(form, dat, regul=NULL)
  lm(form, data = dat)

penalized_linear_mle <- function(form, dat, regul=1e-2) {
  form <- validate_response_formula(form)
  fit_penalized <- function(f) {
    all <- extract_terms_response(f)
    response <- dat[, as.character(all$response)]
    ## no group terms allowed for penalized MLE
    covariates <- c(all$individual_terms, all$int_terms)
    if (length(covariates) == 0)
      linear_mle(f, dat)
    else {
      covariates_form <- as.formula(paste0("~ ", paste(covariates, collapse=" + ")))
      penalized::penalized(response=response, penalized=covariates_form,
                           lambda2=regul, data=dat, model="linear",
                           standardize=TRUE) # standard l2 regularisation
    }
  }
  if (class(form) == "formula")
    fit_penalized(form)
  else if (class(form) == "list")
    lapply(form, function(f) fit_penalized(f))
  else
    stop("The provided formula is neither a formula nor a list")
}

multilevel_mle <- function(form, dat, regul=NULL) {
  form <- validate_response_formula(form)
  fit_lmer_callback <- function(f) {
    tryCatch(lme4::lmer(f, data = dat),
             error=function(e) {
               if (!is.na(str_match(as.character(e),
                                    "No random effects")))
                 lm(f, data = dat)
               else if (!is.na(str_match(as.character(e),
                                         "not positive definite")))
                 lme4::lmer(f, data = dat,
                            control=lmerControl(optimizer="optimx",
                                                optCtrl=list( method="nlminb" )))
               else
                 e
             })
  }
  if (class(form) == "formula")
    fit_lmer_callback(form)
  else if (class(form) == "list")
    lapply(form, function(f) (fit_lmer_callback(f)))
  else
    stop("The provided formula is neither a formula object nor a list")
}

multilevel_proj_predfun <- function(fit, newdata=NULL) {
  if ("list" %in% class(fit)) {
    do.call(cbind, lapply(fit, function(fit)
      predict(fit, newdata=newdata, allow.new.levels=TRUE)))
  }
  else
    predict(fit, newdata=newdata, allow.new.levels=TRUE)
}

linear_proj_predfun <- function(fit, newdata=NULL) {
  if (!is.null(newdata))
    predict(fit, newdata=newdata)
  else
    predict(fit)
}

penalized_linear_proj_predfun <- function(fit, newdata=NULL) {
  if ("list" %in% class(fit)) {
    do.call(cbind, lapply(fit, function(fit)
      penalized_linear_proj_predfun(fit, newdata=newdata)))
  } else if (class(fit) == "penfit")
    if (is.null(newdata))
      fit@fitted
    else {
      covariates <- as.formula(paste0("~ ", paste(names(fit@penalized), collapse=" + ")))
      matrix(penalized::predict(fit, penalized=covariates, data=newdata), nrow(newdata), 2)[, 1]
    }
  else
    linear_proj_predfun(fit, newdata)
}
