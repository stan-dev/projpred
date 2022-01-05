#' Latent projection predictive inference for generalised linear models
#'
#' Perform latent projection predictive variable selection for generalised
#' linear models.

# utility function to perform kfolds validation over posterior samples
#' @export
.latent_cvfun <- function(folds, ...) {
  cvres <- brms::kfold(
    object,
    K = max(folds),
    save_fits = TRUE,
    folds = folds,
    ...
  )
  fits <- cvres$fits[, "fit"]
  return(fits)
}

#' @export
.latent_nlist <- function(...) {
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names)
    FALSE
  else
    nzchar(names(dots))
  if (all(has_name)) {
    return(dots)
  }
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

#' @export
.extract_latent_model_data <- function(object,
                               nlist = projpred:::.latent_nlist,
                               newdata = NULL,
                               wrhs = NULL,
                               orhs = NULL,
                               extract_y = TRUE) {
  if (!extract_y) {
    resp_form <- NULL
  } else {
    resp_form <- ~ .y
  }

  if (is.null(newdata)) {
    newdata <- data
  }

  if (is.null(wrhs) && !is.null(object) &&
      !is.null(object$weights) && length(object$weights) != 0) {
    wrhs <- ~ weights
    newdata <- cbind(newdata, weights = object$weights)
  }

  if (is.null(orhs) && !is.null(object) &&
      !is.null(object$offset) && length(object$offset) != 0) {
    orhs <- ~ offset
    newdata <- cbind(newdata, offset = object$offset)
  }

  args <- nlist(object, newdata, wrhs, orhs, resp_form)
  return(do_call(projpred:::.extract_model_data, args))
}

# use the sample predictive function as for the refrence model
#' @export
.latent_ref_predfun <- function(fit, newdata = NULL) {
  return(t(posterior_linpred(
    fit, newdata = newdata, transform = FALSE
  )))
}

#' Extract the posterior draws of the latent predictor
#'
#' @param fit A reference model written
#' @param data The data on which the reference model was trained
#' @param ref_predfun The predictive function of the reference model
#' @return An array of the posterior draws from the latent predictor
#' @export
extract_eta <- function(fit, data) {
  # explicit the predictive function for the refrence model posterior latent
  # predictor
  eta_post_draws <-
    t(brms::posterior_linpred(fit, newdata = data, transform = FALSE))

  # TODO: fit a generalised Pareto distribution to the tail of the latent
  # predictor's posterior distribution, and identify which distribution
  # is optimal for future fitting

  # TODO: investigate performing an inverse-logit transformation on the latent
  # predictor to restrict its support, ensure that there exists a mean, and
  # approximate with a Beta distribution. If this doesn't work, consider
  # approximating with a mixture of Betas.

  recommended_latent_family <- ""
  if (TRUE) {
    recommended_latent_family <- "Gaussian"
    brms_recommended_latent_family_command <-
      "latent_family = brms::brmsfamily('gaussian')"
  } else if (FALSE) {
    recommended_latent_family <- "mixture of Betas"
    brms_recommended_latent_family_command <-
      ""
  }

  message(
    paste(
      "To continue with latent projection predictive inference, kindly pass",
      "the returned `eta_post_draws` object to the `fit_latent` method.\n"
    )
  )
  message(
    paste(
      "The recommended distribution for the latent predictor is",
      recommended_latent_family,
      "since it has infinite support and two finite moments. Kindly pass",
      brms_recommended_latent_family_command,
      "to `fit_latent` method for best results (Catalina et al., 2021)."
    )
  )
  return(eta_post_draws)
}

#' Fit the latent reference model
#'
#' This function takes in a reference model, the posterior draws of its latent
#' predictor, and the data on which the reference model was trained to fit a
#' reference model in terms of the latent predictor using a Gaussianity
#' assumption.
#'
#' @param fit A reference model written
#' @param eta_post_draws The posterior predictive draws of the latent predictor from
#' the reference model
#' @param data The data on which the reference model was trained
#' @param latent_div_minimizer The divergence minimising formula to be used
#' on the latent predictor
#' @param latent_proj_predfun The projection predictive function
#' @param dis_latent An array of the dispertion parameter for the latent predictor
#' @param extract_model_data A function to extract data from the model
#' @param ref_predfun The predictive function of the reference model
#' @param cv_fun A function for cross-validation
#' @return A latent predictor reference model
#' @export
fit_latent <-
  function(fit,
           eta_post_draws,
           data,
           latent_family = brms::brmsfamily("gaussian"),
           latent_div_minimizer = projpred:::linear_mle,
           latent_proj_predfun = projpred:::linear_proj_predfun,
           dis_latent = rep(1, 4000), # fixed and uniform
           extract_model_data = .extract_latent_model_data,
           ref_predfun = .latent_ref_predfun,
           cv_fun = .latent_cvfun) {
    # update the reference model's formula to fit to the latent predictor
    latent_formula <- formula(update(fit$formula, ".y ~ ."))
    # add latent predictor to input data
    data[[".y"]] <- rowMeans(eta_post_draws)
    # define new reference model over latent parameter space with Gaussian family
    latent_ref <- init_refmodel(
      fit,
      data = data,
      formula = latent_formula,
      family = latent_family,
      ref_predfun = ref_predfun,
      div_minimizer = latent_div_minimizer,
      proj_predfun = latent_proj_predfun,
      dis = dis_latent,
      extract_model_data = extract_model_data,
      cvfun = cv_fun
    )
    message(
      paste(
        "To proceed to variable selection in the latent space, kindly pass",
        "the returned `latent_ref` object to the `varsel` or `cv_varsel` method."
      )
    )
    return(latent_ref)
  }
