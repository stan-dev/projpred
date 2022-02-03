# Latent projection predictive inference for generalized linear models
#
# Perform latent projection predictive variable selection for generalized
# linear models.

# TODO: I guess the documentation is not up-to-date anymore.
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
