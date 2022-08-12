# Family-specific helper functions
#
# `extend_family(family)` returns a [`family`] object augmented with auxiliary
# functions that are needed for computing KL-divergence, log predictive density,
# dispersion projection, etc.
#
# Missing: Quasi-families are not implemented. If dis_gamma is the correct shape
# parameter for projected Gamma regression, everything should be OK for gamma.

#' Extend a family
#'
#' This function adds some internally required elements to a [`family`] object.
#' It is called internally by [init_refmodel()], so you will rarely need to call
#' it yourself.
#'
#' @param family A [`family`] object.
#' @param latent A single logical value indicating whether to use the latent
#'   projection (`TRUE`) or not (`FALSE`). Note that setting `latent = TRUE`
#'   causes all augmented-data arguments to be ignored.
#' @param latent_y_unqs Only relevant for a latent projection where the original
#'   response space has finite support (i.e., the original response values may
#'   be regarded as categories), in which case this needs to be the character
#'   vector of unique response values. See also section "Latent projection"
#'   below.
#' @param latent_ilink Only relevant for the latent projection, in which case
#'   this needs to be the inverse-link function (supplied as a function, not a
#'   character string, for example). Can also be `NULL`, but then downstream
#'   functions will have limited functionality (a message thrown by
#'   [extend_family()] will state what exactly won't be available). See also
#'   section "Latent projection" below.
#' @param latent_llOrig Only relevant for the latent projection, in which case
#'   this needs to be the function (supplied as a function, not a character
#'   string, for example) computing response-scale (not latent-scale)
#'   log-likelihood values. If `!is.null(latent_y_unqs)`, then `latent_llOrig`
#'   can be `NULL`, in which case an internal default will be used. Can also be
#'   `NULL` in all other cases, but then downstream functions will have limited
#'   functionality (a message thrown by [extend_family()] will state what
#'   exactly won't be available). See also section "Latent projection" below.
#' @param latent_ppdOrig Only relevant for the latent projection, in which case
#'   this needs to be the function (supplied as a function, not a character
#'   string, for example) sampling response values given latent predictors that
#'   have been transformed to response scale using `latent_ilink`. If
#'   `!is.null(latent_y_unqs)`, then `latent_ppdOrig` can be `NULL`, in which
#'   case an internal default will be used. Can also be `NULL` in all other
#'   cases, but then downstream functions will have limited functionality (a
#'   message thrown by [extend_family()] will state what exactly won't be
#'   available). See also section "Latent projection" below. Note that although
#'   this function has the abbreviation "PPD" in its name (which stands for
#'   "posterior predictive distribution"), \pkg{projpred} currently only uses it
#'   in [proj_predict()], i.e., for sampling from what would better be termed
#'   posterior-projection predictive distribution (PPPD).
#' @param augdat_y_unqs Only relevant for augmented-data projection, in which
#'   case this needs to be the character vector of unique response values (which
#'   will be assigned to `family$cats` internally) or may be left at `NULL` if
#'   `family$cats` is already non-`NULL`. See also section "Augmented-data
#'   projection" below.
#' @param augdat_link Only relevant for augmented-data projection, in which case
#'   this needs to be the link function (supplied as a function, not a character
#'   string, for example). Use `NULL` for the traditional (non-augmented-data)
#'   projection. See also section "Augmented-data projection" below.
#' @param augdat_ilink Only relevant for augmented-data projection, in which
#'   case this needs to be the inverse-link function (supplied as a function,
#'   not a character string, for example). Use `NULL` for the traditional
#'   (non-augmented-data) projection. See also section "Augmented-data
#'   projection" below.
#' @param augdat_args_link Only relevant for augmented-data projection, in which
#'   case this may be a named `list` of arguments to pass to the function
#'   supplied to `augdat_link`.
#' @param augdat_args_ilink Only relevant for augmented-data projection, in
#'   which case this may be a named `list` of arguments to pass to the function
#'   supplied to `augdat_ilink`.
#' @param ... Ignored (exists only to swallow up further arguments which might
#'   be passed to this function).
#'
#' @details
#'
#' # Augmented-data projection
#'
#' In the following, \eqn{N}, \eqn{C_{\mathrm{cat}}}{C_cat},
#' \eqn{C_{\mathrm{lat}}}{C_lat}, \eqn{S_{\mathrm{ref}}}{S_ref}, and
#' \eqn{S_{\mathrm{prj}}}{S_prj} from help topic [refmodel-init-get] are used.
#' Furthermore, let \eqn{S} denote either \eqn{S_{\mathrm{ref}}}{S_ref} or
#' \eqn{S_{\mathrm{prj}}}{S_prj}, whichever is appropriate in the context where
#' it is used.
#'
#' Then, as their first input, the functions supplied to arguments `augdat_link`
#' and `augdat_ilink` have to accept:
#' * For `augdat_link`: an \eqn{S \times N \times C_{\mathrm{cat}}}{S x N x
#' C_cat} array containing the probabilities for the response categories. The
#' order of the response categories is the same as in `family$cats` (see
#' argument `augdat_y_unqs`).
#' * For `augdat_ilink`: an \eqn{S \times N \times C_{\mathrm{lat}}}{S x N x
#' C_lat} array containing the linear predictors.
#'
#' The return value of these functions needs to be:
#' * For `augdat_link`: an \eqn{S \times N \times C_{\mathrm{lat}}}{S x N x
#' C_lat} array containing the linear predictors.
#' * For `augdat_ilink`: an \eqn{S \times N \times C_{\mathrm{cat}}}{S x N x
#' C_cat} array containing the probabilities for the response categories. The
#' order of the response categories has to be the same as in `family$cats` (see
#' argument `augdat_y_unqs`).
#'
#' For the augmented-data projection, the response vector resulting from
#' `extract_model_data` (see [init_refmodel()]) is coerced to a `factor` (using
#' [as.factor()]). Inside of [init_refmodel()], the levels of this `factor` have
#' to be identical to `family$cats`. Everywhere else, they have to be a subset
#' of `family$cats`. See argument `augdat_y_unqs` for how to control
#' `family$cats`.
#'
#' For ordinal \pkg{brms} families, be aware that the submodels (onto which the
#' reference model is projected) currently have the following restrictions:
#' * The discrimination parameter `disc` is not supported (i.e., it is a
#' constant with value 1).
#' * The thresholds are `"flexible"` (see [brms::brmsfamily()]).
#' * The thresholds do not vary across the levels of a `factor`-like variable
#' (see argument `gr` of [brms::resp_thres()]).
#' * The `"probit_approx"` link is replaced by `"probit"`.
#'
#' For the [brms::categorical()] family, be aware that:
#' * For multilevel submodels, the group-level effects are allowed to be
#' correlated between different response categories.
#' * For multilevel submodels, \pkg{mclogit} versions < 0.9.4 may throw the
#' error \code{'a' (<number> x 1) must be square}. Updating \pkg{mclogit} to a
#' version >= 0.9.4 should fix this.
#'
#' # Latent projection
#'
#' **TODO (latent)**: Add details about arguments and return values of the
#' functions required for the latent projection.
#' If `!is.null(latent_y_unqs)`, then `latent_ilink()` needs to return
#' *probabilities* (those for the response categories given in `latent_y_unqs`).
#' The default for `wdraws_ref` needs to be `rep(1, length(cl_ref))`, but
#' \pkg{projpred} tries to specify all arguments (i.e., not to rely on defaults)
#' when calling these user-specified functions. Function [cl_agg()] might be
#' helpful for constructing `latent_ilink`. In fact, argument `wdraws_ref` of
#' `latent_ilink` is needed (i.e., the weights passed to this argument are
#' nonconstant) only in case of [cv_varsel()] with `cv_method = "LOO"` and
#' `validate_search = TRUE`. In that case, the weights passed to this argument
#' are the PSIS-LOO CV weights for one observation. Note that although argument
#' `wdraws_ref` has the suffix `_ref`, `wdraws_ref` does not necessarily obtain
#' weights for the *initial* reference model's posterior draws: In case of
#' [cv_varsel()] with `cv_method = "kfold`, these weights may refer to one of
#' the \eqn{K} reference model re-fits (but in that case, they are constant
#' anyway).
#'
#' @return The [`family`] object extended in the way needed by \pkg{projpred}.
#'
#' @export
extend_family <- function(family,
                          latent = FALSE,
                          latent_y_unqs = NULL,
                          latent_ilink = NULL,
                          latent_llOrig = NULL,
                          latent_ppdOrig = NULL,
                          augdat_y_unqs = NULL,
                          augdat_link = NULL,
                          augdat_ilink = NULL,
                          augdat_args_link = list(),
                          augdat_args_ilink = list(),
                          ...) {
  if (.has_family_extras(family)) {
    # If the family was already extended using this function, then return as-is:
    return(family)
  }
  if (latent) {
    familyOrig_tmp <- family$family
    linkOrig_tmp <- family$link
    family <- gaussian()
    family$familyOrig <- familyOrig_tmp
    family$linkOrig <- linkOrig_tmp
  }
  aug_data <- !is.null(augdat_link) && !is.null(augdat_ilink) && !latent
  if (!aug_data) {
    extend_family_specific <- paste0("extend_family_", tolower(family$family))
    if (!exists(extend_family_specific, mode = "function")) {
      stop("Family '", family$family, "' is not supported by projpred.")
    }
    extend_family_specific <- get(extend_family_specific, mode = "function")
    family <- extend_family_specific(family)

    # If `family$cats` weren't `NULL`, then downstream code in projpred would
    # have to be adapted:
    stopifnot(is.null(family$cats))

    if (latent) {
      if (!is.null(latent_y_unqs)) {
        family$cats <- latent_y_unqs
        if (is.null(latent_llOrig)) {
          latent_llOrig <- latent_llOrig_cats
        }
        if (is.null(latent_ppdOrig)) {
          latent_ppdOrig <- latent_ppdOrig_cats
        }
      }
      if (is.null(latent_ilink)) {
        message("`latent_ilink` is `NULL`, so cv_varsel() with ",
                "`cv_method = \"LOO\"` won't be usable. Furthermore, ",
                "predict.refmodel(), summary.vsel(), print.vsel(), ",
                "plot.vsel(), suggest_size.vsel(), proj_linpred(), and ",
                "proj_predict() won't work on response scale (only on latent ",
                "scale).")
      }
      if (is.null(latent_llOrig)) {
        message("`latent_llOrig` is `NULL`, so cv_varsel() with ",
                "`cv_method = \"LOO\"` won't be usable. Furthermore, some ",
                "features of predict.refmodel(), summary.vsel(), ",
                "print.vsel(), plot.vsel(), suggest_size.vsel(), and ",
                "proj_linpred() won't work on response scale (only on latent ",
                "scale).")
      }
      if (is.null(latent_ppdOrig)) {
        message("`latent_ppdOrig` is `NULL`, so proj_predict() won't work on ",
                "response scale (only on latent scale).")
      }
      family$lat2resp_possible <- !is.null(latent_ilink) &&
        !is.null(latent_llOrig)
      family$ppdOrig_possible <- !is.null(latent_ilink) &&
        !is.null(latent_ppdOrig)
      family$latent_ilink <- latent_ilink
      family$latent_llOrig <- latent_llOrig
      family$latent_ppdOrig <- latent_ppdOrig
    }
    family$for_latent <- latent
    family$for_augdat <- FALSE
  } else {
    if (!is.null(augdat_y_unqs)) {
      family$cats <- augdat_y_unqs
    } else if (is.null(family$cats)) {
      stop("Please supply argument `augdat_y_unqs` or ensure that ",
           "`family$cats` is not `NULL`.")
    }

    # Checks for special 'brms' features (currently only necessary for argument
    # `refcat` of the brms::categorical() family):
    if (!(is.null(family$refcat) ||
          identical(family$refcat, utils::head(family$cats, 1)))) {
      stop("Currently, the first category must be the reference category.")
    }

    family$linkfun <- function(mu) {
      mu_arr <- augmat2arr(mu, margin_draws = 1)
      eta_arr <- do.call(augdat_link, c(list(mu_arr), augdat_args_link))
      return(arr2augmat(eta_arr, margin_draws = 1))
    }
    family$linkinv <- function(eta) {
      eta_arr <- augmat2arr(eta, margin_draws = 1)
      mu_arr <- do.call(augdat_ilink, c(list(eta_arr), augdat_args_ilink))
      return(arr2augmat(mu_arr, margin_draws = 1))
    }
    family$kl_ptwise <- function(mu_ref, mu_sub, w_obs = 1) {
      mu_ref_arr <- augmat2arr(mu_ref)
      mu_sub_arr <- augmat2arr(mu_sub)
      stopifnot(identical(dim(mu_ref_arr), dim(mu_sub_arr)))
      n_obs <- dim(mu_ref_arr)[1]
      n_prjdraws <- dim(mu_ref_arr)[3]
      if (length(w_obs) == 0) {
        w_obs <- rep(1, n_obs)
      } else if (length(w_obs) == 1) {
        w_obs <- rep(w_obs, n_obs)
      } else if (length(w_obs) != n_obs) {
        stop("Argument `w_obs` needs to be of length 0, 1, or `n_obs`.")
      }
      return(do.call(rbind, lapply(seq_len(n_obs), function(i_obs) {
        do.call(c, lapply(seq_len(n_prjdraws), function(i_prjdraws) {
          prbs_ref <- mu_ref_arr[i_obs, , i_prjdraws]
          prbs_sub <- mu_sub_arr[i_obs, , i_prjdraws]
          # Assign some nonzero value to have a finite log() value (the specific
          # value doesn't matter for `prbs_ref` because of the multiplication
          # with zero):
          prbs_ref[prbs_ref == 0] <- .Machine$double.eps
          # Assign some nonzero value to have a finite log() value (for
          # `prbs_sub`, the specific value matters, in contrast to `prbs_ref`):
          prbs_sub[prbs_sub == 0] <- .Machine$double.eps
          # In analogy to `extend_family(binomial())$kl()`,
          # `extend_family(poisson())$kl()`, and
          # `extend_family(gaussian())$kl()`, multiply by the observation
          # weight:
          return(w_obs[i_obs] * sum(prbs_ref * (log(prbs_ref) - log(prbs_sub))))
        }))
      })))
    }
    family$kl <- function(pref, data, psub) {
      data$weights <- data$weights / sum(data$weights)
      return(colSums(
        family$kl_ptwise(mu_ref = pref$mu,
                         mu_sub = psub$mu,
                         w_obs = data$weights)
      ))
    }
    family$dis_fun <- function(pref, psub, wobs = 1) {
      return(rep(NA, ncol(pref$mu)))
    }
    family$predvar <- function(mu, dis, wsample = 1) {
      return(rep(NA, NROW(mu)))
    }
    family$ll_fun <- function(mu, dis = NULL, y, weights = 1) {
      mu_arr <- augmat2arr(mu)
      return(ll_cats(mu_arr, y = y, wobs = weights))
    }
    family$ppd <- function(mu, dis, weights = 1) {
      mu_arr <- augmat2arr(augvec2augmat(mu))
      return(ppd_cats(mu_arr, wobs = weights, return_vec = TRUE))
    }
    family$for_latent <- FALSE
    family$for_augdat <- TRUE
  }
  family$is_extended <- TRUE
  return(family)
}

extend_family_binomial <- function(family) {
  kl_dev <- function(pref, data, psub) {
    data$weights <- data$weights / sum(data$weights)
    data$weights <- rep(data$weights, ncol(pref$mu))
    colSums(family$dev.resids(pref$mu, psub$mu, data$weights)) / 2
  }
  dis_na <- function(pref, psub, wobs = 1) {
    rep(NA, ncol(pref$mu))
  }
  predvar_na <- function(mu, dis, wsample = 1) {
    rep(NA, NROW(mu))
  }
  ll_binom <- function(mu, dis, y, weights = 1) {
    y <- as.matrix(y)
    dbinom(y, weights, mu, log = TRUE)
  }
  dev_binom <- function(mu, y, weights = 1, dis = NULL) {
    if (NCOL(y) < NCOL(mu)) {
      y <- matrix(y, nrow = length(y), ncol = NCOL(mu))
    }
    -2 * weights * (y * log(mu) + (1 - y) * log(1 - mu))
  }
  ppd_binom <- function(mu, dis, weights = 1) {
    rbinom(length(mu), weights, mu)
  }
  initialize_binom <- expression({
    if (NCOL(y) == 1) {
      if (is.factor(y)) {
        y <- y != levels(y)[1L]
      }
      n <- rep.int(1, nobs)
      y[weights == 0] <- 0
      if (any(y < 0 | y > 1)) {
        stop("y values must be 0 <= y <= 1")
      }
      mustart <- (weights * y + 0.5) / (weights + 1)
      m <- weights * y
      if ("binomial" == "binomial" && any(abs(m - round(m)) >
                                          0.001)) {
        ### Deactivated because in general, this will be the case in 'projpred':
        # warning(gettextf("non-integer #successes in a %s glm!",
        #                  "binomial"), domain = NA)
        ###
      }
    }
    else if (NCOL(y) == 2) {
      if ("binomial" == "binomial" && any(abs(y - round(y)) >
                                          0.001)) {
        warning(gettextf("non-integer counts in a %s glm!",
                         "binomial"), domain = NA)
      }
      n <- (y1 <- y[, 1L]) + y[, 2L]
      y <- y1 / n
      if (any(n0 <- n == 0)) {
        y[n0] <- 0
      }
      weights <- weights * n
      mustart <- (n * y + 0.5) / (n + 1)
    } else {
      stop(gettextf(paste("for the '%s' family, y must be a vector of 0 and",
                          "1's\nor a 2 column matrix where col 1 is no.",
                          "successes and col 2 is no. failures"),
                    "binomial"), domain = NA)
    }
  })

  family$initialize <- initialize_binom
  family$kl <- kl_dev
  family$dis_fun <- dis_na
  family$predvar <- predvar_na
  family$ll_fun <- ll_binom
  family$deviance <- dev_binom
  family$ppd <- ppd_binom

  return(family)
}

extend_family_poisson <- function(family) {
  kl_dev <- function(pref, data, psub) {
    data$weights <- data$weights / sum(data$weights)
    data$weights <- rep(data$weights, ncol(pref$mu))
    colSums(family$dev.resids(pref$mu, psub$mu, data$weights)) / 2
  }
  dis_na <- function(pref, psub, wobs = 1) {
    rep(NA, ncol(pref$mu))
  }
  predvar_na <- function(mu, dis, wsample = 1) {
    rep(NA, NROW(mu))
  }
  ll_poiss <- function(mu, dis, y, weights = 1) {
    y <- as.matrix(y)
    weights * dpois(y, mu, log = TRUE)
  }
  dev_poiss <- function(mu, y, weights = 1, dis = NULL) {
    if (NCOL(y) < NCOL(mu)) {
      y <- matrix(y, nrow = length(y), ncol = NCOL(mu))
    }
    -2 * weights * (y * log(mu) - mu)
  }
  ppd_poiss <- function(mu, dis, weights = 1) {
    rpois(length(mu), mu)
  }

  family$kl <- kl_dev
  family$dis_fun <- dis_na
  family$predvar <- predvar_na
  family$ll_fun <- ll_poiss
  family$deviance <- dev_poiss
  family$ppd <- ppd_poiss

  return(family)
}

extend_family_gaussian <- function(family) {
  kl_gauss <- function(pref, data, psub) {
    data$weights <- data$weights / sum(data$weights)
    colSums(data$weights * (psub$mu - pref$mu)^2)
  } # not the actual KL but reasonable surrogate..
  dis_gauss <- function(pref, psub, wobs = 1) {
    sqrt(colSums(wobs / sum(wobs) * (pref$var + (pref$mu - psub$mu)^2)))
  }
  predvar_gauss <- function(mu, dis, wsample = 1) {
    wsample <- wsample / sum(wsample)
    mu_mean <- mu %*% wsample
    mu_var <- mu^2 %*% wsample - mu_mean^2
    as.vector(sum(wsample * dis^2) + mu_var)
  }
  ll_gauss <- function(mu, dis, y, weights = 1) {
    y <- as.matrix(y)
    dis <- matrix(rep(dis, each = length(y)), ncol = NCOL(mu))
    weights * dnorm(y, mu, dis, log = TRUE)
  }
  dev_gauss <- function(mu, y, weights = 1, dis = NULL) {
    if (is.null(dis)) {
      dis <- 1
    } else {
      dis <- matrix(rep(dis, each = length(y)), ncol = NCOL(mu))
    }
    if (NCOL(y) < NCOL(mu)) {
      y <- matrix(y, nrow = length(y), ncol = NCOL(mu))
    }
    -2 * weights * (-0.5 / dis^2 * (y - mu)^2 - log(dis))
  }
  ppd_gauss <- function(mu, dis, weights = 1) {
    rnorm(length(mu), mu, dis)
  }

  family$kl <- kl_gauss
  family$dis_fun <- dis_gauss
  family$predvar <- predvar_gauss
  family$ll_fun <- ll_gauss
  family$deviance <- dev_gauss
  family$ppd <- ppd_gauss

  return(family)
}

extend_family_gamma <- function(family) {
  kl_gamma <- function(pref, data, psub) {
    stop("KL-divergence for gamma not implemented yet.")
    ## mean(data$weights*(
    ##   p_sub$dis*(log(pref$dis)-log(p_sub$dis)+log(psub$mu)-log(pref$mu)) +
    ##     digamma(pref$dis)*(pref$dis - p_sub$dis) - lgamma(pref$dis) +
    ##     lgamma(p_sub$dis) + pref$mu*p_sub$dis/p_sub$mu - pref$dis))
  }
  dis_gamma <- function(pref, psub, wobs = 1) {
    ## TODO (Gamma()), IMPLEMENT THIS
    stop("Projection of dispersion parameter not yet implemented for family",
         " Gamma.")
    ## mean(wobs*((pref$mu - p_sub$mu)/
    ##                      family$mu.eta(family$linkfun(p_sub$mu))^2))
  }
  predvar_gamma <- function(mu, dis, wsample = 1) {
    stop("Family Gamma not implemented yet.")
  }
  ll_gamma <- function(mu, dis, y, weights = 1) {
    y <- as.matrix(y)
    dis <- matrix(rep(dis, each = length(y)), ncol = NCOL(mu))
    weights * dgamma(y, dis, dis / matrix(mu), log = TRUE)
  }
  dev_gamma <- function(mu, dis, y, weights = 1) {
    stop("Loss function not implemented for Gamma-family yet.")
    ## dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    ## weights*dgamma(y, dis, dis/matrix(mu), log= TRUE)
  }
  ppd_gamma <- function(mu, dis, weights = 1) {
    rgamma(length(mu), dis, dis / mu)
  }

  family$kl <- kl_gamma
  family$dis_fun <- dis_gamma
  family$predvar <- predvar_gamma
  family$ll_fun <- ll_gamma
  family$deviance <- dev_gamma
  family$ppd <- ppd_gamma

  return(family)
}

extend_family_student_t <- function(family) {
  kl_student_t <- function(pref, data, psub) {
    log(psub$dis)
  } #- 0.5*log(pref$var) # FIX THIS, NOT CORRECT
  dis_student_t <- function(pref, psub, wobs = 1) {
    s2 <- colSums(psub$w / sum(wobs) *
                    (pref$var + (pref$mu - psub$mu)^2)) # CHECK THIS
    sqrt(s2)
    ## stop('Projection of dispersion not yet implemented for student-t')
  }
  predvar_student_t <- function(mu, dis, wsample = 1) {
    wsample <- wsample / sum(wsample)
    mu_mean <- mu %*% wsample
    mu_var <- mu^2 %*% wsample - mu_mean^2
    as.vector(family$nu / (family$nu - 2) * sum(wsample * dis^2) + mu_var)
  }
  ll_student_t <- function(mu, dis, y, weights = 1) {
    y <- as.matrix(y)
    dis <- matrix(rep(dis, each = length(y)), ncol = NCOL(mu))
    weights * (dt((y - mu) / dis, family$nu, log = TRUE) - log(dis))
  }
  dev_student_t <- function(mu, y, weights = 1, dis = NULL) {
    if (is.null(dis)) {
      dis <- 1
    } else {
      dis <- matrix(rep(dis, each = length(y)), ncol = NCOL(mu))
    }
    if (NCOL(y) < NCOL(mu)) {
      y <- matrix(y, nrow = length(y), ncol = NCOL(mu))
    }
    (-2 * weights * (-0.5 * (family$nu + 1)
                     * log(1 + 1 / family$nu * ((y - mu) / dis)^2) - log(dis)))
  }
  ppd_student_t <- function(mu, dis, weights = 1) {
    rt(length(mu), family$nu) * dis + mu
  }

  family$kl <- kl_student_t
  family$dis_fun <- dis_student_t
  family$predvar <- predvar_student_t
  family$ll_fun <- ll_student_t
  family$deviance <- dev_student_t
  family$ppd <- ppd_student_t

  return(family)
}

.has_dispersion <- function(family) {
  # a function for checking whether the family has a dispersion parameter
  family$family %in% c("gaussian", "Student_t", "Gamma")
}

# A function for checking whether a `family` object has the required extra
# functions, that is, whether it has already been extended (typically by a call
# to extend_family()):
.has_family_extras <- function(family) {
  return(isTRUE(family$is_extended))
}
