# Family-specific helper functions
#
# `extend_family(family)` returns a `family` object augmented with auxiliary
# functions that are needed for computing KL-divergence, log predictive density,
# dispersion projection, etc.
#
# Missing: Quasi-families are not implemented. If dis_gamma is the correct shape
# parameter for projected Gamma regression, everything should be OK for gamma.

#' Extend a family
#'
#' This function adds some internally required elements to an object of class
#' `family` (see, e.g., [family()]). It is called internally by
#' [init_refmodel()], so you will rarely need to call it yourself.
#'
#' @param family An object of class `family`.
#'
#' @return The `family` object extended in the way needed by \pkg{projpred}.
#'
#' @export
extend_family <- function(family) {
  if (.has_family_extras(family)) {
    # If the family was already extended using this function, then return as-is:
    return(family)
  }
  extend_family_specific <- paste0("extend_family_", tolower(family$family))
  if (!exists(extend_family_specific, mode = "function")) {
    stop("Family '", family$family, "' is not supported by projpred.")
  }
  extend_family_specific <- get(extend_family_specific, mode = "function")
  family <- extend_family_specific(family)
  family$is_extended <- TRUE
  return(family)
}

extend_family_binomial <- function(family) {
  # Helper function for calculating the log PMF of the binomial distribution,
  # but (i) modified to be non-zero at `x` not contained in the support and (ii)
  # "reduced" in the sense of lacking the (additive) part
  # `log(choose(size, size * x))`:
  dbinom_log_reduced <- function(x, size, prob) {
    size * (x * log(prob) + (1 - x) * log(1 - prob))
  }

  ce_binom <- function(pref, data, psub) {
    ce_sums <- -colSums(
      dbinom_log_reduced(x = pref$mu, size = data$weights, prob = psub$mu)
    )
    return(ce_sums / sum(data$weights))
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
    -2 * dbinom_log_reduced(x = y, size = weights, prob = mu)
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
  family$ce <- ce_binom
  family$dis_fun <- dis_na
  family$predvar <- predvar_na
  family$ll_fun <- ll_binom
  family$deviance <- dev_binom
  family$ppd <- ppd_binom

  return(family)
}

extend_family_poisson <- function(family) {
  # Helper function for calculating the log PMF of the Poisson distribution,
  # but (i) modified to be non-zero at `x` not contained in the support and (ii)
  # "reduced" in the sense of lacking the (additive) part
  # `- wobs * log(factorial(x))`:
  dpois_log_reduced <- function(x, lamb, wobs) {
    wobs * (x * log(lamb) - lamb)
  }

  ce_poiss <- function(pref, data, psub) {
    ce_sums <- -colSums(
      dpois_log_reduced(x = pref$mu, lamb = psub$mu, wobs = data$weights)
    )
    return(ce_sums / sum(data$weights))
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
    -2 * dpois_log_reduced(x = y, lamb = mu, wobs = weights)
  }
  ppd_poiss <- function(mu, dis, weights = 1) {
    rpois(length(mu), mu)
  }

  family$ce <- ce_poiss
  family$dis_fun <- dis_na
  family$predvar <- predvar_na
  family$ll_fun <- ll_poiss
  family$deviance <- dev_poiss
  family$ppd <- ppd_poiss

  return(family)
}

extend_family_gaussian <- function(family) {
  # ce_gauss() does not give the actual cross-entropy (not even the one which
  # would result from dropping terms which would cancel out when calculating the
  # KL divergence) but a reasonable surrogate. This additional approximation was
  # already made back when this used to be the KL divergence, not the
  # cross-entropy.
  ce_gauss <- function(pref, data, psub) {
    ce_sums <- colSums(data$weights * (-2 * pref$mu * psub$mu + psub$mu^2))
    return(ce_sums / sum(data$weights))
  }
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

  family$ce <- ce_gauss
  family$dis_fun <- dis_gauss
  family$predvar <- predvar_gauss
  family$ll_fun <- ll_gauss
  family$deviance <- dev_gauss
  family$ppd <- ppd_gauss

  return(family)
}

extend_family_gamma <- function(family) {
  ce_gamma <- function(pref, data, psub) {
    stop("Cross-entropy for the Gamma() family not implemented yet.")
    ### TODO (Gamma()): This commented code stems from a time when this was
    ### still the actual KL divergence and not the (possibly reduced)
    ### cross-entropy ("possibly reduced" means: possibly reduced to only those
    ### terms which would not cancel out when calculating the KL divergence):
    ## mean(data$weights*(
    ##   p_sub$dis*(log(pref$dis)-log(p_sub$dis)+log(psub$mu)-log(pref$mu)) +
    ##     digamma(pref$dis)*(pref$dis - p_sub$dis) - lgamma(pref$dis) +
    ##     lgamma(p_sub$dis) + pref$mu*p_sub$dis/p_sub$mu - pref$dis))
    ###
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

  family$ce <- ce_gamma
  family$dis_fun <- dis_gamma
  family$predvar <- predvar_gamma
  family$ll_fun <- ll_gamma
  family$deviance <- dev_gamma
  family$ppd <- ppd_gamma

  return(family)
}

extend_family_student_t <- function(family) {
  ce_student_t <- function(pref, data, psub) {
    stop("Cross-entropy for the Student_t() family not implemented yet.")
    ### TODO (Student_t()): This commented code stems from a time when this was
    ### still the actual KL divergence and not the (possibly reduced)
    ### cross-entropy ("possibly reduced" means: possibly reduced to only those
    ### terms which would not cancel out when calculating the KL divergence):
    # log(psub$dis)
    ##- 0.5*log(pref$var) # FIX THIS, NOT CORRECT
    ###
  }
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

  family$ce <- ce_student_t
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
