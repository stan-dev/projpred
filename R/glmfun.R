# The functions in this file are used to compute the elastic net coefficient
# paths for a GLM. The main function is glm_elnet, other functions are
# auxiliaries. The L1-regularized projection path is computed by replacing the
# actual data y by the fit of the reference model when calling glm_elnet. Uses
# functions in glmfun.cpp.


standardization <- function(x, center = TRUE, scale = TRUE, weights = NULL) {
  #
  # return the shift and scaling for each variable based on data matrix x.
  #
  w <- weights / sum(weights)
  if (center) {
    mx <- colSums(x * w)
  } else {
    mx <- rep(0, ncol(x))
  }
  if (scale) {
    sx <- apply(x, 2, weighted.sd, w)
  } else {
    sx <- rep(1, ncol(x))
  }
  return(list(shift = mx, scale = sx))
}


pseudo_data <- function(f, y, family, offset = rep(0, NROW(f)),
                        weights = rep(1.0, NROW(f)), obsvar = 0, wprev = NULL) {
  #
  # Returns locations z and weights w (inverse-variances) of the Gaussian
  # pseudo-observations based on the linear approximation to the link function
  # at f = eta = x*beta + beta0, as explained in McGullagh and Nelder (1989).
  # Returns also the loss (= negative log likelihood) and its pointwise
  # derivative w.r.t f at the current f.
  #
  mu <- family$linkinv(f + offset)
  dmu_df <- family$mu.eta(f + offset)
  z <- f + (y - mu) / dmu_df

  if (family$family == "Student_t") {
    # Student-t does not belong to the exponential family and thus it has its
    # own way of computing the observation weights
    if (is.null(wprev)) {
      # initialization of the em-iteration; loop recursively until stable
      # initial weights are found
      wprev <- weights
      while (TRUE) {
        wtemp <- pseudo_data(f, y, family, offset = offset, weights = weights,
                             wprev = wprev, obsvar = obsvar)$wobs
        if (max(abs(wtemp - wprev)) < 1e-6) {
          break
        }
        wprev <- wtemp
      }
    }
    # given the weights from the previous em-iteration, update s2 based on the
    # previous weights and mu, and then compute new weights w
    nu <- family$nu
    s2 <- sum(wprev * (obsvar + (y - mu)^2)) / sum(weights)
    wobs <- weights * (nu + 1) / (nu + 1 / s2 * (obsvar + (y - mu)^2))
    loss <- 0.5 * sum(family$deviance(mu, y, weights, sqrt(s2)))
    grad <- weights * (mu - y) / (nu * s2) * (nu + 1) /
      (1 + (y - mu)^2 / (nu * s2)) * dmu_df
  } else if (family$family %in% c("gaussian", "poisson", "binomial")) {
    # exponential family distributions
    wobs <- (weights * dmu_df^2) / family$variance(mu) # 2* because of deviance
    loss <- 0.5 * sum(family$deviance(mu, y, weights))
    grad <- -wobs * (z - f)
  } else {
    stop("Don't know how to compute quadratic approximation and gradients",
         sprintf(" for family '%s'.", family$family))
  }

  return(nlist(z, wobs, loss, grad))
}

lambda_grid <- function(x, y, family, offset, weights, intercept, penalty,
                        obsvar = 0, alpha = 1.0, lambda_min_ratio = 1e-2,
                        nlam = 100) {
  #
  # Standard lambda sequence as described in Friedman et al. (2009), section
  # 2.5. The grid will have nlam values, evenly spaced in the log-space between
  # lambda_max and lambda_min. lambda_max is the smallest value for which all
  # the regression coefficients will be zero (assuming alpha > 0, alpha = 0 will
  # be initialized as if alpha = 0.01). Returns also the initial solution
  # corresponding to the largest lambda (intercept and the unpenalized variables
  # will be nonzero).

  n <- dim(x)[1]

  if (alpha == 0) {
    # initialize ridge as if alpha = 0.01
    alpha <- 0.01
  }

  # find the initial solution, that is, values for the intercept (if included)
  # and those covariates that have penalty=0 (those which are always included,
  # if such exist)
  init <- glm_ridge(x[, penalty == 0, drop = FALSE], y,
                    family = family, lambda = 0, weights = weights,
                    offset = offset, obsvar = obsvar, intercept = intercept
  )
  f0 <- init$beta0 * rep(1, n)
  if (length(init$beta) > 0) {
    f0 <- f0 + as.vector(x[, penalty == 0, drop = FALSE] %*% init$beta)
  }

  obs <- pseudo_data(f0, y, family, offset, weights, obsvar = obsvar)
  resid <- obs$z - f0 # residual from the initial solution
  lambda_max_cand <- abs(t(x) %*% (resid * obs$wobs)) / (penalty * alpha)
  lambda_max <- max(lambda_max_cand[is.finite(lambda_max_cand)])
  ## to prevent some variable from entering at the first step due to numerical
  ## inaccuracy
  lambda_max <- 1.001 * lambda_max
  lambda_min <- lambda_min_ratio * lambda_max
  loglambda <- seq(log(lambda_min), log(lambda_max), len = nlam)

  beta <- rep(0, ncol(x))
  beta[penalty == 0] <- init$beta
  return(list(lambda = rev(exp(loglambda)), beta = beta,
              beta0 = init$beta0, w0 = obs$wobs))
}

glm_elnet <- function(x, y, family = gaussian(), nlambda = 100,
                      lambda_min_ratio = 1e-3, lambda = NULL, alpha = 1.0,
                      qa_updates_max = ifelse(family$family == "gaussian" &&
                                                family$link == "identity",
                                              1, 100),
                      pmax = dim(as.matrix(x))[2] + 1, pmax_strict = FALSE,
                      weights = NULL, offset = NULL, obsvar = 0,
                      intercept = TRUE, normalize = TRUE, penalty = NULL,
                      thresh = 1e-6) {
  #
  # Fits GLM with elastic net penalty on the regression coefficients.
  # Computes the whole regularization path.
  # Does not handle any dispersion parameters.
  #
  if (!.has_family_extras(family)) {
    family <- extend_family(family)
  }

  # ensure x is in matrix form and fill in missing weights and offsets
  x <- as.matrix(x)
  if (is.null(weights)) {
    weights <- rep(1.0, nrow(x))
  }
  if (is.null(offset)) {
    offset <- rep(0.0, nrow(x))
  }
  if (is.null(penalty)) {
    penalty <- rep(1.0, ncol(x))
  } else if (length(penalty) != ncol(x)) {
    stop(paste0("Incorrect length of penalty vector (should be ",
                ncol(x), ")."))
  }

  # standardize the features (notice that the variables are centered only if
  # intercept is used because otherwise the intercept would become nonzero
  # unintentionally)
  transf <- standardization(x, center = intercept, scale = normalize,
                            weights = weights)
  penalty[transf$scale == 0] <- Inf # ignore variables with zero variance
  transf$scale[transf$scale == 0] <- 1
  x <- t((t(x) - transf$shift) / transf$scale)

  # default lambda-sequence, including optimal start point
  if (is.null(lambda)) {
    temp <- lambda_grid(x, y, family, offset, weights, intercept, penalty,
                        alpha = alpha, obsvar = obsvar, nlam = nlambda,
                        lambda_min_ratio = lambda_min_ratio
    )
    lambda <- temp$lambda
    w0 <- temp$w0
    beta <- temp$beta
    beta0 <- temp$beta0
  } else {
    beta <- rep(0, ncol(x))
    beta0 <- 0
    w0 <- weights
  }

  # call the C++-function that serves as the workhorse
  pseudo_obs <- function(f, wprev) {
    pseudo_data(f, y, family, offset = offset, weights = weights,
                obsvar = obsvar, wprev = wprev)
  }
  out <- glm_elnet_c(
    x, pseudo_obs, lambda, alpha, intercept, penalty,
    thresh, qa_updates_max, pmax, pmax_strict, beta, beta0, w0
  )
  beta <- out[[1]]
  beta0 <- as.vector(out[[2]])

  # return the intercept and the coefficients on the original scale
  beta <- beta / transf$scale
  beta0 <- beta0 - colSums(transf$shift * beta)

  return(nlist(
    beta, beta0, w = out[[3]], lambda = lambda[seq_len(ncol(beta))],
    npasses = out[[4]], updates_qa = as.vector(out[[5]]),
    updates_as = as.vector(out[[6]])
  ))
}

glm_ridge <- function(x, y, family = gaussian(), lambda = 0, thresh = 1e-7,
                      qa_updates_max = NULL, weights = NULL, offset = NULL,
                      obsvar = 0, intercept = TRUE, penalty = NULL,
                      normalize = TRUE, la_approx = FALSE, beta_init = NULL,
                      beta0_init = NULL, ls_iter_max = 30) {
  #
  # Fits GLM with ridge penalty on the regression coefficients.
  # Does not handle any dispersion parameters.
  #
  if (is.null(x)) {
    x <- matrix(ncol = 0, nrow = length(y))
  }
  if (!.has_family_extras(family)) {
    family <- extend_family(family)
  }
  if (family$family == "gaussian" && family$link == "identity") {
    qa_updates_max <- 1
    ls_iter_max <- 1
  } else if (is.null(qa_updates_max)) {
    qa_updates_max <- 100
  }

  if (is.null(weights)) {
    weights <- rep(1.0, length(y))
  }
  if (is.null(offset)) {
    offset <- rep(0.0, length(y))
  }
  if (is.null(beta0_init)) {
    beta0_init <- 0
  }
  if (is.null(beta_init)) {
    beta_init <- rep(0, NCOL(x))
  }
  if (intercept) {
    beta_start <- c(beta0_init, beta_init)
  } else {
    beta_start <- beta_init
  }
  if (is.null(penalty)) {
    penalty <- rep(1.0, NCOL(x))
  }


  if (length(x) == 0) {
    if (intercept) {
      # model with intercept only (fit like model with no intercept but with one
      # constant predictor)
      x <- matrix(rep(1, length(y)), ncol = 1)
      w0 <- weights
      pseudo_obs <- function(f, wprev)
        pseudo_data(f, y, family, offset = offset, weights = weights,
                    obsvar = obsvar, wprev = wprev)
      out <- glm_ridge_c(x, pseudo_obs, lambda, FALSE, 1, beta_start, w0,
                         thresh, qa_updates_max, ls_iter_max)
      return(list(beta = matrix(integer(length = 0)),
                  beta0 = as.vector(out[[1]]), w = out[[3]], loss = out[[4]],
                  qa_updates = out[[5]]))
    } else {
      # null model with no predictors and no intercept
      pseudo_obs <- function(f, wprev)
        pseudo_data(f, y, family, offset = offset, weights = weights,
                    obsvar = obsvar, wprev = wprev)
      pobs <- pseudo_obs(rep(0, length(y)), weights)
      return(list(beta = matrix(integer(length = 0)), beta0 = 0, w = pobs$wobs,
                  qa_updates = 0))
    }
  }

  # normal case, at least one predictor
  x <- as.matrix(x) # ensure x is a matrix

  # standardize the features (notice that the variables are centered only if
  # intercept is used because otherwise the intercept would become nonzero
  # unintentionally)
  transf <- standardization(x, center = intercept, scale = normalize,
                            weights = weights)
  penalty[transf$scale == 0] <- Inf # ignore variables with zero variance
  transf$scale[transf$scale == 0] <- 1
  x <- t((t(x) - transf$shift) / transf$scale)

  # compute the solution
  w0 <- weights
  pseudo_obs <- function(f, wprev)
    pseudo_data(f, y, family, offset = offset, weights = weights,
                obsvar = obsvar, wprev = wprev)
  out <- glm_ridge_c(x, pseudo_obs, lambda, intercept, penalty, beta_start, w0,
                     thresh, qa_updates_max, ls_iter_max)
  beta <- out[[1]]
  beta0 <- as.vector(out[[2]])
  w <- out[[3]]
  loss <- out[[4]]

  # return the intercept and the coefficients on the original scale
  beta_orig <- beta / transf$scale
  beta0_orig <- beta0 - sum(transf$shift * beta_orig)

  out <- nlist(beta = beta_orig, beta0 = beta0_orig, w,
               qa_updates = out[[5]])
  return(out)
}

glm_forward <- function(x, y, family = gaussian(), lambda = 0, thresh = 1e-7,
                        qa_updates_max = NULL, weights = NULL, offset = NULL,
                        obsvar = 0, intercept = TRUE, penalty = NULL,
                        normalize = TRUE, pmax = dim(as.matrix(x))[2]) {
  #
  # Runs forward stepwise regression. Does not handle any dispersion parameters.
  #
  if (is.null(x)) {
    x <- matrix(ncol = 0, nrow = length(y))
  }
  if (!.has_family_extras(family)) {
    family <- extend_family(family)
  }
  if (family$family == "gaussian" && family$link == "identity") {
    qa_updates_max <- 1
  } else if (is.null(qa_updates_max)) {
    qa_updates_max <- 100
  }
  if (is.null(penalty)) {
    penalty <- rep(1.0, ncol(x))
  }


  # compute the null model
  out <- glm_ridge(NULL, y, family = family, lambda = lambda, thresh = thresh,
                   qa_updates_max = qa_updates_max, weights = weights,
                   offset = offset, obsvar = obsvar, intercept = intercept,
                   penalty = penalty
  )
  nullmodel <- list(beta = out$beta, beta0 = out$beta0,
                    varorder = integer(length = 0), w = out$w)

  if (length(x) == 0) {
    # return only the null model
    nullmodel$varorder <- integer(length = 0)
    return(nullmodel)
  }

  # normal case, at least one predictor
  x <- as.matrix(x)
  if (is.null(weights)) {
    weights <- rep(1.0, nrow(x))
  }
  if (is.null(offset)) {
    offset <- rep(0.0, nrow(x))
  }

  # standardize the features (notice that the variables are centered only if
  # intercept is used because otherwise the intercept would become nonzero
  # unintentionally)
  transf <- standardization(x, center = intercept, scale = normalize,
                            weights = weights)
  penalty[transf$scale == 0] <- Inf # ignore variables with zero variance
  transf$scale[transf$scale == 0] <- 1
  x <- t((t(x) - transf$shift) / transf$scale)

  # forward search (use the C++ function)
  w0 <- weights
  pseudo_obs <- function(f, wprev)
    pseudo_data(f, y, family, offset = offset, weights = weights,
                obsvar = obsvar, wprev = wprev)
  path <- glm_forward_c(x, pseudo_obs, lambda, intercept, penalty,
                        thresh, qa_updates_max, pmax, w0)
  beta <- cbind(rep(0, ncol(x)), path[[1]])
  beta0 <- c(nullmodel$beta0, as.vector(path[[2]]))

  # return the intercept and the coefficients on the original scale
  beta <- beta / transf$scale
  beta0 <- beta0 - colSums(transf$shift * beta)

  return(nlist(beta, beta0, varorder = as.vector(path[[3]]) + 1,
               w = cbind(nullmodel$w, path[[4]])))
}
