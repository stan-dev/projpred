context("ridge")

# tests for ridge regression, currently untested combinations
# - gaussian with inverse-link
# - binomial with log or cloglog-link
# - poisson with sqrt or id-link
# - Gamma with inverse or id-link
# - everything except gaussian with id-link for ridge penalty

set.seed(1235)
n <- 40
nterms <- 10
nterms_fit <- nterms - 5
x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
b <- runif(nterms) - 0.5
dis <- runif(1, 1, 2)
x_tr <- x[, 1:nterms_fit]
weights <- sample(1:4, n, replace = TRUE)
offset <- rnorm(n, 0, 1)

tol <- 3e-03
# some link-functions seem to need higher thresh-argument for glm_ridge
# (gaussian-log, binomial-cauchit, Gamma-log)
extra_thresh <- 1e-10

test_that(paste("glmfun: gradients should give the same results as finite",
                "differences"), {
  fdiffu <- function(f, x, h = 1e-3, order = 1) {
    ## function for computing derivative of univariate function f at x using
    ## finite difference.
    if (order != 1 && order != 2) {
      stop("Order must be either 1 or 2.")
    }
    n <- length(x)
    df <- rep(0, n)
    for (i in 1:n) {
      if (order == 1) {
        df[i] <- (f(x[i] + h) - f(x[i])) / h
      } else if (order == 2) {
        df[i] <- (f(x[i] + h) - 2 * f(x[i]) + f(x[i] - h)) / h^2
      }
    }
    return(df)
  }

  fams <- list(
    extend_family(gaussian(link = "identity")),
    extend_family(gaussian(link = "log")),
    extend_family(binomial(link = "logit")),
    extend_family(binomial(link = "probit")),
    extend_family(binomial(link = "cauchit")),
    extend_family(poisson(link = "log")),
    extend_family(Student_t(nu = 3, link = "identity")),
    extend_family(Student_t(nu = 4, link = "log")),
    extend_family(Student_t(nu = 7, link = "inverse"))
  )

  n <- 10
  weights <- sample(1:4, n, replace = TRUE)
  offset <- rnorm(n, 0, 1)

  for (i in seq_along(fams)) {
    fam <- fams[[i]]
    if (fam$family == "gaussian" || fam$family == "Student_t") {
      y <- rnorm(n)
    } else if (fam$family == "binomial") {
      y <- rbinom(n, 1, 0.6)
    } else if (fam$family == "poisson") {
      y <- rpois(n, 1)
    }

    devfun <- function(f)
      projpred:::pseudo_data(f, y, fam, weights = weights, offset = offset)$loss
    zfun <- function(f)
      projpred:::pseudo_data(f, y, fam, weights = weights, offset = offset)$z
    wfun <- function(f)
      projpred:::pseudo_data(f, y, fam, weights = weights, offset = offset)$w
    gradan <- function(f)
      sum(projpred:::pseudo_data(f, y, fam, weights = weights,
                                 offset = offset)$grad)
    gradfd <- function(f)
      fdiffu(devfun, f, h = 1e-5) # finite difference

    # compare analytic and finite difference gradients
    fval <- seq(-5, 5, len = 100)
    gan <- sapply(fval, gradan)
    gfd <- sapply(fval, gradfd)
    expect_equal(gan, gfd, tol = 1e-4 * max(abs(gan), abs(gfd)))
  }
})

test_that("glm_ridge: gaussian, id-link, intercept, lambda = 0", {
  fam <- extend_family(gaussian(link = "identity"))
  y <- rnorm(n, x %*% b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = TRUE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
    tolerance = tol
  )
})

test_that("glm_ridge: gaussian, id-link, no intercept, lambda = 0", {
  fam <- extend_family(gaussian(link = "identity"))
  y <- rnorm(n, x %*% b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = FALSE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: gaussian, id-link, intercept, lambda = 0.5", {
  fam <- extend_family(gaussian(link = "identity"))
  y <- rnorm(n, x %*% b, dis)
  lambda <- 0.5

  ridgefit <- glm_ridge(x_tr, y,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = TRUE
  )
  # analytic solution, no penalty on the intercept term
  penalty <- 0.5 * diag(c(0, rep(lambda, nterms_fit)))
  exp_beta <- c(solve(
    crossprod(cbind(1, x_tr) * sqrt(weights)) + penalty,
    crossprod(cbind(1, x_tr) * weights, y - offset)
  ))

  expect_equal(exp_beta, c(ridgefit$beta0, ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: gaussian, log-link, intercept, lambda = 0", {
  fam <- extend_family(gaussian(link = "log"))
  # intercept of 4 to ensure that y are positive
  y <- rnorm(n, fam$linkinv(x %*% b + 4), dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
    family = fam, lambda = lambda, weights = weights,
    offset = offset, intercept = TRUE,
    thresh = extra_thresh
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
    tolerance = tol
  )
})

test_that("glm_ridge: binomial, logit-link, intercept, lambda = 0", {
  fam <- extend_family(binomial(link = "logit"))
  y <- rbinom(n, weights, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights - y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y / weights,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = TRUE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
    tolerance = tol
  )
})

test_that("glm_ridge: binomial, logit-link, no intercept, lambda = 0", {
  fam <- extend_family(binomial(link = "logit"))
  y <- rbinom(n, weights, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights - y) ~ x_tr - 1, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y / weights,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = FALSE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: binomial, probit-link, intercept, lambda = 0", {
  fam <- extend_family(binomial(link = "probit"))
  y <- rbinom(n, weights, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights - y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y / weights,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = TRUE,
    thresh = extra_thresh
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
    tolerance = tol
  )
})

test_that("glm_ridge: binomial, cauchit-link, intercept, lambda = 0", {
  fam <- extend_family(binomial(link = "cauchit"))
  y <- rbinom(n, weights, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights - y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y / weights,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = TRUE,
    thresh = extra_thresh
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
    tolerance = tol
  )
})

test_that("glm_ridge: poisson, log-link, intercept, lambda = 0", {
  fam <- extend_family(poisson(link = "log"))
  y <- rpois(n, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = TRUE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
    tolerance = tol
  )
})

test_that("glm_ridge: poisson, log-link, no intercept, lambda = 0", {
  fam <- extend_family(poisson(link = "log"))
  y <- rpois(n, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
    family = fam, lambda = lambda,
    weights = weights, offset = offset, intercept = FALSE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

# test_that("glm_ridge: Gamma, log-link, intercept, lambda = 0", {
#   fam <- extend_family(Gamma(link = 'log'))
#   y <- rgamma(n, fam$linkinv(x%*%b + 1))
#   lambda <- 0
#
#   glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
#   ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda,
#                         weights = weights, offset = offset, intercept = TRUE,
#                         thresh = extra_thresh)
#
#   expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
#                tolerance = tol)
# })
