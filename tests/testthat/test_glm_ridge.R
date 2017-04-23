# tests for ridge regression, currently untested combinations
# - gaussian with inverse-link
# - binomial with log or cloglog-link
# - poisson with sqrt or id-link
# - Gamma with inverse or id-link
# - everything except gaussian with id-link for ridge penalty

set.seed(1235)
n <- 40
nv <- 10
nv_fit <- nv - 5
x <- matrix(rnorm(n*nv, 0, 1), n, nv)
b <- runif(nv)-0.5
dis <- runif(1, 1, 2)
x_tr <- x[,1:nv_fit]
weights <- sample(1:4, n, replace = T)
offset <- rnorm(n, 0, 1)

tol <- 1e-04
# some link-functions seem to need higher thresh-argument for glm_ridge
# (gaussian-log, binomial-cauchit, Gamma-log)
extra_thresh <- 1e-10


context("ridge")
test_that("glm_ridge: gaussian, id-link, intercept, lambda = 0", {
  fam <- gaussian(link = 'identity')
  y <- rnorm(n, x%*%b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol)
})

test_that("glm_ridge: gaussian, id-link, no intercept, lambda = 0", {
  fam <- gaussian(link = 'identity')
  y <- rnorm(n, x%*%b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = FALSE)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: gaussian, id-link, intercept, lambda = 0.5", {
  fam <- gaussian(link = 'identity')
  y <- rnorm(n, x%*%b, dis)
  lambda <- 0.5

  ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE)
  # analytic solution, penalty on the intercept term?
  penalty <- diag(c(lambda, rep(lambda, nv_fit)))
  exp_beta <- c(solve(crossprod(cbind(1, x_tr) * sqrt(weights)) + penalty,
                      crossprod(cbind(1, x_tr) * weights, y - offset)))

  expect_equal(exp_beta, c(ridgefit$beta0, ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: gaussian, log-link, intercept, lambda = 0", {
  fam <- gaussian(link = 'log')
  # intercept of 4 to ensure that y are positive
  y <- rnorm(n, fam$linkinv(x%*%b+4), dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda, weights = weights,
                        offset = offset, intercept = TRUE, qa_updates_max = 300,
                        thresh = extra_thresh)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol)
})

test_that("glm_ridge: binomial, logit-link, intercept, lambda = 0", {
  fam <- binomial(link = 'logit')
  y <- rbinom(n, weights, fam$linkinv(x%*%b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights-y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y/weights, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol)
})

test_that("glm_ridge: binomial, logit-link, no intercept, lambda = 0", {
  fam <- binomial(link = 'logit')
  y <- rbinom(n, weights, fam$linkinv(x%*%b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights-y) ~ x_tr - 1, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y/weights, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = FALSE)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: binomial, probit-link, intercept, lambda = 0", {
  fam <- binomial(link = 'probit')
  y <- rbinom(n, weights, fam$linkinv(x%*%b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights-y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y/weights, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE,
                        thresh = extra_thresh)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol)
})

test_that("glm_ridge: binomial, cauchit-link, intercept, lambda = 0", {
  fam <- binomial(link = 'cauchit')
  y <- rbinom(n, weights, fam$linkinv(x%*%b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights-y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y/weights, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE,
                        thresh = extra_thresh)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol)
})

test_that("glm_ridge: poisson, log-link, intercept, lambda = 0", {
  fam <- poisson(link = 'log')
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol)
})

test_that("glm_ridge: poisson, log-link, no intercept, lambda = 0", {
  fam <- poisson(link = 'log')
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0

  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = FALSE)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: Gamma, log-link, intercept, lambda = 0", {
  fam <- Gamma(link = 'log')
  y <- rgamma(n, fam$linkinv(x%*%b + 1))
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE,
                        thresh = extra_thresh)

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol)
})

