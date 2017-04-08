# tests for glm_elnet, currently untested combinations:
# - gaussian with non-id-link
# - binomial with non-logit-link
# - binomial with weights
# - poisson with non-log-link
# - ridge penalty
# - gamma with any link

set.seed(1235)
n <- 40
nv <- 10
nv_fit <- nv - 5
x <- mvrnorm(n, rep(0, nv), diag(rep(1, nv)))
b <- runif(nv)-0.5
dis <- runif(1, 1, 2)
x_tr <- x[,1:nv_fit]
weights <- sample(1:4, n, replace = T)
weights_norm <- weights / sum(weights) * n
# change this to something else once offsets work
offset <- rep(0, n)

tol <- 1e-04
extra_thresh <- 1e-10

context("lasso")
test_that("glm_elnet: gaussian, id-link, intercept, lambda = 0", {
  fam <- gaussian(link = 'identity')
  y <- rnorm(n, x%*%b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights,
                offset = offset)
  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, intercept = TRUE,
                        normalize = FALSE, weights = weights, offset = offset,
                        thresh = extra_thresh)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta0, lassofit$beta),
               tolerance = tol)
})

test_that("glm_elnet: gaussian, id-link, no intercept, lambda = 0", {
  fam <- gaussian(link = 'identity')
  y <- rnorm(n, x%*%b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)
  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, intercept = FALSE,
                        normalize = FALSE, weights = weights, offset = offset,
                        thresh = extra_thresh)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: gaussian, id-link, intercept, lambda = 0.5", {
  fam <- gaussian(link = 'identity')
  y <- rnorm(n, x%*%b, dis)
  lambda <- 0.5

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        weights = weights_norm, offset = offset,
                        intercept = TRUE, normalize = FALSE)
  # Calculated with:
  # glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
  #                             weights = weights_norm, offset = offset,
  #                             lambda = lambda/n,
  #                             intercept = TRUE, standardize = FALSE)
  # exp_beta <- unname(c(glmnetfit$a0, as.matrix(glmnetfit$beta)))
  exp_beta <- c(
    -0.3309525, -0.4230152, 0.2262128, -0.2064023, -0.2174287, -0.2842894
  )

  expect_equal(exp_beta, c(lassofit$beta0, lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: gaussian, id-link, no intercept, lambda = 0.5", {
  fam <- gaussian(link = 'identity')
  y <- rnorm(n, x%*%b, dis)
  lambda <- 0.5

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        weights = weights_norm, offset = offset,
                        intercept = FALSE, normalize = FALSE)
  # Calculated with:
  # glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
  #                             weights = weights_norm, offset = offset,
  #                             lambda = lambda/n,
  #                             intercept = FALSE, standardize = FALSE)
  # exp_beta <- c(as.matrix(glmnetfit$beta))
  exp_beta <- c(
    -0.1515209, 0.1867172, 0.05876578, -0.4579532, -0.4606212
  )

  expect_equal(exp_beta, c(lassofit$beta), tolerance = tol)
})


test_that("glm_elnet: binomial, logit-link, intercept, lambda = 0", {
  fam <- binomial(link = 'logit')
  y <- rbinom(n, weights, fam$linkinv(x%*%b))
  lambda <- 0

  lassofit <- glm_elnet(x_tr, y/weights, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights,
                        intercept = TRUE, normalize = FALSE)
  glmfit <- glm(cbind(y, weights-y) ~ x_tr, family = fam, offset = offset)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta0, lassofit$beta),
               tolerance = tol)
})

test_that("glm_elnet: binomial, logit-link, no intercept, lambda = 0", {
  fam <- binomial(link = 'logit')
  y <- rbinom(n, weights, fam$linkinv(x%*%b))
  lambda <- 0

  lassofit <- glm_elnet(x_tr, y/weights, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights,
                        intercept = FALSE, normalize = FALSE)
  glmfit <- glm(cbind(y, weights-y) ~ x_tr - 1, family = fam, offset = offset)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta), tolerance = tol)
})


test_that("glm_elnet: binomial, logit-link, intercept, lambda = 0.5", {
  fam <- binomial(link = 'logit')
  y <- rbinom(n, 1, fam$linkinv(x%*%b))
  lambda <- 0.5

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, thresh = extra_thresh,
                        intercept = TRUE, normalize = FALSE)
  # Calculated with:
  # glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
  #                             offset = offset, lambda = lambda/n,
  #                             intercept = TRUE, standardize = FALSE)
  # exp_beta <- unname(c(glmnetfit$a0, as.matrix(glmnetfit$beta)))
  exp_beta <- c(
    -0.4117375, -0.2879587, 0.7805274, 0.1626663, -0.135424, 0
  )

  expect_equal(exp_beta, c(lassofit$beta0, lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: binomial, logit-link, no intercept, lambda = 0.5", {
  fam <- binomial(link = 'logit')
  y <- rbinom(n, 1, fam$linkinv(x%*%b))
  lambda <- 0.5

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset,
                        intercept = FALSE, normalize = FALSE)
  # Calculated with:
  # glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
  #                             offset = offset, lambda = lambda/n,
  #                             intercept = FALSE, standardize = FALSE)
  # exp_beta <- c(as.matrix(glmnetfit$beta))
  exp_beta <- c(
    -0.3031029, 0.209245, 0, 0, -0.3442062
  )

  expect_equal(exp_beta, c(lassofit$beta), tolerance = tol)
})


test_that("glm_elnet: poisson, log-link, intercept, lambda = 0", {
  fam <- poisson(link = 'log')
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights,
                        intercept = TRUE, normalize = FALSE)
  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta0, lassofit$beta),
               tolerance = tol)
})

test_that("glm_elnet: poisson, log-link, no intercept, lambda = 0", {
  fam <- poisson(link = 'log')
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights,
                        intercept = FALSE, normalize = FALSE)
  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: poisson, log-link, intercept, lambda = 0.5", {
  fam <- poisson(link = 'log')
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0.5

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights_norm,
                        intercept = TRUE, normalize = FALSE, thresh = extra_thresh)
  # Calculated with:
  # glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
  #                             offset = offset,  weights = weights_norm, lambda = lambda/n,
  #                             intercept = TRUE, standardize = FALSE, thresh = extra_thresh)
  # exp_beta <- unname(c(glmnetfit$a0, as.matrix(glmnetfit$beta)))
  exp_beta <- c(
    -0.3980256, -0.8064176, 0.3061166, 0.4252627, 0.3223456, -0.2962351
  )

  expect_equal(exp_beta, c(lassofit$beta0, lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: poisson, log-link, no intercept, lambda = 0.5", {
  fam <- poisson(link = 'log')
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0.5

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights_norm,
                        intercept = FALSE, normalize = FALSE)
  # Calculated with:
  # glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
  #                             offset = offset,  weights = weights_norm, lambda = lambda/n,
  #                             intercept = FALSE, standardize = FALSE)
  # exp_beta <- c(as.matrix(glmnetfit$beta))
  exp_beta <- c(
    -0.2162512, 0.09895655, 0.1547091, -0.2937255, -0.2399947
  )

  expect_equal(exp_beta, c(lassofit$beta), tolerance = tol)
})
