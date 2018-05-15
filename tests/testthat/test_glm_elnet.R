# tests for glm_elnet, currently untested combinations:
# - gaussian with non-id-link
# - binomial with non-logit-link
# - binomial with weights
# - poisson with non-log-link
# - ridge penalty
# - gamma with any link

if (!requireNamespace('glmnet', quietly = TRUE)) {
  stop('glmnet needed for this function to work. Please install it.',
       call. = FALSE)
}

set.seed(1235)
n <- 40
nv <- 10
nv_fit <- nv - 5
x <- matrix(rnorm(n*nv, 0, 1), n, nv)
b <- runif(nv)-0.5
dis <- runif(1, 1, 2)
x_tr <- x[,1:nv_fit]
weights <- sample(1:4, n, replace = T)
weights_norm <- weights / sum(weights) * n
offset <- rnorm(n)
penalty <- runif(ncol(x_tr))+0.5
penalty <- penalty/sum(penalty)*ncol(x_tr) # must scale the penalties to be comparable to glmnet

tol <- 1e-04
extra_thresh <- 1e-10

context("elnet")
test_that("glm_elnet: gaussian, id-link, intercept, lambda = 0", {
  fam <- kl_helpers(gaussian(link = 'identity'))
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
  fam <- kl_helpers(gaussian(link = 'identity'))
  y <- rnorm(n, x%*%b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)
  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, intercept = FALSE,
                        normalize = FALSE, weights = weights, offset = offset,
                        thresh = extra_thresh)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: gaussian, id-link, intercept, lambda = 7", {
	fam <- kl_helpers(gaussian(link = 'identity'))
  y <- rnorm(n, x%*%b, dis)
  lambda <- 7.0

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        weights = weights_norm, offset = offset,
                        intercept = TRUE, normalize = FALSE, penalty = penalty)
  glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
                              weights = weights_norm, offset = offset,
                              lambda = lambda/n,
                              intercept = TRUE, standardize = FALSE,
                              penalty.factor = penalty)
  exp_beta <- unname(c(glmnetfit$a0, as.matrix(glmnetfit$beta)))

  expect_equal(exp_beta, c(lassofit$beta0, lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: gaussian, id-link, no intercept, lambda = 5", {
	fam <- kl_helpers(gaussian(link = 'identity'))
  y <- rnorm(n, x%*%b, dis)
  lambda <- 5

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        weights = weights_norm, offset = offset,
                        intercept = FALSE, normalize = FALSE, penalty = penalty)
  glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
                              weights = weights_norm, offset = offset,
                              lambda = lambda/n,
                              intercept = FALSE, standardize = FALSE,
                              penalty.factor = penalty)
  exp_beta <- c(as.matrix(glmnetfit$beta))

  
  expect_equal(exp_beta, c(lassofit$beta), tolerance = tol)
})


test_that("glm_elnet: binomial, logit-link, intercept, lambda = 0", {
  fam <- kl_helpers(binomial(link = 'logit'))
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
	fam <- kl_helpers(binomial(link = 'logit'))
  y <- rbinom(n, weights, fam$linkinv(x%*%b))
  lambda <- 0

  lassofit <- glm_elnet(x_tr, y/weights, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights,
                        intercept = FALSE, normalize = FALSE)
  glmfit <- glm(cbind(y, weights-y) ~ x_tr - 1, family = fam, offset = offset)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta), tolerance = tol)
})


test_that("glm_elnet: binomial, logit-link, intercept, lambda = 3", {
	fam <- kl_helpers(binomial(link = 'logit'))
  y <- rbinom(n, 1, fam$linkinv(x%*%b))
  lambda <- 3

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, thresh = extra_thresh,
                        intercept = TRUE, normalize = FALSE, penalty = penalty)

  glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
                              offset = offset, lambda = lambda/n,
                              intercept = TRUE, standardize = FALSE,
                              penalty.factor = penalty)
  exp_beta <- unname(c(glmnetfit$a0, as.matrix(glmnetfit$beta)))
  
  expect_equal(exp_beta, c(lassofit$beta0, lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: binomial, logit-link, no intercept, lambda = 2", {
	fam <- kl_helpers(binomial(link = 'logit'))
  y <- rbinom(n, 1, fam$linkinv(x%*%b))
  lambda <- 2

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset,
                        intercept = FALSE, normalize = FALSE, penalty = penalty)
  glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
                              offset = offset, lambda = lambda/n,
                              intercept = FALSE, standardize = FALSE,
                              penalty.factor = penalty)
  exp_beta <- c(as.matrix(glmnetfit$beta))
  
  expect_equal(exp_beta, c(lassofit$beta), tolerance = tol)
})


test_that("glm_elnet: poisson, log-link, intercept, lambda = 0", {
  fam <- kl_helpers(poisson(link = 'log'))
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights, intercept = TRUE,
                        normalize = FALSE, thresh = extra_thresh)
  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta0, lassofit$beta),
               tolerance = tol)
})

test_that("glm_elnet: poisson, log-link, no intercept, lambda = 0", {
	fam <- kl_helpers(poisson(link = 'log'))
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights,
                        intercept = FALSE, normalize = FALSE)
  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)

  expect_equal(unname(coef(glmfit)), c(lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: poisson, log-link, intercept, lambda = 0.5", {
	fam <- kl_helpers(poisson(link = 'log'))
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 0.5

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights_norm,
                        intercept = TRUE, normalize = FALSE, thresh = extra_thresh,
                        penalty = penalty)

  glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
                              offset = offset,  weights = weights_norm, lambda = lambda/n,
                              intercept = TRUE, standardize = FALSE, thresh = extra_thresh,
                              penalty.factor = penalty)
  exp_beta <- unname(c(glmnetfit$a0, as.matrix(glmnetfit$beta)))

  expect_equal(exp_beta, c(lassofit$beta0, lassofit$beta), tolerance = tol)
})

test_that("glm_elnet: poisson, log-link, no intercept, lambda = 3", {
	fam <- kl_helpers(poisson(link = 'log'))
  y <- rpois(n, fam$linkinv(x%*%b))
  lambda <- 3

  lassofit <- glm_elnet(x_tr, y, family = fam, lambda = lambda, alpha = 1,
                        offset = offset, weights = weights_norm,
                        intercept = FALSE, normalize = FALSE, penalty = penalty)
  glmnetfit <- glmnet::glmnet(x_tr, y, family = fam$family, alpha = 1,
                              offset = offset,  weights = weights_norm, lambda = lambda/n,
                              intercept = FALSE, standardize = FALSE, penalty.factor = penalty)
  exp_beta <- c(as.matrix(glmnetfit$beta))
  
  expect_equal(exp_beta, c(lassofit$beta), tolerance = tol)
})


test_that("glm_elnet: poisson, log-link, normalization should not affect the maximum likelihood solution", {
  fam <- kl_helpers(poisson(link = 'log'))
  y <- rpois(n, fam$linkinv(x%*%b))
  
  nlam <- 100
  elnetfit1 <- glm_elnet(x_tr, y, family = fam, nlambda=nlam, lambda_min_ratio=1e-7,
                        offset = offset, weights = weights_norm,
                        intercept = TRUE, normalize = FALSE)
  elnetfit2 <- glm_elnet(x_tr, y, family = fam, nlambda=nlam, lambda_min_ratio=1e-7,
                         offset = offset, weights = weights_norm,
                         intercept = TRUE, normalize = TRUE)
  
  expect_equal(c(elnetfit1$beta0[nlam], elnetfit1$beta[,nlam]), 
               c(elnetfit2$beta0[nlam], elnetfit2$beta[,nlam]), tolerance = tol)
})



test_that("glm_elnet with alpha=0 and glm_ridge give the same result.", {
  
  for (famstr in c('gaussian', 'binomial', 'poisson')) {
    
    if (famstr == 'gaussian') {
      fam <- kl_helpers(gaussian(link = 'identity'))
      y <- rnorm(n, x%*%b, 0.5)
    } else if (famstr == 'binomial') {
      fam <- kl_helpers(binomial(link = 'logit'))
      y <- rbinom(n, weights, fam$linkinv(x%*%b))
    } else if (famstr == 'poisson') {
      fam <- kl_helpers(poisson(link = 'log'))
      y <- rpois(n, fam$linkinv(x%*%b))
    }
    
    for (lambda in c(0.1, 1, 10)) {
      for (intercept in c(T,F)) {
        for (normalize in c(T,F)) {
          elnetfit <- projpred:::glm_elnet(x_tr, y/weights, family = fam, lambda = lambda, alpha = 0,
                                           offset = offset, weights = weights, penalty = penalty,
                                           intercept = intercept, normalize = normalize, thresh = 1e-12)
          ridgefit <- projpred:::glm_ridge(x_tr, y/weights, family = fam, lambda = lambda, 
                                           offset = offset, weights = weights, penalty = penalty,
                                           intercept = intercept, normalize = normalize, thresh = 1e-12)
          
          expect_equal( c(ridgefit$beta0, ridgefit$beta), c(elnetfit$beta0, elnetfit$beta),
                        tolerance = 1e-5)
        }
      }
    }  
  }
})
