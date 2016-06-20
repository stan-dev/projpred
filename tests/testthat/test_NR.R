library(glmproj)

context("NR (gaussian)")
test_that("NR finds the correct weight vector and the dispersion parameter", {
  set.seed(1337)
  n <- 40
  d_p <- 10
  d_q <- d_p - 5
  family <- gaussian()
  x <- MASS::mvrnorm(n, rep(0, d_p), diag(rep(1, d_p)))
  x_q <- x[,1:d_q]
  b_p <- runif(d_p)-0.5
  mu_p <- family$linkinv(x%*%b_p)
  dis_p <- 1
  b0 <- rep(0,d_q)
  w <- rep(1,n)
  funs <- family_kls(family)
  nrfit <- NR(mu_p, x_q, b0, w, dis_p, funs)
  glmfit <- glm(mu_p~x_q-1, family = family)

  exp_beta <- unname(coef(glmfit))
  exp_dis <- sqrt(summary(glmfit)$dispersion*(n-d_q)/n + dis_p^2)

  expect_equal(exp_beta, drop(nrfit$b), tolerance = 0.001)
  expect_equal(exp_dis, nrfit$dis, tolerance = 0.001)
})

context("NR (binomial, logit)")
test_that("NR finds the correct weight vector", {
  set.seed(1337)
  n <- 40
  d_p <- 10
  d_q <- d_p - 5
  family <- binomial(link = 'logit')
  x <- MASS::mvrnorm(n, rep(0, d_p), diag(rep(1, d_p)))
  x_q <- x[,1:d_q]
  b_p <- runif(d_p)-0.5
  w <- sample(1:4, n, replace = T)
  mu_p <- family$linkinv(x%*%b_p)
  b0 <- rep(0 ,d_q)
  funs <- family_kls(family)
  nrfit <- NR(mu_p, x_q, b0, w, NA, funs)
  glmfit <- suppressWarnings(glm(mu_p~x_q-1, family = family, weights = w))

  exp_beta <- unname(coef(glmfit))

  expect_equal(exp_beta, drop(nrfit$b), tolerance = 0.001)

})

context("NR (binomial, probit)")
test_that("NR finds the correct weight vector", {
  set.seed(1337)
  n <- 40
  d_p <- 10
  d_q <- d_p - 5
  family <- binomial(link = 'probit')
  x <- MASS::mvrnorm(n, rep(0, d_p), diag(rep(1, d_p)))
  x_q <- x[,1:d_q]
  b_p <- runif(d_p)-0.5
  w <- sample(1:4, n, replace = T)
  mu_p <- family$linkinv(x%*%b_p)
  b0 <- rep(0 ,d_q)
  funs <- family_kls(family)
  nrfit <- NR(mu_p, x_q, b0, w, NA, funs)
  glmfit <- suppressWarnings(glm(mu_p~x_q-1, family = family, weights = w))

  exp_beta <- unname(coef(glmfit))

  expect_equal(exp_beta, drop(nrfit$b), tolerance = 0.001)

})

context("NR (poisson)")
test_that("NR finds the correct weight vector", {
  set.seed(1337)
  n <- 40
  d_p <- 10
  d_q <- d_p - 5
  family <- poisson()
  x <- MASS::mvrnorm(n, rep(0, d_p), diag(rep(1, d_p)))
  x_q <- x[, 1:d_q]
  b_p <- runif(d_p)-0.5
  w <- rep(1, n)
  mu_p <- family$linkinv(x%*%b_p)
  b0 <- rep(0 ,d_q)
  funs <- family_kls(family)
  nrfit <- NR(mu_p, x_q, b0, w, NA, family_kls(family))
  glmfit <- suppressWarnings(glm(mu_p~x_q-1, family = family))

  exp_beta <- unname(coef(glmfit))

  expect_equal(exp_beta, drop(nrfit$b), tolerance = 0.001)

})
