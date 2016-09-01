library(glmproj, quietly = T)

set.seed(1234)
n <- 40
nv <- 10
nv_fit <- nv - 5
x <- MASS::mvrnorm(n, rep(0, nv), diag(rep(1, nv)))
d <- list(x = x[,1:nv_fit], w = rep(1, n), offset = rep(0, n))
p <- list(b = runif(nv)-0.5, dis = 2)

b0 <- rep(0, nv_fit)
d_binom <- within(d, w <- sample(1:4, n, replace = T))
tol <- 1.5e-06

context("NR")
test_that("NR works for gaussian model", {
  family_kl <- kl_helpers(gaussian())
  p <- within(p, mu <- family_kl$linkinv(x%*%p$b))
  nrfit <- NR(p, d, b0, family_kl)
  glmfit <- glm(p$mu ~ d$x - 1, family = family_kl)

  exp_beta <- unname(coef(glmfit))
  exp_dis <- sqrt(summary(glmfit)$dispersion*(n-nv_fit)/n + p$dis^2)
  beta <- drop(nrfit$b)
  dis <- nrfit$dis

  expect_equal(exp_beta, beta, tolerance = tol)
  expect_equal(exp_dis, dis, tolerance = tol)
})

test_that("NR works for binomial model with logit link", {
  family_kl <- kl_helpers(binomial(link = 'logit'))
  p <- within(p, mu <- family_kl$linkinv(x%*%p$b))
  nrfit <- NR(p, d, b0, family_kl)
  glmfit <- suppressWarnings(glm(p$mu ~ d$x - 1, family= family_kl))

  exp_beta <- unname(coef(glmfit))
  beta <- drop(nrfit$b)

  expect_equal(exp_beta, beta, tolerance = tol)
})

test_that("NR works for binomial model with probit link", {
  family_kl <- kl_helpers(binomial(link = 'probit'))
  p <- within(p, mu <- family_kl$linkinv(x%*%p$b))
  nrfit <- NR(p, d, b0, family_kl)
  glmfit <- suppressWarnings(glm(p$mu ~ d$x - 1, family= family_kl))

  exp_beta <- unname(coef(glmfit))
  beta <- drop(nrfit$b)

  expect_equal(exp_beta, beta, tolerance = tol)
})

test_that("NR works for poisson model", {
  family_kl <- kl_helpers(poisson())
  p <- within(p, mu <- family_kl$linkinv(x%*%p$b))
  nrfit <- NR(p, d, b0, family_kl)
  glmfit <- suppressWarnings(glm(p$mu ~ d$x - 1, family= family_kl))

  exp_beta <- unname(coef(glmfit))
  beta <- drop(nrfit$b)

  expect_equal(exp_beta, beta, tolerance = tol)
})
