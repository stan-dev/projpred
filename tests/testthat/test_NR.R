library(glmproj, quietly = T)

set.seed(1234)
n <- 40
d_p <- 10
d_q <- d_p - 5
x <- MASS::mvrnorm(n, rep(0, d_p), diag(rep(1, d_p)))
x_q <- x[,1:d_q]
b_p <- runif(d_p)-0.5
dis_p <- 2
b0 <- rep(0,d_q)
w <- rep(1,n)
w_binom <- sample(1:4, n, replace = T)

context("NR")
test_that("NR works for gaussian model", {
  family <- gaussian()
  mu_p <- family$linkinv(x%*%b_p)
  funs <- kl_helpers(family)
  nrfit <- NR(mu_p, x_q, b0, w, dis_p, funs)
  glmfit <- glm(mu_p ~ x_q - 1, family = family)

  exp_beta <- unname(coef(glmfit))
  exp_dis <- sqrt(summary(glmfit)$dispersion*(n-d_q)/n + dis_p^2)
  beta <- drop(nrfit$b)
  dis <- nrfit$dis

  expect_equal(exp_beta, beta, tolerance = 0.001)
  expect_equal(exp_dis, dis, tolerance = 0.001)
})

test_that("NR works for binomial model with logit link", {
  family <- binomial(link = 'logit')
  mu_p <- family$linkinv(x%*%b_p)
  funs <- kl_helpers(family)
  nrfit <- NR(mu_p, x_q, b0, w_binom, NA, funs)
  glmfit <- suppressWarnings(glm(mu_p ~ x_q - 1, family = family, weights = w_binom))

  exp_beta <- unname(coef(glmfit))
  beta <- drop(nrfit$b)

  expect_equal(exp_beta, beta, tolerance = 0.001)
})

test_that("NR works for binomial model with logit link", {
  family <- binomial(link = 'probit')
  mu_p <- family$linkinv(x%*%b_p)
  funs <- kl_helpers(family)
  nrfit <- NR(mu_p, x_q, b0, w_binom, NA, funs)
  glmfit <- suppressWarnings(glm(mu_p~x_q-1, family = family, weights = w_binom))

  exp_beta <- unname(coef(glmfit))
  beta <- drop(nrfit$b)

  expect_equal(exp_beta, beta, tolerance = 0.001)
})

test_that("NR works for poisson model", {
  family <- poisson()
  mu_p <- family$linkinv(x%*%b_p)
  funs <- kl_helpers(family)
  nrfit <- NR(mu_p, x_q, b0, w, NA, kl_helpers(family))
  glmfit <- suppressWarnings(glm(mu_p ~ x_q - 1, family = family))

  exp_beta <- unname(coef(glmfit))
  beta <- drop(nrfit$b)

  expect_equal(exp_beta, beta, tolerance = 0.001)
})
