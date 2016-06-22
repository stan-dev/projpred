library(glmproj)

set.seed(1234)
n <- 80
d <- 5
s <- 4
family <- gaussian()
funs <- kl_helpers(family)
x <- MASS::mvrnorm(n, rep(0, d), diag(rep(1, d)))
b_p <- t(MASS::mvrnorm(s, d:1, diag(rep(1,d))))
mu_p <- x%*%b_p
dis_p <- sqrt(colMeans((rnorm(n,mu_p)-mu_p)^2))
chosen <- 1:(d-1)
helperf <- function(ind) coef(lm(mu_p[, ind] ~ x[, 1:(d-1)] - 1))
exp_beta <- unname(sapply(1:s, helperf))

context("Parameter projection")
test_that("Check that projection works for 4 samples with 1 core", {
  cores <- 1
  beta <- proj_params(mu_p, x, b_p, w, dis_p, funs, chosen, cores, F)$b[[d-1]]

  expect_equal(exp_beta, beta, tolerance = 0.001)
})

test_that("Check that projection works for 4 samples and multiple cores", {
  cores <- min(2, parallel::detectCores())
  beta <- proj_params(mu_p, x, b_p, w, dis_p, funs, chosen, cores, F)$b[[d-1]]

  expect_equal(exp_beta, beta, tolerance = 0.001)
})
