library(glmproj)
set.seed(1234)
n <- 80
d <- 10
s <- 40
family <- gaussian()
funs <- kl_helpers(family)
x <- MASS::mvrnorm(n, rep(0, d), diag(rep(1, d)))
b_p <- t(MASS::mvrnorm(s, sample(d), diag(rep(0.1, d))))
mu_p <- x%*%b_p
dis_p <- sqrt(colMeans((rnorm(n, mu_p) - mu_p)^2))

exp_chosen <- c(10, 6, 8, 7, 4, 9, 2, 5, 1)

context("Forward selection")
test_that("Forward selection returns the correct sequence
          when avg = F and using a single core", {
  avg <- F
  cores <- 1
  chosen <- fsel(mu_p, x, b_p, w, dis_p, funs, avg, d - 1, cores)

  expect_equal(chosen, exp_chosen)
})

test_that("Forward selection returns the correct sequence
          when avg = F and using multiple cores", {
  avg <- F
  cores <- min(2, parallel::detectCores())
  chosen <- fsel(mu_p, x, b_p, w, dis_p, funs, avg, d - 1, cores)

  expect_equal(chosen, exp_chosen)
})

test_that("Forward selection returns the correct sequence
          when avg = T", {
  avg <- T
  cores <- 1
  chosen <- fsel(mu_p, x, b_p, w, dis_p, funs, avg, d - 1, cores)

  expect_equal(chosen, exp_chosen)
})



