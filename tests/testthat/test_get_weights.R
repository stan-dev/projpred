library(glmproj, quietly = T)
suppressPackageStartupMessages(library(rstanarm, quietly = T))

set.seed(1234)
n <- 40
s <- 250
family <- binomial()
x <- rnorm(n)
exp_w <- sample(1:5, n, replace = T)
y <- rbinom(n, exp_w, family$linkinv(x))
# suppress messages from stan
if (.Platform$OS.type == 'windows') {
  sink("NUL")
} else {
  sink("/dev/null")
}
fit1 <- stan_glm(y/exp_w ~ x, weights = exp_w,
                 family = family, chains = 1, iter = s, seed = 1234)
fit2 <- stan_glm(cbind(y, exp_w - y) ~ x,
                 family = family, chains = 1, iter = s, seed = 1234)
sink()

context("Get weights")
test_that("Weights are correctly extracted from a binomial fit
          when output is specified as percentage.", {
  w1 <- get_weights(fit1)
  expect_equal(w1, exp_w)
})
test_that("Weights are correctly extracted from a binomial fit
          when output is specified as successes and failures.", {
  w2 <- get_weights(fit2)
  expect_equal(w2, exp_w)
})
