library(glmproj)

set.seed(1234)
n <- 100
s <- 800
d <- 10
x <- cbind(MASS::mvrnorm(n, rep(0, d), diag(rep(1, d))))
trb <- sample(1:d)
y <- rnorm(n, x%*%trb)
y_bin <- rbinom(n, 1, binomial()$linkinv(x%*%trb))

# suppress messages from stan
if (.Platform$OS.type == 'windows') {
  sink("NUL")
} else {
  sink("/dev/null")
}
fit <- stan_glm(y ~ x, chains = 1, iter = s, seed = 1234)
fit_bin <- stan_glm(y_bin ~ x, chains = 1, iter = s, seed = 1234, family = binomial())
sink()

context("glm_proj")
test_that("glm_proj returns non-increasing kld for gaussian models", {
  kl <- glm_proj(fit, avg = T, n_out = 200)$kl
  expect_true(all(kl[-1] <= kl[-length(kl)]))
})

test_that("glm_proj returns non-increasing kld for gaussian models", {
  kl_bin <- glm_proj(fit_bin, avg = T, n_out = 200)$kl
  expect_true(all(kl_bin[-1] <= kl_bin[-length(kl_bin)]))
})

test_that("changing d affects the output as expected", {
  exp_d <- 5
  d1 <- length(glm_proj(fit, d = exp_d, avg = T, n_out = 200)$kl)
  expect_equal(exp_d, d1)
})

test_that("changing n_out affects the output as expected", {
  exp_n_out <- 100
  n_out <- length(glm_proj(fit, n_out = exp_n_out, d = 2, avg = T)$b[[1]])
  expect_equal(exp_n_out, n_out)
})

