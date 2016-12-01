# set.seed(1234)
# n <- 40
# nv <- 10
# nv_fit <- nv - 5
# x <- MASS::mvrnorm(n, rep(0, nv), diag(rep(1, nv)))
# b = runif(nv)-0.5
# d_train <- list(x = x[,1:nv_fit], weights = sample(1:4, n, replace = T), offset = rep(0, n))
# p_full <- list(dis = 2)
#
# b0 <- rep(0, nv_fit)
# tol <- 1.5e-06
#
# context("IRLS")
# test_that("IRLS works for a gaussian model with inverse link", {
#   fam <- gaussian(link = 'log')
#   args <- .init_args(list(), d_train, fam)
#   p_full$mu <- args$family_kl$linkinv(x%*%b + d_train$offset)
#   irlsfit <- IRLS(p_full, d_train, b0, args)
#   glmfit <- glm(p_full$mu ~ d_train$x - 1, family = fam,
#                 weights = d_train$weights, offset = d_train$offset)
#
#   exp_beta <- unname(coef(glmfit))
#   exp_dis <- sqrt(summary(glmfit)$dispersion*(n-nv_fit)/n + p_full$dis^2)
#   beta <- drop(irlsfit$b)
#   dis <- irlsfit$dis
#
#   expect_equal(exp_beta, beta, tolerance = tol)
#   expect_equal(exp_dis, dis, tolerance = tol)
# })
#
# test_that("IRLS works for binomial model with logit link", {
#   fam <- binomial(link = 'logit')
#   args <- .init_args(list(), d_train, fam)
#   p_full$mu <- args$family_kl$linkinv(x%*%b + d_train$offset)
#   irlsfit <- IRLS(p_full, d_train, b0, args)
#   glmfit <- suppressWarnings(
#     glm(p_full$mu ~ d_train$x - 1, family = fam,
#         weights = d_train$weights, offset = d_train$offset))
#
#   exp_beta <- unname(coef(glmfit))
#   beta <- drop(irlsfit$b)
#
#   expect_equal(exp_beta, beta, tolerance = tol)
# })
#
# test_that("IRLS works for poisson model with log link", {
#   fam <- poisson(link = 'log')
#   args <- .init_args(list(), d_train, fam)
#   p_full$mu <- args$family_kl$linkinv(x%*%b + d_train$offset)
#   irlsfit <- IRLS(p_full, d_train, b0, args)
#   glmfit <- suppressWarnings(
#     glm(p_full$mu ~ d_train$x - 1, family = fam,
#         weights = d_train$weights, offset = d_train$offset))
#
#   exp_beta <- unname(coef(glmfit))
#   beta <- drop(irlsfit$b)
#
#   expect_equal(exp_beta, beta, tolerance = tol)
# })
