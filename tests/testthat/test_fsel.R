set.seed(1234)
n <- 80
nv <- 10
ns <- 10
d_train <- list(x = MASS::mvrnorm(n, c(1,-1)*rep(1, nv), diag(rep(1, nv))),
                y = rnorm(n), weights = rep(1, n), offset = rep(0, n))
args <- .init_args(list(verbose = F), d_train, kl_helpers(gaussian()))

b <- t(MASS::mvrnorm(ns, (nv:1)^2/15, diag(rep(1,nv))))
p_full <- list(mu = args$family_kl$linkinv(d_train$x%*%b),
               cluster_w = rep(1/ns, ns))
p_full$dis <- sqrt(colMeans((rnorm(n, p_full$mu)-p_full$mu)^2))

chosen <- 1:(nv-1)
b0 <- matrix(rowMeans(b))

context("Forward selection")
test_that("Forward selection returns a sensible sequence
          with a gaussian likelihood.", {
  sel <- fsel(p_full, d_train, d_train, p_clust = NULL, b0, args)
  expect_equal(sel$chosen[1:4], chosen[1:4])
})

