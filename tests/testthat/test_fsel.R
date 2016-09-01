library(glmproj, quietly = T)

set.seed(1234)
n <- 80
nv <- 10
ns <- 10
family_kl <- kl_helpers(gaussian())
d <- list(x = MASS::mvrnorm(n, c(1,-1)*rep(1, nv), diag(rep(1, nv))),
          w = rep(1, n), y = rnorm(n),
          offset = rep(0, n))
p <- list(b = t(MASS::mvrnorm(ns, (nv:1)^2/15, diag(rep(1,nv)))),
          cluster_w = rep(1, ns))
p$mu <- family_kl$linkinv(d$x%*%p$b)
p$dis <- sqrt(colMeans((rnorm(n, p$mu)-p$mu)^2))
chosen <- 1:(nv-1)
b0 <- matrix(rowMeans(p$b))
args <- list(family_kl = family_kl, clust = F, nv = nv-1, intercept = F,
             verbose = F, ns = ns, rank_x = Matrix::rankMatrix(d$x))

context("Forward selection")
test_that("Forward selection returns a sensible sequence
          with a gaussian likelihood.", {
  sel <- fsel(p, d, d, p_clust = NULL, b0, args)
  expect_equal(sel$chosen[1:4], chosen[1:4])
})

