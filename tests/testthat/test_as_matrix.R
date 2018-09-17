context('as.matrix.projection')
library(rstanarm)

# tests for as_matrix


set.seed(1235)
n <- 40
nv <- 5
x <- matrix(rnorm(n*nv, 0, 1), n, nv)
b <- runif(nv)-0.5
dis <- runif(1, 1, 2)
weights <- sample(1:4, n, replace = T)
offset <- rnorm(n)
chains <- 2
seed <- 1235
iter <- 500
source(file.path('helpers', 'SW.R'))


f_gauss <- gaussian()
df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x%*%b), dis), x = x)
f_binom <- binomial()
df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x%*%b)), x = x)

SW(
  fit_gauss <- stan_glm(y ~ x, family = f_gauss, data = df_gauss, QR = T,
                        weights = weights, offset = offset,
                        chains = chains, seed = seed, iter = iter)
)
SW(
  fit_binom <- stan_glm(cbind(y, weights-y) ~ x, family = f_binom, QR = T,
                        data = df_binom, weights = weights, offset = offset,
                        chains = chains, seed = seed, iter = iter)
)

vs_gauss <- varsel(fit_gauss)
vs_binom <- varsel(fit_binom)
vind <- c(2,3)
ns <- 100
p_gauss <- project(vs_gauss, vind = vind, ns = ns)
p_binom <- project(vs_binom, vind = vind, ns = ns)



test_that("as.matrix.projection returns the relevant variables for gaussian", {
  m <- as.matrix(p_gauss)
  expect_equal(colnames(m), c(names(coef(fit_gauss))[c(1, vind + 1)], 'sigma'))
  expect_equal(dim(m), c(ns, length(vind) + 2))
})

test_that("as.matrix.projection returns the relevant variables for binomial", {
  m <- as.matrix(p_binom)
  expect_equal(colnames(m), names(coef(fit_binom))[c(1, vind + 1)])
  expect_equal(dim(m), c(ns, length(vind) + 1))
})

test_that("as.matrix.projection works as expected without an intercept", {
  p_nointercept <- project(vs_gauss, vind = vind, ns = ns, intercept = FALSE)
  m <- as.matrix(p_nointercept)
  expect_equal(colnames(m), c(names(coef(fit_gauss))[vind + 1], 'sigma'))
  expect_equal(dim(m), c(ns, length(vind) + 1))
})

test_that("as.matrix.projection works as expected with zero variables", {
  p_novars <- project(vs_gauss, nv = 0, ns = ns, intercept = F)
  m <- as.matrix(p_novars)
  expect_equal(colnames(m), 'sigma')
  expect_equal(dim(m), c(ns, 1))
})



test_that("as.matrix.projection gives a warning but works with clustering", {
  nc <- 3
  p_clust <- project(vs_gauss, vind = vind, nc = nc)
  expect_warning(m <- as.matrix(p_clust))
  expect_equal(colnames(m), c(names(coef(fit_gauss))[c(1, vind + 1)], 'sigma'))
  expect_equal(dim(m), c(nc, length(vind) + 2))
})

