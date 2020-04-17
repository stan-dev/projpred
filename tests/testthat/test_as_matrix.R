context("as.matrix.projection")

# tests for as_matrix

if (require(rstanarm) && require(brms)) {
  set.seed(1235)
  n <- 40
  nv <- 5
  x <- matrix(rnorm(n * nv, 0, 1), n, nv)
  b <- runif(nv) - 0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = TRUE)
  offset <- rnorm(n)
  chains <- 2
  seed <- 1235
  iter <- 500
  source(file.path("helpers", "SW.R"))

  f_gauss <- gaussian()
  df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x %*% b), dis), x = x)
  f_binom <- binomial()
  df_binom <- data.frame(
    y = rbinom(n, weights, f_binom$linkinv(x %*% b)), x = x,
    weights = weights, offset = offset
  )

  SW(
    fit_gauss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_gauss, data = df_gauss,
      chains = chains, seed = seed, iter = iter
    )
  )
  SW(
    fit_binom <- stan_glm(cbind(y, weights - y) ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_binom,
      data = df_binom, chains = chains, seed = seed, iter = iter
    )
  )

  vs_gauss <- varsel(fit_gauss)
  vs_binom <- varsel(fit_binom)
  solution_terms <- c(2, 3)
  number_samples <- 100
  p_gauss <- project(vs_gauss, solution_terms = solution_terms,
                     number_samples = number_samples)
  p_binom <- project(vs_binom, solution_terms = solution_terms,
                     number_samples = number_samples)

  ## test_that("as.matrix.projection returns the relevant variables for gaussian", {
  ##   m <- as.matrix(p_gauss)
  ##   expect_length(setdiff(colnames(m),
  ##                 c('Intercept', vs_gauss$solution_terms[solution_terms],
  ##                   'sigma')), 0)
  ##   expect_equal(dim(m), c(number_samples, length(solution_terms) + 2))
  ## })

  ## test_that("as.matrix.projection returns the relevant variables for binomial", {
  ##   m <- as.matrix(p_binom)
  ##   expect_length(setdiff(colnames(m), c("Intercept", vs_binom$solution_terms[solution_terms])), 0)
  ##   expect_equal(dim(m), c(number_samples, length(solution_terms) + 1))
  ## })

  ## test_that("as.matrix.projection works as expected with zero variables", {
  ##   p_novars <- project(vs_gauss, nv = 0, number_samples = number_samples)
  ##   m <- as.matrix(p_novars)
  ##   expect_length(setdiff(colnames(m), c('Intercept', 'sigma')), 0)
  ##   expect_equal(dim(m), c(number_samples, 2))
  ## })

  ## test_that("as.matrix.projection works with clustering", {
  ##   number_clusters <- 3
  ##   p_clust <- project(vs_gauss, solution_terms = solution_terms, number_clusters = number_clusters)
  ##   m <- as.matrix(p_clust)
  ##   expect_length(setdiff(colnames(m), c("Intercept", vs_gauss$solution_terms[solution_terms], "sigma")), 0)
  ##   expect_equal(dim(m), c(number_clusters, length(solution_terms) + 2))
  ## })
}
