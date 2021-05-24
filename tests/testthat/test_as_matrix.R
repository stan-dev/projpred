context("as.matrix.projection")

# Gaussian and binomial reference models without multilevel or additive terms:
if (require(rstanarm)) {
  set.seed(1235)
  n <- 40
  nterms <- 5
  x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
  b <- runif(nterms) - 0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = TRUE)
  offset <- rnorm(n)
  chains <- 2
  seed <- 1235
  iter <- 500
  source(testthat::test_path("helpers", "SW.R"))

  f_gauss <- gaussian()
  df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x %*% b), dis), x = x)
  f_binom <- binomial()
  df_binom <- data.frame(
    y = rbinom(n, weights, f_binom$linkinv(x %*% b)), x = x,
    weights = weights, offset = offset
  )

  SW({
    fit_gauss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          data = df_gauss, family = f_gauss,
                          chains = chains, seed = seed, iter = iter)
    fit_binom <- stan_glm(cbind(y, weights - y) ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          data = df_binom, family = f_binom,
                          weights = weights,
                          chains = chains, seed = seed, iter = iter)
  })

  solution_terms <- c("x.3", "x.5")
  ndraws <- 100
  p_gauss <- project(fit_gauss,
                     solution_terms = solution_terms,
                     ndraws = ndraws)
  SW(p_binom <- project(fit_binom,
                        solution_terms = solution_terms,
                        ndraws = ndraws))

  test_that(paste(
    "as.matrix.projection returns the relevant variables for",
    "gaussian"
  ), {
    m <- as.matrix(p_gauss)
    expect_length(setdiff(
      colnames(m),
      c(
        paste0("b_", c("Intercept", solution_terms)),
        "sigma"
      )
    ), 0)
    expect_equal(dim(m), c(ndraws, length(solution_terms) + 2))
  })

  test_that(paste(
    "as.matrix.projection returns the relevant variables for",
    "binomial"
  ), {
    m <- as.matrix(p_binom)
    expect_length(setdiff(colnames(m),
                          paste0("b_",
                                 c("Intercept",
                                   solution_terms))),
                  0)
    expect_equal(dim(m), c(ndraws, length(solution_terms) + 1))
  })

  test_that("as.matrix.projection works as expected with zero variables", {
    p_novars <- project(fit_gauss,
                        solution_terms = character(),
                        ndraws = ndraws)
    m <- as.matrix(p_novars)
    expect_length(setdiff(colnames(m), c("b_Intercept", "sigma")), 0)
    expect_equal(dim(m), c(ndraws, 2))
  })

  test_that("as.matrix.projection works with clustering", {
    nclusters <- 3
    p_clust <- project(fit_gauss,
                       solution_terms = solution_terms,
                       nclusters = nclusters)
    SW(m <- as.matrix(p_clust))
    expect_length(
      setdiff(
        colnames(m),
        c(
          paste0("b_", c("Intercept", solution_terms)),
          "sigma"
        )
      ),
      0
    )
    expect_equal(dim(m), c(nclusters, length(solution_terms) + 2))
  })
}
