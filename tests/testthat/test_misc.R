context("miscellaneous")

# miscellaneous tests

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
    weights = weights
  )
  f_poiss <- poisson()
  df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x %*% b)), x = x)

  SW(
    fit_gauss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_gauss, data = df_gauss,
      chains = chains, seed = seed, iter = iter
    )
  )
  SW(
    fit_binom <- stan_glm(cbind(y, weights - y) ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_binom, data = df_binom, weights = weights,
      chains = chains, seed = seed, iter = iter
    )
  )
  SW(
    fit_poiss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_poiss, data = df_poiss,
      chains = chains, seed = seed, iter = iter
    )
  )
  fit_list <- list(
    gauss = fit_gauss, binom = fit_binom, poiss = fit_poiss
  )

  test_that(paste(
    "check that the main function calls do not return the",
    "same RNG state every time"
  ), {
    s <- 5

    for (seed in c(130927, NULL)) {
      for (i in seq_along(fit_list)) {
        fit <- fit_list[[i]]

        # varsel
        foo <- varsel(fit, seed = seed)
        r1 <- rnorm(s)
        foo <- varsel(fit, seed = seed)
        r2 <- rnorm(s)
        expect_true(any(r1 != r2))

        # cv_varsel
        SW(foo <- cv_varsel(fit, seed = seed))
        r1 <- rnorm(s)
        SW(foo <- cv_varsel(fit, seed = seed))
        r2 <- rnorm(s)
        expect_true(any(r1 != r2))

        # project
        solution_terms <- c(1, 2)
        foo <- project(fit,
          solution_terms = solution_terms,
          ndraws = 100, seed = seed
        )
        r1 <- rnorm(s)
        foo <- project(fit,
          solution_terms = solution_terms,
          ndraws = 100, seed = seed
        )
        r2 <- rnorm(s)
        expect_true(any(r1 != r2))

        # proj_linpred
        solution_terms <- c(1, 3)
        frame <- cbind(data.frame(x = x)[, solution_terms], weights = weights)
        foo <- proj_linpred(fit, frame,
          solution_terms = solution_terms,
          seed = seed, weightsnew = ~weights
        )
        r1 <- rnorm(s)
        foo <- proj_linpred(fit, frame,
          solution_terms = solution_terms,
          seed = seed, weightsnew = ~weights
        )
        r2 <- rnorm(s)
        expect_true(any(r1 != r2))

        # proj_predict
        solution_terms <- c(1, 3)
        frame <- cbind(data.frame(x = x)[, solution_terms], weights = weights)
        foo <- proj_predict(fit, frame,
          solution_terms = solution_terms, seed = seed
        )
        r1 <- rnorm(s)
        foo <- proj_predict(fit, frame,
          solution_terms = solution_terms, seed = seed
        )
        r2 <- rnorm(s)
        expect_true(any(r1 != r2))
      }
    }
  })


  test_that("check that providing seed has the expected effect", {
    for (seed in c(130927, 1524542)) {
      for (i in seq_along(fit_list)) {
        fit <- fit_list[[i]]

        # varsel
        foo <- varsel(fit, seed = seed)
        bar <- varsel(fit, seed = seed)
        expect_equal(foo, bar)

        # cv_varsel
        SW(foo <- cv_varsel(fit, seed = seed))
        SW(bar <- cv_varsel(fit, seed = seed))
        expect_equal(foo, bar)

        # project
        solution_terms <- c(1, 2)
        foo <- project(fit,
          solution_terms = solution_terms,
          nclusters = 10, seed = seed
        )
        bar <- project(fit,
          solution_terms = solution_terms,
          nclusters = 10, seed = seed
        )
        expect_equal(foo, bar)


        # proj_linpred
        solution_terms <- c(1, 3)
        frame <- cbind(data.frame(x = x)[, solution_terms], weights = weights)
        foo <- proj_linpred(fit, frame,
          solution_terms = solution_terms,
          seed = seed
        )
        bar <- proj_linpred(fit, frame,
          solution_terms = solution_terms,
          seed = seed
        )
        expect_equal(foo, bar)

        # proj_predict
        solution_terms <- c(1, 3)
        frame <- cbind(data.frame(x = x)[, solution_terms], weights = weights)
        foo <- proj_predict(fit, frame,
          solution_terms = solution_terms,
          seed = seed
        )
        bar <- proj_predict(fit, frame,
          solution_terms = solution_terms,
          seed = seed
        )
        expect_equal(foo, bar)
      }
    }
  })
}
