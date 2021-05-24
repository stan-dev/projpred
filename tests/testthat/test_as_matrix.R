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

  settings_list <- list(
    gauss = list(
      fitobj = fit_gauss,
      solution_terms_list = list(character(), c("x.3", "x.5")),
      ndraws_list = list(100, 3, 1)
    ),
    binom = list(
      fitobj = fit_binom,
      solution_terms_list = list(c("x.3", "x.5")),
      ndraws_list = list(100)
    )
  )

  for (fam_type in settings_list) {
    for (solution_terms in fam_type$solution_terms_list) {
      for (ndraws in fam_type$ndraws_list) {
        # Expected warning (more precisely: regexp which is matched against the
        # warning; NA means no warning) for project() and family-specific
        # parameters:
        if (family(fam_type$fitobj)$family == "gaussian") {
          warn_prj_expect <- NA
          npars_fam <- "sigma"
        } else if (family(fam_type$fitobj)$family == "binomial") {
          # For the binomial family with > 1 trials, we expect a warning (see
          # GitHub issue #136):
          warn_prj_expect <- paste("Using formula\\(x\\) is deprecated when x",
                                   "is a character vector of length > 1")
          npars_fam <- character()
        }

        expect_warning(prj <- project(fam_type$fitobj,
                                      solution_terms = solution_terms,
                                      ndraws = ndraws),
                       warn_prj_expect)

        # Expected warning (more precisely: regexp which is matched against the
        # warning; NA means no warning) for as.matrix.projection():
        if (ndraws > 20) {
          warn_prjmat_expect <- NA
        } else {
          # Clustered projection, so we expect a warning:
          warn_prjmat_expect <- "the clusters might have different weights"
        }

        expect_warning(m <- as.matrix(prj), warn_prjmat_expect)

        test_that("as.matrix.projection()'s output structure is correct", {
          expect_equal(
            dim(m),
            c(ndraws, length(solution_terms) + 1 + length(npars_fam))
          )
          expect_identical(
            colnames(m),
            c(paste0("b_", c("Intercept", solution_terms)), npars_fam)
          )
        })
      }
    }
  }
}
