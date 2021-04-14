context("project")


# tests for project

if (require(rstanarm)) {
  seed <- 1235
  set.seed(seed)
  n <- 40
  nterms <- 5
  ndraws <- 1
  ndraws_pred <- 5
  x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
  b <- runif(nterms) - 0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = TRUE)
  offset <- rnorm(n)
  chains <- 2
  iter <- 500
  source(testthat::test_path("helpers", "SW.R"))

  f_gauss <- gaussian()
  df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x %*% b), dis), x = x)
  f_binom <- binomial()
  df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x %*% b)),
                         x = x, weights = weights)
  f_poiss <- poisson()
  df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x %*% b)), x = x)

  SW({
    fit_gauss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          family = f_gauss, data = df_gauss, QR = TRUE,
                          weights = weights, offset = offset,
                          chains = chains, seed = seed, iter = iter
    )
    fit_binom <- stan_glm(cbind(y, weights - y) ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          family = f_binom, weights = weights,
                          data = df_binom, chains = chains, seed = seed,
                          iter = iter
    )
    fit_poiss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          family = f_poiss, data = df_poiss,
                          chains = chains, seed = seed, iter = iter
    )
    fit_list <- list( # fit_gauss,
      fit_binom, fit_poiss
    )
    vs_list <- lapply(fit_list, varsel,
                      nterms_max = nterms + 1,
                      verbose = FALSE
    )
  })

  test_that("object returned by project contains the relevant fields", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nterms = 0:nterms)
      expect_type(p, "list")
      expect_length(p, nterms + 1)

      for (j in 1:length(p)) {
        expect_s3_class(p[[j]], "projection")
        expect_named(p[[j]], c(
          "kl", "weights", "dis", "sub_fit", "solution_terms",
          "p_type", "family", "intercept", "extract_model_data",
          "refmodel"
        ),
        ignore.order = TRUE, info = i_inf
        )
        # number of ndraws should equal to the number of draw weights
        ndraws <- length(p[[j]]$weights)
        ## expect_equal(NROW(as.matrix(p[[j]])), ndraws)
        expect_length(p[[j]]$dis, ndraws)
        # j:th element should have j variables, including the intercept
        expect_length(p[[j]]$solution_terms, max(j - 1, 1))
        # family kl
        expect_equal(p[[j]]$family, vs_list[[i]]$family,
                     info = i_inf
        )
      }
      # kl should be non-increasing on training data
      klseq <- sapply(p, function(x) sum(x$kl))
      # remove intercept from the comparison
      expect_true(all(diff(klseq)[-1] - 1e-1 < 0), info = i_inf)
      ## expect_equal(klseq, cummin(klseq), info = i_inf)

      # all submodels should use the same clustering
      expect_equal(p[[1]]$weights, p[[nterms]]$weights, info = i_inf)
    }
  })

  test_that(paste(
    "project: error when varsel has not been performed for the",
    "object"
  ), {
    expect_error(
      project(1, newdata = x),
      "is not a variable selection -object"
    )
    expect_error(
      project(fit_gauss, newdata = x),
      "is not a variable selection -object"
    )
  })

  test_that("project: nterms is checked", {
    expect_error(
      project(vs_list[[1]], nterms = 1000),
      "Cannot perform the projection with 1000 variables"
    )
    expect_error(
      project(vs_list[[1]], nterms = -1),
      "must contain non-negative values"
    )
    expect_error(
      project(vs_list[[1]], nterms = "a"),
      "must contain non-negative values"
    )
    expect_error(
      project(vs_list[[1]], nterms = df_gauss),
      "must contain non-negative values"
    )
  })

  test_that("project: setting nterms = NULL has the expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nterms = NULL)
      # if only one model size is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_equal(count_terms_chosen(p$solution_terms) - 1,
                   vs_list[[i]]$suggested_size, info = i_inf)
      expect_length(p$solution_terms, vs_list[[i]]$suggested_size)
    }
  })

  test_that("project: setting nterms = 0 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nterms <- 0
      p <- project(vs_list[[i]], nterms = nterms)
      # if only one model size is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_equal(count_terms_chosen(p$solution_terms) - 1, nterms,
                   info = i_inf)
      expect_length(p$solution_terms, 1)
    }
  })

  test_that("project: setting nterms = 3 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nterms <- 3
      p <- project(vs_list[[i]], nterms = nterms)
      # if only one model is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_length(p$solution_terms, nterms)
    }
  })

  test_that("project: setting solution_terms to 4 has an expected effect", {
    for (i in 1:length(vs_list)) {
      solution_terms <- 4
      p <- project(vs_list[[i]],
                   solution_terms = vs_list[[i]]$solution_terms[solution_terms])
      expect_equivalent(p$solution_terms,
                        vs_list[[i]]$solution_terms[solution_terms])
    }
  })

  test_that("project: setting solution_terms to 1:2 has an expected effect", {
    for (i in 1:length(vs_list)) {
      solution_terms <- 1:2
      p <- project(vs_list[[i]],
                   solution_terms = vs_list[[i]]$solution_terms[solution_terms])
      expect_equivalent(p$solution_terms,
                        vs_list[[i]]$solution_terms[solution_terms])
    }
  })

  ## test_that(paste(
  ##   "project: setting solution_terms to something nonsensical",
  ##   "returns an error"
  ## ), {
  ##   # variable selection objects
  ##   expect_error(
  ##     project(vs_list[[1]],
  ##             solution_terms = vs_list[[1]]$solution_terms[1:10]),
  ##     "solution_terms contains an index larger than"
  ##   )
  ##
  ##   # fit objects
  ##   expect_error(
  ##     SW(project(fit_list[[1]], solution_terms = 1:10)),
  ##     "solution_terms contains an index larger than"
  ##   )
  ##   expect_error(
  ##     SW(project(fit_list[[1]], solution_terms = 17)),
  ##     "solution_terms contains an index larger than"
  ##   )
  ## })

  test_that("project: setting ndraws to 1 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      ndraws <- 1
      p <- project(vs_list[[i]], ndraws = ndraws, nterms = nterms)
      # expected number of ndraws
      expect_length(p$weights, ndraws)
      expect_equal(NROW(as.matrix(p)), ndraws)
      expect_equal(p$weights, 1, info = i_inf)
    }
  })

  test_that("project: setting ndraws to 40 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      ndraws <- 40
      p <- project(vs_list[[i]], ndraws = ndraws, nterms = nterms)
      # expected number of ndraws
      expect_length(p$weights, ndraws)
      ## expect_equal(NROW(as.matrix(p)), ndraws)
    }
  })

  test_that("project: setting nclusters to 1 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nclusters <- 1
      p <- project(vs_list[[i]], nclusters = nclusters, nterms = nterms)
      # expected number of clusters
      expect_length(p$weights, nclusters)
    }
  })

  test_that("project: setting nclusters to 20 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nclusters <- 20
      p <- project(vs_list[[i]], nclusters = nclusters, nterms = nterms)
      # expected number of ndraws
      expect_length(p$weights, nclusters)
      expect_length(p$sub_fit, nclusters)
    }
  })

  test_that(paste(
    "project: setting ndraws or nclusters too big causes them to be cut off at",
    "the number of posterior draws in the reference model"
  ), {
    p <- project(vs_list[[1]], ndraws = 400000, nterms = nterms)
    expect_length(p$weights, nrow(as.matrix(fit_list[[1]])))
    expect_length(p$sub_fit, nrow(as.matrix(fit_list[[1]])))
    expect_length(p$dis, nrow(as.matrix(fit_list[[1]])))
    p <- project(vs_list[[1]], nclusters = 400000, nterms = nterms)
    expect_length(p$weights, nrow(as.matrix(fit_list[[1]])))
    expect_length(p$sub_fit, nrow(as.matrix(fit_list[[1]])))
    expect_length(p$dis, nrow(as.matrix(fit_list[[1]])))
  })

  test_that("project: specifying the seed does not cause errors", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nterms = nterms, seed = seed)
      expect_named(p, c(
        "kl", "weights", "dis", "sub_fit", "solution_terms",
        "p_type", "family", "intercept", "extract_model_data",
        "refmodel"
      ),
      ignore.order = TRUE, info = i_inf
      )
    }
  })

  test_that(paste(
    "project: projecting full model onto itself does not change",
    "results"
  ), {
    tol <- 1e-3

    for (i in 1:length(fit_list)) {
      fit <- fit_list[[i]]
      draws <- as.data.frame(fit)
      alpha_ref <- draws[, "(Intercept)"]
      beta_ref <- draws[, 1 + seq_len(nterms), drop = FALSE]
      S <- nrow(draws)
      SW(vs <- varsel(fit))
      proj <- project(vs,
                      solution_terms = vs$solution_terms[1:nterms],
                      seed = seed, ndraws = S
      )

      # test alpha and beta
      coefs <- as.matrix(proj)
      dalpha <- abs(mean(coefs[, 1]) - mean(alpha_ref))
      order <- match(colnames(fit_list[[i]]$data), proj$solution_terms)
      order <- order[!is.na(order)]
      dbeta <- max(abs(colMeans(coefs[, -1, drop = FALSE][, order])
                       - colMeans(beta_ref)))
      expect_lt(dalpha, tol)
      expect_lt(dbeta, tol)
    }
  })

  test_that("project: works as expected from a vsel object", {
    SW({
      cvs <- cv_varsel(fit_binom,
                       nterms_max = 3, verbose = FALSE, ndraws = ndraws,
                       ndraws_pred = ndraws_pred
      )
      p <- project(cvs, nterms = 3)
    })
    expect_length(p$solution_terms, 3)
  })
}
