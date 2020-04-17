context("project")


# tests for project

if (require(rstanarm)) {
  seed <- 1235
  set.seed(seed)
  n <- 40
  nv <- 5
  x <- matrix(rnorm(n * nv, 0, 1), n, nv)
  b <- runif(nv) - 0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = TRUE)
  offset <- rnorm(n)
  chains <- 2
  iter <- 500
  source(file.path("helpers", "SW.R"))

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
      family = f_binom,
      data = df_binom, chains = chains, seed = seed, iter = iter
    )
    fit_poiss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_poiss, data = df_poiss,
      chains = chains, seed = seed, iter = iter
    )
  })
  fit_list <- list( # fit_gauss,
    fit_binom, fit_poiss
  )
  vs_list <- lapply(fit_list, varsel, nv_max = nv + 1, verbose = FALSE)


  test_that("object returned by project contains the relevant fields", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nv = 0:nv)
      expect_type(p, "list")
      expect_length(p, nv + 1)

      for (j in 1:length(p)) {
        expect_s3_class(p[[j]], "projection")
        expect_named(p[[j]], c(
          "kl", "weights", "dis", "sub_fit", "solution_terms",
          "p_type", "family", "intercept"
        ),
        ignore.order = TRUE, info = i_inf
        )
        # number of draws should equal to the number of draw weights
        number_samples <- length(p[[j]]$weights)
        ## expect_equal(NROW(as.matrix(p[[j]])), number_samples)
        expect_length(p[[j]]$dis, number_samples)
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
      expect_equal(p[[1]]$weights, p[[nv]]$weights, info = i_inf)
    }
  })

  test_that(paste("project: error when varsel has not been performed for the",
                  "object"), {
    expect_error(
      project(1, xnew = x),
      "is not a variable selection -object"
    )
    expect_error(
      project(fit_gauss, xnew = x),
      "is not a variable selection -object"
    )
  })

  test_that("project: nv is checked", {
    expect_error(
      project(vs_list[[1]], nv = 1000),
      "Cannot perform the projection with 1000 variables"
    )
    expect_error(
      project(vs_list[[1]], nv = -1),
      "must contain non-negative values"
    )
    expect_error(
      project(vs_list[[1]], nv = "a"),
      "must contain non-negative values"
    )
    expect_error(
      project(vs_list[[1]], nv = df_gauss),
      "must contain non-negative values"
    )
  })

  test_that("project: setting nv = NULL has the expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nv = NULL)
      # if only one model size is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_equal(count_terms_chosen(p$solution_terms) - 1,
                   vs_list[[i]]$suggested_size, info = i_inf)
      expect_length(p$solution_terms, vs_list[[i]]$suggested_size)
    }
  })

  test_that("project: setting nv = 0 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nv <- 0
      p <- project(vs_list[[i]], nv = nv)
      # if only one model size is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_equal(count_terms_chosen(p$solution_terms) - 1, nv, info = i_inf)
      expect_length(p$solution_terms, 1)
    }
  })

  test_that("project: setting nv = 3 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nv <- 3
      p <- project(vs_list[[i]], nv = nv)
      # if only one model is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_length(p$solution_terms, nv)
    }
  })

  test_that("project: setting solution_terms to 4 has an expected effect", {
    for (i in 1:length(vs_list)) {
      solution_terms <- 4
      p <- project(vs_list[[i]], solution_terms = solution_terms)
      expect_equivalent(p$solution_terms,
                        vs_list[[i]]$solution_terms[solution_terms])
    }
  })

  test_that("project: setting solution_terms to 1:2 has an expected effect", {
    for (i in 1:length(vs_list)) {
      solution_terms <- 1:2
      p <- project(vs_list[[i]], solution_terms = solution_terms)
      expect_equivalent(p$solution_terms,
                        vs_list[[i]]$solution_terms[solution_terms])
    }
  })

  test_that(paste("project: setting solution_terms to something nonsensical",
                  "returns an error"), {
    # variable selection objects
    expect_error(
      project(vs_list[[1]], solution_terms = 1:10),
      "solution_terms contains an index larger than"
    )
    expect_error(
      project(vs_list[[1]], solution_terms = 17),
      "solution_terms contains an index larger than"
    )

    # fit objects
    expect_error(
      project(fit_list[[1]], solution_terms = 1:10),
      "solution_terms contains an index larger than"
    )
    expect_error(
      project(fit_list[[1]], solution_terms = 17),
      "solution_terms contains an index larger than"
    )
  })

  test_that("project: setting number_samples to 1 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      number_samples <- 1
      p <- project(vs_list[[i]], number_samples = number_samples, nv = nv)
      # expected number of draws
      expect_length(p$weights, number_samples)
      expect_equal(NROW(as.matrix(p)), number_samples)
      expect_equal(p$weights, 1, info = i_inf)
    }
  })

  test_that("project: setting number_samples to 40 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      number_samples <- 40
      p <- project(vs_list[[i]], number_samples = number_samples, nv = nv)
      # expected number of draws
      expect_length(p$weights, number_samples)
      ## expect_equal(NROW(as.matrix(p)), number_samples)
    }
  })

  test_that("project: setting number_clusters to 1 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      number_clusters <- 1
      p <- project(vs_list[[i]], number_clusters = number_clusters, nv = nv)
      # expected number of clusters
      expect_length(p$weights, number_clusters)
      expect_length(p$kl, number_clusters)
    }
  })

  test_that("project: setting number_clusters to 20 has an expected effect", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      number_clusters <- 20
      p <- project(vs_list[[i]], number_clusters = number_clusters, nv = nv)
      # expected number of draws
      expect_length(p$weights, number_clusters)
      expect_length(p$kl, number_clusters)
    }
  })

  test_that(paste("project: setting number_samples or number_clusters to too",
                  "big throws an error"), {
    expect_error(
      project(vs_list[[1]], number_samples = 400000, nv = nv),
      "exceed the number of columns"
    )
    expect_error(
      project(vs_list[[1]], number_clusters = 400000, nv = nv),
      "exceed the number of columns"
    )
  })

  test_that("project: specifying the seed does not cause errors", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nv = nv, seed = seed)
      expect_named(p, c(
        "kl", "weights", "dis", "sub_fit", "solution_terms",
        "p_type", "family", "intercept"
      ),
      ignore.order = TRUE, info = i_inf
      )
    }
  })

  test_that(paste("project: projecting full model onto itself does not change",
                  "results"), {
    tol <- 0.6

    for (i in 1:length(fit_list)) {
      fit <- fit_list[[i]]
      draws <- as.data.frame(fit)
      alpha_ref <- draws$`b_Intercept`
      beta_ref <- draws[, 1 + seq_len(nv), drop = FALSE]
      S <- nrow(draws)
      vs <- varsel(fit)
      proj <- project(vs, solution_terms = 1:nv, seed = seed,
                      number_samples = S)

      # test alpha and beta
      ## coefs <- as.matrix(proj)
      ## dalpha <- max(abs(coefs[, 1] - alpha_ref))
      ## order <- match(colnames(fit_list[[i]]$data), proj$solution_terms)
      ## order <- order[!is.na(order)]
      ## dbeta <- max(abs(coefs[, -1, drop = FALSE][, order] - beta_ref))
      ## expect_lt(dalpha, tol)
      ## expect_lt(dbeta, tol)
    }
  })

  test_that("project: works as expected from a cvsel object", {
    SW({
      cvs <- cv_varsel(fit_binom, nv_max = 3, verbose = FALSE)
      p <- project(cvs, nv = 3)
    })
    expect_length(p$solution_terms, 3)
  })
}
