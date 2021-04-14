context("varsel")

# tests for varsel and cv_varsel

if (require(rstanarm)) {
  seed <- 1235
  set.seed(seed)
  n <- 50
  nterms <- 5
  ndraws <- 1
  ndraws_pred <- 5
  x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
  b <- runif(nterms) - 0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = TRUE)
  chains <- 2
  iter <- 500
  offset <- rnorm(n)
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
                          family = f_gauss, data = df_gauss,
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
    fit_glmer <- stan_glmer(mpg ~ wt + (1 | cyl),
                            data = mtcars,
                            chains = chains, seed = seed, iter = iter
    )
    fit_list <- list(
      gauss = fit_gauss, binom = fit_binom, poiss = fit_poiss
    )

    formula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5
    vsf <- function(x, m)
      varsel(x,
             method = m, nterms_max = nterms + 1, verbose = FALSE,
             ndraws = ndraws, ndraws_pred = ndraws_pred
      )
    vs_list <- list(
      l1 = lapply(fit_list, vsf, "L1"),
      fs = lapply(fit_list, vsf, "forward")
    )
  })

  extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                 orhs = NULL, extract_y = FALSE) {
    if (!is.null(object)) {
      formula <- formula(object)
      tt <- extract_terms_response(formula)
      response_name <- tt$response
    } else {
      response_name <- NULL
    }

    if (is.null(newdata)) {
      newdata <- object$data
    }

    resp_form <- NULL
    if (is.null(object)) {
      if ("weights" %in% colnames(newdata))
        wrhs <- ~ weights
      if ("offset" %in% colnames(newdata))
        orhs <- ~ offset
      if ("y" %in% colnames(newdata))
        resp_form <- ~ y
    }

    args <- nlist(object, newdata, wrhs, orhs, resp_form)
    return(do_call(.extract_model_data, args))
  }

  ref_gauss <- init_refmodel(
    object = NULL, df_gauss, formula,
    family = f_gauss, extract_model_data = extract_model_data
  )
  ref_binom <- init_refmodel(
    object = NULL, df_binom, formula,
    family = f_binom, extract_model_data = extract_model_data
  )
  ref_list <- list(ref_gauss = ref_gauss, ref_binom = ref_binom)
  vsref_list <- list(
    l1 = lapply(ref_list, vsf, "L1"),
    fs = lapply(ref_list, vsf, "forward")
  )

  test_that('varsel returns an object of type "vsel"', {
    for (i in seq_len(length(vs_list))) {
      for (j in seq_len(length(vs_list[[i]]))) {
        expect_s3_class(vs_list[[i]][[j]], "vsel")
      }
    }
  })

  test_that("object returned by varsel contains the relevant fields", {
    for (i in seq_len(length(vs_list))) {
      i_inf <- names(vs_list)[i]
      for (j in seq_len(length(vs_list[[i]]))) {
        j_inf <- names(vs_list[[i]])[j]
        # refmodel seems legit
        expect_s3_class(vs_list[[i]][[j]]$refmodel, "refmodel")
        # solution_terms seems legit
        expect_length(vs_list[[i]][[j]]$solution_terms, nterms)
        expect_true(all(!is.na(match(
          colnames(fit_gauss$data[, -1]),
          vs_list[[i]][[j]]$solution_terms
        ))),
        info = paste(i_inf, j_inf)
        )
        # kl seems legit
        expect_length(vs_list[[i]][[j]]$kl, nterms + 1)
        # decreasing
        expect_equal(vs_list[[i]][[j]]$kl,
                     cummin(vs_list[[i]][[j]]$kl),
                     tolerance = 2e-2,
                     info = paste(i_inf, j_inf)
        )
        # summaries seems legit
        expect_named(vs_list[[i]][[j]]$summaries, c("sub", "ref"),
                     info = paste(i_inf, j_inf)
        )
        expect_length(vs_list[[i]][[j]]$summaries$sub, nterms + 1)
        expect_named(vs_list[[i]][[j]]$summaries$sub[[1]],
                     c("mu", "lppd", "draws"),
                    ignore.order = TRUE,
                    info = paste(i_inf, j_inf)
        )
        expect_named(vs_list[[i]][[j]]$summaries$ref,
                     c("mu", "lppd", "draws"),
                     ignore.order = TRUE,
                     info = paste(i_inf, j_inf)
        )
        # family seems legit
        expect_equal(vs_list[[i]][[j]]$family$family,
                     vs_list[[i]][[j]]$family$family,
                     info = paste(i_inf, j_inf)
        )
        expect_equal(vs_list[[i]][[j]]$family$link,
                     vs_list[[i]][[j]]$family$link,
                     info = paste(i_inf, j_inf)
        )
        expect_true(length(vs_list[[i]][[j]]$family) >=
                      length(vs_list[[i]][[j]]$family$family),
                    info = paste(i_inf, j_inf)
        )
      }
    }
  })

  test_that("search method is valid", {
    expect_error(
      varsel(fit_gauss, method = "k-fold"),
      "Unknown search method"
    )
  })

  test_that("nterms_max has an effect on varsel for gaussian models", {
    vs1 <- varsel(fit_gauss, method = "forward", nterms_max = 3,
                  verbose = FALSE)
    expect_length(vs1$solution_terms, 3)
  })

  test_that("nterms_max has an effect on varsel for non-gaussian models", {
    SW(vs1 <- varsel(fit_binom, method = "forward", nterms_max = 3,
                     verbose = FALSE))
    expect_length(vs1$solution_terms, 3)
  })

  test_that("specifying the number of clusters has an expected effect", {
    SW(vs <- varsel(fit_binom, method = "forward", nterms_max = 3,
                    nclusters = 10))
    expect_length(vs$solution_terms, 3)
  })

  test_that("specifying d_test has the expected effect", {
    refmodel_ <- vs_list[[1]][[1]]$refmodel
    d_test <- list(
      data = fit_gauss$data, y = refmodel_$y,
      test_points = seq_along(refmodel_$y),
      weights = refmodel_$wobs,
      type = "test"
    )
    vs <- varsel(fit_gauss, d_test = d_test, nterms_max = 3)
    expect_length(vs$solution_terms, 3)
  })

  test_that("Having something else than stan_glm as the fit throws an error", {
    expect_error(varsel(rnorm(5), verbose = FALSE),
                 regexp = "no applicable method")
  })


  test_that("varsel: adding more regularization has an expected effect", {
    regul <- c(1e-6, 1e-3, 1e-1, 1e1, 1e4)
    for (i in 1:length(fit_list)) {
      nonzeros <- rep(0, length(regul))
      msize <- 3
      for (j in 1:length(regul)) {
        SW(
          vsel <- varsel(fit_list[[i]], regul = regul[j], nterms_max = 6)
        )
        x <- vsel$search_path$sub_fits[[6]]
        sol <- rbind(x$alpha, x$beta)
        nonzeros[j] <- length(which(sol != 0))
      }
      for (j in 1:(length(regul) - 1)) {
        expect_gte(nonzeros[j], nonzeros[j + 1])
      }
    }
  })

  test_that("varsel: length of the penalty vector is checked", {
    vsf <- function(obj, penalty) {
      varsel(obj,
             method = "L1", nterms_max = nterms + 1,
             verbose = FALSE, penalty = penalty,
             ndraws = ndraws, ndraws_pred = ndraws_pred
      )
    }
    expect_error(vsf(fit_list$gauss, rep(1, nterms + 10)))
    expect_error(vsf(fit_list$gauss, 1))
  })

  test_that(paste(
    "varsel: specifying penalties for variables has an expected",
    "effect"
  ), {
    penalty <- rep(1, nterms)
    ind_zeropen <- c(3, 5) # a few variables without cost
    ind_infpen <- c(2) # one variable with infinite penalty
    penalty[ind_zeropen] <- 0
    penalty[ind_infpen] <- Inf
    vsf <- function(obj)
      varsel(obj,
             method = "L1", nterms_max = nterms, verbose = FALSE,
             penalty = penalty, ndraws = ndraws, ndraws_pred = ndraws_pred
      )
    SW(vs_list_pen <- lapply(fit_list, vsf))
    for (i in seq_along(vs_list_pen)) {
      # check that the variables with no cost are selected first and the ones
      # with inf penalty last
      sub_fits <- vs_list_pen[[i]]$search_path$sub_fits
      sdiff <- length(which(sub_fits[[nterms + 1]]$beta == 0))
      expect_gte(sdiff, 1)
    }
  })


  # -------------------------------------------------------------
  context("cv_varsel")

  cvsf <- function(x, m, cvm, K = NULL, ...) {
    cv_varsel(x, method = m, cv_method = cvm, nterms_max = nterms, K = K,
              ndraws = ndraws, ndraws_pred = ndraws_pred, verbose = FALSE,
              ...)
  }

  if (Sys.getenv("NOT_CRAN") == "true") {
    SW({
      cvs_list <- list(
        l1 = lapply(fit_list, cvsf, "L1", "LOO"),
        fs = lapply(fit_list, cvsf, "forward", "LOO", validate_search = FALSE)
      )

      # without weights/offset because kfold does not support them currently
      # test only with one family to make the tests faster

      # the chains, seed and iter arguments to the rstanarm functions here must
      # be specified directly rather than through a variable (eg, seed = 1235
      # instead of seed = seed), otherwise when the calls are evaluated in
      # refmodel$cvfun() they may not be found in the evaluation frame of the
      # calling function, causing the test to fail
      glm_simp <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
                           family = poisson(), data = df_poiss,
                           chains = 2, seed = 1235, iter = 400
      )
      lm_simp <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          data = df_gauss, family = gaussian(),
                          chains = 2, seed = 1235, iter = 400
      )
      simp_list <- list(glm = lm_simp)

      cv_kf_list <- list(
        l1 = lapply(simp_list, cvsf, "L1", "kfold", K = 2),
        fs = lapply(simp_list, cvsf, "forward", "kfold", K = 2)
      )

      # LOO cannot be performed without a genuine probabilistic model
      cvsref_list <- list(
        l1 = lapply(ref_list, cvsf, "L1", "kfold"),
        fs = lapply(ref_list, cvsf, "forward", "kfold")
      )
    })

    test_that('cv_varsel returns an object of type "vsel"', {
      for (i in seq_len(length(cvs_list))) {
        for (j in seq_len(length(cvs_list[[i]]))) {
          expect_s3_class(cvs_list[[i]][[j]], "vsel")
        }
      }
    })

    test_that("object returned by cv_varsel contains the relevant fields", {
      for (i in seq_len(length(cvs_list))) {
        i_inf <- names(cvs_list)[i]
        for (j in seq_len(length(cvs_list[[i]]))) {
          j_inf <- names(cvs_list[[i]])[j]
          # refmodel seems legit
          expect_s3_class(cvs_list[[i]][[j]]$refmodel, "refmodel")
          # solution_terms seems legit
          expect_length(cvs_list[[i]][[j]]$solution_terms, nterms)
          expect_true(all(!is.na(match(
            colnames(fit_gauss$data[, -1]),
            cvs_list[[i]][[j]]$solution_terms
          ))),
          info = paste(i_inf, j_inf)
          )
          # kl seems legit
          expect_length(cvs_list[[i]][[j]]$kl, nterms + 1)
          # decreasing
          expect_equal(cvs_list[[i]][[j]]$kl,
                       cummin(cvs_list[[i]][[j]]$kl),
                       tolerance = 23e-2,
                       info = paste(i_inf, j_inf)
          )
          # summaries seems legit
          expect_named(cvs_list[[i]][[j]]$summaries, c("sub", "ref"),
                       info = paste(i_inf, j_inf)
          )
          expect_length(cvs_list[[i]][[j]]$summaries$sub, nterms + 1)
          expect_named(cvs_list[[i]][[j]]$summaries$sub[[1]],
                       c("lppd", "mu", "w", "draws"),
                       info = paste(i_inf, j_inf)
          )
          expect_named(cvs_list[[i]][[j]]$summaries$ref, c("lppd", "mu", "draws"),
                       info = paste(i_inf, j_inf)
          )
          # family seems legit
          expect_equal(cvs_list[[i]][[j]]$family$family,
                       cvs_list[[i]][[j]]$family$family,
                       info = paste(i_inf, j_inf)
          )
          expect_equal(cvs_list[[i]][[j]]$family$link,
                       cvs_list[[i]][[j]]$family$link,
                       info = paste(i_inf, j_inf)
          )
          expect_true(length(cvs_list[[i]][[j]]$family) >=
                        length(cvs_list[[i]][[j]]$family$family),
                      info = paste(i_inf, j_inf)
          )
        }
      }
    })

    test_that("nterms_max has an effect on cv_varsel for gaussian models", {
      suppressWarnings(
        vs1 <- cv_varsel(fit_gauss,
                         method = "forward", nterms_max = 3,
                         verbose = FALSE, ndraws = ndraws,
                         ndraws_pred = ndraws_pred,
                         validate_search = FALSE
        )
      )
      expect_length(vs1$solution_terms, 3)
    })

    test_that("nterms_max has an effect on cv_varsel for non-gaussian models", {
      suppressWarnings(
        vs1 <- cv_varsel(fit_binom,
                         method = "forward", nterms_max = 3,
                         verbose = FALSE, ndraws = ndraws,
                         ndraws_pred = ndraws_pred,
                         validate_search = FALSE
        )
      )
      expect_length(vs1$solution_terms, 3)
    })

    test_that("nloo works as expected", {
      expect_error(
        SW(
          cv_varsel(fit_gauss,
                    cv_method = "LOO", nloo = -1, ndraws = ndraws,
                    ndraws_pred = ndraws_pred, validate_search = FALSE
          )
        ),
        "must be at least 1"
      )
      SW({
        expect_equal(
          cv_varsel(fit_gauss,
                    cv_method = "LOO", nterms_max = nterms, seed = seed,
                    nloo = NULL, ndraws = ndraws, ndraws_pred = ndraws_pred
          ),
          cv_varsel(fit_gauss,
                    cv_method = "LOO", nterms_max = nterms, seed = seed,
                    nloo = 1000, ndraws = ndraws, ndraws_pred = ndraws_pred
          )
        )

        # nloo less than number of observations
        out <- cv_varsel(fit_gauss,
                         cv_method = "LOO", nloo = 20, verbose = FALSE,
                         ndraws = ndraws, ndraws_pred = ndraws_pred
        )
        expect_equal(sum(!is.na(out$summaries$sub[[1]]$lppd)), 20)
      })
    })

    test_that("the validate_search option works as expected", {
      SW({
        vs1 <- cv_varsel(fit_gauss,
                         validate_search = FALSE,
                         ndraws = ndraws, ndraws_pred = ndraws_pred
        )
        vs2 <- cv_varsel(fit_gauss,
                         validate_search = TRUE,
                         ndraws = ndraws, ndraws_pred = ndraws_pred
        )
      })
      expect_true(all(summary(vs1)$selection$elpd >=
                        summary(vs2)$selection$elpd))
    })

    test_that(paste(
      "Having something else than stan_glm as the fit throws an error"
    ), {
      expect_error(cv_varsel(rnorm(5), verbose = FALSE),
                   regexp = "no applicable method"
      )
    })

    test_that(paste(
      "object returned by cv_varsel, kfold contains the relevant",
      "fields"
    ), {
      for (i in seq_len(length(cv_kf_list))) {
        i_inf <- names(cv_kf_list)[i]
        for (j in seq_len(length(cv_kf_list[[i]]))) {
          j_inf <- names(cv_kf_list[[i]])[j]
          # solution_terms seems legit
          expect_length(cv_kf_list[[i]][[j]]$solution_terms, nterms)
          expect_true(all(!is.na(match(
            colnames(fit_gauss$data[, -1]),
            cv_kf_list[[i]][[j]]$solution_terms
          ))),
          info = paste(i_inf, j_inf)
          )
          # kl seems legit
          expect_length(cv_kf_list[[i]][[j]]$kl, nterms + 1)

          # decreasing
          expect_equal(cv_kf_list[[i]][[j]]$kl[-1],
                       cummin(cv_kf_list[[i]][[j]]$kl[-1]),
                       info = paste(i_inf, j_inf),
                       tolerance = 24e-2
          )

          # summaries seems legit
          expect_named(cv_kf_list[[i]][[j]]$summaries, c("sub", "ref"),
                       info = paste(i_inf, j_inf)
          )
          expect_length(cv_kf_list[[i]][[j]]$summaries$sub, nterms + 1)
          expect_named(cv_kf_list[[i]][[j]]$summaries$sub[[1]],
                       c("mu", "lppd", "w", "draws"),
                       ignore.order = TRUE, info = paste(i_inf, j_inf)
          )
          expect_named(cv_kf_list[[i]][[j]]$summaries$ref, c("mu", "lppd", "draws"),
                       ignore.order = TRUE, info = paste(i_inf, j_inf)
          )
          # family seems legit
          expect_equal(cv_kf_list[[i]][[j]]$family$family,
                       cv_kf_list[[i]][[j]]$family$family,
                       info = paste(i_inf, j_inf)
          )
          expect_equal(cv_kf_list[[i]][[j]]$family$link,
                       cv_kf_list[[i]][[j]]$family$link,
                       info = paste(i_inf, j_inf)
          )
          expect_true(length(cv_kf_list[[i]][[j]]$family) >=
                        length(cv_kf_list[[i]][[j]]$family$family),
                      info = paste(i_inf, j_inf)
          )
          # pct_solution_terms_cv seems legit
          expect_equal(dim(cv_kf_list[[i]][[j]]$pct_solution_terms_cv),
                       c(nterms, nterms + 1),
                       info = paste(i_inf, j_inf)
          )
          expect_true(all(
            cv_kf_list[[i]][[j]]$pct_solution_terms_cv[, -1] <= 1 &
              cv_kf_list[[i]][[j]]$pct_solution_terms_cv[, -1] >= 0
          ),
          info = paste(i_inf, j_inf)
          )
          expect_equal(cv_kf_list[[i]][[j]]$pct_solution_terms_cv[, 1],
                       1:nterms,
                       info = paste(i_inf, j_inf)
          )
          expect_equal(colnames(cv_kf_list[[i]][[j]]$pct_solution_terms_cv),
                       c("size", cv_kf_list[[i]][[j]]$solution_terms),
                       info = paste(i_inf, j_inf)
          )
        }
      }
    })

    test_that("cross-validation method is valid", {
      expect_error(
        cv_varsel(fit_gauss, cv_method = "k-fold"),
        "Unknown cross-validation method"
      )
    })

    test_that("K is valid for cv_method='kfold'", {
      expect_error(
        cv_varsel(glm_simp, cv_method = "kfold", K = 1),
        "must be at least 2"
      )
      expect_error(
        cv_varsel(glm_simp, cv_method = "kfold", K = 1000),
        "cannot exceed n"
      )
      expect_error(
        cv_varsel(glm_simp, cv_method = "kfold", K = c(4, 9)),
        "a single integer value"
      )
      expect_error(
        cv_varsel(glm_simp, cv_method = "kfold", K = "a"),
        "a single integer value"
      )
      expect_error(
        cv_varsel(glm_simp, cv_method = "kfold", K = df_poiss),
        "a single integer value"
      )
    })

    test_that("providing k_fold works", {
      out <- SW({
        k_fold <- kfold(glm_simp, K = 2, save_fits = TRUE)
        folds <- seq_len(nrow(glm_simp$data))
        for (K in seq_len(2)) {
          folds[as.numeric(rownames(k_fold$fit[[K]]$data))] <- K
        }
        attr(k_fold, "folds") <- folds
        fit_cv <- cv_varsel(glm_simp,
                            cv_method = "kfold", cvfits = k_fold,
                            ndraws = ndraws, ndraws_pred = ndraws_pred,
                            verbose = FALSE
        )
      })
      expect_false(any(grepl("k_fold not provided", out)))
      expect_length(fit_cv$solution_terms, nterms)

      # kl seems legit
      expect_length(fit_cv$kl, nterms + 1)

      # decreasing
      expect_equal(fit_cv$kl, cummin(fit_cv$kl), tolerance = 1e-3)

      # summaries seems legit
      expect_named(fit_cv$summaries, c("sub", "ref"))
      expect_length(fit_cv$summaries$sub, nterms + 1)
      expect_named(fit_cv$summaries$sub[[1]], c("mu", "lppd", "w", "draws"),
                   ignore.order = TRUE
      )
      expect_named(fit_cv$summaries$ref, c("mu", "lppd", "draws"),
                   ignore.order = TRUE
      )
      # family seems legit
      expect_equal(
        fit_cv$family$family,
        fit_cv$family$family
      )
      expect_equal(fit_cv$family$link, fit_cv$family$link)
      expect_true(length(fit_cv$family) >= length(fit_cv$family$family))
      # pct_solution_terms_cv seems legit
      expect_equal(dim(fit_cv$pct_solution_terms_cv), c(nterms, nterms + 1))
      expect_true(all(fit_cv$pct_solution_terms_cv[, -1] <= 1 &
                        fit_cv$pct_solution_terms_cv[, -1] >= 0))

      expect_equal(fit_cv$pct_solution_terms_cv[, 1], 1:nterms)
      expect_equal(
        colnames(fit_cv$pct_solution_terms_cv),
        c("size", fit_cv$solution_terms)
      )
    })
  }

  # -------------------------------------------------------------
  context("summary")

  valid_stats_all <- c("elpd", "mlpd")
  valid_stats_gauss_only <- c("mse", "rmse")
  valid_stats_binom_only <- c("acc", "auc")
  valid_stats_gauss <- c(valid_stats_all, valid_stats_gauss_only)
  valid_stats_binom <- c(valid_stats_all, valid_stats_binom_only)
  vs_funs <- c(summary, plot, suggest_size)

  ## test_that("invalid objects are rejected", {
  ##   for (fun in vs_funs) {
  ##     expect_error(fun(NULL), "is not a variable selection object")
  ##     expect_error(fun(fit_gauss), "is not a variable selection object")
  ##   }
  ## })

  ## test_that("invalid stats are rejected", {
  ##   for (fun in vs_funs) {
  ##     expect_error(fun(vs_list[[1]][["gauss"]], stat = NULL),
  ##                  "specified as NULL")
  ##     expect_error(fun(vs_list[[1]][["gauss"]], stat = NA),
  ##                  "not recognized")
  ##     expect_error(fun(vs_list[[1]][["gauss"]], stat = "zzz"),
  ##                  "not recognized")
  ##     expect_error(fun(vs_list[[1]][["gauss"]], stat = "acc"),
  ##                  "available only for the binomial family")
  ##     expect_error(
  ##       fun(vs_list[[1]][["gauss"]], stat = "auc"),
  ##       "available only for the binomial family"
  ##     )
  ##   }
  ## })

  ## test_that("invalid 'baseline' arguments are rejected", {
  ##   expect_error(
  ##     summary(vs_list[[1]][["gauss"]], baseline = "zzz"),
  ##     "Argument 'baseline' must be either 'ref' or 'best'"
  ##   )
  ## })

  test_that("summary output seems legit", {
    skip_on_cran()
    for (i in seq_along(cvs_list)) {
      for (j in seq_along(cvs_list[[i]])) {
        cvs <- cvs_list[[i]][[j]]
        if (cvs$family$family == "gaussian") {
          stats_str <- valid_stats_gauss
        } else if (cvs$family$family == "binomial") {
          stats_str <- valid_stats_binom
        } else {
          stats_str <- valid_stats_all
        }
        cv_method <- cvs_list[[i]][[j]]$cv_method
        stats <- summary(cvs,
                         stats = stats_str,
                         type = c("mean", "lower", "upper", "se")
        )$selection
        expect_true(nrow(stats) == nterms + 1)
        expect_true(all(c(
          "size", "solution_terms", paste0(stats_str, "_", tolower(cv_method)),
          paste0(stats_str, "_", c("se", "upper", "lower"))
        ) %in% names(stats)))
        expect_true(all(stats[, paste0("mlpd_", tolower(cv_method))] >
                          stats[, "mlpd_lower"]))
        expect_true(all(stats[, paste0("mlpd_", tolower(cv_method))] <
                          stats[, "mlpd_upper"]))
      }
    }
  })

  test_that("summary works with reference models", {
    for (i in seq_along(vsref_list)) {
      for (j in seq_along(vsref_list[[i]])) {
        vs <- vsref_list[[i]][[j]]
        if (vs$family$family == "gaussian") {
          stats_str <- valid_stats_gauss
        } else {
          stats_str <- valid_stats_binom
        }
        stats <- summary(vs, stats = stats_str)$selection
        expect_true(is.data.frame(stats))
      }
    }
  })

  test_that("print works as expected", {

    skip_on_cran()
    # default rounding
    expect_output(out <- print(vs_list[[1]][[1]]))
    expect_equal(out$selection$elpd, round(out$selection$elpd, 2),
                 tolerance = 1e-3
    )
    expect_output(out <- print(cvs_list[[1]][[1]]))
    expect_equal(out$selection$elpd, round(out$selection$elpd, 2),
                 tolerance = 1e-3
    )

    # rounding to 4 decimal places
    expect_output(out <- print(vs_list[[1]][[1]], digits = 4))
    expect_equal(out$selection$elpd, round(out$selection$elpd, 4),
                 tolerance = 1e-3
    )
    expect_output(out <- print(cvs_list[[1]][[1]], digits = 4))
    expect_equal(out$selection$elpd, round(out$selection$elpd, 4),
                 tolerance = 1e-3
    )
    # options to summary
    expect_output(out <- print(vs_list[[1]][[1]],
                               nterms_max = 3,
                               stats = "mse"
    ))
    expect_equal(nrow(out$selection) - 1, 3)
    expect_named(out$selection, c(
      "size", "solution_terms",
      "mse", "se",
      "diff", "diff_se"
    ))

    expect_output(out <- print(cvs_list[[1]][[1]],
                               nterms_max = 3,
                               stats = "mse"
    ))
    expect_equal(nrow(out$selection) - 1, 3)
    expect_named(out$selection, c(
      "size", "solution_terms",
      paste0("mse.", tolower(out$cv_method)), "se",
      "diff", "diff_se"
      # "pct_solution_terms_cv"
    ))
  })


  # -------------------------------------------------------------
  context("plots")

  test_that("plotting works", {
    expect_s3_class(plot(vs_list[[1]][[1]]), "ggplot")
    expect_visible(plot(vs_list[[1]][[1]], nterms_max = 3))
  })

  test_that("invalid 'baseline' arguments are rejected", {
    expect_error(
      plot(vs_list[[1]][[1]], baseline = "zzz"),
      "Argument 'baseline' must be either 'ref' or 'best'"
    )
  })

  test_that("the value of nterms_max is valid", {
    expect_error(
      plot(vs_list[[1]][[1]], nterms_max = 0),
      "nterms_max must be at least 1"
    )
  })

  ## test_that("nterms_max is capped to the largest model size", {
  ##   expect_equal(
  ##     plot(vs_list[[1]][[1]]),
  ##     plot(vs_list[[1]][[1]], nterms_max = 1000)
  ##   )
  ## })


  # -------------------------------------------------------------
  context("suggest_size")

  test_that("suggest_size checks the length of stat", {
    expect_error(suggest_size(vs_list[[1]][["gauss"]], stat = valid_stats_all),
                 "Only one statistic")
  })

  test_that("suggest_size works on all stats", {
    for (stat in valid_stats_gauss) {
      suggested_size <- suggest_size(vs_list[[1]][["gauss"]], stat = stat)
      expect_true(!is.na(suggested_size))
      expect_true(suggested_size >= 0)
    }
    for (stat in valid_stats_binom) {
      suggested_size <- suggest_size(vs_list[[1]][["binom"]], stat = stat)
      expect_true(!is.na(suggested_size))
      expect_true(suggested_size >= 0)
    }
  })
}
