context("proj_linpred")


# tests for proj_linpred and proj_predict

if (require(rstanarm) && Sys.getenv("NOT_CRAN") == "true") {
  seed <- 1235
  set.seed(seed)
  n <- 40
  nterms <- 5
  x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
  b <- runif(nterms) - 0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = TRUE)
  offset <- rnorm(n)
  chains <- 2
  iter <- 500
  source(file.path("helpers", "SW.R"))

  f_gauss <- gaussian()
  df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x %*% b), dis), x = x)
  f_binom <- binomial()
  df_binom <- data.frame(
    y = rbinom(n, weights, f_binom$linkinv(x %*% b)), x = x,
    weights = weights
  )
  f_poiss <- poisson()
  df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x %*% b)), x = x)
  ys <- list()
  ys[[1]] <- df_gauss$y
  ys[[2]] <- df_binom$y
  ys[[3]] <- df_poiss$y

  SW({
    fit_gauss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_gauss, data = df_gauss,
      chains = chains, seed = seed, iter = iter
    )
    fit_binom <- stan_glm(cbind(y, weights - y) ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_binom, data = df_binom, weights = weights,
      chains = chains, seed = seed, iter = iter
    )
    fit_poiss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_poiss, data = df_poiss,
      chains = chains, seed = seed, iter = iter
    )
    fit_list <- list(
      gauss = fit_gauss,
      binom = fit_binom,
      poiss = fit_poiss
    )
    vs_list <- lapply(fit_list, varsel,
      nterms_max = nterms + 1,
      verbose = FALSE
    )
    proj_solution_terms_list <- lapply(vs_list, project,
      solution_terms = c(2, 3),
      seed = seed
    )
    proj_all_list <- lapply(vs_list, project,
      seed = seed,
      nterms = 0:nterms
    )
  })

  test_that("proj_linpred: newdata is specified correctly", {
    ## expect_error(
    ##   proj_linpred(proj_solution_terms_list),
    ##   'argument "newdata" is missing, with no default'
    ## )
    ## expect_error(
    ##   proj_linpred(proj_solution_terms_list, newdata = NULL),
    ##   "must be a data.frame or a matrix"
    ## )
    expect_error(
      proj_linpred(proj_solution_terms_list, newdata = x[, 1]),
      "must be a data.frame or a matrix"
    )
    expect_error(
      proj_linpred(proj_solution_terms_list, newdata = data.frame(x = x),
                   solution_terms = 1:10000),
      "number of columns in newdata does not match"
    )
    expect_error(
      proj_linpred(proj_solution_terms_list, newdata = data.frame(x = x)[, 1:2],
                   solution_terms = 1:3),
      "number of columns in newdata does not match"
    )
  })

  test_that("output of proj_linpred is sensible with fit-object as input", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      y <- vs_list[[i]]$refmodel$y
      pl <- proj_linpred(vs_list[[i]], newdata = data.frame(y = y, x = x),
                         nterms = 0:nterms)
      expect_length(pl, nterms + 1)
    }
  })

  test_that("output of proj_linpred is sensible with project-object as input", {
    for (i in 1:length(proj_solution_terms_list)) {
      i_inf <- names(proj_solution_terms_list)[i]
      y <- proj_solution_terms_list[[i]]$refmodel$y
      pl <- proj_linpred(proj_solution_terms_list[[i]],
        newdata = data.frame(y = y, x = x)
      )
    }
    for (i in 1:length(proj_all_list)) {
      i_inf <- names(proj_all_list)[i]
      y <- proj_all_list[[i]][[1]]$refmodel$y
      pl <- proj_linpred(proj_all_list[[i]],
        newdata = data.frame(y = y, x = x)
      )
      expect_length(pl, nterms + 1)
    }
  })

  test_that(paste("proj_linpred: error when varsel has not been performed on",
                  "the object (and 'solution_terms' is provided neither)"), {
    expect_error(
      proj_linpred(1, newdata = data.frame(x = x)),
      "is not a variable selection -object"
    )
    expect_error(
      proj_linpred(fit_gauss, newdata = data.frame(x = x)),
      "is not a variable selection -object"
    )
    expect_error(
      proj_linpred(c(proj_solution_terms_list, list(x)), newdata = x),
      "Invalid object supplied to argument 'object'\\."
    )
  })

  ## test_that("proj_linpred: specifying ynew incorrectly produces an error", {
  ##   expect_error(
  ##     proj_linpred(vs_list[["gauss"]], newdata = data.frame(x = x),
  ##                  ynew = x[, 1:3]),
  ##     "y cannot have more than two columns"
  ##   )
  ##   expect_error(
  ##     proj_linpred(vs_list[["gauss"]], newdata = data.frame(x = x),
  ##                  ynew = factor(ys[[1]])),
  ##     "cannot be a factor"
  ##   )
  ##   expect_error(
  ##     proj_linpred(vs_list[["poiss"]], newdata = data.frame(x = x),
  ##                  ynew = factor(ys[[3]])),
  ##     "cannot be a factor"
  ##   )
  ##   expect_error(
  ##     proj_linpred(vs_list[["binom"]], newdata = data.frame(x = x),
  ##                  ynew = factor(ys[[1]])),
  ##     "y cannot contain more than two classes"
  ##   )
  ## })

  ## test_that("proj_linpred: specifying ynew has an expected effect", {
  ##   for (i in 1:length(vs_list)) {
  ##     i_inf <- names(vs_list)[i]
  ##     pl <- proj_linpred(vs_list[[i]],
  ##       newdata = df_binom, ynew = ys[[i]],
  ##       weightsnew = ~weights, nterms = 0:nterms
  ##     )
  ##     pl2 <- proj_linpred(vs_list[[i]],
  ##       newdata = data.frame(x = x, weights = weights),
  ##       weightsnew = ~weights, nterms = 0:nterms
  ##     )
  ##     for (j in 1:length(pl)) {
  ##       expect_named(pl[[j]], c("pred", "lpd"))
  ##       expect_equal(ncol(pl[[j]]$pred), n, info = i_inf)
  ##       expect_equal(nrow(pl[[j]]$lpd), n, info = i_inf)
  ##     }
  ##   }
  ## })

  ## test_that(paste("proj_linpred: specifying ynew as a factor works in a",
  ##                 "binomial model"), {
  ##   yfactor <- factor(rbinom(n, 1, 0.5))
  ##   pl <- proj_linpred(vs_list[["binom"]], newdata = data.frame(x = x),
  ##                      ynew = yfactor)
  ##   expect_named(pl, c("pred", "lpd"))
  ##   expect_equal(ncol(pl$pred), n)
  ##   expect_equal(nrow(pl$lpd), n)
  ## })

  test_that("proj_linpred: specifying weights has an expected effect", {
    for (i in 1:length(proj_solution_terms_list)) {
      # for binomial models weights have to be specified
      if (proj_solution_terms_list[[i]]$family$family != "binomial") {
        i_inf <- names(proj_solution_terms_list)[i]
        weightsnew <- sample(1:4, n, replace = TRUE)
        plw <- proj_linpred(proj_solution_terms_list[[i]],
          newdata = data.frame(y = ys[[i]], x = x, weights = weightsnew),
          weightsnew = ~weights
        )
        pl <- proj_linpred(proj_solution_terms_list[[i]],
          newdata = data.frame(y = ys[[i]], x = x, weights = weights),
          weightsnew = ~weights
        )
        expect_named(plw, c("pred", "lpd"))
        expect_equal(ncol(plw$pred), n, info = i_inf)
        expect_equal(nrow(plw$lpd), n, info = i_inf)
        expect_false(all(plw$lpd == pl$lpd))
      }
    }
  })

  test_that("proj_linpred: specifying offset has an expected effect", {
    for (i in 1:length(proj_solution_terms_list)) {
      i_inf <- names(proj_solution_terms_list)[i]
      plo <- proj_linpred(proj_solution_terms_list[[i]],
        newdata = data.frame(
          y = ys[[i]], x = x, weights = weights,
          offset = offset
        ),
        weightsnew = ~weights, offsetnew = ~offset
      )
      pl <- proj_linpred(proj_solution_terms_list[[i]],
        newdata = data.frame(y = ys[[i]], x = x, weights = weights),
        weightsnew = ~weights
      )
      expect_named(plo, c("pred", "lpd"))
      expect_equal(ncol(plo$pred), n, info = i_inf)
      expect_equal(nrow(plo$lpd), n, info = i_inf)
      expect_equal(t(plo$pred) - offset, t(pl$pred), tol = 1e-8)
    }
  })

  test_that("proj_linpred: specifying transform has an expected effect", {
    for (i in 1:length(proj_solution_terms_list)) {
      i_inf <- names(proj_solution_terms_list)[i]
      y <- proj_solution_terms_list[[i]]$refmodel$y
      plt <- proj_linpred(proj_solution_terms_list[[i]],
                          newdata = data.frame(y = y, x = x), transform = TRUE)
      plf <- proj_linpred(proj_solution_terms_list[[i]],
                          newdata = data.frame(y = y, x = x), transform = FALSE)
      expect_equal(proj_solution_terms_list[[i]]$family$linkinv(plf$pred),
                   plt$pred, info = i_inf)
    }
  })

  test_that("proj_linpred: specifying integrated has an expected effect", {
    for (i in 1:length(proj_solution_terms_list)) {
      i_inf <- names(proj_solution_terms_list)[i]
      y <- proj_solution_terms_list[[i]]$refmodel$y
      plt <- proj_linpred(proj_solution_terms_list[[i]],
                          newdata = data.frame(y = y, x = x), integrated = TRUE)
      plf <- proj_linpred(proj_solution_terms_list[[i]],
                          newdata = data.frame(y = y, x = x), integrated = FALSE)
      expect_equal(as.vector(proj_solution_terms_list[[i]]$weights %*%
                             plf$pred),
                   plt$pred, info = i_inf)
    }
  })

  test_that("proj_linpred: adding more regularization has an expected effect", {
    regul <- c(1e-6, 1e-1, 1e2)
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      norms <- rep(0, length(regul))
      for (j in 1:length(regul)) {
        y <- vs_list[[i]]$refmodel$y
        pred <- proj_linpred(vs_list[[i]],
          newdata = data.frame(y = y, x = x), nterms = 2, transform = FALSE,
          integrated = TRUE, regul = regul[j]
        )
        norms[j] <- sum(pred$pred^2)
      }
      for (j in 1:(length(regul) - 1)) {
        expect_true(all(norms[j] >= norms[j + 1]), info = i_inf)
      }
    }
  })


  test_that("proj_linpred: arguments passed to project work accordingly", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      y <- vs_list[[i]]$refmodel$y
      SW(pr <- project(vs_list[[i]],
        nterms = c(2, 4), nclusters = 2, ndraws = 20,
        intercept = FALSE, regul = 1e-8, seed = 12
      ))
      prl1 <- proj_linpred(pr, newdata = data.frame(y = y, x = x))
      SW(prl2 <- proj_linpred(vs_list[[i]],
        newdata = data.frame(y = y, x = x), nterms = c(2, 4), nclusters = 2,
        ndraws = 20, intercept = FALSE, regul = 1e-8, seed = 12
      ))
      expect_equal(prl1$pred, prl2$pred, info = i_inf)
    }
  })

  test_that("proj_linpred: providing newdata as a data frame works as expected",
  {
    for (i in 1:length(proj_solution_terms_list)) {
      i_inf <- names(proj_solution_terms_list)[i]
      y <- proj_solution_terms_list[[i]]$refmodel$y
      pl <- proj_predict(proj_solution_terms_list[[i]],
        newdata = data.frame(y = y, x = x)
      )
      expect_equal(ncol(pl), n, info = i_inf)
    }
    SW(
      fit_form <- stan_glm(mpg ~ (drat + wt)^2,
        data = mtcars, QR = TRUE,
        chains = chains, seed = seed, iter = iter
      )
    )
    vs_form <- varsel(fit_form)
    p1 <- proj_linpred(vs_form, newdata = mtcars, nterms = 3, seed = 2)
    x <- get_x(fit_form)[, -1]
    newdata <- data.frame(mpg = get_y(fit_form), x)
    p2 <- proj_linpred(vs_form,
      newdata = newdata, nterms = 3,
      seed = 2
    )
    expect_equal(p1$pred, p2$pred)
  })


  # -------------------------------------------------------------
  context("proj_predict")

  test_that("proj_predict: newdata is specified correctly", {
    ## expect_error(
    ##   proj_predict(proj_solution_terms_list),
    ##   'argument "newdata" is missing, with no default'
    ## )
    ## expect_error(
    ##   proj_predict(proj_solution_terms_list, newdata = NULL),
    ##   "must be a data.frame or a matrix"
    ## )
    expect_error(
      proj_predict(proj_solution_terms_list, newdata = x[, 1]),
      "must be a data.frame or a matrix"
    )
    expect_error(
      proj_predict(proj_solution_terms_list, newdata = data.frame(x = x),
                   solution_terms = 1:1000),
      "number of columns in newdata does not match"
    )
    expect_error(
      proj_predict(proj_solution_terms_list,
        newdata = data.frame(x = x)[, 1:2],
        solution_terms = 1:3
      ),
      "number of columns in newdata does not match"
    )
  })

  test_that("output of proj_predict is sensible with fit-object as input", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      pl <- proj_predict(vs_list[[i]], newdata = data.frame(x = x),
                         nterms = 0:nterms)
      expect_length(pl, nterms + 1)
      for (j in 1:length(pl)) {
        expect_equal(ncol(pl[[j]]), n, info = i_inf)
      }
    }
  })

  test_that("output of proj_predict is sensible with project-object as input", {
    for (i in 1:length(proj_solution_terms_list)) {
      i_inf <- names(proj_solution_terms_list)[i]
      y <- proj_solution_terms_list[[i]]$refmodel$y
      pl <- proj_predict(proj_solution_terms_list[[i]],
                         newdata = data.frame(y = y, x = x))
      expect_equal(ncol(pl), n, info = i_inf)
    }
    for (i in 1:length(proj_all_list)) {
      i_inf <- names(proj_all_list)[i]
      pl <- proj_predict(proj_all_list[[i]], newdata = data.frame(x = x))
      expect_length(pl, nterms + 1)
      for (j in 1:length(pl)) {
        expect_equal(ncol(pl[[j]]), n, info = i_inf)
      }
    }
  })

  test_that(paste("proj_predict: error when varsel has not been performed on",
                  "the object (and 'solution_terms' is provided neither)"), {
    expect_error(
      proj_predict(1, newdata = data.frame(x = x)),
      "is not a variable selection -object"
    )
    expect_error(
      proj_predict(fit_gauss, newdata = data.frame(x = x)),
      "is not a variable selection -object"
    )
    expect_error(
      proj_predict(c(proj_solution_terms_list, list(x)),
                   newdata = data.frame(x = x)),
      "Invalid object supplied to argument 'object'\\."
    )
  })

  ## test_that("proj_predict: specifying ynew has an expected effect", {
  ##   for (i in seq_along(vs_list)) {
  ##     pl <- proj_predict(vs_list[[i]], newdata = data.frame(x = x),
  ##                        ynew = ys[[i]], nterms = 0:3)
  ##     pl2 <- proj_predict(vs_list[[i]], newdata = data.frame(x = x),
  ##                         nterms = 0:3)
  ##     for (j in seq_len(length(pl))) {
  ##       expect_equal(dim(pl[[j]]), dim(pl2[[j]]))
  ##     }
  ##   }
  ## })

  ## test_that(paste("proj_predict: specifying ynew as a factor works in a",
  ##                 "binomial model"), {
  ##   yfactor <- factor(rbinom(n, 1, 0.5))
  ##   pl <- proj_predict(vs_list[["binom"]], newdata = data.frame(x = x),
  ##                      ynew = yfactor)
  ##   expect_equal(ncol(pl), n)
  ##   expect_true(all(pl %in% c(0, 1)))
  ## })

  test_that("proj_predict: specifying weightsnew has an expected effect", {
    pl <- proj_predict(proj_solution_terms_list[["binom"]],
      newdata = data.frame(x = x, weights = rep(1, NROW(x))),
      seed = seed
    )
    plw <- proj_predict(proj_solution_terms_list[["binom"]],
      newdata = data.frame(x = x, weights = weights), seed = seed,
      weightsnew = ~weights
    )
    expect_true(sum(pl != plw) > 0)
  })

  test_that("proj_predict: specifying offsetnew has an expected effect", {
    for (i in seq_len(length(proj_solution_terms_list))) {
      i_inf <- names(proj_solution_terms_list)[i]
      pl <- proj_predict(proj_solution_terms_list[[i]],
        newdata = data.frame(x = x), ndraws = iter,
        seed = seed
      )
      plo <- proj_predict(proj_solution_terms_list[[i]],
        newdata = data.frame(x = x, offset = offset), ndraws = iter,
        seed = seed, offsetnew = ~offset
      )
      expect_true(sum(pl != plo) > 0, info = i_inf)
    }
  })

  test_that("proj_predict: specifying ndraws has an expected effect", {
    for (i in 1:length(proj_solution_terms_list)) {
      i_inf <- names(proj_solution_terms_list)[i]
      pl <- proj_predict(proj_solution_terms_list[[i]],
                         newdata = data.frame(x = x), ndraws = iter)
      expect_equal(dim(pl), c(iter, n))
    }
  })

  test_that("proj_predict: specifying seed_sam has an expected effect", {
    for (i in 1:length(proj_solution_terms_list)) {
      i_inf <- names(proj_solution_terms_list)[i]
      pl1 <- proj_predict(proj_solution_terms_list[[i]],
                          newdata = data.frame(x = x), seed = seed)
      pl2 <- proj_predict(proj_solution_terms_list[[i]],
                          newdata = data.frame(x = x), seed = seed)
      expect_equal(pl1, pl2, info = i_inf)
    }
  })

  test_that("proj_predict: arguments passed to project work accordingly", {
    for (i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      prp1 <- proj_predict(vs_list[[i]],
        newdata = data.frame(x = x), ndraws = 100,
        seed = 12, nterms = c(2, 4), nclusters = 2,
        regul = 1e-08
      )
      prp2 <- proj_predict(vs_list[[i]],
        newdata = data.frame(x = x), ndraws = 100,
        nterms = c(2, 4), nclusters = 2, regul = 1e-8,
        seed = 12
      )
      prp3 <- proj_predict(vs_list[[i]],
        newdata = data.frame(x = x), ndraws = 100,
        seed = 120, nterms = c(2, 4), nclusters = 2,
        regul = 1e-08
      )
      expect_equal(prp1, prp2, info = i_inf)
      expect_false(all(unlist(lapply(seq_along(prp1), function(i) {
        all(prp1[[i]] == prp3[[i]])
      }))),
      info = i_inf
      )
    }
  })
}
