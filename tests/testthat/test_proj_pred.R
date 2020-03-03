context('proj_linpred')


# tests for proj_linpred and proj_predict

if (require(brms) && require(rstanarm)) {

  seed <- 1235
  set.seed(seed)
  n <- 40
  nv <- 5
  x <- matrix(rnorm(n*nv, 0, 1), n, nv)
  b <- runif(nv)-0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = T)
  offset <- rnorm(n)
  chains <- 2
  iter <- 500
  source(file.path('helpers', 'SW.R'))


  f_gauss <- gaussian()
  df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x%*%b), dis), x = x)
  f_binom <- binomial()
  df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x%*%b)), x = x, weights=weights)
  f_poiss <- poisson()
  df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x%*%b)), x = x)
  ys <- list()
  ys[[1]] <- df_gauss$y
  ys[[2]] <- df_binom$y/weights
  ys[[3]] <- df_poiss$y

  SW({
    fit_gauss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_gauss, data = df_gauss,
      chains = chains, seed = seed, iter = iter
    )
    fit_binom <- brm(y | trials(weights) ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_binom, data = df_binom,
      chains = chains, seed = seed, iter = iter
    )
    fit_poiss <- brm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
      family = f_poiss, data = df_poiss,
      chains = chains, seed = seed, iter = iter
    )
  })
  fit_list <- list(gauss = fit_gauss, binom = fit_binom, poiss = fit_poiss)
  vs_list <- lapply(fit_list, varsel, nv_max = nv + 1, verbose = FALSE)
  proj_vind_list <- lapply(vs_list, project, vind = c(2,3), seed = seed)
  proj_all_list <- lapply(vs_list, project, seed = seed, nv=0:nv)


  test_that("proj_linpred: xnew is specified correctly", {
    expect_error(proj_linpred(proj_vind_list),
                 'argument "xnew" is missing, with no default')
    expect_error(proj_linpred(proj_vind_list, xnew = NULL),
                 'must be a data.frame or a matrix')
    expect_error(proj_linpred(proj_vind_list, xnew = data.frame(x.1=x[, 1])),
                 'must be a data.frame or a matrix')
    expect_error(proj_linpred(proj_vind_list, xnew = data.frame(x=x), vind = 1:1000),
                 'number of columns in xnew does not match')
    expect_error(proj_linpred(proj_vind_list, xnew = data.frame(x=x[, 1:2])),
                 'xnew has 2 columns, but vind expects 3 columns')
  })

  test_that("output of proj_linpred is sensible with fit-object as input", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      pl <- proj_linpred(vs_list[[i]], xnew = data.frame(x=x), nv = 0:nv)
      expect_length(pl, nv + 1)
      for(j in 1:length(pl))
        expect_equal(ncol(pl[[j]]), n, info = i_inf)
    }
  })

  test_that("output of proj_linpred is sensible with project-object as input", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      pl <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x))
      expect_equal(ncol(pl), n, info = i_inf)
    }
    for(i in 1:length(proj_all_list)) {
      i_inf <- names(proj_all_list)[i]
      pl <- proj_linpred(proj_all_list[[i]], xnew = data.frame(x=x))
      expect_length(pl, nv + 1)
      for(j in 1:length(pl))
        expect_equal(ncol(pl[[j]]), n, info = i_inf)
    }
  })

  test_that("proj_linpred: error when varsel has not been performed on the object", {
    expect_error(proj_linpred(1, xnew = data.frame(x=x)),
                 'is not a variable selection object')
    expect_error(proj_linpred(fit_gauss, xnew = data.frame(x=x)),
                 'is not a variable selection object')
    expect_error(proj_linpred(c(proj_vind_list, list(x)), xnew = x),
                 'contains objects not created by varsel')
  })

  test_that("proj_linpred: specifying ynew incorrectly produces an error", {
    expect_error(proj_linpred(vs_list[["gauss"]], xnew = data.frame(x=x), ynew = x[, 1:3]),
                 'y cannot have more than two columns')
    expect_error(proj_linpred(vs_list[["gauss"]], xnew = data.frame(x=x), ynew = factor(ys[[1]])),
                 'cannot be a factor')
    expect_error(proj_linpred(vs_list[["poiss"]], xnew = data.frame(x=x), ynew = factor(ys[[3]])),
                 'cannot be a factor')
    expect_error(proj_linpred(vs_list[["binom"]], xnew = data.frame(x=x), ynew = ys[[1]]),
                 'y values must be 0 <= y <= 1 for the binomial model')
    expect_error(proj_linpred(vs_list[["binom"]], xnew = data.frame(x=x), ynew = factor(ys[[1]])),
                 'y cannot contain more than two classes')
  })

  test_that("proj_linpred: specifying ynew has an expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      pl <- proj_linpred(vs_list[[i]], xnew = data.frame(x=x), ynew = ys[[i]], weightsnew=weights, nv = 0:nv)
      pl2 <- proj_linpred(vs_list[[i]], xnew = data.frame(x=x), weightsnew=weights, nv = 0:nv)
      for(j in 1:length(pl)) {
        expect_named(pl[[j]], c('pred', 'lpd'))
        expect_equal(ncol(pl[[j]]$pred), n, info = i_inf)
        expect_equal(ncol(pl[[j]]$lpd), n, info = i_inf)
        expect_equal(ncol(pl2[[j]]), n, info = i_inf)
      }
    }
  })

  test_that("proj_linpred: specifying ynew as a factor works in a binomial model", {
    yfactor <- factor(rbinom(n, 1, 0.5))
    pl <- proj_linpred(vs_list[["binom"]], xnew = data.frame(x=x), ynew = yfactor)
    expect_named(pl, c('pred', 'lpd'))
    expect_equal(ncol(pl$pred), n)
    expect_equal(ncol(pl$lpd), n)
  })

  test_that("proj_linpred: specifying weights has an expected effect", {
    for(i in 1:length(proj_vind_list)) {
      # for binomial models weights have to be specified
      if (proj_vind_list[[i]]$family$family != 'binomial') {
        i_inf <- names(proj_vind_list)[i]
        plw <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x), ynew = ys[[i]],
                            weightsnew = weights)
        pl <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x), ynew = ys[[i]])
        expect_named(plw, c('pred', 'lpd'))
        expect_equal(ncol(plw$pred), n, info = i_inf)
        expect_equal(ncol(plw$lpd), n, info = i_inf)
        expect_true(sum(plw$lpd != pl$lpd) > 0, info = i_inf)
      }
    }
  })

  test_that("proj_linpred: specifying offset has an expected effect", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      plo <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x), ynew = ys[[i]], weightsnew=weights,
                          offsetnew = offset)
      pl <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x), ynew = ys[[i]], weightsnew=weights)
      expect_named(plo, c('pred', 'lpd'))
      expect_equal(ncol(plo$pred), n, info = i_inf)
      expect_equal(ncol(plo$lpd), n, info = i_inf)
      expect_true(sum(plo$lpd != pl$lpd) > 0, info = i_inf)
    }
  })

  test_that("proj_linpred: specifying transform has an expected effect", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      plt <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x), transform = TRUE)
      plf <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x), transform = FALSE)
      expect_equal(proj_vind_list[[i]]$family$linkinv(plf), plt, info = i_inf)
    }
  })

  test_that("proj_linpred: specifying integrated has an expected effect", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      plt <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x), integrated = TRUE)
      plf <- proj_linpred(proj_vind_list[[i]], xnew = data.frame(x=x), integrated = FALSE)
      expect_equal(drop(proj_vind_list[[i]]$weights%*%plf), plt, info = i_inf)
    }
  })

  test_that("proj_linpred: adding more regularization has an expected effect", {
    regul <- c(1e-6, 1e-1, 1e2)
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      norms <- rep(0, length(regul))
      for (j in 1:length(regul)) {
        pred <- proj_linpred(vs_list[[i]], xnew = data.frame(x=x), nv = 2, transform = FALSE,
                             integrated = TRUE, regul=regul[j])
        norms[j] <- sum(pred^2)
      }
      for (j in 1:(length(regul)-1))
        expect_true(all(norms[j] >= norms[j+1]), info = i_inf)
    }
  })


  test_that("proj_linpred: arguments passed to project work accordingly", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      pr <- project(vs_list[[i]], nv = c(2, 4), nc = 2, ns = 20,
                    intercept = FALSE, regul = 1e-8, seed = 12)
      prl1 <- proj_linpred(pr, xnew = data.frame(x=x))
      prl2 <- proj_linpred(vs_list[[i]], xnew = data.frame(x=x), nv = c(2, 4), nc = 2, ns = 20,
                           intercept = FALSE, regul = 1e-8, seed = 12)
      expect_equal(prl1, prl2, info = i_inf)
    }
  })

  test_that("proj_linpred: providing xnew as a data frame works as expected", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      pl <- proj_predict(proj_vind_list[[i]],
                         xnew = data.frame(x))
      expect_equal(ncol(pl), n, info = i_inf)
    }
    SW(
      fit_form <- stan_glm(mpg~(drat + wt)^2, data = mtcars, QR = T,
                           chains = chains, seed = seed, iter = iter)
    )
    vs_form <- varsel(fit_form)
    p1 <- proj_linpred(vs_form, xnew = mtcars, nv = 3, seed = 2)
    p2 <- proj_linpred(vs_form, xnew = get_x(fit_form)[,-1], nv = 3, seed = 2)
    expect_equal(p1, p2)
  })


  # -------------------------------------------------------------
  context('proj_predict')

  test_that("proj_predict: xnew is specified correctly", {
    expect_error(proj_predict(proj_vind_list),
                 'argument "xnew" is missing, with no default')
    expect_error(proj_predict(proj_vind_list, xnew = NULL),
                 'must be a data.frame or a matrix')
    expect_error(proj_predict(proj_vind_list, xnew = data.frame(x.1=x[, 1])),
                 'must be a data.frame or a matrix')
    expect_error(proj_predict(proj_vind_list, xnew = data.frame(x=x), vind = 1:1000),
                 'number of columns in xnew does not match')
    expect_error(proj_predict(proj_vind_list, xnew = data.frame(x=x[, 1:2])),
                 'xnew has 2 columns, but vind expects 3 columns')
  })

  test_that("output of proj_predict is sensible with fit-object as input", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      pl <- proj_predict(vs_list[[i]], xnew = data.frame(x=x), nv = 0:nv)
      expect_length(pl, nv + 1)
      for(j in 1:length(pl))
        expect_equal(ncol(pl[[j]]), n, info = i_inf)
    }
  })

  test_that("output of proj_predict is sensible with project-object as input", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      pl <- proj_predict(proj_vind_list[[i]], xnew = data.frame(x=x))
      expect_equal(ncol(pl), n, info = i_inf)
    }
    for(i in 1:length(proj_all_list)) {
      i_inf <- names(proj_all_list)[i]
      pl <- proj_predict(proj_all_list[[i]], xnew = data.frame(x=x))
      expect_length(pl, nv + 1)
      for(j in 1:length(pl))
        expect_equal(ncol(pl[[j]]), n, info = i_inf)
    }
  })

  test_that("proj_predict: error when varsel has not been performed on the object", {
    expect_error(proj_predict(1, xnew = data.frame(x=x)),
                 'is not a variable selection object')
    expect_error(proj_predict(fit_gauss, xnew = data.frame(x=x)),
                 'is not a variable selection object')
    expect_error(proj_predict(c(proj_vind_list, list(x)), xnew = data.frame(x=x)),
                 'contains objects not created by varsel')
  })

  test_that("proj_predict: specifying ynew has an expected effect", {
    for (i in 1:length(vs_list)) {
      pl <- proj_predict(vs_list[[i]], xnew = data.frame(x=x), ynew = ys[[i]], nv = 0:3)
      pl2 <- proj_predict(vs_list[[i]], xnew = data.frame(x=x), nv = 0:3)
      for (j in 1:length(pl)) {
        expect_equal(dim(pl[[j]]), dim(pl2[[j]]))
      }
    }
  })

  test_that("proj_predict: specifying ynew as a factor works in a binomial model", {
    yfactor <- factor(rbinom(n, 1, 0.5))
    pl <- proj_predict(vs_list[["binom"]], xnew = data.frame(x=x), ynew = yfactor)
    expect_equal(ncol(pl), n)
    expect_true(all(pl %in% c(0, 1)))
  })

  test_that("proj_predict: specifying weightsnew has an expected effect", {
    pl <- proj_predict(proj_vind_list[['binom']], xnew = data.frame(x=x), seed = seed)
    plw <- proj_predict(proj_vind_list[['binom']], xnew = data.frame(x=x), seed = seed,
                        weightsnew = weights)
    expect_true(sum(pl != plw)>0)
  })

  test_that("proj_predict: specifying offsetnew has an expected effect", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      pl <- proj_predict(proj_vind_list[[i]], xnew = data.frame(x=x), draws = iter,
                         seed = seed)
      plo <- proj_predict(proj_vind_list[[i]], xnew = data.frame(x=x), draws = iter,
                          seed = seed, offsetnew = offset)
      expect_true(sum(pl != plo) > 0, info = i_inf)
    }
  })

  test_that("proj_predict: specifying draws has an expected effect", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      pl <- proj_predict(proj_vind_list[[i]], xnew = data.frame(x=x), draws = iter)
      expect_equal(dim(pl), c(iter, n))
    }
  })

  test_that("proj_predict: specifying seed_sam has an expected effect", {
    for(i in 1:length(proj_vind_list)) {
      i_inf <- names(proj_vind_list)[i]
      pl1 <- proj_predict(proj_vind_list[[i]], xnew = data.frame(x=x), seed = seed)
      pl2 <- proj_predict(proj_vind_list[[i]], xnew = data.frame(x=x), seed = seed)
      expect_equal(pl1, pl2, info = i_inf)
    }
  })

  test_that("proj_predict: arguments passed to project work accordingly", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      pr1 <- project(vs_list[[i]], nv = c(2, 4), nc = 2, ns = 20,
                     intercept = FALSE, regul = 1e-8, seed = 12)
      prp1 <- proj_predict(pr1, xnew = data.frame(x=x), draws = 100, seed = 11)
      prp2 <- proj_predict(vs_list[[i]], xnew = data.frame(x=x), draws = 100, seed = 11,
                           nv = c(2, 4), nc = 2, ns = 20, intercept = FALSE,
                           regul = 1e-8, seed = 12)
      expect_equal(prp1, prp2, info = i_inf)
    }
  })

}
