context("proj_linpred")

test_that("proj_linpred: `newdata` is checked correctly", {
  expect_error(
    proj_linpred(prjs_solterms, newdata = dat[, 1]),
    "must be a data.frame or a matrix"
  )
  expect_error(
    proj_linpred(prjs_solterms,
                 solution_terms = paste0("x.", 1:10000)),
    paste("^The number of solution terms is greater than the number of",
          "columns in newdata\\.$")
  )
  stopifnot(length(solterms_x) > 1)
  expect_error(
    proj_linpred(prjs_solterms$glm.gauss.solterms_x.clust,
                 newdata = dat[, 1, drop = FALSE],
                 solution_terms = solterms_x),
    paste("^The number of solution terms is greater than the number of",
          "columns in newdata\\.$")
  )
})

test_that(paste(
  "proj_linpred: \"refmodel\" object as input leads to correct output",
  "structure"
), {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      tstsetup <- unlist(nlist(mod_nm, fam_nm))
      pl <- proj_linpred(refmods[[mod_nm]][[fam_nm]],
                         solution_terms = solterms_x,
                         nclusters = nclusters_pred_tst)
      expect_named(pl, c("pred", "lpd"), info = tstsetup)
      expect_identical(dim(pl$pred), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
      expect_identical(dim(pl$lpd), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
    }
  }
})

test_that(paste(
  "proj_linpred: \"vsel\" object as input leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    y <- vs_list[[i]]$refmodel$y
    pl <- proj_linpred(vs_list[[i]], nclusters = nclusters_pred_tst,
                       newdata = data.frame(y = y, x = x),
                       nterms = 0:nterms)
    expect_length(pl, nterms + 1)
    for (j in seq_along(pl)) {
      expect_named(pl[[!!j]], c("pred", "lpd"), info = i)
      expect_identical(dim(pl[[!!j]]$pred), c(nclusters_pred_tst, n), info = i)
      expect_identical(dim(pl[[!!j]]$lpd), c(nclusters_pred_tst, n), info = i)
    }
  }
})

test_that(paste(
  "proj_linpred: \"projection\" object as input leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    y <- prjs_solterms[[i]]$refmodel$y
    pl <- proj_linpred(prjs_solterms[[i]],
                       newdata = data.frame(y = y, x = x))
    expect_named(pl, c("pred", "lpd"), info = i)
    expect_identical(dim(pl$pred), c(nclusters_pred_tst, n), info = i)
    expect_identical(dim(pl$lpd), c(nclusters_pred_tst, n), info = i)
  }
})

test_that(paste(
  "proj_linpred: \"proj_list\" object (an informal class) as input leads to",
  "correct output structure"
), {
  for (i in fam_nms) {
    y <- proj_all_list[[i]][[1]]$refmodel$y
    pl <- proj_linpred(proj_all_list[[i]],
                       newdata = data.frame(y = y, x = x))
    expect_length(pl, nterms + 1)
    for (j in seq_along(pl)) {
      expect_named(pl[[!!j]], c("pred", "lpd"), info = i)
      expect_identical(dim(pl[[!!j]]$pred), c(nclusters_pred_tst, n), info = i)
      expect_identical(dim(pl[[!!j]]$lpd), c(nclusters_pred_tst, n), info = i)
    }
  }
})

test_that("proj_linpred: output structure is also correct in edge cases", {
  for (i in fam_nms) {
    y <- refmod_list[[i]]$y
    for (n_crr in c(1L, 12L)) {
      for (nclusters_pred_crr in c(1L, 4L)) {
        for (integrated_crr in c(FALSE, TRUE)) {
          pl <- proj_linpred(
            refmod_list[[i]], nclusters = nclusters_pred_crr,
            newdata = head(data.frame(y = y, x = x), n_crr),
            integrated = integrated_crr,
            solution_terms = c("x.3", "x.5")
          )
          tstsetup <- unlist(nlist(i, n_crr, nclusters_pred_crr,
                                   integrated_crr))
          expect_named(pl, c("pred", "lpd"), info = tstsetup)
          nprjdraws_crr <- ifelse(integrated_crr,
                                  1L, nclusters_pred_crr)
          expect_identical(dim(pl$pred), c(nprjdraws_crr, n_crr),
                           info = tstsetup)
          expect_identical(dim(pl$lpd), c(nprjdraws_crr, n_crr),
                           info = tstsetup)
        }
      }
    }
  }
})

test_that(paste(
  "proj_linpred: error when varsel has not been performed on",
  "the object (and `solution_terms` is provided neither)"
), {
  expect_error(
    proj_linpred(1, newdata = data.frame(x = x)),
    "is not an object of class \"vsel\""
  )
  expect_error(
    proj_linpred(fit_gauss, newdata = data.frame(x = x)),
    "is not an object of class \"vsel\""
  )
  expect_error(
    proj_linpred(c(prjs_solterms, list(x)), newdata = x),
    "Invalid object supplied to argument `object`\\."
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
##   for (i in fam_nms) {
##     pl <- proj_linpred(vs_list[[i]], nclusters = nclusters_pred_tst,
##       newdata = df_binom, ynew = ys[[i]],
##       weightsnew = ~weights, nterms = 0:nterms
##     )
##     pl2 <- proj_linpred(vs_list[[i]], nclusters = nclusters_pred_tst,
##       newdata = data.frame(x = x, weights = weights),
##       weightsnew = ~weights, nterms = 0:nterms
##     )
##     for (j in 1:length(pl)) {
##       expect_named(pl[[j]], c("pred", "lpd"))
##       expect_equal(ncol(pl[[!!j]]$pred), n, info = i)
##       expect_equal(nrow(pl[[!!j]]$lpd), n, info = i)
##     }
##   }
## })

## test_that(paste(
##   "proj_linpred: specifying ynew as a factor works in a",
##   "binomial model"
## ), {
##   yfactor <- factor(rbinom(n, 1, 0.5))
##   pl <- proj_linpred(vs_list[["binom"]], nclusters = nclusters_pred_tst,
##                      newdata = data.frame(x = x),
##                      ynew = yfactor)
##   expect_named(pl, c("pred", "lpd"))
##   expect_equal(ncol(pl$pred), n)
##   expect_equal(nrow(pl$lpd), n)
## })

test_that(paste(
  "proj_linpred: omitting the response causes output element `lpd` to be",
  "`NULL`."
), {
  stopifnot(!exists("y"))
  for (i in fam_nms) {
    i_resampled <- sample.int(nrow(x))
    stopifnot(identical(sort(i_resampled), seq_len(nrow(x))))
    pl <- proj_linpred(prjs_solterms[[i]],
                       newdata = data.frame(
                         x = x[i_resampled, , drop = FALSE]
                       ))
    expect_named(pl, c("pred", "lpd"), info = i)
    expect_identical(dim(pl$pred), c(nclusters_pred_tst, n), info = i)
    expect_null(pl$lpd, info = i)
  }
})

test_that("proj_linpred: specifying weights has an expected effect", {
  for (i in fam_nms) {
    # for binomial models weights have to be specified
    if (prjs_solterms[[i]]$family$family != "binomial") {
      weightsnew <- sample(1:4, n, replace = TRUE)
      plw <- proj_linpred(prjs_solterms[[i]],
                          newdata = data.frame(y = ys[[i]], x = x,
                                               weights = weightsnew),
                          weightsnew = ~weights)
      pl <- proj_linpred(prjs_solterms[[i]],
                         newdata = data.frame(y = ys[[i]], x = x,
                                              weights = weights),
                         weightsnew = ~weights)
      expect_named(plw, c("pred", "lpd"))
      expect_equal(ncol(plw$pred), n, info = i)
      expect_equal(ncol(plw$lpd), n, info = i)
      expect_false(all(plw$lpd == pl$lpd))
    }
  }
})

test_that("proj_linpred: specifying offset has an expected effect", {
  for (i in fam_nms) {
    plo <- proj_linpred(prjs_solterms[[i]],
                        newdata = data.frame(
                          y = ys[[i]], x = x, weights = weights,
                          offset = offset
                        ),
                        weightsnew = ~weights, offsetnew = ~offset)
    pl <- proj_linpred(prjs_solterms[[i]],
                       newdata = data.frame(y = ys[[i]], x = x,
                                            weights = weights),
                       weightsnew = ~weights)
    expect_named(plo, c("pred", "lpd"))
    expect_equal(ncol(plo$pred), n, info = i)
    expect_equal(ncol(plo$lpd), n, info = i)
    expect_equal(t(plo$pred) - offset, t(pl$pred), tol = 1e-8)
  }
})

test_that("proj_linpred: specifying transform has an expected effect", {
  for (i in fam_nms) {
    y <- prjs_solterms[[i]]$refmodel$y
    plt <- proj_linpred(prjs_solterms[[i]],
                        newdata = data.frame(y = y, x = x), transform = TRUE)
    plf <- proj_linpred(prjs_solterms[[i]],
                        newdata = data.frame(y = y, x = x), transform = FALSE)
    expect_equal(prjs_solterms[[!!i]]$family$linkinv(plf$pred),
                 plt$pred)
  }
})

test_that("proj_linpred: specifying integrated has an expected effect", {
  for (i in fam_nms) {
    y <- prjs_solterms[[i]]$refmodel$y
    plt <- proj_linpred(prjs_solterms[[i]],
                        newdata = data.frame(y = y, x = x),
                        integrated = TRUE)
    plf <- proj_linpred(prjs_solterms[[i]],
                        newdata = data.frame(y = y, x = x),
                        integrated = FALSE)
    expect_equal(
      prjs_solterms[[!!i]]$weights %*% plf$pred,
      plt$pred
    )
    expect_length(plt$lpd, length(plt$pred))
  }
})

test_that("proj_linpred: adding more regularization has an expected effect", {
  regul <- c(1e-6, 1e-1, 1e2)
  for (i in fam_nms) {
    norms <- rep(0, length(regul))
    for (j in 1:length(regul)) {
      y <- vs_list[[i]]$refmodel$y
      pred <- proj_linpred(vs_list[[i]],
                           nclusters = nclusters_pred_tst,
                           newdata = data.frame(y = y, x = x), nterms = 2,
                           transform = FALSE,
                           integrated = TRUE, regul = regul[j])
      norms[j] <- sum(pred$pred^2)
    }
    for (j in 1:(length(regul) - 1)) {
      expect_true(all(norms[!!j] >= norms[!!(j + 1)]), info = i)
    }
  }
})


test_that("proj_linpred: arguments passed to project work accordingly", {
  for (i in fam_nms) {
    y <- vs_list[[i]]$refmodel$y
    SW(pr <- project(vs_list[[i]],
                     nterms = c(2, 4), nclusters = nclusters_pred_tst,
                     regul = 1e-8, seed = 12))
    prl1 <- proj_linpred(pr, newdata = data.frame(y = y, x = x))
    SW(prl2 <- proj_linpred(vs_list[[i]],
                            nclusters = nclusters_pred_tst,
                            newdata = data.frame(y = y, x = x),
                            nterms = c(2, 4),
                            regul = 1e-8,
                            seed = 12))
    expect_equal(prl1$pred, prl2$pred, info = i)
  }
})

test_that(paste(
  "proj_linpred: providing newdata as a data frame works as expected"
), {
  SW(
    fit_form <- stan_glm(mpg ~ (drat + wt)^2,
                         data = mtcars, QR = TRUE,
                         chains = chains, seed = seed, iter = iter)
  )
  vs_form <- varsel(fit_form,
                    nclusters = nclusters_tst,
                    nclusters_pred = nclusters_pred_tst)
  p1 <- proj_linpred(vs_form, nclusters = nclusters_pred_tst,
                     newdata = mtcars, nterms = 3, seed = 2)
  x <- rstanarm::get_x(fit_form)[, -1]
  newdata <- data.frame(mpg = rstanarm::get_y(fit_form), x)
  p2 <- proj_linpred(vs_form, nclusters = nclusters_pred_tst,
                     newdata = newdata, nterms = 3,
                     seed = 2)
  expect_equal(p1$pred, p2$pred)
})


# -------------------------------------------------------------
context("proj_predict")

test_that("proj_predict: newdata is specified correctly", {
  ## expect_error(
  ##   proj_predict(prjs_solterms),
  ##   'argument "newdata" is missing, with no default'
  ## )
  ## expect_error(
  ##   proj_predict(prjs_solterms, newdata = NULL),
  ##   "must be a data.frame or a matrix"
  ## )
  expect_error(
    proj_predict(prjs_solterms, newdata = x[, 1]),
    "must be a data.frame or a matrix"
  )
  expect_error(
    proj_predict(prjs_solterms, newdata = data.frame(x = x),
                 solution_terms = paste0("x.", 1:1000)),
    paste("^The number of solution terms is greater than the number of",
          "columns in newdata\\.$")
  )
  expect_error(
    proj_predict(prjs_solterms,
                 newdata = data.frame(x = x)[, 1:2],
                 solution_terms = paste0("x.", 1:3)),
    paste("^The number of solution terms is greater than the number of",
          "columns in newdata\\.$")
  )
})

test_that(paste(
  "proj_predict: \"refmodel\" object as input leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    pl <- proj_predict(refmod_list[[i]],
                       nclusters = nclusters_pred_tst,
                       newdata = data.frame(x = x),
                       solution_terms = c("x.3", "x.5"))
    expect_identical(dim(pl), c(nresample_clusters_default, n), info = i)
  }
})

test_that(paste(
  "proj_predict: \"vsel\" object as input leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    pl <- proj_predict(vs_list[[i]],
                       nclusters = nclusters_pred_tst,
                       newdata = data.frame(x = x),
                       nterms = 0:nterms)
    expect_length(pl, nterms + 1)
    for (j in seq_along(pl)) {
      expect_identical(dim(pl[[!!j]]), c(nresample_clusters_default, n),
                       info = i)
    }
  }
})

test_that(paste(
  "proj_predict: \"projection\" object as input leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    pl <- proj_predict(prjs_solterms[[i]],
                       newdata = data.frame(x = x))
    expect_identical(dim(pl), c(nresample_clusters_default, n), info = i)
  }
})

test_that(paste(
  "proj_predict: \"proj_list\" object (an informal class) as input leads to",
  "correct output structure"
), {
  for (i in fam_nms) {
    pl <- proj_predict(proj_all_list[[i]], newdata = data.frame(x = x))
    expect_length(pl, nterms + 1)
    for (j in seq_along(pl)) {
      expect_identical(dim(pl[[!!j]]), c(nresample_clusters_default, n),
                       info = i)
    }
  }
})

test_that(paste(
  "proj_predict: output structure is also correct in edge cases",
  "(using `nclusters`)"
), {
  for (i in fam_nms) {
    for (n_crr in c(1L, 12L)) {
      for (nclusters_pred_crr in c(1L, 4L, 24L)) {
        for (nresample_clusters_crr in c(1L, 8L)) {
          pl <- proj_predict(
            refmod_list[[i]], nclusters = nclusters_pred_crr,
            newdata = head(data.frame(x = x), n_crr),
            nresample_clusters = nresample_clusters_crr,
            .seed = seed + 1,
            solution_terms = c("x.3", "x.5")
          )
          tstsetup <- unlist(nlist(i, n_crr, nclusters_pred_crr,
                                   nresample_clusters_crr))
          expect_identical(dim(pl), c(nresample_clusters_crr, n_crr),
                           info = tstsetup)
        }
      }
    }
  }
})

test_that(paste(
  "proj_predict: output structure is also correct in edge cases",
  "(using `ndraws`)"
), {
  for (i in fam_nms) {
    for (n_crr in c(1L, 12L)) {
      for (ndraws_pred_crr in c(1L, 4L, 24L)) {
        for (nresample_clusters_crr in c(1L, 8L)) {
          pl <- proj_predict(
            refmod_list[[i]], ndraws = ndraws_pred_crr,
            newdata = head(data.frame(x = x), n_crr),
            nresample_clusters = nresample_clusters_crr,
            .seed = seed + 1,
            solution_terms = c("x.3", "x.5")
          )
          tstsetup <- unlist(nlist(i, n_crr, ndraws_pred_crr,
                                   nresample_clusters_crr))
          nprjdraws_crr <- ifelse(ndraws_pred_crr <= 20,
                                  nresample_clusters_crr,
                                  ndraws_pred_crr)
          expect_identical(dim(pl), c(nprjdraws_crr, n_crr),
                           info = tstsetup)
        }
      }
    }
  }
})

test_that(paste(
  "proj_predict: error when varsel has not been performed on",
  "the object (and `solution_terms` is provided neither)"
), {
  expect_error(
    proj_predict(1, newdata = data.frame(x = x)),
    "is not an object of class \"vsel\""
  )
  expect_error(
    proj_predict(fit_gauss, newdata = data.frame(x = x)),
    "is not an object of class \"vsel\""
  )
  expect_error(
    proj_predict(c(prjs_solterms, list(x)),
                 newdata = data.frame(x = x)),
    "Invalid object supplied to argument `object`\\."
  )
})

## test_that("proj_predict: specifying ynew has an expected effect", {
##   for (i in seq_along(vs_list)) {
##     pl <- proj_predict(vs_list[[i]],
##                        nclusters = nclusters_pred_tst,
##                        newdata = data.frame(x = x),
##                        ynew = ys[[i]], nterms = 0:3)
##     pl2 <- proj_predict(vs_list[[i]],
##                         nclusters = nclusters_pred_tst,
##                         newdata = data.frame(x = x),
##                         nterms = 0:3)
##     for (j in seq_len(length(pl))) {
##       expect_equal(dim(pl[[j]]), dim(pl2[[j]]))
##     }
##   }
## })

## test_that(paste(
##   "proj_predict: specifying ynew as a factor works in a",
##   "binomial model"
## ), {
##   yfactor <- factor(rbinom(n, 1, 0.5))
##   pl <- proj_predict(vs_list[["binom"]],
##                      nclusters = nclusters_pred_tst,
##                      newdata = data.frame(x = x),
##                      ynew = yfactor)
##   expect_equal(ncol(pl), n)
##   expect_true(all(pl %in% c(0, 1)))
## })

test_that("proj_predict: specifying weightsnew has an expected effect", {
  pl <- proj_predict(prjs_solterms[["binom"]],
                     newdata = data.frame(x = x, weights = rep(1, NROW(x))),
                     seed = seed, .seed = seed)
  plw <- proj_predict(prjs_solterms[["binom"]],
                      newdata = data.frame(x = x, weights = weights),
                      seed = seed, .seed = seed,
                      weightsnew = ~weights)
  expect_true(sum(pl != plw) > 0)
})

test_that("proj_predict: specifying offsetnew has an expected effect", {
  for (i in seq_len(length(prjs_solterms))) {
    pl <- proj_predict(prjs_solterms[[i]],
                       newdata = data.frame(x = x),
                       seed = seed, .seed = seed)
    plo <- proj_predict(prjs_solterms[[i]],
                        newdata = data.frame(x = x, offset = offset),
                        seed = seed, .seed = seed, offsetnew = ~offset)
    expect_true(sum(pl != plo) > 0, info = i)
  }
})

test_that(paste(
  "proj_predict: specifying nresample_clusters has an expected effect"
), {
  for (i in fam_nms) {
    pl <- proj_predict(prjs_solterms[[i]],
                       nresample_clusters = nresample_clusters_tst,
                       newdata = data.frame(x = x))
    expect_equal(dim(pl), c(nresample_clusters_tst, n))
  }
})

test_that(paste(
  "proj_predict: specifying seed and .seed has an expected",
  "effect"
), {
  for (i in fam_nms) {
    pl1 <- proj_predict(prjs_solterms[[i]],
                        newdata = data.frame(x = x),
                        seed = seed, .seed = seed)
    pl2 <- proj_predict(prjs_solterms[[i]],
                        newdata = data.frame(x = x),
                        seed = seed, .seed = seed)
    expect_equal(pl1, pl2, info = i)
  }
})

test_that("proj_predict: arguments passed to project work accordingly", {
  for (i in fam_nms) {
    prp1 <- proj_predict(vs_list[[i]],
                         newdata = data.frame(x = x),
                         nresample_clusters = nresample_clusters_tst,
                         seed = 12, .seed = 12, nterms = c(2, 4),
                         nclusters = nclusters_pred_tst,
                         regul = 1e-08)
    prp2 <- proj_predict(vs_list[[i]],
                         newdata = data.frame(x = x),
                         nresample_clusters = nresample_clusters_tst,
                         nterms = c(2, 4),
                         nclusters = nclusters_pred_tst, regul = 1e-8,
                         seed = 12, .seed = 12)
    prp3 <- proj_predict(vs_list[[i]],
                         newdata = data.frame(x = x),
                         nresample_clusters = nresample_clusters_tst,
                         seed = 120, .seed = 120, nterms = c(2, 4),
                         nclusters = nclusters_pred_tst,
                         regul = 1e-08)
    expect_equal(prp1, prp2, info = i)
    expect_false(all(unlist(lapply(seq_along(prp1), function(i) {
      all(prp1[[i]] == prp3[[i]])
    }))),
    info = i
    )
  }
})
