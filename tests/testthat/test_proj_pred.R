# proj_linpred() ----------------------------------------------------------

context("proj_linpred()")

test_that("proj_linpred(): passing arguments to project() works correctly", {
  for (tstsetup in names(prjs_solterms)) {
    pl_from_prj <- proj_linpred(prjs_solterms[[tstsetup]])
    args_prj_i <- args_prj[[tstsetup]]
    pl_direct <- do.call(proj_linpred, c(
      list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]]),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    expect_equal(pl_from_prj, pl_direct, info = tstsetup)
  }
})

test_that(paste(
  "proj_linpred(): `object` of class \"refmodel\" leads to correct output",
  "structure"
), {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      tstsetup <- unlist(nlist(mod_nm, fam_nm))
      pl <- proj_linpred(refmods[[mod_nm]][[fam_nm]],
                         solution_terms = solterms_x,
                         nclusters = nclusters_pred_tst,
                         seed = seed_tst)
      expect_named(pl, c("pred", "lpd"), info = tstsetup)
      expect_identical(dim(pl$pred), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
      expect_identical(dim(pl$lpd), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
    }
  }
})

test_that(paste(
  "proj_linpred(): `object` of class \"vsel\" leads to correct output",
  "structure"
), {
  for (tstsetup in grep("^glm\\.gauss", names(vss), value = TRUE)) {
    nterms_max_crr <- args_vs[[tstsetup]]$nterms_max
    pl <- proj_linpred(vss[[tstsetup]],
                       nterms = 0:nterms_max_crr,
                       nclusters = nclusters_pred_tst,
                       seed = seed_tst)
    expect_length(pl, nterms_max_crr + 1)
    for (j in seq_along(pl)) {
      expect_named(pl[[!!j]], c("pred", "lpd"), info = tstsetup)
      expect_identical(dim(pl[[!!j]]$pred), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
      expect_identical(dim(pl[[!!j]]$lpd), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
    }
  }
})

test_that(paste(
  "proj_linpred(): `object` of class \"projection\" leads to correct output",
  "structure"
), {
  for (tstsetup in names(prjs_solterms)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    stopifnot(length(ndr_ncl_nm) == 1)
    nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    pl <- proj_linpred(prjs_solterms[[tstsetup]])
    expect_named(pl, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl$lpd), c(nprjdraws, n_tst), info = tstsetup)
  }
})

test_that(paste(
  "proj_linpred(): `object` of (informal) class \"proj_list\" leads to correct",
  "output structure"
), {
  pl <- proj_linpred(prj_nterms)
  expect_length(pl, nterms_max_tst + 1)
  for (j in seq_along(pl)) {
    expect_named(pl[[!!j]], c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl[[!!j]]$pred), c(nclusters_pred_tst, n_tst),
                     info = tstsetup)
    expect_identical(dim(pl[[!!j]]$lpd), c(nclusters_pred_tst, n_tst),
                     info = tstsetup)
  }
})

test_that(paste(
  "proj_linpred(): error if `object` is not of class \"vsel\" (and",
  "`solution_terms` is provided neither)"
), {
  expect_error(proj_linpred(1), "is not an object of class \"vsel\"")
  expect_error(proj_linpred(fits$glm$gauss),
               "is not an object of class \"vsel\"")
  expect_error(proj_linpred(c(prjs_solterms, list(dat))),
               "Invalid object supplied to argument `object`\\.")
})

test_that("proj_linpred(): `newdata` is checked correctly", {
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
    proj_linpred(prjs_solterms[[grep("^glm\\.gauss", names(prjs_solterms))[1]]],
                 newdata = dat[, 1, drop = FALSE],
                 solution_terms = solterms_x),
    paste("^The number of solution terms is greater than the number of",
          "columns in newdata\\.$")
  )
})

test_that(paste(
  "proj_linpred(): `newdata` and `integrated` lead to correct output structure",
  "(even in edge cases)"
), {
  for (tstsetup in names(prjs_solterms)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    stopifnot(length(ndr_ncl_nm) == 1)
    nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    for (n_crr in c(1L, 12L)) {
      for (integrated_crr in c(FALSE, TRUE)) {
        tstsetup_crr <- unlist(nlist(tstsetup, n_crr, integrated_crr))
        pl <- proj_linpred(prjs_solterms[[tstsetup]],
                           newdata = head(dat, n_crr),
                           integrated = integrated_crr)
        expect_named(pl, c("pred", "lpd"), info = tstsetup_crr)
        nprjdraws_crr <- ifelse(integrated_crr, 1L, nprjdraws)
        expect_identical(dim(pl$pred), c(nprjdraws_crr, n_crr),
                         info = tstsetup_crr)
        expect_identical(dim(pl$lpd), c(nprjdraws_crr, n_crr),
                         info = tstsetup_crr)
      }
    }
  }
})

test_that(paste(
  "proj_linpred(): `newdata` set to the original dataset doesn't change results"
), {
  for (tstsetup in names(prjs_solterms)) {
    pl_orig <- proj_linpred(prjs_solterms[[tstsetup]])
    pl_newdata <- proj_linpred(prjs_solterms[[tstsetup]], newdata = dat)
    expect_equal(pl_orig, pl_newdata, info = tstsetup)
  }
})

test_that(paste(
  "proj_linpred(): omitting the response in `newdata` causes output element",
  "`lpd` to be `NULL`"
), {
  for (tstsetup in names(prjs_solterms)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    stopifnot(length(ndr_ncl_nm) == 1)
    nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    resp_nm <- extract_terms_response(
      prjs_solterms[[tstsetup]]$refmodel$formula
    )$response
    stopifnot(!exists(resp_nm))
    pl <- proj_linpred(prjs_solterms[[tstsetup]],
                       newdata = dat[, setdiff(names(dat), resp_nm)])
    expect_named(pl, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_null(pl$lpd, info = tstsetup)
  }
})

test_that("proj_linpred(): `weightsnew` has an expected effect", {
  dat_ones <- within(dat, {
    wobs_col <- NULL
    wobs_col_ones <- rep_len(1, length.out = n_tst)
  })
  dat_new <- within(dat, {
    wobs_col <- NULL
    wobs_col_new <- rep_len(2:5, length.out = n_tst)
  })
  for (tstsetup in names(prjs_solterms)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    stopifnot(length(ndr_ncl_nm) == 1)
    nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]

    pl_orig <- proj_linpred(prjs_solterms[[tstsetup]])
    expect_named(pl_orig, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl_orig$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl_orig$lpd), c(nprjdraws, n_tst), info = tstsetup)

    pl_ones <- proj_linpred(prjs_solterms[[tstsetup]],
                            newdata = dat_ones,
                            weightsnew = ~ wobs_col_ones)
    expect_named(pl_ones, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl_ones$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl_ones$lpd), c(nprjdraws, n_tst), info = tstsetup)

    pl <- proj_linpred(prjs_solterms[[tstsetup]],
                       newdata = dat,
                       weightsnew = ~ wobs_col)
    expect_named(pl, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl$lpd), c(nprjdraws, n_tst), info = tstsetup)

    plw <- proj_linpred(prjs_solterms[[tstsetup]],
                        newdata = dat_new,
                        weightsnew = ~ wobs_col_new)
    expect_named(plw, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(plw$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(plw$lpd), c(nprjdraws, n_tst), info = tstsetup)

    expect_equal(pl_orig$pred, pl_ones$pred, info = tstsetup)
    expect_equal(pl_orig$pred, pl$pred, info = tstsetup)
    expect_equal(pl_orig$pred, plw$pred, info = tstsetup)
    ### Note: This equivalence might in fact be undesired:
    expect_equal(pl_orig$lpd, pl_ones$lpd, info = tstsetup)
    ###
    ### Note: This inequality might in fact be undesired:
    expect_false(isTRUE(all.equal(pl_orig$lpd, pl$lpd)), info = tstsetup)
    ###
    expect_false(isTRUE(all.equal(pl_orig$lpd, plw$lpd)), info = tstsetup)
    expect_false(isTRUE(all.equal(pl$lpd, plw$lpd)), info = tstsetup)
  }
})

test_that("proj_linpred(): `offsetnew` has an expected effect", {
  dat_zeros <- within(dat, {
    offs_col <- NULL
    offs_col_zeros <- rep_len(0, length.out = n_tst)
  })
  dat_new <- within(dat, {
    offs_col <- NULL
    offs_col_new <- seq(-2, 2, length.out = n_tst)
  })
  for (tstsetup in names(prjs_solterms)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    stopifnot(length(ndr_ncl_nm) == 1)
    nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]

    pl_orig <- proj_linpred(prjs_solterms[[tstsetup]])
    expect_named(pl_orig, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl_orig$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl_orig$lpd), c(nprjdraws, n_tst), info = tstsetup)

    pl_zeros <- proj_linpred(prjs_solterms[[tstsetup]],
                             newdata = dat_zeros,
                             offsetnew = ~ offs_col_zeros)
    expect_named(pl_zeros, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl_zeros$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl_zeros$lpd), c(nprjdraws, n_tst), info = tstsetup)

    pl <- proj_linpred(prjs_solterms[[tstsetup]],
                       newdata = dat,
                       offsetnew = ~ offs_col)
    expect_named(pl, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl$lpd), c(nprjdraws, n_tst), info = tstsetup)

    plo <- proj_linpred(prjs_solterms[[tstsetup]],
                        newdata = dat_new,
                        offsetnew = ~ offs_col_new)
    expect_named(plo, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(plo$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(plo$lpd), c(nprjdraws, n_tst), info = tstsetup)

    ### Note: This equivalence might in fact be undesired:
    expect_equal(pl_orig, pl_zeros, info = tstsetup)
    ###
    ### Note: This inequality might in fact be undesired:
    expect_false(isTRUE(all.equal(pl_orig, pl)), info = tstsetup)
    ###
    expect_equal(t(pl_orig$pred), t(pl$pred) - dat$offs_col)
    expect_equal(t(pl_orig$pred), t(plo$pred) - dat_new$offs_col_new)
    expect_false(isTRUE(all.equal(pl_orig$lpd, pl$lpd)), info = tstsetup)
    expect_false(isTRUE(all.equal(pl_orig$lpd, plo$lpd)), info = tstsetup)
    expect_false(isTRUE(all.equal(pl$lpd, plo$lpd)), info = tstsetup)
  }
})

test_that("proj_linpred(): `transform` has an expected effect", {
  for (tstsetup in names(prjs_solterms)) {
    plt <- proj_linpred(prjs_solterms[[tstsetup]], transform = TRUE)
    plf <- proj_linpred(prjs_solterms[[tstsetup]], transform = FALSE)
    expect_equal(prjs_solterms[[!!tstsetup]]$family$linkinv(plf$pred), plt$pred)
  }
})

test_that("proj_linpred(): `integrated` has an expected effect", {
  for (tstsetup in names(prjs_solterms)) {
    plt <- proj_linpred(prjs_solterms[[tstsetup]], integrated = TRUE)
    plf <- proj_linpred(prjs_solterms[[tstsetup]], integrated = FALSE)
    expect_equal(prjs_solterms[[!!tstsetup]]$weights %*% plf$pred, plt$pred)
  }
})

test_that("proj_linpred(): `regul` has an expected effect", {
  regul_tst <- c(1e-6, 1e-1, 1e2)
  stopifnot(identical(regul_tst, sort(regul_tst)))
  for (fam_nm in fam_nms) {
    norms <- sapply(regul_tst, function(regul_crr) {
      pl <- proj_linpred(refmods$glm[[fam_nm]],
                         integrated = TRUE,
                         solution_terms = solterms_x,
                         nclusters = nclusters_pred_tst,
                         seed = seed_tst,
                         regul = regul_crr)
      return(sum(pl$pred^2))
    })
    for (j in head(seq_along(regul_tst), -1)) {
      expect_true(all(norms[!!j] >= norms[!!(j + 1)]), info = tstsetup)
    }
  }
})

# proj_predict() ----------------------------------------------------------

context("proj_predict()")

test_that("proj_predict(): `newdata` is checked correctly", {
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
  "proj_predict(): `object` of class \"refmodel\" leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    pl <- proj_predict(refmod_list[[i]],
                       nclusters = nclusters_pred_tst,
                       newdata = data.frame(x = x),
                       solution_terms = c("x.3", "x.5"))
    expect_identical(dim(pl), c(nresample_clusters_default, n_tst), info = tstsetup)
  }
})

test_that(paste(
  "proj_predict(): `object` of class \"vsel\" leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    pl <- proj_predict(vs_list[[i]],
                       nclusters = nclusters_pred_tst,
                       newdata = data.frame(x = x),
                       nterms = 0:nterms)
    expect_length(pl, nterms + 1)
    for (j in seq_along(pl)) {
      expect_identical(dim(pl[[!!j]]), c(nresample_clusters_default, n_tst),
                       info = i)
    }
  }
})

test_that(paste(
  "proj_predict(): `object` of class \"projection\" leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    pl <- proj_predict(prjs_solterms[[i]],
                       newdata = data.frame(x = x))
    expect_identical(dim(pl), c(nresample_clusters_default, n_tst), info = tstsetup)
  }
})

test_that(paste(
  "proj_predict(): `object` of (informal) class \"proj_list\" leads to correct",
  "output structure"
), {
  for (i in fam_nms) {
    pl <- proj_predict(proj_all_list[[i]], newdata = data.frame(x = x))
    expect_length(pl, nterms + 1)
    for (j in seq_along(pl)) {
      expect_identical(dim(pl[[!!j]]), c(nresample_clusters_default, n_tst),
                       info = i)
    }
  }
})

test_that(paste(
  "proj_predict(): `newdata` and `nresample_clusters` lead to correct output",
  "structure (even in edge cases)",
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
  "proj_predict(): `newdata` and `nresample_clusters` lead to correct output",
  "structure (even in edge cases)",
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
  "proj_predict(): error if `object` is not of class \"vsel\" (and",
  "`solution_terms` is provided neither)"
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

## test_that("proj_predict(): specifying ynew has an expected effect", {
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
##   "proj_predict(): specifying ynew as a factor works in a",
##   "binomial model"
## ), {
##   yfactor <- factor(rbinom(n, 1, 0.5))
##   pl <- proj_predict(vs_list[["binom"]],
##                      nclusters = nclusters_pred_tst,
##                      newdata = data.frame(x = x),
##                      ynew = yfactor)
##   expect_equal(ncol(pl), n_tst)
##   expect_true(all(pl %in% c(0, 1)))
## })

test_that("proj_predict(): `weightsnew` has an expected effect", {
  pl <- proj_predict(prjs_solterms[["binom"]],
                     newdata = data.frame(x = x, weights = rep(1, NROW(x))),
                     seed = seed, .seed = seed)
  plw <- proj_predict(prjs_solterms[["binom"]],
                      newdata = data.frame(x = x, weights = weights),
                      seed = seed, .seed = seed,
                      weightsnew = ~weights)
  expect_true(sum(pl != plw) > 0)
})

test_that("proj_predict(): `offsetnew` has an expected effect", {
  for (i in seq_len(length(prjs_solterms))) {
    pl <- proj_predict(prjs_solterms[[i]],
                       newdata = data.frame(x = x),
                       seed = seed, .seed = seed)
    plo <- proj_predict(prjs_solterms[[i]],
                        newdata = data.frame(x = x, offset = offset),
                        seed = seed, .seed = seed, offsetnew = ~offset)
    expect_true(sum(pl != plo) > 0, info = tstsetup)
  }
})

test_that("proj_predict(): `nresample_clusters` has an expected effect", {
  for (i in fam_nms) {
    pl <- proj_predict(prjs_solterms[[i]],
                       nresample_clusters = nresample_clusters_tst,
                       newdata = data.frame(x = x))
    expect_equal(dim(pl), c(nresample_clusters_tst, n_tst))
  }
})

test_that("proj_predict(): `seed` and `.seed` have an expected effect", {
  for (i in fam_nms) {
    pl1 <- proj_predict(prjs_solterms[[i]],
                        newdata = data.frame(x = x),
                        seed = seed, .seed = seed)
    pl2 <- proj_predict(prjs_solterms[[i]],
                        newdata = data.frame(x = x),
                        seed = seed, .seed = seed)
    expect_equal(pl1, pl2, info = tstsetup)
  }
})

test_that("proj_predict(): passing arguments to project() works correctly", {
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
    expect_equal(prp1, prp2, info = tstsetup)
    expect_false(all(unlist(lapply(seq_along(prp1), function(i) {
      all(prp1[[i]] == prp3[[i]])
    }))),
    info = i
    )
  }
})
