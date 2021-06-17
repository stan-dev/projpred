# proj_linpred() ----------------------------------------------------------

context("proj_linpred()")

test_that("proj_linpred(): passing arguments to project() works correctly", {
  for (tstsetup in names(prjs)) {
    pl_from_prj <- proj_linpred(prjs[[tstsetup]])
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
  "proj_linpred(): `object` of class \"vsel\" (created by varsel()) leads",
  "to correct output structure"
), {
  skip_if_not(exists("vss"))
  tstsetups <- grep("^glm\\.gauss", names(vss), value = TRUE)
  stopifnot(length(tstsetups) == 1)
  for (tstsetup in tstsetups) {
    nterms_crr <- nterms_avail$subvec
    pl <- proj_linpred(vss[[tstsetup]],
                       nterms = nterms_crr,
                       nclusters = nclusters_pred_tst,
                       seed = seed_tst)
    expect_length(pl, length(nterms_crr))
    for (j in seq_along(pl)) {
      expect_named(pl[[!!j]], c("pred", "lpd"), info = tstsetup)
      expect_identical(dim(pl[[!!j]]$pred), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
      expect_identical(dim(pl[[!!j]]$lpd), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
    }
    expect_equal(pl, proj_linpred(prjs_vs$glm.gauss.subvec))
  }
})

test_that(paste(
  "proj_linpred(): `object` of class \"vsel\" (created by cv_varsel()) leads",
  "to correct output structure"
), {
  skip_if_not(exists("cvvss"))
  tstsetups <- grep("^glm\\.gauss", names(cvvss), value = TRUE)
  stopifnot(length(tstsetups) == 1)
  for (tstsetup in tstsetups) {
    nterms_crr <- nterms_avail$subvec
    pl <- proj_linpred(cvvss[[tstsetup]],
                       nterms = nterms_crr,
                       nclusters = nclusters_pred_tst,
                       seed = seed_tst)
    expect_length(pl, length(nterms_crr))
    for (j in seq_along(pl)) {
      expect_named(pl[[!!j]], c("pred", "lpd"), info = tstsetup)
      expect_identical(dim(pl[[!!j]]$pred), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
      expect_identical(dim(pl[[!!j]]$lpd), c(nclusters_pred_tst, n_tst),
                       info = tstsetup)
    }
    expect_equal(pl, proj_linpred(prjs_cvvs$glm.gauss.subvec))
  }
})

test_that(paste(
  "proj_linpred(): `object` of class \"projection\" leads to correct output",
  "structure"
), {
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    pl <- proj_linpred(prjs[[tstsetup]])
    expect_named(pl, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl$lpd), c(nprjdraws, n_tst), info = tstsetup)
  }
})

test_that(paste(
  "proj_linpred(): `object` of (informal) class \"proj_list\" (created by",
  "varsel()) leads to correct output structure"
), {
  skip_if_not(exists("prjs_vs"))
  for (tstsetup in names(prjs_vs)) {
    pl <- proj_linpred(prjs_vs[[tstsetup]])
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      tstsetup_vs <- grep(
        paste0("^", args_prj_vs[[tstsetup]]$mod_nm,
               "\\.", args_prj_vs[[tstsetup]]$fam_nm),
        names(vss),
        value = TRUE
      )
      stopifnot(length(tstsetup_vs) == 1)
      # Subtract 1L for the intercept:
      nterms_crr <- vss[[tstsetup_vs]]$suggested_size - 1L
    }
    if (length(nterms_crr) == 1) {
      # In fact, we don't have a "proj_list" object in this case, but since
      # incorporating this case is so easy, we create one:
      pl <- list(pl)
    }
    expect_length(pl, length(nterms_crr))
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
  "proj_linpred(): `object` of (informal) class \"proj_list\" (created by",
  "cv_varsel()) leads to correct output structure"
), {
  skip_if_not(exists("prjs_cvvs"))
  for (tstsetup in names(prjs_cvvs)) {
    pl <- proj_linpred(prjs_cvvs[[tstsetup]])
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      tstsetup_cvvs <- grep(
        paste0("^", args_prj_cvvs[[tstsetup]]$mod_nm,
               "\\.", args_prj_cvvs[[tstsetup]]$fam_nm),
        names(cvvss),
        value = TRUE
      )
      stopifnot(length(tstsetup_cvvs) == 1)
      # Subtract 1L for the intercept:
      nterms_crr <- cvvss[[tstsetup_cvvs]]$suggested_size - 1L
    }
    if (length(nterms_crr) == 1) {
      # In fact, we don't have a "proj_list" object in this case, but since
      # incorporating this case is so easy, we create one:
      pl <- list(pl)
    }
    expect_length(pl, length(nterms_crr))
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
  "proj_linpred(): error if `object` is not of class \"vsel\" and",
  "`solution_terms` is provided neither"
), {
  expect_error(proj_linpred(1), "is not an object of class \"vsel\"")
  expect_error(proj_linpred(fits$glm$gauss),
               "is not an object of class \"vsel\"")
  expect_error(proj_linpred(c(prjs, list(dat))),
               "Invalid object supplied to argument `object`\\.")
})

test_that("proj_linpred(): `newdata` is checked correctly", {
  expect_error(
    proj_linpred(prjs, newdata = dat[, 1]),
    "must be a data.frame or a matrix"
  )
  expect_error(
    proj_linpred(prjs,
                 solution_terms = rep_len(solterms_x, length.out = 1e4)),
    paste("^The number of solution terms is greater than the number of",
          "columns in newdata\\.$")
  )
  stopifnot(length(solterms_x) > 1)
  expect_error(
    proj_linpred(prjs[[grep("^glm\\.gauss", names(prjs))[1]]],
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
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    for (n_crr in c(1L, 12L)) {
      for (integrated_crr in c(FALSE, TRUE)) {
        tstsetup_crr <- unlist(nlist(tstsetup, n_crr, integrated_crr))
        pl <- proj_linpred(prjs[[tstsetup]],
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
  for (tstsetup in names(prjs)) {
    pl_orig <- proj_linpred(prjs[[tstsetup]])
    pl_newdata <- proj_linpred(prjs[[tstsetup]], newdata = dat)
    expect_equal(pl_orig, pl_newdata, info = tstsetup)
  }
})

test_that(paste(
  "proj_linpred(): omitting the response in `newdata` causes output element",
  "`lpd` to be `NULL`"
), {
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    resp_nm <- extract_terms_response(
      prjs[[tstsetup]]$refmodel$formula
    )$response
    stopifnot(!exists(resp_nm))
    pl <- proj_linpred(prjs[[tstsetup]],
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
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }

    pl_orig <- proj_linpred(prjs[[tstsetup]])
    expect_named(pl_orig, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl_orig$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl_orig$lpd), c(nprjdraws, n_tst), info = tstsetup)

    pl_ones <- proj_linpred(prjs[[tstsetup]],
                            newdata = dat_ones,
                            weightsnew = ~ wobs_col_ones)
    expect_named(pl_ones, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl_ones$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl_ones$lpd), c(nprjdraws, n_tst), info = tstsetup)

    pl <- proj_linpred(prjs[[tstsetup]],
                       newdata = dat,
                       weightsnew = ~ wobs_col)
    expect_named(pl, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl$lpd), c(nprjdraws, n_tst), info = tstsetup)

    plw <- proj_linpred(prjs[[tstsetup]],
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
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }

    pl_orig <- proj_linpred(prjs[[tstsetup]])
    expect_named(pl_orig, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl_orig$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl_orig$lpd), c(nprjdraws, n_tst), info = tstsetup)

    pl_zeros <- proj_linpred(prjs[[tstsetup]],
                             newdata = dat_zeros,
                             offsetnew = ~ offs_col_zeros)
    expect_named(pl_zeros, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl_zeros$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl_zeros$lpd), c(nprjdraws, n_tst), info = tstsetup)

    pl <- proj_linpred(prjs[[tstsetup]],
                       newdata = dat,
                       offsetnew = ~ offs_col)
    expect_named(pl, c("pred", "lpd"), info = tstsetup)
    expect_identical(dim(pl$pred), c(nprjdraws, n_tst), info = tstsetup)
    expect_identical(dim(pl$lpd), c(nprjdraws, n_tst), info = tstsetup)

    plo <- proj_linpred(prjs[[tstsetup]],
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
    expect_equal(t(pl_orig$pred), t(pl$pred) - dat$offs_col,
                 info = tstsetup)
    expect_equal(t(pl_orig$pred), t(plo$pred) - dat_new$offs_col_new,
                 info = tstsetup)
    expect_false(isTRUE(all.equal(pl_orig$lpd, pl$lpd)), info = tstsetup)
    expect_false(isTRUE(all.equal(pl_orig$lpd, plo$lpd)), info = tstsetup)
    expect_false(isTRUE(all.equal(pl$lpd, plo$lpd)), info = tstsetup)
  }
})

test_that("proj_linpred(): `transform` has an expected effect", {
  for (tstsetup in names(prjs)) {
    plt <- proj_linpred(prjs[[tstsetup]], transform = TRUE)
    plf <- proj_linpred(prjs[[tstsetup]], transform = FALSE)
    expect_equal(prjs[[!!tstsetup]]$family$linkinv(plf$pred), plt$pred,
                 info = tstsetup)
  }
})

test_that("proj_linpred(): `integrated` has an expected effect", {
  for (tstsetup in names(prjs)) {
    plt <- proj_linpred(prjs[[tstsetup]], integrated = TRUE)
    plf <- proj_linpred(prjs[[tstsetup]], integrated = FALSE)
    expect_equal(prjs[[!!tstsetup]]$weights %*% plf$pred, plt$pred,
                 info = tstsetup)
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
      expect_true(all(norms[!!j] >= norms[!!(j + 1)]), info = fam_nm)
    }
  }
})

test_that(paste(
  "proj_linpred(): `filter_nterms` works correctly (for an `object` of class",
  "\"projection\")"
), {
  pl_orig <- proj_linpred(prjs$glm.gauss.solterms_x.clust)
  nterms_avail_x <- length(solterms_x)
  nterms_unavail_x <- c(0L, nterms_avail_x + 130L)
  stopifnot(!nterms_avail_x %in% nterms_unavail_x)
  for (filter_nterms_crr in nterms_unavail_x) {
    expect_error(proj_linpred(prjs$glm.gauss.solterms_x.clust,
                              filter_nterms = !!filter_nterms_crr),
                 "subscript out of bounds")
  }
  pl <- proj_linpred(prjs$glm.gauss.solterms_x.clust,
                     filter_nterms = nterms_avail_x)
  expect_equal(pl_orig, pl)
})

test_that(paste(
  "proj_linpred(): `filter_nterms` works correctly (for an `object` of",
  "(informal) class \"proj_list\")"
), {
  skip_if_not(exists("prjs_vs"))
  prjs_vs_crr <- prjs_vs$glm.gauss.full
  # The special case of all possible numbers of terms:
  pl_orig <- proj_linpred(prjs_vs_crr)
  # Unavailable number(s) of terms:
  for (filter_nterms_crr in nterms_unavail) {
    expect_error(proj_linpred(prjs_vs_crr,
                              filter_nterms = !!filter_nterms_crr),
                 "subscript out of bounds")
  }
  # Available number(s) of terms:
  nterms_avail_filter <- c(
    nterms_avail,
    list(partvec = c(nterms_max_tst %/% 2L, nterms_max_tst + 130L))
  )
  for (filter_nterms_crr in nterms_avail_filter) {
    tstsetup_crr <- paste(filter_nterms_crr, collapse = ", ")
    pl_crr <- proj_linpred(prjs_vs_crr,
                           filter_nterms = filter_nterms_crr)
    if (is.null(filter_nterms_crr)) filter_nterms_crr <- 0:nterms_max_tst
    nhits_nterms <- sum(filter_nterms_crr <= nterms_max_tst)
    if (nhits_nterms == 1) pl_crr <- list(pl_crr)
    expect_length(pl_crr, nhits_nterms)
    for (j in seq_along(pl_crr)) {
      expect_named(pl_crr[[!!j]], c("pred", "lpd"), info = tstsetup_crr)
      expect_identical(dim(pl_crr[[!!j]]$pred), c(nclusters_pred_tst, n_tst),
                       info = tstsetup_crr)
      expect_identical(dim(pl_crr[[!!j]]$lpd), c(nclusters_pred_tst, n_tst),
                       info = tstsetup_crr)
    }
    if (identical(filter_nterms_crr, 0:nterms_max_tst)) {
      expect_equal(pl_orig, pl_crr)
    }
  }
})

# proj_predict() ----------------------------------------------------------

context("proj_predict()")

test_that("proj_predict(): `.seed` has an expected effect", {
  for (tstsetup in names(prjs)) {
    pp1 <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    pp2 <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    pp3 <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst + 1L)
    pp4 <- proj_predict(prjs[[tstsetup]])
    expect_equal(pp1, pp2, info = tstsetup)
    expect_false(isTRUE(all.equal(pp1, pp3)), info = tstsetup)
    expect_false(isTRUE(all.equal(pp1, pp4)), info = tstsetup)
    expect_false(isTRUE(all.equal(pp3, pp4)), info = tstsetup)
  }
})

test_that("proj_predict(): passing arguments to project() works correctly", {
  for (tstsetup in names(prjs)) {
    pp_from_prj <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    args_prj_i <- args_prj[[tstsetup]]
    pp_direct <- do.call(proj_predict, c(
      list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]],
           .seed = seed2_tst),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    expect_equal(pp_from_prj, pp_direct, info = tstsetup)
  }
})

test_that(paste(
  "proj_predict(): `object` of class \"refmodel\" leads to correct output",
  "structure"
), {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      tstsetup <- unlist(nlist(mod_nm, fam_nm))
      pp <- proj_predict(refmods[[mod_nm]][[fam_nm]],
                         .seed = seed2_tst,
                         solution_terms = solterms_x,
                         nclusters = nclusters_pred_tst,
                         seed = seed_tst)
      expect_identical(dim(pp), c(nresample_clusters_default, n_tst),
                       info = tstsetup)
    }
  }
})

test_that(paste(
  "proj_predict(): `object` of class \"vsel\" (created by varsel()) leads",
  "to correct output structure"
), {
  skip_if_not(exists("vss"))
  tstsetups <- grep("^glm\\.gauss", names(vss), value = TRUE)
  stopifnot(length(tstsetups) == 1)
  for (tstsetup in tstsetups) {
    nterms_crr <- nterms_avail$subvec
    pp <- proj_predict(vss[[tstsetup]],
                       .seed = seed2_tst,
                       nterms = nterms_crr,
                       nclusters = nclusters_pred_tst,
                       seed = seed_tst)
    expect_length(pp, length(nterms_crr))
    for (j in seq_along(pp)) {
      expect_identical(dim(pp[[!!j]]), c(nresample_clusters_default, n_tst),
                       info = tstsetup)
    }
    expect_equal(pp, proj_predict(prjs_vs$glm.gauss.subvec, .seed = seed2_tst))
  }
})

test_that(paste(
  "proj_predict(): `object` of class \"vsel\" (created by cv_varsel()) leads",
  "to correct output structure"
), {
  skip_if_not(exists("cvvss"))
  tstsetups <- grep("^glm\\.gauss", names(cvvss), value = TRUE)
  stopifnot(length(tstsetups) == 1)
  for (tstsetup in tstsetups) {
    nterms_crr <- nterms_avail$subvec
    pp <- proj_predict(cvvss[[tstsetup]],
                       .seed = seed2_tst,
                       nterms = nterms_crr,
                       nclusters = nclusters_pred_tst,
                       seed = seed_tst)
    expect_length(pp, length(nterms_crr))
    for (j in seq_along(pp)) {
      expect_identical(dim(pp[[!!j]]), c(nresample_clusters_default, n_tst),
                       info = tstsetup)
    }
    expect_equal(pp, proj_predict(prjs_cvvs$glm.gauss.subvec, .seed = seed2_tst))
  }
})

test_that(paste(
  "proj_predict(): `object` of class \"projection\" leads to correct output",
  "structure"
), {
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    if (ndr_ncl_nm == "nclusters" || nprjdraws <= 20) {
      nprjdraws_out <- nresample_clusters_default
    } else {
      nprjdraws_out <- nprjdraws
    }
    pp <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    expect_identical(dim(pp), c(nprjdraws_out, n_tst), info = tstsetup)
  }
})

test_that(paste(
  "proj_predict(): `object` of (informal) class \"proj_list\" (created by",
  "varsel()) leads to correct output structure"
), {
  skip_if_not(exists("prjs_vs"))
  for (tstsetup in names(prjs_vs)) {
    pp <- proj_predict(prjs_vs[[tstsetup]], .seed = seed2_tst)
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      tstsetup_vs <- grep(
        paste0("^", args_prj_vs[[tstsetup]]$mod_nm,
               "\\.", args_prj_vs[[tstsetup]]$fam_nm),
        names(vss),
        value = TRUE
      )
      stopifnot(length(tstsetup_vs) == 1)
      # Subtract 1L for the intercept:
      nterms_crr <- vss[[tstsetup_vs]]$suggested_size - 1L
    }
    if (length(nterms_crr) == 1) {
      # In fact, we don't have a "proj_list" object in this case, but since
      # incorporating this case is so easy, we create one:
      pp <- list(pp)
    }
    expect_length(pp, length(nterms_crr))
    for (j in seq_along(pp)) {
      expect_identical(dim(pp[[!!j]]), c(nresample_clusters_default, n_tst),
                       info = tstsetup)
    }
  }
})

test_that(paste(
  "proj_predict(): `object` of (informal) class \"proj_list\" (created by",
  "cv_varsel()) leads to correct output structure"
), {
  skip_if_not(exists("prjs_cvvs"))
  for (tstsetup in names(prjs_cvvs)) {
    pp <- proj_predict(prjs_cvvs[[tstsetup]], .seed = seed2_tst)
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      tstsetup_cvvs <- grep(
        paste0("^", args_prj_cvvs[[tstsetup]]$mod_nm,
               "\\.", args_prj_cvvs[[tstsetup]]$fam_nm),
        names(cvvss),
        value = TRUE
      )
      stopifnot(length(tstsetup_cvvs) == 1)
      # Subtract 1L for the intercept:
      nterms_crr <- cvvss[[tstsetup_cvvs]]$suggested_size - 1L
    }
    if (length(nterms_crr) == 1) {
      # In fact, we don't have a "proj_list" object in this case, but since
      # incorporating this case is so easy, we create one:
      pp <- list(pp)
    }
    expect_length(pp, length(nterms_crr))
    for (j in seq_along(pp)) {
      expect_identical(dim(pp[[!!j]]), c(nresample_clusters_default, n_tst),
                       info = tstsetup)
    }
  }
})

test_that(paste(
  "proj_predict(): error if `object` is not of class \"vsel\" and",
  "`solution_terms` is provided neither"
), {
  expect_error(proj_predict(1, .seed = seed2_tst),
               "is not an object of class \"vsel\"")
  expect_error(proj_predict(fits$glm$gauss, .seed = seed2_tst),
               "is not an object of class \"vsel\"")
  expect_error(proj_predict(c(prjs, list(dat)), .seed = seed2_tst),
               "Invalid object supplied to argument `object`\\.")
})

test_that("proj_predict(): `newdata` is checked correctly", {
  expect_error(
    proj_predict(prjs, newdata = dat[, 1], .seed = seed2_tst),
    "must be a data.frame or a matrix"
  )
  expect_error(
    proj_predict(prjs,
                 .seed = seed2_tst,
                 solution_terms = rep_len(solterms_x, length.out = 1e4)),
    paste("^The number of solution terms is greater than the number of",
          "columns in newdata\\.$")
  )
  stopifnot(length(solterms_x) > 1)
  expect_error(
    proj_predict(prjs[[grep("^glm\\.gauss", names(prjs))[1]]],
                 newdata = dat[, 1, drop = FALSE],
                 .seed = seed2_tst,
                 solution_terms = solterms_x),
    paste("^The number of solution terms is greater than the number of",
          "columns in newdata\\.$")
  )
})

test_that(paste(
  "proj_predict(): `newdata` and `nresample_clusters` lead to correct output",
  "structure (even in edge cases)"
), {
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    for (n_crr in c(1L, 12L)) {
      for (nresample_clusters_crr in c(1L, 100L)) {
        tstsetup_crr <- unlist(nlist(tstsetup, n_crr, nresample_clusters_crr))
        pp <- proj_predict(prjs[[tstsetup]],
                           newdata = head(dat, n_crr),
                           nresample_clusters = nresample_clusters_crr,
                           .seed = seed2_tst)
        if (ndr_ncl_nm == "nclusters" || nprjdraws <= 20) {
          nprjdraws_crr <- nresample_clusters_crr
        } else {
          nprjdraws_crr <- nprjdraws
        }
        expect_identical(dim(pp), c(nprjdraws_crr, n_crr), info = tstsetup_crr)
      }
    }
  }
})

test_that(paste(
  "proj_predict(): `newdata` set to the original dataset doesn't change results"
), {
  for (tstsetup in names(prjs)) {
    pp_orig <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    pp_newdata <- proj_predict(prjs[[tstsetup]],
                               newdata = dat,
                               .seed = seed2_tst)
    expect_equal(pp_orig, pp_newdata, info = tstsetup)
  }
})

test_that(paste(
  "proj_predict(): omitting the response in `newdata` doesn't change results"
), {
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    pp_orig <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    resp_nm <- extract_terms_response(
      prjs[[tstsetup]]$refmodel$formula
    )$response
    stopifnot(!exists(resp_nm))
    pp_noresp <- proj_predict(prjs[[tstsetup]],
                              newdata = dat[, setdiff(names(dat), resp_nm)],
                              .seed = seed2_tst)
    expect_equal(pp_orig, pp_noresp, info = tstsetup)
  }
})

test_that("proj_predict(): `weightsnew` has an expected effect", {
  dat_ones <- within(dat, {
    wobs_col <- NULL
    wobs_col_ones <- rep_len(1, length.out = n_tst)
  })
  dat_new <- within(dat, {
    wobs_col <- NULL
    wobs_col_new <- rep_len(2:5, length.out = n_tst)
  })
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    if (ndr_ncl_nm == "nclusters" || nprjdraws <= 20) {
      nprjdraws_out <- nresample_clusters_default
    } else {
      nprjdraws_out <- nprjdraws
    }

    pp_orig <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    expect_identical(dim(pp_orig), c(nprjdraws_out, n_tst), info = tstsetup)

    pp_ones <- proj_predict(prjs[[tstsetup]],
                            newdata = dat_ones,
                            weightsnew = ~ wobs_col_ones,
                            .seed = seed2_tst)
    expect_identical(dim(pp_ones), c(nprjdraws_out, n_tst), info = tstsetup)

    pp <- proj_predict(prjs[[tstsetup]],
                       newdata = dat,
                       weightsnew = ~ wobs_col,
                       .seed = seed2_tst)
    expect_identical(dim(pp), c(nprjdraws_out, n_tst), info = tstsetup)

    ppw <- proj_predict(prjs[[tstsetup]],
                        newdata = dat_new,
                        weightsnew = ~ wobs_col_new,
                        .seed = seed2_tst)
    expect_identical(dim(ppw), c(nprjdraws_out, n_tst), info = tstsetup)

    # Weights are only relevant for the binomial() family:
    if (args_prj[[tstsetup]]$fam_nm != "binom") {
      expect_equal(pp_orig, pp_ones, info = tstsetup)
      expect_equal(pp_orig, pp, info = tstsetup)
      expect_equal(pp_orig, ppw, info = tstsetup)
    } else {
      ### Note: This equivalence might in fact be undesired:
      expect_equal(pp_orig, pp_ones, info = tstsetup)
      ###
      ### Note: This inequality might in fact be undesired:
      expect_false(isTRUE(all.equal(pp_orig, pp)), info = tstsetup)
      ###
      expect_false(isTRUE(all.equal(pp_orig, ppw)), info = tstsetup)
      expect_false(isTRUE(all.equal(pp, ppw)), info = tstsetup)
    }
  }
})

test_that("proj_predict(): `offsetnew` has an expected effect", {
  dat_zeros <- within(dat, {
    offs_col <- NULL
    offs_col_zeros <- rep_len(0, length.out = n_tst)
  })
  dat_new <- within(dat, {
    offs_col <- NULL
    offs_col_new <- seq(-2, 2, length.out = n_tst)
  })
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      nprjdraws <- ndraws_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    if (ndr_ncl_nm == "nclusters" || nprjdraws <= 20) {
      nprjdraws_out <- nresample_clusters_default
    } else {
      nprjdraws_out <- nprjdraws
    }

    pp_orig <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    expect_identical(dim(pp_orig), c(nprjdraws_out, n_tst), info = tstsetup)

    pp_zeros <- proj_predict(prjs[[tstsetup]],
                             newdata = dat_zeros,
                             offsetnew = ~ offs_col_zeros,
                             .seed = seed2_tst)
    expect_identical(dim(pp_zeros), c(nprjdraws_out, n_tst),
                     info = tstsetup)

    pp <- proj_predict(prjs[[tstsetup]],
                       newdata = dat,
                       offsetnew = ~ offs_col,
                       .seed = seed2_tst)
    expect_identical(dim(pp), c(nprjdraws_out, n_tst), info = tstsetup)

    ppo <- proj_predict(prjs[[tstsetup]],
                        newdata = dat_new,
                        offsetnew = ~ offs_col_new,
                        .seed = seed2_tst)
    expect_identical(dim(ppo), c(nprjdraws_out, n_tst), info = tstsetup)

    ### Note: This equivalence might in fact be undesired:
    expect_equal(pp_orig, pp_zeros, info = tstsetup)
    ###
    ### Note: This inequality might in fact be undesired:
    expect_false(isTRUE(all.equal(pp_orig, pp)), info = tstsetup)
    ###
    # For the gaussian() family, we can perform an easy check (because of the
    # identity link):
    if (args_prj[[tstsetup]]$fam_nm == "gauss") {
      expect_equal(t(pp_orig), t(pp) - dat$offs_col, info = tstsetup)
      expect_equal(t(pp_orig), t(ppo) - dat_new$offs_col_new, info = tstsetup)
    } else {
      expect_false(isTRUE(all.equal(pp_orig, ppo)), info = tstsetup)
      expect_false(isTRUE(all.equal(pp, ppo)), info = tstsetup)
    }
  }
})

test_that(paste(
  "proj_predict(): `filter_nterms` works correctly (for an `object` of class",
  "\"projection\")"
), {
  pp_orig <- proj_predict(prjs$glm.gauss.solterms_x.clust,
                          .seed = seed2_tst)
  nterms_avail_x <- length(solterms_x)
  nterms_unavail_x <- c(0L, nterms_avail_x + 130L)
  stopifnot(!nterms_avail_x %in% nterms_unavail_x)
  for (filter_nterms_crr in nterms_unavail_x) {
    expect_error(proj_predict(prjs$glm.gauss.solterms_x.clust,
                              filter_nterms = !!filter_nterms_crr,
                              .seed = seed2_tst),
                 "subscript out of bounds")
  }
  pp <- proj_predict(prjs$glm.gauss.solterms_x.clust,
                     filter_nterms = nterms_avail_x,
                     .seed = seed2_tst)
  expect_equal(pp_orig, pp)
})

test_that(paste(
  "proj_predict(): `filter_nterms` works correctly (for an `object` of",
  "(informal) class \"proj_list\")"
), {
  skip_if_not(exists("prjs_vs"))
  prjs_vs_crr <- prjs_vs$glm.gauss.full
  # The special case of all possible numbers of terms:
  pp_orig <- proj_predict(prjs_vs_crr, .seed = seed2_tst)
  # Unavailable number(s) of terms:
  for (filter_nterms_crr in nterms_unavail) {
    expect_error(proj_predict(prjs_vs_crr,
                              filter_nterms = !!filter_nterms_crr,
                              .seed = seed2_tst),
                 "subscript out of bounds")
  }
  # Available number(s) of terms:
  nterms_avail_filter <- c(
    nterms_avail,
    list(partvec = c(nterms_max_tst %/% 2L, nterms_max_tst + 130L))
  )
  for (filter_nterms_crr in nterms_avail_filter) {
    tstsetup_crr <- paste(filter_nterms_crr, collapse = ", ")
    pp_crr <- proj_predict(prjs_vs_crr,
                           filter_nterms = filter_nterms_crr,
                           .seed = seed2_tst)
    if (is.null(filter_nterms_crr)) filter_nterms_crr <- 0:nterms_max_tst
    nhits_nterms <- sum(filter_nterms_crr <= nterms_max_tst)
    if (nhits_nterms == 1) pp_crr <- list(pp_crr)
    expect_length(pp_crr, nhits_nterms)
    for (j in seq_along(pp_crr)) {
      expect_identical(dim(pp_crr[[!!j]]), c(nresample_clusters_default, n_tst),
                       info = tstsetup_crr)
    }
    if (identical(filter_nterms_crr, 0:nterms_max_tst)) {
      expect_equal(pp_orig, pp_crr)
    }
  }
})
