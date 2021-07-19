# proj_linpred() ----------------------------------------------------------

context("proj_linpred()")

## object -----------------------------------------------------------------

test_that("`object` of class \"projection\" works", {
  for (tstsetup in names(prjs)) {
    pl_tester(pls[[tstsetup]],
              nprjdraws_expected = ndr_ncl_dtls(args_prj[[tstsetup]])$nprjdraws,
              info_str = tstsetup)
  }
})

test_that(paste(
  "`object` of (informal) class \"proj_list\" (created by varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs)) {
    tstsetup_vs <- args_prj_vs[[tstsetup]]$tstsetup
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- vss[[tstsetup_vs]]$suggested_size
    }
    pl_tester(pls_vs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_expected =
                ndr_ncl_dtls(args_prj_vs[[tstsetup]])$nprjdraws,
              info_str = tstsetup)
  }
})

test_that(paste(
  "`object` of (informal) class \"proj_list\" (created by cv_varsel()) works"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(prjs_cvvs)) {
    tstsetup_cvvs <- args_prj_cvvs[[tstsetup]]$tstsetup
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- cvvss[[tstsetup_cvvs]]$suggested_size
    }
    pl_tester(pls_cvvs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_expected =
                ndr_ncl_dtls(args_prj_cvvs[[tstsetup]])$nprjdraws,
              info_str = tstsetup)
  }
})

test_that(paste(
  "`object` of (informal) class \"proj_list\" (created manually) works"
), {
  tstsetups <- grep("^glm\\.gauss.*clust1", names(prjs), value = TRUE)
  stopifnot(length(tstsetups) > 1)
  pl <- proj_linpred(prjs[tstsetups])
  pl_tester(pl,
            len_expected = length(tstsetups),
            nprjdraws_expected = 1L,
            info_str = paste(tstsetups, collapse = ","))
})

test_that(paste(
  "`object` of class \"refmodel\" and passing arguments to project() works"
), {
  tstsetups <- grep("^glm\\.gauss\\.solterms_x\\.clust", names(prjs),
                    value = TRUE)[1]
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    pl_from_refmod <- do.call(proj_linpred, c(
      list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]]),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    pl_from_prj <- pls[[tstsetup]]
    expect_equal(pl_from_refmod, pl_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"stanreg\" and passing arguments to project() works"
), {
  tstsetups <- grep("^glm\\.gauss\\.solterms_x\\.clust", names(prjs),
                    value = TRUE)[1]
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    pl_from_fit <- do.call(proj_linpred, c(
      list(object = fits[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]]),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    pl_from_prj <- pls[[tstsetup]]
    expect_equal(pl_from_fit, pl_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by varsel()) and passing arguments",
  "to project() works"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("^glm\\.gauss\\.default_meth\\.subvec", names(prjs_vs),
                    value = TRUE)[1]
  for (tstsetup in tstsetups) {
    args_prj_vs_i <- args_prj_vs[[tstsetup]]
    pl_from_vsel <- do.call(proj_linpred, c(
      list(object = vss[[args_prj_vs_i$tstsetup]]),
      args_prj_vs_i[setdiff(names(args_prj_vs_i), c("tstsetup"))]
    ))
    pl_from_prj <- pls_vs[[tstsetup]]
    expect_equal(pl_from_vsel, pl_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by cv_varsel()) and passing arguments",
  "to project() works"
), {
  skip_if_not(run_cvvs)
  tstsetups <- grep("^glm\\.gauss\\.default_meth\\.default_cvmeth\\.subvec",
                    names(prjs_cvvs), value = TRUE)[1]
  for (tstsetup in tstsetups) {
    args_prj_cvvs_i <- args_prj_cvvs[[tstsetup]]
    pl_from_vsel <- do.call(proj_linpred, c(
      list(object = cvvss[[args_prj_cvvs_i$tstsetup]]),
      args_prj_cvvs_i[setdiff(names(args_prj_cvvs_i), c("tstsetup"))]
    ))
    pl_from_prj <- pls_cvvs[[tstsetup]]
    expect_equal(pl_from_vsel, pl_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "error if `object` is not of class \"vsel\" and `solution_terms` is provided",
  "neither"
), {
  expect_error(proj_linpred(1), "is not an object of class \"vsel\"")
  expect_error(proj_linpred(fits$glm$gauss),
               "is not an object of class \"vsel\"")
  expect_error(proj_linpred(refmods$glm$gauss),
               "is not an object of class \"vsel\"")
  expect_error(proj_linpred(c(prjs, list(dat))),
               "Invalid object supplied to argument `object`\\.")
})

## newdata and integrated -------------------------------------------------

test_that("incorrect `newdata` fails", {
  expect_error(
    proj_linpred(prjs, newdata = dat[, 1]),
    "must be a data\\.frame or a matrix"
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

test_that("`newdata` and `integrated` work (even in edge cases)", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    for (nobsv_crr in nobsv_tst) {
      pl_false <- proj_linpred(prjs[[tstsetup]],
                               newdata = head(dat, nobsv_crr))
      pl_tester(pl_false,
                nprjdraws_expected = ndr_ncl$nprjdraws,
                nobsv_expected = nobsv_crr,
                info_str = paste(tstsetup, nobsv_crr, sep = "__"))
      pl_true <- proj_linpred(prjs[[tstsetup]],
                              newdata = head(dat, nobsv_crr),
                              integrated = TRUE)
      pl_tester(pl_true,
                nprjdraws_expected = 1L,
                nobsv_expected = nobsv_crr,
                info_str = paste(tstsetup, nobsv_crr, "integrated", sep = "__"))
      expect_equal(prjs[[!!tstsetup]]$weights %*% pl_false$pred, pl_true$pred,
                   info = nobsv_crr)
    }
  }
})

test_that("`newdata` set to the original dataset doesn't change results", {
  for (tstsetup in names(prjs)) {
    pl_newdata <- proj_linpred(prjs[[tstsetup]], newdata = dat)
    pl_orig <- pls[[tstsetup]]
    expect_equal(pl_newdata, pl_orig, info = tstsetup)
  }
})

test_that(paste(
  "omitting the response in `newdata` causes output element `lpd` to be `NULL`"
), {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    resp_nm <- extract_terms_response(
      prjs[[tstsetup]]$refmodel$formula
    )$response
    stopifnot(!exists(resp_nm))
    pl <- proj_linpred(prjs[[tstsetup]],
                       newdata = dat[, setdiff(names(dat), resp_nm)])
    pl_tester(pl,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              lpd_null_expected = TRUE,
              info_str = tstsetup)
  }
})

## weightsnew -------------------------------------------------------------

test_that("`weightsnew` works", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    pl_orig <- pls[[tstsetup]]
    pl_ones <- proj_linpred(prjs[[tstsetup]],
                            newdata = dat_wobs_ones,
                            weightsnew = ~ wobs_col_ones)
    pl_tester(pl_ones,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    pl <- proj_linpred(prjs[[tstsetup]],
                       newdata = dat,
                       weightsnew = ~ wobs_col)
    pl_tester(pl,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    plw <- proj_linpred(prjs[[tstsetup]],
                        newdata = dat_wobs_new,
                        weightsnew = ~ wobs_col_new)
    pl_tester(plw,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    expect_equal(pl_ones$pred, pl_orig$pred, info = tstsetup)
    expect_equal(pl$pred, pl_orig$pred, info = tstsetup)
    expect_equal(plw$pred, pl_orig$pred, info = tstsetup)
    ### Note: This equivalence might in fact be undesired:
    expect_equal(pl_ones$lpd, pl_orig$lpd, info = tstsetup)
    ###
    ### Note: This inequality might in fact be undesired:
    expect_false(isTRUE(all.equal(pl$lpd, pl_orig$lpd)), info = tstsetup)
    ###
    expect_false(isTRUE(all.equal(plw$lpd, pl_orig$lpd)), info = tstsetup)
    expect_false(isTRUE(all.equal(plw$lpd, pl$lpd)), info = tstsetup)
  }
})

## offsetnew --------------------------------------------------------------

test_that("`offsetnew` works", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    pl_orig <- pls[[tstsetup]]
    pl_zeros <- proj_linpred(prjs[[tstsetup]],
                             newdata = dat_offs_zeros,
                             offsetnew = ~ offs_col_zeros)
    pl_tester(pl_zeros,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    pl <- proj_linpred(prjs[[tstsetup]],
                       newdata = dat,
                       offsetnew = ~ offs_col)
    pl_tester(pl,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    plo <- proj_linpred(prjs[[tstsetup]],
                        newdata = dat_offs_new,
                        offsetnew = ~ offs_col_new)
    pl_tester(plo,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    ### Note: This equivalence might in fact be undesired:
    expect_equal(pl_zeros, pl_orig, info = tstsetup)
    ###
    ### Note: This inequality might in fact be undesired:
    expect_false(isTRUE(all.equal(pl, pl_orig)), info = tstsetup)
    ###
    expect_equal(t(pl$pred) - dat$offs_col, t(pl_orig$pred),
                 info = tstsetup)
    expect_equal(t(plo$pred) - dat_offs_new$offs_col_new, t(pl_orig$pred),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(pl$lpd, pl_orig$lpd)), info = tstsetup)
    expect_false(isTRUE(all.equal(plo$lpd, pl_orig$lpd)), info = tstsetup)
    expect_false(isTRUE(all.equal(plo$lpd, pl$lpd)), info = tstsetup)
  }
})

## transform --------------------------------------------------------------

test_that("`transform` works", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    pl_false <- pls[[tstsetup]]
    pl_true <- proj_linpred(prjs[[tstsetup]], transform = TRUE)
    pl_tester(pl_true,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    expect_equal(prjs[[!!tstsetup]]$family$linkinv(pl_false$pred), pl_true$pred,
                 info = tstsetup)
  }
})

## regul ------------------------------------------------------------------

test_that("`regul` works", {
  regul_tst <- c(1e-6, 1e-1, 1e2)
  stopifnot(identical(regul_tst, sort(regul_tst)))
  tstsetups <- grep("^glm\\..*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    norms <- sapply(regul_tst, function(regul_crr) {
      pl <- do.call(proj_linpred, c(
        list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]],
             integrated = TRUE,
             regul = regul_crr),
        args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
      ))
      pl_tester(pl,
                nprjdraws_expected = 1L,
                info_str = tstsetup)
      return(sum(pl$pred^2))
    })
    for (j in head(seq_along(regul_tst), -1)) {
      expect_true(all(norms[!!j] >= norms[!!(j + 1)]), info = tstsetup)
    }
  }
})

## filter_nterms ----------------------------------------------------------

test_that(paste(
  "`filter_nterms` works correctly (for an `object` of class \"projection\")"
), {
  tstsetups <- grep("^glm\\.gauss\\.solterms_x\\.clust", names(prjs),
                    value = TRUE)[1]
  for (tstsetup in tstsetups) {
    nterms_avail_crr <- length(args_prj[[tstsetup]]$solution_terms)
    nterms_unavail_crr <- c(0L, nterms_avail_crr + 130L)
    stopifnot(!nterms_avail_crr %in% nterms_unavail_crr)
    for (filter_nterms_crr in nterms_unavail_crr) {
      expect_error(proj_linpred(prjs[[tstsetup]],
                                filter_nterms = filter_nterms_crr),
                   "subscript out of bounds",
                   info = paste(tstsetup, filter_nterms_crr, sep = "__"))
    }
    pl <- proj_linpred(prjs[[tstsetup]],
                       filter_nterms = nterms_avail_crr)
    pl_orig <- pls[[tstsetup]]
    expect_equal(pl, pl_orig, info = tstsetup)
  }
})

test_that(paste(
  "`filter_nterms` works correctly (for an `object` of (informal) class",
  "\"proj_list\")"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("^glm\\.gauss\\.default_meth\\.full$", names(prjs_vs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    # Unavailable number(s) of terms:
    for (filter_nterms_crr in nterms_unavail) {
      expect_error(proj_linpred(prjs_vs[[tstsetup]],
                                filter_nterms = filter_nterms_crr),
                   "subscript out of bounds",
                   info = paste(tstsetup,
                                paste(filter_nterms_crr, collapse = ","),
                                sep = "__"))
    }
    # Available number(s) of terms:
    nterms_avail_filter <- c(
      nterms_avail,
      list(partvec = c(nterms_max_tst %/% 2L, nterms_max_tst + 130L))
    )
    for (filter_nterms_crr in nterms_avail_filter) {
      pl_crr <- proj_linpred(prjs_vs[[tstsetup]],
                             filter_nterms = filter_nterms_crr)
      if (is.null(filter_nterms_crr)) filter_nterms_crr <- 0:nterms_max_tst
      nhits_nterms <- sum(filter_nterms_crr <= nterms_max_tst)
      pl_tester(pl_crr,
                len_expected = nhits_nterms,
                info_str = paste(tstsetup,
                                 paste(filter_nterms_crr, collapse = ","),
                                 sep = "__"))
      if (identical(filter_nterms_crr, 0:nterms_max_tst)) {
        # The special case of all possible numbers of terms:
        pl_orig <- pls_vs[[tstsetup]]
        expect_equal(pl_crr, pl_orig, info = tstsetup)
      }
    }
  }
})

# proj_predict() ----------------------------------------------------------

context("proj_predict()")

## seed -------------------------------------------------------------------

test_that("`.seed` works", {
  for (tstsetup in names(prjs)) {
    .Random.seed_orig1 <- .Random.seed
    pp_orig <- pps[[tstsetup]]
    .Random.seed_orig2 <- .Random.seed
    rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_new1 <- .Random.seed
    pp_new <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst + 1L)
    .Random.seed_new2 <- .Random.seed
    rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_repr1 <- .Random.seed
    pp_repr <- proj_predict(prjs[[tstsetup]], .seed = seed2_tst)
    .Random.seed_repr2 <- .Random.seed
    rand_repr <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_null1 <- .Random.seed
    pp_null <- proj_predict(prjs[[tstsetup]])
    .Random.seed_null2 <- .Random.seed

    expect_equal(pp_orig, pp_repr, info = tstsetup)
    expect_false(isTRUE(all.equal(pp_orig, pp_new)), info = tstsetup)
    expect_false(isTRUE(all.equal(pp_orig, pp_null)), info = tstsetup)
    expect_false(isTRUE(all.equal(pp_new, pp_null)), info = tstsetup)

    expect_equal(.Random.seed_orig2, .Random.seed_orig1, info = tstsetup)
    expect_equal(.Random.seed_new2, .Random.seed_new1, info = tstsetup)
    expect_equal(.Random.seed_repr2, .Random.seed_repr1, info = tstsetup)
    expect_equal(.Random.seed_null2, .Random.seed_null1, info = tstsetup)

    expect_false(isTRUE(all.equal(rand_new, rand_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(rand_repr, rand_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(rand_repr, rand_new)), info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_new2, .Random.seed_orig2)),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_repr2, .Random.seed_orig2)),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_repr2, .Random.seed_new2)),
                 info = tstsetup)
  }
})

## object -----------------------------------------------------------------

test_that("`object` of class \"projection\" works", {
  for (tstsetup in names(prjs)) {
    pp_tester(pps[[tstsetup]],
              nprjdraws_out_expected =
                ndr_ncl_dtls(args_prj[[tstsetup]])$nprjdraws_out,
              info_str = tstsetup)
  }
})

test_that(paste(
  "`object` of (informal) class \"proj_list\" (created by varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs)) {
    tstsetup_vs <- args_prj_vs[[tstsetup]]$tstsetup
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- vss[[tstsetup_vs]]$suggested_size
    }
    pp_tester(pps_vs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_out_expected =
                ndr_ncl_dtls(args_prj_vs[[tstsetup]])$nprjdraws_out,
              info_str = tstsetup)
  }
})

test_that(paste(
  "`object` of (informal) class \"proj_list\" (created by cv_varsel()) works"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(prjs_cvvs)) {
    tstsetup_cvvs <- args_prj_cvvs[[tstsetup]]$tstsetup
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- cvvss[[tstsetup_cvvs]]$suggested_size
    }
    pp_tester(pps_cvvs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_out_expected =
                ndr_ncl_dtls(args_prj_cvvs[[tstsetup]])$nprjdraws_out,
              info_str = tstsetup)
  }
})

test_that(paste(
  "`object` of (informal) class \"proj_list\" (created manually) works"
), {
  tstsetups <- grep("^glm\\.gauss.*clust1", names(prjs), value = TRUE)
  stopifnot(length(tstsetups) > 1)
  pp <- proj_predict(prjs[tstsetups], .seed = seed2_tst)
  pp_tester(pp,
            len_expected = length(tstsetups),
            info_str = paste(tstsetups, collapse = ","))
})

test_that(paste(
  "`object` of class \"refmodel\" and passing arguments to project() works"
), {
  tstsetups <- grep("^glm\\.gauss\\.solterms_x\\.clust", names(prjs),
                    value = TRUE)[1]
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    pp_from_refmod <- do.call(proj_predict, c(
      list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]],
           .seed = seed2_tst),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    pp_from_prj <- pps[[tstsetup]]
    expect_equal(pp_from_refmod, pp_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"stanreg\" and passing arguments to project() works"
), {
  tstsetups <- grep("^glm\\.gauss\\.solterms_x\\.clust", names(prjs),
                    value = TRUE)[1]
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    pp_from_fit <- do.call(proj_predict, c(
      list(object = fits[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]],
           .seed = seed2_tst),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    pp_from_prj <- pps[[tstsetup]]
    expect_equal(pp_from_fit, pp_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by varsel()) and passing arguments",
  "to project() works"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("^glm\\.gauss\\.default_meth\\.subvec", names(prjs_vs),
                    value = TRUE)[1]
  for (tstsetup in tstsetups) {
    args_prj_vs_i <- args_prj_vs[[tstsetup]]
    pp_from_vsel <- do.call(proj_predict, c(
      list(object = vss[[args_prj_vs_i$tstsetup]],
           .seed = seed2_tst),
      args_prj_vs_i[setdiff(names(args_prj_vs_i), c("tstsetup"))]
    ))
    pp_from_prj <- pps_vs[[tstsetup]]
    expect_equal(pp_from_vsel, pp_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by cv_varsel()) and passing arguments",
  "to project() works"
), {
  skip_if_not(run_cvvs)
  tstsetups <- grep("^glm\\.gauss\\.default_meth\\.default_cvmeth\\.subvec",
                    names(prjs_cvvs), value = TRUE)[1]
  for (tstsetup in tstsetups) {
    args_prj_cvvs_i <- args_prj_cvvs[[tstsetup]]
    pp_from_vsel <- do.call(proj_predict, c(
      list(object = cvvss[[args_prj_cvvs_i$tstsetup]],
           .seed = seed2_tst),
      args_prj_cvvs_i[setdiff(names(args_prj_cvvs_i), c("tstsetup"))]
    ))
    pp_from_prj <- pps_cvvs[[tstsetup]]
    expect_equal(pp_from_vsel, pp_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "error if `object` is not of class \"vsel\" and `solution_terms` is provided",
  "neither"
), {
  expect_error(proj_predict(1, .seed = seed2_tst),
               "is not an object of class \"vsel\"")
  expect_error(proj_predict(fits$glm$gauss, .seed = seed2_tst),
               "is not an object of class \"vsel\"")
  expect_error(proj_predict(refmods$glm$gauss, .seed = seed2_tst),
               "is not an object of class \"vsel\"")
  expect_error(proj_predict(c(prjs, list(dat)), .seed = seed2_tst),
               "Invalid object supplied to argument `object`\\.")
})

## newdata and nresample_clusters -----------------------------------------

test_that("incorrect `newdata` fails", {
  expect_error(
    proj_predict(prjs, newdata = dat[, 1], .seed = seed2_tst),
    "must be a data\\.frame or a matrix"
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

test_that("`newdata` and `nresample_clusters` work (even in edge cases)", {
  for (tstsetup in names(prjs)) {
    for (nobsv_crr in nobsv_tst) {
      for (nresample_clusters_crr in nresample_clusters_tst) {
        pp <- proj_predict(prjs[[tstsetup]],
                           newdata = head(dat, nobsv_crr),
                           nresample_clusters = nresample_clusters_crr,
                           .seed = seed2_tst)
        ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]],
                                nresample_clusters_crr = nresample_clusters_crr)
        pp_tester(pp,
                  nprjdraws_out_expected = ndr_ncl$nprjdraws_out,
                  nobsv_expected = nobsv_crr,
                  info_str = paste(tstsetup, nobsv_crr, nresample_clusters_crr,
                                   sep = "__"))
      }
    }
  }
})

test_that("`newdata` set to the original dataset doesn't change results", {
  for (tstsetup in names(prjs)) {
    pp_newdata <- proj_predict(prjs[[tstsetup]],
                               newdata = dat,
                               .seed = seed2_tst)
    pp_orig <- pps[[tstsetup]]
    expect_equal(pp_newdata, pp_orig, info = tstsetup)
  }
})

test_that("omitting the response in `newdata` doesn't change results", {
  for (tstsetup in names(prjs)) {
    resp_nm <- extract_terms_response(
      prjs[[tstsetup]]$refmodel$formula
    )$response
    stopifnot(!exists(resp_nm))
    pp_noresp <- proj_predict(prjs[[tstsetup]],
                              newdata = dat[, setdiff(names(dat), resp_nm)],
                              .seed = seed2_tst)
    pp_orig <- pps[[tstsetup]]
    expect_equal(pp_noresp, pp_orig, info = tstsetup)
  }
})

## weightsnew -------------------------------------------------------------

test_that("`weightsnew` works", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])

    pp_orig <- pps[[tstsetup]]
    expect_identical(dim(pp_orig), c(ndr_ncl$nprjdraws_out, nobsv),
                     info = tstsetup)

    pp_ones <- proj_predict(prjs[[tstsetup]],
                            newdata = dat_wobs_ones,
                            weightsnew = ~ wobs_col_ones,
                            .seed = seed2_tst)
    expect_identical(dim(pp_ones), c(ndr_ncl$nprjdraws_out, nobsv),
                     info = tstsetup)

    pp <- proj_predict(prjs[[tstsetup]],
                       newdata = dat,
                       weightsnew = ~ wobs_col,
                       .seed = seed2_tst)
    expect_identical(dim(pp), c(ndr_ncl$nprjdraws_out, nobsv), info = tstsetup)

    ppw <- proj_predict(prjs[[tstsetup]],
                        newdata = dat_wobs_new,
                        weightsnew = ~ wobs_col_new,
                        .seed = seed2_tst)
    expect_identical(dim(ppw), c(ndr_ncl$nprjdraws_out, nobsv), info = tstsetup)

    # Weights are only relevant for the binomial() family:
    if (args_prj[[tstsetup]]$fam_nm != "binom") {
      expect_equal(pp_ones, pp_orig, info = tstsetup)
      expect_equal(pp, pp_orig, info = tstsetup)
      expect_equal(ppw, pp_orig, info = tstsetup)
    } else {
      ### Note: This equivalence might in fact be undesired:
      expect_equal(pp_ones, pp_orig, info = tstsetup)
      ###
      ### Note: This inequality might in fact be undesired:
      expect_false(isTRUE(all.equal(pp, pp_orig)), info = tstsetup)
      ###
      expect_false(isTRUE(all.equal(ppw, pp_orig)), info = tstsetup)
      expect_false(isTRUE(all.equal(ppw, pp)), info = tstsetup)
    }
  }
})

test_that("`offsetnew` works", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])

    pp_orig <- pps[[tstsetup]]
    expect_identical(dim(pp_orig), c(ndr_ncl$nprjdraws_out, nobsv),
                     info = tstsetup)

    pp_zeros <- proj_predict(prjs[[tstsetup]],
                             newdata = dat_offs_zeros,
                             offsetnew = ~ offs_col_zeros,
                             .seed = seed2_tst)
    expect_identical(dim(pp_zeros), c(ndr_ncl$nprjdraws_out, nobsv),
                     info = tstsetup)

    pp <- proj_predict(prjs[[tstsetup]],
                       newdata = dat,
                       offsetnew = ~ offs_col,
                       .seed = seed2_tst)
    expect_identical(dim(pp), c(ndr_ncl$nprjdraws_out, nobsv), info = tstsetup)

    ppo <- proj_predict(prjs[[tstsetup]],
                        newdata = dat_offs_new,
                        offsetnew = ~ offs_col_new,
                        .seed = seed2_tst)
    expect_identical(dim(ppo), c(ndr_ncl$nprjdraws_out, nobsv), info = tstsetup)

    ### Note: This equivalence might in fact be undesired:
    expect_equal(pp_zeros, pp_orig, info = tstsetup)
    ###
    ### Note: This inequality might in fact be undesired:
    expect_false(isTRUE(all.equal(pp, pp_orig)), info = tstsetup)
    ###
    # For the gaussian() family, we can perform an easy check (because of the
    # identity link):
    if (args_prj[[tstsetup]]$fam_nm == "gauss") {
      expect_equal(t(pp) - dat$offs_col, t(pp_orig), info = tstsetup)
      expect_equal(t(ppo) - dat_offs_new$offs_col_new, t(pp_orig), info = tstsetup)
    } else {
      expect_false(isTRUE(all.equal(ppo, pp_orig)), info = tstsetup)
      expect_false(isTRUE(all.equal(ppo, pp)), info = tstsetup)
    }
  }
})

test_that(paste(
  "`filter_nterms` works correctly (for an `object` of class",
  "\"projection\")"
), {
  nterms_avail_crr <- length(solterms_x)
  nterms_unavail_crr <- c(0L, nterms_avail_crr + 130L)
  stopifnot(!nterms_avail_crr %in% nterms_unavail_crr)
  for (filter_nterms_crr in nterms_unavail_crr) {
    expect_error(proj_predict(prjs$glm.gauss.solterms_x.clust,
                              filter_nterms = !!filter_nterms_crr,
                              .seed = seed2_tst),
                 "subscript out of bounds")
  }
  pp <- proj_predict(prjs$glm.gauss.solterms_x.clust,
                     filter_nterms = nterms_avail_crr,
                     .seed = seed2_tst)
  pp_orig <- proj_predict(prjs$glm.gauss.solterms_x.clust,
                          .seed = seed2_tst)
  expect_equal(pp, pp_orig)
})

test_that(paste(
  "`filter_nterms` works correctly (for an `object` of",
  "(informal) class \"proj_list\")"
), {
  skip_if_not(run_vs)
  prjs_vs_crr <- prjs_vs$glm.gauss.default_meth.full
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
      expect_identical(dim(pp_crr[[!!j]]), c(nresample_clusters_default, nobsv),
                       info = tstsetup_crr)
    }
    if (identical(filter_nterms_crr, 0:nterms_max_tst)) {
      # The special case of all possible numbers of terms:
      pp_orig <- proj_predict(prjs_vs_crr, .seed = seed2_tst)
      expect_equal(pp_crr, pp_orig)
    }
  }
})
