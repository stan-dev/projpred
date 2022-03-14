# proj_linpred() ----------------------------------------------------------

context("proj_linpred()")

## object -----------------------------------------------------------------

test_that("pl: `object` of class \"projection\" works", {
  for (tstsetup in names(prjs)) {
    pl_tester(pls[[tstsetup]],
              nprjdraws_expected = ndr_ncl_dtls(args_prj[[tstsetup]])$nprjdraws,
              info_str = tstsetup)
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(pls[[tstsetup]]))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "pl: `object` of (informal) class \"proj_list\" (based on varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs)) {
    tstsetup_vs <- args_prj_vs[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- vss[[tstsetup_vs]]$suggested_size
    }
    pl_tester(pls_vs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_expected =
                ndr_ncl_dtls(args_prj_vs[[tstsetup]])$nprjdraws,
              info_str = tstsetup)
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(pls_vs[[tstsetup]]))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "pl: `object` of (informal) class \"proj_list\" (based on cv_varsel()) works"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(prjs_cvvs)) {
    tstsetup_cvvs <- args_prj_cvvs[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- cvvss[[tstsetup_cvvs]]$suggested_size
    }
    pl_tester(pls_cvvs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_expected =
                ndr_ncl_dtls(args_prj_cvvs[[tstsetup]])$nprjdraws,
              info_str = tstsetup)
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(pls_cvvs[[tstsetup]]))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "`object` of (informal) class \"proj_list\" (created manually) works"
), {
  tstsetups <- grep("\\.glm\\.gauss\\..*\\.clust$", names(prjs), value = TRUE)
  stopifnot(length(tstsetups) > 1)
  pl <- proj_linpred(prjs[tstsetups])
  pl_tester(pl,
            len_expected = length(tstsetups),
            info_str = paste(tstsetups, collapse = ","))
})

test_that(paste(
  "`object` of class \"refmodel\" and passing arguments to project() works"
), {
  tstsetups <- grep("\\.glm\\.gauss.*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    pl_from_refmod <- do.call(proj_linpred, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]]),
      excl_nonargs(args_prj_i)
    ))
    pl_from_prj <- pls[[tstsetup]]
    expect_equal(pl_from_refmod, pl_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"stanreg\" or \"brmsfit\" and passing arguments to",
  "project() works"
), {
  tstsetups <- grep("\\.glm\\.gauss.*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    pl_from_fit <- do.call(proj_linpred, c(
      list(object = fits[[args_prj_i$tstsetup_fit]]),
      excl_nonargs(args_prj_i)
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
  tstsetups <- head(grep("\\.glm\\.gauss.*\\.default_meth\\.subvec",
                         names(prjs_vs), value = TRUE),
                    1)
  for (tstsetup in tstsetups) {
    args_prj_vs_i <- args_prj_vs[[tstsetup]]
    pl_from_vsel <- do.call(proj_linpred, c(
      list(object = vss[[args_prj_vs_i$tstsetup_vsel]]),
      excl_nonargs(args_prj_vs_i)
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
  tstsetups <- head(
    grep("\\.glm\\.gauss.*\\.default_meth\\.default_cvmeth\\.subvec",
         names(prjs_cvvs), value = TRUE),
    1
  )
  for (tstsetup in tstsetups) {
    args_prj_cvvs_i <- args_prj_cvvs[[tstsetup]]
    pl_from_vsel <- do.call(proj_linpred, c(
      list(object = cvvss[[args_prj_cvvs_i$tstsetup_vsel]]),
      excl_nonargs(args_prj_cvvs_i)
    ))
    pl_from_prj <- pls_cvvs[[tstsetup]]
    expect_equal(pl_from_vsel, pl_from_prj, info = tstsetup)
  }
})

test_that("`object` not of class \"vsel\" and missing `solution_terms` fails", {
  expect_error(
    proj_linpred(1),
    paste("^Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.$")
  )
  expect_error(
    proj_linpred(fits[[1]]),
    paste("^Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.$")
  )
  expect_error(
    proj_linpred(refmods[[1]]),
    paste("^Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.$")
  )
  expect_error(
    proj_linpred(c(prjs, list(dat))),
    paste("Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.")
  )
})

## newdata and integrated -------------------------------------------------

test_that("invalid `newdata` fails", {
  expect_error(
    proj_linpred(prjs, newdata = dat[, 1]),
    "must be a data\\.frame or a matrix"
  )
  stopifnot(length(solterms_x) > 1)
  expect_error(
    proj_linpred(prjs[[head(grep("\\.glm\\.gauss.*\\.solterms_x", names(prjs)),
                            1)]],
                 newdata = dat[, head(solterms_x, -1), drop = FALSE]),
    "^object '.*' not found$"
  )
})

test_that("`newdata` and `integrated` work (even in edge cases)", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    dat_crr <- get_dat(tstsetup)
    for (nobsv_crr in nobsv_tst) {
      if (args_prj[[tstsetup]]$mod_nm == "gamm") {
        # TODO (GAMMs): Fix this.
        next
      }
      pl_false <- proj_linpred(prjs[[tstsetup]],
                               newdata = head(dat_crr, nobsv_crr))
      pl_tester(pl_false,
                nprjdraws_expected = ndr_ncl$nprjdraws,
                nobsv_expected = nobsv_crr,
                info_str = paste(tstsetup, nobsv_crr, sep = "__"))
      pl_true <- proj_linpred(prjs[[tstsetup]],
                              newdata = head(dat_crr, nobsv_crr),
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
    dat_crr <- get_dat(tstsetup)
    pl_newdata <- proj_linpred(prjs[[tstsetup]], newdata = dat_crr)
    pl_orig <- pls[[tstsetup]]
    expect_equal(pl_newdata, pl_orig, info = tstsetup)
  }
})

test_that(paste(
  "omitting the response in `newdata` (not possible for `\"brmsfit\"`-based",
  "reference models) causes output element `lpd` to be `NULL` but doesn't",
  "change results otherwise"
), {
  for (tstsetup in names(prjs)) {
    if (args_prj[[tstsetup]]$pkg_nm == "brms") {
      next
    }
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    resp_nm <- extract_terms_response(
      prjs[[tstsetup]]$refmodel$formula
    )$response
    stopifnot(!exists(resp_nm))
    dat_crr <- get_dat(tstsetup)
    pl_noresp <- proj_linpred(
      prjs[[tstsetup]],
      newdata = dat_crr[, setdiff(names(dat_crr), resp_nm)]
    )
    pl_tester(pl_noresp,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              lpd_null_expected = TRUE,
              info_str = tstsetup)
    pl_orig <- pls[[tstsetup]]
    expect_equal(pl_noresp$pred, pl_orig$pred, info = tstsetup)
  }
})

## weightsnew -------------------------------------------------------------

test_that("`weightsnew` works", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    pl_orig <- pls[[tstsetup]]
    if (!(args_prj[[tstsetup]]$pkg_nm == "brms" &&
          args_prj[[tstsetup]]$fam_nm == "binom")) {
      # TODO (brms): Fix or document why this doesn't work for "brmsfit"s.
      pl_ones <- proj_linpred(prjs[[tstsetup]],
                              newdata = get_dat(tstsetup, dat_wobs_ones),
                              weightsnew = ~ wobs_col_ones)
      pl_tester(pl_ones,
                nprjdraws_expected = ndr_ncl$nprjdraws,
                info_str = tstsetup)
    }
    pl <- proj_linpred(prjs[[tstsetup]],
                       newdata = get_dat(tstsetup, dat),
                       weightsnew = ~ wobs_col)
    pl_tester(pl,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    if (!(args_prj[[tstsetup]]$pkg_nm == "brms" &&
          args_prj[[tstsetup]]$fam_nm == "binom")) {
      # TODO (brms): Fix or document why this doesn't work for "brmsfit"s.
      plw <- proj_linpred(prjs[[tstsetup]],
                          newdata = get_dat(tstsetup, dat_wobs_new),
                          weightsnew = ~ wobs_col_new)
      pl_tester(plw,
                nprjdraws_expected = ndr_ncl$nprjdraws,
                info_str = tstsetup)
    }
    if (!(args_prj[[tstsetup]]$pkg_nm == "brms" &&
          args_prj[[tstsetup]]$fam_nm == "binom")) {
      expect_equal(pl_ones$pred, pl_orig$pred, info = tstsetup)
    }
    expect_equal(pl$pred, pl_orig$pred, info = tstsetup)
    if (!(args_prj[[tstsetup]]$pkg_nm == "brms" &&
          args_prj[[tstsetup]]$fam_nm == "binom")) {
      expect_equal(plw$pred, pl_orig$pred, info = tstsetup)
    }
    if (grepl("\\.with_wobs", tstsetup)) {
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        ### TODO: This might in fact be undesired:
        expect_equal(pl_ones$lpd, pl_orig$lpd, info = tstsetup)
        expect_false(isTRUE(all.equal(pl$lpd, pl_orig$lpd)), info = tstsetup)
        ###
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        expect_false(isTRUE(all.equal(pl_ones$lpd, pl_orig$lpd)),
                     info = tstsetup)
        expect_equal(pl$lpd, pl_orig$lpd, info = tstsetup)
      }
      expect_false(isTRUE(all.equal(plw$lpd, pl_orig$lpd)), info = tstsetup)
      expect_false(isTRUE(all.equal(pl$lpd, pl_ones$lpd)), info = tstsetup)
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        expect_false(isTRUE(all.equal(plw$lpd, pl_ones$lpd)), info = tstsetup)
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        expect_equal(plw$lpd, pl_ones$lpd, info = tstsetup)
      }
      expect_false(isTRUE(all.equal(plw$lpd, pl$lpd)), info = tstsetup)
    } else {
      if (!(args_prj[[tstsetup]]$pkg_nm == "brms" &&
            args_prj[[tstsetup]]$fam_nm == "binom")) {
        expect_equal(pl_ones$lpd, pl_orig$lpd, info = tstsetup)
      }
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        expect_false(isTRUE(all.equal(pl$lpd, pl_orig$lpd)), info = tstsetup)
        expect_false(isTRUE(all.equal(plw$lpd, pl_orig$lpd)), info = tstsetup)
        expect_false(isTRUE(all.equal(pl$lpd, pl_ones$lpd)), info = tstsetup)
        expect_false(isTRUE(all.equal(plw$lpd, pl_ones$lpd)), info = tstsetup)
        expect_false(isTRUE(all.equal(plw$lpd, pl$lpd)), info = tstsetup)
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        ### TODO: This might in fact be undesired:
        expect_equal(pl$lpd, pl_orig$lpd, info = tstsetup)
        if (args_prj[[tstsetup]]$fam_nm != "binom") {
          expect_equal(plw$lpd, pl_orig$lpd, info = tstsetup)
          expect_equal(pl$lpd, pl_ones$lpd, info = tstsetup)
          expect_equal(plw$lpd, pl_ones$lpd, info = tstsetup)
          expect_equal(plw$lpd, pl$lpd, info = tstsetup)
        }
        ###
      }
    }
  }
})

## offsetnew --------------------------------------------------------------

test_that("`offsetnew` works", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    pl_orig <- pls[[tstsetup]]
    if (args_prj[[tstsetup]]$pkg_nm != "brms") {
      # TODO (brms): Fix or document why this doesn't work for "brmsfit"s.
      pl_zeros <- proj_linpred(prjs[[tstsetup]],
                               newdata = get_dat(tstsetup, dat_offs_zeros),
                               offsetnew = ~ offs_col_zeros)
      pl_tester(pl_zeros,
                nprjdraws_expected = ndr_ncl$nprjdraws,
                info_str = tstsetup)
    }
    pl <- proj_linpred(prjs[[tstsetup]],
                       newdata = get_dat(tstsetup, dat),
                       offsetnew = ~ offs_col)
    pl_tester(pl,
              nprjdraws_expected = ndr_ncl$nprjdraws,
              info_str = tstsetup)
    if (args_prj[[tstsetup]]$pkg_nm != "brms") {
      # TODO (brms): Fix or document why this doesn't work for "brmsfit"s.
      plo <- proj_linpred(prjs[[tstsetup]],
                          newdata = get_dat(tstsetup, dat_offs_new),
                          offsetnew = ~ offs_col_new)
      pl_tester(plo,
                nprjdraws_expected = ndr_ncl$nprjdraws,
                info_str = tstsetup)
    }
    if (grepl("\\.with_offs", tstsetup)) {
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        ### TODO: This might in fact be undesired:
        expect_equal(pl_zeros, pl_orig, info = tstsetup)
        expect_false(isTRUE(all.equal(pl, pl_orig)), info = tstsetup)
        expect_equal(t(pl$pred) - dat$offs_col, t(pl_orig$pred),
                     info = tstsetup)
        expect_equal(t(plo$pred) - dat_offs_new$offs_col_new, t(pl_orig$pred),
                     info = tstsetup)
        expect_false(isTRUE(all.equal(pl$lpd, pl_orig$lpd)), info = tstsetup)
        ###
        expect_false(isTRUE(all.equal(plo$lpd, pl_orig$lpd)), info = tstsetup)
        expect_false(isTRUE(all.equal(plo$lpd, pl$lpd)), info = tstsetup)
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        expect_equal(pl, pl_orig, info = tstsetup)
      }
    } else {
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        expect_equal(pl_zeros, pl_orig, info = tstsetup)
        if (args_prj[[tstsetup]]$fam_nm %in% c("brnll", "binom")) {
          # To avoid failing tests due to numerical inaccuracies for extreme
          # values:
          is_extreme <- which(abs(pl_orig$pred) > f_binom$linkfun(1 - 1e-12),
                              arr.ind = TRUE)
          pl_orig$pred[is_extreme] <- NA
          pl$pred[is_extreme] <- NA
          plo$pred[is_extreme] <- NA
        }
        expect_equal(t(pl$pred) - dat$offs_col, t(pl_orig$pred),
                     info = tstsetup)
        expect_equal(t(plo$pred) - dat_offs_new$offs_col_new, t(pl_orig$pred),
                     info = tstsetup)
        expect_false(isTRUE(all.equal(pl$lpd, pl_orig$lpd)), info = tstsetup)
        expect_false(isTRUE(all.equal(plo$lpd, pl_orig$lpd)), info = tstsetup)
        expect_false(isTRUE(all.equal(plo$lpd, pl$lpd)), info = tstsetup)
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        ### TODO: This might in fact be undesired:
        expect_equal(pl, pl_orig, info = tstsetup)
        ###
      }
    }
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
    expect_equal(prjs[[!!tstsetup]]$refmodel$family$linkinv(pl_false$pred),
                 pl_true$pred,
                 info = tstsetup)
  }
})

## regul ------------------------------------------------------------------

test_that("`regul` works", {
  regul_tst <- c(1e-6, 1e-1, 1e2)
  stopifnot(identical(regul_tst, sort(regul_tst)))
  tstsetups <- grep("\\.glm\\..*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    norms <- sapply(regul_tst, function(regul_crr) {
      pl <- do.call(proj_linpred, c(
        list(object = refmods[[args_prj_i$tstsetup_ref]],
             integrated = TRUE,
             regul = regul_crr),
        excl_nonargs(args_prj_i)
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

test_that("`filter_nterms` works (for an `object` of class \"projection\")", {
  tstsetups <- grep("\\.glm\\.gauss.*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    nterms_avail_crr <- length(args_prj[[tstsetup]]$solution_terms)
    nterms_unavail_crr <- c(0L, nterms_avail_crr + 130L)
    stopifnot(!nterms_avail_crr %in% nterms_unavail_crr)
    for (filter_nterms_crr in nterms_unavail_crr) {
      expect_error(proj_linpred(prjs[[tstsetup]],
                                filter_nterms = filter_nterms_crr),
                   "Invalid `filter_nterms`\\.",
                   info = paste(tstsetup, filter_nterms_crr, sep = "__"))
    }
    pl <- proj_linpred(prjs[[tstsetup]],
                       filter_nterms = nterms_avail_crr)
    pl_orig <- pls[[tstsetup]]
    expect_equal(pl, pl_orig, info = tstsetup)
  }
})

test_that(paste(
  "`filter_nterms` works (for an `object` of (informal) class \"proj_list\")"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.glm\\.gauss.*\\.default_meth\\.full$", names(prjs_vs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    # Unavailable number(s) of terms:
    for (filter_nterms_crr in nterms_unavail) {
      expect_error(proj_linpred(prjs_vs[[tstsetup]],
                                filter_nterms = filter_nterms_crr),
                   "Invalid `filter_nterms`\\.",
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

## Single observation, single draw ----------------------------------------

test_that(paste(
  "a single observation and a single draw work (which implicitly tests",
  "this edge case for family$ll_fun(), too)"
), {
  for (tstsetup in names(prjs)) {
    if (args_prj[[tstsetup]]$mod_nm == "gamm") {
      # TODO (GAMMs): Fix this.
      next
    }
    pl1 <- proj_linpred(refmods[[args_prj[[tstsetup]]$tstsetup_ref]],
                        newdata = head(get_dat(tstsetup), 1),
                        solution_terms = args_prj[[tstsetup]]$solution_terms,
                        nclusters = 1L,
                        seed = seed_tst)
    pl_tester(pl1,
              nprjdraws_expected = 1L,
              nobsv_expected = 1L,
              info_str = tstsetup)
  }
})

# proj_predict() ----------------------------------------------------------

context("proj_predict()")

## seed -------------------------------------------------------------------

test_that("`.seed` works (and restores the RNG state afterwards)", {
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

test_that("pp: `object` of class \"projection\" works", {
  for (tstsetup in names(prjs)) {
    pp_tester(pps[[tstsetup]],
              nprjdraws_out_expected =
                ndr_ncl_dtls(args_prj[[tstsetup]])$nprjdraws_out,
              info_str = tstsetup)
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(pps[[tstsetup]]))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "pp: `object` of (informal) class \"proj_list\" (based on varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs)) {
    tstsetup_vs <- args_prj_vs[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- vss[[tstsetup_vs]]$suggested_size
    }
    pp_tester(pps_vs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_out_expected =
                ndr_ncl_dtls(args_prj_vs[[tstsetup]])$nprjdraws_out,
              info_str = tstsetup)
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(pps_vs[[tstsetup]]))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "pp: `object` of (informal) class \"proj_list\" (based on cv_varsel()) works"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(prjs_cvvs)) {
    tstsetup_cvvs <- args_prj_cvvs[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- cvvss[[tstsetup_cvvs]]$suggested_size
    }
    pp_tester(pps_cvvs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_out_expected =
                ndr_ncl_dtls(args_prj_cvvs[[tstsetup]])$nprjdraws_out,
              info_str = tstsetup)
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(pps_cvvs[[tstsetup]]))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "`object` of (informal) class \"proj_list\" (created manually) works"
), {
  tstsetups <- grep("\\.glm\\.gauss\\..*\\.clust$", names(prjs), value = TRUE)
  stopifnot(length(tstsetups) > 1)
  pp <- proj_predict(prjs[tstsetups], .seed = seed2_tst)
  pp_tester(pp,
            len_expected = length(tstsetups),
            info_str = paste(tstsetups, collapse = ","))
})

test_that(paste(
  "`object` of class \"refmodel\" and passing arguments to project() works"
), {
  tstsetups <- grep("\\.glm\\.gauss.*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    pp_from_refmod <- do.call(proj_predict, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]],
           .seed = seed2_tst),
      excl_nonargs(args_prj_i)
    ))
    pp_from_prj <- pps[[tstsetup]]
    expect_equal(pp_from_refmod, pp_from_prj, info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"stanreg\" or \"brmsfit\" and passing arguments to",
  "project() works"
), {
  tstsetups <- grep("\\.glm\\.gauss.*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    pp_from_fit <- do.call(proj_predict, c(
      list(object = fits[[args_prj_i$tstsetup_fit]],
           .seed = seed2_tst),
      excl_nonargs(args_prj_i)
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
  tstsetups <- head(grep("\\.glm\\.gauss.*\\.default_meth\\.subvec",
                         names(prjs_vs), value = TRUE),
                    1)
  for (tstsetup in tstsetups) {
    args_prj_vs_i <- args_prj_vs[[tstsetup]]
    pp_from_vsel <- do.call(proj_predict, c(
      list(object = vss[[args_prj_vs_i$tstsetup_vsel]],
           .seed = seed2_tst),
      excl_nonargs(args_prj_vs_i)
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
  tstsetups <- head(
    grep("\\.glm\\.gauss.*\\.default_meth\\.default_cvmeth\\.subvec",
         names(prjs_cvvs), value = TRUE),
    1
  )
  for (tstsetup in tstsetups) {
    args_prj_cvvs_i <- args_prj_cvvs[[tstsetup]]
    pp_from_vsel <- do.call(proj_predict, c(
      list(object = cvvss[[args_prj_cvvs_i$tstsetup_vsel]],
           .seed = seed2_tst),
      excl_nonargs(args_prj_cvvs_i)
    ))
    pp_from_prj <- pps_cvvs[[tstsetup]]
    expect_equal(pp_from_vsel, pp_from_prj, info = tstsetup)
  }
})

test_that("`object` not of class \"vsel\" and missing `solution_terms` fails", {
  expect_error(
    proj_predict(1, .seed = seed2_tst),
    paste("^Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.$")
  )
  expect_error(
    proj_predict(fits[[1]], .seed = seed2_tst),
    paste("^Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.$")
  )
  expect_error(
    proj_predict(refmods[[1]], .seed = seed2_tst),
    paste("^Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.$")
  )
  expect_error(
    proj_predict(c(prjs, list(dat)), .seed = seed2_tst),
    paste("Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.")
  )
})

## newdata and nresample_clusters -----------------------------------------

test_that("invalid `newdata` fails", {
  expect_error(
    proj_predict(prjs, newdata = dat[, 1], .seed = seed2_tst),
    "must be a data\\.frame or a matrix"
  )
  stopifnot(length(solterms_x) > 1)
  expect_error(
    proj_predict(prjs[[head(grep("\\.glm\\.gauss.*\\.solterms_x", names(prjs)),
                            1)]],
                 newdata = dat[, head(solterms_x, -1), drop = FALSE],
                 .seed = seed2_tst,
                 solution_terms = solterms_x),
    "^object '.*' not found$"
  )
})

test_that("`newdata` and `nresample_clusters` work (even in edge cases)", {
  for (tstsetup in names(prjs)) {
    for (nobsv_crr in nobsv_tst) {
      if (args_prj[[tstsetup]]$mod_nm == "gamm") {
        # TODO (GAMMs): Fix this.
        next
      }
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

test_that(paste(
  "omitting the response in `newdata` (not possible for `\"brmsfit\"`-based",
  "reference models) doesn't change results"
), {
  for (tstsetup in names(prjs)) {
    if (args_prj[[tstsetup]]$pkg_nm == "brms") {
      next
    }
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
    if (!(args_prj[[tstsetup]]$pkg_nm == "brms" &&
          args_prj[[tstsetup]]$fam_nm == "binom")) {
      # TODO (brms): Fix or document why this doesn't work for "brmsfit"s.
      pp_ones <- proj_predict(prjs[[tstsetup]],
                              newdata = dat_wobs_ones,
                              weightsnew = ~ wobs_col_ones,
                              .seed = seed2_tst)
      pp_tester(pp_ones,
                nprjdraws_out_expected = ndr_ncl$nprjdraws_out,
                info_str = tstsetup)
    }
    pp <- proj_predict(prjs[[tstsetup]],
                       newdata = dat,
                       weightsnew = ~ wobs_col,
                       .seed = seed2_tst)
    pp_tester(pp,
              nprjdraws_out_expected = ndr_ncl$nprjdraws_out,
              info_str = tstsetup)
    if (!(args_prj[[tstsetup]]$pkg_nm == "brms" &&
          args_prj[[tstsetup]]$fam_nm == "binom")) {
      # TODO (brms): Fix or document why this doesn't work for "brmsfit"s.
      ppw <- proj_predict(prjs[[tstsetup]],
                          newdata = dat_wobs_new,
                          weightsnew = ~ wobs_col_new,
                          .seed = seed2_tst)
      pp_tester(ppw,
                nprjdraws_out_expected = ndr_ncl$nprjdraws_out,
                info_str = tstsetup)
    }
    # Weights are only relevant for the binomial() family:
    if (!args_prj[[tstsetup]]$fam_nm %in% c("brnll", "binom")) {
      expect_equal(pp_ones, pp_orig, info = tstsetup)
      expect_equal(pp, pp_orig, info = tstsetup)
      expect_equal(ppw, pp_orig, info = tstsetup)
    } else if (args_prj[[tstsetup]]$fam_nm == "brnll") {
      expect_equal(pp_ones, pp_orig, info = tstsetup)
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        expect_false(isTRUE(all.equal(pp, pp_orig)), info = tstsetup)
        expect_false(isTRUE(all.equal(ppw, pp_orig)), info = tstsetup)
        expect_false(isTRUE(all.equal(ppw, pp)), info = tstsetup)
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        ### TODO: This might in fact be undesired:
        expect_equal(pp, pp_orig, info = tstsetup)
        expect_equal(ppw, pp_orig, info = tstsetup)
        expect_equal(ppw, pp, info = tstsetup)
        ###
      }
    } else if (args_prj[[tstsetup]]$fam_nm == "binom") {
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        ### TODO: This might in fact be undesired:
        expect_equal(pp_ones, pp_orig, info = tstsetup)
        expect_false(isTRUE(all.equal(pp, pp_orig)), info = tstsetup)
        ###
        expect_false(isTRUE(all.equal(ppw, pp_orig)), info = tstsetup)
        expect_false(isTRUE(all.equal(ppw, pp)), info = tstsetup)
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        expect_equal(pp, pp_orig, info = tstsetup)
      }
    }
  }
})

## offsetnew --------------------------------------------------------------

test_that("`offsetnew` works", {
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    pp_orig <- pps[[tstsetup]]
    if (args_prj[[tstsetup]]$pkg_nm != "brms") {
      # TODO (brms): Fix or document why this doesn't work for "brmsfit"s.
      pp_zeros <- proj_predict(prjs[[tstsetup]],
                               newdata = dat_offs_zeros,
                               offsetnew = ~ offs_col_zeros,
                               .seed = seed2_tst)
      pp_tester(pp_zeros,
                nprjdraws_out_expected = ndr_ncl$nprjdraws_out,
                info_str = tstsetup)
    }
    pp <- proj_predict(prjs[[tstsetup]],
                       newdata = dat,
                       offsetnew = ~ offs_col,
                       .seed = seed2_tst)
    pp_tester(pp,
              nprjdraws_out_expected = ndr_ncl$nprjdraws_out,
              info_str = tstsetup)
    if (args_prj[[tstsetup]]$pkg_nm != "brms") {
      # TODO (brms): Fix or document why this doesn't work for "brmsfit"s.
      ppo <- proj_predict(prjs[[tstsetup]],
                          newdata = dat_offs_new,
                          offsetnew = ~ offs_col_new,
                          .seed = seed2_tst)
      pp_tester(ppo,
                nprjdraws_out_expected = ndr_ncl$nprjdraws_out,
                info_str = tstsetup)
    }
    if (grepl("\\.with_offs", tstsetup)) {
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        ### TODO: This might in fact be undesired:
        expect_equal(pp_zeros, pp_orig, info = tstsetup)
        expect_false(isTRUE(all.equal(pp, pp_orig)), info = tstsetup)
        ###
        # For the gaussian() family, we can perform an easy check (because of
        # the identity link):
        if (args_prj[[tstsetup]]$fam_nm == "gauss") {
          ### TODO: This might in fact be undesired (see above):
          expect_equal(t(pp) - dat$offs_col, t(pp_orig), info = tstsetup)
          expect_equal(t(ppo) - dat_offs_new$offs_col_new, t(pp_orig),
                       info = tstsetup)
          ###
        } else {
          expect_false(isTRUE(all.equal(ppo, pp_orig)), info = tstsetup)
          expect_false(isTRUE(all.equal(ppo, pp)), info = tstsetup)
        }
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        expect_equal(pp, pp_orig, info = tstsetup)
      }
    } else {
      if (args_prj[[tstsetup]]$pkg_nm == "rstanarm") {
        expect_equal(pp_zeros, pp_orig, info = tstsetup)
        expect_false(isTRUE(all.equal(pp, pp_orig)), info = tstsetup)
        # For the gaussian() family, we can perform an easy check (because of
        # the identity link):
        if (args_prj[[tstsetup]]$fam_nm == "gauss") {
          expect_equal(t(pp) - dat$offs_col, t(pp_orig), info = tstsetup)
          expect_equal(t(ppo) - dat_offs_new$offs_col_new, t(pp_orig),
                       info = tstsetup)
        } else {
          expect_false(isTRUE(all.equal(ppo, pp_orig)), info = tstsetup)
          expect_false(isTRUE(all.equal(ppo, pp)), info = tstsetup)
        }
      } else if (args_prj[[tstsetup]]$pkg_nm == "brms") {
        ### TODO: This might in fact be undesired:
        expect_equal(pp, pp_orig, info = tstsetup)
        ###
      }
    }
  }
})

## filter_nterms ----------------------------------------------------------

test_that("`filter_nterms` works (for an `object` of class \"projection\")", {
  tstsetups <- grep("\\.glm\\.gauss.*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    nterms_avail_crr <- length(args_prj[[tstsetup]]$solution_terms)
    nterms_unavail_crr <- c(0L, nterms_avail_crr + 130L)
    stopifnot(!nterms_avail_crr %in% nterms_unavail_crr)
    for (filter_nterms_crr in nterms_unavail_crr) {
      expect_error(proj_predict(prjs[[tstsetup]],
                                filter_nterms = filter_nterms_crr,
                                .seed = seed2_tst),
                   "Invalid `filter_nterms`\\.",
                   info = paste(tstsetup, filter_nterms_crr, sep = "__"))
    }
    pp <- proj_predict(prjs[[tstsetup]],
                       filter_nterms = nterms_avail_crr,
                       .seed = seed2_tst)
    pp_orig <- pps[[tstsetup]]
    expect_equal(pp, pp_orig, info = tstsetup)
  }
})

test_that(paste(
  "`filter_nterms` works (for an `object` of (informal) class \"proj_list\")"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.glm\\.gauss.*\\.default_meth\\.full$", names(prjs_vs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    # Unavailable number(s) of terms:
    for (filter_nterms_crr in nterms_unavail) {
      expect_error(proj_predict(prjs_vs[[tstsetup]],
                                filter_nterms = filter_nterms_crr,
                                .seed = seed2_tst),
                   "Invalid `filter_nterms`\\.",
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
      pp_crr <- proj_predict(prjs_vs[[tstsetup]],
                             filter_nterms = filter_nterms_crr,
                             .seed = seed2_tst)
      if (is.null(filter_nterms_crr)) filter_nterms_crr <- 0:nterms_max_tst
      nhits_nterms <- sum(filter_nterms_crr <= nterms_max_tst)
      pp_tester(pp_crr,
                len_expected = nhits_nterms,
                info_str = paste(tstsetup,
                                 paste(filter_nterms_crr, collapse = ","),
                                 sep = "__"))
      if (identical(filter_nterms_crr, 0:nterms_max_tst)) {
        # The special case of all possible numbers of terms:
        pp_orig <- pps_vs[[tstsetup]]
        expect_equal(pp_crr, pp_orig, info = tstsetup)
      }
    }
  }
})

## Single observation, single draw ----------------------------------------

test_that(paste(
  "a single observation and a single draw work (which implicitly tests",
  "this edge case for family$ppd(), too)"
), {
  for (tstsetup in names(prjs)) {
    if (args_prj[[tstsetup]]$mod_nm == "gamm") {
      # TODO (GAMMs): Fix this.
      next
    }
    pp1 <- proj_predict(refmods[[args_prj[[tstsetup]]$tstsetup_ref]],
                        newdata = head(get_dat(tstsetup), 1),
                        nresample_clusters = 1L,
                        .seed = seed2_tst,
                        solution_terms = args_prj[[tstsetup]]$solution_terms,
                        nclusters = 1L,
                        seed = seed_tst)
    pp_tester(pp1,
              nprjdraws_out_expected = 1L,
              nobsv_expected = 1L,
              info_str = tstsetup)
  }
})
