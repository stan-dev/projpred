context("miscellaneous")

# get_refdist() ----------------------------------------------------------

## ndraws and nclusters ---------------------------------------------------

test_that(paste(
  "get_refdist(): `ndraws = NULL` and `nclusters = NULL` leads to",
  "`ndraws = S` (and `nclusters = NULL`)"
), {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  set.seed(seed2_tst)
  for (tstsetup in names(refmods)) {
    refdist <- get_refdist(refmods[[tstsetup]])
    if (args_ref[[tstsetup]]$prj_nm == "augdat") {
      nobsv_fac <- length(refmods[[tstsetup]]$family$cats)
    } else {
      nobsv_fac <- 1L
    }
    # The following refdist_tester() call runs more expectations than necessary
    # for this test (only the one for `refdist$clust_used` and the dim() test
    # for `refdist$mu` are actually necessary):
    refdist_tester(
      refdist,
      nobsv_expected = nobsv * nobsv_fac,
      nprjdraws_expected = nrefdraws,
      clust_expected = FALSE,
      fam_expected = refmods[[tstsetup]]$family$family,
      info_str = tstsetup
    )
  }
  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})

test_that(paste(
  "`ndraws` and/or `nclusters` too big causes them to be cut off at the number",
  "of posterior draws in the reference model"
), {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  set.seed(seed2_tst)
  for (tstsetup in names(refmods)) {
    for (ndraws_crr in list(nrefdraws + 1L)) {
      for (nclusters_crr in list(NULL, nrefdraws + 1L)) {
        refdist <- get_refdist(refmods[[tstsetup]],
                               ndraws = ndraws_crr,
                               nclusters = nclusters_crr)
        if (args_ref[[tstsetup]]$prj_nm == "augdat") {
          nobsv_fac <- length(refmods[[tstsetup]]$family$cats)
        } else {
          nobsv_fac <- 1L
        }
        refdist_tester(
          refdist,
          nobsv_expected = nobsv * nobsv_fac,
          nprjdraws_expected = nrefdraws,
          clust_expected = FALSE,
          fam_expected = refmods[[tstsetup]]$family$family,
          info_str = paste(tstsetup, ndraws_crr, nclusters_crr, sep = "__")
        )
      }
    }
  }
  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})

# force_search_terms() ----------------------------------------------------

test_that("force_search_terms() works", {
  expect_identical(
    force_search_terms(
      forced_terms = paste0("X", 1:2),
      optional_terms = paste0("X", 3:5)
    ),
    c("X1 + X2", "X1 + X2 + X3", "X1 + X2 + X4", "X1 + X2 + X5"),
    info = "two forced, three optional terms"
  )
  expect_identical(
    force_search_terms(
      forced_terms = paste0("X", 1),
      optional_terms = paste0("X", 3:5)
    ),
    c("X1", "X1 + X3", "X1 + X4", "X1 + X5"),
    info = "one forced, three optional terms"
  )
  expect_identical(
    force_search_terms(
      forced_terms = paste0("X", 1:2),
      optional_terms = paste0("X", 3)
    ),
    c("X1 + X2", "X1 + X2 + X3"),
    info = "two forced, one optional term"
  )
  expect_identical(
    force_search_terms(
      forced_terms = paste0("X", 1),
      optional_terms = paste0("X", 3)
    ),
    c("X1", "X1 + X3"),
    info = "one forced, one optional term"
  )
  expect_error(
    force_search_terms(
      forced_terms = character(),
      optional_terms = paste0("X", 3)
    ),
    "length\\(forced_terms\\) > 0 is not TRUE",
    info = "zero forced, one optional term"
  )
  expect_error(
    force_search_terms(
      forced_terms = paste0("X", 1),
      optional_terms = character()
    ),
    "length\\(optional_terms\\) > 0 is not TRUE",
    info = "one forced, zero optional terms"
  )
})

# rstanarm: special formulas ----------------------------------------------

test_that("rstanarm: special formulas work", {
  # Skip this on CRAN to avoid depending too strongly on rstanarm internals:
  skip_on_cran()
  # Note: This test only tests that rstanarm handles special formulas correctly.
  # Within projpred, arithmetic expressions on the left-hand side of a formula
  # are handled by get_refmodel() and init_refmodel(); arithmetic expressions on
  # the right-hand side of a formula are handled by the `div_minimizer`.
  tstsetups <- grep("^rstanarm.*\\.spclformul", names(fits), value = TRUE)
  # Compare the "special formula" fit with the corresponding "standard formula"
  # fit (which does not have the special formula elements):
  for (tstsetup in tstsetups) {
    mf_spclformul <- fits[[tstsetup]]$model
    if (grepl("\\.glmm\\.", tstsetup)) {
      expect_null(mf_spclformul, info = tstsetup)
      mf_spclformul <- fits[[tstsetup]]$glmod$fr
    } else {
      expect_false(is.null(mf_spclformul), info = tstsetup)
    }
    nms_spclformul <- grep("y_|xco", names(mf_spclformul), value = TRUE)

    tstsetup_stdformul <- sub("\\.spclformul", ".stdformul", tstsetup)
    stopifnot(tstsetup_stdformul != tstsetup)
    if (tstsetup_stdformul %in% names(fits)) {
      mf_stdformul <- fits[[tstsetup_stdformul]]$model
      if (grepl("\\.glmm\\.", tstsetup_stdformul)) {
        expect_null(mf_stdformul, info = tstsetup_stdformul)
        mf_stdformul <- fits[[tstsetup_stdformul]]$glmod$fr
      } else {
        expect_false(is.null(mf_stdformul), info = tstsetup_stdformul)
      }
      nms_stdformul <- grep("y_|xco", names(mf_stdformul), value = TRUE)
      expect_equal(mf_spclformul[, setdiff(names(mf_spclformul),
                                           nms_spclformul)],
                   mf_stdformul[, setdiff(names(mf_stdformul),
                                          nms_stdformul)],
                   info = tstsetup)
    }
    # Check arithmetic expressions:
    for (nm_spclformul in nms_spclformul) {
      expect_equal(mf_spclformul[, nm_spclformul],
                   eval(str2lang(nm_spclformul), envir = dat),
                   info = paste(tstsetup, nm_spclformul, sep = "__"))
    }
  }
})

# Other internal functions ------------------------------------------------

test_that(paste(
  "`pseudo_data(f = 0, [...], family = extend_family(gaussian()), [...])` is",
  "essentially an identity function (apart from offsets)."
), {
  mu_crr <- matrix(1:3, ncol = 1)
  wobs_crr <- c(2.5, 6, 4)
  offs_crr <- c(-4.2, 3.5, 1.1)
  psdat <- pseudo_data(
    f = 0,
    y = mu_crr,
    family = extend_family(gaussian()),
    weights = wobs_crr,
    offset = offs_crr
  )
  expect_equal(psdat$z, mu_crr - offs_crr)
  expect_false(isTRUE(all.equal(psdat$z, mu_crr)))
  expect_equal(psdat$wobs, wobs_crr)
})

test_that(".srs_diff_est_w() propagates input `NA`s to its output", {
  nloo_tst <- nobsv %/% 5L
  loo_inds_tst <- ceiling(seq(1L, nobsv, length.out = nloo_tst))
  skip_if_not(length(unique(loo_inds_tst)) == nloo_tst)
  srs_diff_res_NAs <- list(
    m = nloo_tst,
    N = nobsv,
    y_hat = NA_real_,
    v_y_hat = NA_real_,
    hat_v_y = NA_real_
  )
  expect_identical(
    .srs_diff_est_w(y_approx = rep(NA, nobsv),
                    y = rep(NA, nloo_tst),
                    y_idx = loo_inds_tst,
                    wobs = rep(NA, nobsv)),
    srs_diff_res_NAs,
    info = "all inputs (except `y_idx`) are all-`NA`s"
  )
  expect_identical(
    .srs_diff_est_w(y_approx = rep(-0.8, nobsv),
                    y = rep(NA, nloo_tst),
                    y_idx = loo_inds_tst,
                    wobs = rep_len(c(2, 3), length.out = nobsv)),
    srs_diff_res_NAs,
    info = "`y_approx` without `NA`s, `y` only `NA`s, `wobs` without `NA`s"
  )
  expect_identical(
    .srs_diff_est_w(y_approx = rep(-0.8, nobsv),
                    y = c(rep(-0.7, nloo_tst - 1L), NA),
                    y_idx = loo_inds_tst,
                    wobs = rep_len(c(2, 3), length.out = nobsv)),
    srs_diff_res_NAs,
    info = paste0("`y_approx` without `NA`s, `y` with a single `NA`, ",
                  "`wobs` without `NA`s")
  )
  expect_identical(
    .srs_diff_est_w(y_approx = c(NA, rep(-0.8, nobsv - 1L)),
                    y = c(rep(NA, nloo_tst - 1L), -0.7),
                    y_idx = loo_inds_tst,
                    wobs = rep_len(c(2, 3), length.out = nobsv)),
    srs_diff_res_NAs,
    info = paste0("`y_approx` with a single `NA`, `y` with a single `NA`, ",
                  "`wobs` without `NA`s")
  )
})

test_that(".weighted_sd() with `na.rm = FALSE` propagates input `NA`s to its output", {
  expect_identical(
    .weighted_sd(x = rep(NA, nobsv),
                 w = rep(NA, nobsv)),
    NA_real_,
    info = "all inputs are all-`NA`s"
  )
  expect_identical(
    .weighted_sd(x = rep(-0.8, nobsv),
                 w = c(rep_len(c(2, 3), length.out = nobsv - 1L), NA)),
    NA_real_,
    info = "`x` without `NA`s, `w` with a single `NA`"
  )
  expect_identical(
    .weighted_sd(x = c(rep(-0.8, nobsv - 1L), NA),
                 w = rep_len(c(2, 3), length.out = nobsv)),
    NA_real_,
    info = "`x` with a single `NA`, `w` without `NA`s"
  )
  expect_identical(
    .weighted_sd(x = rep(-0.8, nobsv),
                 w = rep(NA, nobsv)),
    NA_real_,
    info = "`x` without `NA`s, `w` all-`NA`s"
  )
  expect_identical(
    .weighted_sd(x = rep(NA, nobsv),
                 w = rep_len(c(2, 3), length.out = nobsv)),
    NA_real_,
    info = "`x` all-`NA`s, `w` without `NA`s"
  )
})

test_that(".auc() works", {
  nobsv_auc <- 19L
  expect_equal(
    .auc(cbind(rep_len(c(0, 1), length.out = nobsv_auc),
               imperfect_alternation(c(0.3, 0.8), length.out = nobsv_auc),
               rep(1, nobsv_auc))),
    0.6833333333333333333333,
    info = "`wobs` column only `1`s"
  )
  expect_equal(
    .auc(cbind(rep_len(c(0, 1), length.out = nobsv_auc),
               imperfect_alternation(c(0.3, 0.8), length.out = nobsv_auc),
               imperfect_alternation(c(1, 2), n_tail = 11L, length.out = nobsv_auc))),
    0.717948717948718062587,
    info = "`wobs` column with imperfect alternation of `1`s and `2`s"
  )
})

test_that(".auc() propagates input `NA`s to its output", {
  nobsv_auc <- 19L
  expect_identical(
    .auc(cbind(rep(NA, nobsv_auc),
               rep(NA, nobsv_auc),
               rep(NA, nobsv_auc))),
    NA_real_,
    info = "all inputs are all-`NA`s"
  )
  expect_identical(
    .auc(cbind(rep(1, nobsv_auc),
               rep(NA, nobsv_auc),
               rep(1, nobsv_auc))),
    NA_real_,
    info = paste0("`resp` column without `NA`s, ",
                  "`pred` column only `NA`s, ",
                  "`wobs` column without `NA`s")
  )
  expect_identical(
    .auc(cbind(rep(1, nobsv_auc),
               rep(NA, nobsv_auc),
               c(rep(1, nobsv_auc - 1L), NA))),
    NA_real_,
    info = paste0("`resp` column without `NA`s, ",
                  "`pred` column only `NA`s, ",
                  "`wobs` column with a single `NA`")
  )
  expect_identical(
    .auc(cbind(c(rep_len(c(0, 1), length.out = nobsv_auc - 1L), NA),
               rep(0.7, nobsv_auc),
               rep(1, nobsv_auc))),
    NA_real_,
    info = paste0("`resp` column with a single `NA`, ",
                  "`pred` column without `NA`s, ",
                  "`wobs` column without `NA`s")
  )
  expect_identical(
    .auc(cbind(rep(NA, nobsv_auc),
               rep(0.7, nobsv_auc),
               rep(1, nobsv_auc))),
    NA_real_,
    info = paste0("`resp` column only `NA`s, ",
                  "`pred` column without `NA`s, ",
                  "`wobs` column without `NA`s")
  )
  expect_identical(
    .auc(cbind(c(rep_len(c(0, 1), length.out = nobsv_auc - 1L), NA),
               c(rep(0.7, nobsv_auc - 1L), NA),
               rep(1, nobsv_auc))),
    NA_real_,
    info = paste0("`resp` column with a single `NA`, ",
                  "`pred` column with a single `NA`, ",
                  "`wobs` column without `NA`s")
  )
  expect_identical(
    .auc(cbind(c(rep_len(c(0, 1), length.out = nobsv_auc - 1L), NA),
               c(NA, rep(0.7, nobsv_auc - 1L)),
               rep(1, nobsv_auc))),
    NA_real_,
    info = paste0("`resp` column with a single `NA`, ",
                  "`pred` column with a single `NA` (at different position), ",
                  "`wobs` column without `NA`s")
  )
  expect_identical(
    .auc(cbind(c(rep_len(c(0, 1), length.out = nobsv_auc - 1L), NA),
               rep(NA, nobsv_auc),
               rep(1, nobsv_auc))),
    NA_real_,
    info = paste0("`resp` column with a single `NA`, ",
                  "`pred` column only `NA`s, ",
                  "`wobs` column without `NA`s")
  )
  expect_identical(
    .auc(cbind(rep(NA, nobsv_auc),
               rep(0.7, nobsv_auc),
               c(rep(1, nobsv_auc - 1L), NA))),
    NA_real_,
    info = paste0("`resp` column only `NA`s, ",
                  "`pred` column without `NA`s, ",
                  "`wobs` column with a single `NA`")
  )
  expect_identical(
    .auc(cbind(rep_len(c(0, 1), length.out = nobsv_auc),
               imperfect_alternation(c(0.3, 0.8), length.out = nobsv_auc),
               rep(NA, nobsv_auc))),
    NA_real_,
    info = paste0("`resp` column without `NA`s, ",
                  "`pred` column without `NA`s, ",
                  "`wobs` column only `NA`s")
  )
})
