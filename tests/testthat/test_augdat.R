context("augmented-data projection")

# Setup -------------------------------------------------------------------

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

# augdat.R ----------------------------------------------------------------

nobs_orig_tst <- 4L
ncat_tst <- 2L
ndraws_tst <- 5L
xtst <- seq_len(nobs_orig_tst * ncat_tst * ndraws_tst)
arrtst <- array(xtst, dim = c(nobs_orig_tst, ncat_tst, ndraws_tst))
augmtst <- arr2augmat(arrtst)

## arr2augmat(), augmat2arr() ---------------------------------------------

test_that("arr2augmat(), augmat2arr()", {
  expect_identical(augmtst, structure(matrix(xtst, ncol = ndraws_tst),
                                      nobs_orig = nobs_orig_tst,
                                      class = "augmat"))

  # With "margin_draws = 3":
  xtst_m3 <- 0.1
  stopifnot(length(xtst_m3) == 1)
  arr_m3 <- array(xtst_m3, dim = c(1, 1, 1))
  augmat_m3 <- arr2augmat(arr_m3)
  augmat_m3_ch <- structure(matrix(xtst_m3),
                            nobs_orig = 1L,
                            class = "augmat")
  expect_identical(augmat_m3, augmat_m3_ch)
  arr_m3_ch <- augmat2arr(augmat_m3)
  expect_identical(arr_m3_ch, arr_m3)

  # With "margin_draws = 1":
  xtst_m1 <- seq(-1, 1, by = 0.1)
  stopifnot(length(xtst_m1) == 21)
  arr_m1 <- array(xtst_m1, dim = c(7, 3, 1))
  augmat_m1 <- arr2augmat(arr_m1, margin_draws = 1)
  augmat_m1_ch <- structure(matrix(xtst_m1, ncol = 7, byrow = TRUE),
                            nobs_orig = 3L,
                            class = "augmat")
  expect_identical(augmat_m1, augmat_m1_ch)
  arr_m1_ch <- augmat2arr(augmat_m1, margin_draws = 1)
  expect_identical(arr_m1_ch, arr_m1)
})

## t.augmat() -------------------------------------------------------------

test_that("t.augmat()", {
  expect_identical(t(augmtst), matrix(xtst, nrow = ndraws_tst, byrow = TRUE))
})

## `[.augmat`() -----------------------------------------------------------

test_that("`[.augmat`()", {
  expect_identical(augmtst[], augmtst)
  # Subsetting columns:
  expect_identical(augmtst[, 1:2],
                   structure(matrix(head(xtst, nobs_orig_tst * ncat_tst * 2),
                                    ncol = 2),
                             nobs_orig = nobs_orig_tst,
                             class = "augmat"))
  expect_identical(augmtst[, 1, drop = FALSE],
                   structure(matrix(head(xtst, nobs_orig_tst * ncat_tst)),
                             nobs_orig = nobs_orig_tst,
                             class = "augmat"))
  expect_identical(augmtst[, 1],
                   structure(head(xtst, nobs_orig_tst * ncat_tst),
                             nobs_orig = nobs_orig_tst,
                             class = "augvec"))
  # Subsetting rows (in fact, this is not of interest, at least currently):
  xrow1 <- xtst[nobs_orig_tst * ncat_tst * (seq_len(ndraws_tst) - 1L) + 1L]
  expect_identical(augmtst[1, , drop = FALSE],
                   structure(t(xrow1),
                             nobs_orig = nobs_orig_tst,
                             class = "augmat"))
  expect_identical(augmtst[1, ],
                   structure(xrow1,
                             nobs_orig = nobs_orig_tst,
                             class = "augvec"))
  # Subsetting rows and columns:
  expect_identical(augmtst[1:2, 1:2],
                   structure(arrtst[1:2, 1, 1:2],
                             nobs_orig = nobs_orig_tst,
                             class = "augmat"))
  expect_identical(augmtst[1, 1],
                   structure(head(xtst, 1),
                             nobs_orig = nobs_orig_tst,
                             class = "augvec"))

  # Assigning:
  augmtst[1, 1] <- 0
  expect_identical(augmtst, structure(matrix(c(0, xtst[-1]), ncol = ndraws_tst),
                                      nobs_orig = nobs_orig_tst,
                                      class = "augmat"))
})

## augmat2augvec(), augvec2augmat() ---------------------------------------

augmtst_1col <- augmtst[, 1, drop = FALSE]
augvtst <- augmat2augvec(augmtst_1col)

test_that("augmat2augvec(), augvec2augmat()", {
  # The structure of `augmtst[, 1]` has already been tested above, so simply use
  # `augmtst[, 1]` as expectation here:
  expect_identical(augvtst, augmtst[, 1])
  expect_identical(augvec2augmat(augvtst), augmtst_1col)
})

## t.augvec() -------------------------------------------------------------

test_that("t.augvec()", {
  expect_identical(t(augvtst),
                   matrix(head(xtst, nobs_orig_tst * ncat_tst), nrow = 1))
})

## `[.augvec`() -----------------------------------------------------------

test_that("`[.augvec`()", {
  expect_identical(augvtst[], augvtst)
  expect_identical(augvtst[1:2],
                   structure(head(xtst, 2),
                             nobs_orig = nobs_orig_tst,
                             class = "augvec"))
  expect_identical(augvtst[1],
                   structure(head(xtst, 1),
                             nobs_orig = nobs_orig_tst,
                             class = "augvec"))
  expect_identical(augvtst[integer()],
                   structure(head(xtst, 0),
                             nobs_orig = nobs_orig_tst,
                             class = "augvec"))

  # Assigning:
  augvtst[1] <- 0
  expect_identical(augvtst,
                   structure(c(0, head(xtst, nobs_orig_tst * ncat_tst)[-1]),
                             nobs_orig = nobs_orig_tst,
                             class = "augvec"))
})

## catmaxprb() ------------------------------------------------------------

test_that("catmaxprb()", {
  stopifnot(ncat_tst == 2)
  expect_identical(catmaxprb(augvtst, lvls = c("0", "1")),
                   factor(rep("1", nobs_orig_tst), levels = c("0", "1")))

  # With probabilities (and a single observation):
  augvtst_pr <- structure(c(0.7, 0.3),
                          nobs_orig = 1L,
                          class = "augvec")
  expect_identical(catmaxprb(augvtst_pr, lvls = c("lvl1", "lvl2")),
                   factor("lvl1", levels = c("lvl1", "lvl2")))
})

## augdat_link_binom(), augdat_ilink_binom() ------------------------------

test_that("augdat_link_binom(), augdat_ilink_binom()", {
  arr_pr <- array(c(0.7, 0.3), dim = c(1, 1, 2))
  arr_lat <- augdat_link_binom(arr_pr)
  expect_equal(arr_lat, array(log(0.3 / 0.7), dim = c(1, 1, 1)))
  expect_equal(augdat_ilink_binom(arr_lat), arr_pr)
})

# Comparison with traditional projection ----------------------------------

## Clustering -------------------------------------------------------------

test_that(paste(
  "clustering `<refmodel_object>$mu` in case of the `brnll` family is the same",
  "no matter whether `<refmodel_object>` was set up for the augmented-data or",
  "the traditional projection"
), {
  skip_if_not(run_prj)
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  tstsetups <- grep("\\.brnll\\..*\\.augdat\\.", names(prjs), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(prjs)) next

    args_prj_i <- args_prj[[tstsetup]]
    args_prj_i_trad <- args_prj[[tstsetup_trad]]
    refmod_crr <- refmods[[args_prj_i$tstsetup_ref]]
    refmod_crr_trad <- refmods[[args_prj_i_trad$tstsetup_ref]]
    set.seed(args_prj_i$seed)
    pref_aug <- .get_refdist(refmod_crr,
                             ndraws = args_prj_i$ndraws,
                             nclusters = args_prj_i$nclusters)
    set.seed(args_prj_i_trad$seed)
    pref_trad <- .get_refdist(refmod_crr_trad,
                              ndraws = args_prj_i_trad$ndraws,
                              nclusters = args_prj_i_trad$nclusters)

    eta_aug <- refmod_crr$family$linkfun(pref_aug$mu)
    eta_trad <- refmod_crr_trad$family$linkfun(pref_trad$mu)
    expect_identical(structure(unclass(eta_aug), nobs_orig = NULL), eta_trad,
                     info = tstsetup)
  }
  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})

## Projection -------------------------------------------------------------

test_that(paste(
  "augmented-data and traditional projection give similar results (in case of",
  "the `brnll` family and the default `regul` for submodels fitted by",
  "fit_glm_ridge_callback())"
), {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.brnll\\..*\\.augdat\\.", names(prjs), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(prjs)) next

    prjmat <- suppressWarnings(as.matrix(prjs[[tstsetup]]))
    prjmat_trad <- suppressWarnings(as.matrix(prjs[[tstsetup_trad]]))

    tol_coefs <- ifelse(
      args_prj[[tstsetup]]$mod_nm == "glmm" &&
        any(grepl("^\\(.*\\)$", args_prj[[tstsetup]]$solution_terms)),
      1e-3, 1e-6
    )
    expect_equal(prjmat, prjmat_trad, tolerance = tol_coefs, info = tstsetup)
  }
})

test_that(paste(
  "augmented-data and traditional projection give (almost) the same results",
  "for submodels fitted by fit_glm_ridge_callback() if `regul = 0` is used (in",
  "case of the `brnll` family)"
), {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.glm\\.brnll\\..*\\.augdat\\.", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(prjs)) next

    args_prj_i <- args_prj[[tstsetup]]
    args_prj_i_trad <- args_prj[[tstsetup_trad]]
    prj <- do.call(project, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]], regul = 0),
      excl_nonargs(args_prj_i)
    ))
    prj_trad <- do.call(project, c(
      list(object = refmods[[args_prj_i_trad$tstsetup_ref]], regul = 0),
      excl_nonargs(args_prj_i_trad)
    ))
    prjmat <- suppressWarnings(as.matrix(prj))
    prjmat_trad <- suppressWarnings(as.matrix(prj_trad))

    tol_coefs <- ifelse(ndr_ncl_dtls(args_prj_i)$clust_used,
                        1e-9, 1e-14)
    expect_equal(prjmat, prjmat_trad, tolerance = tol_coefs, info = tstsetup)
  }
})

## Prediction -------------------------------------------------------------

test_that(paste(
  "proj_linpred() gives similar results for the augmented-data and the",
  "traditional projection (in case of the `brnll` family and the default",
  "`regul` for submodels fitted by fit_glm_ridge_callback())"
), {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.brnll\\..*\\.augdat\\.", names(pls), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(pls)) next

    pl <- pls[[tstsetup]]
    expect_length(dim(pl$pred), 3)
    expect_identical(dim(pl$pred)[3], 1L)
    pl$pred <- matrix(pl$pred,
                      nrow = dim(pl$pred)[1],
                      ncol = dim(pl$pred)[2])
    dimnames(pl$pred) <- list(NULL, as.character(seq_len(ncol(pl$pred))))
    dimnames(pl$lpd) <- list(NULL, as.character(seq_len(ncol(pl$lpd))))
    pl_trad <- pls[[tstsetup_trad]]

    tol_lpreds <- ifelse(
      args_prj[[tstsetup]]$mod_nm == "glmm" &&
        any(grepl("^\\(.*\\)$", args_prj[[tstsetup]]$solution_terms)),
      1e-4, 1e-6
    )
    expect_equal(pl, pl_trad, tolerance = tol_lpreds, info = tstsetup)
  }
})

## Variable selection -----------------------------------------------------

test_that(paste(
  "varsel() gives similar results for the augmented-data and the traditional",
  "projection (in case of the `brnll` family and the default `regul` for",
  "submodels fitted by fit_glm_ridge_callback())"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.brnll\\..*\\.augdat\\.", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.default_meth", ".trad_compare.forward",
                         tstsetup)
    if (!tstsetup_trad %in% names(vss)) next

    vs <- vss[[tstsetup]]
    vs_trad <- vss[[tstsetup_trad]]

    summs_sub <- vs$summaries$sub
    summs_sub_mu <- do.call(cbind, lapply(summs_sub, "[[", "mu"))[
      (nobsv + 1):(2 * nobsv), , drop = FALSE
    ]
    summs_sub_trad <- vs_trad$summaries$sub
    summs_sub_mu_trad <- do.call(cbind, lapply(summs_sub_trad, "[[", "mu"))
    expect_equal(summs_sub_mu, summs_sub_mu_trad, tolerance = 1e-6,
                 info = tstsetup)
    summs_sub_lppd <- do.call(cbind, lapply(summs_sub, "[[", "lppd"))
    summs_sub_lppd_trad <- unname(
      do.call(cbind, lapply(summs_sub_trad, "[[", "lppd"))
    )
    expect_equal(summs_sub_lppd, summs_sub_lppd_trad, tolerance = 1e-6,
                 info = tstsetup)
    summs_ref <- vs$summaries$ref
    summs_ref$mu <- structure(unclass(summs_ref$mu), nobs_orig = NULL)
    summs_ref$mu <- summs_ref$mu[(nobsv + 1):(2 * nobsv)]
    tol_ref <- .Machine$double.eps
    if (args_vs[[tstsetup]]$mod_nm == "glmm") {
      tol_ref <- 1e3 * tol_ref
    }
    expect_equal(summs_ref, vs_trad$summaries$ref, tolerance = tol_ref,
                 info = tstsetup)
  }
})

test_that(paste(
  "summary.vsel() applied to varsel() output gives similar results for the",
  "augmented-data and the traditional projection (in case of the `brnll`",
  "family and the default `regul` for submodels fitted by",
  "fit_glm_ridge_callback())"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.brnll\\..*\\.augdat\\.", names(smmrys_vs), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.default_meth", ".trad_compare.forward",
                         tstsetup)
    if (!tstsetup_trad %in% names(smmrys_vs)) next

    smmry_vs <- smmrys_vs[[tstsetup]]
    smmry_vs_trad <- smmrys_vs[[tstsetup_trad]]

    expect_equal(
      smmry_vs[setdiff(names(smmry_vs), c("family", "selection"))],
      smmry_vs_trad[setdiff(names(smmry_vs_trad), c("family", "selection"))],
      info = tstsetup
    )
    expect_identical(
      smmry_vs$selection[, c("size", "solution_terms")],
      smmry_vs_trad$selection[, c("size", "solution_terms")],
      info = tstsetup
    )
    expect_equal(smmry_vs$selection, smmry_vs_trad$selection, tolerance = 1e-6,
                 info = tstsetup)
  }
})

# The `stats = "acc"` case is excluded in the test above, so test it separately
# here:
test_that(paste(
  "summary.vsel() applied to varsel() output gives similar results for the",
  "augmented-data and the traditional projection when using `stats = \"acc\"`",
  "(in case of the `brnll` family and the default `regul` for submodels fitted",
  "by fit_glm_ridge_callback())"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.brnll\\..*\\.augdat_stats\\.", names(smmrys_vs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.default_meth", ".trad_compare.forward",
                         tstsetup)
    tstsetup_trad <- sub("\\.augdat_stats\\.", ".binom_stats.", tstsetup_trad)
    if (!tstsetup_trad %in% names(smmrys_vs)) next

    smmry_vs <- smmrys_vs[[tstsetup]]
    smmry_vs_trad <- smmrys_vs[[tstsetup_trad]]

    expect_equal(
      smmry_vs[setdiff(names(smmry_vs), c("family", "selection"))],
      smmry_vs_trad[setdiff(names(smmry_vs_trad),
                            c("family", "selection"))],
      info = tstsetup
    )
    expect_identical(
      smmry_vs$selection[, c("size", "solution_terms")],
      smmry_vs_trad$selection[, c("size", "solution_terms")],
      info = tstsetup
    )
    # Exclude statistics which are not supported for the augmented-data
    # projection:
    smmry_vs_trad$selection <- smmry_vs_trad$selection[
      , -grep("mse|auc", names(smmry_vs_trad$selection)), drop = FALSE
    ]
    expect_equal(smmry_vs$selection, smmry_vs_trad$selection,
                 tolerance = 1e-6, info = tstsetup)
  }
})

test_that(paste(
  "cv_varsel() gives similar results for the augmented-data and the",
  "traditional projection (in case of the `brnll` family and the default",
  "`regul` for submodels fitted by fit_glm_ridge_callback())"
), {
  skip_if_not(run_cvvs)
  tstsetups <- grep("\\.brnll\\..*\\.augdat\\.", names(cvvss), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.default_meth", ".trad_compare.forward",
                         tstsetup)
    if (!tstsetup_trad %in% names(cvvss)) next

    cvvs <- cvvss[[tstsetup]]
    cvvs_trad <- cvvss[[tstsetup_trad]]

    is_kfold <- identical(args_cvvs[[tstsetup]]$cv_method, "kfold")
    if (is_kfold) {
      is_glmm <- args_cvvs[[tstsetup]]$mod_nm == "glmm"
      if (is_glmm) {
        tol_sub <- 1e-4
      } else {
        tol_sub <- 1e-5
      }
    } else {
      tol_sub <- 1e-6
    }

    summs_sub <- cvvs$summaries$sub
    summs_sub_mu <- do.call(cbind, lapply(summs_sub, "[[", "mu"))[
      (nobsv + 1):(2 * nobsv), , drop = FALSE
    ]
    summs_sub_trad <- cvvs_trad$summaries$sub
    summs_sub_mu_trad <- do.call(cbind, lapply(summs_sub_trad, "[[", "mu"))
    expect_equal(summs_sub_mu, summs_sub_mu_trad, tolerance = tol_sub,
                 info = tstsetup)
    summs_sub_lppd <- do.call(cbind, lapply(summs_sub, "[[", "lppd"))
    summs_sub_lppd_trad <- unname(
      do.call(cbind, lapply(summs_sub_trad, "[[", "lppd"))
    )
    # Sometimes, there are seemingly large differences on log scale which are
    # probably due to numerical overflow of probabilities (i.e., on exp scale)
    # towards 1 or underflow towards zero. Thus, compare on exp scale if the
    # comparison on log scale fails:
    tryCatch(expect_equal(summs_sub_lppd, summs_sub_lppd_trad,
                          tolerance = tol_sub, info = tstsetup),
             error = function(e) {
               stopifnot(is_kfold)
               expect_equal(exp(summs_sub_lppd), exp(summs_sub_lppd_trad),
                            tolerance = tol_sub, info = tstsetup)
             })
    summs_ref <- cvvs$summaries$ref
    summs_ref$mu <- structure(unclass(summs_ref$mu), nobs_orig = NULL)
    summs_ref$mu <- summs_ref$mu[(nobsv + 1):(2 * nobsv)]
    expect_equal(summs_ref, cvvs_trad$summaries$ref, tolerance = 1e-12,
                 info = tstsetup)
  }
})

test_that(paste(
  "summary.vsel() applied to cv_varsel() output gives similar results for the",
  "augmented-data and the traditional projection (in case of the `brnll`",
  "family and the default `regul` for submodels fitted by",
  "fit_glm_ridge_callback())"
), {
  skip_if_not(run_cvvs)
  tstsetups <- grep("\\.brnll\\..*\\.augdat\\.", names(smmrys_cvvs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.default_meth", ".trad_compare.forward",
                         tstsetup)
    if (!tstsetup_trad %in% names(smmrys_cvvs)) next

    smmry_cvvs <- smmrys_cvvs[[tstsetup]]
    smmry_cvvs_trad <- smmrys_cvvs[[tstsetup_trad]]

    expect_equal(
      smmry_cvvs[setdiff(names(smmry_cvvs), c("family", "selection"))],
      smmry_cvvs_trad[setdiff(names(smmry_cvvs_trad),
                              c("family", "selection"))],
      info = tstsetup
    )
    expect_identical(
      smmry_cvvs$selection[, c("size", "solution_terms")],
      smmry_cvvs_trad$selection[, c("size", "solution_terms")],
      info = tstsetup
    )
    is_kfold <- identical(
      args_cvvs[[args_smmry_cvvs[[tstsetup]]$tstsetup_vsel]]$cv_method,
      "kfold"
    )
    if (is_kfold) {
      tol_smmry <- 1e-5
    } else {
      tol_smmry <- 1e-6
    }
    compare_exp <- function(e) {
      stopifnot(is_kfold)
      # Check that we have the default `stats` in this case, meaning only the
      # ELPD:
      stopifnot(is.null(args_smmry_cvvs[[tstsetup]]$stats))

      smmry_pd <- smmry_cvvs$selection
      smmry_pd[[grep("elpd", names(smmry_pd), value = TRUE)]] <- exp(
        smmry_pd[[grep("elpd", names(smmry_pd), value = TRUE)]]
      )
      smmry_pd$lower <- exp(smmry_pd$lower)
      smmry_pd$upper <- exp(smmry_pd$upper)

      smmry_pd_trad <- smmry_cvvs_trad$selection
      smmry_pd_trad[[grep("elpd", names(smmry_pd_trad), value = TRUE)]] <- exp(
        smmry_pd_trad[[grep("elpd", names(smmry_pd_trad), value = TRUE)]]
      )
      smmry_pd_trad$lower <- exp(smmry_pd_trad$lower)
      smmry_pd_trad$upper <- exp(smmry_pd_trad$upper)

      expect_equal(smmry_pd[, setdiff(names(smmry_pd), "se")],
                   smmry_pd_trad[, setdiff(names(smmry_pd_trad), "se")],
                   tolerance = 1e-10, info = tstsetup)
      expect_equal(smmry_pd$se, smmry_pd_trad$se, tolerance = 1e-1,
                   info = tstsetup)
    }
    tryCatch(expect_equal(smmry_cvvs$selection, smmry_cvvs_trad$selection,
                          tolerance = tol_smmry, info = tstsetup),
             error = compare_exp)
  }
})

# The `stats = "acc"` case is excluded in the test above, so test it separately
# here:
test_that(paste(
  "summary.vsel() applied to cv_varsel() output gives similar results for the",
  "augmented-data and the traditional projection when using `stats = \"acc\"`",
  "(in case of the `brnll` family and the default `regul` for submodels fitted",
  "by fit_glm_ridge_callback())"
), {
  skip_if_not(run_cvvs)
  tstsetups <- grep("\\.brnll\\..*\\.augdat_stats\\.", names(smmrys_cvvs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.augdat\\.default_meth", ".trad_compare.forward",
                         tstsetup)
    tstsetup_trad <- sub("\\.augdat_stats\\.", ".binom_stats.", tstsetup_trad)
    if (!tstsetup_trad %in% names(smmrys_cvvs)) next

    smmry_cvvs <- smmrys_cvvs[[tstsetup]]
    smmry_cvvs_trad <- smmrys_cvvs[[tstsetup_trad]]

    expect_equal(
      smmry_cvvs[setdiff(names(smmry_cvvs), c("family", "selection"))],
      smmry_cvvs_trad[setdiff(names(smmry_cvvs_trad),
                              c("family", "selection"))],
      info = tstsetup
    )
    expect_identical(
      smmry_cvvs$selection[, c("size", "solution_terms")],
      smmry_cvvs_trad$selection[, c("size", "solution_terms")],
      info = tstsetup
    )
    # Exclude statistics which are not supported for the augmented-data
    # projection:
    smmry_cvvs_trad$selection <- smmry_cvvs_trad$selection[
      , -grep("mse\\.|auc\\.", names(smmry_cvvs_trad$selection)), drop = FALSE
    ]
    is_kfold <- identical(
      args_cvvs[[args_smmry_cvvs[[tstsetup]]$tstsetup_vsel]]$cv_method,
      "kfold"
    )
    if (is_kfold) {
      tol_smmry <- 1e-5
    } else {
      tol_smmry <- 1e-6
    }
    compare_exp <- function(e) {
      stopifnot(is_kfold)

      smmry_pd <- smmry_cvvs$selection
      se_cols <- grep("lpd\\.se$", names(smmry_pd), value = TRUE)
      elpd_cols <- setdiff(grep("elpd", names(smmry_pd), value = TRUE),
                           se_cols)
      smmry_pd[elpd_cols] <- exp(smmry_pd[elpd_cols])
      mlpd_cols <- setdiff(grep("mlpd", names(smmry_pd), value = TRUE),
                           se_cols)

      smmry_pd_trad <- smmry_cvvs_trad$selection
      se_cols_trad <- grep("lpd\\.se$", names(smmry_pd_trad), value = TRUE)
      elpd_cols_trad <- setdiff(grep("elpd", names(smmry_pd_trad), value = TRUE),
                                se_cols_trad)
      smmry_pd_trad[elpd_cols_trad] <- exp(smmry_pd_trad[elpd_cols_trad])
      mlpd_cols_trad <- setdiff(grep("mlpd", names(smmry_pd_trad), value = TRUE),
                                se_cols_trad)

      expect_equal(smmry_pd[, setdiff(names(smmry_pd),
                                      c(mlpd_cols, se_cols))],
                   smmry_pd_trad[, setdiff(names(smmry_pd_trad),
                                           c(mlpd_cols_trad, se_cols_trad))],
                   tolerance = 1e-10, info = tstsetup)
      expect_equal(smmry_pd[, c(mlpd_cols, se_cols)],
                   smmry_pd_trad[, c(mlpd_cols_trad, se_cols_trad)],
                   tolerance = 1e-3, info = tstsetup)
    }
    tryCatch(expect_equal(smmry_cvvs$selection, smmry_cvvs_trad$selection,
                          tolerance = tol_smmry, info = tstsetup),
             error = compare_exp)
  }
})

# Teardown ----------------------------------------------------------------

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
