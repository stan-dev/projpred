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
                                      ndiscrete = ncat_tst, class = "augmat"))

  # With "margin_draws = 3":
  xtst_m3 <- 0.1
  stopifnot(length(xtst_m3) == 1)
  arr_m3 <- array(xtst_m3, dim = c(1, 1, 1))
  augmat_m3 <- arr2augmat(arr_m3)
  augmat_m3_ch <- structure(matrix(xtst_m3), ndiscrete = 1L, class = "augmat")
  expect_identical(augmat_m3, augmat_m3_ch)
  arr_m3_ch <- augmat2arr(augmat_m3)
  expect_identical(arr_m3_ch, arr_m3)

  # With "margin_draws = 1":
  xtst_m1 <- seq(-1, 1, by = 0.1)
  stopifnot(length(xtst_m1) == 21)
  arr_m1 <- array(xtst_m1, dim = c(7, 3, 1))
  augmat_m1 <- arr2augmat(arr_m1, margin_draws = 1)
  augmat_m1_ch <- structure(matrix(xtst_m1, ncol = 7, byrow = TRUE),
                            ndiscrete = 1L, class = "augmat")
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
  # Subsetting:
  # Subsetting with selection of all elements:
  expect_identical(augmtst[], augmtst)
  # Subsetting columns:
  expect_identical(augmtst[, 1:2],
                   structure(matrix(head(xtst, nobs_orig_tst * ncat_tst * 2),
                                    ncol = 2),
                             ndiscrete = ncat_tst, class = "augmat"))
  expect_identical(augmtst[, 1, drop = FALSE],
                   structure(matrix(head(xtst, nobs_orig_tst * ncat_tst)),
                             ndiscrete = ncat_tst, class = "augmat"))
  expect_identical(augmtst[, 1],
                   structure(head(xtst, nobs_orig_tst * ncat_tst),
                             ndiscrete = ncat_tst, class = "augvec"))
  # Subsetting rows:
  xrows1 <- xtst[as.vector(t(sapply(
    nobs_orig_tst * (seq_len(ncat_tst) - 1L) + 1L,
    function(idx) {
      nobs_orig_tst * ncat_tst * (seq_len(ndraws_tst) - 1L) + idx
    }
  )))]
  xrows1 <- matrix(xrows1, nrow = ncat_tst)
  expect_identical(augmtst[nobs_orig_tst * (seq_len(ncat_tst) - 1L) + 1L, ,
                           drop = FALSE],
                   structure(xrows1, ndiscrete = ncat_tst, class = "augmat"))
  # "Illegally" subsetting not only observations, but also the (possibly latent)
  # response categories (within projpred, this should never be used, but testing
  # here anyway for the sake of completeness):
  xrow1 <- xtst[nobs_orig_tst * ncat_tst * (seq_len(ndraws_tst) - 1L) + 1L]
  expect_identical(augmtst[1, ],
                   structure(xrow1, ndiscrete = ncat_tst, class = "augvec"))
  expect_identical(augmtst[1, , drop = FALSE],
                   structure(t(xrow1), ndiscrete = ncat_tst, class = "augmat"))
  # Subsetting rows and columns:
  expect_identical(
    augmtst[sort(c(nobs_orig_tst * (seq_len(ncat_tst) - 1L) + 1L,
                   nobs_orig_tst * (seq_len(ncat_tst) - 1L) + 2L)), 1:2],
    structure(cbind(as.vector(arrtst[1:2, , 1]), as.vector(arrtst[1:2, , 2])),
              ndiscrete = ncat_tst, class = "augmat")
  )
  expect_identical(augmtst[nobs_orig_tst * (seq_len(ncat_tst) - 1L) + 1L, 1],
                   structure(xrows1[, 1], ndiscrete = ncat_tst,
                             class = "augvec"))
  # "Illegally" subsetting not only observations, but also the (possibly latent)
  # response categories (within projpred, this should never be used, but testing
  # here anyway for the sake of completeness):
  expect_identical(augmtst[1:3, 1:2],
                   structure(arrtst[1:3, 1, 1:2], ndiscrete = ncat_tst,
                             class = "augmat"))
  expect_identical(augmtst[1, 1],
                   structure(head(xtst, 1), ndiscrete = ncat_tst,
                             class = "augvec"))
  expect_identical(augmtst[1, 1, drop = FALSE],
                   structure(matrix(head(xtst, 1)), ndiscrete = ncat_tst,
                             class = "augmat"))

  # Assigning:
  augmtst[1, 1] <- 0
  expect_identical(augmtst, structure(matrix(c(0, xtst[-1]), ncol = ndraws_tst),
                                      ndiscrete = ncat_tst, class = "augmat"))
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
  # Subsetting:
  expect_identical(augvtst[], augvtst)
  xrows1_1 <- xtst[nobs_orig_tst * (seq_len(ncat_tst) - 1L) + 1L]
  expect_identical(augvtst[nobs_orig_tst * (seq_len(ncat_tst) - 1L) + 1L],
                   structure(xrows1_1, ndiscrete = ncat_tst, class = "augvec"))
  expect_identical(augvtst[1:2],
                   structure(head(xtst, 2), ndiscrete = ncat_tst,
                             class = "augvec"))
  expect_identical(augvtst[1],
                   structure(head(xtst, 1), ndiscrete = ncat_tst,
                             class = "augvec"))
  expect_identical(augvtst[integer()],
                   structure(integer(), ndiscrete = ncat_tst, class = "augvec"))

  # Assigning:
  augvtst[1] <- 0
  expect_identical(augvtst,
                   structure(c(0, head(xtst, nobs_orig_tst * ncat_tst)[-1]),
                             ndiscrete = ncat_tst, class = "augvec"))
})

## catmaxprb() ------------------------------------------------------------

test_that("catmaxprb()", {
  stopifnot(ncat_tst == 2)
  expect_identical(catmaxprb(augvtst, lvls = c("0", "1")),
                   factor(rep("1", nobs_orig_tst), levels = c("0", "1")))

  # With probabilities (and a single observation):
  augvtst_pr <- structure(c(0.7, 0.3), ndiscrete = 2L, class = "augvec")
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
    pref_aug <- get_refdist(refmod_crr,
                            ndraws = args_prj_i$ndraws,
                            nclusters = args_prj_i$nclusters)
    set.seed(args_prj_i_trad$seed)
    pref_trad <- get_refdist(refmod_crr_trad,
                             ndraws = args_prj_i_trad$ndraws,
                             nclusters = args_prj_i_trad$nclusters)

    eta_aug <- refmod_crr$family$linkfun(pref_aug$mu)
    eta_trad <- refmod_crr_trad$family$linkfun(pref_trad$mu)
    expect_identical(structure(unclass(eta_aug), ndiscrete = NULL), eta_trad,
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

    cl_used_gt1 <- ndr_ncl_dtls(args_prj[[tstsetup]])$clust_used_gt1
    prjmat <- as.matrix(prjs[[tstsetup]],
                        allow_nonconst_wdraws_prj = cl_used_gt1)
    prjmat_trad <- as.matrix(prjs[[tstsetup_trad]],
                             allow_nonconst_wdraws_prj = cl_used_gt1)

    tol_coefs <- ifelse(
      args_prj[[tstsetup]]$mod_nm == "glmm" &&
        any(grepl("^\\(.*\\)$", args_prj[[tstsetup]]$predictor_terms)),
      1e-3, 1e-5
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
    cl_used_gt1 <- ndr_ncl_dtls(args_prj_i)$clust_used_gt1
    prjmat <- as.matrix(prj, allow_nonconst_wdraws_prj = cl_used_gt1)
    prjmat_trad <- as.matrix(prj_trad, allow_nonconst_wdraws_prj = cl_used_gt1)

    tol_coefs <- ifelse(cl_used_gt1, 1e-9, 1e-14)
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
    pl$pred <- structure(matrix(pl$pred,
                                nrow = dim(pl$pred)[1],
                                ncol = dim(pl$pred)[2]),
                         wdraws_prj = attr(pl$pred, "wdraws_prj"))
    dimnames(pl$pred) <- list(NULL, as.character(seq_len(ncol(pl$pred))))
    dimnames(pl$lpd) <- list(NULL, as.character(seq_len(ncol(pl$lpd))))
    pl_trad <- pls[[tstsetup_trad]]

    tol_lpreds <- ifelse(
      args_prj[[tstsetup]]$mod_nm == "glmm" &&
        any(grepl("^\\(.*\\)$", args_prj[[tstsetup]]$predictor_terms)),
      1e-4, 1e-5
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
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
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
    summs_ref$mu <- structure(unclass(summs_ref$mu), ndiscrete = NULL)
    summs_ref$mu <- summs_ref$mu[(nobsv + 1):(2 * nobsv)]
    tol_ref <- 1e1 * .Machine$double.eps
    if (args_vs[[tstsetup]]$mod_nm == "glmm" ||
        (args_vs[[tstsetup]]$mod_nm == "glm" &&
         args_vs[[tstsetup]]$fam_nm %in% c("brnll", "binom"))) {
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
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(smmrys_vs)) next

    smmry_vs <- smmrys_vs[[tstsetup]]
    smmry_vs_trad <- smmrys_vs[[tstsetup_trad]]

    expect_equal(
      smmry_vs[setdiff(names(smmry_vs), c("family", "perf_sub"))],
      smmry_vs_trad[setdiff(names(smmry_vs_trad), c("family", "perf_sub"))],
      info = tstsetup
    )
    expect_identical(
      smmry_vs$perf_sub[, c("size", "ranking_fulldata")],
      smmry_vs_trad$perf_sub[, c("size", "ranking_fulldata")],
      info = tstsetup
    )
    expect_equal(smmry_vs$perf_sub, smmry_vs_trad$perf_sub, tolerance = 1e-6,
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
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    tstsetup_trad <- sub("\\.augdat_stats\\.", ".binom_stats.", tstsetup_trad)
    if (!tstsetup_trad %in% names(smmrys_vs)) next

    smmry_vs <- smmrys_vs[[tstsetup]]
    smmry_vs_trad <- smmrys_vs[[tstsetup_trad]]

    expect_equal(
      smmry_vs[setdiff(names(smmry_vs), c("family", "perf_sub", "perf_ref"))],
      smmry_vs_trad[setdiff(names(smmry_vs_trad),
                            c("family", "perf_sub", "perf_ref"))],
      info = tstsetup
    )
    expect_identical(
      smmry_vs$perf_sub[, c("size", "ranking_fulldata", "cv_proportions_diag")],
      smmry_vs_trad$perf_sub[
        , c("size", "ranking_fulldata", "cv_proportions_diag")
      ],
      info = tstsetup
    )
    # Exclude statistics which are not supported for the augmented-data
    # projection:
    smmry_vs_trad$perf_sub <- smmry_vs_trad$perf_sub[
      , -grep("mse|R2|auc", names(smmry_vs_trad$perf_sub)), drop = FALSE
    ]
    smmry_vs_trad$perf_ref <- smmry_vs_trad$perf_ref[
      -grep("mse|R2|auc", names(smmry_vs_trad$perf_ref))
    ]
    expect_equal(smmry_vs$perf_sub, smmry_vs_trad$perf_sub,
                 tolerance = 1e-6, info = tstsetup)
    expect_equal(smmry_vs$perf_ref, smmry_vs_trad$perf_ref, info = tstsetup)
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
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(cvvss)) next

    cvvs <- cvvss[[tstsetup]]
    cvvs_trad <- cvvss[[tstsetup_trad]]

    is_kfold <- identical(args_cvvs[[tstsetup]]$cv_method, "kfold")

    summs_sub <- cvvs$summaries$sub
    summs_sub_mu <- do.call(cbind, lapply(summs_sub, "[[", "mu"))[
      (nobsv + 1):(2 * nobsv), , drop = FALSE
    ]
    summs_sub_trad <- cvvs_trad$summaries$sub
    summs_sub_mu_trad <- do.call(cbind, lapply(summs_sub_trad, "[[", "mu"))
    # Sometimes, there are seemingly large differences on logit scale which are
    # probably due to numerical overflow of probabilities (i.e., on probability
    # scale) towards 1 or underflow towards zero. Thus, compare on probability
    # scale if the comparison on logit scale fails:
    tryCatch(expect_equal(summs_sub_mu, summs_sub_mu_trad,
                          tolerance = 1e-6, info = tstsetup),
             error = function(e) {
               stopifnot(is_kfold)
               expect_equal(ilinkfun_raw(summs_sub_mu, link_nm = "logit"),
                            ilinkfun_raw(summs_sub_mu_trad, link_nm = "logit"),
                            tolerance = 1e-5, info = tstsetup)
             })
    summs_sub_lppd <- do.call(cbind, lapply(summs_sub, "[[", "lppd"))
    summs_sub_lppd_trad <- unname(
      do.call(cbind, lapply(summs_sub_trad, "[[", "lppd"))
    )
    # Sometimes, there are seemingly large differences on log scale which are
    # probably due to numerical overflow of probabilities (i.e., on exp scale)
    # towards 1 or underflow towards zero. Thus, compare on exp scale if the
    # comparison on log scale fails:
    tryCatch(expect_equal(summs_sub_lppd, summs_sub_lppd_trad,
                          tolerance = 1e-6, info = tstsetup),
             error = function(e) {
               stopifnot(is_kfold)
               expect_equal(exp(summs_sub_lppd), exp(summs_sub_lppd_trad),
                            tolerance = 1e-5, info = tstsetup)
             })
    summs_ref <- cvvs$summaries$ref
    summs_ref$mu <- structure(unclass(summs_ref$mu), ndiscrete = NULL)
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
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(smmrys_cvvs)) next

    smmry_cvvs <- smmrys_cvvs[[tstsetup]]
    smmry_cvvs_trad <- smmrys_cvvs[[tstsetup_trad]]

    expect_equal(
      smmry_cvvs[setdiff(names(smmry_cvvs), c("family", "perf_sub"))],
      smmry_cvvs_trad[setdiff(names(smmry_cvvs_trad),
                              c("family", "perf_sub"))],
      info = tstsetup
    )
    expect_identical(
      smmry_cvvs$perf_sub[
        , c("size", "ranking_fulldata", "cv_proportions_diag")
      ],
      smmry_cvvs_trad$perf_sub[
        , c("size", "ranking_fulldata", "cv_proportions_diag")
      ],
      info = tstsetup
    )
    is_kfold <- identical(
      args_cvvs[[args_smmry_cvvs[[tstsetup]]$tstsetup_vsel]]$cv_method,
      "kfold"
    )
    if (is_kfold) {
      tol_smmry <- 1e-4
    } else {
      tol_smmry <- 1e-6
    }
    compare_exp <- function(e) {
      stopifnot(is_kfold)
      # Check that we have the default `stats` in this case, meaning only the
      # ELPD:
      stopifnot(is.null(args_smmry_cvvs[[tstsetup]]$stats))

      smmry_pd <- smmry_cvvs$perf_sub
      smmry_pd$elpd <- exp(smmry_pd$elpd)
      smmry_pd$elpd.lower <- exp(smmry_pd$elpd.lower)
      smmry_pd$elpd.upper <- exp(smmry_pd$elpd.upper)

      smmry_pd_trad <- smmry_cvvs_trad$perf_sub
      smmry_pd_trad$elpd <- exp(smmry_pd_trad$elpd)
      smmry_pd_trad$elpd.lower <- exp(smmry_pd_trad$elpd.lower)
      smmry_pd_trad$elpd.upper <- exp(smmry_pd_trad$elpd.upper)

      expect_equal(smmry_pd[, setdiff(names(smmry_pd), "elpd.se")],
                   smmry_pd_trad[, setdiff(names(smmry_pd_trad), "elpd.se")],
                   tolerance = tol_smmry, info = tstsetup)
      expect_equal(smmry_pd$elpd.se, smmry_pd_trad$elpd.se, tolerance = 1e-3,
                   info = tstsetup)
    }
    tryCatch(expect_equal(smmry_cvvs$perf_sub, smmry_cvvs_trad$perf_sub,
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
    tstsetup_trad <- sub("\\.augdat\\.", ".trad_compare.", tstsetup)
    tstsetup_trad <- sub("\\.augdat_stats\\.", ".binom_stats.", tstsetup_trad)
    if (!tstsetup_trad %in% names(smmrys_cvvs)) next

    smmry_cvvs <- smmrys_cvvs[[tstsetup]]
    smmry_cvvs_trad <- smmrys_cvvs[[tstsetup_trad]]

    expect_equal(
      smmry_cvvs[setdiff(names(smmry_cvvs),
                         c("family", "perf_sub", "perf_ref"))],
      smmry_cvvs_trad[setdiff(names(smmry_cvvs_trad),
                              c("family", "perf_sub", "perf_ref"))],
      info = tstsetup
    )
    expect_identical(
      smmry_cvvs$perf_sub[
        , c("size", "ranking_fulldata", "cv_proportions_diag")
      ],
      smmry_cvvs_trad$perf_sub[
        , c("size", "ranking_fulldata", "cv_proportions_diag")
      ],
      info = tstsetup
    )
    # Exclude statistics which are not supported for the augmented-data
    # projection:
    smmry_cvvs_trad$perf_sub <- smmry_cvvs_trad$perf_sub[
      , -grep("mse|R2|auc", names(smmry_cvvs_trad$perf_sub)), drop = FALSE
    ]
    smmry_cvvs_trad$perf_ref <- smmry_cvvs_trad$perf_ref[
      -grep("mse|R2|auc", names(smmry_cvvs_trad$perf_ref))
    ]
    is_kfold <- identical(
      args_cvvs[[args_smmry_cvvs[[tstsetup]]$tstsetup_vsel]]$cv_method,
      "kfold"
    )
    if (is_kfold) {
      tol_smmry <- 1e-4
    } else {
      tol_smmry <- 1e-6
    }
    compare_exp <- function(e) {
      stopifnot(is_kfold)

      smmry_pd <- smmry_cvvs$perf_sub
      se_cols <- grep("lpd\\.se$", names(smmry_pd), value = TRUE)
      elpd_cols <- setdiff(grep("elpd", names(smmry_pd), value = TRUE),
                           se_cols)
      smmry_pd[elpd_cols] <- exp(smmry_pd[elpd_cols])
      mlpd_cols <- setdiff(grep("mlpd", names(smmry_pd), value = TRUE),
                           se_cols)
      smmry_pd[mlpd_cols] <- exp(smmry_pd[mlpd_cols])

      smmry_pd_trad <- smmry_cvvs_trad$perf_sub
      se_cols_trad <- grep("lpd\\.se$", names(smmry_pd_trad), value = TRUE)
      elpd_cols_trad <- setdiff(
        grep("elpd", names(smmry_pd_trad), value = TRUE),
        se_cols_trad
      )
      smmry_pd_trad[elpd_cols_trad] <- exp(smmry_pd_trad[elpd_cols_trad])
      mlpd_cols_trad <- setdiff(
        grep("mlpd", names(smmry_pd_trad), value = TRUE),
        se_cols_trad
      )
      smmry_pd_trad[mlpd_cols_trad] <- exp(smmry_pd_trad[mlpd_cols_trad])

      expect_equal(smmry_pd[setdiff(names(smmry_pd), se_cols)],
                   smmry_pd_trad[setdiff(names(smmry_pd_trad), se_cols_trad)],
                   tolerance = tol_smmry, info = tstsetup)
      expect_equal(smmry_pd[se_cols],
                   smmry_pd_trad[se_cols_trad],
                   tolerance = 1e-3, info = tstsetup)
    }
    tryCatch(expect_equal(smmry_cvvs$perf_sub, smmry_cvvs_trad$perf_sub,
                          tolerance = tol_smmry, info = tstsetup),
             error = compare_exp)
    expect_equal(smmry_cvvs$perf_ref, smmry_cvvs_trad$perf_ref,
                 tolerance = tol_smmry, info = tstsetup)
  }
})

# Teardown ----------------------------------------------------------------

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
