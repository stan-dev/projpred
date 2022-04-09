context("miscellaneous")

# .get_refdist() ----------------------------------------------------------

## ndraws and nclusters ---------------------------------------------------

test_that(paste(
  ".get_refdist(): `ndraws = NULL` and `nclusters = NULL` leads to",
  "`ndraws = S` (and `nclusters = NULL`)"
), {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  set.seed(seed2_tst)
  for (tstsetup in names(refmods)) {
    refdist <- .get_refdist(refmods[[tstsetup]])
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
        refdist <- .get_refdist(refmods[[tstsetup]],
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
          info_str = paste(tstsetup, ndraws_crr, nclusters_crr, sep = "__")
        )
      }
    }
  }
  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
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
    nms_spclformul <- setdiff(
      grep("y_|xco", names(mf_spclformul), value = TRUE),
      "xco.1"
    )

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
      nms_stdformul <- setdiff(
        grep("y_|xco", names(mf_stdformul), value = TRUE),
        "xco.1"
      )
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

# pseudo_data() -----------------------------------------------------------

test_that(paste(
  "`pseudo_data(f = 0, [...], family = extend_family(gaussian()), [...])` is",
  "essentially an identity function (apart from offsets)."
), {
  mu_crr <- matrix(1:12, nrow = 3)
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
