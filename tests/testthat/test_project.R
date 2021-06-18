context("project")

test_that(paste(
  "`object` of class \"refmodel\", correctly specified `solution_terms`, and",
  "correctly specified `ndraws` or `nclusters` lead to correct output structure"
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
    projection_tester(
      prjs[[tstsetup]],
      solterms_expected = args_prj[[tstsetup]]$solution_terms,
      nprjdraws_expected = nprjdraws,
      info_str = tstsetup
    )
  }
})

test_that(paste(
  "specifying `solution_terms` incorrectly leads to a warning or an error"
), {
  for (mod_nm in mod_nms["glm"]) {
    for (fam_nm in fam_nms["gauss"]) {
      expect_error(project(refmods[[!!mod_nm]][[!!fam_nm]],
                           solution_terms = NULL),
                   "is not an object of class \"vsel\"")
      for (solterms_crr in list(2, 1:3, "1", list(solterms_x, solterms_x))) {
        tstsetup_crr <- paste(
          c(mod_nm, fam_nm, paste(solterms_crr, collapse = ", ")),
          collapse = "."
        )
        expect_warning(
          p <- project(refmods[[mod_nm]][[fam_nm]],
                       nclusters = nclusters_pred_tst,
                       solution_terms = solterms_crr,
                       seed = seed_tst),
          paste("At least one element of `solution_terms` could not be found",
                "among the terms in the reference model"),
          info = tstsetup_crr
        )
        projection_tester(p,
                          solterms_expected = character(),
                          nprjdraws_expected = nclusters_pred_tst,
                          info_str = tstsetup_crr)
      }
    }
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by varsel()) and correctly specified",
  "`nterms` lead to correct output structure"
), {
  skip_if_not(exists("prjs_vs"))
  for (tstsetup in names(prjs_vs)) {
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    tstsetup_vs <- grep(
      paste0("^", args_prj_vs[[tstsetup]]$mod_nm,
             "\\.", args_prj_vs[[tstsetup]]$fam_nm),
      names(vss),
      value = TRUE
    )
    stopifnot(length(tstsetup_vs) == 1)
    if (is.null(nterms_crr)) {
      # Subtract 1L for the intercept:
      nterms_crr <- vss[[tstsetup_vs]]$suggested_size - 1L
    }
    if (length(nterms_crr) == 1) {
      if (nterms_crr > 0) {
        solterms_expected_crr <- vss[[tstsetup_vs]]$solution_terms[
          seq_len(nterms_crr)
        ]
      } else {
        solterms_expected_crr <- "1"
      }
      projection_tester(
        prjs_vs[[tstsetup]],
        solterms_expected = solterms_expected_crr,
        nprjdraws_expected = args_prj_vs[[tstsetup]]$nclusters,
        info_str = tstsetup
      )
    } else {
      proj_list_tester(
        prjs_vs[[tstsetup]],
        len_expected = length(nterms_crr),
        is_seq = all(diff(nterms_crr) == 1),
        info_str = tstsetup,
        nprjdraws_expected = args_prj_vs[[tstsetup]]$nclusters,
        fam_expected = prjs_vs[[tstsetup]][[1]]$family,
        prjdraw_weights_expected = prjs_vs[[tstsetup]][[1]]$weights
      )
    }
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by cv_varsel()) and correctly specified",
  "`nterms` lead to correct output structure"
), {
  skip_if_not(exists("prjs_cvvs"))
  for (tstsetup in names(prjs_cvvs)) {
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    tstsetup_cvvs <- grep(
      paste0("^", args_prj_cvvs[[tstsetup]]$mod_nm,
             "\\.", args_prj_cvvs[[tstsetup]]$fam_nm),
      names(cvvss),
      value = TRUE
    )
    stopifnot(length(tstsetup_cvvs) == 1)
    if (is.null(nterms_crr)) {
      # Subtract 1L for the intercept:
      nterms_crr <- cvvss[[tstsetup_cvvs]]$suggested_size - 1L
    }
    if (length(nterms_crr) == 1) {
      if (nterms_crr > 0) {
        solterms_expected_crr <- cvvss[[tstsetup_cvvs]]$solution_terms[
          seq_len(nterms_crr)
        ]
      } else {
        solterms_expected_crr <- "1"
      }
      projection_tester(
        prjs_cvvs[[tstsetup]],
        solterms_expected = solterms_expected_crr,
        nprjdraws_expected = args_prj_cvvs[[tstsetup]]$nclusters,
        info_str = tstsetup
      )
    } else {
      proj_list_tester(
        prjs_cvvs[[tstsetup]],
        len_expected = length(nterms_crr),
        is_seq = all(diff(nterms_crr) == 1),
        info_str = tstsetup,
        nprjdraws_expected = args_prj_cvvs[[tstsetup]]$nclusters,
        fam_expected = prjs_cvvs[[tstsetup]][[1]]$family,
        prjdraw_weights_expected = prjs_cvvs[[tstsetup]][[1]]$weights
      )
    }
  }
})

test_that(paste(
  "error if `object` is not of class \"vsel\" and `solution_terms` is provided",
  "neither"
), {
  expect_error(project(fits$glm$gauss), "is not an object of class \"vsel\"")
})

test_that("specifying `nterms` incorrectly leads to an error", {
  skip_if_not(exists("vss"))
  for (tstsetup in grep("^glm\\.gauss", names(vss), value = TRUE)[1]) {
    for (nterms_crr in nterms_unavail) {
      expect_error(project(vss[[!!tstsetup]], nterms = !!nterms_crr),
                   paste("Cannot perform the projection with", max(nterms_crr),
                         "variables"))
    }
    for (nterms_crr in list(neg = -1, char = "a", dafr = dat)) {
      expect_error(project(vss[[!!tstsetup]], nterms = !!nterms_crr),
                   "must contain non-negative values")
    }
  }
})

test_that("specifying `ndraws` incorrectly leads to an error", {
  for (mod_nm in mod_nms["glm"]) {
    for (fam_nm in fam_nms["gauss"]) {
      expect_error(project(refmods[[!!mod_nm]][[!!fam_nm]],
                           ndraws = NULL,
                           solution_terms = solterms_x),
                   "^!is\\.null\\(ndraws\\) is not TRUE$")
    }
  }
})

test_that(paste(
  "specifying `ndraws` and/or `nclusters` too big causes them to be cut off",
  "at the number of posterior draws in the reference model"
), {
  for (mod_nm in mod_nms["glm"]) {
    for (fam_nm in fam_nms["gauss"]) {
      S <- nrow(as.matrix(fits[[mod_nm]][[fam_nm]]))
      for (ndraws_crr in list(S + 1L)) {
        for (nclusters_crr in list(NULL, S + 1L)) {
          tstsetup <- paste(c(mod_nm, fam_nm, ndraws_crr, nclusters_crr),
                            collapse = ".")
          p <- project(refmods[[mod_nm]][[fam_nm]],
                       ndraws = ndraws_crr,
                       nclusters = nclusters_crr,
                       solution_terms = solterms_x,
                       seed = seed_tst)
          projection_tester(p,
                            solterms_expected = solterms_x,
                            nprjdraws_expected = S,
                            info_str = tstsetup)
        }
      }
    }
  }
})

test_that("specifying `seed` correctly leads to reproducible results", {
  for (mod_nm in mod_nms["glm"]) {
    for (fam_nm in fam_nms["gauss"]) {
      tstsetup <- paste(c(mod_nm, fam_nm), collapse = ".")
      p1 <- project(refmods[[mod_nm]][[fam_nm]],
                    nclusters = nclusters_pred_tst,
                    solution_terms = solterms_x,
                    seed = seed_tst)
      p2 <- project(refmods[[mod_nm]][[fam_nm]],
                    nclusters = nclusters_pred_tst,
                    solution_terms = solterms_x,
                    seed = seed_tst + 1L)
      p3 <- project(refmods[[mod_nm]][[fam_nm]],
                    nclusters = nclusters_pred_tst,
                    solution_terms = solterms_x,
                    seed = seed_tst)
      p4 <- project(refmods[[mod_nm]][[fam_nm]],
                    nclusters = nclusters_pred_tst,
                    solution_terms = solterms_x)

      # Expected equality:
      expect_true(isTRUE(all.equal(p1, p3)), info = tstsetup)
      # The resulting objects are even identical when ignoring the environments of
      # functions:
      expect_identical(p1, p3, info = tstsetup, ignore.environment = TRUE)

      # Expected inequality:
      expect_false(isTRUE(all.equal(p1, p2)), info = tstsetup)
      expect_false(isTRUE(all.equal(p1, p4)), info = tstsetup)
      expect_false(isTRUE(all.equal(p2, p3)), info = tstsetup)
      expect_false(isTRUE(all.equal(p2, p4)), info = tstsetup)
      expect_false(isTRUE(all.equal(p3, p4)), info = tstsetup)
    }
  }
})

test_that(paste(
  "projecting the reference model onto the full model (i.e.,",
  "itself) does not change results on average (even though this is not",
  "guaranteed; see the comments)"
), {
  # NOTE: Projecting the reference model onto the full model (i.e., itself)
  # does not necessarily have to give results close to the reference model's
  # since in contrast to the reference model, the projection is "fitting to
  # the fit" of the reference model, not to the observed response. The fact
  # that the tolerance for the Gaussian reference model needs to be
  # increased here (see below) might be an indicator for this inequality.
  tol <- setNames(rep(1e-3, length(fam_nms)), fam_nms)
  tol["gauss"] <- 0.25

  for (i in fam_nms) {
    draws <- as.data.frame(fit_list[[i]])
    alpha_ref <- draws[, "(Intercept)"]
    beta_ref <- draws[, 1 + seq_len(nterms), drop = FALSE]
    S <- nrow(draws)
    proj <- project(refmods[[mod_nm]][[fam_nm]],
                    solution_terms = paste0("x.", 1:nterms),
                    ndraws = S)

    # test alpha and beta
    coefs <- as.matrix(proj)
    dalpha <- abs(mean(coefs[, 1]) - mean(alpha_ref))
    order <- match(colnames(fit_list[[i]]$data), proj$solution_terms)
    order <- order[!is.na(order)]
    dbeta <- max(abs(colMeans(coefs[, -1, drop = FALSE][, order]) -
                       colMeans(beta_ref)))
    expect_lt(dalpha, tol[!!i])
    expect_lt(dbeta, tol[!!i])
  }
})
