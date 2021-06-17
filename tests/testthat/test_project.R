context("project")

test_that(paste(
  "`object` of class \"refmodel\" leads to correct output structure"
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

test_that(paste(
  "specifying `solution_terms` incorrectly leads to a warning or an error"
), {
  i <- "gauss"
  expect_error(project(refmods[[mod_nm]][[fam_nm]], nclusters = nclusters_pred_tst,
                       solution_terms = NULL),
               "is not an object of class \"vsel\"")
  for (solterms_crr in list(2, 1:3, "1", list(c("x.3", "x.5"),
                                              c("x.2", "x.4")))) {
    expect_warning(
      p <- project(refmods[[mod_nm]][[fam_nm]], nclusters = nclusters_pred_tst,
                   solution_terms = solterms_crr),
      paste("At least one element of `solution_terms` could not be found",
            "among the terms in the reference model"),
      info = as.character(solterms_crr)
    )
    expect_s3_class(p, "projection")
    expect_named(p, projection_nms, info = solterms_crr)
    expect_length(p$sub_fit, nclusters_pred_tst)
    expect_length(p$weights, nclusters_pred_tst)
    expect_length(p$dis, nclusters_pred_tst)
    SW(nprjdraws <- NROW(as.matrix(p)))
    expect_identical(nprjdraws, nclusters_pred_tst, info = solterms_crr)
    expect_identical(p$solution_terms, "1")
  }
})

test_that(paste(
  "specifying `solution_terms` correctly leads to correct output structure"
), {
  for (i in fam_nms) {
    for (solterms_crr in list(character(), "x.3", c("x.2", "x.4"))) {
      p <- project(refmods[[mod_nm]][[fam_nm]], nclusters = nclusters_pred_tst,
                   solution_terms = solterms_crr)
      expect_s3_class(p, "projection")
      expect_named(p, projection_nms, info = tstsetup)
      expect_length(p$sub_fit, nclusters_pred_tst)
      expect_length(p$weights, nclusters_pred_tst)
      expect_length(p$dis, nclusters_pred_tst)
      SW(nprjdraws <- NROW(as.matrix(p)))
      expect_identical(nprjdraws, nclusters_pred_tst, info = tstsetup)
      solterms_out <- if (length(solterms_crr) == 0) {
        "1"
      } else {
        solterms_crr
      }
      expect_identical(p$solution_terms, solterms_out)
    }
  }
})

test_that("specifying `ndraws` incorrectly leads to an error", {
  i <- "gauss"
  expect_error(project(refmods[[mod_nm]][[fam_nm]],
                       ndraws = NULL,
                       solution_terms = solterms_tst),
               "^!is\\.null\\(ndraws\\) is not TRUE$", info = tstsetup)
})

test_that(paste(
  "specifying `ndraws` and/or `nclusters` too big causes them to be cut off",
  "at the number of posterior draws in the reference model"
), {
  i <- "gauss"
  S <- nrow(as.matrix(fit_list[[i]]))
  for (ndraws_crr in list(S + 1L)) {
    for (nclusters_crr in list(NULL, S + 1L)) {
      tstsetup <- unlist(nlist(i, ndraws_crr, nclusters_crr))
      p <- project(refmods[[mod_nm]][[fam_nm]],
                   ndraws = ndraws_crr,
                   nclusters = nclusters_crr,
                   solution_terms = solterms_tst)
      expect_s3_class(p, "projection")
      expect_named(p, projection_nms, info = tstsetup)
      expect_length(p$sub_fit, S)
      expect_length(p$weights, S)
      expect_length(p$dis, S)
      SW(nprjdraws <- NROW(as.matrix(p)))
      expect_identical(nprjdraws, S, info = tstsetup)
      solterms_out <- if (length(solterms_tst) == 0) "1" else solterms_tst
      expect_identical(p$solution_terms, solterms_out)
    }
  }
})

test_that(paste(
  "specifying `ndraws` and/or `nclusters` correctly leads to correct output",
  "structure"
), {
  for (i in fam_nms) {
    for (ndraws_crr in list(1L, 20L, 21L)) {
      for (nclusters_crr in list(NULL, 1L, 2L, 3L)) {
        tstsetup <- unlist(nlist(i, ndraws_crr, nclusters_crr))
        p <- project(refmods[[mod_nm]][[fam_nm]],
                     ndraws = ndraws_crr,
                     nclusters = nclusters_crr,
                     solution_terms = solterms_tst)
        expect_s3_class(p, "projection")
        expect_named(p, projection_nms, info = tstsetup)
        nprjdraws_out <- if (!is.null(nclusters_crr)) {
          nclusters_crr
        } else {
          ndraws_crr
        }
        nprjdraws_sub_fit <- if (nprjdraws_out == 1) {
          length(sub_fit_nms)
        } else {
          nprjdraws_out
        }
        expect_length(p$sub_fit, nprjdraws_sub_fit)
        expect_length(p$weights, nprjdraws_out)
        expect_length(p$dis, nprjdraws_out)
        SW(nprjdraws <- NROW(as.matrix(p)))
        expect_identical(nprjdraws, nprjdraws_out, info = tstsetup)
        solterms_out <- if (length(solterms_tst) == 0) "1" else solterms_tst
        expect_identical(p$solution_terms, solterms_out)
        if (nprjdraws_out == 1) {
          expect_identical(p$weights, 1, info = tstsetup)
        }
      }
    }
  }
})

test_that("specifying `seed` correctly leads to reproducible results", {
  i <- "gauss"
  p1 <- project(refmods[[mod_nm]][[fam_nm]],
                nclusters = nclusters_pred_tst,
                solution_terms = solterms_tst,
                seed = seed_tst)
  p2 <- project(refmods[[mod_nm]][[fam_nm]],
                nclusters = nclusters_pred_tst,
                solution_terms = solterms_tst,
                seed = seed_tst + 1L)
  p3 <- project(refmods[[mod_nm]][[fam_nm]],
                nclusters = nclusters_pred_tst,
                solution_terms = solterms_tst,
                seed = seed_tst)
  p4 <- project(refmods[[mod_nm]][[fam_nm]],
                nclusters = nclusters_pred_tst,
                solution_terms = solterms_tst)

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
