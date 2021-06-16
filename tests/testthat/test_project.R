context("project")

test_that(paste(
  "error if `object` is not of class \"vsel\" and `solution_terms` is provided",
  "neither"
), {
  expect_error(project(fits$glm$gauss), "is not an object of class \"vsel\"")
})

test_that(paste(
  "`object` of class \"vsel\" (created by varsel()) leads to correct output",
  "structure"
), {
  skip_if_not(exists("prj_nterms_vs"))
  prj_vsel_str_tester(prj_nterms_vs, fam_expected = prj_nterms_vs[[1]]$family)
})

test_that(paste(
  "`object` of class \"vsel\" (created by cv_varsel()) leads to correct output",
  "structure"
), {
  skip_if_not(exists("prj_nterms_cvvs"))
  prj_vsel_str_tester(prj_nterms_cvvs,
                      fam_expected = prj_nterms_cvvs[[1]]$family)
})

test_that(paste(
  "`object` of class \"refmodel\" leads to correct output structure"
), {
  for (i in fam_nms) {
    p <- project(refmods[[mod_nm]][[fam_nm]], solution_terms = solterms_tst,
                 nclusters = nclusters_pred_tst)
    expect_s3_class(p, "projection")
    expect_named(p, projection_nms, info = tstsetup)
    expect_length(p$sub_fit, nclusters_pred_tst)
    expect_length(p$weights, nclusters_pred_tst)
    expect_length(p$dis, nclusters_pred_tst)
    SW(nprjdraws <- NROW(as.matrix(p)))
    expect_identical(nprjdraws, nclusters_pred_tst, info = tstsetup)
    solterms_out <- if (length(solterms_tst) == 0) "1" else solterms_tst
    expect_identical(p$solution_terms, solterms_out)
  }
})

test_that(paste(
  "a fitted model `object` leads to correct output structure and the default",
  "`ndraws` is as expected"
), {
  for (i in fam_nms) {
    if (i == "binom") {
      # For the binomial family with > 1 trials, we expect a warning (see
      # GitHub issue #136):
      warn_prj_expect <- paste("Using formula\\(x\\) is deprecated when x",
                               "is a character vector of length > 1")
    } else {
      warn_prj_expect <- NA
    }
    expect_warning(p <- project(fits$glm$gauss, solution_terms = solterms_tst),
                   warn_prj_expect, info = tstsetup)
    expect_s3_class(p, "projection")
    expect_named(p, projection_nms, info = tstsetup)
    expect_length(p$sub_fit, ndraws_default)
    expect_length(p$weights, ndraws_default)
    expect_length(p$dis, ndraws_default)
    expect_identical(NROW(as.matrix(p)), ndraws_default, info = tstsetup)
    solterms_out <- if (length(solterms_tst) == 0) "1" else solterms_tst
    expect_identical(p$solution_terms, solterms_out)
  }
})

test_that("specifying `nterms` incorrectly leads to an error", {
  i <- "gauss"
  expect_error(project(vss[[tstsetup]], nterms = 1000),
               "Cannot perform the projection with 1000 variables")
  expect_error(project(vss[[tstsetup]], nterms = -1),
               "must contain non-negative values")
  expect_error(project(vss[[tstsetup]], nterms = "a"),
               "must contain non-negative values")
  expect_error(project(vss[[tstsetup]], nterms = df_gauss),
               "must contain non-negative values")
})

test_that("specifying `nterms` correctly leads to correct output structure", {
  for (i in fam_nms) {
    for (nterms_tst in list(NULL, 0, 3, c(1, 3))) {
      p <- project(vss[[tstsetup]], nclusters = nclusters_pred_tst,
                   nterms = nterms_tst)
      out_size <- if (is.null(nterms_tst)) {
        suggest_size(vss[[tstsetup]])
      } else {
        nterms_tst
      }
      if (length(out_size) == 1) {
        expect_s3_class(p, "projection")
        expect_named(p, projection_nms, info = tstsetup)
        expect_length(p$sub_fit, nclusters_pred_tst)
        expect_length(p$weights, nclusters_pred_tst)
        expect_length(p$dis, nclusters_pred_tst)
        SW(nprjdraws <- NROW(as.matrix(p)))
        expect_identical(nprjdraws, nclusters_pred_tst, info = tstsetup)
        expect_length(p$solution_terms, max(out_size, 1))
        expect_equal(count_terms_chosen(p$solution_terms) - 1, out_size,
                     info = i)
      } else {
        expect_type(p, "list")
        expect_length(p, length(out_size))
        expect_true(.is_proj_list(p), info = tstsetup)

        prjdraw_weights <- p[[1]]$weights
        for (j in seq_along(p)) {
          expect_s3_class(p[[!!j]], "projection")
          expect_named(p[[!!j]], projection_nms, info = tstsetup)
          expect_length(p[[!!j]]$sub_fit, nclusters_pred_tst)
          expect_length(p[[!!j]]$weights, nclusters_pred_tst)
          expect_length(p[[!!j]]$dis, nclusters_pred_tst)
          SW(nprjdraws <- NROW(as.matrix(p[[!!j]])))
          expect_identical(nprjdraws, nclusters_pred_tst, info = tstsetup)
          expect_length(p[[!!j]]$solution_terms, max(out_size[j], 1))
          expect_equal(count_terms_chosen(p[[!!j]]$solution_terms) - 1,
                       out_size[!!j], info = tstsetup)
          expect_identical(p[[!!j]]$family, vss[[!!tstsetup]]$family, info = tstsetup)
          expect_identical(p[[!!j]]$weights, prjdraw_weights, info = tstsetup)
        }
      }
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
