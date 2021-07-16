context("project")

test_that(paste(
  "`object` of class \"refmodel\", correctly specified `solution_terms`, and",
  "correctly specified `ndraws` or `nclusters` lead to correct output structure"
), {
  for (tstsetup in names(prjs)) {
    ndr_ncl_nm <- intersect(names(args_prj[[tstsetup]]),
                            c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      ndr_ncl_nm <- "ndraws"
      nprjdraws <- ndraws_pred_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj[[tstsetup]][[ndr_ncl_nm]]
    }
    projection_tester(
      prjs[[tstsetup]],
      solterms_expected = args_prj[[tstsetup]]$solution_terms,
      nprjdraws_expected = nprjdraws,
      p_type_expected = (ndr_ncl_nm == "nclusters" || nprjdraws <= 20),
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
                          p_type_expected = TRUE,
                          info_str = tstsetup_crr)
      }
    }
  }
})

test_that("a fitted model `object` leads to correct output structure", {
  tstsetups <- grep("^glm\\.gauss\\.solterms_x\\.clust", names(prjs),
                    value = TRUE)[1]
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    p_fit <- do.call(project, c(
      list(object = fits[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]]),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    expect_identical(p_fit, prjs[[tstsetup]], ignore.environment = TRUE,
                     info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by varsel()) and correctly specified",
  "`nterms` lead to correct output structure"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs)) {
    tstsetup_vs <- args_prj_vs[[tstsetup]]$tstsetup
    stopifnot(length(tstsetup_vs) > 0)
    mod_crr <- args_vs[[tstsetup_vs]]$mod_nm
    fam_crr <- args_vs[[tstsetup_vs]]$fam_nm
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- vss[[tstsetup_vs]]$suggested_size
    }
    if (length(nterms_crr) == 1) {
      solterms_expected_crr <- vss[[tstsetup_vs]]$solution_terms[
        seq_len(nterms_crr)
      ]
      projection_tester(
        prjs_vs[[tstsetup]],
        solterms_expected = solterms_expected_crr,
        nprjdraws_expected = args_prj_vs[[tstsetup]]$nclusters,
        p_type_expected = TRUE,
        info_str = tstsetup
      )
      # Check that projecting from the "vsel" object onto a single submodel
      # gives the same output as projecting the reference model onto that
      # submodel directly:
      tstsetup_tries <- grep(paste0("^", mod_crr, ".", fam_crr, ".*\\.clust$"),
                             names(prjs), value = TRUE)
      match_prj <- sapply(tstsetup_tries, function(tstsetup_try) {
        setequal(solterms_expected_crr, prjs[[tstsetup_try]]$solution_terms)
      })
      tstsetup_match_prj <- tstsetup_tries[match_prj]
      if (length(tstsetup_match_prj) == 1) {
        expect_identical(prjs_vs[[tstsetup]], prjs[[tstsetup_match_prj]],
                         ignore.environment = TRUE, info = tstsetup)
      } else if (length(tstsetup_match_prj) > 1) {
        stop("Unexpected number of matches.")
      }
    } else {
      proj_list_tester(
        prjs_vs[[tstsetup]],
        len_expected = length(nterms_crr),
        is_seq = all(diff(nterms_crr) == 1),
        info_str = tstsetup,
        nprjdraws_expected = args_prj_vs[[tstsetup]]$nclusters,
        p_type_expected = TRUE,
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
    tstsetup_cvvs <- args_prj_cvvs[[tstsetup]]$tstsetup
    stopifnot(length(tstsetup_cvvs) > 0)
    mod_crr <- args_cvvs[[tstsetup_cvvs]]$mod_nm
    fam_crr <- args_cvvs[[tstsetup_cvvs]]$fam_nm
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- cvvss[[tstsetup_cvvs]]$suggested_size
    }
    if (length(nterms_crr) == 1) {
      solterms_expected_crr <- cvvss[[tstsetup_cvvs]]$solution_terms[
        seq_len(nterms_crr)
      ]
      projection_tester(
        prjs_cvvs[[tstsetup]],
        solterms_expected = solterms_expected_crr,
        nprjdraws_expected = args_prj_cvvs[[tstsetup]]$nclusters,
        p_type_expected = TRUE,
        info_str = tstsetup
      )
      # Check that projecting from the "vsel" object onto a single submodel
      # gives the same output as projecting the reference model onto that
      # submodel directly:
      tstsetup_tries <- grep(paste0("^", mod_crr, ".", fam_crr, ".*\\.clust$"),
                             names(prjs), value = TRUE)
      match_prj <- sapply(tstsetup_tries, function(tstsetup_try) {
        setequal(solterms_expected_crr, prjs[[tstsetup_try]]$solution_terms)
      })
      tstsetup_match_prj <- tstsetup_tries[match_prj]
      if (length(tstsetup_match_prj) == 1) {
        expect_identical(prjs_cvvs[[tstsetup]], prjs[[tstsetup_match_prj]],
                         ignore.environment = TRUE, info = tstsetup)
      } else if (length(tstsetup_match_prj) > 1) {
        stop("Unexpected number of matches.")
      }
    } else {
      proj_list_tester(
        prjs_cvvs[[tstsetup]],
        len_expected = length(nterms_crr),
        is_seq = all(diff(nterms_crr) == 1),
        info_str = tstsetup,
        nprjdraws_expected = args_prj_cvvs[[tstsetup]]$nclusters,
        p_type_expected = TRUE,
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
  skip_if_not(run_vs)
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
                            p_type_expected = !is.null(nclusters_crr),
                            info_str = tstsetup)
        }
      }
    }
  }
})

test_that(paste(
  "specifying `seed` correctly leads to reproducible results (and restores the",
  "RNG state afterwards)"
), {
  # Note: Extensive tests for reproducibility may be found among the tests for
  # .get_refdist().
  for (mod_nm in mod_nms[1]) {
    for (fam_nm in fam_nms[1]) {
      tstsetup <- paste(c(mod_nm, fam_nm), collapse = ".")
      tstsetup_prj <- tail(
        grep(paste0("^", mod_nm, "\\.", fam_nm, ".*", "\\.clust$"),
             names(args_prj), value = TRUE),
        1
      )
      args_prj_i <- args_prj[[tstsetup_prj]]
      ndr_ncl_nm <- intersect(names(args_prj_i), c("ndraws", "nclusters"))
      stopifnot(identical(ndr_ncl_nm, "nclusters"))
      p_orig <- prjs[[tstsetup_prj]]
      rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
      .Random.seed_new1 <- .Random.seed
      p_new <- do.call(project, c(
        list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]],
             seed = args_prj_i$seed + 1L),
        args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm", "seed"))]
      ))
      .Random.seed_new2 <- .Random.seed
      rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
      .Random.seed_repr1 <- .Random.seed
      p_repr <- do.call(project, c(
        list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]]),
        args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
      ))
      .Random.seed_repr2 <- .Random.seed
      # Expected equality:
      expect_equal(p_repr, p_orig, info = tstsetup)
      expect_equal(.Random.seed_new2, .Random.seed_new1, info = tstsetup)
      expect_equal(.Random.seed_repr2, .Random.seed_repr1, info = tstsetup)
      # Expected inequality:
      expect_false(isTRUE(all.equal(p_new, p_orig)), info = tstsetup)
      expect_false(isTRUE(all.equal(rand_new, rand_orig)), info = tstsetup)
      expect_false(isTRUE(all.equal(.Random.seed_repr2, .Random.seed_new2)),
                   info = tstsetup)
    }
  }
})

test_that("for non-GLMs, `regul` has no effect", {
  regul_tst <- 1e-1
  for (mod_crr in setdiff(mod_nms, "glm")) {
    tstsetups <- grep(paste0("^", mod_crr, "\\.gauss\\.solterms_x\\.clust"),
                      names(prjs), value = TRUE)[1]
    stopifnot(length(tstsetups) > 0)
    for (tstsetup in tstsetups) {
      args_prj_i <- args_prj[[tstsetup]]
      p_regul <- do.call(project, c(
        list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]],
             regul = regul_tst),
        args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
      ))
      expect_equal(p_regul, prjs[[tstsetup]], info = tstsetup)
    }
  }
})

test_that("for GLMs, `regul` has an expected effect", {
  regul_tst <- c(regul_default, 1e-1, 1e2)
  stopifnot(regul_tst[1] == regul_default)
  stopifnot(all(diff(regul_tst) > 0))
  tstsetups <- grep("^glm\\..*\\.clust$", names(prjs), value = TRUE)
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    ndr_ncl_nm <- intersect(names(args_prj_i), c("ndraws", "nclusters"))
    if (length(ndr_ncl_nm) == 0) {
      ndr_ncl_nm <- "ndraws"
      nprjdraws <- ndraws_pred_default
    } else {
      stopifnot(length(ndr_ncl_nm) == 1)
      nprjdraws <- args_prj_i[[ndr_ncl_nm]]
    }
    ssq_regul_alpha <- rep(NA, length(regul_tst))
    ssq_regul_beta <- rep(NA, length(regul_tst))
    for (j in seq_along(regul_tst)) {
      # Run project() if necessary:
      if (regul_tst[j] == regul_default) {
        prj_regul <- prjs[[tstsetup]]
      } else {
        prj_regul <- do.call(project, c(
          list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]],
               regul = regul_tst[j]),
          args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
        ))
        projection_tester(
          prj_regul,
          solterms_expected = args_prj_i$solution_terms,
          nprjdraws_expected = nprjdraws,
          p_type_expected = (ndr_ncl_nm == "nclusters" || nprjdraws <= 20),
          info_str = tstsetup
        )
      }

      # Run as.matrix.projection():
      if (ndr_ncl_nm == "nclusters" || nprjdraws <= 20) {
        # Clustered projection, so we expect a warning:
        warn_prjmat_expect <- "the clusters might have different weights"
      } else {
        warn_prjmat_expect <- NA
      }
      expect_warning(prjmat <- as.matrix(prj_regul),
                     warn_prjmat_expect, info = tstsetup)

      # Reduce to only those columns which are necessary here:
      prjmat <- prjmat[, grep("^b_", colnames(prjmat)), drop = FALSE]

      # Posterior means:
      prjmat_mean <- colMeans(prjmat)

      # Calculate the Euclidean norm for intercept and coefficients:
      ssq_regul_alpha[j] <- prjmat_mean["b_Intercept"]^2
      coef_colnms <- setdiff(names(prjmat_mean), "b_Intercept")
      if (length(coef_colnms) > 0) {
        ssq_regul_beta[j] <- sum(prjmat_mean[coef_colnms]^2)
      }
    }
    if (length(args_prj_i$solution_terms) == 0) {
      # For an intercept-only model:
      expect_length(unique(ssq_regul_alpha), 1)
      stopifnot(all(is.na(ssq_regul_beta)))
    } else {
      # All other (i.e., not intercept-only) models:
      for (j in seq_along(ssq_regul_alpha)[-1]) {
        expect_equal(ssq_regul_alpha[!!j], ssq_regul_alpha[j - 1],
                     tolerance = 1e-1, info = tstsetup)
      }
      for (j in seq_along(ssq_regul_beta)[-1]) {
        expect_lt(ssq_regul_beta[!!j], ssq_regul_beta[j - 1])
      }
    }
  }
})
