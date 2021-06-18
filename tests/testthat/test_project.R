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
      nprjdraws <- ndraws_default
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
  stopifnot(length(tstsetups) == 1)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    p_fit <- do.call(project, c(
      list(object = fits[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]]),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    expect_identical(prjs[[tstsetup]], p_fit, ignore.environment = TRUE,
                     info = tstsetup)
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
        p_type_expected = TRUE,
        info_str = tstsetup
      )
      tstsetup_tries <- grep(
        paste0("^", args_prj_vs[[tstsetup]]$mod_nm,
               ".", args_prj_vs[[tstsetup]]$fam_nm, ".*\\.clust$"),
        names(prjs), value = TRUE
      )
      match_prj <- sapply(tstsetup_tries, function(tstsetup_try) {
        setequal(solterms_expected_crr, prjs[[tstsetup_try]]$solution_terms)
      })
      tstsetup_match_prj <- tstsetup_tries[match_prj]
      if (length(tstsetup_match_prj) == 1) {
        # cat("Found match:", tstsetup, "and", tstsetup_match_prj, "\n")
        expect_equal(prjs_vs[[tstsetup]], prjs[[tstsetup_match_prj]],
                     info = tstsetup)
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
        p_type_expected = TRUE,
        info_str = tstsetup
      )
      tstsetup_tries <- grep(
        paste0("^", args_prj_cvvs[[tstsetup]]$mod_nm,
               ".", args_prj_cvvs[[tstsetup]]$fam_nm, ".*\\.clust$"),
        names(prjs), value = TRUE
      )
      match_prj <- sapply(tstsetup_tries, function(tstsetup_try) {
        setequal(solterms_expected_crr, prjs[[tstsetup_try]]$solution_terms)
      })
      tstsetup_match_prj <- tstsetup_tries[match_prj]
      if (length(tstsetup_match_prj) == 1) {
        # cat("Found match:", tstsetup, "and", tstsetup_match_prj, "\n")
        expect_equal(prjs_cvvs[[tstsetup]], prjs[[tstsetup_match_prj]],
                     info = tstsetup)
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
                            p_type_expected = !is.null(nclusters_crr),
                            info_str = tstsetup)
        }
      }
    }
  }
})

test_that("specifying `seed` correctly leads to reproducible results", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      tstsetup <- paste(c(mod_nm, fam_nm), collapse = ".")
      tstsetup_prj <- tail(
        grep(paste0("^", mod_nm, "\\.", fam_nm, ".*", "\\.clust$"),
             names(args_prj), value = TRUE),
        1
      )
      ndr_ncl_nm <- intersect(names(args_prj[[tstsetup_prj]]),
                              c("ndraws", "nclusters"))
      stopifnot(identical(ndr_ncl_nm, "nclusters"))
      p_orig <- prjs[[tstsetup_prj]]
      p_new <- project(
        refmods[[mod_nm]][[fam_nm]],
        nclusters = args_prj[[tstsetup_prj]][[ndr_ncl_nm]],
        solution_terms = args_prj[[tstsetup_prj]]$solution_terms,
        seed = args_prj[[tstsetup_prj]]$seed + 1L
      )
      p_repr <- project(
        refmods[[mod_nm]][[fam_nm]],
        nclusters = args_prj[[tstsetup_prj]][[ndr_ncl_nm]],
        solution_terms = args_prj[[tstsetup_prj]]$solution_terms,
        seed = args_prj[[tstsetup_prj]]$seed
      )
      # Expected equality:
      expect_equal(p_orig, p_repr, info = tstsetup)
      # Expected inequality:
      expect_false(isTRUE(all.equal(p_orig, p_new)), info = tstsetup)
    }
  }
})

test_that("for non-GLMs, `regul` has no effect", {
  for (tstsetup in grep("^glm\\.", names(prjs), value = TRUE, invert = TRUE)) {
    args_prj_i <- args_prj[[tstsetup]]
    p_regul <- do.call(project, c(
      list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]],
           regul = 1e-1),
      args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
    ))
    expect_equal(prjs[[tstsetup]], p_regul, info = tstsetup)
  }
})
