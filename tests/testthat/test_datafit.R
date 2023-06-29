context("datafit")

# Setup -------------------------------------------------------------------

# Note: Since PR #351, offsets are not supported anymore for `datafit`s. Here,
# we use the data as generated in `setup.R`, i.e., sometimes with offsets. So
# far, this doesn't seem to cause problems in the submodel fitting routines, but
# it might do in the future. Then, either the scenarios including offsets have
# to be excluded or new data has to be generated without offsets.

.extrmoddat_datafit <- function(object, newdata = NULL, wrhs = NULL,
                                orhs = NULL, resp_form = NULL) {
  if (is.null(newdata)) {
    newdata <- object$data
  }

  if (inherits(wrhs, "formula")) {
    weights <- eval_rhs(wrhs, newdata)
  } else if (is.null(wrhs)) {
    weights <- newdata$wobs_col
  } else {
    weights <- wrhs
  }

  offset <- rep(0, nrow(newdata))

  if (inherits(resp_form, "formula")) {
    y <- eval_el2(resp_form, newdata)
  } else {
    y <- NULL
  }

  return(nlist(y, weights, offset))
}

## Reference model --------------------------------------------------------
## (actually "datafit"s)

# Exclude brms reference models (since `datafit`s don't make use of a reference
# model fit, it doesn't make a difference if rstanarm or brms is used as the
# basis here for retrieving the formula, data, and family):
args_datafit <- lapply(setNames(
  nm = grep("^brms\\.", names(refmods), value = TRUE, invert = TRUE)
), function(tstsetup_ref) {
  tstsetup_fit <- args_ref[[tstsetup_ref]]$tstsetup_fit
  c(nlist(tstsetup_fit), only_nonargs(args_fit[[tstsetup_fit]]))
})
# The latent and the augmented-data projection are not supported for `datafit`s:
args_datafit <- args_datafit[!grepl("\\.(latent|augdat)", names(args_datafit))]

datafits <- lapply(args_datafit, function(args_datafit_i) {
  formul_crr <- args_fit[[args_datafit_i$tstsetup_fit]]$formula
  formul_crr <- rm_addresp(formul_crr)
  if (!is.null(args_fit[[args_datafit_i$tstsetup_fit]]$random)) {
    formul_crr <- update(formul_crr, paste(
      ". ~ . + ",
      tail(as.character(args_fit[[args_datafit_i$tstsetup_fit]]$random), 1)
    ))
  }
  extrmoddat <- function(object, newdata = NULL, wrhs = NULL, orhs = NULL,
                         extract_y = TRUE) {
    resp_form <- if (!extract_y) NULL else lhs(formul_crr)
    if (is.null(newdata)) {
      newdata <- dat
    }
    if (args_datafit_i$fam_nm == "brnll") {
      newdata$wobs_col <- 1
    }
    args <- nlist(object, newdata, wrhs, orhs, resp_form)
    return(do.call(.extrmoddat_datafit, args))
  }
  return(init_refmodel(
    object = NULL,
    data = dat,
    formula = formul_crr,
    family = get(paste0("f_", args_datafit_i$fam_nm)),
    extract_model_data = extrmoddat
  ))
})

## Variable selection -----------------------------------------------------

### varsel() --------------------------------------------------------------

if (run_vs) {
  stopifnot(all(names(args_datafit) %in% names(args_ref)))
  args_vs_datafit <- args_vs[
    sapply(args_vs, "[[", "tstsetup_ref") %in% names(datafits)
  ]
  args_vs_datafit <- lapply(args_vs_datafit, function(args_vs_i) {
    names(args_vs_i)[names(args_vs_i) == "tstsetup_ref"] <- "tstsetup_datafit"
    return(args_vs_i)
  })
  # For `"datafit"`s, we always have 1 cluster by default, so omit related
  # arguments:
  args_vs_datafit <- lapply(args_vs_datafit, function(args_vs_i) {
    return(args_vs_i[setdiff(names(args_vs_i),
                             c("ndraws", "nclusters",
                               "ndraws_pred", "nclusters_pred"))])
  })

  vss_datafit <- lapply(args_vs_datafit, function(args_vs_i) {
    do.call(varsel, c(
      list(object = datafits[[args_vs_i$tstsetup_datafit]]),
      excl_nonargs(args_vs_i)
    ))
  })
}

### cv_varsel() -----------------------------------------------------------

if (run_cvvs) {
  args_cvvs_datafit <- args_cvvs[
    sapply(args_cvvs, "[[", "tstsetup_ref") %in% names(datafits)
  ]
  args_cvvs_datafit <- lapply(args_cvvs_datafit, function(args_cvvs_i) {
    names(args_cvvs_i)[names(args_cvvs_i) == "tstsetup_ref"] <-
      "tstsetup_datafit"
    return(args_cvvs_i)
  })
  # (PSIS-)LOO CV is not possible for `"datafit"`s, so only use K-fold CV:
  args_cvvs_datafit <- lapply(args_cvvs_datafit, function(args_cvvs_i) {
    args_cvvs_i$cv_method <- NULL
    args_cvvs_i$K <- NULL
    args_cvvs_i$validate_search <- TRUE
    return(c(args_cvvs_i, list(cv_method = "kfold", K = K_tst)))
  })
  names(args_cvvs_datafit) <- gsub("default_cvmeth", "kfold",
                                   names(args_cvvs_datafit))
  args_cvvs_datafit <- args_cvvs_datafit[unique(names(args_cvvs_datafit))]
  # For `"datafit"`s, we always have 1 cluster by default, so omit related
  # arguments:
  args_cvvs_datafit <- lapply(args_cvvs_datafit, function(args_cvvs_i) {
    return(args_cvvs_i[setdiff(names(args_cvvs_i),
                               c("ndraws", "nclusters",
                                 "ndraws_pred", "nclusters_pred"))])
  })

  cvvss_datafit <- lapply(args_cvvs_datafit, function(args_cvvs_i) {
    do.call(cv_varsel, c(
      list(object = datafits[[args_cvvs_i$tstsetup_datafit]]),
      excl_nonargs(args_cvvs_i)
    ))
  })
}

## Projection -------------------------------------------------------------

### From varsel() ---------------------------------------------------------

if (run_vs) {
  args_prj_vs_datafit <- args_prj_vs[
    sapply(args_prj_vs, "[[", "tstsetup_ref") %in% names(datafits) &
      sapply(args_prj_vs, "[[", "tstsetup_vsel") %in% names(vss_datafit)
  ]
  args_prj_vs_datafit <- lapply(args_prj_vs_datafit, function(args_prj_vs_i) {
    names(args_prj_vs_i)[names(args_prj_vs_i) == "tstsetup_ref"] <-
      "tstsetup_datafit"
    return(args_prj_vs_i)
  })
  # For `"datafit"`s, we always have 1 cluster by default, so omit related
  # arguments:
  args_prj_vs_datafit <- lapply(args_prj_vs_datafit, function(args_prj_vs_i) {
    return(args_prj_vs_i[setdiff(names(args_prj_vs_i),
                                 c("ndraws", "nclusters"))])
  })

  prjs_vs_datafit <- lapply(args_prj_vs_datafit, function(args_prj_vs_i) {
    args_prj_vs_i$refit_prj <- FALSE
    do.call(project, c(
      list(object = vss_datafit[[args_prj_vs_i$tstsetup_vsel]]),
      excl_nonargs(args_prj_vs_i)
    ))
  })
}

## Prediction -------------------------------------------------------------

### From "proj_list" ------------------------------------------------------

if (run_vs) {
  pls_vs_datafit <- lapply(prjs_vs_datafit, proj_linpred, .seed = seed2_tst)
  pps_vs_datafit <- lapply(prjs_vs_datafit, proj_predict, .seed = seed2_tst)
}

# Tests (projpred only) ---------------------------------------------------

## Reference model --------------------------------------------------------

test_that("init_refmodel(): `object` of class \"datafit\" works", {
  for (tstsetup in names(datafits)) {
    tstsetup_fit <- args_datafit[[tstsetup]]$tstsetup_fit
    with_spclformul_crr <- grepl("\\.spclformul", tstsetup)
    if (args_datafit[[tstsetup]]$fam_nm == "binom" ||
        grepl("\\.with_wobs", tstsetup)) {
      wobs_expected_crr <- wobs_tst
    } else {
      wobs_expected_crr <- rep(1, nobsv)
    }
    refmodel_tester(
      datafits[[tstsetup]],
      is_datafit = TRUE,
      pkg_nm = args_datafit[[tstsetup]]$pkg_nm,
      fit_expected = NULL,
      formul_expected = get_formul_from_fit(fits[[tstsetup_fit]]),
      with_spclformul = with_spclformul_crr,
      wobs_expected = wobs_expected_crr,
      offs_expected = rep(0, nobsv),
      nrefdraws_expected = 1L,
      fam_orig = get(paste0("f_", args_datafit[[tstsetup]]$fam_nm)),
      mod_nm = args_datafit[[tstsetup]]$mod_nm,
      fam_nm = args_datafit[[tstsetup]]$fam_nm,
      info_str = tstsetup
    )
  }
})

test_that("predict.refmodel(): `object` of class \"datafit\" fails", {
  for (tstsetup in names(datafits)) {
    expect_error(
      predict(datafits[[tstsetup]], newdata = dat),
      "^Cannot make predictions for an `object` of class \"datafit\"\\.$",
      info = tstsetup
    )
  }
})

## Variable selection -----------------------------------------------------

test_that(paste(
  "varsel(): `object` of class \"datafit\", `method`, and `nterms_max` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(vss_datafit)) {
    mod_crr <- args_vs_datafit[[tstsetup]]$mod_nm
    meth_exp_crr <- args_vs_datafit[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    extra_tol_crr <- 1.5
    if (any(grepl(":", ranking(vss_datafit[[tstsetup]])[["fulldata"]]))) {
      ### Testing for non-increasing element `ce` (for increasing model size)
      ### doesn't make sense if the ranking of predictors involved in
      ### interactions has been changed, so we choose a higher `extra_tol`:
      extra_tol_crr <- 3
      ###
    }
    vsel_tester(
      vss_datafit[[tstsetup]],
      from_datafit = TRUE,
      refmod_expected =
        datafits[[args_vs_datafit[[tstsetup]]$tstsetup_datafit]],
      solterms_len_expected = args_vs_datafit[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      search_trms_empty_size =
        length(args_vs_datafit[[tstsetup]]$search_terms) &&
        all(grepl("\\+", args_vs_datafit[[tstsetup]]$search_terms)),
      extra_tol = extra_tol_crr,
      info_str = tstsetup
    )
  }
})

test_that(paste(
  "cv_varsel(): `object` of class \"datafit\", `method`, `cv_method`, and",
  "`nterms_max` work"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(cvvss_datafit)) {
    mod_crr <- args_cvvs_datafit[[tstsetup]]$mod_nm
    meth_exp_crr <- args_cvvs_datafit[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    vsel_tester(
      cvvss_datafit[[tstsetup]],
      with_cv = TRUE,
      from_datafit = TRUE,
      refmod_expected =
        datafits[[args_cvvs_datafit[[tstsetup]]$tstsetup_datafit]],
      solterms_len_expected = args_cvvs_datafit[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "kfold",
      valsearch_expected = args_cvvs_datafit[[tstsetup]]$validate_search,
      search_trms_empty_size =
        length(args_cvvs_datafit[[tstsetup]]$search_terms) &&
        all(grepl("\\+", args_cvvs_datafit[[tstsetup]]$search_terms)),
      extra_tol = 1.2,
      info_str = tstsetup
    )
  }
})

## Projection -------------------------------------------------------------

test_that("project(): `object` of class \"datafit\" fails", {
  skip_if_not(run_prj)
  # A prerequisite for this project() test (otherwise, it would have to be
  # adapted):
  stopifnot(all(names(args_datafit) %in% names(args_ref)))

  tstsetups <- grep("\\.solterms_x.*\\.clust$", names(args_prj), value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    if (!args_prj_i$tstsetup_ref %in% names(datafits)) next
    args_prj_i$refit_prj <- FALSE
    expect_error(
      do.call(project, c(
        list(object = datafits[[args_prj_i$tstsetup_ref]]),
        excl_nonargs(args_prj_i)
      )),
      paste("^project\\(\\) does not support an `object` of class",
            "\"datafit\"\\.$"),
      info = tstsetup
    )
  }
})

test_that(paste(
  "project(): `object` of class \"vsel\" (created by varsel() applied to an",
  "`object` of class \"datafit\"), `nclusters`, and `nterms` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs_datafit)) {
    tstsetup_vs <- args_prj_vs_datafit[[tstsetup]]$tstsetup_vsel
    stopifnot(length(tstsetup_vs) > 0)
    nterms_crr <- args_prj_vs_datafit[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- suggest_size(vss_datafit[[tstsetup_vs]], warnings = FALSE)
    }
    with_L1 <- (args_vs_datafit[[tstsetup_vs]]$mod_nm == "glm" &&
                  is.null(args_vs_datafit[[tstsetup_vs]]$method)) ||
      identical(args_vs_datafit[[tstsetup_vs]]$method, "L1")
    if (length(nterms_crr) == 1) {
      solterms_expected_crr <- vss_datafit[[tstsetup_vs]]$solution_terms[
        seq_len(nterms_crr)
      ]
      projection_tester(
        prjs_vs_datafit[[tstsetup]],
        refmod_expected =
          datafits[[args_prj_vs_datafit[[tstsetup]]$tstsetup_datafit]],
        solterms_expected = solterms_expected_crr,
        nprjdraws_expected = 1L,
        with_clusters = TRUE,
        const_wdraws_prj_expected = TRUE,
        from_vsel_L1_search = with_L1,
        info_str = tstsetup
      )
      ### TODO: Currently, the as.matrix() call below is not possible for
      ### `datafit`s. Fix this.
      # if (run_snaps) {
      #   if (testthat_ed_max2) local_edition(3)
      #   suppressWarnings(m <- as.matrix(prjs_vs_datafit[[tstsetup]]))
      #   expect_snapshot({
      #     print(tstsetup)
      #     print(rlang::hash(m)) # cat(m)
      #   })
      #   if (testthat_ed_max2) local_edition(2)
      # }
      ###
    } else {
      proj_list_tester(
        prjs_vs_datafit[[tstsetup]],
        len_expected = length(nterms_crr),
        is_seq = all(diff(nterms_crr) == 1),
        info_str = tstsetup,
        refmod_expected =
          datafits[[args_prj_vs_datafit[[tstsetup]]$tstsetup_datafit]],
        nprjdraws_expected = 1L,
        with_clusters = TRUE,
        const_wdraws_prj_expected = TRUE,
        prjdraw_weights_expected = prjs_vs_datafit[[tstsetup]][[1]]$wdraws_prj,
        from_vsel_L1_search = with_L1
      )
      ### TODO: Currently, the as.matrix() call below is not possible for
      ### `datafit`s. Fix this.
      # if (run_snaps) {
      #   if (testthat_ed_max2) local_edition(3)
      #   res_vs <- lapply(prjs_vs_datafit[[tstsetup]], function(prjs_vs_i) {
      #     suppressWarnings(m <- as.matrix(prjs_vs_i))
      #     expect_snapshot({
      #       print(tstsetup)
      #       print(prjs_vs_i$solution_terms)
      #       print(rlang::hash(m)) # cat(m)
      #     })
      #     return(invisible(TRUE))
      #   })
      #   if (testthat_ed_max2) local_edition(2)
      # }
      ###
    }
  }
})

## Prediction -------------------------------------------------------------

test_that(paste(
  "proj_linpred(): `object` of (informal) class \"proj_list\" (based on",
  "varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs_datafit)) {
    tstsetup_vs <- args_prj_vs_datafit[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_vs_datafit[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- suggest_size(vss_datafit[[tstsetup_vs]], warnings = FALSE)
    }
    pl_tester(pls_vs_datafit[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_expected = 1L,
              info_str = tstsetup)
    pl_with_args <- proj_linpred(
      prjs_vs_datafit[[tstsetup]],
      newdata = head(
        get_dat_formul(
          args_fit[[args_prj_vs_datafit[[tstsetup]]$tstsetup_fit]]$formula,
          needs_adj = grepl("\\.spclformul", tstsetup)
        ),
        tail(nobsv_tst, 1)
      ),
      weightsnew = ~ wobs_col,
      filter_nterms = nterms_crr[1],
      .seed = seed2_tst
    )
    pl_tester(pl_with_args,
              len_expected = 1L,
              nprjdraws_expected = 1L,
              nobsv_expected = tail(nobsv_tst, 1),
              info_str = paste(tstsetup, "with_args", sep = "__"))
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(pl_with_args))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "proj_predict(): `object` of (informal) class \"proj_list\" (based on",
  "varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs_datafit)) {
    tstsetup_vs <- args_prj_vs_datafit[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_vs_datafit[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- suggest_size(vss_datafit[[tstsetup_vs]], warnings = FALSE)
    }
    pp_tester(pps_vs_datafit[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_out_expected = 1L,
              info_str = tstsetup)
    pp_with_args <- proj_predict(
      prjs_vs_datafit[[tstsetup]],
      newdata = head(dat, tail(nobsv_tst, 1)),
      weightsnew = ~ wobs_col,
      filter_nterms = nterms_crr[1],
      nresample_clusters = tail(nresample_clusters_tst, 1),
      .seed = seed2_tst
    )
    pp_tester(pp_with_args,
              len_expected = 1L,
              nprjdraws_out_expected = 1L,
              nobsv_expected = tail(nobsv_tst, 1),
              info_str = paste(tstsetup, "with_args", sep = "__"))
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(pp_with_args))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

## summary.vsel() ---------------------------------------------------------

test_that("summary.vsel(): `object` of class \"datafit\" fails", {
  for (tstsetup in names(datafits)) {
    expect_error(
      summary.vsel(datafits[[tstsetup]]),
      paste("^The object is not a variable selection object\\. Run variable",
            "selection first$"),
      info = tstsetup
    )
  }
})

test_that("summary.vsel(): `baseline = \"ref\"` and `deltas = TRUE` fails", {
  skip_if_not(run_vs)
  for (tstsetup in head(names(vss_datafit), 1)) {
    expect_error(
      summary(vss_datafit[[tstsetup]], baseline = "ref", deltas = TRUE),
      paste("^Cannot use deltas = TRUE and baseline = 'ref' when there is no",
            "reference model\\.$"),
      info = tstsetup
    )
  }
})

test_that(paste(
  "summary.vsel(): `object` of class \"vsel\" (created by varsel() applied to",
  "an `object` of class \"datafit\"), `stats`, and `type` work"
), {
  skip_if_not(run_vs)
  tstsetups <- unname(unlist(lapply(mod_nms, function(mod_nm) {
    unlist(lapply(fam_nms, function(fam_nm) {
      grep(paste0("\\.", mod_nm, "\\.", fam_nm), names(vss_datafit),
           value = TRUE)
    }))
  })))
  for (tstsetup in tstsetups) {
    smmry <- summary(vss_datafit[[tstsetup]],
                     stats = stats_common,
                     type = type_tst,
                     seed = seed3_tst)
    smmry_tester(
      smmry,
      vsel_expected = vss_datafit[[tstsetup]],
      search_trms_empty_size =
        length(args_vs_datafit[[tstsetup]]$search_terms) &&
        all(grepl("\\+", args_vs_datafit[[tstsetup]]$search_terms)),
      info_str = tstsetup,
      stats_expected = stats_common,
      type_expected = type_tst,
      solterms_expected = vss_datafit[[tstsetup]]$solution_terms
    )
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(smmry, digits = 6)
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "summary.vsel(): `object` of class \"vsel\" (created by cv_varsel() applied",
  "to an `object` of class \"datafit\"), `stats`, and `type` work"
), {
  skip_if_not(run_cvvs)
  tstsetups <- unname(unlist(lapply(mod_nms, function(mod_nm) {
    unlist(lapply(fam_nms, function(fam_nm) {
      grep(paste0("\\.", mod_nm, "\\.", fam_nm), names(cvvss_datafit),
           value = TRUE)
    }))
  })))
  for (tstsetup in tstsetups) {
    smmry <- summary(cvvss_datafit[[tstsetup]],
                     stats = stats_common,
                     type = type_tst,
                     seed = seed3_tst)
    smmry_tester(
      smmry,
      vsel_expected = cvvss_datafit[[tstsetup]],
      search_trms_empty_size =
        length(args_cvvs_datafit[[tstsetup]]$search_terms) &&
        all(grepl("\\+", args_cvvs_datafit[[tstsetup]]$search_terms)),
      info_str = tstsetup,
      stats_expected = stats_common,
      type_expected = type_tst,
      cv_method_expected =
        args_cvvs_datafit[[tstsetup]]$cv_method %||% "LOO",
      solterms_expected = cvvss_datafit[[tstsetup]]$solution_terms
    )
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(smmry, digits = 6)
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

# Comparison with glmnet --------------------------------------------------

# below are some tests that check Lasso solution computed with varsel is the
# same as that of glmnet. (notice that glm_ridge and glm_elnet are already
# tested separately, so these would only check that the results do not change
# due to varsel/cv_varsel etc.)

test_that(paste(
  "L1-projection with data reference gives the same results as",
  "Lasso from glmnet."
), {
  skip_if_not_installed("glmnet")
  # This test sometimes behaves inpredictably when run in `R CMD check`, so skip
  # it on CRAN:
  skip_on_cran()
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1235)
  n <- 100
  nterms <- 10
  x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
  b <- seq(0, 1, length.out = nterms)
  dis <- runif(1, 0.3, 0.5)
  weights <- sample(1:4, n, replace = TRUE)

  fams <- list(gaussian(), binomial(), poisson())
  x_list <- lapply(fams, function(fam) x)
  y_list <- lapply(fams, function(fam) {
    if (fam$family == "gaussian") {
      y <- rnorm(n, x %*% b, 0.5)
      weights <- NULL
      y_glmnet <- y
    } else if (fam$family == "binomial") {
      y <- rbinom(n, weights, fam$linkinv(x %*% b))
      ## y <- y / weights
      ## different way of specifying binomial y for glmnet
      y_glmnet <- cbind(1 - y / weights, y / weights)
      weights <- weights
    } else if (fam$family == "poisson") {
      y <- rpois(n, fam$linkinv(x %*% b))
      y_glmnet <- y
      weights <- NULL
    }
    nlist(y, y_glmnet, weights)
  })

  extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                 orhs = NULL, extract_y = FALSE) {
    if (!is.null(object)) {
      formula <- formula(object)
      tt <- extract_terms_response(formula)
      response_name <- tt$response
    } else {
      response_name <- NULL
    }

    if (is.null(newdata)) {
      newdata <- object$data
    }

    resp_form <- NULL
    if (is.null(object)) {
      if ("weights" %in% colnames(newdata)) {
        wrhs <- ~ weights
      }
      if ("y" %in% colnames(newdata)) {
        resp_form <- ~ y
      }
    }

    args <- nlist(object, newdata, wrhs, orhs, resp_form)
    return(do_call(.extract_model_data, args))
  }

  for (i in seq_along(fams)) {
    x <- x_list[[i]]
    y <- y_list[[i]]$y
    y_glmnet <- y_list[[i]]$y_glmnet
    fam <- fams[[i]]
    weights <- y_list[[i]]$weights
    if (is.null(weights)) {
      weights <- rep(1, NROW(y))
    }

    lambda_min_ratio <- 1e-7
    nlambda <- 1500

    df <- data.frame(y = y, x = x, weights = weights)
    formula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10
    # Lasso solution with projpred
    ref <- init_refmodel(
      object = NULL, data = df, formula = formula,
      family = fam, extract_model_data = extract_model_data
    )
    vs <- suppressWarnings(varsel(
      ref,
      method = "L1", lambda_min_ratio = lambda_min_ratio,
      nlambda = nlambda, thresh = 1e-12, verbose = FALSE
    ))
    pred1 <- proj_linpred(vs,
                          newdata = data.frame(x = x, weights = weights),
                          transform = FALSE, .seed = seed2_tst,
                          nterms = 0:nterms, refit_prj = FALSE)

    # compute the results for the Lasso
    lasso <- glmnet::glmnet(x, y_glmnet,
                            family = fam$family, weights = weights,
                            lambda.min.ratio = lambda_min_ratio,
                            nlambda = nlambda, thresh = 1e-12)
    solution_terms <- predict(lasso, type = "nonzero", s = lasso$lambda)
    nselected <- sapply(solution_terms, function(e) length(e))
    lambdainds <- sapply(unique(nselected), function(nterms) {
      max(which(nselected == nterms))
    })
    lambdaval <- lasso$lambda[lambdainds]
    pred2 <- predict(lasso, newx = x, type = "link", s = lambdaval)
    solution_terms_lasso <- integer()
    lasso_coefs <- as.matrix(lasso$beta[, tail(lambdainds, -1), drop = FALSE])
    for (ii in seq_len(ncol(lasso_coefs))) {
      solution_terms_lasso <- c(
        solution_terms_lasso,
        setdiff(which(lasso_coefs[, ii] > 0), solution_terms_lasso)
      )
    }

    # check that the predictions agree
    pred1 <- unname(sapply(pred1, "[[", "pred"))
    pred2 <- unname(pred2)
    # Sometimes, glmnet terminates the coefficient path computation too early
    # for some reason:
    if (ncol(pred1) > ncol(pred2)) {
      pred1 <- pred1[, seq_len(ncol(pred2)), drop = FALSE]
    }
    expect_equal(pred1, pred2, tolerance = 1e-2, info = as.character(i))

    # check that the coefficients are similar
    ind <- match(vs$solution_terms, setdiff(split_formula(formula), "1"))
    betas <- sapply(vs$search_path$outdmins, function(x) x[[1]]$beta %||% 0)
    delta <- sapply(seq_len(length(lambdainds) - 1), function(i) {
      abs(t(betas[[i + 1]]) - lasso$beta[ind[1:i], lambdainds[i + 1]])
    })
    expect_true(median(unlist(delta)) < 6e-2)
    expect_true(median(abs(
      sapply(head(vs$search_path$outdmins, length(lambdainds)),
             function(x) {
               x[[1]]$alpha
             }) -
        lasso$a0[lambdainds])
    ) < 1.5e-1)
    # Sometimes, glmnet terminates the coefficient path computation too early
    # for some reason:
    if (length(ind) > length(solution_terms_lasso)) {
      ind <- ind[seq_along(solution_terms_lasso)]
    }
    expect_identical(ind, solution_terms_lasso, info = as.character(i))
  }
  RNGversion(getRversion())
  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})
