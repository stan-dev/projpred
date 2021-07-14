# varsel() ----------------------------------------------------------------

context("varsel()")

test_that(paste(
  "`object` of class \"refmodel\", correctly specified `method`, `nterms_max`,",
  "`nclusters`, and `nclusters_pred` lead to correct output structure"
), {
  skip_if_not(exists("vss"))
  for (tstsetup in names(vss)) {
    mod_crr <- args_vs[[tstsetup]]$mod
    fam_crr <- args_vs[[tstsetup]]$fam
    meth_exp_crr <- args_vs[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    vsel_tester(
      vss[[tstsetup]],
      refmod_expected = refmods[[mod_crr]][[fam_crr]],
      solterms_len_expected = args_vs[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      nclusters_expected = args_vs[[tstsetup]]$nclusters,
      nclusters_pred_expected = args_vs[[tstsetup]]$nclusters_pred,
      info_str = tstsetup
    )
  }
})

test_that("specifying `object` incorrectly leads to an error", {
  expect_error(varsel(rnorm(5), verbose = FALSE),
               "no applicable method")
})

test_that("specifying `method` incorrectly leads to an error", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      expect_error(varsel(refmods[[!!mod_nm]][[!!fam_nm]], method = "k-fold"),
                   "Unknown search method")
      if (mod_nm == "glmm") {
        expect_error(varsel(refmods[[!!mod_nm]][[!!fam_nm]], method = "L1"),
                     "^L1 search is not supported for multilevel models\\.$")
      } else if (!mod_nm %in% c("glm", "glmm")) {
        ### TODO:
        stop("Still to-do.")
        # expect_error(varsel(refmods[[!!mod_nm]][[!!fam_nm]], method = "L1"),
        #              "ENTER EXPECTED TEXT HERE")
        ###
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
  skip_if_not(exists("vss"))
  for (tstsetup in names(vss)[1]) {
    args_vs_i <- args_vs[[tstsetup]]
    vs_orig <- vss[[tstsetup]]
    rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_new1 <- .Random.seed
    vs_new <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
           seed = args_vs_i$seed + 1L),
      args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm", "seed"))]
    ))
    .Random.seed_new2 <- .Random.seed
    rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_repr1 <- .Random.seed
    vs_repr <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]]),
      args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
    ))
    .Random.seed_repr2 <- .Random.seed
    # Expected equality:
    expect_equal(vs_repr, vs_orig, info = tstsetup)
    expect_equal(.Random.seed_new2, .Random.seed_new1, info = tstsetup)
    expect_equal(.Random.seed_repr2, .Random.seed_repr1, info = tstsetup)
    # Expected inequality:
    expect_false(isTRUE(all.equal(vs_new, vs_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(rand_new, rand_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_repr2, .Random.seed_new2)),
                 info = tstsetup)
  }
})

test_that("specifying d_test has an expected effect", {
  skip_if_not(exists("vss"))
  tstsetups <- grep("^glm\\.gauss", names(vss), value = TRUE)[1]
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    mod_crr <- args_vs_i$mod_nm
    fam_crr <- args_vs_i$fam_nm
    refmod_crr <- refmods[[mod_crr]][[fam_crr]]
    d_test_crr <- list(
      y = refmod_crr$y,
      test_points = seq_along(refmod_crr$y),
      data = refmod_crr$fit$data,
      weights = refmod_crr$wobs,
      type = "test",
      offset = refmod_crr$offset
    )
    # We expect a warning which in fact should be suppressed, though (see
    # issue #162):
    expect_warning(
      vs_repr <- do.call(varsel, c(
        list(object = refmod_crr, d_test = d_test_crr),
        args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
      )),
      paste("^'offset' argument is NULL but it looks like you estimated the",
            "model using an offset term\\.$")
    )
    meth_exp_crr <- args_vs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    vsel_tester(
      vs_repr,
      refmod_expected = refmod_crr,
      dtest_expected = d_test_crr,
      solterms_len_expected = args_vs_i$nterms_max,
      method_expected = meth_exp_crr,
      nclusters_expected = args_vs_i$nclusters,
      nclusters_pred_expected = args_vs_i$nclusters_pred,
      info_str = tstsetup
    )
    expect_identical(vs_repr$d_test, d_test_crr, info = tstsetup)
    expect_equal(vs_repr[setdiff(names(vs_repr), vsel_nms_dtest)],
                 vss[[tstsetup]][setdiff(names(vss[[tstsetup]]),
                                         vsel_nms_dtest)],
                 info = tstsetup)
  }
})

## Regularization ---------------------------------------------------------

# In fact, `regul` is already checked in `test_project.R`, so the `regul` tests
# could be omitted here since varsel() and cv_varsel() also pass `regul` to
# project_submodel() (usually via .get_submodels(), just like project()). This
# doesn't hold for L1 search, though. So for L1 search, the `regul` tests are
# still needed.

test_that("for non-GLMs, `regul` has no effect", {
  skip_if_not(exists("vss"))
  regul_tst <- 1e-1
  for (mod_crr in setdiff(mod_nms, "glm")) {
    tstsetups <- grep(paste0("^", mod_crr, "\\.gauss"), names(vss),
                      value = TRUE)[1]
    stopifnot(length(tstsetups) > 0)
    for (tstsetup in tstsetups) {
      args_vs_i <- args_vs[[tstsetup]]
      vs_regul <- do.call(varsel, c(
        list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
             regul = regul_tst),
        args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
      ))
      expect_equal(vs_regul, vss[[tstsetup]], info = tstsetup)
    }
  }
})

test_that(paste(
  "for GLMs with L1 search, `regul` only has an effect on prediction, not on",
  "selection"
), {
  skip_if_not(exists("vss"))
  regul_tst <- c(regul_default, 1e-1, 1e2)
  stopifnot(regul_tst[1] == regul_default)
  stopifnot(all(diff(regul_tst) > 0))
  tstsetups <- setdiff(grep("^glm\\.", names(vss), value = TRUE),
                       grep("^glm\\..*\\.forward", names(vss), value = TRUE))
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    m_max <- args_vs_i$nterms_max + 1L
    ssq_regul_prd <- array(dim = c(length(regul_tst), m_max))
    for (j in seq_along(regul_tst)) {
      if (regul_tst[j] == regul_default) {
        vs_regul <- vss[[tstsetup]]
      } else {
        vs_regul <- do.call(varsel, c(
          list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
               regul = regul_tst[j]),
          args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
        ))
        vsel_tester(
          vs_regul,
          refmod_expected = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
          solterms_len_expected = args_vs_i$nterms_max,
          method_expected = "L1",
          nclusters_expected = args_vs_i$nclusters,
          nclusters_pred_expected = args_vs_i$nclusters_pred,
          info_str = tstsetup
        )
        # Expect equality for all components not related to prediction:
        expect_equal(vs_regul[vsel_nms_nonpred],
                     vss[[tstsetup]][vsel_nms_nonpred],
                     info = paste(tstsetup, j, sep = "__"))
        ### Excluded for the sake of speed (and because the inequality of the
        ### prediction components is checked below in detail):
        # # Expect inequality when taking only the components related to
        # # prediction:
        # expect_false(isTRUE(all.equal(vs_regul[vsel_nms_pred],
        #                               vss[[tstsetup]][vsel_nms_pred])),
        #              info = paste(tstsetup, j, sep = "__"))
        ###
      }
      # Check the inequality of the prediction components in detail: Expect a
      # reduction of the sum of the squared coefficients (excluding the
      # intercept) for increasing `regul`:
      for (m in seq_len(m_max)) {
        # Since varsel() doesn't output object `p_sub`, use the linear predictor
        # `$summaries$sub[[m]]$mu` here (instead of the coefficients themselves,
        # which would only be accessible from `p_sub`):
        mu_jm_regul <- vs_regul$summaries$sub[[m]]$mu
        # In fact, `sum((mu - intercept)^2)` would make more sense than
        # `var(mu) = sum((mu - mean(mu))^2)` but since varsel() doesn't output
        # object `p_sub`, the intercept from the prediction is not accessible
        # here.
        ssq_regul_prd[j, m] <- var(mu_jm_regul)
      }
    }
    # For the intercept-only model, the linear predictor consists only
    # of the intercept, so we expect no variation in `mu_jm_regul`:
    expect_true(all(ssq_regul_prd[, 1] == 0), info = tstsetup)
    # All other (i.e., not intercept-only) models:
    for (j in seq_len(dim(ssq_regul_prd)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_prd)[2])[-1]) {
        expect_lt(ssq_regul_prd[!!j, !!m], ssq_regul_prd[j - 1, m])
      }
    }
  }
})

test_that(paste(
  "for GLMs with forward search, `regul` has an expected effect on selection",
  "as well as on prediction"
), {
  skip_if_not(exists("vss"))
  regul_tst <- c(regul_default, 1e-1, 1e2)
  stopifnot(regul_tst[1] == regul_default)
  stopifnot(all(diff(regul_tst) > 0))
  tstsetups <- grep("^glm\\..*\\.forward", names(vss), value = TRUE)
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    m_max <- args_vs_i$nterms_max + 1L
    if (identical(args_vs_i$method, "forward")) {
      ncl_crr <- args_vs_i$nclusters
    } else {
      ncl_crr <- 1L
    }
    ssq_regul_sel_alpha <- array(dim = c(length(regul_tst), m_max, ncl_crr))
    ssq_regul_sel_beta <- array(dim = c(length(regul_tst), m_max, ncl_crr))
    ssq_regul_prd <- array(dim = c(length(regul_tst), m_max))
    for (j in seq_along(regul_tst)) {
      if (regul_tst[j] == regul_default) {
        vs_regul <- vss[[tstsetup]]
      } else {
        vs_regul <- do.call(varsel, c(
          list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
               regul = regul_tst[j]),
          args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
        ))
        vsel_tester(
          vs_regul,
          refmod_expected = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
          solterms_len_expected = args_vs_i$nterms_max,
          method_expected = "forward",
          nclusters_expected = args_vs_i$nclusters,
          nclusters_pred_expected = args_vs_i$nclusters_pred,
          info_str = tstsetup
        )
      }
      for (m in seq_len(m_max)) {
        # Selection:
        subfits_jm_regul <- vs_regul$search_path$sub_fits[[m]]
        if (ncl_crr == 1) {
          subfits_jm_regul <- list(subfits_jm_regul)
        } else {
          stopifnot(identical(ncl_crr, length(subfits_jm_regul)))
        }
        for (nn in seq_len(ncl_crr)) {
          stopifnot(length(subfits_jm_regul[[nn]]$alpha) == 1)
          ssq_regul_sel_alpha[j, m, nn] <- subfits_jm_regul[[nn]]$alpha^2
          if (length(subfits_jm_regul[[nn]]$beta) > 0) {
            ssq_regul_sel_beta[j, m, nn] <- sum(subfits_jm_regul[[nn]]$beta^2)
          }
        }
        # Prediction:
        # Since varsel() doesn't output object `p_sub`, use the linear predictor
        # `$summaries$sub[[m]]$mu` here (instead of the coefficients themselves,
        # which would only be accessible from `p_sub`):
        mu_jm_regul <- vs_regul$summaries$sub[[m]]$mu
        # In fact, `sum((mu - intercept)^2)` would make more sense than
        # `var(mu) = sum((mu - mean(mu))^2)` but since varsel() doesn't output
        # object `p_sub`, the intercept from the prediction is not accessible
        # here.
        ssq_regul_prd[j, m] <- var(mu_jm_regul)
      }
    }
    # Selection:
    ### Excluded because of issue #169:
    # # For the intercept-only model:
    # for (nn in seq_len(dim(ssq_regul_sel_alpha)[3])) {
    #   expect_length(unique(ssq_regul_sel_alpha[, 1, !!n]), 1)
    # }
    # # All other (i.e., not intercept-only) models:
    # for (j in seq_len(dim(ssq_regul_sel_alpha)[1])[-1]) {
    #   for (m in seq_len(dim(ssq_regul_sel_alpha)[2])[-1]) {
    #     for (nn in seq_len(dim(ssq_regul_sel_alpha)[3])) {
    #       expect_equal(ssq_regul_sel_alpha[!!j, !!m, !!nn],
    #                    ssq_regul_sel_alpha[j - 1, m, nn],
    #                    tolerance = 1e-2)
    #     }
    #   }
    # }
    ###
    # For the intercept-only model:
    expect_true(all(is.na(ssq_regul_sel_beta[, 1, ])), info = tstsetup)
    # All other (i.e., not intercept-only) models:
    for (j in seq_len(dim(ssq_regul_sel_beta)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_sel_beta)[2])[-1]) {
        for (nn in seq_len(dim(ssq_regul_sel_beta)[3])) {
          expect_lt(ssq_regul_sel_beta[!!j, !!m, !!nn],
                    ssq_regul_sel_beta[j - 1, m, nn])
        }
      }
    }
    # Prediction:
    # For the intercept-only model, the linear predictor consists only
    # of the intercept, so we expect no variation in `mu_jm_regul`:
    expect_true(all(ssq_regul_prd[, 1] == 0), info = tstsetup)
    # All other (i.e., not intercept-only) models:
    for (j in seq_len(dim(ssq_regul_prd)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_prd)[2])[-1]) {
        expect_lt(ssq_regul_prd[!!j, !!m], ssq_regul_prd[j - 1, m])
      }
    }
  }
})

## Penalty ----------------------------------------------------------------

test_that("`penalty` of incorrect length causes an error", {
  skip_if_not(exists("args_vs"))
  tstsetups <- setdiff(
    grep("^glm\\.", names(args_vs), value = TRUE),
    grep("^glm\\..*\\.forward", names(args_vs), value = TRUE)
  )
  stopifnot(length(tstsetups) > 0)
  # Note: As mentioned in issue #149, the reference level of a categorical
  # predictor actually should not have its own coefficient:
  len_penal <- sum(grepl("^xco\\.", names(dat))) +
    sum(sapply(dat[, grep("^xca\\.", names(dat))], nlevels))
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    penal_tst <- list(rep(1, args_vs_i$nterms_max + 10),
                      1)
    for (penal_crr in penal_tst) {
      expect_error(
        do.call(varsel, c(
          list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
               penalty = penal_crr),
          args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
        )),
        paste0("^Incorrect length of penalty vector \\(should be ",
               len_penal, "\\)\\.$")
      )
    }
  }
})

test_that("for forward search, `penalty` has no effect", {
  skip_if_not(exists("vss"))
  penal_tst <- 2
  tstsetups <- union(grep("\\.forward", names(vss), value = TRUE),
                     grep("^glm\\.", names(vss), value = TRUE, invert = TRUE))
  tstsetups <- tstsetups[1] # For the sake of speed.
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    vs_penal <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
           penalty = penal_tst),
      args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
    ))
    expect_equal(vs_penal, vss[[tstsetup]], info = tstsetup)
  }
})

test_that("for L1 search, `penalty` has an expected effect", {
  skip_if_not(exists("vss"))
  tstsetups <- setdiff(grep("^glm\\.", names(vss), value = TRUE),
                       grep("^glm\\..*\\.forward", names(vss), value = TRUE))
  stopifnot(length(tstsetups) > 0)
  # Note: As mentioned in issue #149, the reference level of a categorical
  # predictor actually should not have its own coefficient:
  len_penal <- sum(grepl("^xco\\.", names(dat))) +
    sum(sapply(dat[, grep("^xca\\.", names(dat))], nlevels))
  penal_crr <- rep(1, len_penal)
  stopifnot(len_penal >= 5)
  idx_penal_0 <- c(1, 2) # A few variables without cost.
  idx_penal_Inf <- c(3) # One variable with infinite penalty.
  penal_crr[idx_penal_0] <- 0
  penal_crr[idx_penal_Inf] <- Inf
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    vs_penal <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
           penalty = penal_crr),
      args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
    ))
    vsel_tester(
      vs_penal,
      refmod_expected = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]],
      solterms_len_expected = args_vs_i$nterms_max,
      method_expected = "L1",
      nclusters_expected = args_vs_i$nclusters,
      nclusters_pred_expected = args_vs_i$nclusters_pred,
      info_str = tstsetup
    )
    # Check that the variables with no cost are selected first and the ones
    # with infinite penalty last:
    formula_crr <- refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]]$formula
    solterms_orig <- setdiff(split_formula(formula_crr), "1")
    solterms_penal <- vs_penal$solution_terms
    # Note: This test probably needs to be adopted properly to categorical
    # predictors.
    stopifnot(all(grep("^xca\\.", solterms_orig) >= max(c(idx_penal_0,
                                                          idx_penal_Inf))))
    expect_identical(solterms_penal[seq_along(idx_penal_0)],
                     solterms_orig[idx_penal_0],
                     info = tstsetup)
    expect_identical(rev(solterms_penal)[seq_along(idx_penal_Inf)],
                     rev(solterms_orig[idx_penal_Inf]),
                     info = tstsetup)
  }
})

# cv_varsel() -------------------------------------------------------------

context("cv_varsel()")

test_that(paste(
  "`object` of class \"refmodel\", correctly specified `method`, `nterms_max`,",
  "`nclusters`, and `nclusters_pred` lead to correct output structure"
), {
  skip_if_not(exists("cvvss"))
  for (tstsetup in names(cvvss)) {
    mod_crr <- args_cvvs[[tstsetup]]$mod
    fam_crr <- args_cvvs[[tstsetup]]$fam
    meth_exp_crr <- args_cvvs[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    vsel_tester(
      cvvss[[tstsetup]],
      with_cv = TRUE,
      refmod_expected = refmods[[mod_crr]][[fam_crr]],
      solterms_len_expected = args_cvvs[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "LOO",
      valsearch_expected = args_cvvs[[tstsetup]]$validate_search,
      nclusters_expected = args_cvvs[[tstsetup]]$nclusters,
      nclusters_pred_expected = args_cvvs[[tstsetup]]$nclusters_pred,
      info_str = tstsetup
    )
  }
})

test_that("specifying `object` incorrectly leads to an error", {
  expect_error(cv_varsel(rnorm(5)),
               "^no applicable method for")
})

test_that("specifying `method` incorrectly leads to an error", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      expect_error(cv_varsel(refmods[[!!mod_nm]][[!!fam_nm]],
                             method = "k-fold"),
                   "^Unknown search method$")
      if (mod_nm == "glmm") {
        ### Excluded because of issue #171:
        # expect_error(cv_varsel(refmods[[!!mod_nm]][[!!fam_nm]], method = "L1"),
        #              "^L1 search is not supported for multilevel models\\.$")
        ###
      } else if (!mod_nm %in% c("glm", "glmm")) {
        ### TODO:
        stop("Still to-do.")
        # expect_error(cv_varsel(refmods[[!!mod_nm]][[!!fam_nm]], method = "L1"),
        #              "ENTER EXPECTED TEXT HERE")
        ###
      }
    }
  }
})

test_that("specifying `cv_method` incorrectly leads to an error", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      expect_error(cv_varsel(refmods[[!!mod_nm]][[!!fam_nm]],
                             cv_method = "k-fold"),
                   "^Unknown cross-validation method$")
    }
  }
})

test_that("specifying `nloo` incorrectly leads to an error", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      # Use SW() because of occasional warnings concerning Pareto k diagnostics:
      expect_error(SW(cv_varsel(refmods[[!!mod_nm]][[!!fam_nm]],
                                nloo = -1)),
                   "^nloo must be at least 1$",
                   info = paste(mod_nm, fam_nm, sep = "__"))
    }
  }
})

test_that(paste(
  "setting `nloo` at least as large as the number of observations doesn't",
  "change results"
), {
  skip_if_not(exists("cvvss"))
  nloo_tst <- n_tst + 1L
  # To save time: Pick only a single scenario with a GLM and a forward search
  # (the latter because of `validate_search = FALSE`):
  tstsetups <- grep("^glm\\..*\\.forward", names(cvvss), value = TRUE)[1]
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    cvvs_nloo <- do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$mod_nm]][[args_cvvs_i$fam_nm]],
           nloo = nloo_tst),
      args_cvvs_i[setdiff(names(args_cvvs_i), c("tstsetup"))]
    ))
    expect_equal(cvvs_nloo, cvvss[[tstsetup]], info = tstsetup)
  }
})

test_that(paste(
  "setting `nloo` smaller than the number of observations leads to correct",
  "output structure"
), {
  skip_if_not(exists("cvvss"))
  nloo_tst <- n_tst - 1L
  # To save time: Pick only a single scenario with a GLM and a forward search
  # (the latter because of `validate_search = FALSE`):
  tstsetups <- grep("^glm\\..*\\.forward", names(cvvss), value = TRUE)[1]
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    mod_crr <- args_cvvs_i$mod
    fam_crr <- args_cvvs_i$fam
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    # Expect a warning (see issue #172):
    expect_warning(
      cvvs_nloo <- do.call(cv_varsel, c(
        list(object = refmods[[args_cvvs_i$mod_nm]][[args_cvvs_i$fam_nm]],
             nloo = nloo_tst),
        args_cvvs_i[setdiff(names(args_cvvs_i), c("tstsetup"))]
      )),
      "longer object length is not a multiple of shorter object length",
      info = tstsetup
    )
    vsel_tester(
      cvvs_nloo,
      with_cv = TRUE,
      refmod_expected = refmods[[mod_crr]][[fam_crr]],
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "LOO",
      valsearch_expected = args_cvvs_i$validate_search,
      nclusters_expected = args_cvvs_i$nclusters,
      nclusters_pred_expected = args_cvvs_i$nclusters_pred,
      expect_const_obs_w = FALSE,
      info_str = tstsetup
    )
    expect_false(isTRUE(all.equal(cvvs_nloo, cvvss[[tstsetup]])),
                 info = tstsetup)
    # TODO (extend):
    for (j in seq_along(cvvs_nloo$summaries$sub)) {
      expect_equal(sum(!is.na(cvvs_nloo$summaries$sub[[!!j]]$lppd)), nloo_tst,
                   info = tstsetup)
    }
  }
})

test_that("`validate_search` works", {
  skip_if_not(exists("cvvss"))
  tstsetups <- grep("^glm\\..*\\.default_meth", names(cvvss), value = TRUE)
  stopifnot(length(tstsetups) > 0)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    stopifnot(is.null(args_cvvs_i$validate_search) ||
                isTRUE(args_cvvs_i$validate_search))
    mod_crr <- args_cvvs_i$mod
    fam_crr <- args_cvvs_i$fam
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    # Use SW() because of occasional warnings concerning Pareto k diagnostics:
    SW(cvvs_valsearch <- do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$mod_nm]][[args_cvvs_i$fam_nm]],
           validate_search = FALSE),
      args_cvvs_i[setdiff(names(args_cvvs_i), c("tstsetup"))]
    )))
    vsel_tester(
      cvvs_valsearch,
      with_cv = TRUE,
      refmod_expected = refmods[[mod_crr]][[fam_crr]],
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "LOO",
      valsearch_expected = FALSE,
      nclusters_expected = args_cvvs_i$nclusters,
      nclusters_pred_expected = args_cvvs_i$nclusters_pred,
      info_str = tstsetup
    )
    expect_false(isTRUE(all.equal(cvvs_valsearch, cvvss[[tstsetup]])),
                 info = tstsetup)
    # TODO (extend):
    expect_true(all(summary(cvvs_valsearch)$selection$elpd.loo >=
                      summary(cvvss[[tstsetup]])$selection$elpd.loo),
                info = tstsetup)
  }
})

# TODO:
test_that(paste(
  "object returned by cv_varsel, kfold contains the relevant",
  "fields"
), {
  for (i in seq_len(length(cv_kf_list))) {
    i_inf <- names(cv_kf_list)[i]
    for (j in seq_len(length(cv_kf_list[[i]]))) {
      j_inf <- names(cv_kf_list[[i]])[j]
      # solution_terms seems legit
      expect_length(cv_kf_list[[i]][[j]]$solution_terms, nterms)
      expect_true(all(!is.na(match(
        colnames(fit_gauss$data[, -1]),
        cv_kf_list[[i]][[j]]$solution_terms
      ))),
      info = paste(i_inf, j_inf)
      )
      # kl seems legit
      expect_length(cv_kf_list[[i]][[j]]$kl, nterms + 1)

      # decreasing
      expect_equal(cv_kf_list[[i]][[j]]$kl[-1],
                   cummin(cv_kf_list[[i]][[j]]$kl[-1]),
                   info = paste(i_inf, j_inf),
                   tolerance = 24e-2
      )

      # summaries seems legit
      expect_named(cv_kf_list[[i]][[j]]$summaries, c("sub", "ref"),
                   info = paste(i_inf, j_inf)
      )
      expect_length(cv_kf_list[[i]][[j]]$summaries$sub, nterms + 1)
      expect_named(cv_kf_list[[i]][[j]]$summaries$sub[[1]],
                   c("mu", "lppd", "w"),
                   ignore.order = TRUE, info = paste(i_inf, j_inf)
      )
      expect_named(cv_kf_list[[i]][[j]]$summaries$ref, c("mu", "lppd"),
                   ignore.order = TRUE, info = paste(i_inf, j_inf)
      )
      # family seems legit
      expect_equal(cv_kf_list[[i]][[j]]$family$family,
                   cv_kf_list[[i]][[j]]$family$family,
                   info = paste(i_inf, j_inf)
      )
      expect_equal(cv_kf_list[[i]][[j]]$family$link,
                   cv_kf_list[[i]][[j]]$family$link,
                   info = paste(i_inf, j_inf)
      )
      expect_true(length(cv_kf_list[[i]][[j]]$family) >=
                    length(cv_kf_list[[i]][[j]]$family$family),
                  info = paste(i_inf, j_inf)
      )
      # pct_solution_terms_cv seems legit
      expect_equal(dim(cv_kf_list[[i]][[j]]$pct_solution_terms_cv),
                   c(nterms, nterms + 1),
                   info = paste(i_inf, j_inf)
      )
      expect_true(all(
        cv_kf_list[[i]][[j]]$pct_solution_terms_cv[, -1] <= 1 &
          cv_kf_list[[i]][[j]]$pct_solution_terms_cv[, -1] >= 0
      ),
      info = paste(i_inf, j_inf)
      )
      expect_equal(cv_kf_list[[i]][[j]]$pct_solution_terms_cv[, 1],
                   1:nterms,
                   info = paste(i_inf, j_inf)
      )
      expect_equal(colnames(cv_kf_list[[i]][[j]]$pct_solution_terms_cv),
                   c("size", cv_kf_list[[i]][[j]]$solution_terms),
                   info = paste(i_inf, j_inf)
      )
    }
  }
})

test_that("K is valid for cv_method='kfold'", {
  # the chains, seed and iter arguments to the rstanarm functions here must
  # be specified directly rather than through a variable (eg, seed = 1235
  # instead of seed = seed), otherwise when the calls are evaluated in
  # refmodel$cvfun() they may not be found in the evaluation frame of the
  # calling function, causing the test to fail
  glm_simp <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
                       family = poisson(), data = df_poiss,
                       chains = 2, seed = 1235, iter = 400)

  expect_error(
    cv_varsel(glm_simp, cv_method = "kfold", K = 1),
    "must be at least 2"
  )
  expect_error(
    cv_varsel(glm_simp, cv_method = "kfold", K = 1000),
    "cannot exceed n"
  )
  expect_error(
    cv_varsel(glm_simp, cv_method = "kfold", K = c(4, 9)),
    "a single integer value"
  )
  expect_error(
    cv_varsel(glm_simp, cv_method = "kfold", K = "a"),
    "a single integer value"
  )
  expect_error(
    cv_varsel(glm_simp, cv_method = "kfold", K = df_poiss),
    "a single integer value"
  )
})

test_that("providing `cvfits` works", {
  out <- SW({
    k_fold <- kfold(glm_simp, K = 2, save_fits = TRUE)
    folds <- seq_len(nrow(glm_simp$data))
    for (K in seq_len(2)) {
      folds[as.numeric(rownames(k_fold$fit[[K]]$data))] <- K
    }
    attr(k_fold, "folds") <- folds
    fit_cv <- cv_varsel(glm_simp,
                        cv_method = "kfold", cvfits = k_fold,
                        ndraws = ndraws, ndraws_pred = ndraws_pred,
                        verbose = FALSE)
  })
  expect_false(any(grepl("k_fold not provided", out)))
  expect_length(fit_cv$solution_terms, nterms)

  # kl seems legit
  expect_length(fit_cv$kl, nterms + 1)

  # decreasing
  expect_equal(fit_cv$kl, cummin(fit_cv$kl), tolerance = 1e-3)

  # summaries seems legit
  expect_named(fit_cv$summaries, c("sub", "ref"))
  expect_length(fit_cv$summaries$sub, nterms + 1)
  expect_named(fit_cv$summaries$sub[[1]], c("mu", "lppd", "w"),
               ignore.order = TRUE
  )
  expect_named(fit_cv$summaries$ref, c("mu", "lppd"),
               ignore.order = TRUE
  )
  # family seems legit
  expect_equal(
    fit_cv$family$family,
    fit_cv$family$family
  )
  expect_equal(fit_cv$family$link, fit_cv$family$link)
  expect_true(length(fit_cv$family) >= length(fit_cv$family$family))
  # pct_solution_terms_cv seems legit
  expect_equal(dim(fit_cv$pct_solution_terms_cv), c(nterms, nterms + 1))
  expect_true(all(fit_cv$pct_solution_terms_cv[, -1] <= 1 &
                    fit_cv$pct_solution_terms_cv[, -1] >= 0))

  expect_equal(fit_cv$pct_solution_terms_cv[, 1], 1:nterms)
  expect_equal(
    colnames(fit_cv$pct_solution_terms_cv),
    c("size", fit_cv$solution_terms)
  )
})

# summary() ---------------------------------------------------------------

context("summary()")

valid_stats_all <- c("elpd", "mlpd")
valid_stats_gauss_only <- c("mse", "rmse")
valid_stats_binom_only <- c("acc", "auc")
valid_stats_gauss <- c(valid_stats_all, valid_stats_gauss_only)
valid_stats_binom <- c(valid_stats_all, valid_stats_binom_only)
vs_funs <- c(summary, plot, suggest_size)

## test_that("invalid objects are rejected", {
##   for (fun in vs_funs) {
##     expect_error(fun(NULL), "is not a variable selection object")
##     expect_error(fun(fit_gauss), "is not a variable selection object")
##   }
## })

## test_that("invalid stats are rejected", {
##   for (fun in vs_funs) {
##     expect_error(fun(vs_list[[1]][["gauss"]], stat = NULL),
##                  "specified as NULL")
##     expect_error(fun(vs_list[[1]][["gauss"]], stat = NA),
##                  "not recognized")
##     expect_error(fun(vs_list[[1]][["gauss"]], stat = "zzz"),
##                  "not recognized")
##     expect_error(fun(vs_list[[1]][["gauss"]], stat = "acc"),
##                  "available only for the binomial family")
##     expect_error(
##       fun(vs_list[[1]][["gauss"]], stat = "auc"),
##       "available only for the binomial family"
##     )
##   }
## })

## test_that("invalid 'baseline' arguments are rejected", {
##   expect_error(
##     summary(vs_list[[1]][["gauss"]], baseline = "zzz"),
##     "Argument 'baseline' must be either 'ref' or 'best'"
##   )
## })

test_that("summary output seems legit", {
  skip_on_cran()
  for (i in seq_along(cvs_list)) {
    for (j in seq_along(cvs_list[[i]])) {
      cvs <- cvs_list[[i]][[j]]
      if (cvs$family$family == "gaussian") {
        stats_str <- valid_stats_gauss
      } else if (cvs$family$family == "binomial") {
        stats_str <- valid_stats_binom
      } else {
        stats_str <- valid_stats_all
      }
      cv_method <- cvs_list[[i]][[j]]$cv_method
      stats <- summary(cvs,
                       stats = stats_str,
                       type = c("mean", "lower", "upper", "se")
      )$selection
      expect_true(nrow(stats) == nterms + 1)
      expect_true(all(c(
        "size", "solution_terms", paste0(stats_str, ".", tolower(cv_method)),
        paste0(stats_str, ".", c("se", "upper", "lower"))
      ) %in% names(stats)))
      expect_true(all(stats[, paste0("mlpd.", tolower(cv_method))] >
                        stats[, "mlpd.lower"]))
      expect_true(all(stats[, paste0("mlpd.", tolower(cv_method))] <
                        stats[, "mlpd.upper"]))
    }
  }
})

test_that("summary works with reference models", {
  for (i in seq_along(vsref_list)) {
    for (j in seq_along(vsref_list[[i]])) {
      vs <- vsref_list[[i]][[j]]
      if (vs$family$family == "gaussian") {
        stats_str <- valid_stats_gauss
      } else {
        stats_str <- valid_stats_binom
      }
      stats <- summary(vs, stats = stats_str)$selection
      expect_true(is.data.frame(stats))
    }
  }
})

# print() -----------------------------------------------------------------

context("print()")

test_that("print() works as expected", {

  skip_on_cran()
  # default rounding
  expect_output(out <- print(vs_list[[1]][[1]]))
  expect_equal(out$selection$elpd, round(out$selection$elpd, 2),
               tolerance = 1e-3
  )
  expect_output(out <- print(cvs_list[[1]][[1]]))
  expect_equal(out$selection$elpd, round(out$selection$elpd, 2),
               tolerance = 1e-3
  )

  # rounding to 4 decimal places
  expect_output(out <- print(vs_list[[1]][[1]], digits = 4))
  expect_equal(out$selection$elpd, round(out$selection$elpd, 4),
               tolerance = 1e-3
  )
  expect_output(out <- print(cvs_list[[1]][[1]], digits = 4))
  expect_equal(out$selection$elpd, round(out$selection$elpd, 4),
               tolerance = 1e-3
  )
  # options to summary
  expect_output(out <- print(vs_list[[1]][[1]],
                             nterms_max = 3,
                             stats = "mse"
  ))
  expect_equal(nrow(out$selection) - 1, 3)
  expect_named(out$selection, c(
    "size", "solution_terms",
    "mse", "se",
    "diff", "diff.se"
  ))

  expect_output(out <- print(cvs_list[[1]][[1]],
                             nterms_max = 3,
                             stats = "mse"
  ))
  expect_equal(nrow(out$selection) - 1, 3)
  expect_named(out$selection, c(
    "size", "solution_terms",
    paste0("mse.", tolower(out$cv_method)), "se",
    "diff", "diff.se"
    # "pct_solution_terms_cv"
  ))
})

# plot() ------------------------------------------------------------------

context("plot()")

test_that("plotting works", {
  expect_s3_class(plot(vs_list[[1]][[1]]), "ggplot")
  expect_visible(plot(vs_list[[1]][[1]], nterms_max = 3))
})

test_that("invalid 'baseline' arguments are rejected", {
  expect_error(
    plot(vs_list[[1]][[1]], baseline = "zzz"),
    "Argument 'baseline' must be either 'ref' or 'best'"
  )
})

test_that("the value of nterms_max is valid", {
  expect_error(
    plot(vs_list[[1]][[1]], nterms_max = 0),
    "nterms_max must be at least 1"
  )
})

## test_that("nterms_max is capped to the largest model size", {
##   expect_equal(
##     plot(vs_list[[1]][[1]]),
##     plot(vs_list[[1]][[1]], nterms_max = 1000)
##   )
## })

# suggest_size() ----------------------------------------------------------

context("suggest_size()")

test_that("suggest_size() checks the length of stat", {
  expect_error(suggest_size(vs_list[[1]][["gauss"]], stat = valid_stats_all),
               "Only one statistic")
})

test_that("suggest_size() works on all stats", {
  for (stat in valid_stats_gauss) {
    suggested_size <- suggest_size(vs_list[[1]][["gauss"]], stat = stat)
    expect_true(!is.na(suggested_size))
    expect_true(suggested_size >= 0)
  }
  for (stat in valid_stats_binom) {
    suggested_size <- suggest_size(vs_list[[1]][["binom"]], stat = stat)
    expect_true(!is.na(suggested_size))
    expect_true(suggested_size >= 0)
  }
})
