# varsel() ----------------------------------------------------------------

context("varsel()")

test_that(paste(
  "`object` of class \"refmodel\", `method`, `nterms_max`, `nclusters`, and",
  "`nclusters_pred` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(vss)) {
    tstsetup_ref <- args_vs[[tstsetup]]$tstsetup_ref
    mod_crr <- args_vs[[tstsetup]]$mod_nm
    fam_crr <- args_vs[[tstsetup]]$fam_nm
    meth_exp_crr <- args_vs[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    vsel_tester(
      vss[[tstsetup]],
      refmod_expected = refmods[[tstsetup_ref]],
      solterms_len_expected = args_vs[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      nclusters_expected = args_vs[[tstsetup]]$nclusters,
      nclusters_pred_expected = args_vs[[tstsetup]]$nclusters_pred,
      info_str = tstsetup
    )
  }
})

test_that("invalid `object` fails", {
  expect_error(varsel(rnorm(5), verbose = FALSE),
               "no applicable method")
})

test_that("invalid `method` fails", {
  for (tstsetup in names(refmods)) {
    expect_error(varsel(refmods[[tstsetup]], method = "k-fold"),
                 "Unknown search method",
                 info = tstsetup)
    if (args_ref[[tstsetup]]$mod_nm != "glm") {
      expect_error(varsel(refmods[[tstsetup]], method = "L1"),
                   "^L1 search is only supported for GLMs\\.$",
                   info = tstsetup)
    }
  }
})

test_that("`seed` works (and restores the RNG state afterwards)", {
  # Note: Extensive tests for reproducibility may be found among the tests for
  # .get_refdist().
  skip_if_not(run_vs)
  # To save time:
  tstsetups <- grep("^glm\\.gauss\\.", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    vs_orig <- vss[[tstsetup]]
    rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_new1 <- .Random.seed
    vs_new <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$tstsetup_ref]],
           seed = args_vs_i$seed + 1L),
      excl_nonargs(args_vs_i, nms_excl_add = "seed")
    ))
    .Random.seed_new2 <- .Random.seed
    rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_repr1 <- .Random.seed
    vs_repr <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$tstsetup_ref]]),
      excl_nonargs(args_vs_i)
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

test_that("`d_test` works", {
  skip_if_not(run_vs)
  tstsetups <- names(vss)
  ### Alternative with less test setups:
  # tstsetups <- grep("^glm\\.", names(vss), value = TRUE)
  ###
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    tstsetup_ref <- args_vs_i$tstsetup_ref
    mod_crr <- args_vs_i$mod_nm
    fam_crr <- args_vs_i$fam_nm
    refmod_crr <- refmods[[tstsetup_ref]]
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
    warn_expected <- if (mod_crr == "glm") {
      paste("^'offset' argument is NULL but it looks like you estimated the",
            "model using an offset term\\.$")
    } else {
      NA
    }
    expect_warning(
      vs_repr <- do.call(varsel, c(
        list(object = refmod_crr, d_test = d_test_crr),
        excl_nonargs(args_vs_i)
      )),
      warn_expected,
      info = tstsetup
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
  skip_if_not(run_vs)
  regul_tst <- 1e-1
  for (mod_crr in setdiff(mod_nms, "glm")) {
    tstsetups <- head(grep(paste0("^", mod_crr, "\\.gauss"), names(vss),
                           value = TRUE),
                      1)
    for (tstsetup in tstsetups) {
      args_vs_i <- args_vs[[tstsetup]]
      vs_regul <- do.call(varsel, c(
        list(object = refmods[[args_vs_i$tstsetup_ref]],
             regul = regul_tst),
        excl_nonargs(args_vs_i)
      ))
      expect_equal(vs_regul, vss[[tstsetup]], info = tstsetup)
    }
  }
})

test_that(paste(
  "for GLMs with L1 search, `regul` only has an effect on prediction, not on",
  "selection"
), {
  skip_if_not(run_vs)
  regul_tst <- c(regul_default, 1e-1, 1e2)
  stopifnot(regul_tst[1] == regul_default)
  stopifnot(all(diff(regul_tst) > 0))
  tstsetups <- setdiff(grep("^glm\\.", names(vss), value = TRUE),
                       grep("^glm\\..*\\.forward", names(vss), value = TRUE))
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    m_max <- args_vs_i$nterms_max + 1L
    ssq_regul_prd <- array(dim = c(length(regul_tst), m_max))
    for (j in seq_along(regul_tst)) {
      if (regul_tst[j] == regul_default) {
        vs_regul <- vss[[tstsetup]]
      } else {
        vs_regul <- do.call(varsel, c(
          list(object = refmods[[args_vs_i$tstsetup_ref]],
               regul = regul_tst[j]),
          excl_nonargs(args_vs_i)
        ))
        vsel_tester(
          vs_regul,
          refmod_expected = refmods[[args_vs_i$tstsetup_ref]],
          solterms_len_expected = args_vs_i$nterms_max,
          method_expected = "L1",
          nclusters_expected = args_vs_i$nclusters,
          nclusters_pred_expected = args_vs_i$nclusters_pred,
          info_str = tstsetup
        )
        # Expect equality for all components not related to prediction:
        expect_equal(vs_regul[setdiff(vsel_nms, vsel_nms_pred)],
                     vss[[tstsetup]][setdiff(vsel_nms, vsel_nms_pred)],
                     info = paste(tstsetup, j, sep = "__"))
        # Expect inequality for the components related to prediction (but note
        # that the components from `vsel_nms_pred_opt` can be, but don't need to
        # be differing):
        for (vsel_nm in setdiff(vsel_nms_pred, vsel_nms_pred_opt)) {
          expect_false(isTRUE(all.equal(vs_regul[[vsel_nm]],
                                        vss[[tstsetup]][[vsel_nm]])),
                       info = paste(tstsetup, j, vsel_nm, sep = "__"))
        }
      }
      # Check the inequality of the prediction components in detail: Expect a
      # reduction of the sum of the squared coefficients (excluding the
      # intercept) for increasing `regul`:
      for (m in seq_len(m_max)) {
        # Since varsel() doesn't output object `p_sub`, use the linear predictor
        # here (instead of the coefficients themselves, which would only be
        # accessible from `p_sub`):
        mu_jm_regul <- vs_regul$family$linkfun(vs_regul$summaries$sub[[m]]$mu) -
          offs_tst
        # In fact, `sum((mu - offset - intercept)^2)` would make more sense than
        # `var(mu - offset) = sum((mu - offset - mean(mu - offset))^2)` but
        # since varsel() doesn't output object `p_sub`, the intercept from the
        # prediction is not accessible here.
        ssq_regul_prd[j, m] <- var(mu_jm_regul)
      }
    }
    # For the intercept-only model, the linear predictor consists only
    # of the intercept, so we expect no variation in `mu_jm_regul`:
    expect_true(all(ssq_regul_prd[, 1] <= 1e-5), info = tstsetup)
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
  skip_if_not(run_vs)
  regul_tst <- c(regul_default, 1e-1, 1e2)
  stopifnot(regul_tst[1] == regul_default)
  stopifnot(all(diff(regul_tst) > 0))
  tstsetups <- grep("^glm\\..*\\.forward", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    m_max <- args_vs_i$nterms_max + 1L
    if (identical(args_vs_i$method, "forward")) {
      ncl_crr <- args_vs_i$nclusters
    } else {
      ncl_crr <- 1L
    }
    if (!grepl("\\.spclformul", tstsetup)) {
      tol_alpha <- 3e-1
    } else {
      tol_alpha <- 5e-1
    }
    ssq_regul_sel_alpha <- array(dim = c(length(regul_tst), m_max, ncl_crr))
    ssq_regul_sel_beta <- array(dim = c(length(regul_tst), m_max, ncl_crr))
    ssq_regul_prd <- array(dim = c(length(regul_tst), m_max))
    for (j in seq_along(regul_tst)) {
      if (regul_tst[j] == regul_default) {
        vs_regul <- vss[[tstsetup]]
      } else {
        vs_regul <- do.call(varsel, c(
          list(object = refmods[[args_vs_i$tstsetup_ref]],
               regul = regul_tst[j]),
          excl_nonargs(args_vs_i)
        ))
        vsel_tester(
          vs_regul,
          refmod_expected = refmods[[args_vs_i$tstsetup_ref]],
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
        # here (instead of the coefficients themselves, which would only be
        # accessible from `p_sub`):
        mu_jm_regul <- vs_regul$family$linkfun(vs_regul$summaries$sub[[m]]$mu) -
          offs_tst
        # In fact, `sum((mu - offset - intercept)^2)` would make more sense than
        # `var(mu - offset) = sum((mu - offset - mean(mu - offset))^2)` but
        # since varsel() doesn't output object `p_sub`, the intercept from the
        # prediction is not accessible here.
        ssq_regul_prd[j, m] <- var(mu_jm_regul)
      }
    }
    # Selection:
    # For the intercept-only model:
    for (nn in seq_len(dim(ssq_regul_sel_alpha)[3])) {
      expect_length(unique(ssq_regul_sel_alpha[, 1, !!nn]), 1)
    }
    # All other (i.e., not intercept-only) models:
    for (j in seq_len(dim(ssq_regul_sel_alpha)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_sel_alpha)[2])[-1]) {
        for (nn in seq_len(dim(ssq_regul_sel_alpha)[3])) {
          expect_equal(ssq_regul_sel_alpha[!!j, !!m, !!nn],
                       ssq_regul_sel_alpha[j - 1, m, nn],
                       tolerance = tol_alpha)
        }
      }
    }
    # For the intercept-only model:
    expect_true(all(is.na(ssq_regul_sel_beta[, 1, ])), info = tstsetup)
    # All other (i.e., not intercept-only) models:
    ssq_regul_sel_beta_cond <- array(
      dim = dim(ssq_regul_sel_beta) + c(-1L, -1L, 0L)
    )
    for (j in seq_len(dim(ssq_regul_sel_beta)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_sel_beta)[2])[-1]) {
        for (nn in seq_len(dim(ssq_regul_sel_beta)[3])) {
          ssq_regul_sel_beta_cond[j - 1, m - 1, nn] <-
            ssq_regul_sel_beta[j, m, nn] < ssq_regul_sel_beta[j - 1, m, nn]
        }
      }
    }
    sum_as_unexpected <- 0L
    expect_true(sum(!ssq_regul_sel_beta_cond) <= sum_as_unexpected,
                info = tstsetup)
    # Prediction:
    # For the intercept-only model, the linear predictor consists only
    # of the intercept, so we expect no variation in `mu_jm_regul`:
    expect_true(all(ssq_regul_prd[, 1] <= 1e-12), info = tstsetup)
    # All other (i.e., not intercept-only) models:
    for (j in seq_len(dim(ssq_regul_prd)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_prd)[2])[-1]) {
        expect_lt(ssq_regul_prd[!!j, !!m], ssq_regul_prd[j - 1, m])
      }
    }
  }
})

## Penalty ----------------------------------------------------------------

test_that("`penalty` of invalid length fails", {
  skip_if_not(run_vs)
  tstsetups <- setdiff(
    grep("^glm\\.", names(args_vs), value = TRUE),
    grep("^glm\\..*\\.forward", names(args_vs), value = TRUE)
  )
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    penal_possbl <- get_penal_possbl(fits[[args_vs_i$tstsetup_fit]]$formula)
    len_penal <- length(penal_possbl)
    # The `penalty` objects to be tested:
    penal_tst <- list(rep(1, len_penal + 1), rep(1, len_penal - 1))
    for (penal_crr in penal_tst) {
      expect_error(
        do.call(varsel, c(
          list(object = refmods[[args_vs_i$tstsetup_ref]],
               penalty = penal_crr),
          excl_nonargs(args_vs_i)
        )),
        paste0("^Incorrect length of penalty vector \\(should be ",
               len_penal, "\\)\\.$"),
        info = paste(tstsetup, which(sapply(penal_tst, identical, penal_crr)),
                     sep = "__")
      )
    }
  }
})

test_that("for forward search, `penalty` has no effect", {
  skip_if_not(run_vs)
  penal_tst <- 2
  tstsetups <- union(grep("\\.forward", names(vss), value = TRUE),
                     grep("^glm\\.", names(vss), value = TRUE, invert = TRUE))
  # To save time:
  tstsetups <- head(tstsetups, 1)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    vs_penal <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$tstsetup_ref]],
           penalty = penal_tst),
      excl_nonargs(args_vs_i)
    ))
    expect_equal(vs_penal, vss[[tstsetup]], info = tstsetup)
  }
})

test_that("for L1 search, `penalty` has an expected effect", {
  skip_if_not(run_vs)
  tstsetups <- setdiff(grep("^glm\\.", names(vss), value = TRUE),
                       grep("^glm\\..*\\.forward", names(vss), value = TRUE))
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]

    penal_possbl <- get_penal_possbl(fits[[args_vs_i$tstsetup_fit]]$formula)
    len_penal <- length(penal_possbl)
    penal_crr <- rep(1, len_penal)
    stopifnot(len_penal >= 3)
    idx_penal_0 <- c(1, 2) # A few variables without cost.
    idx_penal_Inf <- c(3) # One variable with infinite penalty.
    penal_crr[idx_penal_0] <- 0
    penal_crr[idx_penal_Inf] <- Inf
    # Note: This test probably needs to be adopted properly to categorical
    # predictors. As a workaround, perform the following check:
    stopifnot(all(grep("^xca\\.", penal_possbl) >= max(c(idx_penal_0,
                                                         idx_penal_Inf))))

    vs_penal <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$tstsetup_ref]],
           penalty = penal_crr),
      excl_nonargs(args_vs_i, nms_excl_add = "nterms_max")
    ))
    nterms_max_crr <- count_terms_in_formula(
      fits[[args_vs_i$tstsetup_fit]]$formula
    ) - 1L
    vsel_tester(
      vs_penal,
      refmod_expected = refmods[[args_vs_i$tstsetup_ref]],
      solterms_len_expected = nterms_max_crr,
      method_expected = "L1",
      nclusters_expected = args_vs_i$nclusters,
      nclusters_pred_expected = args_vs_i$nclusters_pred,
      info_str = tstsetup
    )
    # Check that the variables with no cost are selected first and the ones
    # with infinite penalty last:
    solterms_penal <- vs_penal$solution_terms
    expect_identical(solterms_penal[seq_along(idx_penal_0)],
                     penal_possbl[idx_penal_0],
                     info = tstsetup)
    expect_identical(rev(solterms_penal)[seq_along(idx_penal_Inf)],
                     rev(penal_possbl[idx_penal_Inf]),
                     info = tstsetup)
  }
})

# cv_varsel() -------------------------------------------------------------

context("cv_varsel()")

test_that(paste(
  "`object` of class \"refmodel\", `method`, `cv_method`, `nterms_max`,",
  "`nclusters`, and `nclusters_pred` work"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(cvvss)) {
    mod_crr <- args_cvvs[[tstsetup]]$mod_nm
    fam_crr <- args_cvvs[[tstsetup]]$fam_nm
    meth_exp_crr <- args_cvvs[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    vsel_tester(
      cvvss[[tstsetup]],
      with_cv = TRUE,
      refmod_expected = refmods[[args_cvvs[[tstsetup]]$tstsetup_ref]],
      solterms_len_expected = args_cvvs[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = args_cvvs[[tstsetup]]$cv_method,
      valsearch_expected = args_cvvs[[tstsetup]]$validate_search,
      nclusters_expected = args_cvvs[[tstsetup]]$nclusters,
      nclusters_pred_expected = args_cvvs[[tstsetup]]$nclusters_pred,
      info_str = tstsetup
    )
  }
})

test_that("invalid `object` fails", {
  expect_error(cv_varsel(rnorm(5)),
               "^no applicable method for")
})

test_that("invalid `method` fails", {
  for (tstsetup in names(refmods)) {
    expect_error(cv_varsel(refmods[[tstsetup]], method = "k-fold"),
                 "^Unknown search method$",
                 info = tstsetup)
    if (args_ref[[tstsetup]]$mod_nm != "glm") {
      expect_error(cv_varsel(refmods[[tstsetup]], method = "L1"),
                   "^L1 search is only supported for GLMs\\.$",
                   info = tstsetup)
    }
  }
})

test_that("invalid `cv_method` fails", {
  for (tstsetup in names(refmods)) {
    expect_error(cv_varsel(refmods[[tstsetup]], cv_method = "k-fold"),
                 "^Unknown cross-validation method$",
                 info = tstsetup)
  }
})

test_that("`seed` works (and restores the RNG state afterwards)", {
  # Note: Extensive tests for reproducibility may be found among the tests for
  # .get_refdist().
  skip_if_not(run_cvvs)
  # To save time:
  tstsetups <- grep("^glm\\.gauss", names(cvvss), value = TRUE)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    cvvs_orig <- cvvss[[tstsetup]]
    rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_new1 <- .Random.seed
    # Use SW() because of occasional warnings concerning Pareto k diagnostics:
    SW(cvvs_new <- do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]],
           seed = args_cvvs_i$seed + 1L),
      excl_nonargs(args_cvvs_i, nms_excl_add = "seed")
    )))
    .Random.seed_new2 <- .Random.seed
    rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_repr1 <- .Random.seed
    SW(cvvs_repr <- do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]]),
      excl_nonargs(args_cvvs_i)
    )))
    .Random.seed_repr2 <- .Random.seed
    # Expected equality:
    expect_equal(cvvs_repr, cvvs_orig, info = tstsetup)
    if (!identical(args_cvvs_i$cv_method, "kfold")) {
      expect_equal(.Random.seed_new2, .Random.seed_new1, info = tstsetup)
      expect_equal(.Random.seed_repr2, .Random.seed_repr1, info = tstsetup)
    }
    # Expected inequality:
    expect_false(isTRUE(all.equal(cvvs_new, cvvs_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(rand_new, rand_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_repr2, .Random.seed_new2)),
                 info = tstsetup)
  }
})

test_that("invalid `nloo` fails", {
  for (tstsetup in names(refmods)) {
    # Use SW() because of occasional warnings concerning Pareto k diagnostics:
    expect_error(SW(cv_varsel(refmods[[tstsetup]], nloo = -1)),
                 "^nloo must be at least 1$",
                 info = tstsetup)
  }
})

test_that(paste(
  "setting `nloo` at least as large as the number of observations doesn't",
  "change results"
), {
  skip_if_not(run_cvvs)
  nloo_tst <- nobsv + 1L
  tstsetups <- grep("^glm\\.gauss\\..*\\.default_cvmeth", names(cvvss),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    # Use SW() because of occasional warnings concerning Pareto k diagnostics:
    SW(cvvs_nloo <- do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]],
           nloo = nloo_tst),
      excl_nonargs(args_cvvs_i)
    )))
    expect_equal(cvvs_nloo, cvvss[[tstsetup]], info = tstsetup)
  }
})

test_that("setting `nloo` smaller than the number of observations works", {
  skip_if_not(run_cvvs)
  nloo_tst <- nobsv - 1L
  tstsetups <- grep("^glm\\.gauss\\..*\\.default_cvmeth", names(cvvss),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    tstsetup_ref <- args_cvvs_i$tstsetup_ref
    mod_crr <- args_cvvs_i$mod_nm
    fam_crr <- args_cvvs_i$fam_nm
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    # Use SW() because of occasional warnings concerning Pareto k diagnostics:
    SW(cvvs_nloo <- do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]],
           nloo = nloo_tst),
      excl_nonargs(args_cvvs_i)
    )))
    vsel_tester(
      cvvs_nloo,
      with_cv = TRUE,
      refmod_expected = refmods[[tstsetup_ref]],
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "LOO",
      valsearch_expected = args_cvvs_i$validate_search,
      nclusters_expected = args_cvvs_i$nclusters,
      nclusters_pred_expected = args_cvvs_i$nclusters_pred,
      nloo_expected = nloo_tst,
      info_str = tstsetup
    )
    # Expected equality for most components with a few exceptions:
    expect_equal(cvvs_nloo[setdiff(vsel_nms_cv, vsel_nms_cv_nloo)],
                 cvvss[[tstsetup]][setdiff(vsel_nms_cv, vsel_nms_cv_nloo)],
                 info = tstsetup)
    # Expected inequality for the exceptions (but note that the components from
    # `vsel_nms_cv_nloo_opt` can be, but don't need to be differing):
    for (vsel_nm in setdiff(vsel_nms_cv_nloo, vsel_nms_cv_nloo_opt)) {
      expect_false(isTRUE(all.equal(cvvs_nloo[[vsel_nm]],
                                    cvvss[[tstsetup]][[vsel_nm]])),
                   info = paste(tstsetup, vsel_nm, sep = "__"))
    }
  }
})

test_that("`validate_search` works", {
  skip_if_not(run_cvvs)
  tstsetups <- grep("\\.default_cvmeth", names(cvvss), value = TRUE)
  if (!run_valsearch_always) {
    tstsetups <- grep("^glm\\.", tstsetups, value = TRUE)
    tstsetups <- grep("\\.forward", tstsetups, value = TRUE, invert = TRUE)
  }
  suggsize_cond <- setNames(rep(NA, length(tstsetups)), nm = tstsetups)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    stopifnot(is.null(args_cvvs_i$validate_search) ||
                isTRUE(args_cvvs_i$validate_search))
    tstsetup_ref <- args_cvvs_i$tstsetup_ref
    mod_crr <- args_cvvs_i$mod_nm
    fam_crr <- args_cvvs_i$fam_nm
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    # Use SW() because of occasional warnings concerning Pareto k diagnostics:
    SW(cvvs_valsearch <- do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]],
           validate_search = FALSE),
      excl_nonargs(args_cvvs_i)
    )))
    vsel_tester(
      cvvs_valsearch,
      with_cv = TRUE,
      refmod_expected = refmods[[tstsetup_ref]],
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "LOO",
      valsearch_expected = FALSE,
      nclusters_expected = args_cvvs_i$nclusters,
      nclusters_pred_expected = args_cvvs_i$nclusters_pred,
      info_str = tstsetup
    )
    # Expected equality for most components with a few exceptions:
    expect_equal(cvvs_valsearch[setdiff(vsel_nms_cv, vsel_nms_cv_valsearch)],
                 cvvss[[tstsetup]][setdiff(vsel_nms_cv, vsel_nms_cv_valsearch)],
                 info = tstsetup)
    expect_identical(cvvs_valsearch$summaries$ref,
                     cvvss[[tstsetup]]$summaries$ref,
                     info = tstsetup)
    # Expected inequality for the exceptions (but note that the components from
    # `vsel_nms_cv_valsearch_opt` can be, but don't need to be differing):
    for (vsel_nm in setdiff(vsel_nms_cv_valsearch, vsel_nms_cv_valsearch_opt)) {
      expect_false(isTRUE(all.equal(cvvs_valsearch[[vsel_nm]],
                                    cvvss[[tstsetup]][[vsel_nm]])),
                   info = paste(tstsetup, vsel_nm, sep = "__"))
    }
    # Check the expected inequalities more specifically:
    # Without a validated search, we expect increased LPPDs (and consequently
    # also an increased ELPD) in the submodels (since the hold-out fold was
    # included in the dataset for fitting the submodels):
    tol_crr <- 5e-2
    # Allow for just a small proportion of extreme differences:
    prop_as_expected <- 0.9
    for (j in seq_along(cvvs_valsearch$summaries$sub)) {
      expect_true(mean(cvvs_valsearch$summaries$sub[[j]]$lppd >=
                         cvvss[[tstsetup]]$summaries$sub[[j]]$lppd - tol_crr) >=
                    prop_as_expected,
                  info = paste(tstsetup, j, sep = "__"))
    }
    expect_true(all(cvvs_valsearch$summary$elpd.loo >=
                      cvvss[[tstsetup]]$summary$elpd.loo),
                info = tstsetup)
    # Without a validated search, we expect overfitting in the suggested model
    # size:
    if (!is.na(cvvs_valsearch$suggested_size) &
        !is.na(cvvss[[tstsetup]]$suggested_size)) {
      suggsize_cond[tstsetup] <- cvvs_valsearch$suggested_size >=
        cvvss[[tstsetup]]$suggested_size
    }
  }
  sum_as_unexpected <- 0L
  expect_true(sum(!suggsize_cond, na.rm = TRUE) <= sum_as_unexpected)
})

test_that("invalid `K` fails", {
  skip_if_not(run_cvvs_kfold)
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = 1),
               "^K must be at least 2$")
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = 1000),
               "^K cannot exceed n$")
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = c(4, 9)),
               "^K must be a single integer value$")
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = "a"),
               "^K must be a single integer value$")
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = dat),
               "^K must be a single integer value$")
})

test_that("`cvfits` (actually passed to init_refmodel()) works", {
  skip_if_not(run_cvvs_kfold)
  tstsetups <- grep("kfold", names(cvvss), value = TRUE)
  if (!run_cvfits_all) {
    tstsetups <- head(tstsetups, 1)
  }
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    tstsetup_fit <- args_cvvs_i$tstsetup_fit
    mod_crr <- args_cvvs_i$mod_nm
    fam_crr <- args_cvvs_i$fam_nm
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    fit_crr <- fits[[tstsetup_fit]]
    K_crr <- args_cvvs_i$K

    # Refit `K_crr` times:
    # rstanarm::kfold() lacks an argument for setting the seed:
    set.seed(seed_tst)
    # Additionally to SW(), suppressMessages() could be used here (but is not
    # necessary since messages seem to be suppressed within test_that()'s
    # `code`):
    SW(kfold_obj <- rstanarm::kfold(fit_crr, K = K_crr, save_fits = TRUE))

    # Create the folds vector:
    folds_vec <- rep(NA, nobsv)
    for (k_crr in seq_len(K_crr)) {
      idcs_fold <- kfold_obj$fits[, "omitted"][[k_crr]]
      stopifnot(identical(
        idcs_fold,
        setdiff(seq_len(nobsv),
                as.integer(rownames(kfold_obj$fits[, "fit"][[k_crr]]$data)))
      ))
      folds_vec[idcs_fold] <- k_crr
    }
    stopifnot(all(!is.na(folds_vec)))
    attr(kfold_obj, "folds") <- folds_vec

    # Create `"refmodel"` object with `cvfits`:
    # For the binomial family with > 1 trials, we currently expect the warning
    # "Using formula(x) is deprecated when x is a character vector of length >
    # 1" (see GitHub issue #136):
    warn_expected_binom <- switch(
      fam_crr,
      "binom" = paste("Using formula\\(x\\) is deprecated when x is a",
                      "character vector of length > 1"),
      NA
    )
    expect_warning(
      refmod_crr <- get_refmodel(fit_crr, cvfits = kfold_obj),
      warn_expected_binom,
      info = tstsetup
    )

    # Run cv_varsel():
    # We expect a warning which in fact should be suppressed, though (see
    # issue #162):
    warn_expected <- switch(
      mod_crr,
      "glm" = paste("^'offset' argument is NULL but it looks like you",
                    "estimated the model using an offset term\\.$"),
      NA
    )
    if (fam_crr == "binom") {
      warn_expected <- paste(c(warn_expected_binom, warn_expected),
                             collapse = "|")
    }
    expect_warning(
      cvvs_cvfits <- do.call(cv_varsel, c(
        list(object = refmod_crr),
        excl_nonargs(args_cvvs_i, nms_excl_add = "K")
      )),
      warn_expected,
      info = tstsetup
    )

    # Checks:
    vsel_tester(
      cvvs_cvfits,
      with_cv = TRUE,
      refmod_expected = refmod_crr,
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "kfold",
      valsearch_expected = args_cvvs_i$validate_search,
      nclusters_expected = args_cvvs_i$nclusters,
      nclusters_pred_expected = args_cvvs_i$nclusters_pred,
      info_str = tstsetup
    )
    # Note: Unfortunately, it is currently not possible to always ensure exactly
    # the same seed when performing K-fold CV with `cvfits` or without `cvfits`.
    # Therefore, the following checks for equality/inequality are quite
    # restricted.
    # Expected equality for some components:
    expect_equal(cvvs_cvfits[setdiff(vsel_nms_cv, vsel_nms_cv_cvfits)],
                 cvvss[[tstsetup]][setdiff(vsel_nms_cv, vsel_nms_cv_cvfits)],
                 info = tstsetup)
    # Expected inequality for the remaining components (but note that the
    # components from `vsel_nms_cv_cvfits_opt` can be, but don't need to be
    # differing):
    for (vsel_nm in setdiff(vsel_nms_cv_cvfits, vsel_nms_cv_cvfits_opt)) {
      expect_false(isTRUE(all.equal(cvvs_cvfits[[vsel_nm]],
                                    cvvss[[tstsetup]][[vsel_nm]])),
                   info = paste(tstsetup, vsel_nm, sep = "__"))
    }
  }
})
