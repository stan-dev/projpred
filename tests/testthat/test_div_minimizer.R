context("div_minimizer")

# Setup -------------------------------------------------------------------

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

# divmin() ----------------------------------------------------------------

test_that("divmin() works", {
  skip_on_cran()
  # For comparison with the divmin_augdat() test:
  divmin_res_brnll_tmp <- list()

  for (tstsetup in names(fits)) {
    args_fit_i <- args_fit[[tstsetup]]
    pkg_crr <- args_fit_i$pkg_nm
    mod_crr <- args_fit_i$mod_nm
    fam_crr <- args_fit_i$fam_nm

    if (fam_crr == "gauss") {
      mu_crr <- posterior_linpred(fits[[tstsetup]])
      # Crude approximations for the predictive variances:
      var_crr <- mean(as.matrix(fits[[tstsetup]])[, "sigma"]^2) +
        (colMeans(mu_crr^2) - colMeans(mu_crr)^2)
    } else {
      var_crr <- rep(NA, nobsv)
    }
    args_fit_i$projpred_var <- matrix(var_crr)
    args_fit_i$projpred_regul <- regul_default

    if (pkg_crr == "brms" && grepl("\\.with_wobs", tstsetup)) {
      args_fit_i$formula <- rm_addresp(args_fit_i$formula)
      args_fit_i$weights <- eval(args_fit_i$data)$wobs_col
    }

    if (fam_crr == "binom") {
      # Use the proportion of successes as (1-column) response and the number of
      # trials as `weights`:
      ybinprop_nm <- paste("ybinprop", mod_crr, sep = "_")
      args_fit_i$data <- within(dat, {
        assign(ybinprop_nm,
               get(paste("y", mod_crr, fam_crr, sep = "_")) / wobs_col)
      })
      args_fit_i$formula <- update(
        args_fit_i$formula,
        as.formula(paste(ybinprop_nm, "~ ."))
      )
      args_fit_i$weights <- wobs_tst
    } else if (fam_crr == "brnll") {
      if (pkg_crr == "brms") {
        args_fit_i$family <- f_binom
      }
    }

    if ("random" %in% names(args_fit_i)) {
      args_fit_i$formula <- update(
        args_fit_i$formula,
        as.formula(paste(". ~ . +", tail(as.character(args_fit_i$random), 1)))
      )
    }

    divmin_res <- do.call(
      divmin,
      args_fit_i[intersect(c("formula", "data", "family", "weights",
                             "projpred_var", "projpred_regul"),
                           names(args_fit_i))]
    )
    if (fam_crr == "brnll") {
      divmin_res_brnll_tmp[[tstsetup]] <- divmin_res
    }

    if (fam_crr == "binom" || grepl("\\.with_wobs", tstsetup)) {
      wobs_expected_crr <- wobs_tst
    } else {
      wobs_expected_crr <- NULL
    }
    submodl_tester(divmin_res,
                   nprjdraws_expected = 1L,
                   sub_formul = list(args_fit_i$formula),
                   sub_data = eval(args_fit_i$data),
                   sub_fam = eval(args_fit_i$family)$family,
                   wobs_expected = wobs_expected_crr,
                   with_offs = grepl("\\.with_offs", tstsetup),
                   info_str = tstsetup)
  }
  assign("divmin_res_brnll_trad", divmin_res_brnll_tmp, envir = .GlobalEnv)
})

# divmin_augdat() ---------------------------------------------------------

test_that("divmin_augdat() works", {
  skip_on_cran()
  # For comparison with the divmin() test:
  divmin_res_brnll_tmp <- list()

  tstsetups_ref <- names(args_ref[sapply(args_ref, "[[", "prj_nm") == "augdat"])
  for (tstsetup in tstsetups_ref) {
    args_ref_i <- args_ref[[tstsetup]]
    args_fit_i <- args_fit[[args_ref_i$tstsetup_fit]]
    pkg_crr <- args_fit_i$pkg_nm
    mod_crr <- args_fit_i$mod_nm
    fam_crr <- args_fit_i$fam_nm
    y_nm <- as.character(args_fit_i$formula)[2]
    y_unqs <- refmods[[tstsetup]]$family$cats
    ncats_crr <- length(y_unqs)

    # Augmented-data weights, here indicating the observed response category
    # with value 1 and all others with 0:
    args_fit_i$projpred_ws_aug <- lapply(y_unqs, function(y_unq) {
      as.integer(eval(args_fit_i$data)[[y_nm]] == as.numeric(y_unq))
    })
    args_fit_i$projpred_ws_aug <- do.call(cbind, args_fit_i$projpred_ws_aug)
    stopifnot(all(rowSums(args_fit_i$projpred_ws_aug) == 1))
    args_fit_i$projpred_ws_aug <- matrix(args_fit_i$projpred_ws_aug)

    # Create the augmented dataset:
    args_fit_i$data <- do.call(rbind, lapply(y_unqs, function(y_unq) {
      dat_j <- eval(args_fit_i$data)
      dat_j[[y_nm]] <- y_unq
      return(dat_j)
    }))
    args_fit_i$data[[y_nm]] <- factor(args_fit_i$data[[y_nm]], levels = y_unqs)

    # Other arguments:
    args_fit_i$weights <- "currently unused"
    args_fit_i$projpred_var <- matrix(nrow = nobsv * ncats_crr)
    args_fit_i$projpred_regul <- regul_default
    if (fam_crr == "brnll" && pkg_crr == "brms") {
      args_fit_i$family <- f_binom
    }

    divmin_res <- do.call(
      divmin_augdat,
      args_fit_i[intersect(c("formula", "data", "family", "weights",
                             "projpred_var", "projpred_regul",
                             "projpred_ws_aug"),
                           names(args_fit_i))]
    )
    if (fam_crr == "brnll") {
      divmin_res_brnll_tmp[[tstsetup]] <- divmin_res
    }

    submodl_tester(divmin_res,
                   nprjdraws_expected = 1L,
                   sub_formul = list(args_fit_i$formula),
                   sub_data = args_fit_i$data,
                   sub_fam = eval(args_fit_i$family)$family,
                   wobs_expected = args_fit_i$projpred_ws_aug,
                   with_offs = grepl("\\.with_offs", tstsetup),
                   augdat_cats = y_unqs,
                   allow_w_zero = TRUE,
                   check_y_from_resp = FALSE,
                   info_str = tstsetup)
  }
  assign("divmin_res_brnll_augdat", divmin_res_brnll_tmp, envir = .GlobalEnv)
})

# Comparison of divmin() and divmin_augdat() ------------------------------

test_that(paste(
  "divmin() and divmin_augdat() give (almost) the same results (in case of the",
  "`brnll` family)"
), {
  skip_on_cran()
  for (tstsetup in names(divmin_res_brnll_augdat)) {
    args_ref_i <- args_ref[[tstsetup]]
    submodl_augdat <- divmin_res_brnll_augdat[[tstsetup]]
    submodl_trad <- divmin_res_brnll_trad[[args_ref_i$tstsetup_fit]]
    expect_length(submodl_augdat, 1)
    expect_length(submodl_trad, 1)
    tol_coefs <- ifelse(args_ref_i$mod_nm == "glmm", 1e-4, 1e-5)
    expect_equal(get_subparams(submodl_augdat[[1]], nm_scheme = "rstanarm"),
                 get_subparams(submodl_trad[[1]], nm_scheme = "rstanarm"),
                 tolerance = tol_coefs, info = tstsetup)
  }
})

# Teardown ----------------------------------------------------------------

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
