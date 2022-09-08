context("div_minimizer")

# Setup -------------------------------------------------------------------

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

# divmin() ----------------------------------------------------------------

test_that("divmin() works", {
  skip_on_cran()
  # For comparison with the divmin_augdat() test:
  divmin_res_brnll_tmp <- list()

  for (tstsetup in grep(fam_nms_aug_regex, names(fits), value = TRUE,
                        invert = TRUE)) {
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
        # Remove some terms which lead to extreme coefficients:
        args_fit_i$formula <- update(args_fit_i$formula, . ~ . - xca.1)
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
  tstsetups_ref <- grep(fam_nms_unsupp_regex, tstsetups_ref, value = TRUE,
                        invert = TRUE)
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
      if (fam_crr == "brnll") {
        y_unq <- as.numeric(y_unq)
      }
      as.integer(eval(args_fit_i$data)[[y_nm]] == y_unq)
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
      # Remove some terms which lead to extreme coefficients:
      args_fit_i$formula <- update(args_fit_i$formula, . ~ . - xca.1)
    } else if (fam_crr == "cumul" && pkg_crr == "rstanarm") {
      args_fit_i$family <- structure(
        list(family = "cumulative_rstanarm",
             link = refmods[[tstsetup]]$family$link,
             cats = y_unqs),
        class = "family"
      )
    } else if (fam_crr == "cumul" && pkg_crr == "brms") {
      fam_crr_long <- get_fam_long(fam_crr)
      args_fit_i$family <- substitute(
        get(fam_crr_long_subst, envir = asNamespace("brms"))(),
        list(fam_crr_long_subst = fam_crr_long)
      )
    } else if (fam_crr == "categ" && mod_crr == "glmm") {
      # Quick-and-dirty solution to get some working results (it's probably due
      # to unfortunate test data simulated here that convergence at the default
      # settings is not given):
      args_fit_i <- c(args_fit_i, list(avoid.increase = TRUE))
    }

    if (fam_crr == "cumul" && mod_crr %in% c("glmm", "gamm")) {
      warn_expected <- paste(
        "^Using formula\\(x\\) is deprecated when x is a character vector of",
        "length > 1"
      )
    } else if (fam_crr == "categ" && mod_crr == "glmm") {
      warn_expected <- paste0(
        "^step size truncated due to possible divergence$|",
        "^Algorithm stopped due to false convergence$"
      )
    } else {
      warn_expected <- NA
    }
    expect_warning(
      divmin_res <- do.call(
        divmin_augdat,
        args_fit_i[intersect(c("formula", "data", "family", "weights",
                               "projpred_var", "projpred_regul",
                               "projpred_ws_aug", "epsilon", "avoid.increase"),
                             names(args_fit_i))]
      ),
      warn_expected
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
    coefs_aug <- get_subparams(submodl_augdat[[1]], nm_scheme = "rstanarm")
    coefs_trad <- get_subparams(submodl_trad[[1]], nm_scheme = "rstanarm")
    if (!identical(is.na(coefs_aug), is.na(coefs_trad))) {
      # There is at least one case where the traditional projection is not able
      # to estimate a correlation of -1:
      expect_true(all(coefs_aug[is.na(coefs_trad)] %in% c(-1, 1)),
                  info = tstsetup)
      coefs_aug <- coefs_aug[!is.na(coefs_trad)]
      coefs_trad <- coefs_trad[!is.na(coefs_trad)]
    }
    expect_equal(coefs_aug, coefs_trad, tolerance = tol_coefs,
                 info = tstsetup)
  }
})

# Teardown ----------------------------------------------------------------

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
