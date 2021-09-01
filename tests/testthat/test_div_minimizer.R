context("div_minimizer")

test_that("all div_minimizer()s work", {
  for (tstsetup in names(fits)) {
    args_fit_i <- args_fit[[tstsetup]]
    mod_crr <- args_fit_i$mod_nm
    fam_crr <- args_fit_i$fam_nm

    if (mod_crr == "glm") {
      divmin_fun <- "linear_mle"
    } else if (mod_crr == "glmm") {
      divmin_fun <- "linear_multilevel_mle"
    } else if (mod_crr %in% c("gam", "gamm")) {
      divmin_fun <- "additive_mle"
    }

    if (fam_crr == "gauss") {
      mu_crr <- posterior_linpred(fits[[tstsetup]])
      # Crude approximations for the predictive variances:
      var_crr <- mean(as.matrix(fits[[tstsetup]])[, "sigma"]^2) +
        (colMeans(mu_crr^2) - colMeans(mu_crr)^2)
      var_crr <- median(var_crr)
    } else {
      var_crr <- 0
    }

    formul_chr <- as.character(args_fit_i$formula)
    ### Note: Here, the offset term is *not* removed (as done in projpred when
    ### fitting submodels) because then, the response would have to be adapted,
    ### too. The code for omitting the offsets would be:
    # if (grepl("[[:blank:]]*\\+[[:blank:]]*offset\\(.*\\)", formul_chr[3])) {
    #   formul_chr[3] <- gsub("[[:blank:]]*\\+[[:blank:]]*offset\\(.*\\)", "",
    #                         formul_chr[3])
    # }
    ###
    # For "brmsfit"s, remove additional response information:
    if (args_fit_i$pkg_nm == "brms" &&
        grepl("[[:blank:]]*\\|[[:blank:]]*weights\\(wobs_col\\)$",
              formul_chr[2])) {
      formul_chr[2] <- gsub("[[:blank:]]*\\|[[:blank:]]*weights\\(wobs_col\\)$",
                            "", formul_chr[2])
      args_fit_i$weights <- eval(args_fit_i$data)$wobs_col
    }
    args_fit_i$formula <- update(args_fit_i$formula,
                                 paste(formul_chr[c(2, 1, 3)], collapse = " "))

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
    }

    if ("random" %in% names(args_fit_i)) {
      args_fit_i$formula <- update(
        args_fit_i$formula,
        as.formula(paste(". ~ . +", tail(as.character(args_fit_i$random), 1)))
      )
    }

    divmin <- do.call(divmin_fun, c(
      args_fit_i[intersect(names(args_fit_i),
                           c("formula", "data", "family", "weights"))],
      list(regul = regul_default, var = var_crr)
    ))

    if (fam_crr == "binom" || grepl("\\.with_wobs", tstsetup)) {
      wobs_expected_crr <- wobs_tst
    } else {
      wobs_expected_crr <- NULL
    }
    sub_fit_tester(divmin,
                   nprjdraws_expected = 1L,
                   sub_formul = list(args_fit_i$formula),
                   sub_data = eval(args_fit_i$data),
                   sub_fam = eval(args_fit_i$family)$family,
                   wobs_expected = wobs_expected_crr,
                   with_offs = TRUE,
                   info_str = tstsetup)
  }
})
