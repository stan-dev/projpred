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

    if (fam_crr == "binom") {
      ybinprop_nm <- paste("ybinprop", mod_crr, sep = "_")
      args_fit_i$formula <- update(
        args_fit_i$formula,
        as.formula(paste(ybinprop_nm, "~ ."))
      )
    }

    divmin <- do.call(divmin_fun, c(
      args_fit_i[c("formula", "data", "family", "weights")],
      list(regul = regul_default, var = var_crr)
    ))

    sub_fit_tester(divmin,
                   nprjdraws_expected = 1L,
                   sub_formul = args_fit_i$formula,
                   sub_data = eval(args_fit_i$data),
                   sub_fam = eval(args_fit_i$family)$family,
                   info_str = tstsetup)
    # TODO: Add more expectations for special formulas.
  }
})
