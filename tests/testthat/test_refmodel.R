# get_refmodel() ----------------------------------------------------------

context("get_refmodel()")

test_that("`object` of class \"stanreg\" or \"brmsfit\" works", {
  for (tstsetup in names(refmods)) {
    tstsetup_fit <- args_ref[[tstsetup]]$tstsetup_fit
    with_spclformul_crr <- grepl("\\.spclformul", tstsetup)
    if (args_ref[[tstsetup]]$fam_nm == "binom" ||
        grepl("\\.with_wobs", tstsetup)) {
      wobs_expected_crr <- wobs_tst
    } else {
      wobs_expected_crr <- rep(1, nobsv)
    }
    if (grepl("\\.with_offs", tstsetup)) {
      offs_expected_crr <- offs_tst
    } else {
      offs_expected_crr <- rep(0, nobsv)
    }
    if (args_ref[[tstsetup]]$prj_nm == "latent") {
      fam_orig_expected <- f_gauss
    } else if (args_ref[[tstsetup]]$pkg_nm == "brms" &&
               args_ref[[tstsetup]]$fam_nm %in% fam_nms_aug) {
      fam_orig_expected <- eval(args_fit[[tstsetup_fit]]$family)
    } else {
      if (args_ref[[tstsetup]]$pkg_nm == "rstanarm" &&
          args_ref[[tstsetup]]$fam_nm == "cumul") {
        f_cumul <- get_f_cumul()
      }
      fam_orig_expected <- get(paste0("f_", args_fit[[tstsetup_fit]]$fam_nm))
    }
    refmodel_tester(
      refmods[[tstsetup]],
      pkg_nm = args_ref[[tstsetup]]$pkg_nm,
      fit_expected = fits[[tstsetup_fit]],
      with_spclformul = with_spclformul_crr,
      wobs_expected = wobs_expected_crr,
      offs_expected = offs_expected_crr,
      fam_orig = fam_orig_expected,
      mod_nm = args_ref[[tstsetup]]$mod_nm,
      fam_nm = args_ref[[tstsetup]]$fam_nm,
      augdat_expected = args_ref[[tstsetup]]$prj_nm == "augdat",
      latent_expected = args_ref[[tstsetup]]$prj_nm == "latent",
      info_str = tstsetup
    )
  }
})

test_that("missing `data` fails", {
  fit_nodata <- suppressWarnings(rstanarm::stan_glm(
    dat$y_glm_gauss ~ dat$xco.1 + dat$xco.2 + dat$xco.3 +
      dat$xca.1 + dat$xca.2 + offset(dat$offs_col),
    family = f_gauss,
    weights = dat$wobs_col,
    chains = chains_tst, seed = seed_fit, iter = iter_tst, QR = TRUE,
    refresh = 0
  ))
  expect_error(
    get_refmodel(fit_nodata),
    paste("^is\\.data\\.frame\\(data\\) is not TRUE$")
  )
})

test_that("`formula` as a character string fails", {
  # If `formula` is a character string, rstanarm::stan_glm() is not able to find
  # objects supplied to arguments `weights` or `offset`, at least when using
  # devtools::test():
  fit_str <- suppressWarnings(rstanarm::stan_glm(
    "y_glm_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2",
    family = f_gauss, data = dat,
    chains = chains_tst, seed = seed_fit, iter = iter_tst, QR = TRUE,
    refresh = 0
  ))
  expect_error(get_refmodel(fit_str),
               "^inherits\\(formula, \"formula\"\\) is not TRUE$")
})

test_that("offsets specified via argument `offset` work", {
  args_fit_i <- args_fit$rstanarm.glm.gauss.stdformul.with_wobs.with_offs
  skip_if_not(!is.null(args_fit_i))
  fit_fun_nm <- switch(args_fit_i$pkg_nm,
                       "rstanarm" = switch(args_fit_i$mod_nm,
                                           "glm" = "stan_glm",
                                           "glmm" = "stan_glmer",
                                           "stan_gamm4"),
                       "brms" = "brm",
                       stop("Unknown `pkg_nm`."))
  upd_no_offs <- paste(". ~", sub(" \\+ offset\\(offs_col\\)", "",
                                  as.character(args_fit_i$formula[3])))
  fit_offs_arg <- suppressWarnings(do.call(
    get(fit_fun_nm, asNamespace(args_fit_i$pkg_nm)),
    c(list(formula = update(args_fit_i$formula, upd_no_offs),
           offset = offs_tst),
      excl_nonargs(args_fit_i, nms_excl_add = "formula"))
  ))
  refmod_offs_arg <- get_refmodel(fit_offs_arg)
  expect_equal(
    as.matrix(fit_offs_arg),
    as.matrix(fits$rstanarm.glm.gauss.stdformul.with_wobs.with_offs),
    tolerance = 1e-12,
    info = "rstanarm.glm.gauss.stdformul.with_wobs.with_offs"
  )
  nms_compare <- c("eta", "mu", "mu_offs", "dis", "y", "intercept", "wobs",
                   "wsample", "offset", "y_oscale")
  expect_equal(
    refmod_offs_arg[nms_compare],
    refmods$rstanarm.glm.gauss.stdformul.with_wobs.with_offs[nms_compare],
    tolerance = 1e-12,
    info = "rstanarm.glm.gauss.stdformul.with_wobs.with_offs"
  )
})

test_that(paste(
  "binomial family with 1-column response and weights which are not all ones",
  "errors"
), {
  dat_prop <- within(dat, {
    ybinprop_glm <- y_glm_binom / wobs_col
  })
  fit_binom_1col_wobs <- suppressWarnings(rstanarm::stan_glm(
    ybinprop_glm ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 + offset(offs_col),
    family = f_binom, data = dat_prop,
    weights = wobs_tst,
    chains = chains_tst, seed = seed_fit, iter = iter_tst, QR = TRUE,
    refresh = 0
  ))
  if ("rstanarm.glm.binom.stdformul.without_wobs.with_offs" %in% names(fits)) {
    expect_equal(
      as.matrix(fit_binom_1col_wobs),
      as.matrix(fits$rstanarm.glm.binom.stdformul.without_wobs.with_offs)
    )
  }
  expect_error(get_refmodel(fit_binom_1col_wobs),
               "response must contain numbers of successes")
})

test_that("get_refmodel() is idempotent", {
  for (tstsetup in names(refmods)) {
    expect_identical(get_refmodel(refmods[[tstsetup]]),
                     refmods[[tstsetup]],
                     info = tstsetup)
  }
})

# predict.refmodel() ------------------------------------------------------

context("predict.refmodel()")

test_that("invalid `type` fails", {
  expect_error(predict(refmods[[1]], dat, type = "zzz"),
               "^type should be one of")
})

test_that("invalid `ynew` fails", {
  expect_error(predict(refmods[[1]], dat, ynew = dat),
               "^Argument `ynew` must be a numeric vector\\.$")
})

test_that(paste(
  "`object` of class `\"refmodel\"`, `newdata`, `ynew`, and `type` work"
), {
  for (tstsetup in names(refmods)) {
    pkg_crr <- args_ref[[tstsetup]]$pkg_nm
    mod_crr <- args_ref[[tstsetup]]$mod_nm
    fam_crr <- args_ref[[tstsetup]]$fam_nm
    prj_crr <- args_ref[[tstsetup]]$prj_nm

    y_crr <- dat[, paste("y", mod_crr, fam_crr, sep = "_")]
    if (prj_crr == "latent") {
      dat_crr <- dat
      if (pkg_crr == "rstanarm" && grepl("\\.with_offs\\.", tstsetup)) {
        dat_crr$projpred_internal_offs_stanreg <- 0
      }
      y_nm <- stdize_lhs(refmods[[tstsetup]]$formula)$y_nm
      # Use `ref_predfun_usr` here (instead of `ref_predfun`) to include
      # offsets:
      refprd_with_offs <- get(
        "ref_predfun_usr",
        envir = environment(refmods[[tstsetup]]$ref_predfun)
      )
      y_crr_link <- rowMeans(unname(
        refprd_with_offs(fit = refmods[[tstsetup]]$fit, newdata = dat_crr)
      ))
    } else {
      y_crr_link <- y_crr
    }

    # Without `ynew`:
    predref_resp <- predict(refmods[[tstsetup]], dat, type = "response")
    predref_link <- predict(refmods[[tstsetup]], dat, type = "link")

    # With `ynew`:
    predref_ynew_resp <- predict(refmods[[tstsetup]], dat, ynew = y_crr,
                                 type = "response")
    predref_ynew_link <- predict(refmods[[tstsetup]], dat, ynew = y_crr_link,
                                 type = "link")

    # Checks without `ynew`:
    if (prj_crr %in% c("latent", "augdat")) {
      if (prj_crr == "augdat" || !is.null(refmods[[tstsetup]]$family$cats)) {
        expect_identical(dim(predref_resp),
                         c(nobsv, length(refmods[[tstsetup]]$family$cats)),
                         info = tstsetup)
        expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                    info = tstsetup)
      } else {
        expect_true(is.vector(predref_resp, "double"), info = tstsetup)
        expect_length(predref_resp, nobsv)
        if (fam_crr %in% c("brnll", "binom")) {
          expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                      info = tstsetup)
        }
      }
      if (prj_crr == "augdat") {
        expect_identical(dim(predref_link),
                         c(nobsv, length(refmods[[tstsetup]]$family$cats) - 1L),
                         info = tstsetup)
      } else if (prj_crr == "latent") {
        expect_true(is.vector(predref_link, "double"), info = tstsetup)
        expect_length(predref_link, nobsv)
      }
    } else {
      expect_true(is.vector(predref_resp, "double"), info = tstsetup)
      expect_length(predref_resp, nobsv)
      if (fam_crr %in% c("brnll", "binom")) {
        expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                    info = tstsetup)
      }
      expect_true(is.vector(predref_link, "double"), info = tstsetup)
      expect_length(predref_link, nobsv)
      if (fam_crr == "gauss") {
        expect_equal(predref_resp, predref_link, info = tstsetup)
      }
    }

    # Checks with `ynew`:
    if (prj_crr != "latent") {
      expect_equal(predref_ynew_resp, predref_ynew_link, info = tstsetup)
    } else {
      expect_false(isTRUE(all.equal(predref_ynew_resp, predref_ynew_link)),
                   info = tstsetup)
    }
    expect_true(is.vector(predref_ynew_resp, "double"), info = tstsetup)
    expect_length(predref_ynew_resp, nobsv)
    expect_false(isTRUE(all.equal(predref_ynew_resp, predref_resp)),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(predref_ynew_resp, predref_link)),
                 info = tstsetup)

    # Snapshots:
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(rlang::hash(predref_resp))
        print(rlang::hash(predref_link))
        print(rlang::hash(predref_ynew_resp))
        print(rlang::hash(predref_ynew_link))
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})
