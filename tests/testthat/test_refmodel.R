# get_refmodel() ----------------------------------------------------------

context("get_refmodel()")

test_that("`object` of class \"stanreg\" works", {
  for (tstsetup in names(refmods)) {
    tstsetup_fit <- args_ref[[tstsetup]]$tstsetup_fit
    if (!grepl("\\.spclformul", tstsetup)) {
      needs_y_overwrite_crr <- FALSE
    } else {
      needs_y_overwrite_crr <- TRUE
    }
    if (args_ref[[tstsetup]]$fam_nm == "binom" ||
        grepl("\\.with_wobs", tstsetup)) {
      wobs_expected_crr <- wobs_tst
    } else {
      wobs_expected_crr <- rep(1, nobsv)
    }
    refmodel_tester(
      refmods[[tstsetup]],
      fit_expected = fits[[tstsetup_fit]],
      needs_y_overwrite = needs_y_overwrite_crr,
      wobs_expected = wobs_expected_crr,
      info_str = tstsetup,
      fam_orig = eval(args_fit[[tstsetup_fit]]$family)
    )
  }
})

test_that("missing `data` fails", {
  SW(fit_nodata <- rstanarm::stan_glm(
    dat$y_glm_gauss ~ dat$xco.1 + dat$xco.2 + dat$xco.3 + dat$xca.1 + dat$xca.2,
    family = f_gauss,
    weights = dat$wobs_col, offset = dat$offs_col,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  ))
  expect_error(get_refmodel(fit_nodata),
               "^object of type 'environment' is not subsettable$")
})

test_that("`formula` as a character string fails", {
  # If `formula` is a character string, rstanarm::stan_glm() is not able to find
  # objects supplied to arguments `weights` or `offset`, at least when using
  # devtools::test():
  SW(fit_str <- rstanarm::stan_glm(
    "y_glm_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2",
    family = f_gauss, data = dat,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  ))
  expect_error(get_refmodel(fit_str),
               "^inherits\\(formula, \"formula\"\\) is not TRUE$")
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
               "^ynew must be a numerical vector$")
})

test_that(paste(
  "`object` of class `\"refmodel\"`, `newdata`, `ynew`, and `type` work"
), {
  for (tstsetup in names(refmods)) {
    mod_crr <- args_ref[[tstsetup]]$mod_nm
    fam_crr <- args_ref[[tstsetup]]$fam_nm

    # We expect a warning which in fact should be suppressed, though (see
    # issue #162):
    warn_expected <- switch(
      mod_crr,
      "glm" = paste("^'offset' argument is NULL but it looks like you",
                    "estimated the model using an offset term\\.$"),
      NA
    )

    y_crr <- dat[, paste("y", mod_crr, fam_crr, sep = "_")]
    if (is.integer(y_crr)) {
      y_crr <- as.numeric(y_crr)
    }

    # Without `ynew`:
    expect_warning(
      predref_resp <- predict(refmods[[tstsetup]], dat, type = "response"),
      warn_expected,
      info = tstsetup
    )
    expect_warning(
      predref_link <- predict(refmods[[tstsetup]], dat, type = "link"),
      warn_expected,
      info = tstsetup
    )

    # With `ynew`:
    expect_warning(
      predref_ynew_resp <- predict(refmods[[tstsetup]], dat, ynew = y_crr,
                                   type = "response"),
      warn_expected,
      info = tstsetup
    )
    expect_warning(
      predref_ynew_link <- predict(refmods[[tstsetup]], dat, ynew = y_crr,
                                   type = "link"),
      warn_expected,
      info = tstsetup
    )

    # Checks without `ynew`:
    expect_true(is.vector(predref_resp, "double"), info = tstsetup)
    expect_length(predref_resp, nobsv)
    if (fam_crr == "binom") {
      expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                  info = tstsetup)
    }
    expect_true(is.vector(predref_link, "double"), info = tstsetup)
    expect_length(predref_link, nobsv)
    if (fam_crr == "gauss") {
      expect_equal(predref_resp, predref_link, info = tstsetup)
    }

    # Checks with `ynew`:
    expect_equal(predref_ynew_resp, predref_ynew_link, info = tstsetup)
    expect_true(is.vector(predref_ynew_resp, "double"), info = tstsetup)
    expect_length(predref_ynew_resp, nobsv)
    expect_false(isTRUE(all.equal(predref_ynew_resp, predref_resp)),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(predref_ynew_resp, predref_link)),
                 info = tstsetup)
  }
})
