# get_refmodel() ----------------------------------------------------------

context("get_refmodel()")

test_that("`object` of class \"stanreg\" works", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      refmodel_tester(
        refmod = refmods[[mod_nm]][[fam_nm]],
        fit_expected = fits[[mod_nm]][[fam_nm]],
        nobsv_expected = nobsv,
        nrefdraws_expected = chains_tst * (iter_tst %/% 2L),
        info_str = paste(mod_nm, fam_nm, sep = "__"),
        fam_orig = get(paste0("f_", fam_nm))
      )
    }
  }
  if (run_cvvs_kfold) {
    for (mod_nm in names(refmods_kfold)) {
      for (fam_nm in names(refmods_kfold[[mod_nm]])) {
        refmodel_tester(
          refmod = refmods_kfold[[mod_nm]][[fam_nm]],
          fit_expected = fits_kfold[[mod_nm]][[fam_nm]],
          nobsv_expected = nobsv,
          nrefdraws_expected = chains_tst * (iter_tst %/% 2L),
          info_str = paste(mod_nm, fam_nm, "kfold", sep = "__"),
          fam_orig = get(paste0("f_", fam_nm))
        )
      }
    }
  }
})

test_that("error if `data` is missing", {
  SW(fit_nodata <- rstanarm::stan_glm(
    dat$y_glm_gauss ~ dat$xco.1 + dat$xco.2 + dat$xco.3 + dat$xca.1 + dat$xca.2,
    family = f_gauss,
    weights = dat$wobs_tst, offset = dat$offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  ))
  expect_error(get_refmodel(fit_nodata),
               "^object of type 'environment' is not subsettable$")
})

test_that("error if `formula` is a character string", {
  SW(fit_str <- rstanarm::stan_glm(
    "y_glm_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2",
    family = f_gauss, data = dat,
    weights = wobs_tst, offset = offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  ))
  expect_error(get_refmodel(fit_str),
               "^inherits\\(formula, \"formula\"\\) is not TRUE$")
})

test_that("get_refmodel() is idempotent", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      expect_identical(get_refmodel(refmods[[mod_nm]][[fam_nm]]),
                       refmods[[mod_nm]][[fam_nm]],
                       info = paste(mod_nm, fam_nm, sep = "__"))
    }
  }
})

# predict.refmodel() ------------------------------------------------------

context("predict.refmodel()")

test_that("error if `type` is invalid", {
  expect_error(predict(refmods$glm$gauss, dat, type = "zzz"),
               "^type should be one of")
})

test_that("error if `ynew` is invalid", {
  expect_error(predict(refmods$glm$gauss, dat, ynew = dat),
               "^ynew must be a numerical vector$")
})

test_that("`object` of class `\"refmodel\"`, `newdata`, and `type` work", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      tstsetup <- paste(mod_nm, fam_nm, sep = "__")
      # We expect a warning which in fact should be suppressed, though (see
      # issue #162):
      warn_expected <- switch(
        mod_nm,
        "glm" = paste("^'offset' argument is NULL but it looks like you",
                      "estimated the model using an offset term\\.$"),
        NA
      )
      expect_warning(
        predref_resp <- predict(refmods[[mod_nm]][[fam_nm]], dat,
                                type = "response"),
        warn_expected,
        info = tstsetup
      )
      expect_warning(
        predref_link <- predict(refmods[[mod_nm]][[fam_nm]], dat,
                                type = "link"),
        warn_expected,
        info = tstsetup
      )
      expect_true(is.vector(predref_resp, "double"), info = tstsetup)
      expect_length(predref_resp, nobsv)
      if (fam_nm == "binom") {
        expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                    info = tstsetup)
      }
      expect_true(is.vector(predref_link, "double"), info = tstsetup)
      expect_length(predref_link, nobsv)
      if (fam_nm == "gauss") {
        expect_equal(predref_resp, predref_link, info = tstsetup)
      }
    }
  }
  if (run_cvvs_kfold) {
    for (mod_nm in names(refmods_kfold)) {
      for (fam_nm in names(refmods_kfold[[mod_nm]])) {
        tstsetup <- paste(mod_nm, fam_nm, "kfold", sep = "__")
        # We expect a warning which in fact should be suppressed, though (see
        # issue #162):
        warn_expected <- switch(
          mod_nm,
          "glm" = paste("^'offset' argument is NULL but it looks like you",
                        "estimated the model using an offset term\\.$"),
          NA
        )
        expect_warning(
          predref_resp <- predict(refmods_kfold[[mod_nm]][[fam_nm]], dat,
                                  type = "response"),
          warn_expected,
          info = tstsetup
        )
        expect_warning(
          predref_link <- predict(refmods_kfold[[mod_nm]][[fam_nm]], dat,
                                  type = "link"),
          warn_expected,
          info = tstsetup
        )
        expect_true(is.vector(predref_resp, "double"), info = tstsetup)
        expect_length(predref_resp, nobsv)
        if (fam_nm == "binom") {
          expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                      info = tstsetup)
        }
        expect_true(is.vector(predref_link, "double"), info = tstsetup)
        expect_length(predref_link, nobsv)
        if (fam_nm == "gauss") {
          expect_equal(predref_resp, predref_link, info = tstsetup)
        }
      }
    }
  }
})

test_that("`ynew` works", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      tstsetup <- paste(mod_nm, fam_nm, sep = "__")
      y_crr <- dat[, paste("y", mod_nm, fam_nm, sep = "_")]
      if (is.integer(y_crr)) {
        y_crr <- as.numeric(y_crr)
      }
      # We expect a warning which in fact should be suppressed, though (see
      # issue #162):
      warn_expected <- switch(
        mod_nm,
        "glm" = paste("^'offset' argument is NULL but it looks like you",
                      "estimated the model using an offset term\\.$"),
        NA
      )
      expect_warning(
        predref_resp <- predict(refmods[[mod_nm]][[fam_nm]], dat, ynew = y_crr,
                                type = "response"),
        warn_expected,
        info = tstsetup
      )
      expect_warning(
        predref_link <- predict(refmods[[mod_nm]][[fam_nm]], dat, ynew = y_crr,
                                type = "link"),
        warn_expected,
        info = tstsetup
      )
      expect_equal(predref_resp, predref_link, info = tstsetup)
      expect_true(is.vector(predref_resp, "double"), info = tstsetup)
      expect_length(predref_resp, nobsv)
    }
  }
  if (run_cvvs_kfold) {
    for (mod_nm in names(refmods_kfold)) {
      for (fam_nm in names(refmods_kfold[[mod_nm]])) {
        tstsetup <- paste(mod_nm, fam_nm, "kfold", sep = "__")
        y_crr <- dat[, paste("y", mod_nm, fam_nm, sep = "_")]
        if (is.integer(y_crr)) {
          y_crr <- as.numeric(y_crr)
        }
        # We expect a warning which in fact should be suppressed, though (see
        # issue #162):
        warn_expected <- switch(
          mod_nm,
          "glm" = paste("^'offset' argument is NULL but it looks like you",
                        "estimated the model using an offset term\\.$"),
          NA
        )
        expect_warning(
          predref_resp <- predict(refmods_kfold[[mod_nm]][[fam_nm]], dat,
                                  ynew = y_crr, type = "response"),
          warn_expected,
          info = tstsetup
        )
        expect_warning(
          predref_link <- predict(refmods_kfold[[mod_nm]][[fam_nm]], dat,
                                  ynew = y_crr, type = "link"),
          warn_expected,
          info = tstsetup
        )
        expect_equal(predref_resp, predref_link, info = tstsetup)
        expect_true(is.vector(predref_resp, "double"), info = tstsetup)
        expect_length(predref_resp, nobsv)
      }
    }
  }
})
