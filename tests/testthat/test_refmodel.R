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

# TODO

test_that("error if `type` is invalid", {
  expect_error(predict(refmods$glm$gauss, dat, type = "zzz"),
               "^type should be one of")
})

test_that("predict produces sensible results for gaussian models", {
  out.resp <- predict(ref_gauss, df_gauss, type = "response")
  expect_vector(out.resp)
  expect_length(out.resp, nrow(df_gauss))

  out.link <- predict(ref_gauss, df_gauss, type = "link")
  expect_equal(out.resp, out.link)
})

test_that("predict produces sensible results for binomial models", {
  out.resp <- predict(ref_binom, df_binom, type = "response")
  expect_vector(out.resp)
  expect_length(out.resp, nrow(df_binom))
  expect_true(all(out.resp >= 0 & out.resp <= 1))

  out.link <- predict(ref_binom, df_binom, type = "link")
  expect_length(out.resp, nrow(df_binom))
})

test_that("predict produces sensible results when specifying ynew", {
  out <- predict(ref_gauss, df_gauss, ynew = df_gauss$y)
  expect_vector(out)
  expect_length(out, length(df_gauss$y))

  expect_error(
    predict(ref_gauss, df_gauss, ynew = df_gauss),
    "must be a numerical vector"
  )
})
