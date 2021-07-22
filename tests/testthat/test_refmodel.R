# get_refmodel() ----------------------------------------------------------

context("get_refmodel()")

test_that("`object` of class \"stanreg\" works", {
  refmod_nms <- c(
    "fit", "formula", "div_minimizer", "family", "mu", "dis", "y", "loglik",
    "intercept", "proj_predfun", "fetch_data", "wobs", "wsample", "offset",
    "folds", "cvfun", "cvfits", "extract_model_data", "ref_predfun"
  )
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      refmod <- refmods[[mod_nm]][[fam_nm]]
      fit_expected <- fits[[mod_nm]][[fam_nm]]
      fam_fullnm_expected <- get(paste0("f_", fam_nm))$family
      info_str <- paste(mod_nm, fam_nm, sep = "__")

      expect_s3_class(refmod, "refmodel", exact = TRUE)
      expect_type(refmod, "list")
      expect_named(refmod, refmod_nms, info = info_str)

      expect_identical(refmod$fit, fit_expected, info = info_str)
      if (fam_fullnm_expected == "binomial") {
        formul_expected_chr <- as.character(fit_expected$formula)
        stopifnot(length(formul_expected_chr) == 3)
        y_expected_chr <- sub("^cbind\\(", "", formul_expected_chr[2])
        y_expected_chr <- sub(",.*\\)$", "", y_expected_chr)
        formul_expected <- as.formula(paste(
          y_expected_chr,
          formul_expected_chr[1],
          formul_expected_chr[3]
        ))
        expect_identical(refmod$formula, formul_expected, info = info_str)
      } else {
        expect_identical(refmod$formula, fit_expected$formula, info = info_str)
      }
      expect_type(refmod$div_minimizer, "closure")
      expect_s3_class(refmod$family, "family")
      expect_identical(refmod$family$family, fam_fullnm_expected,
                       info = info_str)
      expect_identical(refmod$family$family, fit_expected$family$family,
                       info = info_str)
      # extfam_tester(refmod$family)
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

test_that("predict checks the 'type' argument", {
  expect_error(
    predict(ref_gauss, df_gauss, type = "zzz"),
    "type should be one of"
  )
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
