test_that("extend_family() works for the traditional projection", {
  fam_nms_tst <- c("gauss", "binom", "poiss")
  ### TODO: Also test other families extensible by projpred:
  # fam_nms_tst <- c(fam_nms_tst, "Gamma", "Student_t")
  ###
  fam_nms_tst <- setNames(nm = fam_nms_tst)
  for (fam_nm in fam_nms_tst) {
    fam_orig_crr <- get(paste0("f_", fam_nm))
    extfam_crr <- extend_family(fam_orig_crr)
    extfam_tester(extfam = extfam_crr,
                  fam_orig = fam_orig_crr,
                  info_str = fam_nm)
  }
})

test_that("extend_family() works for the augmented-data projection", {
  fam_nms_tst <- c("binom_augdat", "cumul_augdat_rstanarm")
  # For families brms::categorical(), brms::cumulative(), brms::cratio(),
  # brms::sratio(), and brms::acat(), extend_family() can't be tested easily
  # without receiving such a family from a `brmsfit`, so skip them here (they
  # are tested by refmodel_tester() in `test_refmodel.R` anyway).
  fam_nms_tst <- setNames(nm = fam_nms_tst)
  for (fam_nm in fam_nms_tst) {
    if (fam_nm == "binom_augdat") {
      fam_orig_crr <- binomial()
      extfam_crr <- extend_family(fam_orig_crr,
                                  augdat_y_unqs = c("0", "1"),
                                  augdat_link = augdat_link_binom,
                                  augdat_ilink = augdat_ilink_binom)
    } else if (fam_nm == "cumul_augdat_rstanarm") {
      fam_orig_crr <- get_f_cumul()
      extfam_crr <- extend_family(fam_orig_crr,
                                  augdat_y_unqs = paste0("resplvl", 1:3),
                                  augdat_link = augdat_link_cumul,
                                  augdat_ilink = augdat_ilink_cumul,
                                  augdat_args_link = list(link = "logit"),
                                  augdat_args_ilink = list(link = "logit"))
    }
    extfam_tester(extfam = extfam_crr,
                  fam_orig = fam_orig_crr,
                  augdat_expected = TRUE,
                  info_str = fam_nm)
  }
})

test_that("extend_family() works for the latent projection", {
  fam_nms_tst <- c("binom_latent", "cumul_latent_rstanarm")
  # For family brms::cumulative(), extend_family() can't be tested easily
  # without receiving such a family from a `brmsfit`, so skip it here (it is
  # tested by refmodel_tester() in `test_refmodel.R` anyway).
  fam_nms_tst <- setNames(nm = fam_nms_tst)
  ilink_dummy <- identity
  for (fam_nm in fam_nms_tst) {
    if (fam_nm == "binom_latent") {
      extfam_crr <- extend_family(
        binomial(),
        latent = TRUE,
        # latent_y_unqs = c("0", "1"),
        latent_ilink = ilink_dummy,
        latent_ll_oscale = latent_ll_oscale_binom_nocats,
        latent_ppd_oscale = latent_ppd_oscale_binom_nocats
      )
    } else if (fam_nm == "cumul_latent_rstanarm") {
      extfam_crr <- extend_family(get_f_cumul(),
                                  latent = TRUE,
                                  latent_y_unqs = paste0("resplvl", 1:3),
                                  latent_ilink = ilink_dummy,
                                  latent_ll_oscale = latent_ll_oscale_cats,
                                  latent_ppd_oscale = latent_ppd_oscale_cats)
    }
    extfam_tester(extfam = extfam_crr,
                  fam_orig = gaussian(),
                  latent_expected = TRUE,
                  info_str = fam_nm)
  }
})

test_that(paste(
  "the latent projection stops for brms families with several linear",
  "predictors"
), {
  skip_if_not_installed("brms")
  msg_categ <- paste("support the latent projection for the brms family",
                     "`categorical()`")
  expect_error(extend_family(brms::categorical(), latent = TRUE), msg_categ,
               fixed = TRUE)
  for (latent in c(FALSE, TRUE)) {
    expect_error(extend_family(brms::multinomial(), latent = latent),
                 "does not support the brms family `multinomial()`",
                 fixed = TRUE, info = paste("latent:", latent))
  }
  tstsetups <- grep("^brms\\.(glm|glmm)\\.categ\\.", names(fits), value = TRUE)
  for (tstsetup in tstsetups) {
    expect_error(get_refmodel(fits[[tstsetup]], latent = TRUE), msg_categ,
                 fixed = TRUE, info = tstsetup)
  }
})

test_that(paste(
  "extend_family() with `latent = TRUE` gives the `latent_ilink` message only",
  "for families whose `linkinv` does not give the mean"
), {
  fam_tst <- function(fam_nm, link, linkinv) {
    structure(list(family = fam_nm, link = link, linkinv = linkinv),
              class = "family")
  }
  fams_mean <- list(
    neg_binomial_2 = fam_tst("neg_binomial_2", "log", exp),
    negbinomial = fam_tst("negbinomial", "log", exp),
    Gamma = Gamma(link = "log"),
    weibull = fam_tst("weibull", "log", exp)
  )
  for (fam_nm in names(fams_mean)) {
    expect_no_message(
      extend_family(fams_mean[[fam_nm]], latent = TRUE),
      message = "latent_ilink"
    )
  }
  expect_message(
    extend_family(fam_tst("lognormal", "identity", identity), latent = TRUE),
    "latent_ilink",
    info = "lognormal"
  )
})
