context("extend_family()")

test_that("`family` works", {
  fam_nms_tst <- c("gauss", "binom", "poiss", "binom_augdat",
                   "cumul_augdat_rstanarm")
  fam_nms_brms_augdat <- c(
    "categ_augdat", "cumul_augdat", "cratio_augdat", "sratio_augdat",
    "acat_augdat"
  )
  if (run_brms && packageVersion("brms") >= "2.16.6") {
    fam_nms_tst <- c(fam_nms_tst, fam_nms_brms_augdat)
  }
  ### TODO: Also test other families extensible by projpred:
  # fam_nms_tst <- c(fam_nms_tst, "Gamma", "Student_t")
  ###
  fam_nms_tst <- setNames(nm = fam_nms_tst)
  for (fam_nm in fam_nms_tst) {
    if (fam_nm == "binom_augdat") {
      fam_orig_crr <- binomial()
      extfam_crr <- extend_family(fam_orig_crr,
                                  augdat_y_unqs = c("0", "1"),
                                  augdat_link = augdat_link_binom,
                                  augdat_ilink = augdat_ilink_binom)
      augdat_expected_crr <- TRUE
    } else if (fam_nm == "cumul_augdat_rstanarm") {
      fam_orig_crr <- structure(list(family = "cumulative_rstanarm",
                                     link = "logit"),
                                class = "family")
      extfam_crr <- extend_family(fam_orig_crr,
                                  augdat_y_unqs = paste0("resplvl", 1:3),
                                  augdat_link = augdat_link_cumul,
                                  augdat_ilink = augdat_ilink_cumul,
                                  augdat_args_link = list(link = "logit"),
                                  augdat_args_ilink = list(link = "logit"))
      augdat_expected_crr <- TRUE
    } else if (fam_nm %in% fam_nms_brms_augdat) {
      fam_orig_crr <- switch(fam_nm,
                             categ_augdat = brms::categorical(),
                             cumul_augdat = brms::cumulative(),
                             cratio_augdat = brms::cratio(),
                             sratio_augdat = brms::sratio(),
                             acat_augdat = brms::acat())
      # These families are most easily tested when part of a `refmodel` created
      # by `brms:::get_refmodel.brmsfit()`:
      bfit_dummy <- suppressWarnings(brms::brm(
        formula = y_dummy ~ 1,
        data = data.frame(y_dummy = rep(1:3, 2)),
        family = fam_orig_crr,
        chains = chains_tst,
        iter = iter_tst,
        file = file.path(file_pth, paste0("bfit_dummy_", fam_nm)),
        file_refit = "on_change",
        seed = seed_fit,
        refresh = 0
      ))
      extfam_crr <- get_refmodel(bfit_dummy, brms_seed = seed2_tst)$family
      augdat_expected_crr <- TRUE
    } else {
      fam_orig_crr <- get(paste0("f_", fam_nm))
      extfam_crr <- extend_family(fam_orig_crr)
      augdat_expected_crr <- FALSE
    }
    if (fam_nm %in% fam_nms_brms_augdat) {
      extfam_nms_add2_crr <- "mu_fun"
    } else {
      extfam_nms_add2_crr <- character()
    }
    from_brms_crr <- fam_nm %in% fam_nms_brms_augdat
    extfam_tester(extfam = extfam_crr,
                  fam_orig = fam_orig_crr,
                  extfam_nms_add2 = extfam_nms_add2_crr,
                  from_brms = from_brms_crr,
                  augdat_expected = augdat_expected_crr,
                  info_str = fam_nm)
  }
})
