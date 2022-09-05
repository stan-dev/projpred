context("extend_family()")

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
      fam_orig_crr <- structure(list(family = "cumulative_rstanarm",
                                     link = "logit"),
                                class = "family")
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
