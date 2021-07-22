context("extend_family()")

test_that("`family` works", {
  fam_nms_tst <- fam_nms
  if (!"poiss" %in% fam_nms) {
    fam_nms_tst <- setNames(nm = c(fam_nms_tst, "poiss"))
  }
  ### TODO: Also test other families extensible by projpred:
  # fam_nms_tst <- setNames(nm = c(fam_nms_tst, "Gamma", "Student_t"))
  ###
  for (fam_nm in fam_nms_tst) {
    fam_orig_crr <- get(paste0("f_", fam_nm))
    extfam_tester(
      extend_family(fam_orig_crr),
      fam_orig = fam_orig_crr,
      info_str = fam_nm
    )
  }
})
