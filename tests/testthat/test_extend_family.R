context("extend_family()")

test_that("`family` works", {
  for (fam_nm in fam_nms) {
    fam_orig_crr <- get(paste0("f_", fam_nm))
    extfam_tester(
      extend_family(fam_orig_crr),
      fam_orig = fam_orig_crr,
      info_str = fam_nm
    )
  }
})
