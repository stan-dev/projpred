context("as.matrix.projection")

### TODO:
mod_nms <- setdiff(mod_nms, c("gam", "gamm"))
###

settings <- list(
  glm = list(
    gauss = list(
      solterms_list = solterms_glm,
      ndraws_list = list(25L, 2L, 1L)
    ),
    binom = list(
      solterms_list = solterms_glm["somecomb_xco"],
      ndraws_list = list(25L)
    )
  ),
  glmm = list(
    gauss = list(
      solterms_list = solterms_glmm,
      ndraws_list = list(25L, 2L, 1L)
    ),
    binom = list(
      solterms_list = solterms_glmm["somecomb_z"],
      ndraws_list = list(25L)
    )
  )#,
  ### TODO:
  # gam = list(
  #   gauss = list(
  #     solterms_list = list(character(),
  #                          "s(s.1)",
  #                          c("s(s.1)", "s(s.2)")),
  #     ndraws_list = list(25L, 2L, 1L)
  #   ),
  #   binom = list(
  #     solterms_list = list("s(s.1)"),
  #     ndraws_list = list(25L)
  #   )
  # ),
  # gamm = list(
  #   gauss = list(
  #     solterms_list = list(character(),
  #                          "s(s.1)",
  #                          c("s(s.1)", "s(s.2)")),
  #     ndraws_list = list(25L, 2L, 1L)
  #   ),
  #   binom = list(
  #     solterms_list = list("s(s.1)"),
  #     ndraws_list = list(25L)
  #   )
  # )
  ###
)

for (mod_nm in mod_nms) {
  for (fam_nm in fam_nms) {
    refmod <- refmods[[mod_nm]][[fam_nm]]
    # par_nms_orig <- colnames(as.matrix(refmod$fit))
    settings_crr <- settings[[mod_nm]][[fam_nm]]
    for (solterms_crr in settings_crr$solterms_list) {
      for (ndraws_crr in settings_crr$ndraws_list) {
        tstsetup <- unlist(nlist(mod_nm, fam_nm, solterms_crr, ndraws_crr))
        prj <- project(refmod,
                       solution_terms = solterms_crr,
                       ndraws = ndraws_crr)

        # Expected warning (more precisely: regexp which is matched against the
        # warning; NA means no warning) for as.matrix.projection():
        if (ndraws_crr > 20) {
          warn_prjmat_expect <- NA
        } else {
          # Clustered projection, so we expect a warning:
          warn_prjmat_expect <- "the clusters might have different weights"
        }

        expect_warning(m <- as.matrix(prj), warn_prjmat_expect, info = tstsetup)

        if (fam_nm == "gauss") {
          npars_fam <- "sigma"
        } else if (fam_nm == "binom") {
          npars_fam <- character()
        }
        test_that("as.matrix.projection()'s output structure is correct", {
          colnms_prjmat_expect <- c(
            "Intercept",
            grep("^x(co|ca)\\.[[:digit:]]$", solterms_crr,
                 value = TRUE)
          )
          xca_idxs <- as.integer(
            sub("^xca\\.", "",
                grep("^xca\\.", colnms_prjmat_expect, value = TRUE))
          )
          for (xca_idx in xca_idxs) {
            colnms_prjmat_expect <- grep(paste0("^xca\\.", xca_idx, "$"),
                                         colnms_prjmat_expect,
                                         value = TRUE, invert = TRUE)
            colnms_prjmat_expect <- c(
              colnms_prjmat_expect,
              paste0("xca.", xca_idx, "lvl", seq_len(nlvl_fix[xca_idx])[-1])
            )
          }
          if ("xco.1 + (xco.1 | z.1)" %in% solterms_crr) {
            if (!"xco.1" %in% colnms_prjmat_expect) {
              colnms_prjmat_expect <- c(colnms_prjmat_expect, "xco.1")
            }
          }
          colnms_prjmat_expect <- paste0("b_", colnms_prjmat_expect)
          if ("(1 | z.1)" %in% solterms_crr) {
            colnms_prjmat_expect <- c(
              colnms_prjmat_expect,
              "sd_z.1__Intercept"
            )
          }
          if ("xco.1 + (xco.1 | z.1)" %in% solterms_crr) {
            colnms_prjmat_expect <- c(
              colnms_prjmat_expect,
              "sd_z.1__xco.1"
            )
          }
          if (all(c("(1 | z.1)", "xco.1 + (xco.1 | z.1)") %in% solterms_crr)) {
            colnms_prjmat_expect <- c(
              colnms_prjmat_expect,
              "cor_z.1__Intercept__xco.1"
            )
          }
          if ("(1 | z.1)" %in% solterms_crr) {
            colnms_prjmat_expect <- c(
              colnms_prjmat_expect,
              paste0("r_z.1[lvl", seq_len(nlvl_ran[1]), ",Intercept]")
            )
          }
          if ("xco.1 + (xco.1 | z.1)" %in% solterms_crr) {
            colnms_prjmat_expect <- c(
              colnms_prjmat_expect,
              paste0("r_z.1[lvl", seq_len(nlvl_ran[1]), ",xco.1]")
            )
          }
          colnms_prjmat_expect <- c(colnms_prjmat_expect, npars_fam)

          expect_identical(dim(m), c(ndraws_crr, length(colnms_prjmat_expect)),
                           info = tstsetup)
          expect_identical(colnames(m), colnms_prjmat_expect, info = tstsetup)
        })
      }
    }
  }
}
