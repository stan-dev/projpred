context("as.matrix.projection")

# Gaussian and binomial reference models without multilevel or additive terms:
settings_list <- list(
  gauss = list(
    refmod = refmods_glm$gauss,
    solution_terms_list = list(character(), solterms_tst),
    ndraws_list = list(25L, 2L, 1L)
  ),
  binom = list(
    refmod = refmods_glm$binom,
    solution_terms_list = list(solterms_tst),
    ndraws_list = list(25L)
  )
)

for (fam_type in settings_list) {
  for (solution_terms in fam_type$solution_terms_list) {
    for (ndraws in fam_type$ndraws_list) {
      prj <- project(fam_type$refmod,
                     solution_terms = solution_terms,
                     ndraws = ndraws)

      # Expected warning (more precisely: regexp which is matched against the
      # warning; NA means no warning) for as.matrix.projection():
      if (ndraws > 20) {
        warn_prjmat_expect <- NA
      } else {
        # Clustered projection, so we expect a warning:
        warn_prjmat_expect <- "the clusters might have different weights"
      }

      expect_warning(m <- as.matrix(prj), warn_prjmat_expect)

      if (fam_type$refmod$family$family == "gaussian") {
        npars_fam <- "sigma"
      } else if (fam_type$refmod$family$family == "binomial") {
        npars_fam <- character()
      }
      test_that("as.matrix.projection()'s output structure is correct", {
        expect_equal(
          dim(m),
          c(ndraws, length(solution_terms) + 1 + length(npars_fam))
        )
        expect_identical(
          colnames(m),
          c(paste0("b_", c("Intercept", solution_terms)), npars_fam)
        )
      })
    }
  }
}

# GLMMs -------------------------------------------------------------------

settings_list <- list(
  gauss = list(
    refmod = refmods_glmm$gauss,
    solution_terms_list = list(character(),
                               solterms_tst,
                               c("x.3", "(1 | xgr)", "x.1 + (x.1 | xgr)")),
    ndraws_list = list(25L, 2L, 1L)
  ),
  binom = list(
    refmod = refmods_glmm$binom,
    solution_terms_list = list(solterms_tst),
    ndraws_list = list(25L)
  )
)

for (fam_type in settings_list) {
  for (solution_terms in fam_type$solution_terms_list) {
    for (ndraws in fam_type$ndraws_list) {
      prj <- project(fam_type$refmod,
                     solution_terms = solution_terms,
                     ndraws = ndraws)

      # Expected warning (more precisely: regexp which is matched against the
      # warning; NA means no warning) for as.matrix.projection():
      if (ndraws > 20) {
        warn_prjmat_expect <- NA
      } else {
        # Clustered projection, so we expect a warning:
        warn_prjmat_expect <- "the clusters might have different weights"
      }

      expect_warning(m <- as.matrix(prj), warn_prjmat_expect)

      if (fam_type$refmod$family$family == "gaussian") {
        npars_fam <- "sigma"
      } else if (fam_type$refmod$family$family == "binomial") {
        npars_fam <- character()
      }
      test_that("as.matrix.projection()'s output structure is correct", {
        colnms_prjmat_expect <- c(
          "Intercept",
          grep("^x\\.[[:digit:]]$", solution_terms,
               value = TRUE)
        )
        if ("x.1 + (x.1 | xgr)" %in% solution_terms) {
          colnms_prjmat_expect <- c(colnms_prjmat_expect, "x.1")
        }
        colnms_prjmat_expect <- paste0("b_", colnms_prjmat_expect)
        if ("(1 | xgr)" %in% solution_terms) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "sd_xgr__Intercept"
          )
        }
        if ("x.1 + (x.1 | xgr)" %in% solution_terms) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "sd_xgr__x.1"
          )
        }
        if (all(c("(1 | xgr)", "x.1 + (x.1 | xgr)") %in% solution_terms)) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "cor_xgr__Intercept__x.1"
          )
        }
        if ("(1 | xgr)" %in% solution_terms) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            paste0("r_xgr[gr", seq_len(ngr), ",Intercept]")
          )
        }
        if ("x.1 + (x.1 | xgr)" %in% solution_terms) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            paste0("r_xgr[gr", seq_len(ngr), ",x.1]")
          )
        }
        colnms_prjmat_expect <- c(colnms_prjmat_expect, npars_fam)

        expect_equal(dim(m), c(ndraws, length(colnms_prjmat_expect)))
        expect_identical(colnames(m), colnms_prjmat_expect)
      })
    }
  }
}
