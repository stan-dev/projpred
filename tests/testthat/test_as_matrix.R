context("as.matrix.projection")

# Gaussian and binomial reference models without multilevel or additive terms:
settings_list <- list(
  gauss = list(
    refmod = refmods_glm$gauss,
    solterms_list = list(character(), solterms_tst),
    ndraws_list = list(25L, 2L, 1L)
  ),
  binom = list(
    refmod = refmods_glm$binom,
    solterms_list = list(solterms_tst),
    ndraws_list = list(25L)
  )
)

for (fam_type in settings_list) {
  for (solterms_tsttmp in fam_type$solterms_list) {
    for (ndraws_tsttmp in fam_type$ndraws_list) {
      tstsetup <- unlist(nlist(fam_nm = fam_type$refmod$family$family,
                               solterms_tsttmp,
                               ndraws_tsttmp))
      prj <- project(fam_type$refmod,
                     solution_terms = solterms_tsttmp,
                     ndraws = ndraws_tsttmp)

      # Expected warning (more precisely: regexp which is matched against the
      # warning; NA means no warning) for as.matrix.projection():
      if (ndraws_tsttmp > 20) {
        warn_prjmat_expect <- NA
      } else {
        # Clustered projection, so we expect a warning:
        warn_prjmat_expect <- "the clusters might have different weights"
      }

      expect_warning(m <- as.matrix(prj), warn_prjmat_expect, info = tstsetup)

      if (fam_type$refmod$family$family == "gaussian") {
        npars_fam <- "sigma"
      } else if (fam_type$refmod$family$family == "binomial") {
        npars_fam <- character()
      }
      test_that("as.matrix.projection()'s output structure is correct", {
        expect_identical(
          dim(m),
          c(ndraws_tsttmp, length(solterms_tsttmp) + 1L + length(npars_fam)),
          info = tstsetup
        )
        expect_identical(
          colnames(m),
          c(paste0("b_", c("Intercept", solterms_tsttmp)), npars_fam),
          info = tstsetup
        )
      })
    }
  }
}

# GLMMs -------------------------------------------------------------------

settings_list <- list(
  gauss = list(
    refmod = refmods_glmm$gauss,
    solterms_list = list(character(),
                         solterms_tst,
                         c("x.3", "(1 | x.gr)", "x.1 + (x.1 | x.gr)")),
    ndraws_list = list(25L, 2L, 1L)
  ),
  binom = list(
    refmod = refmods_glmm$binom,
    solterms_list = list(solterms_tst),
    ndraws_list = list(25L)
  )
)

for (fam_type in settings_list) {
  for (solterms_tsttmp in fam_type$solterms_list) {
    for (ndraws_tsttmp in fam_type$ndraws_list) {
      tstsetup <- unlist(nlist(fam_nm = fam_type$refmod$family$family,
                               solterms_tsttmp,
                               ndraws_tsttmp))
      prj <- project(fam_type$refmod,
                     solution_terms = solterms_tsttmp,
                     ndraws = ndraws_tsttmp)

      # Expected warning (more precisely: regexp which is matched against the
      # warning; NA means no warning) for as.matrix.projection():
      if (ndraws_tsttmp > 20) {
        warn_prjmat_expect <- NA
      } else {
        # Clustered projection, so we expect a warning:
        warn_prjmat_expect <- "the clusters might have different weights"
      }

      expect_warning(m <- as.matrix(prj), warn_prjmat_expect, info = tstsetup)

      if (fam_type$refmod$family$family == "gaussian") {
        npars_fam <- "sigma"
      } else if (fam_type$refmod$family$family == "binomial") {
        npars_fam <- character()
      }
      test_that("as.matrix.projection()'s output structure is correct", {
        colnms_prjmat_expect <- c(
          "Intercept",
          grep("^x\\.[[:digit:]]$", solterms_tsttmp,
               value = TRUE)
        )
        if ("x.1 + (x.1 | x.gr)" %in% solterms_tsttmp) {
          colnms_prjmat_expect <- c(colnms_prjmat_expect, "x.1")
        }
        colnms_prjmat_expect <- paste0("b_", colnms_prjmat_expect)
        if ("(1 | x.gr)" %in% solterms_tsttmp) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "sd_x.gr__Intercept"
          )
        }
        if ("x.1 + (x.1 | x.gr)" %in% solterms_tsttmp) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "sd_x.gr__x.1"
          )
        }
        if (all(c("(1 | x.gr)", "x.1 + (x.1 | x.gr)") %in% solterms_tsttmp)) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "cor_x.gr__Intercept__x.1"
          )
        }
        if ("(1 | x.gr)" %in% solterms_tsttmp) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            paste0("r_x.gr[gr", seq_len(ngr), ",Intercept]")
          )
        }
        if ("x.1 + (x.1 | x.gr)" %in% solterms_tsttmp) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            paste0("r_x.gr[gr", seq_len(ngr), ",x.1]")
          )
        }
        colnms_prjmat_expect <- c(colnms_prjmat_expect, npars_fam)

        expect_identical(dim(m), c(ndraws_tsttmp, length(colnms_prjmat_expect)),
                         info = tstsetup)
        expect_identical(colnames(m), colnms_prjmat_expect, info = tstsetup)
      })
    }
  }
}
