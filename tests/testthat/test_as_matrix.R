context("as.matrix.projection")

# GLMs --------------------------------------------------------------------

settings_list_glm <- list(
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

for (settings_obj in settings_list_glm) {
  for (solterms_crr in settings_obj$solterms_list) {
    for (ndraws_crr in settings_obj$ndraws_list) {
      tstsetup <- unlist(nlist(fam_nm = settings_obj$refmod$family$family,
                               solterms_crr,
                               ndraws_crr))
      prj <- project(settings_obj$refmod,
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

      if (settings_obj$refmod$family$family == "gaussian") {
        npars_fam <- "sigma"
      } else if (settings_obj$refmod$family$family == "binomial") {
        npars_fam <- character()
      }
      test_that("as.matrix.projection()'s output structure is correct", {
        expect_identical(
          dim(m),
          c(ndraws_crr, length(solterms_crr) + 1L + length(npars_fam)),
          info = tstsetup
        )
        expect_identical(
          colnames(m),
          c(paste0("b_", c("Intercept", solterms_crr)), npars_fam),
          info = tstsetup
        )
      })
    }
  }
}

# GLMMs -------------------------------------------------------------------

settings_list_glmm <- list(
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

for (settings_obj in settings_list_glmm) {
  for (solterms_crr in settings_obj$solterms_list) {
    for (ndraws_crr in settings_obj$ndraws_list) {
      tstsetup <- unlist(nlist(fam_nm = settings_obj$refmod$family$family,
                               solterms_crr,
                               ndraws_crr))
      prj <- project(settings_obj$refmod,
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

      if (settings_obj$refmod$family$family == "gaussian") {
        npars_fam <- "sigma"
      } else if (settings_obj$refmod$family$family == "binomial") {
        npars_fam <- character()
      }
      test_that("as.matrix.projection()'s output structure is correct", {
        colnms_prjmat_expect <- c(
          "Intercept",
          grep("^x\\.[[:digit:]]$", solterms_crr,
               value = TRUE)
        )
        if ("x.1 + (x.1 | x.gr)" %in% solterms_crr) {
          colnms_prjmat_expect <- c(colnms_prjmat_expect, "x.1")
        }
        colnms_prjmat_expect <- paste0("b_", colnms_prjmat_expect)
        if ("(1 | x.gr)" %in% solterms_crr) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "sd_x.gr__Intercept"
          )
        }
        if ("x.1 + (x.1 | x.gr)" %in% solterms_crr) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "sd_x.gr__x.1"
          )
        }
        if (all(c("(1 | x.gr)", "x.1 + (x.1 | x.gr)") %in% solterms_crr)) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            "cor_x.gr__Intercept__x.1"
          )
        }
        if ("(1 | x.gr)" %in% solterms_crr) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            paste0("r_x.gr[gr", seq_len(ngr), ",Intercept]")
          )
        }
        if ("x.1 + (x.1 | x.gr)" %in% solterms_crr) {
          colnms_prjmat_expect <- c(
            colnms_prjmat_expect,
            paste0("r_x.gr[gr", seq_len(ngr), ",x.1]")
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

# GAMs --------------------------------------------------------------------

settings_list_gam <- list(
  gauss = list(
    refmod = refmods_gam$gauss,
    solterms_list = list(character(),
                         solterms_tst_gam,
                         c(solterms_tst_gam, "s(x.2)")),
    ndraws_list = list(25L, 2L, 1L)
  ),
  binom = list(
    refmod = refmods_gam$binom,
    solterms_list = list(solterms_tst_gam),
    ndraws_list = list(25L)
  )
)

for (settings_obj in settings_list_gam) {
  for (solterms_crr in settings_obj$solterms_list) {
    for (ndraws_crr in settings_obj$ndraws_list) {
      tstsetup <- unlist(nlist(fam_nm = settings_obj$refmod$family$family,
                               solterms_crr,
                               ndraws_crr))
      prj <- project(settings_obj$refmod,
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

      if (settings_obj$refmod$family$family == "gaussian") {
        npars_fam <- "sigma"
      } else if (settings_obj$refmod$family$family == "binomial") {
        npars_fam <- character()
      }

      par_nms_orig <- colnames(as.matrix(settings_obj$refmod$fit))
      test_that("as.matrix.projection()'s output structure is correct", {
        # TODO
      })
    }
  }
}

# GAMMs -------------------------------------------------------------------

# Currently deactivated (see `setup.R`).
