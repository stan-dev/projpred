# proj_epred() -------------------------------------------------------------

context("proj_epred()")

## object ------------------------------------------------------------------

test_that("pe: `object` of class `projection` works", {
  skip_if_not(run_prj)
  for (tstsetup in names(prjs)) {
    if (!is.null(refmods[[args_prj[[tstsetup]]$tstsetup_ref]]$family$cats)) {
      ncats_nlats_expected_crr <- length(
        refmods[[args_prj[[tstsetup]]$tstsetup_ref]]$family$cats
      )
    } else {
      ncats_nlats_expected_crr <- integer()
    }
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    if (!has_const_wdr_prj(prjs[[tstsetup]])) {
      wdr_crr <- prjs[[tstsetup]][["wdraws_prj"]]
    } else {
      wdr_crr <- NULL
    }
    pe_tester(pes[[tstsetup]],
              nprjdraws_expected = ndr_ncl$nprjdraws,
              wdraws_prj_expected = wdr_crr,
              ncats_nlats_expected = list(ncats_nlats_expected_crr),
              info_str = tstsetup)
  }
})

test_that(paste(
  "pe: `object` of (informal) class `proj_list` (based on varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs)) {
    tstsetup_vs <- args_prj_vs[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- suggest_size(vss[[tstsetup_vs]], warnings = FALSE)
    }
    if (!is.null(
      refmods[[args_prj_vs[[tstsetup]]$tstsetup_ref]]$family$cats
    )) {
      ncats_nlats_expected_crr <- length(
        refmods[[args_prj_vs[[tstsetup]]$tstsetup_ref]]$family$cats
      )
    } else {
      ncats_nlats_expected_crr <- integer()
    }
    ndr_ncl <- ndr_ncl_dtls(args_prj_vs[[tstsetup]])
    if (!has_const_wdr_prj(prjs_vs[[tstsetup]])) {
      if (length(nterms_crr) > 1) {
        wdr_crr <- drop(unique(do.call(rbind, lapply(prjs_vs[[tstsetup]], "[[",
                                                     "wdraws_prj"))))
      } else {
        wdr_crr <- prjs_vs[[tstsetup]][["wdraws_prj"]]
      }
    } else {
      wdr_crr <- NULL
    }
    pe_tester(pes_vs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_expected = ndr_ncl$nprjdraws,
              wdraws_prj_expected = wdr_crr,
              ncats_nlats_expected = replicate(length(nterms_crr),
                                               ncats_nlats_expected_crr,
                                               simplify = FALSE),
              info_str = tstsetup)
  }
})

test_that(paste(
  "pe: `object` of (informal) class `proj_list` (based on cv_varsel()) works"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(prjs_cvvs)) {
    tstsetup_cvvs <- args_prj_cvvs[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- suggest_size(cvvss[[tstsetup_cvvs]], warnings = FALSE)
    }
    if (!is.null(
      refmods[[args_prj_cvvs[[tstsetup]]$tstsetup_ref]]$family$cats
    )) {
      ncats_nlats_expected_crr <- length(
        refmods[[args_prj_cvvs[[tstsetup]]$tstsetup_ref]]$family$cats
      )
    } else {
      ncats_nlats_expected_crr <- integer()
    }
    ndr_ncl <- ndr_ncl_dtls(args_prj_cvvs[[tstsetup]])
    if (!has_const_wdr_prj(prjs_cvvs[[tstsetup]])) {
      if (length(nterms_crr) > 1) {
        wdr_crr <- drop(unique(do.call(rbind, lapply(prjs_cvvs[[tstsetup]],
                                                     "[[", "wdraws_prj"))))
      } else {
        wdr_crr <- prjs_cvvs[[tstsetup]][["wdraws_prj"]]
      }
    } else {
      wdr_crr <- NULL
    }
    pe_tester(pes_cvvs[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_expected = ndr_ncl$nprjdraws,
              wdraws_prj_expected = wdr_crr,
              ncats_nlats_expected = replicate(length(nterms_crr),
                                               ncats_nlats_expected_crr,
                                               simplify = FALSE),
              info_str = tstsetup)
  }
})

## Equivalence with proj_linpred(transform = TRUE) -------------------------

test_that("pe: equivalent to proj_linpred(transform = TRUE)$pred", {
  skip_if_not(run_prj)
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    pl_true <- proj_linpred(prjs[[tstsetup]], transform = TRUE,
                            allow_nonconst_wdraws_prj = ndr_ncl$clust_used_gt1,
                            .seed = seed2_tst)
    pe_crr <- proj_epred(prjs[[tstsetup]],
                         allow_nonconst_wdraws_prj = ndr_ncl$clust_used_gt1,
                         .seed = seed2_tst)
    expect_equal(pe_crr, pl_true$pred, info = tstsetup)
  }
})

## newdata -----------------------------------------------------------------

test_that("pe: `newdata` works", {
  skip_if_not(run_prj)
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    if (!has_const_wdr_prj(prjs[[tstsetup]])) {
      wdr_crr <- prjs[[tstsetup]][["wdraws_prj"]]
    } else {
      wdr_crr <- NULL
    }
    dat_crr <- get_dat(tstsetup)
    for (nobsv_crr in nobsv_tst) {
      if (args_prj[[tstsetup]]$mod_nm == "gamm") {
        # TODO (GAMMs): Fix this.
        next
      }
      if (!is.null(
        refmods[[args_prj[[tstsetup]]$tstsetup_ref]]$family$cats
      )) {
        ncats_nlats_expected_crr <- length(
          refmods[[args_prj[[tstsetup]]$tstsetup_ref]]$family$cats
        )
      } else {
        ncats_nlats_expected_crr <- integer()
      }
      if (grepl("\\.with_wobs|\\.binom", tstsetup)) {
        wobs_crr <- head(prjs[[tstsetup]]$refmodel$wobs, nobsv_crr)
      } else {
        wobs_crr <- NULL
      }
      if (grepl("\\.with_offs", tstsetup)) {
        offs_crr <- head(prjs[[tstsetup]]$refmodel$offset, nobsv_crr)
      } else {
        offs_crr <- NULL
      }
      expect_warning(
        pe_crr <- proj_epred(
          prjs[[tstsetup]],
          newdata = head(dat_crr, nobsv_crr),
          weightsnew = wobs_crr,
          offsetnew = offs_crr,
          allow_nonconst_wdraws_prj = ndr_ncl$clust_used_gt1,
          .seed = seed2_tst
        ),
        get_warn_wrhs_orhs(tstsetup, weightsnew = wobs_crr,
                           offsetnew = offs_crr),
        info = tstsetup
      )
      pe_tester(pe_crr,
                nprjdraws_expected = ndr_ncl$nprjdraws,
                wdraws_prj_expected = wdr_crr,
                nobsv_expected = nobsv_crr,
                ncats_nlats_expected = list(ncats_nlats_expected_crr),
                info_str = paste(tstsetup, nobsv_crr, sep = "__"))
    }
  }
})

## integrated --------------------------------------------------------------

test_that("pe: `integrated` works", {
  skip_if_not(run_prj)
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    if (!is.null(
      refmods[[args_prj[[tstsetup]]$tstsetup_ref]]$family$cats
    )) {
      ncats_nlats_expected_crr <- length(
        refmods[[args_prj[[tstsetup]]$tstsetup_ref]]$family$cats
      )
    } else {
      ncats_nlats_expected_crr <- integer()
    }
    pe_intgr <- proj_epred(prjs[[tstsetup]], integrated = TRUE,
                           .seed = seed2_tst)
    pe_tester(pe_intgr,
              nprjdraws_expected = 1L,
              ncats_nlats_expected = list(ncats_nlats_expected_crr),
              info_str = paste(tstsetup, "integrated", sep = "__"))
  }
})

## return_draws_matrix -----------------------------------------------------

test_that(paste(
  "pe: `return_draws_matrix` causes a conversion of the output type"
), {
  skip_if_not(run_prj)
  skip_if_not_installed("posterior")
  for (tstsetup in names(prjs)) {
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    pe_orig <- pes[[tstsetup]]
    pe_dr <- proj_epred(prjs[[tstsetup]],
                        return_draws_matrix = TRUE,
                        .seed = seed2_tst)
    if (args_prj[[tstsetup]]$prj_nm == "augdat" ||
        (args_prj[[tstsetup]]$prj_nm == "latent" && !is.null(
          refmods[[args_prj[[tstsetup]]$tstsetup_ref]]$family$cats
        ))) {
      pe_orig_flat <- do.call(rbind, apply(pe_orig, 1, as.vector,
                                           simplify = FALSE))
    } else {
      pe_orig_flat <- pe_orig
    }
    pe_dr_repl <- posterior::as_draws_matrix(pe_orig_flat)
    if (!has_const_wdr_prj(prjs[[tstsetup]])) {
      pe_dr_repl <- posterior::weight_draws(
        pe_dr_repl, weights = prjs[[tstsetup]][["wdraws_prj"]]
      )
    }
    expect_equal(pe_dr, pe_dr_repl, info = tstsetup)
  }
})
