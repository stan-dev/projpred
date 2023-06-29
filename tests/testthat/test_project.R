context("project()")

# object, solution_terms, ndraws / nclusters, nterms ----------------------

test_that(paste(
  "`object` of class \"refmodel\", `solution_terms`, and `ndraws` (or",
  "`nclusters`) work"
), {
  skip_if_not(run_prj)
  for (tstsetup in names(prjs)) {
    args_prj_i <- args_prj[[tstsetup]]
    ndr_ncl <- ndr_ncl_dtls(args_prj_i)
    projection_tester(prjs[[tstsetup]],
                      refmod_expected = refmods[[args_prj_i$tstsetup_ref]],
                      solterms_expected = args_prj_i$solution_terms,
                      nprjdraws_expected = ndr_ncl$nprjdraws,
                      with_clusters = ndr_ncl$ndr_ncl_nm == "nclusters",
                      const_wdraws_prj_expected = !ndr_ncl$clust_used,
                      info_str = tstsetup)
  }
})

test_that("invalid `solution_terms` warns or fails", {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.glm\\.gauss.*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]

    # Non-`vsel` object combined with `solution_terms = NULL`:
    expect_error(
      do.call(project, c(
        list(object = refmods[[args_prj_i$tstsetup_ref]],
             solution_terms = NULL),
        excl_nonargs(args_prj_i, nms_excl_add = "solution_terms")
      )),
      paste("^Please provide an `object` of class \"vsel\" or use argument",
            "`solution_terms`\\.$"),
      info = tstsetup
    )

    # Repeating solution terms:
    p_long <- do.call(project, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]],
           solution_terms = rep_len(args_prj_i$solution_terms,
                                    length.out = 1e4)),
      excl_nonargs(args_prj_i, nms_excl_add = "solution_terms")
    ))
    expect_identical(p_long, prjs[[tstsetup]])

    # Invalid type:
    for (solterms_crr in list(2, 1:3, list(solterms_x, solterms_x))) {
      tstsetup_crr <- paste(tstsetup, paste(solterms_crr, collapse = ","),
                            sep = "__")
      expect_error(
        do.call(project, c(
          list(object = refmods[[args_prj_i$tstsetup_ref]],
               solution_terms = solterms_crr),
          excl_nonargs(args_prj_i, nms_excl_add = "solution_terms")
        )),
        paste(
          "^is\\.null\\(solution_terms\\) \\|\\| is\\.vector\\(solution_terms,",
          "\"character\"\\) is not TRUE$"
        ),
        info = tstsetup_crr
      )
    }

    # Should be working, but result in a projection onto the intercept-only
    # submodel:
    for (solterms_crr in list("1",
                              c("some_dummy_string", "another_dummy_string"))) {
      tstsetup_crr <- paste(tstsetup, paste(solterms_crr, collapse = ","),
                            sep = "__")
      expect_warning(
        p <- do.call(project, c(
          list(object = refmods[[args_prj_i$tstsetup_ref]],
               solution_terms = solterms_crr),
          excl_nonargs(args_prj_i, nms_excl_add = "solution_terms")
        )),
        paste("The following element\\(s\\) of `solution_terms` could not be",
              "found in the table of possible solution terms"),
        info = tstsetup_crr
      )
      projection_tester(p,
                        refmod_expected = refmods[[args_prj_i$tstsetup_ref]],
                        solterms_expected = character(),
                        nprjdraws_expected = nclusters_pred_tst,
                        with_clusters = TRUE,
                        const_wdraws_prj_expected = FALSE,
                        info_str = tstsetup_crr)
    }
  }
})

test_that("`object` of class \"stanreg\" or \"brmsfit\" works", {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.brnll\\..*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    p_fit <- do.call(project, c(
      list(object = fits[[args_prj_i$tstsetup_fit]]),
      excl_nonargs(args_prj_i),
      excl_nonargs(args_ref[[args_prj_i$tstsetup_ref]])
    ))
    expect_identical(p_fit, prjs[[tstsetup]], ignore.environment = TRUE,
                     info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by varsel()), `nclusters`, and",
  "`nterms` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs)) {
    tstsetup_vs <- args_prj_vs[[tstsetup]]$tstsetup_vsel
    stopifnot(length(tstsetup_vs) > 0)
    mod_crr <- args_vs[[tstsetup_vs]]$mod_nm
    fam_crr <- args_vs[[tstsetup_vs]]$fam_nm
    nterms_crr <- args_prj_vs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- suggest_size(vss[[tstsetup_vs]], warnings = FALSE)
    }
    if (length(nterms_crr) == 1) {
      solterms_expected_crr <- vss[[tstsetup_vs]]$solution_terms[
        seq_len(nterms_crr)
      ]
      projection_tester(
        prjs_vs[[tstsetup]],
        refmod_expected = refmods[[args_prj_vs[[tstsetup]]$tstsetup_ref]],
        solterms_expected = solterms_expected_crr,
        nprjdraws_expected = args_prj_vs[[tstsetup]]$nclusters,
        with_clusters = TRUE,
        const_wdraws_prj_expected = FALSE,
        info_str = tstsetup
      )
      # Check that projecting from the "vsel" object onto a single submodel
      # gives the same output as projecting the reference model onto that
      # submodel directly:
      tstsetup_tries <- grep(
        paste0("^",
               gsub("\\.", "\\\\.",
                    sub("(with_offs|without_offs).*", "\\1", tstsetup)),
               ".*\\.clust$"),
        names(prjs),
        value = TRUE
      )
      match_prj <- sapply(tstsetup_tries, function(tstsetup_try) {
        setequal(solterms_expected_crr, prjs[[tstsetup_try]]$solution_terms) &&
          args_prj_vs[[tstsetup]]$prj_nm == args_prj[[tstsetup_try]]$prj_nm
      })
      tstsetup_match_prj <- tstsetup_tries[match_prj]
      if (length(tstsetup_match_prj) == 1) {
        if (identical(prjs_vs[[tstsetup]]$solution_terms,
                      prjs[[tstsetup_match_prj]]$solution_terms)) {
          expect_identical(prjs_vs[[tstsetup]], prjs[[tstsetup_match_prj]],
                           ignore.environment = TRUE, info = tstsetup)
        } else {
          expect_setequal(prjs_vs[[tstsetup]]$solution_terms,
                          prjs[[tstsetup_match_prj]]$solution_terms)
          expect_equal(
            lapply(seq_along(prjs_vs[[tstsetup]]$outdmin), function(s_idx) {
              outdmin_s <- prjs_vs[[tstsetup]]$outdmin[[s_idx]]
              outdmin_s$beta <- outdmin_s$beta[
                rownames(prjs[[tstsetup_match_prj]]$outdmin[[s_idx]]$beta), ,
                drop = FALSE
              ]
              outdmin_s$formula <- prjs[[tstsetup_match_prj]]$outdmin[[s_idx]][[
                "formula"
              ]]
              outdmin_s$x <- outdmin_s$x[
                , colnames(prjs[[tstsetup_match_prj]]$outdmin[[s_idx]]$x),
                drop = FALSE
              ]
              return(outdmin_s)
            }),
            prjs[[tstsetup_match_prj]]$outdmin,
            tolerance = 1e1 * .Machine$double.eps, info = tstsetup
          )
          prj_nms <- names(prjs_vs[[tstsetup]])
          expect_identical(prj_nms, names(prjs[[tstsetup_match_prj]]),
                           info = tstsetup)
          prj_el_excl <- !prj_nms %in% c("solution_terms", "outdmin")
          expect_equal(prjs_vs[[tstsetup]][prj_el_excl],
                       prjs[[tstsetup_match_prj]][prj_el_excl],
                       tolerance = .Machine$double.eps, info = tstsetup)
        }
      } else if (length(tstsetup_match_prj) > 1) {
        stop("Unexpected number of matches.")
      }
    } else {
      proj_list_tester(
        prjs_vs[[tstsetup]],
        len_expected = length(nterms_crr),
        is_seq = all(diff(nterms_crr) == 1),
        info_str = tstsetup,
        refmod_expected = refmods[[args_prj_vs[[tstsetup]]$tstsetup_ref]],
        nprjdraws_expected = args_prj_vs[[tstsetup]]$nclusters,
        with_clusters = TRUE,
        const_wdraws_prj_expected = FALSE,
        prjdraw_weights_expected = prjs_vs[[tstsetup]][[1]]$wdraws_prj
      )
    }
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by cv_varsel()), `nclusters`, and",
  "`nterms` work"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(prjs_cvvs)) {
    tstsetup_cvvs <- args_prj_cvvs[[tstsetup]]$tstsetup_vsel
    stopifnot(length(tstsetup_cvvs) > 0)
    mod_crr <- args_cvvs[[tstsetup_cvvs]]$mod_nm
    fam_crr <- args_cvvs[[tstsetup_cvvs]]$fam_nm
    nterms_crr <- args_prj_cvvs[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- suggest_size(cvvss[[tstsetup_cvvs]], warnings = FALSE)
    }
    if (length(nterms_crr) == 1) {
      solterms_expected_crr <- cvvss[[tstsetup_cvvs]]$solution_terms[
        seq_len(nterms_crr)
      ]
      projection_tester(
        prjs_cvvs[[tstsetup]],
        refmod_expected = refmods[[args_prj_cvvs[[tstsetup]]$tstsetup_ref]],
        solterms_expected = solterms_expected_crr,
        nprjdraws_expected = args_prj_cvvs[[tstsetup]]$nclusters,
        with_clusters = TRUE,
        const_wdraws_prj_expected = FALSE,
        info_str = tstsetup
      )
      # Check that projecting from the "vsel" object onto a single submodel
      # gives the same output as projecting the reference model onto that
      # submodel directly:
      tstsetup_tries <- grep(
        paste0("^",
               gsub("\\.", "\\\\.",
                    sub("(with_offs|without_offs).*", "\\1", tstsetup)),
               ".*\\.clust$"),
        names(prjs),
        value = TRUE
      )
      match_prj <- sapply(tstsetup_tries, function(tstsetup_try) {
        setequal(solterms_expected_crr, prjs[[tstsetup_try]]$solution_terms) &&
          args_prj_cvvs[[tstsetup]]$prj_nm == args_prj[[tstsetup_try]]$prj_nm
      })
      tstsetup_match_prj <- tstsetup_tries[match_prj]
      if (length(tstsetup_match_prj) == 1) {
        if (identical(prjs_cvvs[[tstsetup]]$solution_terms,
                      prjs[[tstsetup_match_prj]]$solution_terms)) {
          expect_identical(prjs_cvvs[[tstsetup]], prjs[[tstsetup_match_prj]],
                           ignore.environment = TRUE, info = tstsetup)
        } else {
          expect_setequal(prjs_cvvs[[tstsetup]]$solution_terms,
                          prjs[[tstsetup_match_prj]]$solution_terms)
          expect_equal(
            lapply(seq_along(prjs_cvvs[[tstsetup]]$outdmin), function(s_idx) {
              outdmin_s <- prjs_cvvs[[tstsetup]]$outdmin[[s_idx]]
              outdmin_s$beta <- outdmin_s$beta[
                rownames(prjs[[tstsetup_match_prj]]$outdmin[[s_idx]]$beta), ,
                drop = FALSE
              ]
              outdmin_s$formula <- prjs[[tstsetup_match_prj]]$outdmin[[s_idx]][[
                "formula"
              ]]
              outdmin_s$x <- outdmin_s$x[
                , colnames(prjs[[tstsetup_match_prj]]$outdmin[[s_idx]]$x),
                drop = FALSE
              ]
              return(outdmin_s)
            }),
            prjs[[tstsetup_match_prj]]$outdmin,
            tolerance = 1e1 * .Machine$double.eps, info = tstsetup
          )
          prj_nms <- names(prjs_cvvs[[tstsetup]])
          expect_identical(prj_nms, names(prjs[[tstsetup_match_prj]]),
                           info = tstsetup)
          prj_el_excl <- !prj_nms %in% c("solution_terms", "outdmin")
          expect_equal(prjs_cvvs[[tstsetup]][prj_el_excl],
                       prjs[[tstsetup_match_prj]][prj_el_excl],
                       tolerance = .Machine$double.eps, info = tstsetup)
        }
      } else if (length(tstsetup_match_prj) > 1) {
        stop("Unexpected number of matches.")
      }
    } else {
      proj_list_tester(
        prjs_cvvs[[tstsetup]],
        len_expected = length(nterms_crr),
        is_seq = all(diff(nterms_crr) == 1),
        info_str = tstsetup,
        refmod_expected = refmods[[args_prj_cvvs[[tstsetup]]$tstsetup_ref]],
        nprjdraws_expected = args_prj_cvvs[[tstsetup]]$nclusters,
        with_clusters = TRUE,
        const_wdraws_prj_expected = FALSE,
        prjdraw_weights_expected = prjs_cvvs[[tstsetup]][[1]]$wdraws_prj
      )
    }
  }
})

test_that("`object` not of class \"vsel\" and missing `solution_terms` fails", {
  skip_if_not(length(fits) > 0)
  expect_error(
    project(fits[[1]]),
    paste("^Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.$")
  )
  expect_error(
    project(refmods[[1]]),
    paste("^Please provide an `object` of class \"vsel\" or use argument",
          "`solution_terms`\\.$")
  )
})

# nterms ------------------------------------------------------------------

test_that("invalid `nterms` fails", {
  skip_if_not(run_vs)
  tstsetups <- head(grep("\\.glm\\.gauss", names(vss), value = TRUE), 1)
  for (tstsetup in tstsetups) {
    for (nterms_crr in nterms_unavail) {
      expect_error(project(vss[[!!tstsetup]], nterms = !!nterms_crr),
                   paste("Cannot perform the projection with", max(nterms_crr),
                         "variables"))
    }
    for (nterms_crr in list(neg = -1, char = "a", dafr = dat)) {
      expect_error(project(vss[[!!tstsetup]], nterms = !!nterms_crr),
                   "must contain non-negative values")
    }
  }
})

# seed --------------------------------------------------------------------

test_that("non-clustered projection does not require a seed", {
  skip_if_not(run_prj)
  # This test is important to ensure that we don't have to set a seed where we
  # don't expect it to be necessary.
  tstsetups <- grep("\\.noclust$|\\.default_ndr_ncl$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    p_orig <- prjs[[tstsetup]]
    rand_new1 <- runif(1) # Just to advance `.Random.seed[2]`.
    if (args_prj_i$prj_nm == "augdat" && args_prj_i$fam_nm == "cumul" &&
        !any(grepl("\\|", args_prj_i$solution_terms))) {
      warn_expected <- "non-integer #successes in a binomial glm!"
    } else if (!is.null(args_prj_i$avoid.increase) &&
               any(grepl("\\|", args_prj_i$solution_terms))) {
      warn_expected <- warn_mclogit
    } else {
      warn_expected <- NA
    }
    expect_warning(
      p_new <- do.call(project, c(
        list(object = refmods[[args_prj_i$tstsetup_ref]]),
        excl_nonargs(args_prj_i, nms_excl_add = "seed")
      )),
      warn_expected
    )
    if (args_prj_i$mod_nm %in% c("glmm", "gamm") &&
        any(grepl("\\|", args_prj_i$solution_terms))) {
      if (getOption("projpred.mlvl_pred_new", FALSE)) {
        # In this case, the multilevel submodel fitters (fit_glmer_callback(),
        # fit_gamm_callback(), fit_cumul_mlvl(), fit_categ_mlvl()) should still
        # be deterministic, but the prediction from the fitted submodels is not
        # (because of the group-level effects drawn randomly by repair_re() (for
        # all group levels; here, only the existing ones are relevant)). Thus,
        # we cannot test the whole project() output, but need to restrict
        # ourselves to the output of as.matrix.projection().
        if (args_prj_i$mod_nm == "gamm") {
          # Skipping GAMMs because of issue #131.
          # TODO (GAMMs): Fix this.
          next
        }
        prjmat_orig <- as.matrix(p_orig)
        prjmat_new <- as.matrix(p_new)
        if (args_prj_i$fam_nm == "gauss" || args_prj_i$prj_nm == "latent") {
          # The projected dispersion parameter is affected by the group-level
          # effects drawn randomly by repair_re() (for all group levels):
          prjmat_new[, "sigma"] <- prjmat_orig[, "sigma"]
        }
        expect_equal(prjmat_new, prjmat_orig, info = tstsetup,
                     tolerance = .Machine$double.eps)
        # To facilitate the `if` conditions here:
        p_new <- p_orig
      } else if (args_prj_i$prj_nm == "augdat" &&
                 args_prj_i$fam_nm == "cumul" && args_prj_i$mod_nm == "glmm") {
        for (idx_s in seq_along(p_new$outdmin)) {
          if (!is.null(p_new$outdmin[[idx_s]][["L"]])) {
            # We could also use `"sparseMatrix"` instead of `"Matrix"`:
            expect_equal(as(p_new$outdmin[[idx_s]][["L"]], "Matrix"),
                         as(p_orig$outdmin[[idx_s]][["L"]], "Matrix"),
                         info = tstsetup)
            p_new$outdmin[[idx_s]][["L"]] <- p_orig$outdmin[[idx_s]][["L"]]
          }
        }
      }
    }
    expect_equal(p_new, p_orig, info = tstsetup)
  }
})

test_that("`seed` works (and restores the RNG state afterwards)", {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.glm\\.gauss.*\\.solterms_x\\.clust$", names(prjs),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    p_orig <- prjs[[tstsetup]]
    rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_new1 <- .Random.seed
    p_new <- do.call(project, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]],
           seed = args_prj_i$seed + 10L),
      excl_nonargs(args_prj_i, nms_excl_add = "seed")
    ))
    .Random.seed_new2 <- .Random.seed
    rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_repr1 <- .Random.seed
    p_repr <- do.call(project, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]]),
      excl_nonargs(args_prj_i)
    ))
    .Random.seed_repr2 <- .Random.seed
    # Expected equality:
    expect_equal(p_repr, p_orig, info = tstsetup)
    expect_equal(.Random.seed_new2, .Random.seed_new1, info = tstsetup)
    expect_equal(.Random.seed_repr2, .Random.seed_repr1, info = tstsetup)
    # Expected inequality:
    expect_false(isTRUE(all.equal(p_new, p_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(rand_new, rand_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_repr2, .Random.seed_new2)),
                 info = tstsetup)
  }
})

# regul -------------------------------------------------------------------

test_that("for GLMs, `regul` has an expected effect", {
  skip_if_not(run_prj)
  regul_tst <- c(regul_default, 1e-1, 1e2)
  stopifnot(regul_tst[1] == regul_default)
  stopifnot(all(diff(regul_tst) > 0))
  tstsetups <- grep("\\.glm\\..*\\.clust$", names(prjs), value = TRUE)
  tstsetups <- grep(fam_nms_aug_regex, tstsetups, value = TRUE, invert = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    ndr_ncl <- ndr_ncl_dtls(args_prj_i)

    # Calculate the objects for which to run checks:
    ssq_regul_alpha <- rep(NA, length(regul_tst))
    ssq_regul_beta <- rep(NA, length(regul_tst))
    for (j in seq_along(regul_tst)) {
      # Run project() if necessary:
      if (regul_tst[j] == regul_default) {
        prj_regul <- prjs[[tstsetup]]
      } else {
        prj_regul <- do.call(project, c(
          list(object = refmods[[args_prj_i$tstsetup_ref]],
               regul = regul_tst[j]),
          excl_nonargs(args_prj_i)
        ))
        projection_tester(
          prj_regul,
          refmod_expected = refmods[[args_prj_i$tstsetup_ref]],
          solterms_expected = args_prj_i$solution_terms,
          nprjdraws_expected = ndr_ncl$nprjdraws,
          with_clusters = ndr_ncl$ndr_ncl_nm == "nclusters",
          const_wdraws_prj_expected = !ndr_ncl$clust_used,
          info_str = paste(tstsetup, j, sep = "__")
        )
      }

      # Run as.matrix.projection():
      if (ndr_ncl$clust_used) {
        # Clustered projection, so we expect a warning:
        warn_prjmat_expect <- "The projected draws have different .*weights"
      } else {
        warn_prjmat_expect <- NA
      }
      expect_warning(prjmat <- as.matrix(prj_regul, nm_scheme = "brms"),
                     warn_prjmat_expect, info = tstsetup)

      # Reduce to only those columns which are necessary here:
      prjmat <- prjmat[, grep("^b_", colnames(prjmat)), drop = FALSE]

      # Posterior means:
      prjmat_mean <- colMeans(prjmat)

      # Calculate the Euclidean norm for intercept and coefficients:
      ssq_regul_alpha[j] <- prjmat_mean["b_Intercept"]^2
      coef_colnms <- setdiff(names(prjmat_mean), "b_Intercept")
      if (length(coef_colnms) > 0) {
        ssq_regul_beta[j] <- sum(prjmat_mean[coef_colnms]^2)
      }
    }

    # Checks:
    if (length(args_prj_i$solution_terms) == 0) {
      # For an intercept-only model:
      expect_length(unique(ssq_regul_alpha), 1)
      stopifnot(all(is.na(ssq_regul_beta)))
    } else {
      # All other (i.e., not intercept-only) models (note: as discussed at issue
      # #169, the intercept is not tested here to stay the same):
      for (j in seq_along(ssq_regul_beta)[-1]) {
        expect_lt(ssq_regul_beta[!!j], ssq_regul_beta[j - 1])
      }
    }
  }
})
