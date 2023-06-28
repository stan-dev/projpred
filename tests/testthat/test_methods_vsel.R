# Common tests ------------------------------------------------------------

context("summary(), plot(), suggest_size()")

test_that("invalid `object` fails", {
  objs_invalid <- nlist(
    NULL,
    some_numbers = 1:3,
    some_letters = head(letters, 17)
  )
  if (length(fits)) {
    objs_invalid <- c(
      objs_invalid,
      nlist(NULL,
            fit = fits[[1]],
            refmod = refmods[[1]])
    )
  }
  if (run_prj) {
    objs_invalid <- c(
      objs_invalid,
      list(prj = prjs[[1]])
    )
  }
  if (run_vs) {
    objs_invalid <- c(objs_invalid, list(prj_vs = prjs_vs[[1]]))
  }
  if (run_cvvs) {
    objs_invalid <- c(objs_invalid, list(prj_cvvs = prjs_cvvs[[1]]))
  }
  for (obj_nm in names(objs_invalid)) {
    for (vsel_fun in vsel_funs) {
      expect_error(get(vsel_fun, mode = "function")(objs_invalid[[obj_nm]]),
                   "is not a variable selection object",
                   info = paste(obj_nm, vsel_fun, sep = "__"))
    }
  }
})

test_that("invalid `stats` fails", {
  skip_if_not(run_vs)
  tstsetup <- head(grep("\\.gauss\\.", names(vss), value = TRUE), 1)
  stats_invalid <- nlist(NULL, NA, "zzz", "acc", "auc")
  err_expected <- as.list(c(
    "specified as NULL",
    rep("not recognized", 2),
    rep("available only for the binomial family", 2)
  ))
  names(err_expected) <- names(stats_invalid)
  for (stat_nm in names(stats_invalid)) {
    for (vsel_fun in vsel_funs) {
      expect_error(
        get(vsel_fun, mode = "function")(vss[[tstsetup]],
                                         stat = stats_invalid[[stat_nm]]),
        err_expected[[stat_nm]],
        info = paste(tstsetup, stat_nm, vsel_fun, sep = "__")
      )
    }
  }
})

# summary() ---------------------------------------------------------------

context("summary()")

test_that("invalid `baseline` fails", {
  skip_if_not(run_vs)
  for (tstsetup in head(names(vss), 1)) {
    expect_error(summary(vss[[tstsetup]], baseline = "zzz"),
                 "^Argument 'baseline' must be either 'ref' or 'best'\\.$",
                 info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by varsel()), `nterms_max`, `stats`,",
  "and `type` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(smmrys_vs)) {
    tstsetup_vs <- args_smmry_vs[[tstsetup]]$tstsetup_vsel
    smmry_tester(
      smmrys_vs[[tstsetup]],
      vsel_expected = vss[[tstsetup_vs]],
      resp_oscale_expected = args_smmry_vs[[tstsetup]]$resp_oscale %||% TRUE,
      search_trms_empty_size =
        length(args_vs[[tstsetup_vs]]$search_terms) &&
        all(grepl("\\+", args_vs[[tstsetup_vs]]$search_terms)),
      cumul_expected = args_smmry_vs[[tstsetup]]$cumulate,
      info_str = tstsetup,
      stats_expected = args_smmry_vs[[tstsetup]]$stats,
      type_expected = args_smmry_vs[[tstsetup]]$type,
      nterms_max_expected = args_smmry_vs[[tstsetup]]$nterms_max,
      solterms_expected = vss[[tstsetup_vs]]$solution_terms
    )
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by cv_varsel()), `nterms_max`, `stats`,",
  "and `type` work"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(smmrys_cvvs)) {
    tstsetup_cvvs <- args_smmry_cvvs[[tstsetup]]$tstsetup_vsel
    smmry_tester(
      smmrys_cvvs[[tstsetup]],
      vsel_expected = cvvss[[tstsetup_cvvs]],
      resp_oscale_expected = args_smmry_cvvs[[tstsetup]]$resp_oscale %||% TRUE,
      search_trms_empty_size =
        length(args_cvvs[[tstsetup_cvvs]]$search_terms) &&
        all(grepl("\\+", args_cvvs[[tstsetup_cvvs]]$search_terms)),
      cumul_expected = args_smmry_cvvs[[tstsetup]]$cumulate,
      info_str = tstsetup,
      stats_expected = args_smmry_cvvs[[tstsetup]]$stats,
      type_expected = args_smmry_cvvs[[tstsetup]]$type,
      nterms_max_expected = args_smmry_cvvs[[tstsetup]]$nterms_max,
      cv_method_expected =
        args_cvvs[[tstsetup_cvvs]]$cv_method %||% "LOO",
      solterms_expected = cvvss[[tstsetup_cvvs]]$solution_terms
    )
  }
})

# print() -----------------------------------------------------------------

context("print()")

test_that("`x` of class \"vselsummary\" (based on varsel()) works", {
  skip_if_not(run_vs)
  for (tstsetup in names(smmrys_vs)) {
    expect_message(
      expect_output(
        print_obj <- print(smmrys_vs[[tstsetup]]),
        "Family:.*Link function:.*Formula:.*Observations:",
        info = tstsetup
      ),
      NA
    )
    expect_identical(print_obj, smmrys_vs[[tstsetup]], info = tstsetup)
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(smmrys_vs[[tstsetup]], digits = 6)
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that("`x` of class \"vselsummary\" (based on cv_varsel())  works", {
  skip_if_not(run_cvvs)
  for (tstsetup in names(smmrys_cvvs)) {
    args_crr <- args_cvvs[[args_smmry_cvvs[[tstsetup]]$tstsetup_vsel]]
    if (isFALSE(args_crr$validate_search)) {
      mssg_expected <- NA
    } else {
      mssg_expected <- "Column.*contains the full-data predictor ranking"
    }
    expect_message(
      expect_output(
        print_obj <- print(smmrys_cvvs[[tstsetup]]),
        "Family:.*Link function:.*Formula:.*Observations:",
        info = tstsetup
      ),
      mssg_expected
    )
    expect_identical(print_obj, smmrys_cvvs[[tstsetup]], info = tstsetup)
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      width_orig <- options(width = 145)
      expect_snapshot({
        print(tstsetup)
        print(smmrys_cvvs[[tstsetup]], digits = 6)
      })
      options(width_orig)
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

test_that(paste(
  "`x` of class \"vsel\" (created by varsel()) and passing arguments to",
  "summary.vsel() works"
), {
  skip_if_not(run_vs)
  for (tstsetup in head(names(smmrys_vs), 1)) {
    args_smmry_vs_i <- args_smmry_vs[[tstsetup]]
    if (any(c("rmse", "auc") %in% args_smmry_vs_i$stats)) {
      smmry_seed <- list(seed = seed3_tst)
    } else {
      smmry_seed <- list()
    }
    expect_output(
      print_obj <- do.call(print, c(
        list(x = vss[[args_smmry_vs_i$tstsetup_vsel]]),
        excl_nonargs(args_smmry_vs_i),
        smmry_seed
      )),
      "Family:.*Link function:.*Formula:.*Observations:",
      info = tstsetup
    )
    expect_identical(print_obj, smmrys_vs[[tstsetup]], info = tstsetup)
  }
})

test_that(paste(
  "`x` of class \"vsel\" (created by cv_varsel()) and passing arguments to",
  "summary.vsel() works"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in head(names(smmrys_cvvs), 1)) {
    args_smmry_cvvs_i <- args_smmry_cvvs[[tstsetup]]
    if (any(c("rmse", "auc") %in% args_smmry_cvvs_i$stats)) {
      smmry_seed <- list(seed = seed3_tst)
    } else {
      smmry_seed <- list()
    }
    expect_output(
      print_obj <- do.call(print, c(
        list(x = cvvss[[args_smmry_cvvs_i$tstsetup_vsel]]),
        excl_nonargs(args_smmry_cvvs_i),
        smmry_seed
      )),
      "Family:.*Link function:.*Formula:.*Observations:",
      info = tstsetup
    )
    expect_identical(print_obj, smmrys_cvvs[[tstsetup]], info = tstsetup)
  }
})

# plot() ------------------------------------------------------------------

context("plot()")

test_that("`x` of class \"vsel\" (created by varsel()) works", {
  skip_if_not(run_vs)
  common_for_rk_NA <- c("tstsetup_vsel", "text_angle", "nterms_max",
                        "ranking_nterms_max")
  common_for_rk_max <- c("tstsetup_vsel", "text_angle", "nterms_max",
                         "ranking_abbreviate", "ranking_repel")
  for (tstsetup in names(plots_vs)) {
    args_plot_i <- args_plot_vs[[tstsetup]]
    if (identical(args_plot_i$ranking_nterms_max, NA)) {
      matches_tstsetup <- sapply(names(plots_vs), function(tstsetup2) {
        common_for_rk_NA <- intersect(common_for_rk_NA, names(args_plot_i))
        args_plot_i2 <- args_plot_vs[[tstsetup2]]
        common_for_rk_NA2 <- intersect(common_for_rk_NA, names(args_plot_i2))
        if (!setequal(common_for_rk_NA, common_for_rk_NA2)) {
          return(FALSE)
        }
        identical(args_plot_i[common_for_rk_NA], args_plot_i2[common_for_rk_NA])
      })
      if (any(matches_tstsetup)) {
        tstsetup_target <- names(which.max(matches_tstsetup))
      } else {
        tstsetup_target <- tstsetup
      }
    } else if (length(args_plot_i$ranking_nterms_max) &&
               identical(args_plot_i$ranking_nterms_max,
                         args_plot_i$nterms_max)) {
      matches_tstsetup <- sapply(names(plots_vs), function(tstsetup2) {
        common_for_rk_max <- intersect(common_for_rk_max, names(args_plot_i))
        args_plot_i2 <- args_plot_vs[[tstsetup2]]
        common_for_rk_max2 <- intersect(common_for_rk_max, names(args_plot_i2))
        if (!setequal(common_for_rk_max, common_for_rk_max2)) {
          return(FALSE)
        }
        identical(args_plot_i[common_for_rk_max],
                  args_plot_i2[common_for_rk_max])
      })
      if (any(matches_tstsetup)) {
        tstsetup_target <- names(which.max(matches_tstsetup))
      } else {
        tstsetup_target <- tstsetup
      }
    } else if (args_plot_i$ranking_colored || args_plot_i$cumulate) {
      not_common_here <- c("ranking_colored", "cumulate")
      matches_tstsetup <- sapply(names(plots_vs), function(tstsetup2) {
        args_plot_i2 <- args_plot_vs[[tstsetup2]]
        if (!setequal(names(args_plot_i), names(args_plot_i2))) {
          return(FALSE)
        }
        identical(args_plot_i[setdiff(names(args_plot_i), not_common_here)],
                  args_plot_i2[setdiff(names(args_plot_i2), not_common_here)])
      })
      if (any(matches_tstsetup)) {
        tstsetup_target <- names(which.max(matches_tstsetup))
      } else {
        tstsetup_target <- tstsetup
      }
    } else {
      tstsetup_target <- tstsetup
    }
    plot_vsel_tester(
      plots_vs[[tstsetup]],
      nterms_max_expected = args_plot_i$nterms_max %||%
        args_vs[[args_plot_i$tstsetup_vsel]]$nterms_max,
      rk_max_expected = args_plot_i$ranking_nterms_max,
      rk_expected = ranking(vss[[args_plot_i$tstsetup_vsel]])[["fulldata"]],
      abbv_expected = args_plot_i$ranking_abbreviate,
      abbv_args_expected = args_plot_i$ranking_abbreviate_args,
      info_str = tstsetup_target
    )
  }
})

test_that("`x` of class \"vsel\" (created by cv_varsel()) works", {
  skip_if_not(run_cvvs)
  common_for_rk_NA <- c("tstsetup_vsel", "text_angle", "nterms_max",
                        "ranking_nterms_max")
  common_for_rk_max <- c("tstsetup_vsel", "text_angle", "nterms_max",
                         "ranking_abbreviate", "ranking_repel",
                         "ranking_colored", "cumulate")
  for (tstsetup in names(plots_cvvs)) {
    args_plot_i <- args_plot_cvvs[[tstsetup]]
    if (identical(args_plot_i$ranking_nterms_max, NA)) {
      matches_tstsetup <- sapply(names(plots_cvvs), function(tstsetup2) {
        common_for_rk_NA <- intersect(common_for_rk_NA, names(args_plot_i))
        args_plot_i2 <- args_plot_cvvs[[tstsetup2]]
        common_for_rk_NA2 <- intersect(common_for_rk_NA, names(args_plot_i2))
        if (!setequal(common_for_rk_NA, common_for_rk_NA2)) {
          return(FALSE)
        }
        identical(args_plot_i[common_for_rk_NA], args_plot_i2[common_for_rk_NA])
      })
      if (any(matches_tstsetup)) {
        tstsetup_target <- names(which.max(matches_tstsetup))
      } else {
        tstsetup_target <- tstsetup
      }
    } else if (length(args_plot_i$ranking_nterms_max) &&
               identical(args_plot_i$ranking_nterms_max,
                         args_plot_i$nterms_max)) {
      matches_tstsetup <- sapply(names(plots_cvvs), function(tstsetup2) {
        common_for_rk_max <- intersect(common_for_rk_max, names(args_plot_i))
        args_plot_i2 <- args_plot_cvvs[[tstsetup2]]
        common_for_rk_max2 <- intersect(common_for_rk_max, names(args_plot_i2))
        if (!setequal(common_for_rk_max, common_for_rk_max2)) {
          return(FALSE)
        }
        identical(args_plot_i[common_for_rk_max],
                  args_plot_i2[common_for_rk_max])
      })
      if (any(matches_tstsetup)) {
        tstsetup_target <- names(which.max(matches_tstsetup))
      } else {
        tstsetup_target <- tstsetup
      }
    } else {
      tstsetup_target <- tstsetup
    }
    plot_vsel_tester(
      plots_cvvs[[tstsetup]],
      nterms_max_expected = args_plot_i$nterms_max %||%
        args_cvvs[[args_plot_i$tstsetup_vsel]]$nterms_max,
      rk_max_expected = args_plot_i$ranking_nterms_max,
      rk_expected = ranking(cvvss[[args_plot_i$tstsetup_vsel]])[["fulldata"]],
      abbv_expected = args_plot_i$ranking_abbreviate,
      abbv_args_expected = args_plot_i$ranking_abbreviate_args,
      info_str = tstsetup_target
    )
  }
})

test_that("invalid `baseline` fails", {
  skip_if_not(run_vs)
  for (tstsetup in head(names(vss), 1)) {
    expect_error(
      plot(vss[[tstsetup]], baseline = "zzz"),
      "^Argument 'baseline' must be either 'ref' or 'best'\\.$",
      info = tstsetup
    )
  }
})

test_that("invalid `nterms_max` fails", {
  skip_if_not(run_vs)
  for (tstsetup in head(names(vss), 1)) {
    expect_error(
      plot(vss[[tstsetup]], nterms_max = 0),
      "^nterms_max must be at least 1$",
      info = tstsetup
    )
  }
})

test_that(paste(
  "`nterms_max` is capped to the maximum model size (for varsel() output)"
), {
  skip_if_not(run_vs)
  tstsetups <- grep(paste("\\.default_nterms_max_smmry", "\\.default_rk_max",
                          "\\.default_abbv", "\\.default_repel", "\\.colFALSE",
                          "\\.cuFALSE", "\\.default_angle", sep = ".*"),
                    names(plots_vs), value = TRUE)
  for (tstsetup in tstsetups) {
    args_plot_i <- args_plot_vs[[tstsetup]]
    nterms_max_crr <- args_vs[[args_plot_i$tstsetup_vsel]]$nterms_max + 1L
    plot_capped <- plot(vss[[args_plot_i$tstsetup_vsel]],
                        nterms_max = nterms_max_crr)
    plot_vsel_tester(
      plot_capped,
      nterms_max_expected = nterms_max_crr,
      rk_max_expected = args_plot_i$ranking_nterms_max,
      rk_expected = ranking(vss[[args_plot_i$tstsetup_vsel]])[["fulldata"]],
      abbv_expected = args_plot_i$ranking_abbreviate,
      abbv_args_expected = args_plot_i$ranking_abbreviate_args,
      info_str = tstsetup
    )
  }
})

test_that(paste(
  "`nterms_max` is capped to the maximum model size (for cv_varsel() output)"
), {
  skip_if_not(run_cvvs)
  tstsetups <- grep(paste("\\.default_nterms_max_smmry", "\\.default_rk_max",
                          "\\.default_abbv", "\\.default_repel", "\\.colFALSE",
                          "\\.cuFALSE", "\\.default_angle", sep = ".*"),
                    names(plots_cvvs), value = TRUE)
  for (tstsetup in tstsetups) {
    args_plot_i <- args_plot_cvvs[[tstsetup]]
    nterms_max_crr <- args_cvvs[[args_plot_i$tstsetup_vsel]]$nterms_max + 1L
    plot_capped <- plot(cvvss[[args_plot_i$tstsetup_vsel]],
                        nterms_max = nterms_max_crr)
    plot_vsel_tester(
      plot_capped,
      nterms_max_expected = nterms_max_crr,
      rk_max_expected = args_plot_i$ranking_nterms_max,
      rk_expected = ranking(cvvss[[args_plot_i$tstsetup_vsel]])[["fulldata"]],
      abbv_expected = args_plot_i$ranking_abbreviate,
      abbv_args_expected = args_plot_i$ranking_abbreviate_args,
      info_str = tstsetup
    )
  }
})

# suggest_size() ----------------------------------------------------------

context("suggest_size()")

test_that("`stat` of invalid length fails", {
  stopifnot(length(stats_common) > 1)
  skip_if_not(run_vs)
  for (tstsetup in head(names(vss), 1)) {
    expect_error(
      suggest_size(vss[[tstsetup]], stat = stats_common),
      "^Only one statistic can be specified to suggest_size$",
      info = tstsetup
    )
  }
})

test_that("`stat` works", {
  skip_if_not(run_vs)
  tstsetups <- unname(unlist(lapply(mod_nms, function(mod_nm) {
    unlist(lapply(fam_nms, function(fam_nm) {
      head(grep(paste0("\\.", mod_nm, "\\.", fam_nm), names(args_smmry_vs),
                value = TRUE), 1)
    }))
  })))
  tstsetups <- union(tstsetups,
                     grep("\\.augdat\\..*\\.default_stats\\.",
                          names(args_smmry_vs), value = TRUE))
  for (tstsetup in tstsetups) {
    tstsetup_vs <- args_smmry_vs[[tstsetup]]$tstsetup_vsel
    fam_crr <- args_vs[[tstsetup_vs]]$fam_nm
    prj_crr <- args_vs[[tstsetup_vs]]$prj_nm
    stat_crr_nm <- switch(prj_crr,
                          "augdat" = "augdat_stats",
                          switch(fam_crr,
                                 "brnll" = "binom_stats",
                                 "binom" = "binom_stats",
                                 "common_stats"))
    stat_vec <- stats_tst[[stat_crr_nm]]$stats
    for (stat_crr in stat_vec) {
      if (stat_crr %in% c("rmse", "auc")) {
        suggsize_seed <- seed3_tst
      } else {
        suggsize_seed <- NULL
      }
      # Warnings are suppressed, but a suggested size of `NA` (because of a
      # search which was terminated too early) is tested below:
      suggsize <- suppressWarnings(
        suggest_size(vss[[tstsetup_vs]], stat = stat_crr, seed = suggsize_seed)
      )
      expect_length(suggsize, 1)
      if (!is.na(suggsize)) {
        expect_true(is.vector(suggsize, "numeric"),
                    info = paste(tstsetup, stat_crr, sep = "__"))
        expect_true(suggsize >= 0, info = paste(tstsetup, stat_crr, sep = "__"))
      } else {
        expect_identical(suggsize, NA,
                         info = paste(tstsetup, stat_crr, sep = "__"))
        expect_true(
          vss[[tstsetup_vs]]$nterms_max < vss[[tstsetup_vs]]$nterms_all,
          info = paste(tstsetup, stat_crr, sep = "__")
        )
      }
    }
  }
})

# ranking() ---------------------------------------------------------------

context("ranking()")

test_that("`object` of class `vsel` (created by varsel()) works", {
  skip_if_not(run_vs)
  for (tstsetup in names(rks_vs)) {
    tstsetup_vs <- args_rk_vs[[tstsetup]]$tstsetup_vsel
    nterms_max_expected_crr <- args_rk_vs[[tstsetup]][["nterms_max"]]
    ranking_tester(
      rks_vs[[tstsetup]],
      fulldata_expected = vss[[tstsetup_vs]][["solution_terms"]],
      foldwise_expected = NULL,
      nterms_max_expected = nterms_max_expected_crr,
      info_str = tstsetup
    )
  }
})

test_that("`object` of class `vsel` (created by cv_varsel()) works", {
  skip_if_not(run_cvvs)
  for (tstsetup in names(rks_cvvs)) {
    tstsetup_cvvs <- args_rk_cvvs[[tstsetup]]$tstsetup_vsel
    nterms_max_expected_crr <- args_rk_cvvs[[tstsetup]][["nterms_max"]]
    ranking_tester(
      rks_cvvs[[tstsetup]],
      fulldata_expected = cvvss[[tstsetup_cvvs]][["solution_terms"]],
      foldwise_expected = cvvss[[tstsetup_cvvs]][["solution_terms_cv"]],
      nterms_max_expected = nterms_max_expected_crr,
      info_str = tstsetup
    )
  }
})

# cv_proportions() --------------------------------------------------------

context("cv_proportions()")

test_that("`object` of class `ranking` (based on varsel() output) fails", {
  skip_if_not(run_vs)
  expect_length(prs_vs, 0)
})

test_that(paste(
  "`object` of class `ranking` (based on cv_varsel() output) works (if",
  "appropriate)"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(prs_cvvs)) {
    tstsetup_cvvs <- args_pr_cvvs[[tstsetup]]$tstsetup_vsel
    tstsetup_rk <- args_pr_cvvs[[tstsetup]]$tstsetup_rk
    nterms_max_expected_crr <- args_rk_cvvs[[tstsetup_rk]][["nterms_max"]]
    if (is.null(nterms_max_expected_crr)) {
      nterms_max_expected_crr <- args_cvvs[[tstsetup_cvvs]][["nterms_max"]]
    }
    cv_proportions_tester(
      prs_cvvs[[tstsetup]],
      cumulate_expected = args_pr_cvvs[[tstsetup]][["cumulate"]],
      nterms_max_expected = nterms_max_expected_crr,
      cnms_expected = cvvss[[tstsetup_cvvs]][["solution_terms"]],
      info_str = tstsetup
    )
  }
})

test_that("cv_proportions.vsel() is a shortcut", {
  skip_if_not(run_cvvs)
  for (tstsetup in names(prs_cvvs)) {
    args_pr_cvvs_i <- args_pr_cvvs[[tstsetup]]
    args_rk_cvvs_i <- args_rk_cvvs[[args_pr_cvvs_i$tstsetup_rk]]
    pr_from_vsel <- do.call(cv_proportions, c(
      list(object = cvvss[[args_pr_cvvs_i$tstsetup_vsel]]),
      excl_nonargs(args_pr_cvvs_i), excl_nonargs(args_rk_cvvs_i)
    ))
    expect_identical(pr_from_vsel, prs_cvvs[[tstsetup]], info = tstsetup)
  }
})

# Needed to clean up the workspace afterwards:
ls_bu <- ls()

ntrms <- 9L
rk_fdata <- rev(head(letters, ntrms))
rk_fwise <- do.call(rbind, lapply(seq_len(ntrms), function(idx_trm) {
  c(tail(rk_fdata, -idx_trm), head(rk_fdata, idx_trm))
}))
rk <- structure(
  list(fulldata = rk_fdata,
       foldwise = rk_fwise),
  class = "ranking"
)
# With `cumulate = FALSE`:
pr_cF <- cv_proportions(rk)
# With `cumulate = TRUE`:
pr_cT <- cv_proportions(rk, cumulate = TRUE)

test_that("`cumulate = TRUE` works", {
  pr_cT_ch <- structure(apply(pr_cF, 2, cumsum),
                        class = c("cv_proportions_cumul", "cv_proportions"))
  rownames(pr_cT_ch) <- paste0("<=", seq_len(nrow(pr_cT_ch)))
  expect_identical(pr_cT, pr_cT_ch)
})

test_that("ranking proportions are computed correctly", {
  skip_on_cran()
  tol_abs <- sqrt(.Machine$double.eps)
  # With `cumulate = FALSE`:
  cv_proportions_tester(pr_cF, nterms_max_expected = ntrms,
                        cnms_expected = rk_fdata, info_str = "cumulate = FALSE")
  expect_true(all(abs(pr_cF - 1 / ntrms) < tol_abs), info = "cumulate = FALSE")
  expect_true(all(abs(rowSums(pr_cF) - 1) < tol_abs), info = "cumulate = FALSE")
  expect_true(all(abs(colSums(pr_cF) - 1) < tol_abs), info = "cumulate = FALSE")

  # With `cumulate = TRUE`:
  cv_proportions_tester(pr_cT, cumulate_expected = TRUE,
                        nterms_max_expected = ntrms, cnms_expected = rk_fdata,
                        info_str = "cumulate = TRUE")
  pr_cT_expected <- matrix(1:ntrms / ntrms, nrow = ntrms, ncol = ntrms)
  class(pr_cT_expected) <- c("cv_proportions_cumul", "cv_proportions")
  expect_equal(pr_cT, pr_cT_expected, check.attributes = FALSE,
               tolerance = .Machine$double.eps, info = "cumulate = TRUE")
  has_1_colwise <- all(apply(pr_cT, 2, function(x_col) {
    any(abs(x_col - 1) < tol_abs)
  }))
  expect_true(has_1_colwise, info = "cumulate = TRUE")
})

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))

# plot.cv_proportions() ---------------------------------------------------

context("plot.cv_proportions()")

test_that("`x` of class `cv_proportions` works", {
  skip_if_not(run_cvvs)
  for (tstsetup in names(plotprs)) {
    expect_s3_class(plotprs[[tstsetup]], c("gg", "ggplot"))
    expect_visible(plotprs[[tstsetup]], label = tstsetup)
    if (run_snaps) {
      vdiffr::expect_doppelganger(tstsetup, plotprs[[tstsetup]])
    }
  }
})

test_that("plot.ranking() is a shortcut", {
  skip_if_not(run_cvvs)
  for (tstsetup in names(plotprs)) {
    args_plotpr_i <- args_plotpr[[tstsetup]]
    plotpr_from_rk <- do.call(plot, c(
      list(x = rks_cvvs[[args_plotpr_i$tstsetup_rk]]),
      excl_nonargs(args_plotpr_i),
      excl_nonargs(args_pr_cvvs[[args_plotpr_i$tstsetup_pr]])
    ))
    expect_s3_class(plotpr_from_rk, c("gg", "ggplot"))
    expect_visible(plotpr_from_rk, label = tstsetup)
    if (run_snaps) {
      vdiffr::expect_doppelganger(tstsetup, plotpr_from_rk)
    }
  }
})

# Needed to clean up the workspace afterwards:
ls_bu <- ls()

ntrms <- 20L
pr_dummy <- matrix(seq(0, 1, length.out = ntrms^2), nrow = ntrms, ncol = ntrms,
                   dimnames = list("size" = as.character(seq_len(ntrms)),
                                   "predictor" = paste0("x", seq_len(ntrms))))
class(pr_dummy) <- "cv_proportions"
prc_dummy <- pr_dummy
rownames(prc_dummy) <- paste0("<=", rownames(pr_dummy))
class(prc_dummy) <- c("cv_proportions_cumul", class(pr_dummy))

plotpr_dummy <- plot(pr_dummy)
plotprc_dummy <- plot(prc_dummy)

test_that("color gradient behaves as expected", {
  expect_s3_class(plotpr_dummy, c("gg", "ggplot"))
  expect_visible(plotpr_dummy)
  expect_s3_class(plotprc_dummy, c("gg", "ggplot"))
  expect_visible(plotprc_dummy)
  if (run_snaps) {
    vdiffr::expect_doppelganger("plotpr_dummy", plotpr_dummy)
    vdiffr::expect_doppelganger("plotprc_dummy", plotprc_dummy)
  }
})

test_that("`text_angle` works", {
  plotpr_dummy_angle <- plot(pr_dummy, text_angle = 60)
  plotprc_dummy_angle <- plot(prc_dummy, text_angle = 60)
  expect_s3_class(plotpr_dummy_angle, c("gg", "ggplot"))
  expect_visible(plotpr_dummy_angle)
  expect_s3_class(plotprc_dummy_angle, c("gg", "ggplot"))
  expect_visible(plotprc_dummy_angle)
  if (run_snaps) {
    vdiffr::expect_doppelganger("plotpr_dummy_angle", plotpr_dummy_angle)
    vdiffr::expect_doppelganger("plotprc_dummy_angle", plotprc_dummy_angle)
  }
})

test_that("the ggplot can be modified", {
  plotpr_dummy_mod <- plotpr_dummy + theme(legend.position = "none")
  plotprc_dummy_mod <- plotprc_dummy + theme(legend.position = "none")
  expect_s3_class(plotpr_dummy_mod, c("gg", "ggplot"))
  expect_visible(plotpr_dummy_mod)
  expect_s3_class(plotprc_dummy_mod, c("gg", "ggplot"))
  expect_visible(plotprc_dummy_mod)
  if (run_snaps) {
    vdiffr::expect_doppelganger("plotpr_dummy_mod", plotpr_dummy_mod)
    vdiffr::expect_doppelganger("plotprc_dummy_mod", plotprc_dummy_mod)
  }
})

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
