# Common tests ------------------------------------------------------------

context("summary(), plot(), suggest_size()")

test_that("invalid `object` fails", {
  objs_invalid <- nlist(
    NULL,
    fit = fits[[1]],
    refmod = refmods[[1]]
  )
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
      search_trms_empty_size =
        length(args_vs[[tstsetup_vs]]$search_terms) &&
        all(grepl("\\+", args_vs[[tstsetup_vs]]$search_terms)),
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
      search_trms_empty_size =
        length(args_cvvs[[tstsetup_cvvs]]$search_terms) &&
        all(grepl("\\+", args_cvvs[[tstsetup_cvvs]]$search_terms)),
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
    expect_output(
      print_obj <- print(smmrys_vs[[tstsetup]]),
      "Family:.*Link function:.*Formula:.*Observations:",
      info = tstsetup
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
    expect_output(
      print_obj <- print(smmrys_cvvs[[tstsetup]]),
      "Family:.*Link function:.*Formula:.*Observations:",
      info = tstsetup
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
  for (tstsetup in head(names(vss), 1)) {
    plot_obj <- plot(vss[[tstsetup]], nterms_max = nterms_avail$single)
    expect_s3_class(plot_obj, "ggplot")
    expect_visible(plot_obj, label = tstsetup)
  }
})

test_that("`x` of class \"vsel\" (created by cv_varsel()) works", {
  skip_if_not(run_cvvs)
  for (tstsetup in head(names(cvvss), 1)) {
    plot_obj <- plot(cvvss[[tstsetup]], nterms_max = nterms_avail$single)
    expect_s3_class(plot_obj, "ggplot")
    expect_visible(plot_obj, label = tstsetup)
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

test_that("`nterms_max` is capped to the maximum model size", {
  skip_if_not(run_vs)
  for (tstsetup in head(names(vss), 1)) {
    expect_equal(
      plot(vss[[tstsetup]]),
      plot(vss[[tstsetup]], nterms_max = args_vs[[tstsetup]]$nterms_max + 1L),
      info = tstsetup
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
  for (tstsetup in tstsetups) {
    tstsetup_vs <- args_smmry_vs[[tstsetup]]$tstsetup_vsel
    fam_crr <- args_vs[[tstsetup_vs]]$fam_nm
    stat_crr_nm <- switch(fam_crr,
                          "brnll" = "binom_stats",
                          "binom" = "binom_stats",
                          "common_stats")
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
      expect_type(suggsize, "double")
      expect_length(suggsize, 1)
      if (!is.na(suggsize)) {
        expect_true(suggsize >= 0, info = paste(tstsetup, stat_crr, sep = "__"))
        if (stat_crr == "elpd") {
          expect_identical(suggsize, vss[[tstsetup_vs]]$suggested_size,
                           info = paste(tstsetup, stat_crr, sep = "__"))
        }
      } else {
        expect_true(
          vss[[tstsetup_vs]]$nterms_max < vss[[tstsetup_vs]]$nterms_all,
          info = paste(tstsetup, stat_crr, sep = "__")
        )
      }
    }
  }
})
