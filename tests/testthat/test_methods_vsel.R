# Common tests ------------------------------------------------------------

context("summary(), plot(), suggest_size()")

test_that("error if `object` is invalid", {
  objs_invalid <- nlist(
    NULL,
    fit = fits[[1]][[1]],
    refmod = refmods[[1]][[1]],
    prj = prjs[[1]],
    prj_vs = prjs_vs[[1]],
    prj_cvvs = prjs_cvvs[[1]]
  )
  for (obj_nm in names(objs_invalid)) {
    for (vsel_fun in vsel_funs) {
      expect_error(get(vsel_fun, mode = "function")(objs_invalid[[obj_nm]]),
                   "is not a variable selection object",
                   info = paste(obj_nm, vsel_fun, sep = "__"))
    }
  }
})

test_that("error if `stats` is invalid", {
  tstsetup <- grep("\\.gauss\\.", names(vss), value = TRUE)[1]
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

test_that("error if `baseline` is invalid", {
  skip_if_not(run_vs)
  for (tstsetup in names(vss)[1]) {
    expect_error(summary(vss[[tstsetup]], baseline = "zzz"),
                 "^Argument 'baseline' must be either 'ref' or 'best'\\.$",
                 info = tstsetup)
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by varsel()), `nterms_max`, `stats`,",
  "`type`, and `digits` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(smmrys_vs)) {
    tstsetup_vs <- args_smmry_vs[[tstsetup]]$tstsetup_vsel
    smmry_tester(
      smmrys_vs[[tstsetup]],
      vsel_expected = vss[[tstsetup_vs]],
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
  "`type`, and `digits` work"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(smmrys_cvvs)) {
    tstsetup_cvvs <- args_smmry_cvvs[[tstsetup]]$tstsetup_vsel
    smmry_tester(
      smmrys_cvvs[[tstsetup]],
      vsel_expected = cvvss[[tstsetup_cvvs]],
      info_str = tstsetup,
      stats_expected = args_smmry_cvvs[[tstsetup]]$stats,
      type_expected = args_smmry_cvvs[[tstsetup]]$type,
      nterms_max_expected = args_smmry_cvvs[[tstsetup]]$nterms_max,
      cv_method_expected = args_cvvs[[tstsetup]]$cv_method %ORifNULL% "LOO",
      solterms_expected = cvvss[[tstsetup_cvvs]]$solution_terms
    )
  }
})

# print() -----------------------------------------------------------------

context("print()")

test_that(paste(
  "`x` of class \"vsel\" (created by varsel()) and passing arguments to",
  "summary.vsel() works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(smmrys_vs)[1]) {
    args_smmry_vs_i <- args_smmry_vs[[tstsetup]]
    expect_output(
      print_obj <- do.call(print, c(
        list(x = vss[[args_smmry_vs_i$tstsetup_vsel]]),
        args_smmry_vs_i[setdiff(names(args_smmry_vs_i), c("tstsetup_vsel"))]
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
  for (tstsetup in names(smmrys_cvvs)[1]) {
    args_smmry_cvvs_i <- args_smmry_cvvs[[tstsetup]]
    expect_output(
      print_obj <- do.call(print, c(
        list(x = cvvss[[args_smmry_cvvs_i$tstsetup_vsel]]),
        args_smmry_cvvs_i[setdiff(names(args_smmry_cvvs_i), c("tstsetup_vsel"))]
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
  for (tstsetup in names(vss)[1]) {
    plot_obj <- plot(vss[[tstsetup]], nterms_max = nterms_avail$single)
    expect_s3_class(plot_obj, "ggplot")
    expect_visible(plot_obj, label = tstsetup)
  }
})

test_that("`x` of class \"vsel\" (created by cv_varsel()) works", {
  skip_if_not(run_cvvs)
  for (tstsetup in names(cvvss)[1]) {
    plot_obj <- plot(cvvss[[tstsetup]], nterms_max = nterms_avail$single)
    expect_s3_class(plot_obj, "ggplot")
    expect_visible(plot_obj, label = tstsetup)
  }
})

test_that("invalid 'baseline' arguments are rejected", {
  expect_error(
    plot(vs_list[[1]][[1]], baseline = "zzz"),
    "Argument 'baseline' must be either 'ref' or 'best'"
  )
})

test_that("the value of nterms_max is valid", {
  expect_error(
    plot(vs_list[[1]][[1]], nterms_max = 0),
    "nterms_max must be at least 1"
  )
})

## test_that("nterms_max is capped to the largest model size", {
##   expect_equal(
##     plot(vs_list[[1]][[1]]),
##     plot(vs_list[[1]][[1]], nterms_max = 1000)
##   )
## })

# suggest_size() ----------------------------------------------------------

context("suggest_size()")

test_that("suggest_size() checks the length of stat", {
  expect_error(suggest_size(vs_list[[1]][["gauss"]], stat = valid_stats_all),
               "Only one statistic")
})

test_that("suggest_size() works on all stats", {
  for (stat in valid_stats_gauss) {
    suggested_size <- suggest_size(vs_list[[1]][["gauss"]], stat = stat)
    expect_true(!is.na(suggested_size))
    expect_true(suggested_size >= 0)
  }
  for (stat in valid_stats_binom) {
    suggested_size <- suggest_size(vs_list[[1]][["binom"]], stat = stat)
    expect_true(!is.na(suggested_size))
    expect_true(suggested_size >= 0)
  }
})
