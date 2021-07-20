context("summary(), plot(), suggest_size()")

vs_funs <- list(summary.vsel, plot.vsel, suggest_size.vsel)

test_that("error if `object` is incorrect", {
  for (fun in vs_funs) {
    expect_error(fun(NULL), "is not a variable selection object")
    expect_error(fun(fits$glm$gauss), "is not a variable selection object")
    expect_error(fun(refmods$glm$gauss), "is not a variable selection object")
    expect_error(fun(prjs[[1]]), "is not a variable selection object")
    expect_error(fun(prjs_vs[[1]]), "is not a variable selection object")
    expect_error(fun(prjs_cvvs[[1]]), "is not a variable selection object")
  }
})

test_that("error if `stats` is incorrect", {
  for (fun in vs_funs) {
    expect_error(fun(vss[[1]], stat = NULL),
                 "specified as NULL")
    expect_error(fun(vss[[1]], stat = NA),
                 "not recognized")
    expect_error(fun(vss[[1]], stat = "zzz"),
                 "not recognized")
    expect_error(fun(vss[[1]], stat = "acc"),
                 "available only for the binomial family")
    expect_error(fun(vss[[1]], stat = "auc"),
                 "available only for the binomial family")
  }
})

# summary() ---------------------------------------------------------------

context("summary()")

test_that("error if `baseline` is incorrect", {
  skip_if_not(run_vs)
  tstsetups <- grep("^glm\\.gauss\\.default_meth", names(vss), value = TRUE)[1]
  for (tstsetup in tstsetups) {
    expect_error(
      summary(vss[[tstsetup]], baseline = "zzz"),
      "^Argument 'baseline' must be either 'ref' or 'best'\\.$",
      info = tstsetup
    )
  }
})

test_that(paste(
  "`object` of class \"vsel\" (created by cv_varsel()), `stats`, and `type`",
  "work"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(cvvss)) {
    fam_crr <- args_cvvs[[tstsetup]]$fam_nm
    stats_crr <- switch(fam_crr,
                        "gauss" = valid_stats_gauss,
                        "binom" = valid_stats_binom,
                        valid_stats_all)
    smmry <- summary(cvvss[[tstsetup]],
                     stats = stats_crr,
                     type = type_tst)
    smmry_sel_tester(
      smmry$selection,
      stats_expected = stats_crr,
      type_expected = type_tst,
      cv_method_expected = args_cvvs[[tstsetup]]$cv_method %ORifNULL% "LOO",
      solterms_expected = cvvss[[tstsetup]]$solution_terms,
      info_str = tstsetup
    )
  }
})

# TODO:
test_that("summary works with reference models", {
  for (i in seq_along(vsref_list)) {
    for (j in seq_along(vsref_list[[i]])) {
      vs <- vsref_list[[i]][[j]]
      if (vs$family$family == "gaussian") {
        stats_crr <- valid_stats_gauss
      } else {
        stats_crr <- valid_stats_binom
      }
      stats <- summary(vs, stats = stats_crr)$selection
      expect_true(is.data.frame(stats))
    }
  }
})

# print() -----------------------------------------------------------------

context("print()")

test_that("print() works as expected", {

  skip_on_cran()
  # default rounding
  expect_output(out <- print(vs_list[[1]][[1]]))
  expect_equal(out$selection$elpd, round(out$selection$elpd, 2),
               tolerance = 1e-3
  )
  expect_output(out <- print(cvs_list[[1]][[1]]))
  expect_equal(out$selection$elpd, round(out$selection$elpd, 2),
               tolerance = 1e-3
  )

  # rounding to 4 decimal places
  expect_output(out <- print(vs_list[[1]][[1]], digits = 4))
  expect_equal(out$selection$elpd, round(out$selection$elpd, 4),
               tolerance = 1e-3
  )
  expect_output(out <- print(cvs_list[[1]][[1]], digits = 4))
  expect_equal(out$selection$elpd, round(out$selection$elpd, 4),
               tolerance = 1e-3
  )
  # options to summary
  expect_output(out <- print(vs_list[[1]][[1]],
                             nterms_max = 3,
                             stats = "mse"
  ))
  expect_equal(nrow(out$selection) - 1, 3)
  expect_named(out$selection, c(
    "size", "solution_terms",
    "mse", "se",
    "diff", "diff.se"
  ))

  expect_output(out <- print(cvs_list[[1]][[1]],
                             nterms_max = 3,
                             stats = "mse"
  ))
  expect_equal(nrow(out$selection) - 1, 3)
  expect_named(out$selection, c(
    "size", "solution_terms",
    paste0("mse.", tolower(out$cv_method)), "se",
    "diff", "diff.se"
    # "pct_solution_terms_cv"
  ))
})

# plot() ------------------------------------------------------------------

context("plot()")

test_that("plotting works", {
  expect_s3_class(plot(vs_list[[1]][[1]]), "ggplot")
  expect_visible(plot(vs_list[[1]][[1]], nterms_max = 3))
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
