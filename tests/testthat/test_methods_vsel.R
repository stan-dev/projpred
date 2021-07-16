# summary() ---------------------------------------------------------------

context("summary()")

valid_stats_all <- c("elpd", "mlpd")
valid_stats_gauss_only <- c("mse", "rmse")
valid_stats_binom_only <- c("acc", "auc")
valid_stats_gauss <- c(valid_stats_all, valid_stats_gauss_only)
valid_stats_binom <- c(valid_stats_all, valid_stats_binom_only)
vs_funs <- c(summary, plot, suggest_size)

## test_that("invalid objects are rejected", {
##   for (fun in vs_funs) {
##     expect_error(fun(NULL), "is not a variable selection object")
##     expect_error(fun(fit_gauss), "is not a variable selection object")
##   }
## })

## test_that("invalid stats are rejected", {
##   for (fun in vs_funs) {
##     expect_error(fun(vs_list[[1]][["gauss"]], stat = NULL),
##                  "specified as NULL")
##     expect_error(fun(vs_list[[1]][["gauss"]], stat = NA),
##                  "not recognized")
##     expect_error(fun(vs_list[[1]][["gauss"]], stat = "zzz"),
##                  "not recognized")
##     expect_error(fun(vs_list[[1]][["gauss"]], stat = "acc"),
##                  "available only for the binomial family")
##     expect_error(
##       fun(vs_list[[1]][["gauss"]], stat = "auc"),
##       "available only for the binomial family"
##     )
##   }
## })

## test_that("invalid 'baseline' arguments are rejected", {
##   expect_error(
##     summary(vs_list[[1]][["gauss"]], baseline = "zzz"),
##     "Argument 'baseline' must be either 'ref' or 'best'"
##   )
## })

test_that("summary output seems legit", {
  skip_on_cran()
  for (i in seq_along(cvs_list)) {
    for (j in seq_along(cvs_list[[i]])) {
      cvs <- cvs_list[[i]][[j]]
      if (cvs$family$family == "gaussian") {
        stats_str <- valid_stats_gauss
      } else if (cvs$family$family == "binomial") {
        stats_str <- valid_stats_binom
      } else {
        stats_str <- valid_stats_all
      }
      cv_method <- cvs_list[[i]][[j]]$cv_method
      stats <- summary(cvs,
                       stats = stats_str,
                       type = c("mean", "lower", "upper", "se")
      )$selection
      expect_true(nrow(stats) == nterms + 1)
      expect_true(all(c(
        "size", "solution_terms", paste0(stats_str, ".", tolower(cv_method)),
        paste0(stats_str, ".", c("se", "upper", "lower"))
      ) %in% names(stats)))
      expect_true(all(stats[, paste0("mlpd.", tolower(cv_method))] >
                        stats[, "mlpd.lower"]))
      expect_true(all(stats[, paste0("mlpd.", tolower(cv_method))] <
                        stats[, "mlpd.upper"]))
    }
  }
})

test_that("summary works with reference models", {
  for (i in seq_along(vsref_list)) {
    for (j in seq_along(vsref_list[[i]])) {
      vs <- vsref_list[[i]][[j]]
      if (vs$family$family == "gaussian") {
        stats_str <- valid_stats_gauss
      } else {
        stats_str <- valid_stats_binom
      }
      stats <- summary(vs, stats = stats_str)$selection
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
