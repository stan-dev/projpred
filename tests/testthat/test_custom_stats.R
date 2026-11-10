context("custom_stat()")

# Custom reimplementations of built-in statistics:
cstat_mse <- custom_stat(function(y, mu, wobs) (mu - y)^2, name = "c_mse",
                         higher_better = FALSE, label = "Custom MSE")
# Returns a logical vector, which get_stat_custom() converts to numeric:
cstat_acc <- custom_stat(function(y, mu, wobs) round(mu) == y, name = "c_acc")
cstat_auc <- custom_stat(function(y, mu, wobs) .auc(cbind(y, mu, wobs)),
                         name = "c_auc", type = "aggregate")

# Compare the performance columns of the built-in statistic `stat_builtin`
# with those of the custom statistic `cstat`:
expect_perf_equal <- function(vs, stat_builtin, cstat, type, info,
                              els = c("submodels", "reference_model"), ...) {
  perf_builtin <- suppressMessages(
    performances(vs, stats = stat_builtin, type = type, ...)
  )
  perf_custom <- suppressMessages(
    performances(vs, stats = cstat, type = type, ...)
  )
  for (el in els) {
    perf_c <- perf_custom[[el]]
    names(perf_c) <- sub(paste0("^", cstat$name), stat_builtin, names(perf_c))
    expect_equal(perf_c, perf_builtin[[el]], tolerance = 1e-10,
                 info = paste(info, el, sep = "__"))
  }
}

# All `vsel` objects from the setup:
vsel_objs_all <- function() {
  objs <- list()
  if (run_vs) objs <- c(objs, vss)
  if (run_cvvs) objs <- c(objs, cvvss)
  return(objs)
}
# The `vsel` objects for which the custom statistics are available (with the
# default `resp_oscale = TRUE`):
vsel_objs_cstat <- function() {
  Filter(function(vs) {
    !vs$refmodel$family$for_augdat && is.null(vs$refmodel$family$cats)
  }, vsel_objs_all())
}
is_subsampled <- function(vs) {
  !is.null(vs$summaries_fast)
}
is_datafit <- function(vs) {
  inherits(vs$refmodel, "datafit")
}
# Binomial family (on the response scale) with 0/1 responses:
is_binom_unit <- function(vs) {
  fam <- vs$refmodel$family
  if (fam$for_latent) {
    fam_ch <- fam$family_oscale
    y <- vs$y_wobs_test$y_oscale
  } else {
    fam_ch <- fam$family
    y <- vs$y_wobs_test$y
  }
  identical(fam_ch, "binomial") && all(vs$y_wobs_test$wobs == 1) &&
    is.numeric(y) && all(y %in% c(0, 1))
}

test_that("custom_stat() checks its arguments", {
  f <- function(y, mu, wobs) mu
  expect_s3_class(custom_stat(f, name = "a"), "custom_stat")
  expect_error(custom_stat("f", name = "a"), "must be a function")
  expect_error(custom_stat(function(x) x, name = "a"),
               "must have arguments `y`, `mu` and `wobs`")
  for (nm in c("elpd", "auc", "size", "ranking_fulldata",
               "cv_proportions_diag")) {
    expect_error(custom_stat(f, name = nm), "must not be one of", info = nm)
  }
  for (nm in list("a b", "a.se", "_a", "1a", c("a", "b"), NA_character_, 1)) {
    expect_error(custom_stat(f, name = nm), "must be a single character",
                 info = paste(nm, collapse = ","))
  }
  expect_error(custom_stat(f, name = "a", type = "zzz"), "should be one of")
  expect_error(custom_stat(f, name = "a", higher_better = NA),
               "`higher_better` must be")
  expect_error(custom_stat(f, name = "a", label = 1), "`label` must be")
  expect_output(print(cstat_mse),
                "statistic `c_mse` \\(pointwise, lower is better\\)")
})

test_that("parse_stats() accepts strings, custom_stat objects and lists", {
  expect_identical(parse_stats(c("elpd", "mse")),
                   list(nms = c("elpd", "mse"), custom = list()))
  # Other atomic input is left for validate_vsel_object_stats():
  expect_identical(parse_stats(NA)$nms, NA_character_)
  expect_identical(parse_stats(1)$nms, "1")
  ps <- parse_stats(cstat_mse)
  expect_identical(ps$nms, "c_mse")
  expect_identical(ps$custom, list(c_mse = cstat_mse))
  ps <- parse_stats(list("elpd", cstat_acc, "mse"))
  expect_identical(ps$nms, c("elpd", "c_acc", "mse"))
  expect_identical(names(ps$custom), "c_acc")
  expect_error(parse_stats(NULL), "specified as NULL")
  expect_error(parse_stats(list("elpd", 1)), "single character string")
  expect_error(parse_stats(list(cstat_acc, cstat_acc)), "different names")
})

test_that("pointwise custom statistic reproduces \"mse\"", {
  skip_if_not(run_vs || run_cvvs)
  objs <- vsel_objs_cstat()
  expect_gt(length(objs), 0)
  for (tstsetup in names(objs)) {
    vs <- objs[[tstsetup]]
    if (is_subsampled(vs)) {
      # With subsampled LOO-CV, the built-in "mse" estimates the standard error
      # of the difference in another way. Without `deltas`, its uncertainty
      # interval is log-normal:
      expect_perf_equal(vs, "mse", cstat_mse, type = c("mean", "se", "diff"),
                        info = tstsetup)
      next
    }
    expect_perf_equal(vs, "mse", cstat_mse,
                      type = c("mean", "se", "diff", "diff.se", "diff.lower",
                               "diff.upper"),
                      info = tstsetup)
    # For the reference model compared with itself, the built-in "mse" gives a
    # standard error slightly above zero (cancellation), so compare only the
    # submodels:
    expect_perf_equal(vs, "mse", cstat_mse,
                      type = c("mean", "se", "lower", "upper"),
                      els = "submodels", deltas = TRUE,
                      info = paste(tstsetup, "deltas"))
    expect_identical(
      suppressMessages(suggest_size(vs, stat = cstat_mse, warnings = FALSE)),
      suppressMessages(suggest_size(vs, stat = "mse", warnings = FALSE)),
      info = tstsetup
    )
  }
})

test_that("pointwise custom statistic reproduces \"mse\" on latent scale", {
  skip_if_not(run_vs || run_cvvs)
  # This includes the latent projection with a non-`NULL` element
  # `family$cats`:
  objs <- Filter(function(vs) vs$refmodel$family$for_latent, vsel_objs_all())
  skip_if(length(objs) == 0, "no latent-projection `vsel` object in the setup")
  for (tstsetup in names(objs)) {
    type <- if (is_subsampled(objs[[tstsetup]])) {
      c("mean", "se", "diff")
    } else {
      c("mean", "se", "diff", "diff.se")
    }
    expect_perf_equal(objs[[tstsetup]], "mse", cstat_mse, type = type,
                      resp_oscale = FALSE, info = tstsetup)
  }
})

test_that("pointwise custom statistic reproduces \"acc\"", {
  skip_if_not(run_vs || run_cvvs)
  objs <- Filter(is_binom_unit, vsel_objs_cstat())
  expect_gt(length(objs), 0)
  for (tstsetup in names(objs)) {
    vs <- objs[[tstsetup]]
    expect_perf_equal(vs, "acc", cstat_acc,
                      type = c("mean", "se", "lower", "upper", "diff",
                               "diff.se", "diff.lower", "diff.upper"),
                      info = tstsetup)
    if (!is_datafit(vs)) {
      expect_perf_equal(vs, "acc", cstat_acc,
                        type = c("mean", "se", "lower", "upper"),
                        deltas = TRUE, baseline = "ref",
                        info = paste(tstsetup, "deltas"))
    }
    if (!is_subsampled(vs)) {
      expect_perf_equal(vs, "acc", cstat_acc,
                        type = c("mean", "se", "lower", "upper"),
                        deltas = TRUE, baseline = "best",
                        info = paste(tstsetup, "best"))
      expect_identical(
        suppressMessages(suggest_size(vs, stat = cstat_acc, warnings = FALSE)),
        suppressMessages(suggest_size(vs, stat = "acc", warnings = FALSE)),
        info = tstsetup
      )
    }
  }
})

test_that("aggregate custom statistic reproduces \"auc\"", {
  skip_if_not(run_vs || run_cvvs)
  objs <- Filter(function(vs) is_binom_unit(vs) && !is_subsampled(vs),
                 vsel_objs_cstat())
  expect_gt(length(objs), 0)
  for (tstsetup in names(objs)) {
    expect_perf_equal(objs[[tstsetup]], "auc", cstat_auc,
                      type = c("mean", "se", "lower", "upper", "diff",
                               "diff.se", "diff.lower", "diff.upper"),
                      B = 50, seed = 3, info = tstsetup)
  }
})

test_that("custom statistics get `y` as proportions for binomial `wobs > 1`", {
  skip_if_not(run_vs || run_cvvs)
  objs <- Filter(function(vs) {
    identical(vs$refmodel$family$family, "binomial") &&
      !vs$refmodel$family$for_latent && any(vs$y_wobs_test$wobs > 1) &&
      !is_subsampled(vs) && !is_datafit(vs)
  }, vsel_objs_cstat())
  skip_if(length(objs) == 0,
          "no binomial `vsel` object with `wobs > 1` in the setup")
  # Weighted mean of the squared errors of the proportions:
  wmse <- function(y, mu, wobs) sum(wobs * (mu - y)^2) / sum(wobs)
  cstat_wmse <- custom_stat(wmse, name = "c_wmse", type = "aggregate",
                            higher_better = FALSE)
  for (tstsetup in names(objs)) {
    vs <- objs[[tstsetup]]
    y_prop <- vs$y_wobs_test$y / vs$y_wobs_test$wobs
    wobs <- vs$y_wobs_test$wobs
    perf <- performances(vs, stats = list(cstat_wmse, cstat_mse),
                         type = c("mean", "se"), B = 20, seed = 1)
    val_expected <- vapply(vs$summaries$sub, function(s) {
      wmse(y_prop, s$mu, wobs)
    }, numeric(1))
    expect_equal(perf$submodels$c_wmse, val_expected, tolerance = 1e-12,
                 info = tstsetup)
    expect_equal(unname(perf$reference_model["c_wmse"]),
                 wmse(y_prop, vs$summaries$ref$mu, wobs), tolerance = 1e-12,
                 info = tstsetup)
    # The pointwise version gives the same mean:
    expect_equal(perf$submodels$c_mse, val_expected, tolerance = 1e-12,
                 info = tstsetup)
    expect_true(all(perf$submodels$c_wmse.se > 0), info = tstsetup)
  }
})

test_that("aggregate custom statistic fails with subsampled LOO-CV", {
  skip_if_not(run_cvvs)
  objs <- Filter(is_subsampled, vsel_objs_cstat())
  skip_if(length(objs) == 0, "no subsampled-LOO `vsel` object in the setup")
  expect_error(summary(objs[[1]], stats = cstat_auc),
               "aggregate custom statistic is not implemented")
})

test_that("custom statistic with a wrong return value fails", {
  skip_if_not(run_vs || run_cvvs)
  objs <- Filter(function(vs) !is_subsampled(vs), vsel_objs_cstat())
  expect_gt(length(objs), 0)
  cstat_short <- custom_stat(function(y, mu, wobs) mu[-1], name = "c_short")
  expect_error(summary(objs[[1]], stats = cstat_short),
               "must return a numeric vector of length")
  cstat_vec <- custom_stat(function(y, mu, wobs) mu, name = "c_vec",
                           type = "aggregate")
  expect_error(summary(objs[[1]], stats = cstat_vec),
               "must return a single number")
})

test_that("custom statistic works with an all-`NA` reference model", {
  skip_if_not(run_vs || run_cvvs)
  objs <- Filter(function(vs) {
    !is_subsampled(vs) && !vs$refmodel$family$for_latent
  }, vsel_objs_cstat())
  expect_gt(length(objs), 0)
  # A `datafit` has a logical `NA` vector as the reference model `mu`:
  vs <- objs[[1]]
  vs$summaries$ref$mu <- rep(NA, length(vs$summaries$ref$mu))
  # A function that returns a logical vector:
  cstat_close <- custom_stat(function(y, mu, wobs) abs(mu - y) < 1,
                             name = "c_close")
  cstat_close_agg <- custom_stat(function(y, mu, wobs) mean(abs(mu - y) < 1),
                                 name = "c_close_agg", type = "aggregate")
  for (cstat in list(cstat_close, cstat_close_agg)) {
    perf <- performances(vs, stats = cstat, type = c("mean", "se"), B = 20,
                         seed = 1)
    expect_true(all(is.na(perf$reference_model)), info = cstat$name)
    expect_false(anyNA(perf$submodels[[cstat$name]]), info = cstat$name)
  }
})

test_that("custom statistic fails for the augmented-data projection", {
  skip_if_not(run_vs || run_cvvs)
  objs <- Filter(function(vs) vs$refmodel$family$for_augdat, vsel_objs_all())
  skip_if(length(objs) == 0, "no augmented-data `vsel` object in the setup")
  expect_error(summary(objs[[1]], stats = cstat_mse),
               "augmented-data projection may not be combined")
})

test_that("custom statistic fails for the latent projection with `cats`", {
  skip_if_not(run_vs || run_cvvs)
  objs <- Filter(function(vs) {
    vs$refmodel$family$for_latent && !is.null(vs$refmodel$family$cats)
  }, vsel_objs_all())
  skip_if(length(objs) == 0,
          "no latent-projection `vsel` object with `cats` in the setup")
  expect_error(summary(objs[[1]], stats = cstat_mse),
               "a non-`NULL` element `family\\$cats` may not be combined")
})

test_that("plot.vsel() works with custom statistics and uses `label`", {
  skip_if_not(run_vs || run_cvvs)
  # get_strip_labels() exists since ggplot2 3.5.0:
  has_strip_fun <- utils::packageVersion("ggplot2") >= "3.5.0"
  objs <- Filter(function(vs) {
    is_binom_unit(vs) && !is_subsampled(vs) && !is_datafit(vs)
  }, vsel_objs_cstat())
  expect_gt(length(objs), 0)
  pp <- suppressMessages(plot(objs[[1]], stats = list("elpd", cstat_mse)))
  expect_s3_class(pp, "ggplot")
  expect_match(pp$labels$subtitle, "normal-approximation")
  if (has_strip_fun) {
    expect_identical(ggplot2::get_strip_labels(pp)$rows$statistic,
                     c("Custom MSE", "elpd"))
  }
  pp <- suppressMessages(plot(objs[[1]], stats = cstat_auc, deltas = "mixed",
                              B = 20, seed = 1))
  expect_s3_class(pp, "ggplot")
  expect_match(pp$labels$subtitle, "bootstrap")
  if (has_strip_fun) {
    expect_identical(ggplot2::get_strip_labels(pp)$rows$statistic, "c_auc")
  }
})
