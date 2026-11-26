# Helpers of the tests (tests/testthat/helper-*.R and helpers/*.R) ----------

test_that("skip_if_no_setups() skips only for an empty list", {
  expect_condition(skip_if_no_setups(list()), class = "skip")
  expect_no_condition(skip_if_no_setups(list(a = 1)), class = "skip")
})

test_that("expect_warning_all() catches all matching warnings", {
  warn_twice <- function() {
    warning("first match")
    warning("second match")
    "value"
  }
  expect_success(val <- expect_warning_all(warn_twice(), "match"))
  expect_identical(val, "value")
  expect_success(expect_warning_all(warn_twice()))
  expect_failure(expect_warning_all("no warning", "match"))
  # (suppressWarnings() for the warnings that the expectation does not catch.)
  expect_failure(suppressWarnings(expect_warning_all(warn_twice(), "other")))
  expect_success(expect_warning_all("no warning", NA))
  expect_failure(suppressWarnings(expect_warning_all(warn_twice(), NA)))
  # A warning that does not match is not caught:
  expect_warning(
    expect_warning_all({
      warning("match")
      warning("other")
    }, "match"),
    "other"
  )
})

test_that("expect_all_true_tst() works as expect_true(all())", {
  expect_success(expect_all_true_tst(c(TRUE, TRUE)))
  expect_success(expect_all_true_tst(logical(0)))
  expect_success(expect_all_true_tst(matrix(TRUE, 2, 2)))
  expect_failure(expect_all_true_tst(c(TRUE, NA)),
                 "1 of 2 are not TRUE, at 2\\.")
  expect_failure(expect_all_true_tst(c(a = TRUE, b = FALSE, c = FALSE)),
                 "2 of 3 are not TRUE, at b, c\\.")
  expect_failure(expect_all_true_tst(c(a = FALSE, FALSE)), "at a, 2\\.")
  expect_failure(expect_all_true_tst(rep(FALSE, 12)),
                 "12 of 12 are not TRUE, at 1, 2, .*, 10, \\.\\.\\.")
  expect_failure(
    expect_all_true_tst(matrix(c(TRUE, TRUE, FALSE, TRUE), 2, 2)),
    "at \\[1, 2\\]\\."
  )
  expect_failure(expect_all_true_tst(FALSE, info = "xyz"), "xyz")
  expect_failure(expect_all_true_tst(FALSE, label = "lbl"), "lbl")
  expect_error(expect_all_true_tst(c(1, 1)), "must be logical")
  expect_error(expect_all_true_tst("TRUE"), "must be logical")
  x <- list(c(TRUE, TRUE), c(FALSE, TRUE))
  j <- 2
  # `!!` is evaluated before the expression, as in the testthat expectations
  # (without that, `x[[!!j]]` would be `x[[TRUE]]`, the first element):
  expect_failure(expect_all_true_tst(x[[!!j]]), "at 1\\.")
})

test_that("drop_mer_ptrs() keeps the differences between lme4 fits", {
  fml <- Reaction ~ Days + (1 | Subject)
  fit1 <- lme4::lmer(fml, data = lme4::sleepstudy)
  fit2 <- lme4::lmer(fml, data = lme4::sleepstudy)
  dat_other <- lme4::sleepstudy
  dat_other$Reaction[1] <- dat_other$Reaction[1] + 1
  fit_other <- lme4::lmer(fml, data = dat_other)
  expect_failure(expect_equal(fit1, fit2))
  ptr_before <- fit1@resp$Ptr
  expect_success(expect_equal(drop_mer_ptrs(list(a = fit1, b = list(fit1))),
                              drop_mer_ptrs(list(a = fit2, b = list(fit2))),
                              tolerance = 0))
  expect_failure(expect_equal(drop_mer_ptrs(fit1), drop_mer_ptrs(fit_other)))
  x <- structure(list(a = 1, b = list(fit1)), class = "foo")
  expect_s3_class(drop_mer_ptrs(x), "foo", exact = TRUE)
  # The fit itself keeps its pointer:
  expect_identical(fit1@resp$Ptr, ptr_before)
})
