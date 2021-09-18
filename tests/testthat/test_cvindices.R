context("cv-indices")

# tests for cvfolds and cv_ids

test_that("k is valid", {
  cvfuns <- c(cvfolds, cv_ids)
  for (cvfun in cvfuns) {
    expect_error(cvfun(n = 10, K = 1000),
                 "cannot exceed n")
    expect_error(cvfun(n = 10, K = 1),
                 "must be at least 2")
    expect_error(cvfun(n = 10, K = c(4, 9)),
                 "a single integer value")
    expect_error(cvfun(n = 10, K = "a"),
                 "a single integer value")
  }
})

test_that("cvfolds produces sensible results", {
  out <- cvfolds(n = 10, K = 3)
  expect_equal(length(out), 10)
  expect_equal(min(out), 1)
  expect_equal(max(out), 3)
})

test_that("cv_ids checks the \"out\" argument", {
  expect_error(cv_ids(n = 10, K = 3, out = "zzz"),
               "'arg' should be one of")
  expect_error(cv_ids(n = 10, K = 3, out = c("yyy", "zzz")),
               "must be of length 1")
  expect_error(cv_ids(n = 10, K = 3, out = 12),
               "must be NULL or a character vector")
})

test_that("cv_ids produces sensible results with out=\"foldwise\"", {
  out <- cv_ids(n = 10, K = 3, out = "foldwise")
  expect_equal(length(out), 3)
  expect_named(out[[1]], c("tr", "ts"))
})

test_that("cv_ids produces sensible results with out=\"indices\"", {
  out <- cv_ids(n = 10, K = 3, out = "indices")
  expect_equal(length(out), 2)
  expect_named(out, c("tr", "ts"))
  expect_equal(length(out$tr), 3)
  expect_equal(length(out$ts), 3)
})
