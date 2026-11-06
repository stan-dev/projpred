# Helpers of the tests (tests/testthat/helpers/) -----------------------------

test_that("skip_if_no_setups() skips only for an empty list", {
  expect_condition(skip_if_no_setups(list()), class = "skip")
  expect_no_condition(skip_if_no_setups(list(a = 1)), class = "skip")
})
