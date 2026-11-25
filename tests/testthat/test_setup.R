test_that("is_known_setup_err() accepts an error only for its model types", {
  msg_bs <- "singular matrix in 'backsolve'. First zero in diagonal [1]"
  expect_true(is_known_setup_err(msg_bs, mod_nm = "gamm"))
  expect_false(is_known_setup_err(msg_bs, mod_nm = "glmm"))
  expect_false(is_known_setup_err(msg_bs, mod_nm = "gam"))
  expect_false(is_known_setup_err(msg_bs, mod_nm = "gamm", kfold = TRUE))
  msg_pw <- "pwrssUpdate did not converge in (maxit) iterations"
  expect_true(is_known_setup_err(msg_pw, mod_nm = "glmm", kfold = TRUE))
  expect_true(is_known_setup_err(msg_pw, mod_nm = "gamm", kfold = TRUE))
  expect_false(is_known_setup_err(msg_pw, mod_nm = "glm", kfold = TRUE))
  expect_false(is_known_setup_err(msg_pw, mod_nm = "gam", kfold = TRUE))
  expect_false(is_known_setup_err(msg_pw, mod_nm = "glmm"))
  expect_false(is_known_setup_err("other error", mod_nm = "gamm"))
  expect_false(is_known_setup_err(NULL, mod_nm = "gamm"))
  expect_false(is_known_setup_err(msg_bs, mod_nm = NULL))
})

test_that("setup.R reports the dropped test setups", {
  # setup.R drops a test setup only after a known upstream error (see
  # is_known_setup_err()). The skip message lists the dropped test setups.
  if (length(setup_drops) > 0) {
    skip(paste0(
      "setup.R dropped test setups after a known upstream error: ",
      paste(names(setup_drops),
            vapply(setup_drops, paste, character(1), collapse = ", "),
            sep = ": ", collapse = "; ")
    ))
  }
  expect_length(setup_drops, 0)
})
