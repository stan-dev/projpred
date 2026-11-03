# The optimizer of the ordinal::clmm() submodel fits and their convergence
# check. These tests need no reference model fits, so they also run on CRAN.

clmm_test_data <- function() {
  set.seed(5678)
  dat <- data.frame(x = rnorm(150), g = gl(15, 10))
  eta <- 0.6 * dat$x + rnorm(15, sd = 0.7)[dat$g] + rlogis(150)
  dat$y <- factor(cut(eta, c(-Inf, -1, 0, 1, Inf), labels = FALSE),
                  levels = 1:4, ordered = TRUE)
  dat
}

fit_clmm_test <- function(dat, ...) {
  fit_cumul_mlvl(formula = y ~ x + (1 | g), data = dat,
                 family = list(link = "logit"), weights = rep(1, nrow(dat)),
                 ...)
}

test_that("fit_cumul_mlvl() uses ucminf() with forward differences", {
  skip_if_not_installed("ordinal")
  skip_if_not_installed("ucminf")
  dat <- clmm_test_data()
  fit <- fit_clmm_test(dat)
  expect_s3_class(fit, "clmm")
  expect_identical(fit$control$method, "ucminf")
  expect_identical(fit$control$optCtrl$grad, "forward")
  expect_true(check_conv_s(fit))
  # A `method` or `control` from the caller replaces the default:
  fit_nlminb <- fit_clmm_test(dat, method = "nlminb")
  expect_identical(fit_nlminb$control$method, "nlminb")
  expect_true(check_conv_s(fit_nlminb))
  fit_ctrl <- fit_clmm_test(dat, control = ordinal::clmm.control())
  expect_identical(fit_ctrl$control$method, "nlminb")
  # Both optimizers give the same fit within the optimizer tolerance:
  expect_equal(fit$logLik, fit_nlminb$logLik, tolerance = 1e-6)
  expect_equal(fit$coefficients, fit_nlminb$coefficients, tolerance = 1e-3)
})

test_that("fit_cumul_mlvl() refits with nlminb() if ucminf() does not converge", {
  skip_if_not_installed("ordinal")
  skip_if_not_installed("ucminf")
  local_mocked_bindings(check_conv_s = function(fit_s) {
    !identical(fit_s$control$method, "ucminf")
  })
  fit <- fit_clmm_test(clmm_test_data())
  expect_identical(fit$control$method, "nlminb")
})

test_that("check_conv_s() accepts ucminf() fits by code and gradient", {
  skip_if_not_installed("ordinal")
  skip_if_not_installed("ucminf")
  fit <- fit_clmm_test(clmm_test_data())
  fit$gradient[] <- 1e-5
  for (code in c(1, 2, 4)) {
    fit$optRes$convergence <- code
    expect_true(check_conv_s(fit), info = paste("code", code))
  }
  # Evaluation limit or no start:
  for (code in c(3, -2)) {
    fit$optRes$convergence <- code
    expect_false(check_conv_s(fit), info = paste("code", code))
  }
  # A zero step far from the optimum:
  fit$optRes$convergence <- 4
  fit$gradient[1] <- 0.5
  expect_false(check_conv_s(fit))
  fit$gradient[1] <- NA
  expect_false(check_conv_s(fit))
})
