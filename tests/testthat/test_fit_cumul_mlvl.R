# The optimizer of the ordinal::clmm() submodel fits and their convergence
# check. These tests need no reference model fits, so they also run on CRAN.

clmm_test_data <- function() {
  withr::local_seed(5678)
  dat <- data.frame(x = rnorm(150), g = gl(15, 10))
  eta <- 0.6 * dat$x + rnorm(15, sd = 0.7)[dat$g] + rlogis(150)
  dat$y <- factor(cut(eta, c(-Inf, -1, 0, 1, Inf), labels = FALSE),
                  levels = 1:4, ordered = TRUE)
  dat
}

fit_clmm_test <- function(dat, formula = y ~ x + (1 | g), ...) {
  fit_cumul_mlvl(formula = formula, data = dat,
                 family = list(link = "logit"), weights = rep(1, nrow(dat)),
                 ...)
}

test_that("fit_cumul_mlvl() uses ucminf() with forward differences", {
  skip_if_not_installed("ordinal")
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

test_that("fit_cumul_mlvl() uses nlminb() if clmm() cannot use ucminf()", {
  skip_if_not_installed("ordinal")
  # clmm() can use ucminf() only for a single `(1 | <group>)` term:
  expect_no_warning(
    fit <- fit_clmm_test(clmm_test_data(), formula = y ~ x + (x | g)),
    message = "ucminf"
  )
  expect_identical(fit$control$method, "nlminb")
  expect_false(is.null(fit$optRes$objective))
})

test_that("fit_cumul_mlvl() fits the requested link for cloglog and loglog", {
  skip_if_not_installed("ordinal")
  withr::local_seed(4567)
  n <- 400
  dat <- data.frame(x = rnorm(n), g = gl(20, 20))
  # Quantile functions of the latent variable for each link:
  qfuns <- list(cloglog = function(p) log(-log1p(-p)),
                loglog = function(p) -log(-log(p)))
  for (lk in names(qfuns)) {
    lk_other <- setdiff(names(qfuns), lk)
    lat <- dat$x + qfuns[[lk]](runif(n))
    dat$y <- factor(cut(lat, c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE),
                    levels = 1:4, ordered = TRUE)
    # The data have no group effects, so the fit is close to the clm() fit
    # with the same link, and not to the one with the other link:
    fit_clm <- ordinal::clm(y ~ x, data = dat, link = lk)
    fit_clm_other <- ordinal::clm(y ~ x, data = dat, link = lk_other)
    fit <- fit_cumul_mlvl(y ~ x + (1 | g), data = dat,
                          family = list(link = lk), weights = rep(1, n))
    expect_lt(fit$ST[[1]][1, 1], 0.3)
    expect_equal(fit$coefficients, fit_clm$coefficients, tolerance = 0.05,
                 info = lk)
    expect_lt(max(abs(fit$coefficients - fit_clm$coefficients)),
              max(abs(fit$coefficients - fit_clm_other$coefficients)) / 5)
    expect_true(check_conv_s(fit), info = lk)
    # Also with a `control` from the caller:
    fit_ctrl <- fit_cumul_mlvl(y ~ x + (1 | g), data = dat,
                               family = list(link = lk), weights = rep(1, n),
                               control = ordinal::clmm.control())
    expect_equal(fit_ctrl$coefficients, fit_clm$coefficients,
                 tolerance = 0.05, info = lk)
  }
})

test_that("fit_cumul_mlvl() refits with nlminb() if ucminf() does not converge", {
  skip_if_not_installed("ordinal")
  local_mocked_bindings(check_conv_s = function(fit_s) {
    !identical(fit_s$control$method, "ucminf")
  })
  fit <- fit_clmm_test(clmm_test_data())
  expect_identical(fit$control$method, "nlminb")
})

test_that("fit_with_fallback() keeps or drops the first fit and its conditions", {
  fit_first <- function() {
    warning("first warning")
    message("first message")
    1
  }
  fit_second <- function() 2
  # Accepted: the result and the conditions of the first fit.
  expect_message(
    expect_warning(
      res <- fit_with_fallback(fit_first, fit_second, function(x) TRUE),
      "first warning"
    ),
    "first message"
  )
  expect_identical(res, 1)
  # Not accepted: the second fit, without the conditions of the first.
  expect_no_condition(
    res <- fit_with_fallback(fit_first, fit_second, function(x) FALSE)
  )
  expect_identical(res, 2)
  # An error in the first fit:
  expect_no_condition(
    res <- fit_with_fallback(function() stop("failed"), fit_second,
                             function(x) TRUE)
  )
  expect_identical(res, 2)
})

test_that("check_conv_s() accepts ucminf() fits by code and gradient", {
  skip_if_not_installed("ordinal")
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
  # An nlminb() fit is judged by its own code, also if `method` is "ucminf"
  # (clmm() uses nlminb() for models that ucminf() cannot fit):
  fit_nlminb <- fit_clmm_test(clmm_test_data(), method = "nlminb")
  fit_nlminb$control$method <- "ucminf"
  expect_true(check_conv_s(fit_nlminb))
})
