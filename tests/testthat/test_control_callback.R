# control_callback() and the optimizer of lme4::glmer() submodel fits. These
# tests need no reference model fits, so they also run on CRAN.

test_that("control_callback() uses BOBYQA in both stages of glmer()", {
  skip_if_not_installed("lme4")
  expect_identical(control_callback(binomial())$optimizer, c("bobyqa", "bobyqa"))
  # Arguments of the refits after lme4 errors:
  ctrl_tol <- control_callback(poisson(), tolPwrss = 1e-6,
                               optCtrl = list(maxfun = 1e4))
  expect_identical(ctrl_tol$optimizer, c("bobyqa", "bobyqa"))
  expect_identical(ctrl_tol$tolPwrss, 1e-6)
  expect_identical(ctrl_tol$optCtrl$maxfun, 1e4)
  ctrl_optimx <- control_callback(binomial(), optimizer = "optimx",
                                  optCtrl = list(method = "nlminb"))
  expect_identical(ctrl_optimx$optimizer, c("optimx", "optimx"))
  expect_s3_class(control_callback(gaussian()), "lmerControl")
})

test_that("control_callback() is unchanged for the PQL fitter", {
  skip_if_not_installed("nlme")
  opts_old <- options(projpred.PQL = TRUE)
  on.exit(options(opts_old), add = TRUE)
  ctrl <- control_callback(binomial())
  expect_false(inherits(ctrl, "merControl"))
  expect_identical(ctrl, nlme::lmeControl())
})

test_that("fit_glmer_callback() fits with BOBYQA by default", {
  skip_if_not_installed("lme4")
  set.seed(5678)
  dat <- data.frame(x = rnorm(120), g = gl(12, 10))
  dat$y <- rbinom(120, 1, plogis(0.5 * dat$x + rnorm(12, sd = 0.8)[dat$g]))
  fit <- fit_glmer_callback(
    formula = y ~ x + (1 | g),
    projpred_formula_no_random = y ~ x,
    projpred_random = ~ (1 | g),
    family = binomial(),
    data = dat
  )
  expect_s4_class(fit, "glmerMod")
  expect_identical(fit@optinfo$optimizer, "bobyqa")
  # A `control` from the caller replaces the default:
  fit_nm <- fit_glmer_callback(
    formula = y ~ x + (1 | g),
    projpred_formula_no_random = y ~ x,
    projpred_random = ~ (1 | g),
    family = binomial(),
    control = lme4::glmerControl(),
    data = dat
  )
  expect_identical(fit_nm@optinfo$optimizer, "Nelder_Mead")
})
