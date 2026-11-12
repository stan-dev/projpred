# Tests of the automatic refits of fit_glmer_callback() and fit_gamm_callback()
# after errors of the fitter. A mock of the fitter throws the error, records
# the arguments of each call, and returns "fit" for the call that succeeds.
# Except for one check of lme4 itself, the tests do not fit models.

# A copy of `fun` (with the formals of `fun`) that records `control` and
# `nAGQ` of each call and throws the errors in `msgs`: the i-th call throws
# `msgs[i]`, and `NA` returns "fit". After the end of `msgs`, the last entry
# applies. The calls are in `<mock>$calls()`.
mock_fitter <- function(fun, msgs) {
  st <- new.env()
  st$calls <- list()
  mock <- fun
  body(mock) <- quote({
    env <- environment()
    nms <- intersect(c("control", "nAGQ"), names(formals(sys.function())))
    st$calls[[length(st$calls) + 1L]] <- lapply(stats::setNames(nm = nms),
                                                get, envir = env)
    msg <- msgs[min(length(st$calls), length(msgs))]
    if (!is.na(msg)) stop(msg)
    "fit"
  })
  mock_env <- new.env(parent = environment(fun))
  mock_env$st <- st
  mock_env$msgs <- msgs
  environment(mock) <- mock_env
  list(fun = mock, calls = function() st$calls)
}

dat_retry <- function() {
  withr::local_seed(1)
  n <- 40
  data.frame(y = rbinom(n, 1, 0.5), x = rnorm(n), g = gl(4, n / 4),
             h = gl(5, 1, n), f = gl(2, n / 2))
}

fit_glmer_retry <- function(formula, ...) {
  bars <- vapply(reformulas::findbars(formula), deparse1, "")
  random <- stats::as.formula(paste("~", paste0("(", bars, ")",
                                                collapse = " + ")))
  fit_glmer_callback(formula,
                     projpred_formula_no_random = reformulas::nobars(formula),
                     projpred_random = random, family = binomial(),
                     data = dat_retry(), ...)
}

msg_pirls <- "PIRLS step-halvings failed to reduce deviance in pwrssUpdate"
msg_pwrss <- "pwrssUpdate did not converge in (maxit) iterations"

test_that(paste(
  "fit_glmer_callback() refits with optimx and `nlminb` after a `not",
  "positive definite` error and stops if that refit fails too"
), {
  skip_if_not_installed("optimx")
  msg_pd <- "Downdated VtV is not positive definite"
  m <- mock_fitter(lme4::glmer, c(msg_pd, NA))
  testthat::local_mocked_bindings(glmer = m$fun, .package = "lme4")
  expect_identical(fit_glmer_retry(y ~ x + (1 | g)), "fit")
  calls <- m$calls()
  expect_length(calls, 2L)
  expect_identical(calls[[1]]$control$optimizer, c("bobyqa", "bobyqa"))
  expect_identical(calls[[2]]$control$optimizer, c("optimx", "optimx"))
  expect_identical(calls[[2]]$control$optCtrl$method, "nlminb")

  m <- mock_fitter(lme4::glmer, msg_pd)
  testthat::local_mocked_bindings(glmer = m$fun, .package = "lme4")
  expect_error(fit_glmer_retry(y ~ x + (1 | g)),
               "cannot fix this automatically anymore")
  expect_length(m$calls(), 2L)
})

test_that(paste(
  "fit_glmer_callback() refits with `nAGQ` from 20 to 30 after PIRLS errors",
  "for a single scalar group-level term"
), {
  m <- mock_fitter(lme4::glmer, msg_pirls)
  testthat::local_mocked_bindings(glmer = m$fun, .package = "lme4")
  expect_error(fit_glmer_retry(y ~ x + (1 | g)),
               "Encountering a PIRLS error .* cannot fix this automatically")
  nAGQ_calls <- vapply(m$calls(), function(cl) as.integer(cl$nAGQ), 1L)
  expect_identical(nAGQ_calls, c(1L, 20:30))
  # A slope without an intercept of a numeric variable is scalar, too:
  m <- mock_fitter(lme4::glmer, c(msg_pirls, NA))
  testthat::local_mocked_bindings(glmer = m$fun, .package = "lme4")
  expect_identical(fit_glmer_retry(y ~ x + (0 + x | g)), "fit")
  expect_identical(m$calls()[[2]]$nAGQ, 20L)
})

test_that(paste(
  "fit_glmer_callback() stops with the lme4 error after a PIRLS error if the",
  "submodel has group-level slopes or more than one group-level term"
), {
  for (fml in list(y ~ x + (x | g), y ~ x + (1 | g) + (1 | h),
                   y ~ x + (0 + f | g))) {
    m <- mock_fitter(lme4::glmer, msg_pirls)
    testthat::local_mocked_bindings(glmer = m$fun, .package = "lme4")
    expect_error(fit_glmer_retry(fml),
                 paste0("requires a single scalar group-level term.*",
                        "The lme4 error was: PIRLS step-halvings"),
                 info = deparse(fml))
    expect_length(m$calls(), 1L)
  }
})

test_that("glmer_scalar_re() counts the model-matrix columns of the term", {
  d <- dat_retry()
  expect_true(glmer_scalar_re(y ~ x + (1 | g), data = d))
  expect_true(glmer_scalar_re(y ~ x + (0 + x | g), data = d))
  expect_false(glmer_scalar_re(y ~ x + (x | g), data = d))
  expect_false(glmer_scalar_re(y ~ x + (0 + f | g), data = d))
  expect_false(glmer_scalar_re(y ~ x + (1 | g) + (1 | h), data = d))
  expect_false(glmer_scalar_re(y ~ x + (1 + x || g), data = d))
  expect_false(glmer_scalar_re(y ~ x, data = d))
  # A variable that is not found gives `TRUE` (then the refit runs):
  expect_true(glmer_scalar_re(y ~ x + (0 + x_missing | g), data = d))
  # lme4 itself rejects `nAGQ > 1` for a case with one grouping factor:
  fml_two <- y ~ x + (1 | g) + (0 + x | g)
  expect_false(glmer_scalar_re(fml_two, data = d))
  expect_error(suppressWarnings(lme4::glmer(fml_two, data = d,
                                            family = binomial(), nAGQ = 2L)),
               "single.*scalar")
})

test_that(paste(
  "fit_glmer_callback() refits with a larger `tolPwrss` and more optimizer",
  "evaluations after `pwrssUpdate did not converge` errors, then stops"
), {
  m <- mock_fitter(lme4::glmer, msg_pwrss)
  testthat::local_mocked_bindings(glmer = m$fun, .package = "lme4")
  expect_error(fit_glmer_retry(y ~ x + (1 | g)),
               "pwrssUpdate did not converge .* cannot fix this automatically")
  calls <- m$calls()
  tol <- vapply(calls, function(cl) cl$control$tolPwrss, 1)
  expect_equal(tol, 10^(-7:-3))
  maxfun <- vapply(calls[-1], function(cl) cl$control$optCtrl$maxfun, 1)
  expect_equal(maxfun, 10^(4:7))
})

test_that(paste(
  "fit_glmer_callback() with `projpred.PQL = TRUE` refits once with other",
  "MASS::glmmPQL() controls after each of its three known errors"
), {
  # control_callback() calls nlme::lmeControl() with `projpred.PQL = TRUE`:
  skip_if_not_installed("nlme")
  withr::local_options(projpred.PQL = TRUE)
  ctrl_cases <- list(
    list(msg = "iteration limit reached without convergence (10)",
         ctrl = "msMaxIter", val = 100),
    list(msg = "false convergence (8)", ctrl = "niterEM", val = 50),
    list(msg = "fewer observations than random effects in all level 1 groups",
         ctrl = "allow.n.lt.q", val = TRUE)
  )
  for (cs in ctrl_cases) {
    m <- mock_fitter(MASS::glmmPQL, c(cs$msg, NA))
    testthat::local_mocked_bindings(glmmPQL = m$fun, .package = "MASS")
    expect_identical(fit_glmer_retry(y ~ x + (1 | g)), "fit", info = cs$ctrl)
    calls <- m$calls()
    expect_length(calls, 2L)
    expect_equal(calls[[2]]$control[[cs$ctrl]], cs$val, info = cs$ctrl)

    m <- mock_fitter(MASS::glmmPQL, cs$msg)
    testthat::local_mocked_bindings(glmmPQL = m$fun, .package = "MASS")
    expect_error(fit_glmer_retry(y ~ x + (1 | g)),
                 "cannot fix this automatically anymore", info = cs$ctrl)
    expect_length(m$calls(), 2L)
  }
})

fit_gamm_retry <- function() {
  fit_gamm_callback(
    y ~ x + s(x) + (1 | g), projpred_formula_no_random = y ~ x + s(x),
    projpred_random = ~ (1 | g), data = dat_retry(), family = binomial()
  )
}

test_that(paste(
  "fit_gamm_callback() refits with optimx and `nlminb` after a `not",
  "positive definite` error and stops if that refit fails too"
), {
  skip_if_not_installed("optimx")
  m <- mock_fitter(gamm4::gamm4, "Downdated VtV is not positive definite")
  testthat::local_mocked_bindings(gamm4 = m$fun, .package = "gamm4")
  expect_error(fit_gamm_retry(), "cannot fix this automatically anymore")
  calls <- m$calls()
  expect_length(calls, 2L)
  expect_identical(calls[[2]]$control$optimizer, c("optimx", "optimx"))
  expect_identical(calls[[2]]$control$optCtrl$method, "nlminb")
})

test_that(paste(
  "fit_gamm_callback() refits with a larger `tolPwrss` after `pwrssUpdate`",
  "errors, then stops"
), {
  m <- mock_fitter(gamm4::gamm4, msg_pirls)
  testthat::local_mocked_bindings(gamm4 = m$fun, .package = "gamm4")
  expect_error(fit_gamm_retry(), "pwrssUpdate error .* cannot fix this")
  tol <- vapply(m$calls(), function(cl) cl$control$tolPwrss, 1)
  expect_equal(tol, 10^(-7:-4))
})
