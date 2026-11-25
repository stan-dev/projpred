# These tests run on CRAN on purpose, so that the reverse dependency checks of
# lme4, MASS, nlme, mgcv, gamm4, ordinal and mclogit see changes in the
# submodel fits that projpred needs. They call divmin() and divmin_augdat()
# directly with small simulated data (no reference model). They test
# properties of the projection (known solutions, observation weights), not
# the equality with a direct call of the same fitter.

# Data with `n_grp` groups of `n_per` observations. The linear predictor `eta`
# has the group-level effects with standard deviation `sd_grp`.
# The state of the random number generator is restored at the end.
dat_fitters <- function(n_grp = 6, n_per = 8, sd_grp = 0.8, seed = 1) {
  withr::local_seed(seed)
  n <- n_grp * n_per
  g <- factor(paste0("g", rep(seq_len(n_grp), each = n_per)))
  x <- rnorm(n)
  z <- rnorm(n)
  eta <- -0.2 + 0.7 * x + sin(2 * z) + rnorm(n_grp, sd = sd_grp)[g]
  data.frame(x = x, z = z, g = g, eta = eta)
}

# Calls divmin() as proj_to_submodl() does: `mu` has one column per projected
# draw, and the response columns are `.y.1`, `.y.2`, ...
divmin_tst <- function(rhs, data, family, mu, weights = rep(1, nrow(data)),
                       throw_warn = FALSE, check_conv = FALSE, ...) {
  mu <- as.matrix(mu)
  y_nms <- paste0(".y.", seq_len(ncol(mu)))
  data[y_nms] <- mu
  lhs <- if (ncol(mu) > 1) {
    paste0("cbind(", paste(y_nms, collapse = ", "), ")")
  } else {
    y_nms
  }
  divmin(formula = as.formula(paste(lhs, "~", rhs)), data = data,
         family = family, weights = weights,
         projpred_var = matrix(NA_real_, nrow(data), ncol(mu)),
         projpred_ws_aug = mu, verbose_divmin = FALSE,
         throw_warn_sdivmin = throw_warn, do_check_conv = check_conv, ...)
}

# Calls divmin_augdat() with one projected draw: `probs` has one row per
# observation and one column per response category (the augmented-data
# weights).
divmin_augdat_tst <- function(rhs, data, family, probs, throw_warn = FALSE,
                              check_conv = FALSE, ...) {
  data_aug <- do.call(rbind, lapply(family$cats, function(cat_i) {
    data$y <- cat_i
    data
  }))
  data_aug$y <- factor(data_aug$y, levels = family$cats)
  divmin_augdat(formula = as.formula(paste("y ~", rhs)), data = data_aug,
                family = family, weights = "unused",
                projpred_var = matrix(NA_real_, nrow(data_aug), 1),
                projpred_ws_aug = matrix(c(probs)), verbose_divmin = FALSE,
                throw_warn_sdivmin = throw_warn, do_check_conv = check_conv,
                ...)
}

# Inside test_that(), testthat handles messages and warnings before R prints
# them, so capt_mssgs_warns() would capture nothing. While a message sink is
# active, with_r_printing() prints them as R does with `options(warn = 1)`, so
# that the filters of divmin() and divmin_augdat() get the usual lines (R puts
# long messages on a second line; the filters accept both forms).
with_r_printing <- function(expr) {
  is_captured <- function() sink.number(type = "message") != 2L
  withCallingHandlers(
    expr,
    message = function(m) {
      if (is_captured()) {
        cat(conditionMessage(m), file = stderr())
        invokeRestart("muffleMessage")
      }
    },
    warning = function(w) {
      if (is_captured()) {
        call_w <- conditionCall(w)
        head_w <- if (is.null(call_w)) {
          "Warning:"
        } else {
          sprintf("Warning in %s :", deparse(call_w)[1])
        }
        cat(sprintf("%s %s\n", head_w, conditionMessage(w)), file = stderr())
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Observations 1 to 3 get weight 0 and other predictor values:
idxs_w0 <- 1:3
change_w0 <- function(data) {
  data$x[idxs_w0] <- data$x[idxs_w0] + 5
  data$z[idxs_w0] <- data$z[idxs_w0] - 3
  data
}
# Observations 1 to 10 get weight 2, or are duplicated:
idxs_w2 <- 1:10

pars_mer <- function(fit) {
  c(lme4::fixef(fit), lme4::getME(fit, "theta"))
}

# lme4::glmer() -----------------------------------------------------------

# 39 observations in 6 groups. The singular fit (theta = 0) equals the glm()
# solution. The rounded responses separate the data, so a fit to rounded data
# has fixed effects of about +-30 to +-40.
dat_binom_singular <- function() {
  data.frame(
    y = c(0.2565, 0.3552, 0.344, 0.3393, 0.2216, 0.2628, 0.479, 0.2689,
          0.2791, 0.1885, 0.3799, 0.4617, 0.2209, 0.4506, 0.4922, 0.5739,
          0.8292, 0.6354, 0.776, 0.048, 0.5312, 0.0616, 0.1741, 0.9788, 0.917,
          0.3709, 0.5611, 0.111, 0.0479, 0.3284, 0.1255, 0.1129, 0.1322,
          0.2209, 0.1343, 0.1795, 0.0932, 0.1599, 0.3264),
    x = c(0.7458, 0.2496, 0.3323, 0.116, 2.0572, 1.0248, -0.6979, 1.0392,
          1.1827, 2.0272, 0.3178, -0.0954, -2.0399, -1.7392, -1.6437, -0.8398,
          2.3697, 0.0016, 0.1374, -0.769, 0.0934, -0.7415, -0.4593, 1.0565,
          1.3049, 0.1363, 0.9627, 0.0119, -0.3273, 0.6624, -2.0192, -1.7293,
          -0.9185, 0.9509, -1.1285, -0.2256, 1.1649, -0.3001, -0.6718),
    f = factor(rep(c("a", "b", "c"), c(13, 13, 13))),
    g = factor(paste0("g", c(rep(1:6, each = 6), 1, 1, 1)))
  )
}

# The form of the augmented-data projection: each observation gives a row with
# response 1 and weight y, and a row with response 0 and weight 1 - y.
to_augdat_form <- function(d) {
  d1 <- d
  d0 <- d
  d1$w_aug <- d$y
  d0$w_aug <- 1 - d$y
  d1$y <- 1
  d0$y <- 0
  rbind(d1, d0)
}

# The glmer() code path (`mlvl_laplace = FALSE`: no Laplace fit of projpred).
fit_glmer_test <- function(formula, data, family, weights, ...) {
  suppressWarnings(suppressMessages(fit_glmer_callback(
    formula, projpred_formula_no_random = NA, projpred_random = NA,
    family = family, data = data, weights = weights, mlvl_laplace = FALSE, ...
  )))
}

# The `nAGQ` value of a glmer() fit.
nAGQ_of <- function(fit) {
  unname(fit@devcomp$dims[["nAGQ"]])
}

test_that(paste(
  "a binomial glmer() submodel with non-integer data and the default `nAGQ`",
  "is not fitted to rounded data (all lme4 versions)"
), {
  d <- dat_binom_singular()
  fam <- extend_family(binomial())
  coef_glm <- coef(suppressWarnings(glm(y ~ x + f, data = d,
                                        family = binomial())))
  # Traditional projection: non-integer responses.
  fit <- fit_glmer_test(y ~ x + f + (1 | g), d, fam,
                        weights = rep(1, nrow(d)))
  expect_equal(unname(lme4::fixef(fit)), unname(coef_glm), tolerance = 1e-3)
  expect_identical(nAGQ_of(fit), 1L)
  # Augmented-data projection: 0/1 responses with non-integer weights.
  d_aug <- to_augdat_form(d)
  fit_aug <- fit_glmer_test(y ~ x + f + (1 | g), d_aug, fam,
                            weights = d_aug$w_aug)
  expect_equal(unname(lme4::fixef(fit_aug)), unname(coef_glm),
               tolerance = 1e-3)
  expect_identical(nAGQ_of(fit_aug), 1L)
})

test_that(paste(
  "the refit of a binomial glmer() submodel after a PIRLS error uses",
  "`nAGQ = 20`, or the `nAGQ` from `...` plus 1"
), {
  d <- dat_fitters()
  d$y <- plogis(d$eta)
  testthat::local_mocked_bindings(lme4_needs_glmer_family = function() TRUE)
  # A glmer() that records `nAGQ` and the family name and fails in its first
  # call:
  nAGQ_calls <- list()
  fam_calls <- character()
  glmer_mock <- lme4::glmer
  body(glmer_mock) <- quote({
    nAGQ_calls[[length(nAGQ_calls) + 1L]] <<- nAGQ
    fam_calls <<- c(fam_calls, family$family)
    if (length(nAGQ_calls) == 1L) {
      stop("PIRLS step-halvings failed to reduce deviance")
    }
    "fit"
  })
  environment(glmer_mock) <- environment()
  testthat::local_mocked_bindings(glmer = glmer_mock, .package = "lme4")
  fit <- fit_glmer_test(y ~ x + (1 | g), d, extend_family(binomial()),
                        weights = rep(1, nrow(d)))
  expect_identical(fit, "fit")
  expect_identical(nAGQ_calls, list(1L, 20L))
  # The refit also gets the family from glmer_family():
  expect_identical(fam_calls, rep("binomial_projpred", 2))
  nAGQ_calls <- list()
  fit <- fit_glmer_test(y ~ x + (1 | g), d, extend_family(binomial()),
                        weights = rep(1, nrow(d)), nAGQ = 3L)
  expect_identical(fit, "fit")
  expect_identical(nAGQ_calls, list(3L, 4L))
})

test_that(paste(
  "a binomial glmer() submodel with fractional responses gives the glm()",
  "solution for a singular fit"
), {
  d <- dat_binom_singular()
  coef_glm <- coef(suppressWarnings(glm(y ~ x + f, data = d,
                                        family = binomial())))
  for (nAGQ in c(1L, 20L)) {
    fit <- fit_glmer_test(y ~ x + f + (1 | g), d, extend_family(binomial()),
                          weights = rep(1, nrow(d)), nAGQ = nAGQ)
    expect_equal(unname(lme4::fixef(fit)), unname(coef_glm), tolerance = 1e-3,
                 info = paste("nAGQ", nAGQ))
  }
})

test_that(paste(
  "the traditional and the augmented-data form of a binomial glmer()",
  "submodel give the same fit"
), {
  # A singular fit and a fit with theta > 0 (10 groups of 30 observations):
  withr::local_seed(11)
  g <- factor(rep(1:10, each = 30))
  x <- rnorm(300)
  d_theta <- data.frame(y = plogis(-0.3 + 0.8 * x + rnorm(10, sd = 1.5)[g]),
                        x = x, g = g)
  cases <- list(
    singular = list(d = dat_binom_singular(), fml = y ~ x + f + (1 | g)),
    theta = list(d = d_theta, fml = y ~ x + (1 | g))
  )
  for (nm in names(cases)) {
    d <- cases[[nm]]$d
    d_aug <- to_augdat_form(d)
    fit_trad <- fit_glmer_test(cases[[nm]]$fml, d, extend_family(binomial()),
                               weights = rep(1, nrow(d)))
    fit_aug <- fit_glmer_test(cases[[nm]]$fml, d_aug,
                              extend_family(binomial()),
                              weights = d_aug$w_aug)
    expect_equal(lme4::fixef(fit_aug), lme4::fixef(fit_trad),
                 tolerance = 1e-4, info = nm)
    expect_equal(lme4::getME(fit_aug, "theta"),
                 lme4::getME(fit_trad, "theta"), tolerance = 1e-3, info = nm)
    if (nm == "theta") {
      # The fit without rounding has slope 0.780 and theta 1.73; the fit to
      # rounded data has slope 1.69 and theta 3.14:
      expect_equal(unname(lme4::fixef(fit_trad))[2], 0.780, tolerance = 0.01)
      expect_equal(unname(lme4::getME(fit_trad, "theta")), 1.73,
                   tolerance = 0.03)
    }
  }
})

test_that(paste(
  "a Poisson glmer() submodel with non-integer responses estimates the",
  "group-level standard deviation"
), {
  withr::local_seed(4)
  g <- factor(rep(seq_len(8), each = 6))
  x <- rnorm(48)
  b <- rnorm(8, sd = 0.8)
  d <- data.frame(y = exp(0.5 + 0.4 * x + b[g] + rnorm(48, sd = 0.1)),
                  x = x, g = g)
  fit <- fit_glmer_test(y ~ x + (1 | g), d, extend_family(poisson()),
                        weights = rep(1, nrow(d)))
  # With an infinite objective, the optimizer does not move, and theta stays
  # at its start value 1:
  theta <- unname(lme4::getME(fit, "theta"))
  expect_gt(abs(theta - 1), 1e-4)
  expect_gt(theta, 0.2)
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that(paste(
  "glmer_family() renames the binomial and the Poisson family and keeps",
  "their aic for integer data"
), {
  testthat::local_mocked_bindings(lme4_needs_glmer_family = function() TRUE)
  withr::local_seed(5309)
  mu <- plogis(rnorm(20))
  ntr <- rpois(20, 4) + 1
  w1 <- rep(1, 20)
  fam_b <- glmer_family(extend_family(binomial()))
  expect_identical(fam_b$family, "binomial_projpred")
  expect_identical(fam_b$dispersion, 1)
  # Proportions (with the numbers of trials as weights or as `n`), 0/1
  # responses, and zero counts at the bounds of `mu`:
  y <- rbinom(20, ntr, mu) / ntr
  expect_equal(fam_b$aic(y, w1, mu, ntr, 0), binomial()$aic(y, w1, mu, ntr, 0),
               tolerance = 1e-12)
  expect_equal(fam_b$aic(y, ntr, mu, ntr, 0),
               binomial()$aic(y, ntr, mu, ntr, 0), tolerance = 1e-12)
  y01 <- rbinom(20, 1, mu)
  expect_equal(fam_b$aic(y01, w1, mu, w1, 0),
               binomial()$aic(y01, w1, mu, w1, 0), tolerance = 1e-12)
  expect_equal(fam_b$aic(c(0, 1), c(1, 1), c(0, 1), c(1, 1), 0), 0)
  lam <- exp(rnorm(20))
  ycnt <- rpois(20, lam)
  fam_p <- glmer_family(extend_family(poisson()))
  expect_identical(fam_p$family, "poisson_projpred")
  expect_equal(fam_p$aic(ycnt, w1, lam, ntr, 0),
               poisson()$aic(ycnt, w1, lam, ntr, 0), tolerance = 1e-12)
  # Non-integer responses:
  expect_true(is.finite(fam_p$aic(lam, w1, lam, w1, 0)))
  expect_true(is.finite(fam_b$aic(mu, w1, mu, ntr, 0)))
  # Other families, and the other lme4 versions:
  expect_identical(glmer_family(extend_family(gaussian()))$family, "gaussian")
  testthat::with_mocked_bindings(
    expect_identical(glmer_family(extend_family(poisson()))$family, "poisson"),
    lme4_needs_glmer_family = function() FALSE
  )
  # The control (lme4 >= 2.1-0 has `disp_dof_correction`) and the data:
  has_disp_dof <- "disp_dof_correction" %in%
    names(formals(lme4::glmerControl))
  ctrl <- glmer_control(lme4::glmerControl(), family = fam_b)
  expect_identical(isFALSE(ctrl$disp_dof_correction), has_disp_dof)
  ctrl <- glmer_control(list(optimizer = "bobyqa"), family = fam_p)
  expect_identical(names(ctrl),
                   c("optimizer", if (has_disp_dof) "disp_dof_correction"))
  expect_identical(glmer_control(lme4::glmerControl(), family = binomial()),
                   lme4::glmerControl())
  d <- data.frame(y = factor(c("a", "b", "b")), x = 1:3)
  expect_identical(glmer_data(d, y ~ x, family = fam_b)$y, c(0, 1, 1))
  expect_identical(glmer_data(d, relevel(y, "b") ~ x, family = fam_b), d)
  expect_identical(glmer_data(d, y ~ x, family = binomial()), d)
  expect_identical(glmer_data(d, y ~ x, family = fam_p), d)
  d$y <- c(TRUE, FALSE, TRUE)
  expect_identical(glmer_data(d, y ~ x, family = fam_b), d)
})

test_that(paste(
  "the family from glmer_family() gives the fit of lme4's family for integer",
  "data"
), {
  testthat::local_mocked_bindings(lme4_needs_glmer_family = function() TRUE)
  d <- dat_fitters()
  w1 <- rep(1, nrow(d))
  withr::local_seed(2)
  d$y_poisson <- rpois(nrow(d), exp(d$eta))
  d$y_binomial <- rbinom(nrow(d), 1, plogis(d$eta))
  for (nm in c("poisson", "binomial")) {
    fml <- as.formula(paste0("y_", nm, " ~ x + (1 | g)"))
    fam <- extend_family(get(nm)())
    fit <- fit_glmer_test(fml, d, fam, weights = w1)
    expect_identical(fit@resp$family$family, paste0(nm, "_projpred"))
    fit_lme4 <- suppressMessages(lme4::glmer(
      fml, data = d, family = get(nm)(),
      control = lme4::glmerControl(optimizer = "bobyqa")
    ))
    expect_equal(pars_mer(fit), pars_mer(fit_lme4), tolerance = 1e-4,
                 info = nm)
    expect_equal(unname(lme4::VarCorr(fit)$g[1]),
                 unname(lme4::VarCorr(fit_lme4)$g[1]), tolerance = 1e-4,
                 info = nm)
  }
})

test_that(paste(
  "a Poisson glmer() submodel with non-integer responses and without",
  "group-level variation gives the glm() solution"
), {
  withr::local_seed(4)
  g <- factor(rep(seq_len(8), each = 6))
  x <- rnorm(48)
  # (lme4 gives a PIRLS error for an exact fit, so the data have some noise.)
  d <- data.frame(y = exp(0.5 + 0.4 * x + rnorm(48, sd = 0.01)), x = x, g = g)
  coef_glm <- coef(suppressWarnings(glm(y ~ x, data = d, family = poisson())))
  fit <- fit_glmer_test(y ~ x + (1 | g), d, extend_family(poisson()),
                        weights = rep(1, nrow(d)))
  # A singular fit (with lme4's Poisson objective, which is infinite for these
  # data, theta stays at its start value 1):
  expect_lt(unname(lme4::getME(fit, "theta")), 1e-3)
  expect_equal(unname(lme4::fixef(fit)), unname(coef_glm), tolerance = 1e-4)
})

test_that(paste(
  "Poisson glmer() and gamm4() submodels with non-integer responses give no",
  "convergence warning"
), {
  d <- dat_fitters()
  mu <- exp(cbind(d$eta, d$eta + 0.3 * d$z))
  for (rhs in c("x + (1 | g)", "s(z) + x + (1 | g)")) {
    expect_no_warning(
      fits <- with_r_printing(
        divmin_tst(rhs, d, extend_family(poisson()), mu, throw_warn = TRUE,
                   check_conv = TRUE)
      ),
      message = "converge"
    )
    expect_length(fits, 2)
    if (startsWith(rhs, "s(")) {
      # The gamm4() fits are not singular, so gamm_check_conv() checks the
      # gradient:
      for (fit in fits) {
        expect_false(is.null(fit$mer@optinfo$derivs))
      }
    }
  }
})

test_that(paste(
  "a binomial gamm4() submodel with non-integer responses gives no",
  "convergence warning"
), {
  # With fewer data, the fit is singular and has no gradient to check:
  d <- dat_fitters(n_grp = 10, n_per = 20)
  mu <- plogis(cbind(d$eta, d$eta + 0.3 * d$z))
  expect_no_warning(
    fits <- with_r_printing(
      divmin_tst("s(z) + x + (1 | g)", d, extend_family(binomial()), mu,
                 weights = rep(20, nrow(d)), throw_warn = TRUE,
                 check_conv = TRUE)
    ),
    message = "converge"
  )
  expect_length(fits, 2)
  # The fits are not singular, so gamm_check_conv() checks the gradient:
  for (fit in fits) {
    expect_false(is.null(fit$mer@optinfo$derivs))
  }
})

test_that(paste(
  "a Poisson gamm4() submodel with non-integer responses estimates the",
  "group-level standard deviation"
), {
  d <- dat_fitters()
  d$y <- exp(d$eta)
  fam <- extend_family(poisson())
  w1 <- rep(1, nrow(d))
  fit_gamm <- suppressWarnings(suppressMessages(fit_gamm_callback(
    formula = y ~ x + (1 | g), projpred_formula_no_random = y ~ x,
    projpred_random = ~ (1 | g), data = d, family = fam, weights = w1
  )))
  fit_glmer <- fit_glmer_test(y ~ x + (1 | g), d, fam, weights = w1)
  theta_gamm <- unname(lme4::getME(fit_gamm$mer, "theta"))
  # theta is estimated (with lme4's Poisson objective, which is infinite for
  # these data, theta stays at its start value 1):
  expect_gt(abs(theta_gamm - 1), 1e-3)
  expect_equal(theta_gamm, unname(lme4::getME(fit_glmer, "theta")),
               tolerance = 1e-2)
  expect_equal(unname(lme4::fixef(fit_gamm$mer)),
               unname(lme4::fixef(fit_glmer)), tolerance = 1e-3)
  # With a smooth term:
  fit_s <- suppressWarnings(
    divmin_tst("s(z) + x + (1 | g)", d, fam, d$y)
  )[[1]]
  theta_s <- lme4::getME(fit_s$mer, "theta")
  theta_g <- unname(theta_s[grep("^g\\.", names(theta_s))])
  expect_length(theta_g, 1)
  expect_gt(abs(theta_g - 1), 1e-3)
})

test_that(paste(
  "filter_mssgs_warns() removes the warnings about non-integer Poisson",
  "responses and keeps other messages"
), {
  # The lines as capt_mssgs_warns() captures them (with `options(warn = 1)`):
  mssgs_warns <- c(
    "Warning: non-integer response values in a poisson model",
    "Warning in dpois(y, mu, log = TRUE) : non-integer x = 0.152194",
    "Warning in fn(x) :",
    "boundary (singular) fit: see help('isSingular')",
    "Warning: Model failed to converge with max|grad| = 0.01"
  )
  expect_identical(
    filter_mssgs_warns(mssgs_warns),
    c("boundary (singular) fit: see help('isSingular')",
      "Warning: Model failed to converge with max|grad| = 0.01")
  )
})

test_that(paste(
  "with_r_printing() lets capt_mssgs_warns() capture messages and warnings,",
  "and passes on the other ones"
), {
  expect_identical(
    with_r_printing(capt_mssgs_warns({
      message("a message")
      warning("a warning", call. = FALSE)
    })),
    c("a message", "Warning: a warning")
  )
  expect_warning(with_r_printing(warning("not captured")), "^not captured$")
})

test_that(paste(
  "observations with weight 0 do not change a binomial glmer() submodel fit"
), {
  withr::local_options(projpred.mlvl_laplace = FALSE)
  d <- dat_fitters()
  mu <- plogis(d$eta)
  w <- rep(1, nrow(d))
  w[idxs_w0] <- 0
  for (nAGQ in c(1L, 10L)) {
    fit_w0 <- divmin_tst("x + (1 | g)", d, extend_family(binomial()), mu,
                         weights = w, nAGQ = nAGQ)[[1]]
    fit_chg <- divmin_tst("x + (1 | g)", change_w0(d),
                          extend_family(binomial()), mu, weights = w,
                          nAGQ = nAGQ)[[1]]
    fit_rm <- divmin_tst("x + (1 | g)", d[-idxs_w0, ],
                         extend_family(binomial()), mu[-idxs_w0],
                         nAGQ = nAGQ)[[1]]
    expect_equal(pars_mer(fit_chg), pars_mer(fit_w0), tolerance = 1e-6,
                 info = paste("nAGQ", nAGQ))
    expect_equal(pars_mer(fit_rm), pars_mer(fit_w0), tolerance = 1e-6,
                 info = paste("nAGQ", nAGQ))
  }
})

test_that(paste(
  "a glmer() submodel fit with integer weights equals the fit to the data",
  "with duplicated observations"
), {
  withr::local_options(projpred.mlvl_laplace = FALSE)
  d <- dat_fitters()
  idxs_dupl <- c(seq_len(nrow(d)), idxs_w2)
  w <- rep(1, nrow(d))
  w[idxs_w2] <- 2
  cases <- list(
    binomial_nAGQ1 = list(fam = binomial(), mu = plogis(d$eta), nAGQ = 1L),
    binomial_nAGQ10 = list(fam = binomial(), mu = plogis(d$eta), nAGQ = 10L),
    poisson_nAGQ1 = list(fam = poisson(), mu = exp(d$eta), nAGQ = 1L)
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    # For Poisson, divmin() removes the lme4 warning about non-integer
    # responses (the binomial fits are singular and give a message):
    is_pois <- cs$fam$family == "poisson"
    fit_w_fun <- function() {
      with_r_printing(
        divmin_tst("x + (1 | g)", d, extend_family(cs$fam), cs$mu,
                   weights = w, throw_warn = is_pois, nAGQ = cs$nAGQ)[[1]]
      )
    }
    if (is_pois) {
      expect_no_warning(fit_w <- fit_w_fun())
    } else {
      fit_w <- fit_w_fun()
    }
    fit_dupl <- suppressWarnings(
      divmin_tst("x + (1 | g)", d[idxs_dupl, ], extend_family(cs$fam),
                 cs$mu[idxs_dupl], nAGQ = cs$nAGQ)[[1]]
    )
    expect_equal(pars_mer(fit_dupl), pars_mer(fit_w), tolerance = 1e-4,
                 info = nm)
  }
})

# MASS::glmmPQL() ---------------------------------------------------------

test_that(paste(
  "a Gaussian glmmPQL() submodel fit equals the maximum-likelihood lmer() fit"
), {
  skip_if_not_installed("nlme")
  d <- dat_fitters()
  d$y_gauss <- d$eta
  fit_lmer <- lme4::lmer(y_gauss ~ x + (1 | g), data = d, REML = FALSE)
  fit_pql <- withr::with_options(
    list(projpred.PQL = TRUE),
    divmin_tst("x + (1 | g)", d, extend_family(gaussian()), d$eta)[[1]]
  )
  expect_s3_class(fit_pql, "glmmPQL")
  # The differences are about 1e-11:
  expect_equal(nlme::fixef(fit_pql), lme4::fixef(fit_lmer), tolerance = 1e-5)
  expect_equal(fit_pql$sigma, sigma(fit_lmer), tolerance = 1e-5)
})

test_that(paste(
  "glmmPQL() submodel fits give no warnings about messages or convergence"
), {
  skip_if_not_installed("nlme")
  d <- dat_fitters()
  withr::local_options(projpred.PQL = TRUE)
  for (fam in list(gaussian(), binomial(), poisson())) {
    mu <- fam$linkinv(cbind(d$eta, d$eta + 0.3 * d$z))
    expect_no_warning(
      fits <- with_r_printing(
        divmin_tst("x + (1 | g)", d, extend_family(fam), mu,
                   throw_warn = TRUE, check_conv = TRUE)
      )
    )
    expect_length(fits, 2)
  }
})

# mgcv::gam() -------------------------------------------------------------

test_that(paste(
  "gam_outer_conv() reads the result of each smoothing parameter optimizer"
), {
  expect_true(gam_outer_conv(NULL))
  expect_true(gam_outer_conv(list(conv = "full convergence", iter = 5)))
  expect_false(gam_outer_conv(list(conv = "step failed", iter = 5)))
  expect_false(gam_outer_conv(list(conv = "iteration limit reached")))
  # Only the iteration limits (nlm() code 4, optim() code 1) are
  # non-convergence:
  for (code in 1:5) {
    expect_identical(gam_outer_conv(list(minimum = 1, code = code)),
                     code != 4, info = paste("nlm() code", code))
  }
  for (code in c(0, 1, 51, 52)) {
    expect_identical(gam_outer_conv(list(par = 1, convergence = code)),
                     code != 1, info = paste("optim() code", code))
  }
})

# Fits of mgcv::gam() to the data of dat_fitters(). Tests that need a fit that
# did not converge skip if mgcv converges (its stopping rules decide).
fit_gam_tst <- function(formula, ...) {
  d <- dat_fitters()
  d$y_gauss <- d$eta
  d$y_binom <- plogis(d$eta)
  # A linear function of `x` and `z`: the smoothing parameters of s(z) and s(x)
  # go to infinity.
  d$y_lin <- -0.2 + 0.7 * d$x + 0.5 * d$z + 0.3 * sin(3 * seq_len(nrow(d)))
  suppressWarnings(mgcv::gam(formula, data = d, ...))
}

test_that(paste(
  "check_conv_s() accepts converged gam() submodel fits of each smoothing",
  "parameter optimizer"
), {
  fml_gauss <- y_gauss ~ s(z) + x
  # For this linear function, nlm() can stop with code 5:
  fml_lin <- y_lin ~ s(z, k = 5) + s(x, k = 5)
  expect_true(check_conv_s(fit_gam_tst(fml_gauss)))
  for (opt in c("newton", "bfgs", "nlm", "optim")) {
    expect_true(check_conv_s(fit_gam_tst(
      fml_gauss, method = "REML", optimizer = c("outer", opt)
    )), info = opt)
    expect_true(check_conv_s(fit_gam_tst(
      fml_lin, method = "REML", optimizer = c("outer", opt)
    )), info = paste(opt, "linear"))
  }
  expect_true(check_conv_s(fit_gam_tst(fml_gauss, method = "REML",
                                       optimizer = "efs")))
  expect_true(check_conv_s(fit_gam_tst(y_gauss ~ s(z, fx = TRUE) + x)))
})

test_that(paste(
  "check_conv_s() detects a gam() submodel fit at the iteration limit of nlm()"
), {
  fit_nlm <- fit_gam_tst(y_gauss ~ s(z) + x, method = "REML",
                         optimizer = c("outer", "nlm"),
                         control = mgcv::gam.control(nlm = list(iterlim = 1)))
  skip_if_not(identical(fit_nlm$outer.info$code, 4L),
              "nlm() converged in one iteration")
  expect_false(check_conv_s(fit_nlm))
})

test_that(paste(
  "check_conv_s() detects a gam() submodel fit at the iteration limit of the",
  "P-IRLS iterations"
), {
  fit_pirls <- fit_gam_tst(y_binom ~ s(z) + x, family = binomial(),
                           control = mgcv::gam.control(maxit = 1))
  skip_if_not(isFALSE(fit_pirls$converged), "P-IRLS converged")
  expect_false(check_conv_s(fit_pirls))
})

test_that(paste(
  "check_conv_s() detects a gam() submodel fit with a step failure of the",
  "Newton optimizer"
), {
  fit_newton <- fit_gam_tst(y_binom ~ s(z) + x, family = binomial(),
                            method = "REML",
                            control = mgcv::gam.control(
                              newton = list(maxHalf = 0)
                            ))
  skip_if_not(identical(fit_newton$outer.info$conv, "step failed"),
              "no step failure")
  expect_false(check_conv_s(fit_newton))
})

test_that(paste(
  "divmin() warns if the smoothing parameter optimization of gam() submodel",
  "fits stops at the iteration limit without a warning of mgcv"
), {
  d <- dat_fitters()
  mu <- cbind(d$eta, d$eta + 0.3 * d$z)
  divmin_nlm <- function(iterlim, check_conv = TRUE) {
    with_r_printing(divmin_tst(
      "s(z) + x", d, extend_family(gaussian()), mu, throw_warn = TRUE,
      check_conv = check_conv, method = "REML",
      optimizer = c("outer", "nlm"),
      control = mgcv::gam.control(nlm = list(iterlim = iterlim))
    ))
  }
  expect_no_warning(fits <- divmin_nlm(100))
  expect_length(fits, 2)
  fits_lim <- divmin_nlm(1, check_conv = FALSE)
  skip_if_not(all(vapply(fits_lim, function(fit) {
    identical(fit$outer.info$code, 4L)
  }, logical(1))), "nlm() converged in one iteration")
  expect_warning(divmin_nlm(1),
                 "^2 out of 2 submodel fits .* might not have converged")
})

# gamm4::gamm4() ----------------------------------------------------------

test_that(paste(
  "a Gaussian gamm4() submodel fit equals the REML gam() fit with a random",
  "effect smooth"
), {
  d <- dat_fitters()
  d$y_gauss <- d$eta
  fit_gam <- mgcv::gam(y_gauss ~ s(z) + x + s(g, bs = "re"), data = d,
                       method = "REML")
  fit_gamm <- divmin_tst("s(z) + x + (1 | g)", d, extend_family(gaussian()),
                         d$eta)[[1]]
  expect_s3_class(fit_gamm, "gamm4")
  # The differences are about 1e-9:
  expect_equal(unname(fitted(fit_gamm$mer)), unname(fitted(fit_gam)),
               tolerance = 1e-5)
  expect_equal(coef(fit_gamm$gam)[["x"]], coef(fit_gam)[["x"]],
               tolerance = 1e-5)
})

test_that(paste(
  "as.matrix.projection() stops for a gamm4() submodel (its parameter",
  "extraction get_subparams() stops)"
), {
  d <- dat_fitters()
  fit_gamm <- divmin_tst("s(z) + x + (1 | g)", d, extend_family(gaussian()),
                         d$eta)[[1]]
  expect_error(get_subparams(fit_gamm, nm_scheme = "rstanarm"),
               "do not support submodels with smooth and group-level terms")
})

test_that(paste(
  "the traditional and the augmented-data form of a binomial gamm4()",
  "submodel give the same fit"
), {
  d <- dat_fitters()
  mu <- plogis(d$eta)
  fit_trad <- divmin_tst("s(z) + x + (1 | g)", d, extend_family(binomial()),
                         mu)[[1]]
  # Each observation gives a row with response 1 and weight mu, and a row with
  # response 0 and weight 1 - mu:
  d_aug <- rbind(d, d)
  fit_aug <- divmin_tst("s(z) + x + (1 | g)", d_aug,
                        extend_family(binomial()),
                        rep(c(1, 0), each = nrow(d)),
                        weights = c(mu, 1 - mu))[[1]]
  expect_equal(coef(fit_aug$gam), coef(fit_trad$gam), tolerance = 1e-4)
  expect_equal(lme4::getME(fit_aug$mer, "theta"),
               lme4::getME(fit_trad$mer, "theta"), tolerance = 1e-3)
})

test_that(paste(
  "fit_gamm_callback() passes the family from glmer_family() and its",
  "control to gamm4()"
), {
  ctrl_used <- NULL
  fam_used <- NULL
  testthat::local_mocked_bindings(
    gamm4 = function(...) {
      ctrl_used <<- list(...)$control
      fam_used <<- list(...)$family
      stop("mocked gamm4()")
    },
    .package = "gamm4"
  )
  d <- dat_fitters()
  fit_gamm_mocked <- function(family, y) {
    d$y <- y
    expect_error(fit_gamm_callback(
      formula = y ~ s(z) + x + (1 | g),
      projpred_formula_no_random = y ~ s(z) + x,
      projpred_random = ~ (1 | g), data = d, family = family,
      weights = rep(1, nrow(d))
    ), "mocked gamm4")
    return(ctrl_used)
  }
  fam_binom <- extend_family(binomial())
  ctrl <- fit_gamm_mocked(extend_family(poisson()), plogis(d$eta))
  expect_identical(ctrl$optimizer, c("bobyqa", "bobyqa"))
  # The family from glmer_family() and its control:
  suffix <- if (lme4_needs_glmer_family()) "_projpred" else ""
  expect_identical(fam_used$family, paste0("poisson", suffix))
  expect_identical(isFALSE(ctrl$disp_dof_correction),
                   lme4_needs_glmer_family() &&
                     "disp_dof_correction" %in%
                       names(formals(lme4::glmerControl)))
  ctrl <- fit_gamm_mocked(fam_binom, plogis(d$eta))
  expect_identical(fam_used$family, paste0("binomial", suffix))
  ctrl <- fit_gamm_mocked(extend_family(gaussian()), d$eta)
  expect_identical(fam_used$family, "gaussian")
})

test_that(paste(
  "fit_gamm_callback() stops after the optimx refit for a binomial submodel",
  "with non-integer data"
), {
  skip_if_not_installed("optimx")
  optimizers <- list()
  testthat::local_mocked_bindings(
    gamm4 = function(...) {
      optimizers[[length(optimizers) + 1L]] <<- list(...)$control$optimizer
      stop("Downdated VtV is not positive definite")
    },
    .package = "gamm4"
  )
  d <- dat_fitters()
  d$y <- plogis(d$eta)
  expect_error(fit_gamm_callback(
    formula = y ~ s(z) + x + (1 | g),
    projpred_formula_no_random = y ~ s(z) + x,
    projpred_random = ~ (1 | g), data = d, family = extend_family(binomial()),
    weights = rep(1, nrow(d))
  ), "cannot fix this automatically")
  # The first fit and the refit with optimx:
  expect_length(optimizers, 2L)
  expect_identical(optimizers[[2]], c("optimx", "optimx"))
})

# ordinal::clmm() ---------------------------------------------------------

cats_ord <- c("1", "2", "3", "4")
fam_cumul <- structure(list(family = "cumulative_rstanarm", link = "logit",
                            cats = cats_ord),
                       class = "family")
probs_cumul <- function(eta, thres = c(-1, 0, 1)) {
  t(vapply(eta, function(eta_i) diff(c(0, plogis(thres - eta_i), 1)),
           numeric(length(thres) + 1)))
}
pars_clmm <- function(fit) {
  c(fit$coefficients, fit$ST[[1]])
}

test_that(paste(
  "a clmm() submodel fit without group-level variation equals the polr()",
  "fit"
), {
  withr::local_options(projpred.mlvl_laplace = FALSE)
  d <- dat_fitters(sd_grp = 0)
  probs <- probs_cumul(d$eta)
  fit_clmm <- divmin_augdat_tst("x + (1 | g)", d, fam_cumul, probs)[[1]]
  # divmin_augdat() removes the warning of polr() about the non-integer data
  # of the glm() fit for its start values:
  expect_no_warning(fit_polr <- with_r_printing(
    divmin_augdat_tst("x", d, fam_cumul, probs, throw_warn = TRUE)[[1]]
  ))
  expect_s3_class(fit_clmm, "clmm")
  expect_s3_class(fit_polr, "polr")
  expect_lt(fit_clmm$ST[[1]][1, 1], 1e-3)
  expect_equal(unname(fit_clmm$beta), unname(coef(fit_polr)),
               tolerance = 1e-4)
  expect_equal(unname(fit_clmm$alpha), unname(fit_polr$zeta),
               tolerance = 1e-4)
})

test_that(paste(
  "observations with weight 0 do not change a clmm() submodel fit, and",
  "integer weights equal duplicated observations"
), {
  withr::local_options(projpred.mlvl_laplace = FALSE)
  d <- dat_fitters(sd_grp = 2)
  probs <- probs_cumul(d$eta)
  fit <- divmin_augdat_tst("x + (1 | g)", d, fam_cumul, probs)[[1]]
  expect_gt(fit$ST[[1]][1, 1], 0.1)

  probs_w0 <- probs
  probs_w0[idxs_w0, ] <- 0
  fit_w0 <- divmin_augdat_tst("x + (1 | g)", d, fam_cumul, probs_w0)[[1]]
  fit_chg <- divmin_augdat_tst("x + (1 | g)", change_w0(d), fam_cumul,
                               probs_w0)[[1]]
  fit_rm <- divmin_augdat_tst("x + (1 | g)", d[-idxs_w0, ], fam_cumul,
                              probs[-idxs_w0, ])[[1]]
  expect_equal(pars_clmm(fit_chg), pars_clmm(fit_w0), tolerance = 1e-6)
  expect_equal(pars_clmm(fit_rm), pars_clmm(fit_w0), tolerance = 1e-6)

  probs_w2 <- probs
  probs_w2[idxs_w2, ] <- 2 * probs_w2[idxs_w2, ]
  idxs_dupl <- c(seq_len(nrow(d)), idxs_w2)
  fit_w2 <- divmin_augdat_tst("x + (1 | g)", d, fam_cumul, probs_w2)[[1]]
  fit_dupl <- divmin_augdat_tst("x + (1 | g)", d[idxs_dupl, ], fam_cumul,
                                probs[idxs_dupl, ])[[1]]
  expect_equal(pars_clmm(fit_dupl), pars_clmm(fit_w2), tolerance = 1e-5)
})

# mclogit::mblogit() ------------------------------------------------------

cats_cat <- c("1", "2", "3")
fam_categ <- structure(list(family = "categorical", link = "logit",
                            cats = cats_cat),
                       class = "family")
probs_categ <- function(d) {
  eta_mat <- cbind(0, d$eta, -0.5 * d$x)
  exp(eta_mat) / rowSums(exp(eta_mat))
}

test_that(paste(
  "an mblogit() submodel fit without group-level variation equals the",
  "multinom() fit"
), {
  # When the group-level variance goes to 0, the inner iterations of mblogit()
  # warn about a false convergence of nlminb(), but the outer iterations
  # converge. With this seed, mblogit() takes about 0.5 s (several seconds for
  # some other seeds).
  d <- dat_fitters(seed = 4, sd_grp = 0)
  probs <- probs_categ(d)
  fit_mblogit <- suppressWarnings(
    divmin_augdat_tst("x + (1 | g)", d, fam_categ, probs)[[1]]
  )
  fit_multinom <- divmin_augdat_tst("x", d, fam_categ, probs)[[1]]
  expect_s3_class(fit_mblogit, "mmblogit")
  expect_s3_class(fit_multinom, "multinom")
  expect_true(fit_mblogit$converged)
  expect_lt(max(abs(unlist(fit_mblogit$VarCov))), 1e-3)
  coef_multinom <- coef(fit_multinom)
  expect_equal(
    unname(fit_mblogit$coefmat[rownames(coef_multinom),
                               colnames(coef_multinom)]),
    unname(coef_multinom), tolerance = 1e-3
  )
})

# Shared setups of divmin() -----------------------------------------------

# Runs divmin_tst() with two Gaussian draws, with the warnings of divmin() and
# the convergence check, with the global option `opt_nm` set to `shared`.
# `count_fun` is the name of a function of the shared path; its calls with a
# non-`NULL` result are counted. Returns the fits, the warnings and the count.
divmin_shared_tst <- function(opt_nm, shared, count_fun, rhs, d) {
  opts_old <- options(stats::setNames(list(shared), opt_nm))
  on.exit(options(opts_old))
  mu <- cbind(d$eta, d$eta + 0.3 * d$z)
  n_calls <- 0
  fun_orig <- get(count_fun, envir = asNamespace("projpred"))
  fun_count <- function(...) {
    out <- fun_orig(...)
    if (!is.null(out)) n_calls <<- n_calls + 1
    out
  }
  warns <- character()
  fits <- withCallingHandlers(
    with_mocked_bindings(
      with_r_printing(
        divmin_tst(rhs, d, extend_family(gaussian()), mu, throw_warn = TRUE,
                   check_conv = TRUE)
      ),
      !!count_fun := fun_count
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fits = fits, warns = warns, n_calls = n_calls)
}

test_that(paste(
  "the shared lmer() and gam() setups give the warnings and fits of the",
  "draw-wise fits"
), {
  cases <- list(
    lmer = list(opt_nm = "projpred.lmer_shared", count_fun = "lmer_shared_fit",
                rhs = "x + (1 | g)"),
    gam = list(opt_nm = "projpred.gam_shared",
               count_fun = "gam_shared_callback", rhs = "s(z) + x")
  )
  # Without group-level variation, the lmer() fits are singular and lme4 gives
  # a message for each draw:
  for (sd_grp in c(0.8, 0)) {
    d <- dat_fitters(sd_grp = sd_grp)
    for (nm in names(cases)) {
      info <- paste(nm, "sd_grp:", sd_grp)
      args <- cases[[nm]]
      res_shared <- divmin_shared_tst(args$opt_nm, TRUE, args$count_fun,
                                      args$rhs, d)
      res_drawwise <- divmin_shared_tst(args$opt_nm, FALSE, args$count_fun,
                                        args$rhs, d)
      expect_gt(res_shared$n_calls, 0)
      expect_identical(res_drawwise$n_calls, 0)
      if (nm == "lmer" && sd_grp == 0) {
        expect_gt(length(res_drawwise$warns), 0)
      }
      expect_identical(res_shared$warns, res_drawwise$warns, info = info)
      expect_equal(res_shared$fits, res_drawwise$fits, info = info)
    }
  }
})

# A gradient tolerance that no fit meets:
cc_grad_tight <- function() {
  lme4::.makeCC("warning", tol = 1e-12, relTol = NULL)
}

test_that(paste(
  "gamm_check_conv() flags the gradient of a gamm4() fit with only",
  "group-level terms as lme4::checkConv() in glmer() and lmer() does"
), {
  d <- dat_fitters(n_grp = 10, n_per = 20)
  withr::local_seed(2)
  d$yb <- rbinom(nrow(d), 1, plogis(d$eta))
  d$yg <- d$eta + rnorm(nrow(d), sd = 0.5)
  for (fam_nm in c("binomial", "gaussian")) {
    if (fam_nm == "binomial") {
      ctrl <- lme4::glmerControl(optimizer = "bobyqa",
                                 check.conv.grad = cc_grad_tight())
      fit_ref <- suppressWarnings(lme4::glmer(
        yb ~ x + (1 | g), data = d, family = binomial(), control = ctrl
      ))
      fit_gamm <- gamm4::gamm4(yb ~ x, random = ~ (1 | g), data = d,
                               family = binomial(), control = ctrl)
    } else {
      ctrl <- lme4::lmerControl(check.conv.grad = cc_grad_tight())
      fit_ref <- suppressWarnings(lme4::lmer(yg ~ x + (1 | g), data = d,
                                             control = ctrl))
      fit_gamm <- gamm4::gamm4(yg ~ x, random = ~ (1 | g), data = d,
                               control = ctrl)
    }
    class(fit_gamm) <- "gamm4"
    expect_true(check_conv_s(fit_gamm), info = fam_nm)
    expect_warning(fit_gamm <- gamm_check_conv(fit_gamm, control = ctrl),
                   "failed to converge", info = fam_nm)
    cc <- fit_gamm$mer@optinfo$conv$lme4
    cc_ref <- fit_ref@optinfo$conv$lme4
    expect_identical(cc$code, -1L, info = fam_nm)
    expect_identical(cc_ref$code, -1L, info = fam_nm)
    # The component index of lme4 is always 1, so only the gradient value is
    # compared:
    grad_of <- function(msg) {
      as.numeric(sub(".*max\\|grad\\| = ([^ ]+) .*", "\\1", msg[1]))
    }
    expect_equal(grad_of(cc$messages), grad_of(cc_ref$messages),
                 tolerance = 1e-5, info = fam_nm)
    expect_false(check_conv_s(fit_gamm), info = fam_nm)
  }
})

test_that(paste(
  "gamm_check_conv() uses the unscaled gradient for a Hessian that is not",
  "positive definite, and does not check fits at a boundary"
), {
  d <- dat_fitters(n_grp = 10, n_per = 20)
  withr::local_seed(5)
  d$yb <- rbinom(nrow(d), 1, plogis(d$eta))
  ctrl <- lme4::glmerControl(optimizer = "bobyqa")
  fit0 <- gamm4::gamm4(yb ~ x, random = ~ (1 | g), data = d,
                       family = binomial(), control = ctrl)
  class(fit0) <- "gamm4"
  # Parameters: theta, then 2 fixed effects (one of them negative, which needs
  # the lower bound `-Inf`):
  set_derivs <- function(grad, val = c(0.5, -1, 1)) {
    fit <- fit0
    fit$mer@optinfo$val <- val
    fit$mer@optinfo$derivs <- list(gradient = grad, Hessian = -diag(3))
    fit
  }
  tol <- ctrl$checkConv$check.conv.grad$tol
  expect_no_warning(fit <- gamm_check_conv(set_derivs(rep(tol / 2, 3)), ctrl))
  expect_true(check_conv_s(fit))
  expect_warning(fit <- gamm_check_conv(set_derivs(c(0, 2 * tol, 0)), ctrl),
                 "max|grad| = 0.004 (tol = 0.002, component 2)", fixed = TRUE)
  expect_identical(fit$mer@optinfo$conv$lme4$code, -1L)
  expect_false(check_conv_s(fit))
  # As in lme4::checkConv(), no warning for an NA gradient:
  expect_no_warning(fit <- gamm_check_conv(set_derivs(c(NA, 0, 0)), ctrl))
  expect_identical(fit$mer@optinfo$conv$lme4$code, -5L)
  expect_false(check_conv_s(fit))
  # A variance parameter at its lower bound:
  expect_no_warning(fit <- gamm_check_conv(
    set_derivs(rep(1, 3), val = c(0, -1, 1)), ctrl
  ))
  expect_true(check_conv_s(fit))
  # The action of `check.conv.grad`:
  ctrl_act <- ctrl
  ctrl_act$checkConv$check.conv.grad$action <- "ignore"
  expect_no_warning(gamm_check_conv(set_derivs(rep(1, 3)), ctrl_act))
  ctrl_act$checkConv$check.conv.grad$action <- "stop"
  expect_error(gamm_check_conv(set_derivs(rep(1, 3)), ctrl_act),
               "failed to converge")
  # An existing result of lme4::checkConv() is kept:
  fit <- set_derivs(rep(1, 3))
  fit$mer@optinfo$conv$lme4 <- list()
  expect_no_warning(fit <- gamm_check_conv(fit, ctrl))
  expect_identical(fit$mer@optinfo$conv$lme4, list())
})

test_that(paste(
  "divmin() warns about a gamm4() submodel fit that does not meet the",
  "gradient tolerance"
), {
  d <- dat_fitters(n_grp = 10, n_per = 20)
  withr::local_seed(3)
  mu <- d$eta + rnorm(nrow(d), sd = 0.3)
  fam <- extend_family(gaussian())
  expect_warning(
    divmin_tst("s(z) + x + (1 | g)", d, fam, mu, check_conv = TRUE,
               control = lme4::lmerControl(check.conv.grad = cc_grad_tight())),
    "might not have converged"
  )
  expect_no_warning(
    divmin_tst("s(z) + x + (1 | g)", d, fam, mu, check_conv = TRUE)
  )
})

test_that(paste(
  "divmin() does not warn about Gaussian gamm4() submodel fits with a",
  "response without noise"
), {
  # For these data sets, the finite-difference Hessian of lme4::checkConv() is
  # not positive definite (with R 4.5.3, lme4 2.0.6, gamm4 0.3-0), although the
  # fits converged.
  for (seed in c(2, 6, 17)) {
    withr::local_seed(seed)
    n_grp <- 20
    n_per <- 10
    n <- n_grp * n_per
    d <- data.frame(g = factor(rep(seq_len(n_grp), each = n_per)),
                    x = rnorm(n), z = rnorm(n), v = rnorm(n))
    mu <- -0.2 + 0.5 * d$x + sin(2 * d$z) + rnorm(n_grp, sd = 0.7)[d$g]
    expect_no_warning(divmin_tst(
      "s(z) + x + (1 | g)", d, extend_family(gaussian()), mu,
      check_conv = TRUE
    ))
  }
})

test_that(paste(
  "fit_gamm_callback() does not report the smooth term of a linear function",
  "as a singular fit"
), {
  d <- dat_fitters(n_grp = 10, n_per = 20)
  withr::local_seed(4)
  d$y <- d$eta - sin(2 * d$z) + 0.5 * d$z + rnorm(nrow(d), sd = 0.3)
  expect_no_message(fit <- fit_gamm_callback(
    formula = y ~ s(z) + x + (1 | g), projpred_formula_no_random = y ~ s(z) + x,
    projpred_random = ~ (1 | g), data = d, family = extend_family(gaussian())
  ))
  # The variance of the smooth term is at the boundary (the fit is singular).
  # Depending on the machine, the optimizer stops at zero or above zero, but
  # below the singular tolerance 1e-4 of gamm_check_conv():
  expect_lt(lme4::getME(fit$mer, "theta")[["Xr.s(z)"]], 1e-4)
  expect_true(check_conv_s(fit))
})
