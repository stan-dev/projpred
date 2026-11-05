context("submodel fitters from other packages")

# These tests run on CRAN on purpose, so that the reverse dependency checks of
# lme4, MASS, nlme, mgcv, gamm4, ordinal and mclogit see changes in the
# submodel fits that projpred needs. They call divmin() and divmin_augdat()
# directly with small simulated data (no reference model). They test
# properties of the projection (known solutions, observation weights), not
# the equality with a direct call of the same fitter.

# lme4 < 2.1-1 fits binomial glmer() submodels to rounded data and cannot
# estimate the group-level standard deviation of Poisson glmer() submodels.
skip_if_lme4_rounds <- function() {
  skip_if(utils::packageVersion("lme4") < "2.1.1",
          "lme4 < 2.1-1 rounds non-integer binomial and Poisson data")
}

# Data with `n_grp` groups of `n_per` observations. The linear predictor `eta`
# has the group-level effects with standard deviation `sd_grp`.
# An existing state of the random number generator is restored at the end.
dat_fitters <- function(n_grp = 6, n_per = 8, sd_grp = 0.8, seed = 1) {
  if (exists(".Random.seed", envir = globalenv())) {
    seed_old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", seed_old, envir = globalenv()))
  }
  set.seed(seed)
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

fit_glmer_test <- function(formula, data, family, weights, ...) {
  suppressWarnings(suppressMessages(fit_glmer_callback(
    formula, projpred_formula_no_random = NA, projpred_random = NA,
    family = family, data = data, weights = weights, ...
  )))
}

test_that(paste(
  "a binomial glmer() submodel with fractional responses gives the glm()",
  "solution for a singular fit"
), {
  skip_if_lme4_rounds()
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
  skip_if_lme4_rounds()
  # A singular fit and a fit with theta > 0 (10 groups of 30 observations):
  set.seed(11)
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
  skip_if_lme4_rounds()
  set.seed(4)
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
  skip_if_lme4_rounds()
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
  opts_old <- options(projpred.PQL = TRUE)
  fit_pql <- tryCatch(
    divmin_tst("x + (1 | g)", d, extend_family(gaussian()), d$eta)[[1]],
    finally = options(opts_old)
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
  opts_old <- options(projpred.PQL = TRUE)
  tryCatch({
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
  }, finally = options(opts_old))
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
  "the traditional and the augmented-data form of a binomial gamm4()",
  "submodel give the same fit"
), {
  skip_if_lme4_rounds()
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
