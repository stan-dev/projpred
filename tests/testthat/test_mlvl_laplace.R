# projpred's Laplace fit of binomial and Poisson multilevel submodels with
# random intercepts (mlvl_laplace()). These tests need no reference model fits,
# so they also run on CRAN.

# Data with two grouping variables (`g` with 8 levels, `h` with 4 levels, with
# disjoint level names), a linear predictor `eta` and integer responses.
dat_mlvl <- function(n = 160, seed = 3) {
  withr::local_seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(paste0("g", sample(8, n, replace = TRUE))),
                  h = factor(paste0("h", sample(4, n, replace = TRUE))))
  d$eta <- -0.2 + 0.7 * d$x - 0.4 * d$z + rnorm(8, sd = 0.8)[d$g] +
    rnorm(4, sd = 0.4)[d$h]
  d$y_bern <- rbinom(n, 1, plogis(d$eta))
  d$y_pois <- rpois(n, exp(d$eta))
  d
}

fit_mlvl_tst <- function(formula, data, family, weights = rep(1, nrow(data)),
                         ...) {
  fit_glmer_callback(formula, projpred_formula_no_random = NA,
                     projpred_random = NA, family = extend_family(family),
                     data = data, weights = weights, ...)
}

glmer_tst <- function(formula, data, family, weights = rep(1, nrow(data))) {
  data$.w_tst <- weights
  suppressWarnings(suppressMessages(lme4::glmer(
    formula, data = data, family = family, weights = .w_tst,
    control = lme4::glmerControl(optimizer = "bobyqa")
  )))
}

test_that("mlvl_laplace_c() gives the gradient of its objective", {
  d <- dat_mlvl()
  d$y_frac <- plogis(d$eta + 0.3)
  d$y_mean <- exp(d$eta)
  cases <- list(bern = list(y = d$y_bern, fam = 0L),
                frac = list(y = d$y_frac, fam = 0L),
                pois = list(y = d$y_pois, fam = 1L),
                mean = list(y = d$y_mean, fam = 1L))
  x <- cbind(1, d$x, d$z)
  grp <- cbind(as.integer(d$g), as.integer(d$h)) - 1L
  nlev <- c(8L, 4L)
  w <- rep(c(1, 2), length.out = nrow(d))
  par <- c(0.6, 0.3, -0.1, 0.5, -0.3)
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    obj <- function(par) {
      mlvl_laplace_c(par[1:2], par[3:5], x, grp, nlev, cs$y, w,
                     rep(0, nrow(d)), rep(0, 12), cs$fam, FALSE)$D
    }
    ev <- mlvl_laplace_c(par[1:2], par[3:5], x, grp, nlev, cs$y, w,
                         rep(0, nrow(d)), rep(0, 12), cs$fam, TRUE)
    grad_fd <- vapply(seq_along(par), function(j) {
      h <- 1e-5
      e <- replace(numeric(length(par)), j, h)
      (obj(par + e) - obj(par - e)) / (2 * h)
    }, numeric(1))
    expect_equal(as.vector(ev$grad), grad_fd, tolerance = 1e-6, info = nm)
    expect_equal(ev$D, obj(par), info = nm)
  }
})

test_that(paste(
  "mlvl_laplace_c() gives the deviance for linear predictors where mu rounds",
  "to 0 or 1"
), {
  # With theta = 0, the conditional modes are 0, so D is the deviance at
  # eta = offset:
  dev_c <- function(y, eta, fam) {
    n <- length(y)
    mlvl_laplace_c(0, 0, matrix(1, n, 1), matrix(0L, n, 1), 1L, y, rep(1, n),
                   eta, 0, fam, FALSE)$D
  }
  xlogx <- function(x) ifelse(x > 0, x * log(x), 0)
  # Relative error, also for targets near 0 (where expect_equal() compares
  # absolute differences):
  expect_rel <- function(object, expected, info) {
    if (expected == 0) {
      expect_identical(object, 0, info = info)
    } else {
      expect_lt(abs(object / expected - 1), 1e-12, label = info)
    }
  }
  # Binomial: in double precision, plogis() is 0 for eta -800 and 1 for eta 37
  # and 100:
  y <- rep(c(0, 0.3, 1 - 1e-6, 1), 4)
  eta <- rep(c(-800, 2, 37, 100), each = 4)
  dev_ref <- 2 * (xlogx(y) + xlogx(1 - y) -
                    y * plogis(eta, log.p = TRUE) -
                    (1 - y) * plogis(-eta, log.p = TRUE))
  for (i in seq_along(y)) {
    expect_rel(dev_c(y[i], eta[i], 0L), dev_ref[i],
               info = paste("binomial: y", y[i], "eta", eta[i]))
  }
  y2 <- y[eta == 2]
  expect_equal(vapply(y2, dev_c, numeric(1), eta = 2, fam = 0L),
               binomial()$dev.resids(y2, rep(plogis(2), 4), rep(1, 4)),
               tolerance = 1e-12)
  # Poisson: in double precision, exp() is 0 for eta -800 (reference by hand,
  # with log(mu) = eta and mu = 0):
  y <- c(0, 0.5, 3)
  expect_rel(dev_c(y, rep(-800, 3), 1L), 2 * sum(xlogx(y) + 800 * y - y),
             info = "poisson: eta -800")
  expect_equal(dev_c(y, rep(1, 3), 1L),
               sum(poisson()$dev.resids(y, rep(exp(1), 3), rep(1, 3))),
               tolerance = 1e-12)
  # Near the optimum, the terms of each observation nearly cancel.  The
  # deviance keeps its relative accuracy only if each observation adds its
  # terms before the sum:
  withr::local_seed(1)
  n <- 1e4
  y <- rpois(n, 1000)
  eta <- log(y) + rnorm(n, sd = 0.02)
  expect_equal(dev_c(y, eta, 1L),
               sum(poisson()$dev.resids(y, exp(eta), rep(1, n))),
               tolerance = 1e-12)
  eta <- rnorm(n, sd = 2)
  y <- plogis(eta + rnorm(n, sd = 0.01))
  expect_equal(dev_c(y, eta, 0L),
               sum(binomial()$dev.resids(y, plogis(eta), rep(1, n))),
               tolerance = 1e-11)
})

test_that(paste(
  "mlvl_laplace() gives the glmer() fit for integer binomial and Poisson",
  "responses"
), {
  d <- dat_mlvl()
  cases <- list(
    bern = list(fml = y_bern ~ x + z + (1 | g) + (1 | h), fam = binomial()),
    pois = list(fml = y_pois ~ x + z + (1 | g) + (1 | h), fam = poisson()),
    pois_one = list(fml = y_pois ~ x + (1 | h), fam = poisson())
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    fit <- fit_mlvl_tst(cs$fml, d, cs$fam)
    expect_s3_class(fit, "projpred_mlvl")
    expect_true(check_conv_s(fit), info = nm)
    fit_g <- glmer_tst(cs$fml, d, cs$fam)
    expect_equal(fit$coefficients, lme4::fixef(fit_g), tolerance = 1e-3,
                 info = nm)
    # The parameters and the predictions in the names of the lme4 path:
    for (nm_scheme in c("rstanarm", "brms")) {
      expect_equal(get_subparams(fit, nm_scheme = nm_scheme),
                   get_subparams(fit_g, nm_scheme = nm_scheme),
                   tolerance = 1e-3, info = paste(nm, nm_scheme))
    }
    expect_equal(subprd(list(fit, fit), d)[, 2], subprd(list(fit_g), d)[, 1],
                 tolerance = 1e-3, info = nm)
    expect_equal(predict(fit, d), predict(fit_g, d), tolerance = 1e-3,
                 info = nm)
    # For a 0-1 response, the objective is -2 times the log-likelihood of
    # lme4:
    if (nm == "bern") {
      expect_equal(fit$deviance, -2 * as.numeric(logLik(fit_g)),
                   tolerance = 1e-5)
    }
  }
})

test_that(paste(
  "mlvl_laplace() gives the glmer() fit for non-integer responses with",
  "lme4 >= 2.1-1"
), {
  skip_if(utils::packageVersion("lme4") < "2.1.1",
          "lme4 < 2.1-1 rounds non-integer binomial and Poisson data")
  d <- dat_mlvl()
  d$y_frac <- plogis(d$eta + 0.3 * d$z)
  d$y_mean <- exp(d$eta + 0.3 * d$z)
  for (cs in list(list(fml = y_frac ~ x + (1 | g) + (1 | h), fam = binomial()),
                  list(fml = y_mean ~ x + (1 | g) + (1 | h), fam = poisson()))) {
    fit <- fit_mlvl_tst(cs$fml, d, cs$fam)
    fit_g <- glmer_tst(cs$fml, d, cs$fam)
    expect_equal(fit$coefficients, lme4::fixef(fit_g), tolerance = 1e-3)
    expect_equal(unname(fit$sd),
                 unname(vapply(lme4::VarCorr(fit_g)[names(fit$sd)], attr,
                               numeric(1), "stddev")),
                 tolerance = 1e-2)
  }
})

test_that(paste(
  "the traditional and the augmented-data form of a binomial submodel give",
  "the same mlvl_laplace() fit"
), {
  d <- dat_mlvl()
  d$y <- plogis(d$eta + 0.3 * d$z)
  d_aug <- rbind(d, d)
  d_aug$y <- rep(c(1, 0), each = nrow(d))
  fit_trad <- fit_mlvl_tst(y ~ x + (1 | g), d, binomial())
  fit_aug <- fit_mlvl_tst(y ~ x + (1 | g), d_aug, binomial(),
                          weights = c(d$y, 1 - d$y))
  # The two objectives differ by a constant, so the optimizer stops at slightly
  # different points (relative difference 2.3e-5 on one platform):
  expect_equal(fit_aug$coefficients, fit_trad$coefficients, tolerance = 1e-4)
  expect_equal(fit_aug$sd, fit_trad$sd, tolerance = 1e-4)
  # A factor response, as in the augmented-data projection (first level:
  # failure):
  d_aug$y <- factor(d_aug$y, levels = c(0, 1))
  fit_fac <- fit_mlvl_tst(y ~ x + (1 | g), d_aug, binomial(),
                          weights = c(d$y, 1 - d$y))
  expect_s3_class(fit_fac, "projpred_mlvl")
  expect_equal(fit_fac$coefficients, fit_aug$coefficients, tolerance = 1e-10)
})

test_that(paste(
  "a Poisson submodel with non-integer responses estimates the group-level",
  "standard deviation"
), {
  withr::local_seed(4)
  g <- factor(rep(seq_len(8), each = 6))
  x <- rnorm(48)
  b <- rnorm(8, sd = 0.8)
  d <- data.frame(y = exp(0.5 + 0.4 * x + b[g] + rnorm(48, sd = 0.1)),
                  x = x, g = g)
  fit <- fit_mlvl_tst(y ~ x + (1 | g), d, poisson())
  expect_s3_class(fit, "projpred_mlvl")
  expect_gt(abs(fit$sd[["g"]] - 1), 1e-4)
  expect_gt(fit$sd[["g"]], 0.2)
  expect_true(is.finite(fit$deviance))
})

test_that(paste(
  "integer weights in mlvl_laplace() equal duplicated observations, and",
  "divmin() removes observations with weight 0"
), {
  d <- dat_mlvl()
  d$mu_b <- plogis(d$eta)
  d$mu_p <- exp(d$eta)
  idxs_w2 <- 1:10
  w <- rep(1, nrow(d))
  w[idxs_w2] <- 2
  d_dupl <- d[c(seq_len(nrow(d)), idxs_w2), ]
  for (cs in list(list(y = "mu_b", fam = binomial()),
                  list(y = "mu_p", fam = poisson()))) {
    fml <- as.formula(paste(cs$y, "~ x + (1 | g)"))
    fit_w <- fit_mlvl_tst(fml, d, cs$fam, weights = w)
    fit_dupl <- fit_mlvl_tst(fml, d_dupl, cs$fam)
    expect_equal(fit_w$coefficients, fit_dupl$coefficients, tolerance = 1e-6)
    expect_equal(fit_w$sd, fit_dupl$sd, tolerance = 1e-5)
  }
  # Weight 0 (the left-out observation of a fold in LOO CV with validated
  # search), with changed predictor values:
  w0 <- rep(1, nrow(d))
  w0[1:3] <- 0
  d_chg <- d
  d_chg$x[1:3] <- d_chg$x[1:3] + 5
  fit_divmin <- function(data) {
    data$.y.1 <- data$mu_b
    divmin(formula = .y.1 ~ x + (1 | g), data = data,
           family = extend_family(binomial()), weights = w0,
           projpred_var = matrix(NA_real_, nrow(data), 1),
           projpred_ws_aug = matrix(data$mu_b), verbose_divmin = FALSE,
           throw_warn_sdivmin = FALSE, do_check_conv = FALSE)[[1]]
  }
  fit_w0 <- fit_divmin(d)
  fit_chg <- fit_divmin(d_chg)
  expect_s3_class(fit_w0, "projpred_mlvl")
  expect_equal(fit_chg$coefficients, fit_w0$coefficients, tolerance = 1e-10)
  fit_rm <- fit_mlvl_tst(mu_b ~ x + (1 | g), d[-(1:3), ], binomial())
  expect_equal(fit_rm$coefficients, fit_w0$coefficients, tolerance = 1e-6)
})

test_that(paste(
  "fit_glmer_callback() uses glmer() for the cases that mlvl_laplace() does",
  "not handle"
), {
  d <- dat_mlvl()
  expect_glmer <- function(fit) {
    expect_s4_class(fit, "glmerMod")
  }
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  # Random slopes:
  expect_glmer(quiet(fit_mlvl_tst(y_bern ~ x + (x | g), d, binomial())))
  # Another link:
  expect_glmer(quiet(fit_mlvl_tst(y_bern ~ x + (1 | g), d,
                                  binomial(link = "probit"))))
  # A `control` or an `nAGQ` from the caller:
  expect_glmer(quiet(fit_mlvl_tst(y_bern ~ x + (1 | g), d, binomial(),
                                  control = lme4::glmerControl())))
  expect_glmer(quiet(fit_mlvl_tst(y_bern ~ x + (1 | g), d, binomial(),
                                  nAGQ = 5)))
  # An offset:
  d$offs <- 0.1 * d$z
  expect_glmer(quiet(fit_mlvl_tst(y_pois ~ x + offset(offs) + (1 | g), d,
                                  poisson())))
  # A rank-deficient model matrix:
  d$x2 <- 2 * d$x
  expect_glmer(quiet(fit_mlvl_tst(y_bern ~ x + x2 + (1 | g), d, binomial())))
  # The global option:
  withr::local_options(projpred.mlvl_laplace = FALSE)
  expect_glmer(quiet(fit_mlvl_tst(y_bern ~ x + (1 | g), d, binomial())))
})

test_that(paste(
  "fit_glmer_callback() refits with glmer() if mlvl_laplace() does not",
  "converge, and returns the glmer() fit as a `projpred_mlvl` object"
), {
  d <- dat_mlvl()
  fml <- y_bern ~ x + (1 | g)
  local_mocked_bindings(mlvl_laplace_optim = function(des, fam, ...) {
    list(theta = 1, beta = rep(0, ncol(des$x)), u = rep(0, sum(des$nlev)),
         deviance = 1, convergence = 1L, message = "mocked",
         iterations = 0L, n_restart = 0L)
  })
  fit <- suppressWarnings(suppressMessages(fit_mlvl_tst(fml, d, binomial())))
  expect_s3_class(fit, "projpred_mlvl")
  expect_identical(fit$optinfo$message,
                   "lme4::glmer() fit (Laplace fit: mocked)")
  fit_g <- glmer_tst(fml, d, binomial())
  expect_equal(fit$coefficients, lme4::fixef(fit_g), tolerance = 1e-8)
  expect_equal(get_subparams(fit, nm_scheme = "rstanarm"),
               get_subparams(fit_g, nm_scheme = "rstanarm"), tolerance = 1e-8)
  # The deviance on the scale of mlvl_laplace():
  des <- mlvl_laplace_design(fml, d, rep(1, nrow(d)), "g", 0L)
  expect_equal(fit$deviance,
               mlvl_laplace_dev(des, 0L, fit$sd, unname(fit$coefficients)))
  # Non-integer binomial data: the glmer() fallback uses the family of the
  # glmer() code path (see glmer_family()):
  d$y_frac <- plogis(d$eta + 0.3 * d$z)
  fml_frac <- y_frac ~ x + (1 | g)
  fit_frac <- suppressWarnings(suppressMessages(
    fit_mlvl_tst(fml_frac, d, binomial())
  ))
  expect_s3_class(fit_frac, "projpred_mlvl")
  fam_g <- glmer_family(binomial())
  fit_g_frac <- suppressWarnings(suppressMessages(lme4::glmer(
    fml_frac, data = d, family = fam_g,
    control = glmer_control(lme4::glmerControl(optimizer = "bobyqa"), fam_g)
  )))
  expect_equal(fit_frac$coefficients, lme4::fixef(fit_g_frac),
               tolerance = 1e-6)
  expect_identical(fit_frac$convergence, 0L)
  des_frac <- mlvl_laplace_design(fml_frac, d, rep(1, nrow(d)), "g", 0L)
  expect_equal(fit_frac$deviance,
               mlvl_laplace_dev(des_frac, 0L, fit_frac$sd,
                                unname(fit_frac$coefficients)))
})

test_that(paste(
  "as_glmerMod() gives a `glmerMod` object with the parameters and the",
  "predictions of the `projpred_mlvl` fit"
), {
  d <- dat_mlvl()
  d$y_frac <- plogis(d$eta + 0.3 * d$z)
  d$y_mean <- exp(d$eta + 0.3 * d$z)
  w <- rep(c(1, 2), length.out = nrow(d))
  nd <- d[1:6, ]
  nd$g <- factor(c("g_new1", "g1", "g2", "g3", "g_new2", "g4"))
  cases <- list(
    bern = list(fml = y_bern ~ x + z + (1 | g) + (1 | h), fam = binomial()),
    frac = list(fml = y_frac ~ x + (1 | h) + (1 | g), fam = binomial()),
    pois = list(fml = y_pois ~ x + (1 | g), fam = poisson()),
    mean = list(fml = y_mean ~ x + z + (1 | g) + (1 | h), fam = poisson())
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    fit <- fit_mlvl_tst(cs$fml, d, cs$fam, weights = w)
    expect_s3_class(fit, "projpred_mlvl")
    # No warning, also not for the non-integer responses:
    expect_no_warning(fit_g <- as_glmerMod(fit, data = d, weights = w))
    expect_s4_class(fit_g, "glmerMod")
    expect_equal(lme4::fixef(fit_g), fit$coefficients, tolerance = 1e-12,
                 info = nm)
    for (nm_scheme in c("rstanarm", "brms")) {
      expect_equal(get_subparams(fit_g, nm_scheme = nm_scheme),
                   get_subparams(fit, nm_scheme = nm_scheme),
                   tolerance = 1e-6, info = paste(nm, nm_scheme))
    }
    expect_equal(predict(fit_g, d), predict(fit, d), tolerance = 1e-6,
                 info = nm)
    expect_equal(predict(fit_g, nd, allow.new.levels = TRUE),
                 predict(fit, nd), tolerance = 1e-6, info = nm)
  }
  # The fit that mlvl_from_merMod() builds after the glmer() fallback:
  local_mocked_bindings(mlvl_laplace_optim = function(des, fam, ...) {
    list(theta = 1, beta = rep(0, ncol(des$x)), u = rep(0, sum(des$nlev)),
         deviance = 1, convergence = 1L, message = "mocked",
         iterations = 0L, n_restart = 0L)
  })
  fit_fb <- suppressWarnings(suppressMessages(
    fit_mlvl_tst(y_bern ~ x + (1 | g), d, binomial())
  ))
  expect_s3_class(fit_fb, "projpred_mlvl")
  fit_fb_g <- as_glmerMod(fit_fb, data = d)
  expect_equal(predict(fit_fb_g, d), predict(fit_fb, d), tolerance = 1e-6)
})

test_that(paste(
  "as_glmerMod() stores a glmer() call, removes observations with weight 0",
  "and handles the factor response of the augmented-data projection"
), {
  d <- dat_mlvl()
  d$y_frac <- plogis(d$eta + 0.3 * d$z)
  w <- rep(c(1, 2), length.out = nrow(d))
  fit <- fit_mlvl_tst(y_frac ~ x + (1 | g), d, binomial(), weights = w)
  fit_g <- as_glmerMod(fit, data = d, weights = w)
  cl <- getCall(fit_g)
  expect_identical(cl[[1]], quote(lme4::glmer))
  expect_identical(cl$formula, formula(fit))
  expect_identical(cl$data, quote(d))
  expect_identical(cl$weights, quote(w))
  expect_identical(cl$family, quote(binomial(link = "logit")))
  expect_null(getCall(as_glmerMod(
    fit_mlvl_tst(y_bern ~ x + (1 | g), d, binomial()), data = d
  ))$weights)
  # Weight 0: the fit of the other observations, with any values in the rows
  # with weight 0:
  w0 <- w
  w0[1:3] <- 0
  fit_rm <- fit_mlvl_tst(y_frac ~ x + (1 | g), d[-(1:3), ], binomial(),
                         weights = w[-(1:3)])
  d_chg <- d
  d_chg$y_frac[1:3] <- 0.5
  expect_no_warning(fit_rm_g <- as_glmerMod(fit_rm, data = d_chg,
                                            weights = w0))
  expect_identical(nrow(model.frame(fit_rm_g)), nrow(d) - 3L)
  expect_equal(predict(fit_rm_g, d), predict(fit_rm, d), tolerance = 1e-6)
  # The augmented-data form with a factor response (first level: failure):
  d_aug <- rbind(d, d)
  d_aug$y_frac <- factor(rep(c(1, 0), each = nrow(d)), levels = c(0, 1))
  w_aug <- c(d$y_frac, 1 - d$y_frac)
  fit_aug <- fit_mlvl_tst(y_frac ~ x + (1 | g), d_aug, binomial(),
                          weights = w_aug)
  expect_no_warning(fit_aug_g <- as_glmerMod(fit_aug, data = d_aug,
                                             weights = w_aug))
  expect_equal(predict(fit_aug_g, d), predict(fit_aug, d), tolerance = 1e-6)
})

test_that("as_glmerMod() checks its arguments and the data", {
  d <- dat_mlvl()
  d$y_frac <- plogis(d$eta + 0.3 * d$z)
  fit <- fit_mlvl_tst(y_frac ~ x + (1 | g), d, binomial())
  expect_error(as_glmerMod(list(), data = d), "class `projpred_mlvl`")
  expect_error(as_glmerMod(fit, data = as.matrix(d[, c("x", "z")])),
               "`data` must be a `data.frame`")
  expect_error(as_glmerMod(fit, data = d, weights = 1), "`weights` must be")
  expect_error(as_glmerMod(fit, data = cbind(d, .projpred_wobs = 1)),
               "must not have a column named `.projpred_wobs`")
  expect_error(as_glmerMod(fit, data = d[, setdiff(names(d), "z")]), NA)
  expect_error(as_glmerMod(fit, data = d[, setdiff(names(d), "x")]))
  msg_differs <- "converted fit differs from `object`"
  # Another response than that of the fit:
  d_other <- d
  d_other$y_frac <- rev(d$y_frac)
  expect_warning(as_glmerMod(fit, data = d_other), msg_differs)
  # A small change of the response, which changes the group-level intercepts
  # by less than the tolerance of their comparison:
  d_small <- d
  d_small$y_frac <- plogis(qlogis(d$y_frac) + 1e-3 * d$z)
  expect_warning(as_glmerMod(fit, data = d_small), msg_differs)
  # Other weights than those of the fit:
  expect_warning(
    as_glmerMod(fit, data = d, weights = rep(c(1, 20), each = nrow(d) / 2)),
    msg_differs
  )
  # Other group levels than those of the fit:
  d_lvls <- d[d$g != "g1", ]
  expect_warning(as_glmerMod(fit, data = d_lvls), msg_differs)
  # A fit with group-level standard deviation 0 (all group-level intercepts
  # are 0, also for another response):
  fit0 <- fit
  fit0$sd[] <- 0
  fit0$ranef$g[] <- 0
  des <- mlvl_laplace_design(formula(fit), d, rep(1, nrow(d)), "g", 0L)
  fit0$deviance <- mlvl_laplace_dev(des, 0L, 0, unname(fit0$coefficients))
  expect_no_warning(as_glmerMod(fit0, data = d))
  expect_warning(as_glmerMod(fit0, data = d_other), msg_differs)
})

test_that("divmin() passes global option `projpred.mlvl_laplace` on", {
  d <- dat_mlvl()
  d$.y.1 <- plogis(d$eta)
  fit_divmin <- function() {
    divmin(formula = .y.1 ~ x + (1 | g), data = d,
           family = extend_family(binomial()), weights = rep(1, nrow(d)),
           projpred_var = matrix(NA_real_, nrow(d), 1),
           projpred_ws_aug = matrix(d$.y.1), verbose_divmin = FALSE,
           throw_warn_sdivmin = FALSE, do_check_conv = FALSE)[[1]]
  }
  expect_s3_class(fit_divmin(), "projpred_mlvl")
  withr::local_options(projpred.mlvl_laplace = FALSE)
  expect_s4_class(suppressWarnings(suppressMessages(fit_divmin())), "glmerMod")
  # divmin() reads the option itself and passes its value to the fitter, which
  # can run in a parallel worker without the option:
  seen <- NULL
  local_mocked_bindings(fit_glmer_callback = function(...) {
    seen <<- list(...)$mlvl_laplace
    stop("recorded")
  })
  try(fit_divmin(), silent = TRUE)
  expect_identical(seen, FALSE)
  withr::local_options(projpred.mlvl_laplace = TRUE)
  try(fit_divmin(), silent = TRUE)
  expect_identical(seen, TRUE)
})

test_that(paste(
  "subprd() draws the group-level intercepts of new group levels once per",
  "level and keeps those of the existing levels"
), {
  d <- dat_mlvl()
  fit <- fit_mlvl_tst(y_bern ~ x + (1 | g) + (1 | h), d, binomial())
  nd <- d[1:6, ]
  nd$g <- factor(c("g_new1", "g_new1", "g1", "g_new2", "g2", "g_new2"))
  prd_fixed <- predict(fit, nd)
  withr::local_seed(1)
  prd <- subprd(list(fit, fit), nd)
  expect_identical(dim(prd), c(6L, 2L))
  expect_identical(rownames(prd), rownames(nd))
  b_new <- prd[, 1] - prd_fixed
  # Existing levels: no draw; new levels: one draw per level:
  expect_equal(unname(b_new[c(3, 5)]), c(0, 0))
  expect_equal(b_new[[1]], b_new[[2]])
  expect_equal(b_new[[4]], b_new[[6]])
  expect_false(isTRUE(all.equal(b_new[[1]], b_new[[4]])))
  # The draws differ between the draws of the submodel:
  expect_false(isTRUE(all.equal(prd[1, 1], prd[1, 2])))
  # With `projpred.mlvl_pred_new`, all levels are new:
  withr::local_options(projpred.mlvl_pred_new = TRUE)
  prd_new <- subprd(list(fit), nd)
  expect_false(isTRUE(all.equal(prd_new[3, 1], prd_fixed[[3]])))
  # Overlapping level sets of the grouping variables:
  nd$h <- factor(c("g1", "h1", "h1", "h1", "h1", "h1"))
  expect_error(subprd(list(fit), nd), "disjoint level sets")
})

test_that(paste(
  "cv_varsel() with LOO validated search and weight perturbation uses",
  "mlvl_laplace() with non-integer observation weights"
), {
  withr::local_seed(3)
  n <- 60
  S <- 20
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(paste0("g", sample(6, n, replace = TRUE))))
  b_g <- rnorm(6, sd = 0.8)
  eta <- -0.2 + 0.8 * d$x + b_g[d$g]
  # Reference model draws (no MCMC):
  draws <- list(beta = rbind(rnorm(S, -0.2, 0.1), rnorm(S, 0.8, 0.1),
                             rnorm(S, 0, 0.1)),
                bg = matrix(rnorm(6 * S, rep(b_g, S), 0.2), 6, S))
  ref_predfun <- function(fit, newdata = NULL) {
    if (is.null(newdata)) newdata <- d
    cbind(1, newdata$x, newdata$z) %*% fit$beta +
      fit$bg[as.integer(factor(newdata$g, levels = levels(d$g))), ,
             drop = FALSE]
  }
  mlvl_laplace_orig <- mlvl_laplace
  rec <- list()
  local_mocked_bindings(mlvl_laplace = function(formula, data, family,
                                                weights = NULL) {
    fit <- mlvl_laplace_orig(formula, data, family, weights)
    rec[[length(rec) + 1L]] <<- list(w = weights, fit = fit)
    fit
  })
  y_list <- list(binomial = rbinom(n, 1, plogis(eta)),
                 poisson = rpois(n, exp(eta)))
  for (fam_nm in names(y_list)) {
    d$y <- y_list[[fam_nm]]
    ref <- suppressWarnings(init_refmodel(
      object = draws, data = d, formula = y ~ x + z + (1 | g),
      family = get(fam_nm)(), ref_predfun = ref_predfun
    ))
    rec <- list()
    cvvs <- suppressWarnings(suppressMessages(cv_varsel(
      ref, cv_method = "LOO", validate_search = TRUE, method = "forward",
      nloo = 4, ndraws = 4, ndraws_pred = 8, nterms_max = 2, seed = 5,
      verbose = FALSE
    )))
    expect_identical(cvvs$loo_perturb_alpha, 10, info = fam_nm)
    expect_gt(length(rec), 0)
    # The perturbed searches pass non-integer weights:
    w_nonint <- vapply(rec, function(r) {
      any(abs(r$w - round(r$w)) > 1e-8)
    }, logical(1))
    expect_true(any(w_nonint), info = fam_nm)
    expect_true(all(vapply(rec[w_nonint], function(r) {
      inherits(r$fit, "projpred_mlvl")
    }, logical(1))), info = fam_nm)
    expect_true(all(is.finite(summary(cvvs)$perf_sub$elpd)), info = fam_nm)
  }
})
