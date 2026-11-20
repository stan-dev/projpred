# Data with the response columns of several draws, as the submodel fitters get
# them:
n_scd <- 150L
idx_scd <- seq_len(n_scd)
dat_scd <- data.frame(x1 = sin(idx_scd), x2 = cos(1.3 * idx_scd),
                      g = factor(paste0("g", idx_scd %% 6L)))
u_scd <- c(0.8, -0.5, 0.3, -0.9, 0.6, -0.2)
eta_scd <- dat_scd$x1 + 0.5 * dat_scd$x2 + u_scd[as.integer(dat_scd$g)]
n_draws_scd <- 30L
for (s in seq_len(n_draws_scd)) {
  dat_scd[[paste0(".y.", s)]] <- eta_scd +
    0.3 * sin(2.3 * (s + 0.37) * idx_scd + 0.5)
}
dat_scd$yb <- as.integer(eta_scd + 1.2 * sin(2.9 * idx_scd + 0.2) > 0)
dat_scd$yo <- cut(eta_scd + 1.2 * sin(2.9 * idx_scd + 0.2),
                  c(-Inf, -0.5, 0.5, Inf),
                  labels = c("lo", "mid", "hi"), ordered_result = TRUE)
dat_scd$yc <- factor(dat_scd$yo, ordered = FALSE)
# A column that no submodel uses, so only a copy of the whole data holds it:
dat_scd$col_unused <- 1
newdat_scd <- dat_scd[c(3, 50, 77, 120), c("x1", "x2", "g")]

# Whether `x` holds a data frame with column `col` anywhere in its list
# elements, S4 slots and calls (not in environments):
has_df_with <- function(x, col, depth = 0L) {
  if (depth > 8L) {
    return(FALSE)
  }
  if (is.data.frame(x)) {
    return(col %in% names(x))
  }
  parts <- if (isS4(x)) {
    lapply(methods::slotNames(x), function(nm) methods::slot(x, nm))
  } else if (is.list(x) || is.call(x)) {
    as.list(x)
  }
  any(vapply(parts, has_df_with, logical(1), col = col, depth = depth + 1L))
}

test_that(paste(
  "strip_call_data() removes the data frame from the calls of the submodel",
  "fits and keeps their predictions for new data"
), {
  # Each case: the fit (through do.call(), as in the fitters of projpred) and
  # the predictions for new data (`NULL`: no predict() method).
  cases <- list(
    lmer = list(fit = function() {
      do.call(lme4::lmer, list(formula = .y.1 ~ x1 + x2 + (1 | g),
                               data = dat_scd))
    }, pred = function(f) predict(f, newdata = newdat_scd)),
    glmer = list(fit = function() {
      do.call(lme4::glmer, list(formula = yb ~ x1 + (1 | g), data = dat_scd,
                                family = binomial()))
    }, pred = function(f) predict(f, newdata = newdat_scd)),
    gam = list(fit = function() {
      do.call(mgcv::gam, list(formula = .y.1 ~ s(x1) + x2, data = dat_scd))
    }, pred = function(f) predict(f, newdata = newdat_scd)),
    # gamm4() fits keep no data frame in their calls:
    gamm4 = list(fit = function() {
      do.call(gamm4::gamm4, list(formula = .y.1 ~ s(x1) + x2,
                                 random = ~ (1 | g), data = dat_scd))
    }, pred = function(f) predict(f$gam, newdata = newdat_scd),
    no_data = TRUE),
    glmmPQL = list(fit = function() {
      do.call(MASS::glmmPQL, list(fixed = yb ~ x1, random = ~ 1 | g,
                                  family = binomial(), data = dat_scd,
                                  verbose = FALSE))
    }, pred = function(f) predict(f, newdata = newdat_scd, level = 0:1)),
    polr = list(fit = function() {
      do.call(MASS::polr, list(formula = yo ~ x1 + x2, data = dat_scd))
    }, pred = function(f) predict(f, newdata = newdat_scd, type = "probs")),
    clmm = list(fit = function() {
      do.call(ordinal::clmm, list(formula = yo ~ x1 + (1 | g),
                                  data = dat_scd))
    }, pred = NULL),
    multinom = list(fit = function() {
      do.call(nnet::multinom, list(formula = yc ~ x1 + x2, data = dat_scd,
                                   trace = FALSE))
    }, pred = function(f) predict(f, newdata = newdat_scd, type = "probs"))
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    fit <- suppressMessages(suppressWarnings(cs$fit()))
    fit_s <- strip_call_data(fit)
    expect_false(has_df_with(fit_s, "col_unused"), info = nm)
    if (isTRUE(cs$no_data)) {
      expect_identical(fit_s, fit, info = nm)
      next
    }
    expect_true(has_df_with(fit, "col_unused"), info = nm)
    expect_lt(length(serialize(fit_s, NULL)), length(serialize(fit, NULL)),
              label = nm)
    if (!is.null(cs$pred)) {
      expect_identical(cs$pred(fit_s), cs$pred(fit), info = nm)
    }
  }
})

test_that(paste(
  "the call of a stripped fit does not find other data, and lme4 predictions",
  "for the training data use the model frame"
), {
  fit <- suppressMessages(do.call(lme4::lmer, list(
    formula = .y.1 ~ x1 + x2 + (1 | g), data = dat_scd
  )))
  fit_s <- strip_call_data(fit)
  # A data frame named `data` must not replace the data of the fit:
  data <- dat_scd[1:10, ]
  expect_identical(predict(fit_s, re.form = ~ (1 | g)),
                   predict(fit, re.form = ~ (1 | g)))
  expect_error(update(fit_s, . ~ . - x2), "\\.projpred_data_removed")
})

test_that("strip_call_data() leaves other objects as they are", {
  expect_identical(strip_call_data(NULL), NULL)
  expect_identical(strip_call_data(dat_scd), dat_scd)
  x <- list(a = 1, call = quote(f(data = d)), b = list(call = quote(g())))
  expect_identical(strip_call_data(x), x)
})

test_that(paste(
  "divmin() gives multilevel and additive fits without the data in their",
  "calls"
), {
  fml_resp <- paste0("cbind(", paste0(".y.", seq_len(n_draws_scd),
                                      collapse = ", "), ")")
  divmin_scd <- function(rhs) {
    suppressWarnings(divmin(
      formula = as.formula(paste(fml_resp, "~", rhs)), data = dat_scd,
      family = extend_family(gaussian()), weights = rep(1, n_scd),
      projpred_var = matrix(NA, n_scd, n_draws_scd),
      projpred_ws_aug = matrix(NA, n_scd, n_draws_scd), verbose_divmin = FALSE
    ))
  }
  # The shared lme4 and mgcv setups, and the draw-wise lme4 fits:
  for (lmer_shared in c(TRUE, FALSE)) {
    withr::with_options(list(projpred.lmer_shared = lmer_shared), {
      for (rhs in c("x1 + x2 + (1 | g)", "s(x1) + x2")) {
        info <- paste(rhs, "lmer_shared:", lmer_shared)
        fits <- divmin_scd(rhs)
        expect_length(fits, n_draws_scd)
        expect_false(
          any(vapply(fits, has_df_with, logical(1), col = "col_unused")),
          info = info
        )
      }
    })
  }
})

test_that(paste(
  "gamm4_small_envs() removes the data from the environments of a gamm4 fit",
  "and keeps its predictions for new data"
), {
  # A formula environment without the data (the environment of a formula made
  # here has the environment of this file, with `dat_scd`, as an ancestor):
  fml <- as.formula(".y.1 ~ s(x1) + x2", env = new.env(parent = globalenv()))
  fit <- suppressWarnings(do.call(gamm4::gamm4, list(
    formula = fml, random = ~ (1 | g), data = dat_scd
  )))
  # The frame of gamm4() holds the data:
  expect_true(is.data.frame(get0("data", environment(fit$gam$pred.formula),
                                 inherits = FALSE)))
  fit_s <- gamm4_small_envs(fit)
  envs <- list(
    environment(fit_s$gam$pred.formula),
    attr(attr(fit_s$gam$model, "terms"), ".Environment"),
    attr(attr(fit_s$mer@frame, "terms"), ".Environment"),
    environment(attr(fit_s$mer@frame, "formula"))
  )
  for (e in envs) {
    expect_identical(e, environment(fml))
  }
  expect_lt(length(serialize(fit_s, NULL)), length(serialize(fit, NULL)))
  expect_identical(predict(fit_s$gam, newdata = newdat_scd),
                   predict(fit$gam, newdata = newdat_scd))
  expect_identical(predict(fit_s$gam, newdata = newdat_scd, type = "lpmatrix"),
                   predict(fit$gam, newdata = newdat_scd, type = "lpmatrix"))
})

test_that(paste(
  "divmin() gives gamm4 fits with the formula environment in their model",
  "frames and prediction formulas"
), {
  fml <- cbind(.y.1, .y.2, .y.3) ~ s(x1) + x2 + (1 | g)
  fits <- suppressWarnings(divmin(
    formula = fml, data = dat_scd, family = extend_family(gaussian()),
    weights = rep(1, n_scd), projpred_var = matrix(NA, n_scd, 3),
    projpred_ws_aug = matrix(NA, n_scd, 3), verbose_divmin = FALSE
  ))
  expect_length(fits, 3)
  for (fit in fits) {
    expect_s3_class(fit, "gamm4")
    expect_identical(environment(fit$gam$pred.formula), environment(fml))
    expect_identical(attr(attr(fit$mer@frame, "terms"), ".Environment"),
                     environment(fml))
  }
})

test_that("divmin_augdat() gives clmm() fits without the data in their calls", {
  # Without projpred's own Laplace fit, divmin_augdat() uses clmm():
  withr::local_options(projpred.mlvl_laplace = FALSE)
  cats <- c("lo", "mid", "hi")
  d_aug <- do.call(rbind, lapply(cats, function(cat_i) {
    d <- dat_scd[, c("x1", "x2", "g")]
    d$yo <- cat_i
    # A column that only the data of the fitter call holds:
    d$col_unused <- 1
    d
  }))
  d_aug$yo <- factor(d_aug$yo, levels = cats, ordered = TRUE)
  probs <- function(shift) {
    t(vapply(eta_scd + shift, function(eta_i) {
      diff(c(0, plogis(c(-0.5, 0.5) - eta_i), 1))
    }, numeric(3)))
  }
  ws_aug <- cbind(c(probs(0)), c(probs(0.3)))
  fam <- structure(list(family = "cumulative_rstanarm", link = "logit",
                        cats = cats), class = "family")
  fits <- suppressWarnings(divmin_augdat(
    formula = yo ~ x1 + (1 | g), data = d_aug, family = fam,
    weights = "unused", projpred_var = matrix(NA_real_, nrow(d_aug), 2),
    projpred_ws_aug = ws_aug, verbose_divmin = FALSE,
    throw_warn_sdivmin = FALSE, do_check_conv = FALSE
  ))
  expect_length(fits, 2)
  for (fit in fits) {
    expect_s3_class(fit, "clmm")
    expect_false(has_df_with(fit, "col_unused"))
  }
})

test_that(paste(
  "divmin() gives MASS::glmmPQL() fits without the data in their calls and",
  "formula environments"
), {
  withr::local_options(projpred.PQL = TRUE)
  d <- dat_scd[, c("x1", "x2", "g", "col_unused")]
  d$.p.1 <- plogis(eta_scd)
  d$.p.2 <- plogis(eta_scd + 0.3 * sin(1.7 * idx_scd))
  # One group-level term (MASS::glmmPQL() does not find the grouping variables
  # of the nlme::pdBlocked() structure that projpred uses for two terms).
  # A formula environment without the data (see the gamm4 test above):
  fml <- as.formula("cbind(.p.1, .p.2) ~ x1 + (1 | g)",
                    env = new.env(parent = globalenv()))
  fits <- suppressWarnings(divmin(
    formula = fml, data = d, family = extend_family(binomial()),
    weights = rep(1, n_scd), projpred_var = matrix(NA, n_scd, 2),
    projpred_ws_aug = matrix(NA, n_scd, 2), verbose_divmin = FALSE
  ))
  expect_length(fits, 2)
  for (fit in fits) {
    expect_s3_class(fit, "glmmPQL")
    expect_identical(environment(attr(fit$modelStruct$varStruct, "formula")),
                     environment(fml))
    # `fit$data` holds only the columns that the fit uses, so the name of the
    # unused column occurs only in a copy of the whole data:
    expect_length(grepRaw("col_unused", serialize(fit, NULL), fixed = TRUE),
                  0L)
  }
})
