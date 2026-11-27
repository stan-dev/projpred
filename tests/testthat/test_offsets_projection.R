# The projection with offsets: the submodel fits use the offsets, and their
# response is the reference means with offsets (the KL projection; see
# proj_offs_fit()). This file builds its own reference models and does not need
# setup.R.

# A reference model of family `family` with simulated draws of the coefficients
# (no model fit is necessary). `grp`: with group-level intercepts of the
# grouping variable `g`. `off`: the offsets (`NULL` for zeros). `smooth`: the
# formula has the smooth term `s(x1)` instead of `x1` (the draws do not
# change).
refm_offs <- function(family, ndraws = 3, n = 80, seed = 6751, grp = FALSE,
                      off = NULL, smooth = FALSE) {
  withr::local_seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                    g = factor(paste0("g", rep_len(1:8, n))))
  dat$off <- off %||% rep(0, n)
  beta <- c(0.3, 0.6, -0.4, 0.2) +
    matrix(rnorm(4 * ndraws, sd = 0.1), 4, ndraws)
  b_g <- if (grp) matrix(rnorm(8 * ndraws, sd = 0.5), 8, ndraws)
  eta_of <- function(nd, with_off = TRUE) {
    eta <- cbind(1, as.matrix(nd[c("x1", "x2", "x3")])) %*% beta
    if (grp) eta <- eta + b_g[as.integer(nd$g), , drop = FALSE]
    if (with_off) eta <- eta + nd$off
    eta
  }
  eta_mean <- rowMeans(eta_of(dat))
  dat$y <- switch(family$family,
                  poisson = rpois(n, family$linkinv(eta_mean)),
                  binomial = rbinom(n, 1, family$linkinv(eta_mean)),
                  gaussian = rnorm(n, family$linkinv(eta_mean)))
  fml <- if (smooth) y ~ x2 + x3 + s(x1) else y ~ x1 + x2 + x3
  if (grp) fml <- update(fml, . ~ . + (1 | g))
  fit <- list(beta = beta, data = dat)
  # The linear predictor with the offsets (init_refmodel() subtracts them):
  ref_predfun <- function(fit, newdata = NULL) {
    eta_of(if (is.null(newdata)) fit$data else newdata)
  }
  extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                 orhs = NULL, extract_y = TRUE) {
    y_wobs_offs(newdata = newdata %||% object$data, wrhs = wrhs,
                orhs = orhs %||% ~ off,
                resp_form = if (extract_y) ~ y else NULL)
  }
  # K-fold CV: the fold reference models have the draws of the full data (the
  # draws do not depend on the data) and the training data of the fold:
  make_refm <- function(fit_k, from_cv = FALSE) {
    suppressWarnings(init_refmodel(
      fit_k, data = fit_k$data, formula = fml, family = family,
      ref_predfun = ref_predfun, extract_model_data = extract_model_data,
      dis = if (family$family == "gaussian") rep(1, ndraws),
      cvfun = function(folds) rep(list(fit), max(folds)),
      cvrefbuilder = function(cvfit) {
        cvfit$data <- dat[-cvfit$omitted, , drop = FALSE]
        make_refm(cvfit, from_cv = TRUE)
      },
      called_from_cvrefbuilder = from_cv
    ))
  }
  refm <- make_refm(fit)
  list(refm = refm, dat = dat, off = dat$off,
       mu = family$linkinv(eta_of(dat, with_off = FALSE)),
       mu_offs = family$linkinv(eta_of(dat)))
}

# The coefficients of the KL projection of the means `mu` (with offsets) onto
# `fml_rhs`, by stats::glm() with the offsets `off`:
glm_kl <- function(fml_rhs, mu, dat, family, off) {
  quasi <- switch(family$family,
                  poisson = quasipoisson(link = family$link),
                  binomial = quasibinomial(link = family$link))
  d <- dat
  d$mu_ <- mu
  fit <- suppressWarnings(glm(update(fml_rhs, mu_ ~ .), data = d,
                              family = quasi, offset = off,
                              control = glm.control(epsilon = 1e-14,
                                                    maxit = 100)))
  unname(coef(fit))
}

# The coefficients of a `subfit` object:
coef_subfit <- function(fit) unname(c(fit$alpha, fit$beta))

test_that("GLM submodels with offsets give the KL projection", {
  withr::local_seed(6761)
  # Count the fits of the native multi-draw fitter:
  n_multi <- 0L
  multi_orig <- glm_ridge_native_multi
  local_mocked_bindings(glm_ridge_native_multi = function(...) {
    out <- multi_orig(...)
    if (!is.null(out)) n_multi <<- n_multi + 1L
    out
  })
  for (fam in list(poisson(), binomial())) {
    off <- if (fam$family == "poisson") {
      log(runif(80, 0.2, 5))
    } else {
      rnorm(80, sd = 0.7)
    }
    ref <- refm_offs(fam, off = off)
    info <- fam$family
    # The native multi-draw fitter and the draw-wise fitter (R callback):
    for (native in c(TRUE, FALSE)) {
      info_n <- paste(info, "native:", native)
      withr::local_options(projpred.glm_native = native)
      n_multi <- 0L
      prj <- project(ref$refm, predictor_terms = c("x1", "x2"), ndraws = 3,
                     regul = 0, thresh_conv = 1e-15, verbose = FALSE)
      expect_identical(n_multi > 0L, native, info = info_n)
      for (s in 1:3) {
        exp_s <- glm_kl(~ x1 + x2, ref$mu_offs[, s], ref$dat, fam, ref$off)
        expect_equal(coef_subfit(prj$outdmin[[s]]), exp_s, tolerance = 1e-7,
                     info = paste(info_n, "draw", s))
      }
      # The fit to the means without the offsets gives other coefficients:
      exp_no <- glm_kl(~ x1 + x2, ref$mu[, 1], ref$dat, fam, 0)
      expect_gt(max(abs(coef_subfit(prj$outdmin[[1]]) - exp_no)), 1e-3)
      # The predictions add the offsets once:
      pl <- proj_linpred(prj, transform = FALSE)
      x <- cbind(1, as.matrix(ref$dat[c("x1", "x2")]))
      expect_equal(unname(pl$pred[2, ]),
                   drop(x %*% coef_subfit(prj$outdmin[[2]])) + ref$off,
                   tolerance = 1e-12, info = info_n)
    }
  }
})

test_that("the forward search passes the offsets to the cache path", {
  withr::local_seed(6762)
  # The default option `projpred.full_size_ref` avoids the projection onto
  # the full model (see test_full_size_ref.R):
  withr::local_options(projpred.full_size_ref = FALSE)
  off <- log(runif(80, 0.2, 5))
  ref <- refm_offs(poisson(), off = off)
  # The fits of the cache path of proj_to_submodl() call divmin_multi() with
  # the responses `y` (divmin() calls it without `y`):
  n_offs <- 0L
  multi_orig <- divmin_multi
  local_mocked_bindings(
    divmin_multi = function(formulas, shared_x, dot_args, throw_warn_sdivmin,
                            do_check_conv, y = NULL, report = TRUE) {
      if (!is.null(y) && identical(dot_args[["offset"]], ref$off)) {
        n_offs <<- n_offs + 1L
      }
      multi_orig(formulas, shared_x, dot_args, throw_warn_sdivmin,
                 do_check_conv, y, report = report)
    }
  )
  vs <- varsel(ref$refm, method = "forward", nclusters = 1, nterms_max = 3,
               refit_prj = FALSE, regul = 0, thresh_conv = 1e-15,
               verbose = FALSE)
  expect_gt(n_offs, 0L)
  rk <- vs$predictor_ranking
  mu_c <- rowMeans(ref$mu_offs)
  for (k in 0:3) {
    rhs <- if (k == 0) ~ 1 else reformulate(rk[seq_len(k)])
    expect_equal(coef_subfit(vs$search_path$outdmins[[k + 1]][[1]]),
                 glm_kl(rhs, mu_c, ref$dat, poisson(), ref$off),
                 tolerance = 1e-7, info = paste("size", k))
  }
})

test_that("the L1 search and the score screen use the offsets", {
  withr::local_seed(6763)
  off <- log(runif(80, 0.2, 5))
  ref <- refm_offs(poisson(), off = off)
  mu_c <- rowMeans(ref$mu_offs)
  vs <- suppressWarnings(varsel(ref$refm, method = "L1", nterms_max = 3,
                                refit_prj = FALSE, verbose = FALSE))
  # The intercept-only submodel of the L1 path: the weighted means of the
  # reference model and of the submodel (with offsets) agree:
  expect_equal(vs$search_path$outdmins[[1]][[1]]$alpha,
               log(sum(mu_c) / sum(exp(ref$off))), tolerance = 1e-5)
  # The score screen gets the reference means with offsets, and the submodel
  # means with offsets fit them (the score of the intercept is 0):
  args_screen <- NULL
  screen_orig <- screen_gains
  local_mocked_bindings(screen_gains = function(...) {
    args_screen <<- list(...)
    screen_orig(...)
  })
  vs_scr <- varsel(ref$refm, method = "forward", nclusters = 1,
                   nterms_max = 2, search_screen = 1L, regul = 0,
                   thresh_conv = 1e-15, verbose = FALSE)
  expect_false(is.null(args_screen))
  expect_equal(drop(args_screen$mu_ref), unname(mu_c), tolerance = 1e-8)
  expect_equal(sum(args_screen$mu_ref - args_screen$mu_a), 0,
               tolerance = 1e-6)
})

test_that("multilevel submodels with offsets give the KL projection", {
  withr::local_seed(6764)
  off <- log(runif(80, 0.2, 5))
  ref <- refm_offs(poisson(), ndraws = 2, grp = TRUE, off = off)
  fam_g <- glmer_family(poisson())
  fits <- list()
  for (laplace in c(TRUE, FALSE)) {
    withr::local_options(projpred.mlvl_laplace = laplace)
    prj <- suppressWarnings(project(ref$refm,
                                    predictor_terms = c("x1", "(1 | g)"),
                                    ndraws = 2, verbose = FALSE))
    fits[[as.character(laplace)]] <- prj
    if (laplace) {
      expect_s3_class(prj$outdmin[[1]], "projpred_mlvl")
    } else {
      expect_s4_class(prj$outdmin[[1]], "glmerMod")
    }
  }
  for (s in 1:2) {
    d <- ref$dat
    d$mu_ <- ref$mu_offs[, s]
    fit_ref <- suppressWarnings(lme4::glmer(
      mu_ ~ x1 + (1 | g), data = d, family = fam_g, offset = ref$off,
      control = glmer_control(lme4::glmerControl(optimizer = "bobyqa"), fam_g)
    ))
    b_ref <- unname(lme4::fixef(fit_ref))
    re_ref <- lme4::ranef(fit_ref)$g[["(Intercept)"]]
    sub_laplace <- fits[["TRUE"]]$outdmin[[s]]
    sub_glmer <- fits[["FALSE"]]$outdmin[[s]]
    expect_equal(unname(sub_laplace$coefficients), b_ref, tolerance = 1e-4,
                 info = paste("draw", s))
    expect_equal(unname(lme4::fixef(sub_glmer)), b_ref, tolerance = 1e-6,
                 info = paste("draw", s))
    # The predictions add the offsets once:
    for (laplace in c("TRUE", "FALSE")) {
      sub <- fits[[laplace]]$outdmin[[s]]
      if (laplace == "TRUE") {
        b <- unname(sub$coefficients)
        re <- unname(sub$ranef$g[levels(d$g)])
      } else {
        b <- unname(lme4::fixef(sub))
        re <- lme4::ranef(sub)$g[levels(d$g), "(Intercept)"]
      }
      pl <- proj_linpred(fits[[laplace]], newdata = ref$dat,
                         transform = FALSE)
      expect_equal(unname(pl$pred[s, ]),
                   b[1] + b[2] * d$x1 + re[as.integer(d$g)] + ref$off,
                   tolerance = 1e-10, info = paste(laplace, "draw", s))
    }
  }
  # as_glmerMod() with the offsets: lme4 computes the same group-level
  # intercepts from the data and the offsets:
  sub_1 <- fits[["TRUE"]]$outdmin[[1]]
  d <- ref$dat
  d[[all.vars(formula(sub_1))[1]]] <- ref$mu_offs[, 1]
  expect_no_warning(fit_conv <- as_glmerMod(sub_1, data = d,
                                            offset = ref$off))
  expect_equal(lme4::ranef(fit_conv)$g[levels(d$g), "(Intercept)"],
               unname(sub_1$ranef$g[levels(d$g)]), tolerance = 1e-5)
  # Without the offsets, they differ:
  expect_warning(as_glmerMod(sub_1, data = d), "differs from `object`")
})

test_that("additive submodels with offsets give the KL projection", {
  withr::local_seed(6765)
  off <- log(runif(80, 0.2, 5))
  ref <- refm_offs(poisson(), ndraws = 2, off = off, smooth = TRUE)
  # Two draws (shared mgcv setup) and one cluster (fit_gam_callback()):
  for (ndr in c(2, 1)) {
    prj <- suppressWarnings(project(
      ref$refm, predictor_terms = c("x2", "s(x1)"),
      ndraws = 2, nclusters = if (ndr == 1) 1, verbose = FALSE
    ))
    mu_offs <- if (ndr == 1) rowMeans(ref$mu_offs) else ref$mu_offs
    mu_offs <- as.matrix(mu_offs)
    for (s in seq_len(ndr)) {
      d <- ref$dat
      d$mu_ <- mu_offs[, s]
      fit_ref <- suppressWarnings(mgcv::gam(mu_ ~ x2 + s(x1), data = d,
                                            family = poisson(),
                                            offset = ref$off))
      sub <- prj$outdmin[[s]]
      info <- paste("ndraws", ndr, "draw", s)
      expect_equal(unname(coef(sub)), unname(coef(fit_ref)), tolerance = 1e-6,
                   info = info)
      lp <- drop(predict(fit_ref, newdata = d, type = "lpmatrix") %*%
                   coef(fit_ref))
      pl <- proj_linpred(prj, newdata = ref$dat, transform = FALSE)
      pred_s <- if (ndr == 1) pl$pred else pl$pred[s, ]
      expect_equal(unname(drop(pred_s)), unname(lp + ref$off),
                   tolerance = 1e-6, info = info)
    }
  }
})

test_that("divmin() removes the offsets of observations with weight 0", {
  withr::local_seed(6766)
  off <- log(runif(80, 0.2, 5))
  ref <- refm_offs(poisson(), grp = TRUE, off = off)
  d <- ref$dat
  d$mu_ <- ref$mu_offs[, 1]
  w <- rep(1, 80)
  w[7] <- 0
  for (fml in list(mu_ ~ x1 + s(x2), mu_ ~ x1 + (1 | g))) {
    args <- list(formula = fml, family = poisson(), projpred_var =
                   matrix(NA_real_, 80, 1), projpred_ws_aug = NULL,
                 verbose_divmin = FALSE)
    out_w0 <- suppressWarnings(do.call(divmin, c(args, list(
      data = d, weights = w, offset = off
    ))))
    args$projpred_var <- matrix(NA_real_, 79, 1)
    out_rm <- suppressWarnings(do.call(divmin, c(args, list(
      data = d[-7, ], weights = w[-7], offset = off[-7]
    ))))
    expect_equal(coef(out_w0[[1]]), coef(out_rm[[1]]), tolerance = 1e-10,
                 info = deparse(fml))
  }
})

test_that("cv_varsel() with LOO and validated search runs with offsets", {
  withr::local_seed(6767)
  off <- log(runif(80, 0.2, 5))
  ref <- refm_offs(poisson(), grp = TRUE, off = off)
  # The Laplace fits get the offsets of their data rows (the data keep the
  # column `off`); the folds remove the left-out row:
  n_rows <- integer()
  offs_ok <- logical()
  laplace_orig <- mlvl_laplace
  local_mocked_bindings(mlvl_laplace = function(formula, data, family,
                                                weights = NULL,
                                                offset = NULL) {
    n_rows <<- c(n_rows, nrow(data))
    offs_ok <<- c(offs_ok, identical(offset, data$off))
    laplace_orig(formula, data, family, weights, offset)
  })
  cvvs <- suppressWarnings(cv_varsel(
    ref$refm, cv_method = "LOO", validate_search = TRUE, nloo = 3,
    nclusters = 1, nclusters_pred = 2, nterms_max = 2, verbose = FALSE,
    search_terms = c("x1", "(1 | g)", "x1 + (1 | g)")
  ))
  smmry <- summary(cvvs, stats = "elpd")
  expect_true(all(is.finite(smmry$perf_sub$elpd)))
  expect_true(length(offs_ok) > 0 && all(offs_ok))
  expect_true(any(n_rows == 79L))
})

# A regression guard (it passes also without the change):
test_that("the Gaussian family with the identity link ignores the offsets", {
  withr::local_seed(6768)
  ref <- refm_offs(gaussian(), ndraws = 2, off = rnorm(80))
  prj <- project(ref$refm, predictor_terms = c("x1", "x2"), ndraws = 2,
                 regul = 0, verbose = FALSE)
  for (s in 1:2) {
    d <- ref$dat
    d$mu_ <- ref$mu[, s]
    expect_equal(coef_subfit(prj$outdmin[[s]]),
                 unname(coef(lm(mu_ ~ x1 + x2, data = d))), tolerance = 1e-8,
                 info = paste("draw", s))
  }
})

test_that("unsupported cases with offsets stop with a message", {
  withr::local_seed(6769)
  off <- log(runif(80, 0.2, 5))
  ref <- refm_offs(poisson(), grp = TRUE, off = off)
  withr::with_options(list(projpred.glm_fitter = "fit_glm_callback"), {
    expect_error(project(ref$refm, predictor_terms = "x1", ndraws = 1,
                         verbose = FALSE),
                 "must be \"fit_glm_ridge_callback\"")
  })
  withr::with_options(list(projpred.PQL = TRUE), {
    expect_error(project(ref$refm, predictor_terms = c("x1", "(1 | g)"),
                         ndraws = 1, verbose = FALSE),
                 "`projpred.PQL` must be `FALSE`")
  })
  ref_s <- refm_offs(poisson(), grp = TRUE, off = off, smooth = TRUE)
  expect_error(suppressWarnings(project(
    ref_s$refm, predictor_terms = c("s(x1)", "(1 | g)"), ndraws = 1,
    verbose = FALSE
  )), "group-level and additive terms")
  expect_error(refm_offs(poisson(), off = c(-Inf, off[-1])),
               "only if all offsets are finite")
})

# A cumulative reference model (augmented-data projection, no Stan fit) with
# offsets and, for `grp = TRUE`, group-level intercepts.
refm_cumul_offs <- function(ndraws = 2, n = 60, seed = 6752, grp = FALSE) {
  withr::local_seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                    g = factor(paste0("g", rep_len(1:6, n))))
  dat$offs <- runif(n, -1, 1)
  thr <- c(-1, 0.3, 1.2)
  beta <- rbind(rnorm(ndraws, 1.2, 0.2), rnorm(ndraws, -0.8, 0.2))
  b_g <- if (grp) matrix(rnorm(6 * ndraws, sd = 0.6), 6, ndraws)
  thr_draws <- t(sapply(thr, function(t) t + rnorm(ndraws, 0, 0.05)))
  lp_of <- function(nd) {
    lp <- as.matrix(nd[, c("x1", "x2")]) %*% beta + nd$offs
    if (grp) lp <- lp + b_g[as.integer(nd$g), , drop = FALSE]
    lp
  }
  dat$y <- cut(rowMeans(lp_of(dat)) + rlogis(n), c(-Inf, thr, Inf),
               labels = paste0("c", 1:4), ordered_result = TRUE)
  fml <- if (grp) y ~ x1 + x2 + (1 | g) else y ~ x1 + x2
  # The latent predictors `thr_k - lp` (S x N x C_lat), offsets included:
  ref_predfun <- function(fit, newdata = NULL) {
    lp <- t(lp_of(if (is.null(newdata)) dat else newdata))
    simplify2array(lapply(seq_len(nrow(thr_draws)), function(k) {
      thr_draws[k, ] - lp
    }))
  }
  extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                 extract_y = TRUE) {
    y_wobs_offs(newdata = newdata, wrhs = wrhs, orhs = orhs %||% ~ offs,
                resp_form = if (extract_y) ~ y else NULL)
  }
  refm <- suppressWarnings(init_refmodel(
    list(), data = dat, formula = fml,
    family = structure(list(family = "cumulative_rstanarm", link = "logit"),
                       class = "family"),
    ref_predfun = ref_predfun, extract_model_data = extract_model_data,
    augdat_y_unqs = paste0("c", 1:4), augdat_link = augdat_link_cumul,
    augdat_ilink = augdat_ilink_cumul,
    augdat_args_link = list(link = "logit"),
    augdat_args_ilink = list(link = "logit")
  ))
  # The augmented data (one row per observation and category) and the
  # probabilities with offsets as their weights:
  d_aug <- do.call(rbind, lapply(paste0("c", 1:4), function(lv) {
    d <- dat
    d$y <- lv
    d
  }))
  d_aug$y <- factor(d_aug$y, levels = paste0("c", 1:4), ordered = TRUE)
  d_aug$offs <- rep(dat$offs, 4)
  list(refm = refm, dat = dat, d_aug = d_aug,
       ws_offs = unclass(refm$mu_offs), ws = unclass(refm$mu))
}

# The linear predictors alpha_k - x' beta - b_g - offset (N x C_lat) of a
# cumulative submodel fit:
lp_cumul <- function(alpha, beta, x, offs, b = 0) {
  outer(-(drop(x %*% beta) + b + offs), alpha, "+")
}

test_that("cumulative submodels with offsets give the KL projection", {
  ref <- refm_cumul_offs()
  x <- as.matrix(ref$dat["x1"])
  # The native multi-draw fitter and MASS::polr():
  for (native in c(TRUE, FALSE)) {
    withr::local_options(projpred.augdat_native = native)
    prj <- suppressWarnings(project(ref$refm, predictor_terms = "x1",
                                    ndraws = 2, verbose = FALSE))
    for (s in 1:2) {
      d <- ref$d_aug
      d$w_ <- ref$ws_offs[, s]
      fit_ref <- suppressWarnings(MASS::polr(
        y ~ x1 + offset(offs), data = d, weights = w_, method = "logistic",
        control = list(reltol = 1e-14)
      ))
      sub <- prj$outdmin[[s]]
      info <- paste("native:", native, "draw", s)
      # MASS::polr() (BFGS) in projpred has the default tolerance:
      expect_equal(unname(coef(sub)), unname(coef(fit_ref)), tolerance = 1e-4,
                   info = info)
      expect_equal(unname(sub$zeta), unname(fit_ref$zeta), tolerance = 1e-4,
                   info = info)
      # The predictions add the offsets once:
      lp <- ref$refm$mu_fun(prj$outdmin, newdata = ref$dat,
                            offset = ref$dat$offs, transform = FALSE)
      expect_equal(unname(unclass(lp)[, s]),
                   as.vector(lp_cumul(sub$zeta, coef(sub), x, ref$dat$offs)),
                   tolerance = 1e-8, info = info)
    }
    # The fit to the probabilities without offsets gives other coefficients:
    d <- ref$d_aug
    d$w_ <- ref$ws[, 1]
    fit_no <- suppressWarnings(MASS::polr(y ~ x1, data = d, weights = w_,
                                          method = "logistic"))
    expect_gt(max(abs(coef(prj$outdmin[[1]]) - coef(fit_no))), 1e-3)
  }
})

test_that(paste("multilevel cumulative submodels with offsets give the KL",
                "projection"), {
  skip_if_not_installed("ordinal")
  ref <- refm_cumul_offs(grp = TRUE)
  lvls <- levels(ref$dat$g)
  prj <- suppressWarnings(project(ref$refm,
                                  predictor_terms = c("x1", "(1 | g)"),
                                  ndraws = 2, nclusters = 1, verbose = FALSE))
  sub <- prj$outdmin[[1]]
  expect_s3_class(sub, "projpred_mlvl_cumul")
  # The reference: the objective of the Laplace fit minimized by optim().
  # test_mlvl_laplace_cumul.R checks this objective against ordinal::clmm()
  # (without offsets), and the check below the minimization checks its
  # offsets:
  d <- ref$d_aug
  w_ <- rowMeans(ref$ws_offs)
  des <- mlvl_cumul_design(y ~ x1 + (1 | g), d, w_, offset = d$offs)
  link <- mlvl_cumul_link("logit")
  # The objective is even in the standard deviation:
  obj <- function(par) {
    mlvl_cumul_obj(des, link$code, par[1:3], par[4], abs(par[5]))
  }
  opt <- optim(c(-1, 0.3, 1.2, 0, 1), obj,
               control = list(reltol = 1e-15, maxit = 20000))
  opt <- optim(opt$par, obj, control = list(reltol = 1e-15, maxit = 20000))
  opt$par[5] <- abs(opt$par[5])
  expect_equal(unname(c(sub$alpha, sub$coefficients, sub$sd)), opt$par,
               tolerance = 1e-3)
  expect_lte(sub$deviance / 2, opt$value + 1e-8)
  # The offsets in the objective: the same as a predictor column with the
  # coefficient 1:
  des_x <- mlvl_cumul_design(y ~ x1 + offs + (1 | g), d, w_)
  par_t <- c(-0.8, 0.2, 1.1, 0.9, 0.7)
  expect_equal(
    mlvl_cumul_obj(des, link$code, par_t[1:3], par_t[4], par_t[5]),
    mlvl_cumul_obj(des_x, link$code, par_t[1:3], c(par_t[4], 1), par_t[5]),
    tolerance = 1e-12
  )
  # The predictions add the offsets once:
  re <- unname(sub$ranef$g[lvls])
  lp <- ref$refm$mu_fun(prj$outdmin, newdata = ref$dat,
                        offset = ref$dat$offs, transform = FALSE)
  expect_equal(unname(unclass(lp)[, 1]),
               as.vector(lp_cumul(unname(sub$alpha), unname(sub$coefficients),
                                  as.matrix(ref$dat["x1"]), ref$dat$offs,
                                  b = re[as.integer(ref$dat$g)])),
               tolerance = 1e-10)
  # ordinal::clmm() (the Laplace fit switched off). The reference has the
  # offset term before the group-level term: clmm() takes the index of the
  # offset column from the fixed-effects terms (ordinal 2026.7.26).
  withr::local_options(projpred.mlvl_laplace = FALSE)
  prj_c <- suppressWarnings(project(ref$refm,
                                    predictor_terms = c("x1", "(1 | g)"),
                                    ndraws = 2, nclusters = 1,
                                    verbose = FALSE))
  sub_c <- prj_c$outdmin[[1]]
  expect_s3_class(sub_c, "clmm")
  d$w_ <- w_
  fit_ref <- suppressWarnings(ordinal::clmm(
    y ~ offset(offs) + x1 + (1 | g), data = d, weights = w_, link = "logit"
  ))
  expect_equal(unname(c(sub_c$alpha, sub_c$beta)),
               unname(c(fit_ref$alpha, fit_ref$beta)), tolerance = 1e-3)
  expect_equal(unname(c(sub_c$alpha, sub_c$beta)), opt$par[1:4],
               tolerance = 1e-3)
  re_c <- ordinal::ranef(sub_c)$g[lvls, "(Intercept)"]
  lp_c <- ref$refm$mu_fun(prj_c$outdmin, newdata = ref$dat,
                          offset = ref$dat$offs, transform = FALSE)
  expect_equal(unname(unclass(lp_c)[, 1]),
               as.vector(lp_cumul(unname(sub_c$alpha), unname(sub_c$beta),
                                  as.matrix(ref$dat["x1"]), ref$dat$offs,
                                  b = re_c[as.integer(ref$dat$g)])),
               tolerance = 1e-6)
})

test_that("divmin_augdat() gives binomial submodels with offsets", {
  withr::local_seed(6753)
  n <- 50
  dat <- data.frame(x1 = rnorm(n), offs = rnorm(n, sd = 0.7))
  p1 <- plogis(0.2 + 0.7 * dat$x1 + dat$offs + rnorm(n, sd = 0.3))
  d_aug <- rbind(dat, dat)
  d_aug$y <- factor(rep(c("0", "1"), each = n), levels = c("0", "1"))
  out <- divmin_augdat(y ~ x1, data = d_aug, family = binomial(),
                       weights = rep(1, n),
                       projpred_var = matrix(NA_real_, 2 * n, 1),
                       projpred_ws_aug = matrix(c(1 - p1, p1)),
                       verbose_divmin = FALSE, offset = dat$offs,
                       regul = 0, thresh_conv = 1e-15)
  d <- dat
  d$p1 <- p1
  fit_ref <- suppressWarnings(glm(p1 ~ x1, data = d, family = quasibinomial(),
                                  offset = offs,
                                  control = glm.control(epsilon = 1e-14)))
  expect_equal(coef_subfit(out[[1]]), unname(coef(fit_ref)), tolerance = 1e-7)
  # The categorical family stops:
  expect_error(divmin_augdat(
    y ~ x1, data = d_aug,
    family = structure(list(family = "categorical", link = "logit"),
                       class = "family"),
    weights = rep(1, n), projpred_var = matrix(NA_real_, 2 * n, 1),
    projpred_ws_aug = matrix(c(1 - p1, p1)), verbose_divmin = FALSE,
    offset = dat$offs
  ), "brms::categorical")
})

test_that("cv_varsel() with K-fold CV gives the fold projections the offsets", {
  withr::local_seed(6770)
  off <- log(runif(80, 0.2, 5))
  ref <- refm_offs(poisson(), off = off)
  # The fits of the cache path: each response row (the mean over the draws,
  # one cluster) has its own offset; the training data of the folds have fewer
  # rows:
  mu_all <- rowMeans(ref$mu_offs)
  n_rows <- integer()
  offs_ok <- logical()
  multi_orig <- glm_ridge_native_multi
  local_mocked_bindings(
    glm_ridge_native_multi = function(formulas, shared_x, dot_args,
                                      y = NULL) {
      if (!is.null(y) && NCOL(y) == 1) {
        y_v <- as.vector(y)
        idx <- vapply(y_v, function(v) which.min(abs(mu_all - v)), integer(1))
        n_rows <<- c(n_rows, length(y_v))
        offs_ok <<- c(offs_ok,
                      isTRUE(all.equal(y_v, unname(mu_all[idx]),
                                       tolerance = 1e-8)) &&
                        isTRUE(all.equal(unname(dot_args[["offset"]]),
                                         ref$off[idx], tolerance = 0)))
      }
      multi_orig(formulas, shared_x, dot_args, y)
    }
  )
  cvvs <- suppressWarnings(cv_varsel(
    ref$refm, cv_method = "kfold", K = 2, nclusters = 1, nclusters_pred = 1,
    nterms_max = 2, verbose = FALSE
  ))
  expect_true(all(is.finite(summary(cvvs, stats = "elpd")$perf_sub$elpd)))
  expect_true(length(offs_ok) > 0 && all(offs_ok))
  expect_true(any(n_rows < 80L))
})

test_that("the parallel path of divmin_augdat() gets the offsets", {
  ref <- refm_cumul_offs()
  withr::local_options(projpred.augdat_native = FALSE)
  prj_seq <- suppressWarnings(project(ref$refm, predictor_terms = "x1",
                                      ndraws = 2, verbose = FALSE))
  withr::local_options(projpred.parallel_proj_trigger = 1)
  n_par <- 0
  prj_par <- with_mocked_bindings(
    suppressWarnings(project(ref$refm, predictor_terms = "x1", ndraws = 2,
                             verbose = FALSE)),
    par_backend = function(...) "future",
    par_lapply = function(.X, .FUN, ..., .backend, .rng = FALSE,
                          .progress = FALSE) {
      n_par <<- n_par + 1
      par_chunk(list(xs = .X, seeds = NULL),
                list(fun = .FUN, args = list(...), opts = list(),
                     globs = list()))
    }
  )
  expect_gt(n_par, 0)
  for (s in 1:2) {
    expect_equal(coef(prj_par$outdmin[[s]]), coef(prj_seq$outdmin[[s]]),
                 tolerance = 1e-12, info = paste("draw", s))
    expect_equal(prj_par$outdmin[[s]]$zeta, prj_seq$outdmin[[s]]$zeta,
                 tolerance = 1e-12, info = paste("draw", s))
  }
  # The KL projection with offsets (as in the test of the cumulative
  # submodels above):
  d <- ref$d_aug
  d$w_ <- ref$ws_offs[, 1]
  fit_ref <- suppressWarnings(MASS::polr(
    y ~ x1 + offset(offs), data = d, weights = w_, method = "logistic",
    control = list(reltol = 1e-14)
  ))
  expect_equal(unname(coef(prj_par$outdmin[[1]])), unname(coef(fit_ref)),
               tolerance = 1e-4)
})

test_that("a replaced divergence minimizer gets no offsets", {
  withr::local_seed(6772)
  ref <- refm_offs(poisson(), off = log(runif(80, 0.2, 5)))
  refm <- ref$refm
  expect_true(proj_offs_fit(refm))
  n_calls <- 0L
  got_offs <- logical()
  refm$div_minimizer <- function(...) {
    n_calls <<- n_calls + 1L
    got_offs <<- c(got_offs, !is.null(list(...)[["offset"]]))
    divmin(...)
  }
  expect_false(proj_offs_fit(refm))
  prj <- project(refm, predictor_terms = "x1", ndraws = 2, verbose = FALSE)
  expect_gt(n_calls, 0L)
  expect_false(any(got_offs))
})

test_that("offs_fit_family() excludes only Gaussian identity and latent", {
  expect_true(offs_fit_family(extend_family(poisson())))
  expect_true(offs_fit_family(extend_family(binomial())))
  expect_true(offs_fit_family(extend_family(gaussian(link = "log"))))
  expect_true(offs_fit_family(extend_family(Gamma(link = "log"))))
  expect_true(offs_fit_family(list(family = "cumulative", link = "logit",
                                   for_augdat = TRUE)))
  expect_false(offs_fit_family(extend_family(gaussian())))
  expect_false(offs_fit_family(list(family = "binomial", link = "logit",
                                    for_latent = TRUE)))
})

