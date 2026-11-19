# projpred's Laplace fit of cumulative-link multilevel submodels with one
# group-level intercept term (mlvl_laplace_cumul()). These tests need no
# reference model fits, so they also run on CRAN.

links_cumul_tst <- c("logit", "probit", "loglog", "cloglog", "cauchit")

cdf_cumul_tst <- function(z, link) {
  switch(link,
         logit = plogis(z),
         probit = pnorm(z),
         loglog = exp(-exp(-z)),
         cloglog = -expm1(-exp(z)),
         cauchit = pcauchy(z))
}

# Augmented data of a cumulative model with 4 categories: one row per
# observation and category (`y`), and the category probabilities as the
# weights (`w`). `d`: the data of the observations.
dat_cumul_tst <- function(link = "logit", n = 200, n_grp = 10, sd_grp = 0.7,
                          seed = 8) {
  withr::local_seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(paste0("g", sample(n_grp, n, replace = TRUE))))
  eta <- 0.8 * d$x - 0.5 * d$z + rnorm(n_grp, sd = sd_grp)[d$g]
  cum_probs <- cbind(vapply(c(-0.8, 0.2, 1.1), function(a) {
    cdf_cumul_tst(a - eta, link)
  }, numeric(n)), 1)
  probs <- cbind(cum_probs[, 1], cum_probs[, -1] - cum_probs[, -4])
  d_aug <- d[rep(seq_len(n), 4), ]
  d_aug$y <- factor(rep(1:4, each = n))
  list(d = d, d_aug = d_aug, w = as.vector(probs), probs = probs)
}

fam_cumul_tst <- function(link) {
  structure(list(family = "cumulative", link = link, cats = as.character(1:4)),
            class = "family")
}

fit_cumul_tst <- function(dat, link = "logit", formula = y ~ x + z + (1 | g),
                          data = dat$d_aug, weights = dat$w, ...) {
  fit_cumul_mlvl(formula = formula, data = data, family = fam_cumul_tst(link),
                 weights = weights, ...)
}

# The ordinal::clmm() fit. For the links cloglog and loglog, the default code
# of clmm() for a single group-level intercept term swaps the two links, so
# these fits use the other code of clmm() (`useMatrix = TRUE`).
# fit_cumul_mlvl() also sets `useMatrix = TRUE` for these links. The reference
# sets it itself, so that it does not depend on the code under test.
fit_clmm_tst <- function(dat, link = "logit", ...) {
  args <- list(dat, link = link, mlvl_laplace = FALSE, ...)
  if (link %in% c("cloglog", "loglog")) {
    args$control <- ordinal::clmm.control(useMatrix = TRUE)
  }
  suppressWarnings(suppressMessages(do.call(fit_cumul_tst, args)))
}

# The objective of mlvl_laplace_cumul() at the estimates of a clmm() fit:
obj_clmm_tst <- function(fit_clmm, dat, link) {
  des <- mlvl_cumul_design(y ~ x + z + (1 | g), dat$d_aug, dat$w)
  mlvl_cumul_obj(des, mlvl_cumul_link(link)$code, unname(fit_clmm$alpha),
                 unname(fit_clmm$beta[des$x_nms]),
                 unname(attr(ordinal::VarCorr(fit_clmm)$g, "stddev")))
}

test_that("mlvl_cumul_c() gives the gradient of its objective for each link", {
  dat <- dat_cumul_tst(n = 60, n_grp = 6)
  x <- as.matrix(dat$d_aug[, c("x", "z")])
  grp <- as.integer(dat$d_aug$g) - 1L
  y <- as.integer(dat$d_aug$y)
  # Other weights than those of the model, some of them 0:
  w <- dat$w * rep(c(0.5, 1, 2), length.out = length(dat$w))
  w[c(3, 70, 200)] <- 0
  par <- c(-0.9, 0.1, 1.2, 0.4, -0.3, 0.7)
  for (link in links_cumul_tst) {
    code <- mlvl_cumul_link(link)$code
    obj <- function(par) {
      mlvl_cumul_c(par[1:3], par[4:5], par[6], x, grp, 6L, y, w, rep(0, 6),
                   code, FALSE)$L
    }
    ev <- mlvl_cumul_c(par[1:3], par[4:5], par[6], x, grp, 6L, y, w,
                       rep(0, 6), code, TRUE)
    grad_fd <- vapply(seq_along(par), function(j) {
      h <- 1e-5
      e <- replace(numeric(length(par)), j, h)
      (obj(par + e) - obj(par - e)) / (2 * h)
    }, numeric(1))
    expect_equal(as.vector(ev$grad), grad_fd, tolerance = 1e-6, info = link)
    expect_equal(ev$L, obj(par), info = link)
    # Rows with weight 0 do not count:
    keep <- w > 0
    expect_equal(
      mlvl_cumul_c(par[1:3], par[4:5], par[6], x[keep, ], grp[keep], 6L,
                   y[keep], w[keep], rep(0, 6), code, FALSE)$L,
      ev$L, info = link
    )
  }
  # Thresholds that do not increase:
  expect_identical(
    mlvl_cumul_c(c(0.5, 0.1, 1.2), par[4:5], par[6], x, grp, 6L, y, w,
                 rep(0, 6), 1L, FALSE)$L,
    Inf
  )
})

test_that(paste(
  "mlvl_cumul_c() gives `Inf` or a finite objective and gradient for extreme",
  "parameter values"
), {
  dat <- dat_cumul_tst(n = 60, n_grp = 6)
  x <- as.matrix(dat$d_aug[, c("x", "z")])
  grp <- as.integer(dat$d_aug$g) - 1L
  y <- as.integer(dat$d_aug$y)
  # Thresholds far from the linear predictors, with probabilities that are
  # small but not 0 in double precision: a finite objective with the correct
  # gradient:
  for (link in links_cumul_tst) {
    # The distributions of loglog and cloglog have one short tail, where the
    # probabilities are 0 at a smaller distance:
    par <- if (link %in% c("loglog", "cloglog")) {
      c(-2.5, 0, 1.5, 0.5, -0.5, 1)
    } else {
      c(-5, 0, 6, 0.5, -2, 3)
    }
    code <- mlvl_cumul_link(link)$code
    obj <- function(par) {
      mlvl_cumul_c(par[1:3], par[4:5], par[6], x, grp, 6L, y, dat$w,
                   rep(0, 6), code, FALSE)$L
    }
    ev <- mlvl_cumul_c(par[1:3], par[4:5], par[6], x, grp, 6L, y, dat$w,
                       rep(0, 6), code, TRUE)
    expect_true(is.finite(ev$L), info = link)
    grad_fd <- vapply(seq_along(par), function(j) {
      h <- 1e-5
      e <- replace(numeric(length(par)), j, h)
      (obj(par + e) - obj(par - e)) / (2 * h)
    }, numeric(1))
    expect_equal(as.vector(ev$grad), grad_fd, tolerance = 1e-5, info = link)
  }
  alphas <- list(c(-30, -1, 2), c(-8, 0, 8), c(-1, 0, 800), c(-800, 0, 1),
                 c(-800, 0, 800), c(-1e-3, 0, 1e-3))
  for (link in links_cumul_tst) {
    code <- mlvl_cumul_link(link)$code
    for (alpha in alphas) {
      for (tau in c(0, 0.5, 30)) {
        ev <- mlvl_cumul_c(alpha, c(0.1, 20), tau, x, grp, 6L, y, dat$w,
                           rep(0, 6), code, TRUE)
        info_str <- paste(link, paste(alpha, collapse = " "), tau)
        expect_false(is.nan(ev$L), info = info_str)
        if (is.finite(ev$L)) {
          expect_true(all(is.finite(ev$grad)), info = info_str)
          expect_true(all(is.finite(ev$u)), info = info_str)
        } else {
          expect_identical(ev$L, Inf, info = info_str)
        }
      }
    }
  }
})

test_that(paste(
  "mlvl_cumul_c() reports the gradient at the conditional modes, and",
  "mlvl_cumul_optim() uses it for its convergence code"
), {
  skip_if_not_installed("ordinal")
  dat <- dat_cumul_tst(n = 120, n_grp = 6)
  des <- mlvl_cumul_design(y ~ x + z + (1 | g), dat$d_aug, dat$w)
  ev <- function(maxit) {
    mlvl_cumul_c(c(-1, 0, 1), c(0.5, -0.5), 2, des$x, des$grp, des$nlev,
                 des$y, des$w, rep(0, des$nlev), 1L, FALSE, maxit = maxit)
  }
  expect_lt(ev(50)$gmax, 1e-8)
  expect_gt(ev(1)$gmax, 1e-6)
  link <- mlvl_cumul_link("logit")
  expect_identical(mlvl_cumul_optim(des, link)$convergence, 0L)
  opt_tol0 <- mlvl_cumul_optim(des, link, tol_modes = 0)
  expect_identical(opt_tol0$convergence, 2L)
  expect_match(opt_tol0$message, "conditional modes did not converge")
  # fit_cumul_mlvl() then uses clmm():
  mlvl_cumul_c_orig <- mlvl_cumul_c
  local_mocked_bindings(mlvl_cumul_c = function(...) {
    out <- mlvl_cumul_c_orig(...)
    out$gmax <- 1
    out
  })
  fit <- suppressWarnings(suppressMessages(fit_cumul_tst(dat)))
  expect_s3_class(fit, "projpred_mlvl_cumul")
  expect_match(fit$optinfo$message, "ordinal::clmm() fit", fixed = TRUE)
  expect_match(fit$optinfo$message, "conditional modes did not converge")
})

test_that("mlvl_laplace_cumul() gives the clmm() fit for each link", {
  skip_if_not_installed("ordinal")
  for (link in links_cumul_tst) {
    dat <- dat_cumul_tst(link)
    fit <- fit_cumul_tst(dat, link)
    expect_s3_class(fit, c("projpred_mlvl_cumul", "projpred_mlvl"),
                    exact = TRUE)
    expect_identical(fit$link, link)
    expect_true(check_conv_s(fit), info = link)
    fit_c <- fit_clmm_tst(dat, link)
    expect_s3_class(fit_c, "clmm")
    # The objective is minus the log-likelihood of clmm() (the value that
    # clmm() reports differs slightly from the Laplace value at its
    # estimates):
    expect_equal(obj_clmm_tst(fit_c, dat, link), -as.numeric(logLik(fit_c)),
                 tolerance = 1e-4, info = link)
    # The same objective value, or a lower one:
    expect_lt(fit$deviance / 2, obj_clmm_tst(fit_c, dat, link) + 1e-6)
    expect_equal(fit$deviance / 2, obj_clmm_tst(fit_c, dat, link),
                 tolerance = 1e-6, info = link)
    # The parameters in the names of the clmm() path:
    for (nm_scheme in c("rstanarm", "brms")) {
      expect_equal(get_subparams(fit, nm_scheme = nm_scheme),
                   get_subparams(fit_c, nm_scheme = nm_scheme),
                   tolerance = 1e-3, info = paste(link, nm_scheme))
    }
    expect_equal(as.vector(subprd_augdat(list(fit), dat$d)),
                 as.vector(predict(fit, dat$d)), info = link)
    # subprd_augdat() cannot predict from clmm() fits with the loglog link
    # (linkfun_raw() does not have this link):
    if (link == "loglog") next
    # The predictions, with existing and new group levels. For new levels,
    # both paths use the same random numbers. The row names are those of the
    # clmm() path ("1", "2", ...), not those of the new data:
    nd <- dat$d[1:8, ]
    rownames(nd) <- paste0("row", 1:8)
    nd$g <- factor(c("g_new1", "g1", "g_new2", "g_new1", "g2", "g3", "g_new3",
                     "g_new2"))
    for (pred_new in c(FALSE, TRUE)) {
      withr::with_options(list(projpred.mlvl_pred_new = pred_new), {
        prd <- withr::with_seed(1, subprd_augdat(list(fit, fit), nd))
        prd_c <- withr::with_seed(1, subprd_augdat(list(fit_c, fit_c), nd))
      })
      expect_identical(dim(prd), c(8L, 3L, 2L))
      expect_identical(dimnames(prd)[[1]], as.character(1:8))
      expect_equal(prd, prd_c, tolerance = 1e-3,
                   info = paste(link, pred_new))
      # The two draws get different group-level intercepts for new levels:
      expect_false(isTRUE(all.equal(prd[1, , 1], prd[1, , 2])))
    }
  }
  # Link `probit_approx` of brms is fitted as `probit`:
  dat <- dat_cumul_tst("probit")
  expect_equal(fit_cumul_tst(dat, "probit_approx"),
               fit_cumul_tst(dat, "probit"))
})

test_that(paste(
  "an mlvl_laplace_cumul() fit without group-level variation equals the",
  "fit without the group-level term"
), {
  dat <- dat_cumul_tst(sd_grp = 0)
  fit <- fit_cumul_tst(dat)
  fit_polr <- suppressWarnings(MASS::polr(
    y ~ x + z, data = cbind(dat$d_aug, w = dat$w), weights = w
  ))
  expect_lt(fit$sd[["g"]], 1e-3)
  # A standard deviation near 0 starts a second search:
  expect_identical(fit$optinfo$n_restart, 1L)
  expect_equal(fit$coefficients, coef(fit_polr), tolerance = 1e-4)
  expect_equal(fit$alpha, fit_polr$zeta, tolerance = 1e-4)
  expect_equal(fit$deviance, deviance(fit_polr), tolerance = 1e-6)
})

test_that(paste(
  "observations with weight 0 do not change an mlvl_laplace_cumul() fit, and",
  "integer weights equal duplicated observations"
), {
  dat <- dat_cumul_tst()
  n <- nrow(dat$d)
  pars <- function(fit) c(fit$alpha, fit$coefficients, fit$sd)
  idxs_w0 <- 1:3
  probs_w0 <- dat$probs
  probs_w0[idxs_w0, ] <- 0
  d_chg <- dat$d_aug
  d_chg$x[d_chg$x %in% dat$d$x[idxs_w0]] <- 5
  fit_w0 <- fit_cumul_tst(dat, weights = as.vector(probs_w0))
  expect_s3_class(fit_w0, "projpred_mlvl_cumul")
  fit_chg <- fit_cumul_tst(dat, data = d_chg, weights = as.vector(probs_w0))
  keep <- rep(!seq_len(n) %in% idxs_w0, 4)
  fit_rm <- fit_cumul_tst(dat, data = dat$d_aug[keep, ],
                          weights = dat$w[keep])
  expect_equal(pars(fit_chg), pars(fit_w0), tolerance = 1e-10)
  expect_equal(pars(fit_rm), pars(fit_w0), tolerance = 1e-10)

  idxs_w2 <- 1:10
  probs_w2 <- dat$probs
  probs_w2[idxs_w2, ] <- 2 * probs_w2[idxs_w2, ]
  idxs_dupl <- rep(c(seq_len(n), idxs_w2), 4) +
    rep((0:3) * n, each = n + length(idxs_w2))
  fit_w2 <- fit_cumul_tst(dat, weights = as.vector(probs_w2))
  fit_dupl <- fit_cumul_tst(dat, data = dat$d_aug[idxs_dupl, ],
                            weights = dat$w[idxs_dupl])
  expect_equal(pars(fit_dupl), pars(fit_w2), tolerance = 1e-6)
})

test_that(paste(
  "fit_cumul_mlvl() uses clmm() for the cases that mlvl_laplace_cumul() does",
  "not handle"
), {
  skip_if_not_installed("ordinal")
  dat <- dat_cumul_tst(n = 120, n_grp = 6)
  expect_clmm <- function(fit) {
    expect_s3_class(fit, "clmm")
  }
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  # A group-level slope, and two group-level terms:
  expect_clmm(quiet(fit_cumul_tst(dat, formula = y ~ x + (x | g))))
  d2 <- dat$d_aug
  d2$h <- factor(paste0("h", as.integer(d2$g) %% 3))
  expect_clmm(quiet(fit_cumul_tst(dat, formula = y ~ x + (1 | g) + (1 | h),
                                  data = d2)))
  # Arguments for clmm() or its optimizers:
  expect_clmm(quiet(fit_cumul_tst(dat, method = "nlminb")))
  expect_clmm(quiet(fit_cumul_tst(dat, control = ordinal::clmm.control())))
  # An offset:
  d2$offs <- 0.1 * d2$z
  expect_clmm(quiet(fit_cumul_tst(dat, formula = y ~ x + offset(offs) + (1 | g),
                                  data = d2)))
  # A rank-deficient model matrix:
  d2$x2 <- 2 * d2$x
  expect_clmm(quiet(fit_cumul_tst(dat, formula = y ~ x + x2 + (1 | g),
                                  data = d2)))
  # The global option:
  withr::local_options(projpred.mlvl_laplace = FALSE)
  expect_clmm(quiet(fit_cumul_tst(dat)))
})

test_that("mlvl_cumul_design() returns `NULL` for the cases without a fit", {
  dat <- dat_cumul_tst(n = 40, n_grp = 4)
  fml <- y ~ x + (1 | g)
  expect_type(mlvl_cumul_design(fml, dat$d_aug, dat$w), "list")
  # A response category without weight:
  w_cat0 <- dat$w
  w_cat0[dat$d_aug$y == "2"] <- 0
  expect_null(mlvl_cumul_design(fml, dat$d_aug, w_cat0))
  # Missing values, a response that is not a factor, weights of another
  # length:
  d_na <- dat$d_aug
  d_na$x[1] <- NA
  expect_null(mlvl_cumul_design(fml, d_na, dat$w))
  d_num <- dat$d_aug
  d_num$y <- as.integer(d_num$y)
  expect_null(mlvl_cumul_design(fml, d_num, dat$w))
  expect_null(mlvl_cumul_design(fml, dat$d_aug, dat$w[-1]))
  # A link without a fit:
  expect_null(mlvl_laplace_cumul(fml, dat$d_aug, list(link = "softit"),
                                 dat$w))
})

test_that(paste(
  "fit_cumul_mlvl() refits with clmm() if mlvl_laplace_cumul() does not",
  "converge, and returns the clmm() fit as a `projpred_mlvl_cumul` object"
), {
  skip_if_not_installed("ordinal")
  local_mocked_bindings(mlvl_cumul_optim = function(des, link, ...) {
    list(alpha = c(-1, 0, 1), beta = rep(0, ncol(des$x)), tau = 1,
         u = rep(0, des$nlev), objective = 1, convergence = 1L,
         message = "mocked", iterations = 0L, n_restart = 0L)
  })
  # For cloglog and loglog, the clmm() fit must be that of the requested link:
  for (link in c("logit", "cloglog", "loglog")) {
    dat <- dat_cumul_tst(link, n = 120, n_grp = 6)
    fit <- suppressWarnings(suppressMessages(fit_cumul_tst(dat, link)))
    expect_s3_class(fit, c("projpred_mlvl_cumul", "projpred_mlvl"),
                    exact = TRUE)
    expect_identical(fit$link, link)
    expect_identical(fit$optinfo$message,
                     "ordinal::clmm() fit (Laplace fit: mocked)")
    expect_true(check_conv_s(fit), info = link)
    fit_c <- fit_clmm_tst(dat, link)
    expect_equal(get_subparams(fit, nm_scheme = "brms"),
                 get_subparams(fit_c, nm_scheme = "brms"), tolerance = 1e-8,
                 info = link)
    expect_equal(fit$deviance / 2, obj_clmm_tst(fit_c, dat, link), info = link)
    expect_equal(fit$deviance / 2, -as.numeric(logLik(fit_c)),
                 tolerance = 1e-4, info = link)
    # subprd_augdat() cannot predict from clmm() fits with the loglog link:
    if (link == "loglog") next
    nd <- dat$d[1:5, ]
    expect_equal(subprd_augdat(list(fit), nd), subprd_augdat(list(fit_c), nd),
                 tolerance = 1e-8, info = link)
  }
})

test_that(paste(
  "mlvl_laplace_cumul() gives the clmm() fit for factor predictors, a",
  "numeric grouping variable and a submodel without predictors"
), {
  skip_if_not_installed("ordinal")
  fmls <- list(y ~ x + f + (1 | g), y ~ x + (1 | g_num), y ~ 1 + (1 | g),
               y ~ (1 | g))
  cases <- expand.grid(link = c("logit", "cloglog"), fml = seq_along(fmls),
                       stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cases))) {
    link <- cases$link[i]
    fml <- fmls[[cases$fml[i]]]
    dat <- dat_cumul_tst(link, n = 150, n_grp = 6)
    dat$d_aug$f <- factor(c("a", "b", "c")[1 + as.integer(dat$d_aug$z > -0.5) +
                                             as.integer(dat$d_aug$z > 0.5)])
    dat$d_aug$g_num <- as.integer(dat$d_aug$g)
    info_str <- paste(link, paste(deparse(fml), collapse = " "))
    fit <- fit_cumul_tst(dat, link, formula = fml)
    expect_s3_class(fit, "projpred_mlvl_cumul")
    fit_c <- fit_clmm_tst(dat, link, formula = fml)
    expect_s3_class(fit_c, "clmm")
    for (nm_scheme in c("rstanarm", "brms")) {
      expect_equal(get_subparams(fit, nm_scheme = nm_scheme),
                   get_subparams(fit_c, nm_scheme = nm_scheme),
                   tolerance = 1e-3, info = paste(info_str, nm_scheme))
    }
    nd <- dat$d_aug[c(1, 5, 9), ]
    prd <- subprd_augdat(list(fit), nd)
    prd_c <- subprd_augdat(list(fit_c), nd)
    expect_identical(dimnames(prd), dimnames(prd_c), info = info_str)
    expect_true(all(is.finite(prd)), info = info_str)
    # The clmm() path gives `Inf` for a cumulative probability that is 1 in
    # double precision:
    is_fin <- is.finite(prd_c)
    expect_equal(prd[is_fin], prd_c[is_fin], tolerance = 1e-3,
                 info = info_str)
  }
})

test_that(paste(
  "divmin_augdat() gives `projpred_mlvl_cumul` fits, and subprd_augdat()",
  "handles them together with clmm() fits"
), {
  skip_if_not_installed("ordinal")
  dat <- dat_cumul_tst(n = 120, n_grp = 6)
  fit_divmin <- function() {
    divmin_augdat(formula = y ~ x + (1 | g), data = dat$d_aug,
                  family = fam_cumul_tst("logit"), weights = "unused",
                  projpred_var = matrix(NA_real_, nrow(dat$d_aug), 2),
                  projpred_ws_aug = cbind(dat$w, rev(dat$w)),
                  verbose_divmin = FALSE, throw_warn_sdivmin = TRUE,
                  do_check_conv = TRUE)
  }
  expect_no_warning(fits <- fit_divmin())
  expect_length(fits, 2)
  for (fit in fits) {
    expect_s3_class(fit, "projpred_mlvl_cumul")
  }
  withr::local_options(projpred.mlvl_laplace = FALSE)
  fits_c <- suppressWarnings(suppressMessages(fit_divmin()))
  expect_s3_class(fits_c[[1]], "clmm")
  nd <- dat$d[1:5, ]
  prd_mixed <- subprd_augdat(list(fits[[1]], fits_c[[2]]), nd)
  expect_equal(prd_mixed, subprd_augdat(fits_c, nd), tolerance = 1e-3)
  # A single row of new data:
  expect_identical(dim(subprd_augdat(fits, nd[1, ])), c(1L, 3L, 2L))
  expect_identical(dim(subprd_augdat(list(fits[[1]], fits_c[[2]]), nd[1, ])),
                   c(1L, 3L, 2L))
})

test_that("the methods for `projpred_mlvl_cumul` objects work", {
  dat <- dat_cumul_tst(n = 120, n_grp = 6)
  fit <- fit_cumul_tst(dat)
  expect_equal(formula(fit), y ~ x + z + (1 | g))
  prd <- predict(fit, dat$d[1:4, ])
  expect_identical(dimnames(prd), list(as.character(1:4), c("1|2", "2|3",
                                                            "3|4")))
  # The category probabilities of the fit are close to those of the data:
  cum_probs <- cbind(plogis(predict(fit, dat$d)), 1)
  probs <- cbind(cum_probs[, 1], cum_probs[, -1] - cum_probs[, -4])
  expect_lt(mean(abs(probs - dat$probs)), 0.05)
  expect_named(coef(fit), c("1|2", "2|3", "3|4", "x", "z"))
  expect_output(print(fit), "Thresholds")
  expect_error(as_glmerMod(fit, data = dat$d_aug),
               "does not support fits of the cumulative family")
})

test_that(paste(
  "cv_varsel() with LOO validated search and weight perturbation uses",
  "mlvl_laplace_cumul() with perturbed weights and gives the clmm() results"
), {
  skip_if_not_installed("ordinal")
  withr::local_seed(8)
  n <- 60
  S <- 20
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(paste0("g", sample(6, n, replace = TRUE))))
  b_g <- rnorm(6, sd = 0.8)
  eta <- 0.8 * d$x + b_g[d$g]
  d$y <- cut(eta + rlogis(n), c(-Inf, -0.8, 0.2, 1.1, Inf),
             labels = paste0("c", 1:4), ordered_result = TRUE)
  # Reference model draws (no MCMC):
  draws <- list(beta = rbind(rnorm(S, 0.8, 0.1), rnorm(S, 0, 0.1)),
                bg = matrix(rnorm(6 * S, rep(b_g, S), 0.2), 6, S),
                thres = cbind(rnorm(S, -0.8, 0.1), rnorm(S, 0.2, 0.1),
                              rnorm(S, 1.1, 0.1)))
  # An S x N x C_lat array of the thresholds minus the linear predictor:
  ref_predfun <- function(fit, newdata = NULL) {
    if (is.null(newdata)) newdata <- d
    eta <- cbind(newdata$x, newdata$z) %*% fit$beta +
      fit$bg[as.integer(factor(newdata$g, levels = levels(d$g))), ,
             drop = FALSE]
    arr <- vapply(1:3, function(cc) t(fit$thres[, cc] - t(eta)),
                  matrix(0, nrow(eta), S))
    aperm(arr, c(2, 1, 3))
  }
  fam <- extend_family(
    structure(list(family = "cumulative", link = "logit",
                   cats = paste0("c", 1:4)), class = "family"),
    augdat_y_unqs = paste0("c", 1:4), augdat_link = augdat_link_cumul,
    augdat_ilink = augdat_ilink_cumul, augdat_args_link = list(link = "logit"),
    augdat_args_ilink = list(link = "logit")
  )
  ref <- suppressWarnings(init_refmodel(
    object = draws, data = d, formula = y ~ x + z + (1 | g), family = fam,
    ref_predfun = ref_predfun
  ))
  run_cv <- function() {
    suppressWarnings(suppressMessages(cv_varsel(
      ref, cv_method = "LOO", validate_search = TRUE, method = "forward",
      nloo = 4, ndraws = 4, ndraws_pred = 8, nterms_max = 2, seed = 5,
      verbose = FALSE
    )))
  }
  mlvl_laplace_cumul_orig <- mlvl_laplace_cumul
  rec <- list()
  local_mocked_bindings(mlvl_laplace_cumul = function(formula, data, family,
                                                      weights = NULL) {
    fit <- mlvl_laplace_cumul_orig(formula, data, family, weights)
    rec[[length(rec) + 1L]] <<- list(w = weights, fit = fit)
    fit
  })
  cvvs <- run_cv()
  expect_identical(cvvs$loo_perturb_alpha, 10)
  expect_gt(length(rec), 0)
  # The weights of the augmented rows sum to the sum of the observation
  # weights; the perturbed searches give non-integer sums:
  w_sums <- vapply(rec, function(r) sum(r$w), numeric(1))
  w_pert <- abs(w_sums - round(w_sums)) > 1e-6
  expect_true(any(w_pert))
  expect_true(all(vapply(rec[w_pert], function(r) {
    inherits(r$fit, "projpred_mlvl_cumul") && r$fit$convergence == 0
  }, logical(1))))
  # The clmm() fits minimize the same objective:
  withr::local_options(projpred.mlvl_laplace = FALSE)
  cvvs_clmm <- run_cv()
  expect_identical(ranking(cvvs)$fulldata, ranking(cvvs_clmm)$fulldata)
  expect_identical(ranking(cvvs)$foldwise, ranking(cvvs_clmm)$foldwise)
  expect_equal(summary(cvvs)$perf_sub$elpd,
               summary(cvvs_clmm)$perf_sub$elpd, tolerance = 1e-3)
})
