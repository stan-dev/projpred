# Contrasts of factor predictors: the submodels use the `contrasts` attributes
# of the reference model's data.

# Data with a 3-level factor `f`, an ordered 4-level factor `o`, a numeric
# predictor `x` and a grouping variable `g`. `ctr_f` and `ctr_o` are the
# `contrasts` attributes of `f` and `o` (`NULL` for the default contrasts).
dat_contr <- function(n = 60, seed = 4417, ctr_f = NULL, ctr_o = NULL) {
  set.seed(seed)
  dat <- data.frame(
    f = factor(rep_len(c("a", "b", "c"), n)),
    o = factor(rep_len(c("lo", "mid", "hi", "top"), n),
               levels = c("lo", "mid", "hi", "top"), ordered = TRUE),
    x = rnorm(n),
    g = factor(rep_len(paste0("g", 1:6), n))
  )
  dat$y <- c(a = 0, b = 0.8, c = -0.6)[as.character(dat$f)] +
    0.3 * as.integer(dat$o) + 0.5 * dat$x +
    rnorm(6, sd = 0.5)[as.integer(dat$g)] + rnorm(n)
  if (!is.null(ctr_f)) contrasts(dat$f) <- ctr_f
  if (!is.null(ctr_o)) contrasts(dat$o) <- ctr_o
  return(dat)
}

# A Gaussian reference model with simulated draws of the coefficients (no model
# fit is necessary) of the linear model `fml_ref`. Its predictions use the
# default contrasts; they do not depend on the contrasts of `dat`. `fml` is the
# full-model formula for the search.
refm_contr <- function(dat, fml = y ~ f + o + x, fml_ref = y ~ f + o + x,
                       ndraws = 100, seed = 8273) {
  # dat_contr() sets the seed, so evaluate `dat` before set.seed():
  dat_dflt <- strip_contrasts(dat)
  set.seed(seed)
  tt <- delete.response(terms(fml_ref))
  mm <- model.matrix(tt, dat_dflt)
  beta_hat <- drop(solve(crossprod(mm), crossprod(mm, dat$y)))
  beta <- beta_hat + matrix(rnorm(ncol(mm) * ndraws, sd = 0.05), ncol(mm))
  xlev <- .getXlevels(tt, model.frame(tt, dat_dflt))
  fit <- list(beta = beta, data = dat_dflt, tt = tt, xlev = xlev)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- strip_contrasts(if (is.null(newdata)) fit$data else newdata)
    model.matrix(fit$tt, model.frame(fit$tt, nd, xlev = fit$xlev)) %*% fit$beta
  }
  suppressWarnings(init_refmodel(fit, data = dat, formula = fml,
                                 family = gaussian(), ref_predfun = ref_predfun,
                                 dis = rep(1, ndraws)))
}

# Evaluates `expr` and counts the calls of the internal function `fun_nm` that
# returned a result other than `NULL` (a fast path that applied). Returns the
# value of `expr` and the count.
with_return_count <- function(fun_nm, expr) {
  calls <- new.env()
  calls$n <- 0L
  ns <- asNamespace("projpred")
  suppressMessages(trace(
    fun_nm, where = ns, print = FALSE,
    exit = bquote(if (!is.null(returnValue(NULL))) {
      assign("n", get("n", envir = .(calls)) + 1L, envir = .(calls))
    })
  ))
  on.exit(suppressMessages(untrace(fun_nm, where = ns)), add = TRUE)
  value <- expr
  return(list(value = value, n = calls$n))
}

# A custom contrast matrix for `o` (with column names):
ctr_o_cust <- matrix(c(0, 1, 1, 1,
                       0, 0, 1, 1,
                       0, 0, 0, 1), ncol = 3,
                     dimnames = list(NULL, c("_ge_mid", "_ge_hi", "_ge_top")))

test_that("init_refmodel() stores the contrasts of the predictor columns", {
  dat <- dat_contr(ctr_f = contr.sum(3), ctr_o = ctr_o_cust)
  contrasts(dat$g) <- contr.sum(6)
  refm <- refm_contr(dat)
  expect_named(refm$contrasts, c("f", "o"))
  expect_equal(unname(refm$contrasts$f), unname(contr.sum(3)))
  expect_identical(rownames(refm$contrasts$f), levels(dat$f))
  expect_equal(unname(refm$contrasts$o), unname(ctr_o_cust))
  expect_identical(colnames(refm$contrasts$o), colnames(ctr_o_cust))
  # The stored data have no `contrasts` attributes:
  for (nm in names(refm$fetch_data())) {
    expect_null(attr(refm$fetch_data()[[nm]], "contrasts"), info = nm)
  }
  # A character name of a contrast function gives its matrix, so a later change
  # of `options("contrasts")` does not change it:
  dat_nm <- dat_contr()
  contrasts(dat_nm$f) <- "contr.sum"
  refm_nm <- refm_contr(dat_nm)
  opts_old <- options(contrasts = c("contr.helmert", "contr.poly"))
  on.exit(options(opts_old), add = TRUE)
  expect_equal(unname(refm_nm$contrasts$f), unname(contr.sum(3)))
  options(opts_old)
  # Without attributes:
  expect_null(refm_contr(dat_contr())$contrasts)
})

test_that("set_contrasts() keeps the default contrasts for other levels", {
  dat <- dat_contr(ctr_f = contr.sum(3))
  ctrs <- contrasts_from_data(dat, vars = "f")
  dat_dflt <- strip_contrasts(dat)
  expect_equal(attr(set_contrasts(dat_dflt, ctrs)$f, "contrasts"), ctrs$f)
  # A level that no row has:
  dat_sub <- dat_dflt[dat_dflt$f != "c", ]
  expect_null(attr(set_contrasts(dat_sub, ctrs)$f, "contrasts"))
  # Other levels:
  dat_lvl <- dat_dflt
  dat_lvl$f <- factor(dat_lvl$f, levels = c("c", "b", "a"))
  expect_null(attr(set_contrasts(dat_lvl, ctrs)$f, "contrasts"))
  expect_identical(set_contrasts(dat_dflt, NULL), dat_dflt)
})

test_that(paste(
  "the projection uses the contrasts of the reference model, and the",
  "predictions do not depend on them"
), {
  dat_sum <- dat_contr(ctr_f = contr.sum(3), ctr_o = ctr_o_cust)
  refm_sum <- refm_contr(dat_sum)
  refm_dflt <- refm_contr(dat_contr())
  prj_sum <- project(refm_sum, predictor_terms = c("f", "o", "x"),
                     nclusters = 1, verbose = FALSE)
  prj_dflt <- project(refm_dflt, predictor_terms = c("f", "o", "x"),
                      nclusters = 1, verbose = FALSE)
  beta_sum <- prj_sum$outdmin[[1]]$beta
  expect_identical(rownames(beta_sum),
                   c("f1", "f2", "o_ge_mid", "o_ge_hi", "o_ge_top", "x"))
  expect_identical(rownames(prj_dflt$outdmin[[1]]$beta),
                   c("fb", "fc", "o.L", "o.Q", "o.C", "x"))
  # The coefficients agree with a least-squares fit to the mean of the
  # reference model's predictions (`regul` is small):
  mu_ref <- rowMeans(refm_sum$mu)
  dat_lm <- cbind(dat_sum, mu_ref = mu_ref)
  coef_lm <- coef(lm(mu_ref ~ f + o + x, data = dat_lm))
  expect_equal(unname(drop(beta_sum)), unname(coef_lm[-1]), tolerance = 1e-3)
  expect_equal(prj_sum$outdmin[[1]]$alpha, unname(coef_lm[1]),
               tolerance = 1e-3)
  # Predictions for new data without `contrasts` attributes, for new data with
  # them (without a warning), and for the training data:
  nd <- strip_contrasts(dat_sum[c(1, 2, 3, 8), ])
  nd_attr <- dat_sum[c(1, 2, 3, 8), ]
  pred_sum <- proj_linpred(prj_sum, newdata = nd)$pred
  expect_no_warning(
    pred_sum_attr <- proj_linpred(prj_sum, newdata = nd_attr)$pred
  )
  expect_equal(pred_sum_attr, pred_sum, tolerance = 1e-12)
  expect_equal(unname(pred_sum),
               unname(proj_linpred(prj_sum)$pred[, c(1, 2, 3, 8),
                                                 drop = FALSE]),
               tolerance = 1e-10)
  expect_equal(pred_sum, proj_linpred(prj_dflt, newdata = nd)$pred,
               tolerance = 1e-4)
})

test_that(paste(
  "the cross-product search and the draw-wise search use the same contrasts"
), {
  refm <- refm_contr(dat_contr(ctr_f = contr.sum(3), ctr_o = ctr_o_cust))
  des <- gauss_crossprod_design(refm, c("1", "f", "o", "x"), 1e-4)
  expect_equal(des$contrasts$f, refm$contrasts$f)
  args_vs <- list(refm, nclusters = 3, nclusters_pred = 5, nterms_max = 3,
                  seed = 3317, verbose = FALSE)
  res_crossprod <- with_return_count(
    "search_forward_gauss", suppressWarnings(do.call(varsel, args_vs))
  )
  expect_identical(res_crossprod$n, 1L)
  vs_crossprod <- res_crossprod$value
  opts_old <- options(projpred.search_gauss_crossprod = FALSE)
  on.exit(options(opts_old), add = TRUE)
  vs_drawwise <- with_mocked_bindings(
    suppressWarnings(do.call(varsel, args_vs)),
    glm_ridge_gauss_multi = function(...) NULL
  )
  options(opts_old)
  expect_identical(ranking(vs_crossprod)$fulldata,
                   ranking(vs_drawwise)$fulldata)
  expect_equal(performances(vs_crossprod, deltas = FALSE),
               performances(vs_drawwise, deltas = FALSE))
  expect_equal(vs_crossprod$search_path$outdmins,
               vs_drawwise$search_path$outdmins)
  # The LOO CV with validated search (perf_eval_gauss() for the test data):
  args_cv <- list(refm, cv_method = "LOO", validate_search = TRUE, nloo = 5,
                  nclusters = 3, ndraws_pred = 20, nterms_max = 3,
                  seed = 5821, thresh_conv = 1e-15, verbose = FALSE)
  res_cv <- with_return_count(
    "perf_eval_gauss", suppressWarnings(do.call(cv_varsel, args_cv))
  )
  expect_identical(res_cv$n, 5L)
  cv_crossprod <- res_cv$value
  opts_old <- options(projpred.search_gauss_crossprod = FALSE)
  cv_drawwise <- with_mocked_bindings(
    suppressWarnings(do.call(cv_varsel, args_cv)),
    glm_ridge_gauss_multi = function(...) NULL
  )
  options(opts_old)
  expect_equal(cv_crossprod$summaries$sub, cv_drawwise$summaries$sub,
               tolerance = 1e-10)
})

test_that("the L1 search uses the contrasts of the reference model", {
  refm <- refm_contr(dat_contr(ctr_f = contr.sum(3), ctr_o = ctr_o_cust))
  vs <- suppressWarnings(varsel(refm, method = "L1", refit_prj = FALSE,
                                nterms_max = 3, seed = 3906, verbose = FALSE))
  expect_setequal(ranking(vs)$fulldata, c("f", "o", "x"))
  outdmin_full <- vs$search_path$outdmins[[4]][[1]]
  expect_setequal(rownames(outdmin_full$beta),
                  c("f1", "f2", "o_ge_mid", "o_ge_hi", "o_ge_top", "x"))
  prj <- suppressWarnings(project(vs, nterms = 3, refit_prj = FALSE,
                                  verbose = FALSE))
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = refm$fetch_data()[1:10, ])$pred
  expect_equal(unname(drop(pred_new)), unname(drop(pred_train)[1:10]),
               tolerance = 1e-10)
})

test_that("multilevel and additive submodels use the contrasts", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("mgcv")
  dat <- dat_contr(ctr_f = contr.sum(3))
  refm <- refm_contr(dat, fml = y ~ f + x + (1 | g))
  prj_glmm <- suppressWarnings(project(
    refm, predictor_terms = c("f", "(1 | g)"), nclusters = 1, verbose = FALSE
  ))
  expect_identical(names(lme4::fixef(prj_glmm$outdmin[[1]])),
                   c("(Intercept)", "f1", "f2"))
  nd <- strip_contrasts(dat[1:5, ])
  pred_glmm <- proj_linpred(prj_glmm, newdata = nd)$pred
  expect_equal(unname(pred_glmm),
               unname(proj_linpred(prj_glmm)$pred[, 1:5, drop = FALSE]),
               tolerance = 1e-10)

  refm_gam <- refm_contr(dat, fml = y ~ f + s(x))
  prj_gam <- suppressWarnings(project(
    refm_gam, predictor_terms = c("f", "s(x)"), nclusters = 1,
    verbose = FALSE
  ))
  expect_identical(names(coef(prj_gam$outdmin[[1]]))[1:3],
                   c("(Intercept)", "f1", "f2"))
  pred_gam <- proj_linpred(prj_gam, newdata = nd)$pred
  expect_equal(unname(pred_gam),
               unname(proj_linpred(prj_gam)$pred[, 1:5, drop = FALSE]),
               tolerance = 1e-10)
})

test_that(paste(
  "the LOO CV with validated search works with contrasts and a level that",
  "only one observation has"
), {
  dat <- dat_contr(n = 30)
  dat$f <- factor(c("d", as.character(dat$f[-1])))
  contrasts(dat$f) <- contr.sum(4)
  # The GLM path keeps the left-out row with weight 0; the multilevel path
  # removes it and drops the lost level (drop_lost_levels()):
  for (fml in list(y ~ f + x, y ~ f + x + (1 | g))) {
    refm <- refm_contr(dat, fml = fml, fml_ref = y ~ f + x)
    # All observations are left out once, also the one with level "d". The
    # lme4 fits of these small data give convergence warnings; collect the
    # warnings and check that none is about contrasts:
    warns <- character()
    cv <- withCallingHandlers(
      suppressMessages(cv_varsel(
        refm, cv_method = "LOO", validate_search = TRUE, nclusters = 2,
        ndraws_pred = 10, nterms_max = 2, seed = 1, verbose = FALSE
      )),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_false(any(grepl("contrast", warns)), info = deparse(fml))
    expect_true(all(is.finite(
      performances(cv, deltas = FALSE)$submodels$elpd
    )), info = deparse(fml))
  }
})

test_that("the datafit cvrefbuilder() keeps the contrasts", {
  dat <- dat_contr(ctr_f = contr.sum(3))
  # init_refmodel() warns that it uses the internal default `cvrefbuilder`:
  dfit <- suppressWarnings(init_refmodel(object = NULL, data = dat,
                                         formula = y ~ f + x,
                                         family = gaussian()))
  expect_equal(dfit$contrasts, contrasts_from_data(dat, vars = "f"))
  dfit_fold <- suppressWarnings(dfit$cvrefbuilder(list(omitted = 1:5)))
  expect_equal(dfit_fold$contrasts, dfit$contrasts)
})

test_that("get_refmodel() takes the contrasts of an rstanarm fit", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  dat <- dat_contr()
  fit <- suppressWarnings(rstanarm::stan_glm(
    y ~ f + x, family = gaussian(), data = dat,
    contrasts = list(f = "contr.sum"), chains = 1, iter = 500, refresh = 0,
    seed = 7193
  ))
  refm <- get_refmodel(fit)
  expect_named(refm$contrasts, "f")
  expect_equal(unname(refm$contrasts$f), unname(contr.sum(3)))
  expect_null(attr(refm$fetch_data()$f, "contrasts"))
  prj <- project(refm, predictor_terms = c("f", "x"), nclusters = 1,
                 verbose = FALSE)
  expect_identical(rownames(prj$outdmin[[1]]$beta), c("f1", "f2", "x"))
  # The reference model's coefficient names have the same coding:
  expect_true(all(c("f1", "f2") %in% names(coef(fit))))
})

test_that("the augmented-data submodels use the contrasts", {
  cats <- c("1", "2", "3")
  dat <- dat_contr(ctr_f = contr.sum(3))
  ctrs <- contrasts_from_data(dat, vars = "f")
  eta <- c(a = 0, b = 0.8, c = -0.6)[as.character(dat$f)] + 0.5 * dat$x
  fams <- list(
    cumul = structure(list(family = "cumulative_rstanarm", link = "logit",
                           cats = cats), class = "family"),
    categ = structure(list(family = "categorical", link = "logit",
                           cats = cats), class = "family")
  )
  probs <- list(
    cumul = t(vapply(eta, function(eta_i) {
      diff(c(0, plogis(c(-0.5, 0.5) - eta_i), 1))
    }, numeric(3))),
    categ = exp(cbind(0, eta, -eta)) / rowSums(exp(cbind(0, eta, -eta)))
  )
  # The augmented data as in proj_to_submodl() (one projected draw):
  fit_aug <- function(fam_nm, ctrs) {
    d_aug <- do.call(rbind, lapply(cats, function(cat_i) {
      d <- strip_contrasts(dat)
      d$y <- cat_i
      d
    }))
    d_aug$y <- factor(d_aug$y, levels = cats)
    d_aug <- set_contrasts(d_aug, ctrs)
    divmin_augdat(formula = y ~ f + x, data = d_aug, family = fams[[fam_nm]],
                  weights = "unused",
                  projpred_var = matrix(NA_real_, nrow(d_aug), 1),
                  projpred_ws_aug = matrix(c(probs[[fam_nm]])),
                  verbose_divmin = FALSE, throw_warn_sdivmin = FALSE,
                  do_check_conv = FALSE)
  }
  nd <- strip_contrasts(dat[1:5, ])
  for (fam_nm in names(fams)) {
    fits_sum <- suppressWarnings(fit_aug(fam_nm, ctrs))
    fits_dflt <- suppressWarnings(fit_aug(fam_nm, NULL))
    nms_sum <- if (fam_nm == "cumul") {
      names(coef(fits_sum[[1]]))
    } else {
      colnames(coef(fits_sum[[1]]))
    }
    expect_true(all(c("f1", "f2") %in% nms_sum), info = fam_nm)
    expect_false(any(c("fb", "fc") %in% nms_sum), info = fam_nm)
    expect_no_warning(prd_sum <- subprd_augdat(fits_sum, newdata = nd))
    expect_equal(prd_sum, subprd_augdat(fits_dflt, newdata = nd),
                 tolerance = 1e-3, info = fam_nm)
  }
})
