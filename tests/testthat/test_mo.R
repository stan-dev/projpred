context("monotonic terms mo()")

# Helpers -----------------------------------------------------------------

# Imitates brms::mo(), which returns a list that model.frame() cannot use. The
# reference model formulas below have this function in their environment (as
# if brms was attached), so the tests also check that projpred uses its own
# mo().
mo_brms_like <- function(x, id = NA) {
  structure(list(term = deparse(substitute(x))),
            class = c("mo_term", "sp_term"))
}

# Simulated data and draws of a reference model with the brms parametrization
# b * D * sum(zeta[seq_len(j)]) of monotonic effects (`j` is the position of
# the value in the levels, starting at 0): `x` is an integer with levels
# x_min, ..., x_min + D, `f` an ordered factor with 3 levels, `z` and `w` are
# continuous (`w` has no effect), and `g` is a grouping variable.
mo_sim <- function(n = 120, S = 40, D = 4, x_min = 0L, seed = 5602) {
  set.seed(seed)
  dat <- data.frame(
    x = sample(rep_len(x_min + 0:D, n)),
    f = factor(sample(rep_len(c("lo", "mid", "hi"), n)),
               levels = c("lo", "mid", "hi"), ordered = TRUE),
    z = rnorm(n),
    w = rnorm(n),
    g = factor(sample(rep_len(paste0("gr", 1:6), n)))
  )
  zeta_x <- rev(seq_len(D)) / sum(seq_len(D))
  zeta_f <- c(0.7, 0.3)
  b <- cbind(alpha = rnorm(S, 0.2, 0.05), bx = rnorm(S, 0.8, 0.05),
             bf = rnorm(S, -0.6, 0.05), bz = rnorm(S, 0.5, 0.05))
  # Linear predictor (N x S) for data `nd`:
  eta_fun <- function(nd) {
    mo_x <- D * c(0, cumsum(zeta_x))[nd$x - x_min + 1]
    mo_f <- 2 * c(0, cumsum(zeta_f))[as.integer(nd$f)]
    cbind(1, mo_x, mo_f, nd$z) %*% t(b)
  }
  eta_mean <- rowMeans(eta_fun(dat))
  dat$y_gauss <- eta_mean + rnorm(n, sd = 0.5)
  dat$y_binom <- rbinom(n, 1, plogis(eta_mean))
  dat$y_poiss <- rpois(n, exp(eta_mean))
  # The coefficients of the cumulative indicators (S x number of columns):
  coef_true <- cbind(
    "(Intercept)" = b[, "alpha"],
    b[, "bx"] * D * matrix(zeta_x, S, D, byrow = TRUE,
                           dimnames = list(NULL,
                                           paste0("mo(x)", x_min + seq_len(D)))),
    b[, "bf"] * 2 * matrix(zeta_f, S, 2, byrow = TRUE,
                           dimnames = list(NULL, c("mo(f)mid", "mo(f)hi"))),
    z = b[, "bz"]
  )
  return(nlist(dat, eta_fun, coef_true, S, D))
}

# A reference model for `sim` (see mo_sim()). `ref_predfun_arr` turns the N x S
# linear predictor into the output of `ref_predfun` (for the augmented-data
# projection).
mo_refm <- function(sim, fml, family = gaussian(), data = sim$dat,
                    ref_predfun_arr = identity, ...) {
  environment(fml) <- list2env(list(mo = mo_brms_like),
                               parent = environment(fml))
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    ref_predfun_arr(fit$sim$eta_fun(nd))
  }
  dis <- if (family$family == "gaussian") rep(0.5, sim$S) else NULL
  cvfun <- function(folds) {
    lapply(seq_len(max(folds)), function(k) {
      list(sim = sim, data = data[folds != k, , drop = FALSE],
           omitted = which(folds == k))
    })
  }
  dot_args <- list(...)
  cvrefbuilder <- function(cvfit) {
    do.call(init_refmodel, c(
      list(object = cvfit, data = cvfit$data, formula = fml, family = family,
           ref_predfun = ref_predfun, dis = dis, cvfun = cvfun,
           called_from_cvrefbuilder = TRUE),
      dot_args
    ))
  }
  init_refmodel(
    list(sim = sim, data = data), data = data, formula = fml, family = family,
    ref_predfun = ref_predfun, dis = dis, cvfun = cvfun,
    cvrefbuilder = cvrefbuilder, ...
  )
}

# Evaluate `expr` without the warnings whose message contains `msg`.
muffle_warning <- function(expr, msg) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(msg, conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# MASS::polr() warns about non-integer data. projpred filters this warning for
# the draw-wise fits, but the filter misses calls with `:` (here
# `MASS::polr()`), so the warning appears once per draw.
muffle_polr_warning <- function(expr) {
  muffle_warning(expr, "non-integer #successes in a binomial glm!")
}

# The training rows of `nd` (a subset of the training data):
rows_of <- function(nd) as.integer(rownames(nd))

# `a` without attributes other than the dimensions:
strip_attr <- function(a) {
  if (is.null(dim(a))) {
    return(as.vector(a))
  }
  return(array(as.vector(a), dim = dim(a)))
}

# New data that lack some levels of `x` and `f`:
mo_newdata <- function(sim) {
  dat <- sim$dat
  dat[dat$x <= 2 & dat$f != "hi", , drop = FALSE][1:8, ]
}

# mo() --------------------------------------------------------------------

test_that("mo() gives the cumulative indicators of an integer predictor", {
  x <- c(0, 2, 1, 3, 2)
  b <- mo(x)
  expect_identical(dim(b), c(5L, 3L))
  expect_identical(colnames(b), c("1", "2", "3"))
  expect_identical(attr(b, "mo_levels"), 0:3)
  expect_identical(strip_attr(b), outer(x, 1:3, ">=") * 1)
  # With the levels of the training data:
  b_new <- mo(c(1, 0), mo_levels = 0:3)
  expect_identical(strip_attr(b_new), rbind(c(1, 0, 0), c(0, 0, 0)))
  # As in brms, the levels start at the smallest value (also if negative):
  b_shift <- mo(c(3, 1, 2))
  expect_identical(attr(b_shift, "mo_levels"), 1:3)
  expect_identical(strip_attr(b_shift), rbind(c(1, 1), c(0, 0), c(1, 0)))
  b_neg <- mo(c(0, -1, 1))
  expect_identical(colnames(b_neg), c("0", "1"))
  expect_identical(strip_attr(b_neg), rbind(c(1, 0), c(0, 0), c(1, 1)))
  # Missing values give missing rows:
  expect_true(all(is.na(mo(c(1, NA, 2))[2, ])))
})

test_that("mo() gives the cumulative indicators of an ordered factor", {
  f <- factor(c("b", "a", "c", "b"), levels = c("a", "b", "c"), ordered = TRUE)
  b <- mo(f)
  expect_identical(colnames(b), c("b", "c"))
  expect_identical(attr(b, "mo_levels"), c("a", "b", "c"))
  expect_identical(strip_attr(b), cbind(c(1, 0, 1, 1), c(0, 0, 1, 0)))
  # A factor with dropped levels uses the levels of the training data:
  b_new <- mo(droplevels(f[2]), mo_levels = c("a", "b", "c"))
  expect_identical(dim(b_new), c(1L, 2L))
})

test_that("mo() throws errors for unsupported predictors", {
  expect_error(mo(factor(c("a", "b"))), "ordered factor")
  expect_error(mo(c(0, 1.5)), "integer vector")
  expect_error(mo(c("a", "b")), "integer vector")
  expect_error(mo(c(0, 0)), "at least two levels")
  expect_error(mo(c(0, 4), mo_levels = 0:3), "outside the range")
  expect_error(mo(c(-1, 2), mo_levels = 0:3), "outside the range")
  f_new <- factor("d", levels = c("a", "d"), ordered = TRUE)
  expect_error(mo(f_new, mo_levels = c("a", "b", "c")), "do not exist")
  expect_error(mo(c(0, Inf)), "integer vector")
  expect_error(mo(c(0, Inf), mo_levels = 0:3), "integer vector")
})

test_that("mo() matches new data by the type of the training levels", {
  lvls <- c("lo", "mid", "hi")
  b_ref <- strip_attr(mo(factor(c("mid", "hi"), levels = lvls,
                                ordered = TRUE)))
  # Character vectors and unordered factors are matched by their labels:
  expect_identical(strip_attr(mo(c("mid", "hi"), mo_levels = lvls)), b_ref)
  expect_identical(strip_attr(mo(factor(c("mid", "hi")), mo_levels = lvls)),
                   b_ref)
  # Integer codes are not accepted for an ordered factor, and factors are not
  # accepted for an integer:
  expect_error(mo(c(2L, 3L), mo_levels = lvls), "factor or a character")
  expect_error(mo(factor(1), mo_levels = 0:3), "must be an integer in new")
})

test_that("mo() terms in unsupported positions throw an error", {
  sim <- mo_sim(n = 60, S = 4)
  expect_error(mo_refm(sim, y_gauss ~ z + (mo(x) | g)),
               "not supported in group-level terms")
  expect_error(mo_refm(sim, y_gauss ~ mo(x) + s(z) + (1 | g)),
               "not supported in additive multilevel models")
})

test_that(paste(
  "collapse_ranked_predictors() does not change mo() terms by the",
  "replacements for factor terms"
), {
  dat <- data.frame(y = rnorm(6), gbx = c(0:2, 0:2),
                    g = factor(rep(c("a", "b"), 3)))
  fml <- y ~ mo(gbx) + g
  environment(fml) <- list2env(list(mo = mo), parent = globalenv())
  fr <- model.frame(fml, data = dat)
  expect_identical(
    collapse_ranked_predictors(c("mo(gbx)1", "gb", "mo(gbx)2"),
                               formula = fml, data = fr),
    c("mo(gbx)", "g")
  )
})

test_that("the terms of mo() store the training levels in \"predvars\"", {
  dat <- data.frame(x = c(0, 1, 3, 2), y = 1:4)
  mf <- model.frame(y ~ mo(x), data = dat)
  pv <- attr(attr(mf, "terms"), "predvars")
  expect_identical(pv[[3]][[1]], quote(projpred:::mo))
  expect_identical(pv[[3]]$mo_levels, 0:3)
  # The stored call works from the global environment (as for mgcv::gam(),
  # which sets the environment of the terms to the global environment):
  expect_equal(dim(eval(pv[[3]], list(x = 1), globalenv())), c(1L, 3L))
})

test_that("mo_formula_env() makes `mo` refer to projpred's mo()", {
  env_user <- new.env()
  assign("mo", mo_brms_like, envir = env_user)
  assign("user_obj", 1, envir = env_user)
  fml <- y ~ mo(x) + z
  environment(fml) <- env_user
  fml_new <- mo_formula_env(fml)
  expect_identical(get("mo", envir = environment(fml_new)), mo)
  expect_identical(parent.env(environment(fml_new)), env_user)
  expect_identical(get("user_obj", envir = environment(fml_new)), 1)
  # A second call does not add another environment:
  expect_identical(environment(mo_formula_env(fml_new)),
                   environment(fml_new))
  # Formulas without mo() calls (also with a variable called `mo`) stay as
  # they are:
  fml_no <- y ~ mo + z
  environment(fml_no) <- env_user
  expect_identical(environment(mo_formula_env(fml_no)), env_user)
})

# Traditional projection --------------------------------------------------

test_that(paste(
  "the projection onto the full submodel reproduces the monotonic effects of",
  "the reference model draw by draw"
), {
  sim <- mo_sim(x_min = -2L)
  trms <- c("mo(x)", "mo(f)", "z", "w")
  for (fam in list(gaussian(), binomial(), poisson())) {
    y_nm <- switch(fam$family, gaussian = "y_gauss", binomial = "y_binom",
                   poisson = "y_poiss")
    fml <- as.formula(paste(y_nm, "~ mo(x) + mo(f) + z + w"))
    refm <- mo_refm(sim, fml, family = fam)
    prj <- project(refm, predictor_terms = trms, ndraws = sim$S,
                   verbose = FALSE)
    coefs <- do.call(rbind, lapply(prj$outdmin, function(o) {
      drop(coef(o))
    }))
    # The ridge penalty of the submodel fitter (`regul`) has a larger effect
    # for the non-Gaussian families:
    tol <- if (fam$family == "gaussian") 1e-4 else 1e-3
    expect_equal(strip_attr(coefs[, colnames(sim$coef_true)]),
                 strip_attr(sim$coef_true), tolerance = tol,
                 info = fam$family)
    expect_equal(unname(coefs[, "w"]), rep(0, sim$S), tolerance = tol,
                 info = fam$family)
    expect_equal(strip_attr(proj_linpred(prj, transform = FALSE)$pred),
                 t(sim$eta_fun(sim$dat)), tolerance = tol, info = fam$family)
  }
})

test_that(paste(
  "projected predictions for new data with mo() terms equal the predictions",
  "for the same training observations"
), {
  sim <- mo_sim()
  nd <- mo_newdata(sim)
  refm <- mo_refm(sim, y_gauss ~ mo(x) + mo(f) + z + w)
  prj <- project(refm, predictor_terms = c("mo(x)", "mo(f)", "z"),
                 ndraws = 25, seed = 3906, verbose = FALSE)
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = nd)$pred
  expect_equal(unname(pred_new), unname(pred_train[, rows_of(nd)]),
               tolerance = 1e-10)
  # New data with a value of `x` above the training levels:
  nd_out <- nd[1, ]
  nd_out$x <- sim$D + 1L
  expect_error(proj_linpred(prj, newdata = nd_out), "outside the range")
})

test_that("forward search ranks mo() terms as single terms", {
  sim <- mo_sim()
  refm <- mo_refm(sim, y_gauss ~ mo(x) + mo(f) + z + w)
  vs <- varsel(refm, nclusters = 4, nclusters_pred = 8, seed = 3906,
               verbose = FALSE)
  rk <- ranking(vs)$fulldata
  expect_setequal(rk, c("mo(x)", "mo(f)", "z", "w"))
  expect_identical(utils::tail(rk, 1), "w")
})

test_that(paste(
  "L1 search ranks mo() terms as single terms, also with more than 10",
  "levels, and its submodels predict new data"
), {
  sim <- mo_sim(D = 11)
  nd <- sim$dat[sim$dat$x <= 5, , drop = FALSE][1:8, ]
  refm <- mo_refm(sim, y_gauss ~ mo(x) + z + w)
  vs <- varsel(refm, method = "L1", refit_prj = FALSE, seed = 3906,
               verbose = FALSE)
  rk <- ranking(vs)$fulldata
  expect_setequal(rk, c("mo(x)", "z", "w"))
  prj <- project(vs, nterms = 3, refit_prj = FALSE, verbose = FALSE)
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = nd)$pred
  expect_equal(unname(drop(pred_new)), unname(drop(pred_train)[rows_of(nd)]),
               tolerance = 1e-10)
})

test_that("interactions of mo() terms work in forward and L1 search", {
  sim <- mo_sim()
  refm <- mo_refm(sim, y_gauss ~ mo(x) + z + mo(x):z)
  vs_fw <- varsel(refm, nclusters = 4, nclusters_pred = 8, seed = 3906,
                  verbose = FALSE)
  rk_fw <- ranking(vs_fw)$fulldata
  expect_setequal(rk_fw, c("mo(x)", "z", "mo(x):z"))
  expect_identical(rk_fw[3], "mo(x):z")
  vs_l1 <- varsel(refm, method = "L1", seed = 3906, verbose = FALSE)
  expect_setequal(ranking(vs_l1)$fulldata, c("mo(x)", "z", "mo(x):z"))
})

test_that("multilevel and additive submodels with mo() terms work", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("mgcv")
  sim <- mo_sim()
  nd <- mo_newdata(sim)
  # The binomial multilevel fits to non-integer pseudo-data throw lme4
  # convergence warnings; this test only compares the predictions:
  opts_old <- options(projpred.warn_additive_experimental = FALSE,
                      projpred.warn_instable_projections = FALSE,
                      projpred.warn_proj_drawwise = FALSE,
                      projpred.check_convergence = FALSE)
  on.exit(options(opts_old), add = TRUE)
  cases <- list(
    glmm = list(fml = y_gauss ~ mo(x) + z + (1 | g),
                trms = c("mo(x)", "z", "(1 | g)"), family = gaussian()),
    glmm_binom = list(fml = y_binom ~ mo(x) + z + (1 | g),
                      trms = c("mo(x)", "z", "(1 | g)"), family = binomial()),
    gam = list(fml = y_gauss ~ mo(x) + s(z), trms = c("mo(x)", "s(z)"),
               family = gaussian())
  )
  for (nm in names(cases)) {
    refm <- mo_refm(sim, cases[[nm]]$fml, family = cases[[nm]]$family)
    prj <- muffle_warning(muffle_warning(
      project(refm, predictor_terms = cases[[nm]]$trms, ndraws = 25,
              seed = 3906, verbose = FALSE),
      "Model failed to converge"
    ), "Model is nearly unidentifiable")
    pred_train <- proj_linpred(prj)$pred
    pred_new <- proj_linpred(prj, newdata = nd)$pred
    expect_equal(unname(pred_new), unname(pred_train[, rows_of(nd)]),
                 tolerance = 1e-8, info = nm)
  }
})

# Latent and augmented-data projection ------------------------------------

test_that("the latent projection reproduces the monotonic effects", {
  sim <- mo_sim()
  refm <- suppressMessages(mo_refm(sim, y_binom ~ mo(x) + mo(f) + z,
                                   family = binomial(), latent = TRUE))
  prj <- project(refm, predictor_terms = c("mo(x)", "mo(f)", "z"),
                 ndraws = sim$S, verbose = FALSE)
  expect_equal(strip_attr(proj_linpred(prj, transform = FALSE)$pred),
               t(sim$eta_fun(sim$dat)), tolerance = 1e-4)
})

test_that(paste(
  "the augmented-data projection (cumulative family) reproduces the",
  "monotonic effects and predicts new data"
), {
  skip_if_not_installed("MASS")
  # Levels 1, ..., 5 (as in a Likert scale):
  sim <- mo_sim(x_min = 1L)
  thres <- c(-0.5, 0.8)
  dat <- sim$dat
  eta_mean <- rowMeans(sim$eta_fun(dat))
  dat$y_cumul <- factor(
    1L + (runif(nrow(dat)) > plogis(thres[1] - eta_mean)) +
      (runif(nrow(dat)) > plogis(thres[2] - eta_mean)),
    levels = 1:3
  )
  # The output of `ref_predfun`: an S x N x C_lat array of thres_c - eta:
  thres_arr <- function(eta) {
    arr <- vapply(thres, function(th) t(th - eta), t(eta))
    return(arr)
  }
  refm <- mo_refm(
    sim, y_cumul ~ mo(x) + mo(f) + z,
    family = structure(list(family = "cumulative_rstanarm", link = "logit"),
                       class = "family"),
    data = dat, ref_predfun_arr = thres_arr,
    augdat_y_unqs = as.character(1:3),
    augdat_link = augdat_link_cumul, augdat_ilink = augdat_ilink_cumul,
    augdat_args_link = list(link = "logit"),
    augdat_args_ilink = list(link = "logit")
  )
  prj <- muffle_polr_warning(project(
    refm, predictor_terms = c("mo(x)", "mo(f)", "z"), ndraws = sim$S,
    verbose = FALSE
  ))
  pred_train <- strip_attr(proj_linpred(prj, transform = FALSE)$pred)
  expect_equal(pred_train, thres_arr(sim$eta_fun(dat)), tolerance = 1e-3)
  nd <- mo_newdata(sim)
  pred_new <- proj_linpred(prj, newdata = dat[rows_of(nd), ],
                           transform = FALSE)$pred
  expect_equal(strip_attr(pred_new),
               pred_train[, rows_of(nd), , drop = FALSE], tolerance = 1e-10)
})

# Cross-validation --------------------------------------------------------

test_that("cv_varsel() works with mo() terms (LOO and K-fold CV)", {
  # The Pareto k warnings come from the small number of draws:
  opts_old <- options(projpred.warn_psis = FALSE)
  on.exit(options(opts_old), add = TRUE)
  sim <- mo_sim()
  refm <- mo_refm(sim, y_gauss ~ mo(x) + mo(f) + z + w)
  cv_loo <- muffle_warning(cv_varsel(
    refm, cv_method = "LOO", validate_search = TRUE, nloo = 5, nclusters = 3,
    nclusters_pred = 6, seed = 3906, verbose = FALSE
  ), "Some Pareto k diagnostic values are too high")
  expect_setequal(ranking(cv_loo)$fulldata, c("mo(x)", "mo(f)", "z", "w"))
  cv_kf <- cv_varsel(
    refm, cv_method = "kfold", K = 2, nclusters = 3, nclusters_pred = 6,
    seed = 3906, verbose = FALSE
  )
  expect_length(ranking(cv_kf)$foldwise, 2 * 4)
  expect_true(all(is.finite(summary(cv_kf)$perf_sub$elpd)))
})

test_that(paste(
  "K-fold CV throws an error if the training data of a fold lack the largest",
  "value of an integer mo() predictor"
), {
  sim <- mo_sim()
  # Only one observation with the largest value:
  is_max <- sim$dat$x == sim$D
  sim$dat$x[is_max] <- sim$D - 1L
  sim$dat$x[which(is_max)[1]] <- sim$D
  refm <- mo_refm(sim, y_gauss ~ mo(x) + z)
  expect_error(
    cv_varsel(refm, cv_method = "kfold", K = 2, nclusters = 3,
              nclusters_pred = 6, seed = 3906, verbose = FALSE),
    "outside the range of the training data"
  )
})
