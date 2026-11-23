# Projection onto negative binomial (NB) submodels (alternation of the fit of
# the coefficients and the projection of the shape, see nb_alternate()). This
# file builds its own reference models and does not need setup.R.

# An NB reference model (log link) with simulated draws of the coefficients
# and the shape (no model fit is necessary).
refm_nb <- function(ndraws = 4, n = 60, seed = 3611, link = "log",
                    offset = NULL, wobs = NULL, family = NULL) {
  set.seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  beta <- c(0.8, 0.5, -0.4, 0.3) +
    matrix(rnorm(4 * ndraws, sd = 0.15), 4, ndraws)
  shape <- 3 * exp(rnorm(ndraws, sd = 0.3))
  xmat <- cbind(1, as.matrix(dat))
  linv <- MASS::negative.binomial(theta = 3, link = link)$linkinv
  if (link != "log") {
    # A positive linear predictor (for the sqrt link, a negative one gives
    # more than one local optimum of the projection):
    beta <- beta + c(3, rep(0, 3))
  }
  dat$y <- rnbinom(n, size = 3, mu = linv(drop(xmat %*% rowMeans(beta))))
  fit <- list(beta = beta, data = dat)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    cbind(1, as.matrix(nd[c("x1", "x2", "x3")])) %*% fit$beta
  }
  # K-fold CV: the "refits" keep all draws (the draws do not depend on the
  # data):
  cvfun <- function(folds) lapply(seq_len(max(folds)), function(k) fit)
  args <- list(fit, data = dat, formula = y ~ x1 + x2 + x3,
               family = family %||%
                 MASS::negative.binomial(theta = 3, link = link),
               ref_predfun = ref_predfun, dis = shape, cvfun = cvfun)
  if (!is.null(offset) || !is.null(wobs)) {
    dat$off <- offset %||% 0
    dat$wobs <- wobs %||% 1
    args$data <- dat
    args$extract_model_data <- function(object, newdata = NULL,
                                        wrhs = NULL, orhs = NULL,
                                        extract_y = TRUE) {
      nd <- if (is.null(newdata)) dat else newdata
      list(y = if (extract_y) nd$y, weights = nd$wobs, offset = nd$off)
    }
  }
  args_fold <- args
  args$cvrefbuilder <- function(cvfit) {
    suppressWarnings(do.call(init_refmodel, c(
      args_fold, list(called_from_cvrefbuilder = TRUE)
    )))
  }
  refm <- suppressWarnings(do.call(init_refmodel, args))
  list(refm = refm, xmat = xmat, mu = linv(xmat %*% beta), shape = shape)
}

nb_proj_ref_link <- function(x, mu, shape, link) {
  nb_proj_ref(x, mu, shape, 1, link = link)
}

# The projection by an independent method: glm.fit() with the shape fixed and
# optimize() on the expected log density (long sums over the support of the
# mixture of the draws with weights `v`).
nb_proj_ref <- function(x, mu_draws, shape_draws, v, link = "log", w = 1) {
  w <- rep_len(w, nrow(mu_draws))
  n <- nrow(mu_draws)
  pm <- lapply(seq_len(n), function(i) {
    k <- max(qnbinom(1e-15, size = shape_draws, mu = mu_draws[i, ],
                     lower.tail = FALSE))
    y <- 0:k
    list(y = y, p = Reduce(`+`, lapply(seq_along(v), function(s) {
      v[s] * dnbinom(y, size = shape_draws[s], mu = mu_draws[i, s])
    })))
  })
  mu_bar <- drop(mu_draws %*% v)
  phi <- sum(v * shape_draws)
  beta <- NULL
  for (it in 1:200) {
    start <- beta
    if (is.null(start) && link != "log") {
      # The start 0 gives the mean 0 for these links:
      linkfun <- MASS::negative.binomial(theta = 1, link = link)$linkfun
      start <- c(linkfun(mean(mu_bar)), rep(0, ncol(x) - 1))
    }
    fit <- suppressWarnings(glm.fit(
      x, mu_bar, family = MASS::negative.binomial(theta = phi, link = link),
      weights = w, start = start,
      control = glm.control(epsilon = 1e-14, maxit = 100)
    ))
    beta <- fit$coefficients
    mu <- fit$fitted.values
    obj <- function(u) {
      sum(w * vapply(seq_len(n), function(i) {
        sum(pm[[i]]$p * dnbinom(pm[[i]]$y, size = exp(u), mu = mu[i],
                                log = TRUE))
      }, numeric(1)))
    }
    phi_new <- exp(optimize(obj, log(phi) + c(-3, 3), maximum = TRUE,
                            tol = 1e-12)$maximum)
    done <- abs(log(phi_new / phi)) < 1e-10
    phi <- phi_new
    if (done) break
  }
  list(beta = unname(beta), phi = phi)
}

# The projections use `thresh_conv = 1e-15`: with the default 1e-7, glm_ridge()
# stops with coefficients that differ by about 1e-4 from the optimum.

test_that("project() gives the NB projection of each draw", {
  ref <- refm_nb(ndraws = 4)
  prj <- suppressMessages(project(ref$refm, predictor_terms = c("x1", "x2"),
                                  ndraws = 4, regul = 0, thresh_conv = 1e-15))
  x <- ref$xmat[, 1:3]
  for (s in 1:4) {
    exp_s <- nb_proj_ref(x, ref$mu[, s, drop = FALSE], ref$shape[s], 1)
    sub <- prj$outdmin[[s]]
    expect_equal(c(sub$alpha, sub$beta), exp_s$beta, tolerance = 1e-6,
                 info = paste("draw", s))
    expect_equal(prj$dis[s], exp_s$phi, tolerance = 1e-5,
                 info = paste("draw", s))
  }
})

test_that("project() gives the NB projection of a cluster (mixture)", {
  ref <- refm_nb(ndraws = 6)
  prj <- project(ref$refm, predictor_terms = c("x1", "x2"), nclusters = 1,
                 regul = 0, thresh_conv = 1e-15)
  exp_c <- nb_proj_ref(ref$xmat[, 1:3], ref$mu, ref$shape, rep(1 / 6, 6))
  sub <- prj$outdmin[[1]]
  expect_equal(c(sub$alpha, sub$beta), exp_c$beta, tolerance = 1e-6)
  expect_equal(prj$dis, exp_c$phi, tolerance = 1e-5)
  # The exact build (the default is the quadrature) gives the same projection
  # within the error of the quadrature:
  opts_old <- options(projpred.nb_gsum = "exact")
  on.exit(options(opts_old), add = TRUE)
  prj_q <- project(ref$refm, predictor_terms = c("x1", "x2"), nclusters = 1,
                   regul = 0, thresh_conv = 1e-15)
  expect_equal(prj_q$dis, prj$dis, tolerance = 1e-6)
})

test_that("varsel() runs a forward search for an NB reference model", {
  ref <- refm_nb(ndraws = 20)
  vs <- suppressMessages(varsel(ref$refm, method = "forward", nclusters = 3,
                                nclusters_pred = 4, nterms_max = 3,
                                verbose = 0))
  expect_setequal(vs$predictor_ranking, c("x1", "x2", "x3"))
  ce <- vs$ce
  expect_true(all(is.finite(ce)))
  expect_true(all(diff(ce) <= 1e-8))
  smmry <- summary(vs, stats = c("elpd", "mlpd"))
  expect_true(all(is.finite(smmry$perf_sub$elpd)))
})

test_that("NB submodels with group-level terms stop with a message", {
  expect_error(
    divmin(y ~ x + (1 | g), data = data.frame(y = 1, x = 1, g = "a"),
           family = extend_family(MASS::negative.binomial(theta = 2)),
           weights = 1, projpred_var = matrix(NA), projpred_ws_aug = NULL,
           verbose_divmin = FALSE, projpred_nb_phi = 2),
    "only for submodels without group-level"
  )
})

test_that("project() gives the NB projections of two clusters", {
  ref <- refm_nb(ndraws = 12)
  prj <- project(ref$refm, predictor_terms = c("x1", "x2"), nclusters = 2,
                 regul = 0, thresh_conv = 1e-15)
  x <- ref$xmat[, 1:3]
  for (j in 1:2) {
    ind <- which(prj$cl_ref == j)
    exp_j <- nb_proj_ref(x, ref$mu[, ind, drop = FALSE], ref$shape[ind],
                         rep(1 / length(ind), length(ind)))
    sub <- prj$outdmin[[j]]
    expect_equal(c(sub$alpha, sub$beta), exp_j$beta, tolerance = 1e-6,
                 info = paste("cluster", j))
    expect_equal(prj$dis[j], exp_j$phi, tolerance = 1e-5,
                 info = paste("cluster", j))
  }
})

test_that("the NB projection converges with the default fitter settings", {
  ref <- refm_nb(ndraws = 6)
  p_ref <- get_refdist(ref$refm, nclusters = 2)
  rec <- new.env()
  rec$iters <- NULL
  trace("nb_alternate", exit = bquote(assign(
    "iters", c(get("iters", envir = .(rec)), attr(returnValue(), "nb_iter")),
    envir = .(rec)
  )), where = asNamespace("projpred"), print = FALSE)
  on.exit(untrace("nb_alternate", where = asNamespace("projpred")),
          add = TRUE)
  expect_no_warning(
    prj <- project(ref$refm, predictor_terms = c("x1", "x2"), nclusters = 2)
  )
  expect_length(rec$iters, 1)
  expect_lte(rec$iters, 10)
  # The result agrees with the precise fit:
  prj_p <- project(ref$refm, predictor_terms = c("x1", "x2"), nclusters = 2,
                   thresh_conv = 1e-15)
  expect_equal(prj$dis, prj_p$dis, tolerance = 1e-3)
})

test_that("the parallel path of divmin() gives the NB projection", {
  ref <- refm_nb(ndraws = 4)
  prj_seq <- suppressMessages(project(ref$refm, predictor_terms = c("x1", "x2"),
                                      ndraws = 4, regul = 0,
                                      thresh_conv = 1e-15))
  opts_old <- options(projpred.parallel_proj_trigger = 1)
  on.exit(options(opts_old), add = TRUE)
  n_par <- 0
  prj_par <- with_mocked_bindings(
    suppressMessages(project(ref$refm, predictor_terms = c("x1", "x2"),
                             ndraws = 4, regul = 0, thresh_conv = 1e-15)),
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
  expect_equal(prj_par$dis, prj_seq$dis, tolerance = 1e-12)
  for (s in 1:4) {
    expect_equal(prj_par$outdmin[[s]]$beta, prj_seq$outdmin[[s]]$beta,
                 tolerance = 1e-12)
  }
})

test_that("the NB projection works with the identity link", {
  ref <- refm_nb(ndraws = 4, link = "identity")
  prj <- suppressMessages(project(ref$refm, predictor_terms = c("x1", "x2"),
                                  ndraws = 4, regul = 0,
                                  thresh_conv = 1e-15))
  x <- ref$xmat[, 1:3]
  for (s in 1:4) {
    mu_s <- ref$mu[, s, drop = FALSE]
    exp_s <- nb_proj_ref_link(x, mu_s, ref$shape[s], "identity")
    sub <- prj$outdmin[[s]]
    expect_equal(c(sub$alpha, sub$beta), exp_s$beta, tolerance = 1e-5,
                 info = paste("draw", s))
    expect_equal(prj$dis[s], exp_s$phi, tolerance = 1e-5,
                 info = paste("draw", s))
  }
})

test_that("unsupported NB cases stop with a message", {
  expect_error(refm_nb(ndraws = 4, offset = rep(c(0, 0.5), 30)),
               "does not support offsets")
  dat <- data.frame(y = rnbinom(20, size = 2, mu = 3), x1 = rnorm(20))
  expect_error(
    suppressWarnings(init_refmodel(NULL, data = dat, formula = y ~ x1,
                                   family = MASS::negative.binomial(2))),
    "`datafit` with the negative binomial"
  )
})

test_that("cv_varsel() with LOO and validated search runs for NB", {
  ref <- refm_nb(ndraws = 20)
  cvvs <- suppressWarnings(suppressMessages(cv_varsel(
    ref$refm, method = "forward", cv_method = "LOO", nloo = 4,
    validate_search = TRUE, nclusters = 2, nclusters_pred = 3,
    nterms_max = 2, verbose = 0
  )))
  smmry <- summary(cvvs, stats = "elpd")
  expect_true(all(is.finite(smmry$perf_sub$elpd)))
})

test_that("varsel() with L1 search and cv_varsel() with K-fold CV run for NB", {
  ref <- refm_nb(ndraws = 20)
  # (setup.R turns this warning off.)
  opts_old <- options(projpred.warn_instable_projections = TRUE)
  on.exit(options(opts_old), add = TRUE)
  expect_warning(
    vs <- suppressMessages(varsel(ref$refm, method = "L1", nclusters_pred = 4,
                                  verbose = 0)),
    "non-multilevel negative binomial models, an L1 search"
  )
  expect_setequal(vs$predictor_ranking, c("x1", "x2", "x3"))
  expect_true(all(is.finite(vs$ce)))
  cvvs <- suppressWarnings(suppressMessages(cv_varsel(
    ref$refm, method = "forward", cv_method = "kfold", K = 2, nclusters = 2,
    nclusters_pred = 3, nterms_max = 2, verbose = 0
  )))
  smmry <- summary(cvvs, stats = "elpd")
  expect_true(all(is.finite(smmry$perf_sub$elpd)))
})

test_that("init_refmodel() names the NB shape draws in its `dis` message", {
  dat <- data.frame(y = rnbinom(20, size = 2, mu = 3), x1 = rnorm(20))
  fit <- list(beta = matrix(c(1, 0.1), 2, 4))
  fam_brms <- structure(list(family = "negbinomial", link = "log"),
                        class = "family")
  expect_error(
    suppressWarnings(init_refmodel(
      fit, data = dat, formula = y ~ x1, family = fam_brms,
      ref_predfun = function(fit, newdata = NULL) {
        cbind(1, (if (is.null(newdata)) dat else newdata)$x1) %*% fit$beta
      }
    )),
    "draws of the shape parameter"
  )
})

test_that("get_refmodel() takes the NB shape of an rstanarm fit", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  set.seed(5)
  dat <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  dat$y <- rnbinom(80, size = 2, mu = exp(1 + 0.6 * dat$x1 - 0.3 * dat$x2))
  fit <- suppressWarnings(rstanarm::stan_glm(
    y ~ x1 + x2 + x3, family = rstanarm::neg_binomial_2(), data = dat,
    chains = 1, iter = 400, seed = 1140, refresh = 0
  ))
  refm <- get_refmodel(fit)
  expect_identical(refm$family$family, "negbinomial")
  expect_equal(refm$dis, unname(as.matrix(fit)[, "reciprocal_dispersion"]))
  vs <- suppressMessages(varsel(refm, method = "forward", nclusters = 2,
                                nclusters_pred = 4, verbose = 0))
  expect_identical(vs$predictor_ranking[1], "x1")
  expect_true(all(is.finite(vs$ce)))
})

test_that("the NB projection works with the sqrt link and with weights", {
  ref <- refm_nb(ndraws = 3, link = "sqrt")
  prj <- suppressMessages(project(ref$refm, predictor_terms = c("x1", "x2"),
                                  ndraws = 3, regul = 0,
                                  thresh_conv = 1e-15))
  x <- ref$xmat[, 1:3]
  for (s in 1:3) {
    exp_s <- nb_proj_ref(x, ref$mu[, s, drop = FALSE], ref$shape[s], 1,
                         link = "sqrt")
    expect_equal(c(prj$outdmin[[s]]$alpha, prj$outdmin[[s]]$beta),
                 exp_s$beta, tolerance = 1e-5, info = paste("sqrt, draw", s))
    expect_equal(prj$dis[s], exp_s$phi, tolerance = 1e-5,
                 info = paste("sqrt, draw", s))
  }
  set.seed(77)
  wobs <- runif(60, 0.2, 3)
  ref <- refm_nb(ndraws = 3, wobs = wobs)
  prj <- suppressMessages(project(ref$refm, predictor_terms = c("x1", "x2"),
                                  ndraws = 3, regul = 0,
                                  thresh_conv = 1e-15))
  for (s in 1:3) {
    exp_s <- nb_proj_ref(x, ref$mu[, s, drop = FALSE], ref$shape[s], 1,
                         w = wobs)
    expect_equal(c(prj$outdmin[[s]]$alpha, prj$outdmin[[s]]$beta),
                 exp_s$beta, tolerance = 1e-5,
                 info = paste("weights, draw", s))
    expect_equal(prj$dis[s], exp_s$phi, tolerance = 1e-5,
                 info = paste("weights, draw", s))
  }
})

test_that("a bare family object with the brms name works", {
  fam_brms <- structure(list(family = "negbinomial", link = "log"),
                        class = "family")
  ref <- refm_nb(ndraws = 3, family = fam_brms)
  ref_mass <- refm_nb(ndraws = 3)
  prj <- suppressMessages(project(ref$refm, predictor_terms = c("x1", "x2"),
                                  ndraws = 3))
  prj_mass <- suppressMessages(project(ref_mass$refm,
                                       predictor_terms = c("x1", "x2"),
                                       ndraws = 3))
  expect_equal(prj$dis, prj_mass$dis, tolerance = 1e-12)
  expect_equal(prj$outdmin[[2]]$beta, prj_mass$outdmin[[2]]$beta,
               tolerance = 1e-12)
})

test_that("a brms reference model with the shape draws as `dis` works", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(6)
  dat <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  dat$y <- rnbinom(80, size = 2, mu = exp(1 + 0.6 * dat$x1))
  fit <- suppressWarnings(suppressMessages(brms::brm(
    y ~ x1 + x2, family = brms::negbinomial(), data = dat, chains = 1,
    iter = 400, seed = 1141, refresh = 0
  )))
  shape <- as.data.frame(fit)$shape
  refm <- get_refmodel(fit, dis = shape)
  expect_identical(refm$family$family, "negbinomial")
  expect_equal(refm$dis, shape)
  prj <- project(refm, predictor_terms = "x1", nclusters = 2)
  expect_true(all(is.finite(prj$dis) & prj$dis > 0))
})
