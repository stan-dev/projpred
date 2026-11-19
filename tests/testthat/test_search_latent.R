context("latent-screened forward search (argument `search_screen`)")

# Reference models with simulated draws of the parameters (no model fit is
# necessary; the draws are the simulation values plus noise, so they are not a
# posterior). Predictors: 8 numeric columns (x6 to x8 have an effect), a factor
# `f` with 3 levels (no effect) and a grouping factor `g` with 8 levels. The
# predictors with an effect are the last candidates of a step, so a screen
# with a wrong order of the candidates projects x1 to x3 first.

dat_lat <- function(n, seed) {
  set.seed(seed)
  dat <- as.data.frame(matrix(rnorm(n * 8), n, 8,
                              dimnames = list(NULL, paste0("x", 1:8))))
  dat$f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  dat$g <- factor(paste0("g", rep_len(1:8, n)))
  dat$offs <- 0
  return(dat)
}
beta_lat <- c(0, 0, 0, 0, 0, 1.2, -0.9, 0.6)
x_lat <- function(dat) as.matrix(dat[, paste0("x", 1:8)])

# A cumulative reference model (4 categories) for the augmented-data
# projection. `fml` selects the terms; `offset = TRUE` adds an offset column.
refm_latcum <- function(fml = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
                        link = "logit", n = 80, ndraws = 30, seed = 7306,
                        offset = FALSE, family_nm = "cumulative") {
  dat <- dat_lat(n, seed)
  if (offset) dat$offs <- runif(n, -0.5, 0.5)
  thr <- c(-1.2, 0, 1.2)
  cats <- paste0("c", 1:4)
  y_lat <- drop(x_lat(dat) %*% beta_lat) + dat$offs + rlogis(n)
  dat$y <- factor(cats[findInterval(y_lat, thr) + 1L], levels = cats)
  beta <- beta_lat + matrix(rnorm(8 * ndraws, sd = 0.1), 8)
  thr_drws <- apply(thr + matrix(rnorm(3 * ndraws, sd = 0.1), 3), 2, sort)
  # Draws x observations x thresholds (without the offset):
  ref_predfun <- function(fit, newdata = NULL) {
    eta <- x_lat(if (is.null(newdata)) fit$data else newdata) %*% fit$beta
    aperm(vapply(1:3, function(k) sweep(-eta, 2, fit$thr[k, ], "+"), eta),
          c(2, 1, 3))
  }
  extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                 extract_y = TRUE) {
    y_wobs_offs(newdata = newdata, wrhs = wrhs,
                orhs = orhs %||% if (offset) ~ offs else NULL,
                resp_form = if (extract_y) ~ y else NULL)
  }
  suppressWarnings(init_refmodel(
    list(beta = beta, thr = thr_drws, data = dat), data = dat, formula = fml,
    family = structure(list(family = family_nm, link = link, cats = cats),
                       class = "family"),
    ref_predfun = ref_predfun, extract_model_data = extract_model_data,
    augdat_y_unqs = cats, augdat_link = augdat_link_cumul,
    augdat_ilink = augdat_ilink_cumul, augdat_args_link = list(link = link),
    augdat_args_ilink = list(link = link)
  ))
}

# A Bernoulli reference model with group-level intercepts, for the latent
# projection (`latent = TRUE`) or the traditional projection.
refm_latml <- function(fml = y ~ x1 + x2 + x3 + x6 + x7 + x8 + (1 | g),
                       latent = TRUE, n = 96, ndraws = 30, seed = 2954) {
  dat <- dat_lat(n, seed)
  b_g <- rnorm(8, sd = 0.8)
  dat$y <- rbinom(n, 1, plogis(drop(x_lat(dat) %*% beta_lat) +
                                 b_g[as.integer(dat$g)]))
  beta <- beta_lat + matrix(rnorm(8 * ndraws, sd = 0.1), 8)
  b_drws <- b_g + matrix(rnorm(8 * ndraws, sd = 0.1), 8)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    x_lat(nd) %*% fit$beta + fit$b[as.integer(nd$g), , drop = FALSE]
  }
  args <- list(list(beta = beta, b = b_drws, data = dat), data = dat,
               formula = fml, family = binomial(), ref_predfun = ref_predfun)
  if (latent) {
    args <- c(args, list(latent = TRUE, dis = rep(1, ndraws)))
  }
  suppressWarnings(do.call(init_refmodel, args))
}

# Evaluates `expr` and counts the calls of `fun` (a projpred function) with a
# result (`ok`) and with `NULL` (`null`).
count_calls_lat <- function(expr, fun) {
  calls <- new.env()
  calls$ok <- 0L
  calls$null <- 0L
  ns <- asNamespace("projpred")
  suppressMessages(trace(
    fun, where = ns, print = FALSE,
    exit = bquote({
      nm <- if (is.null(returnValue(NULL))) "null" else "ok"
      assign(nm, get(nm, envir = .(calls)) + 1L, envir = .(calls))
    })
  ))
  on.exit(suppressMessages(untrace(fun, where = ns)), add = TRUE)
  value <- expr
  return(list(value = value, ok = calls$ok, null = calls$null))
}

# Muffles the expected warnings of these small reference models (PSIS in the
# LOO CV; fractional responses of the traditional projection of a Bernoulli
# model); other warnings stay visible.
muffle_lat <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Pareto k|PSIS|pareto-k|psis|non-integer #successes",
              conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

# varsel() with the argument `search_screen`.
vs_lat <- function(refm, search_screen, nterms_max = 4, verbose = FALSE, ...) {
  muffle_lat(varsel(refm, search_screen = search_screen, nclusters = 2,
                    ndraws_pred = 10, nterms_max = nterms_max, seed = 5521,
                    verbose = verbose, ...))
}

# Expects that the screen `fun` ran once and gave the plain forward result.
expect_lat_equal_forward <- function(refm, search_screen, fun, info = NULL,
                                     nterms_max = 4) {
  vs_fwd <- vs_lat(refm, NULL, nterms_max = nterms_max)
  res <- count_calls_lat(vs_lat(refm, search_screen, nterms_max = nterms_max),
                         fun = fun)
  expect_identical(res$ok, 1L, info = info)
  expect_identical(res$null, 0L, info = info)
  expect_identical(ranking(res$value)$fulldata, ranking(vs_fwd)$fulldata,
                   info = info)
  expect_equal(res$value$ce, vs_fwd$ce, info = info)
  return(invisible(res$value))
}

# Expects that the screen `fun` returned `NULL` (plain forward search).
expect_lat_fallback <- function(refm, fun, info = NULL, nterms_max = 2) {
  vs_fwd <- vs_lat(refm, NULL, nterms_max = nterms_max)
  res <- count_calls_lat(vs_lat(refm, 2, nterms_max = nterms_max), fun = fun)
  expect_identical(res$ok, 0L, info = info)
  expect_identical(res$null, 1L, info = info)
  expect_identical(ranking(res$value)$fulldata, ranking(vs_fwd)$fulldata,
                   info = info)
  expect_equal(res$value$ce, vs_fwd$ce, info = info)
}

# Cumulative models -------------------------------------------------------

test_that("latent_cum_response() gives the binomial quantities for 2 categories", {
  set.seed(1)
  p1 <- matrix(runif(12, 0.05, 0.95), 6, 2)
  p_arr <- aperm(array(c(p1, 1 - p1), dim = c(6, 2, 2)), c(1, 3, 2))
  wobs <- runif(6, 0.5, 2)
  lat <- latent_cum_response(p_arr, wobs = wobs, link = "logit")
  expect_equal(lat$eta, -qlogis(p1))
  expect_equal(lat$v, wobs * p1 * (1 - p1))
  lat_pr <- latent_cum_response(p_arr, wobs = wobs, link = "probit")
  expect_equal(lat_pr$eta, -qnorm(p1))
  expect_equal(lat_pr$v, wobs * dnorm(qnorm(p1))^2 / (p1 * (1 - p1)))
})

test_that("latent_cum_response() gives the Fisher information for 4 categories", {
  set.seed(4)
  n <- 7
  thr <- c(-1.1, 0.2, 1.5)
  eta <- rnorm(n)
  wobs <- runif(n, 0.5, 2)
  for (link in c("logit", "probit")) {
    p_fun <- if (link == "probit") pnorm else plogis
    prbs <- function(eta) {
      t(apply(cbind(0, p_fun(outer(-eta, thr, "+")), 1), 1, diff))
    }
    p <- prbs(eta)
    lat <- latent_cum_response(array(p, dim = c(n, 4, 1)), wobs = wobs,
                               link = link)
    expect_equal(drop(lat$eta), eta - mean(thr), info = link)
    h <- 1e-5
    dp <- (prbs(eta + h) - prbs(eta - h)) / (2 * h)
    expect_equal(drop(lat$v), wobs * rowSums(dp^2 / p), tolerance = 1e-6,
                 info = link)
  }
})

test_that("latent_cum_gains() gives the reduction of the weighted RSS", {
  set.seed(2)
  n <- 40
  x_c <- scale(matrix(rnorm(n * 6), n, 6), scale = FALSE)
  # Column 6 is in the span of the chosen columns 1 and 2:
  x_c[, 6] <- x_c[, 1] - 2 * x_c[, 2]
  eta_ref <- matrix(rnorm(n * 3), n, 3)
  v_ref <- matrix(runif(n * 3, 0.2, 2), n, 3)
  wdraws <- c(0.5, 0.3, 0.2)
  rss <- function(cols, s) {
    sum(v_ref[, s] * lm.wfit(cbind(1, x_c[, cols, drop = FALSE]), eta_ref[, s],
                             v_ref[, s])$residuals^2)
  }
  gain_exp <- vapply(3:5, function(j) {
    sum(wdraws * vapply(1:3, function(s) rss(1:2, s) - rss(c(1:2, j), s),
                        numeric(1)))
  }, numeric(1))
  gain <- latent_cum_gains(x_c, eta_ref, v_ref, wdraws,
                           cols_cands = as.list(3:6), cols_chosen = 1:2)
  expect_equal(gain[1:3], gain_exp)
  expect_identical(gain[4], 0)
  # Without chosen columns:
  gain0 <- latent_cum_gains(x_c, eta_ref, v_ref, wdraws,
                            cols_cands = list(3L), cols_chosen = integer(0))
  expect_equal(gain0, sum(wdraws * vapply(1:3, function(s) {
    rss(integer(0), s) - rss(3L, s)
  }, numeric(1))))
  # A zero weight matrix has no Cholesky factor:
  expect_null(latent_cum_gains(x_c, eta_ref, 0 * v_ref, wdraws,
                               cols_cands = as.list(3:5), cols_chosen = 1:2))
})

test_that("the cumulative screen gives the plain forward result", {
  skip_on_cran()
  for (link in c("logit", "probit")) {
    refm <- refm_latcum(link = link)
    # All candidates are projected:
    expect_lat_equal_forward(refm, 10, "search_forward_latent_cum",
                             info = paste(link, "m 10"))
    # Three candidates per step are projected (3 predictors have an effect):
    expect_lat_equal_forward(refm, 3, "search_forward_latent_cum",
                             info = paste(link, "m 3"), nterms_max = 3)
  }
})

test_that("the cumulative screen with constant gains gives another path", {
  skip_on_cran()
  refm <- refm_latcum()
  vs_fwd <- vs_lat(refm, NULL, nterms_max = 3)
  expect_identical(ranking(vs_fwd)$fulldata[1], "x6")
  testthat::local_mocked_bindings(
    latent_cum_gains = function(cols_cands, ...) rep(1, length(cols_cands)),
    .package = "projpred"
  )
  vs_const <- vs_lat(refm, 3, nterms_max = 3)
  expect_true(ranking(vs_const)$fulldata[1] %in% c("x1", "x2", "x3"))
})

test_that("the cumulative screen projects `search_screen` candidates per step", {
  skip_on_cran()
  refm <- refm_latcum()
  n_prj <- function(search_screen) {
    count_calls_lat(vs_lat(refm, search_screen, nterms_max = 3,
                           refit_prj = FALSE),
                    fun = "proj_to_submodl")$ok
  }
  # The intercept-only submodel, then m candidates per step:
  expect_identical(n_prj(2), 1L + 2L + 2L + 2L)
  expect_identical(n_prj(NULL), 1L + 8L + 7L + 6L)
})

test_that("the latent screens do not run with `search_screen = NULL`", {
  skip_on_cran()
  refm <- refm_latcum()
  res <- count_calls_lat(vs_lat(refm, NULL, nterms_max = 2),
                         fun = "search_forward_latent")
  expect_identical(res$ok + res$null, 0L)
})

test_that("the cumulative screen returns `NULL` outside its conditions", {
  skip_on_cran()
  fun <- "search_forward_latent_cum"
  expect_lat_fallback(refm_latcum(fml = y ~ x1 + x2 + x3 + f), fun,
                      info = "factor with 3 levels")
  expect_lat_fallback(refm_latcum(fml = y ~ x1 + x2 + x3 + offset(offs),
                                  offset = TRUE), fun, info = "offset")
  expect_lat_fallback(refm_latcum(link = "cloglog"), fun, info = "cloglog")
  expect_lat_fallback(refm_latcum(family_nm = "cumulative_rstanarm"), fun,
                      info = "cumulative_rstanarm")
  expect_lat_fallback(refm_latcum(fml = y ~ x1 + x2 + x3 + (1 | g)), fun,
                      info = "group-level term")
})

test_that("the cumulative screen prints its message and the progress", {
  skip_on_cran()
  refm <- refm_latcum()
  msgs <- character()
  withCallingHandlers(
    vs_lat(refm, 2, nterms_max = 2, verbose = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("latent approximation of the cumulative model", msgs)))
  expect_true(any(grepl("100% of terms selected", msgs, fixed = TRUE)))
  expect_false(any(grepl("cannot be used for this model", msgs)))
})

test_that("cv_varsel() with validated search runs the cumulative screen", {
  skip_on_cran()
  refm <- refm_latcum(n = 40)
  cvvs_lat <- function(search_screen) {
    muffle_lat(cv_varsel(refm, cv_method = "LOO", validate_search = TRUE,
                         search_screen = search_screen, nclusters = 2,
                         ndraws_pred = 10, nterms_max = 2, seed = 5521,
                         verbose = FALSE))
  }
  cv_fwd <- cvvs_lat(NULL)
  res <- count_calls_lat(cvvs_lat(10), fun = "search_forward_latent_cum")
  expect_identical(res$null, 0L)
  # The full-data search and one search per fold:
  expect_identical(res$ok, 41L)
  expect_identical(ranking(res$value), ranking(cv_fwd))
  expect_equal(performances(res$value, deltas = FALSE),
               performances(cv_fwd, deltas = FALSE))
})

# Records the searches of cv_varsel() with LOO CV and validated search for the
# screen `screen_fun` (name of a projpred function) in the returned
# environment: the results of perturb_refmodel() (`prt`, one per fold), the
# get_p_clust() calls with a shift of the draws (`pclust`: the result, and the
# result without the shift) and, for each search, an environment with its
# reference model, `p_ref` and result (`searches`, the full-data search
# first). The recorders of a test add to the environment of the current search
# (cur_search()).
local_fold_recorder <- function(screen_fun, .env = parent.frame()) {
  rec <- new.env()
  rec$prt <- list()
  rec$pclust <- list()
  rec$searches <- list()
  ns <- asNamespace("projpred")
  perturb_orig <- get("perturb_refmodel", ns)
  pclust_orig <- get("get_p_clust", ns)
  screen_orig <- get(screen_fun, ns)
  mocks <- list(
    perturb_refmodel = function(...) {
      out <- perturb_orig(...)
      rec$prt <- c(rec$prt, list(out))
      out
    },
    get_p_clust = function(..., draws_shift = NULL) {
      out <- pclust_orig(..., draws_shift = draws_shift)
      if (!is.null(draws_shift)) {
        rec$pclust <- c(rec$pclust, list(list(
          out = out, unshifted = pclust_orig(..., draws_shift = NULL)
        )))
      }
      out
    }
  )
  mocks[[screen_fun]] <- function(p_ref, refmodel, ...) {
    cur <- new.env()
    cur$refmodel <- refmodel
    cur$p_ref <- p_ref
    rec$searches <- c(rec$searches, list(cur))
    cur$out <- screen_orig(p_ref = p_ref, refmodel = refmodel, ...)
    cur$out
  }
  do.call(testthat::local_mocked_bindings,
          c(mocks, list(.package = "projpred", .env = .env)))
  return(rec)
}

cur_search <- function(rec) rec$searches[[length(rec$searches)]]

# cv_varsel() with LOO CV, validated search in `nloo` folds (sequential) and
# `search_screen = 2`.
cv_lat_folds <- function(refm, nloo) {
  muffle_lat(cv_varsel(refm, cv_method = "LOO", nloo = nloo,
                       validate_search = TRUE, search_screen = 2,
                       nclusters = 2, ndraws_pred = 10, nterms_max = 2,
                       seed = 5521, parallel = FALSE, verbose = FALSE))
}

# Expects that `rec` (see local_fold_recorder()) has the screened full-data
# search with the reference model `refm` and one screened search per fold with
# the perturbed reference model of the fold. `shift` is where this reference
# model has the shift of the draws: in `mu` (shifted draws) or in
# `draws_shift` (get_p_clust() adds it to `p_ref`).
expect_fold_searches_perturbed <- function(rec, refm, nloo, shift) {
  expect_length(rec$prt, nloo)
  expect_false(any(vapply(rec$prt, is.null, logical(1))))
  expect_length(rec$searches, nloo + 1L)
  expect_false(any(vapply(rec$searches, function(cur) is.null(cur$out),
                          logical(1))))
  expect_equal(rec$searches[[1]]$refmodel$wobs, refm$wobs)
  expect_null(rec$searches[[1]]$refmodel$draws_shift)
  expect_length(rec$pclust, if (shift == "draws_shift") nloo else 0L)
  for (k in seq_len(nloo)) {
    cur <- rec$searches[[k + 1L]]
    prt <- rec$prt[[k]]
    info <- paste("fold", k)
    # Weights with a zero for the left-out observation:
    expect_identical(cur$refmodel$wobs, prt$wobs, info = info)
    expect_identical(sum(prt$wobs == 0), 1L, info = info)
    expect_false(isTRUE(all.equal(prt$wobs, refm$wobs)), info = info)
    if (shift == "mu") {
      expect_identical(cur$refmodel$mu, prt$mu, info = info)
      expect_false(isTRUE(all.equal(prt$mu, refm$mu)), info = info)
    } else {
      expect_identical(cur$refmodel$draws_shift, prt$draws_shift,
                       info = info)
      expect_true(any(prt$draws_shift != 0), info = info)
      expect_identical(cur$p_ref, rec$pclust[[k]]$out, info = info)
    }
  }
}

test_that(paste(
  "the cumulative screen in the fold searches of LOO CV uses the perturbed",
  "observation weights and draws"
), {
  skip_on_cran()
  refm <- refm_latcum(n = 40)
  rec <- local_fold_recorder("search_forward_latent_cum")
  resp_orig <- latent_cum_response
  gains_orig <- latent_cum_gains
  testthat::local_mocked_bindings(
    latent_cum_response = function(p_arr, wobs, link) {
      cur <- cur_search(rec)
      cur$p_arr <- p_arr
      cur$wobs <- wobs
      resp_orig(p_arr = p_arr, wobs = wobs, link = link)
    },
    latent_cum_gains = function(...) {
      cur <- cur_search(rec)
      cur$n_gains <- (cur$n_gains %||% 0L) + 1L
      gains_orig(...)
    },
    .package = "projpred"
  )
  nloo <- 4L
  cv_lat_folds(refm, nloo)
  expect_fold_searches_perturbed(rec, refm, nloo, shift = "mu")
  for (k in seq_len(nloo)) {
    cur <- rec$searches[[k + 1L]]
    info <- paste("fold", k)
    # The screen uses the perturbed weights and the probabilities of the
    # search's `p_ref` (observations x 4 categories x projected draws), and
    # its criterion runs:
    expect_identical(cur$wobs, rec$prt[[k]]$wobs, info = info)
    expect_identical(
      cur$p_arr,
      array(unclass(cur$p_ref$mu),
            dim = c(refm$nobs, 4L, length(cur$p_ref$wdraws_prj))),
      info = info
    )
    expect_gt(cur$n_gains %||% 0L, 0)
  }
})

# Latent projection with group-level terms --------------------------------

test_that("latent_mlvl_formula() combines the group-level terms of a factor", {
  fml_base <- .projpred_eta ~ 1
  chr <- function(terms_) {
    paste(deparse(latent_mlvl_formula(terms_, fml_base)[[3]]), collapse = "")
  }
  expect_identical(chr(c("1", "x1")), "x1")
  expect_identical(chr(c("1", "x1", "(1 | g)")), "x1 + (1 | g)")
  expect_identical(chr(c("1", "x1", "(1 | g)", "x2 + (x2 | g)")),
                   "x1 + x2 + (x2 | g)")
  expect_identical(chr(c("1", "x1", "(1 | g)", "(x1 | g)", "x2 + (x2 | g)")),
                   "x1 + x2 + (x1 + x2 | g)")
  expect_identical(chr(c("1", "(1 | g)", "(x1 + x2 | g)", "(1 | h)")),
                   "(x1 + x2 | g) + (1 | h)")
})

test_that("latent_mlvl_rss() gives the weighted RSS of the screen fit", {
  dat <- dat_lat(64, seed = 3)
  dat$.projpred_eta <- dat$x1 - 0.5 * dat$x2 +
    rnorm(8)[as.integer(dat$g)] + rnorm(64, sd = 0.3)
  dat$.projpred_w <- runif(64, 0.5, 2)
  fml_base <- .projpred_eta ~ 1
  # Without group-level terms: weighted least squares.
  fit_lm <- lm(.projpred_eta ~ x1 + f, data = dat, weights = .projpred_w)
  expect_equal(latent_mlvl_rss(c("1", "x1", "f"), dat, fml_base),
               sum(dat$.projpred_w * residuals(fit_lm)^2))
  # With group-level terms: lme4::lmer() with maximum likelihood.
  fit_lmer <- lme4::lmer(.projpred_eta ~ x1 + (x2 | g), data = dat,
                         weights = .projpred_w, REML = FALSE)
  expect_equal(
    latent_mlvl_rss(c("1", "x1", "(1 | g)", "(x2 | g)"), dat, fml_base),
    sum(dat$.projpred_w * (dat$.projpred_eta - fitted(fit_lmer))^2),
    tolerance = 1e-5
  )
  # A failing lme4::lmer() fit gives `NA` (here: as many levels as
  # observations):
  dat$id <- factor(seq_len(nrow(dat)))
  expect_identical(latent_mlvl_rss(c("1", "x1", "(1 | id)"), dat, fml_base),
                   NA_real_)
  # A variable that is not in the data is an error, not a failed fit:
  expect_error(latent_mlvl_rss(c("1", "x_none"), dat, fml_base), "x_none")
  expect_error(latent_mlvl_rss(c("1", "x_none", "(1 | g)"), dat, fml_base),
               "x_none")
  expect_error(latent_mlvl_rss(c("1", "x1", "(1 | g_none)"), dat, fml_base),
               "g_none")
})

test_that("the multilevel screen gives the plain forward result", {
  skip_on_cran()
  refm <- refm_latml()
  expect_lat_equal_forward(refm, 10, "search_forward_latent_mlvl",
                           info = "m 10")
  expect_lat_equal_forward(refm, 3, "search_forward_latent_mlvl",
                           info = "m 3", nterms_max = 3)
})

test_that("the multilevel screen handles group-level slopes", {
  skip_on_cran()
  refm <- refm_latml(fml = y ~ x1 + x2 + x3 + x6 + x7 + (x6 | g))
  expect_lat_equal_forward(refm, 3, "search_forward_latent_mlvl",
                           nterms_max = 3)
})

test_that("the multilevel screen projects the candidates with a failed fit", {
  skip_on_cran()
  refm <- refm_latml()
  testthat::local_mocked_bindings(
    latent_mlvl_rss = function(...) NA_real_, .package = "projpred"
  )
  n_prj <- count_calls_lat(vs_lat(refm, 2, nterms_max = 3, refit_prj = FALSE),
                           fun = "proj_to_submodl")$ok
  # The intercept-only submodel, then all candidates of each step (7 terms):
  expect_identical(n_prj, 1L + 7L + 6L + 5L)
})

test_that("the multilevel screen projects one failed candidate and the m best", {
  skip_on_cran()
  refm <- refm_latml()
  rss_orig <- latent_mlvl_rss
  testthat::local_mocked_bindings(
    latent_mlvl_rss = function(terms_, ...) {
      if (identical(utils::tail(terms_, 1), "x1")) {
        return(NA_real_)
      }
      rss_orig(terms_, ...)
    },
    .package = "projpred"
  )
  res <- count_calls_lat(vs_lat(refm, 2, nterms_max = 2, refit_prj = FALSE),
                         fun = "proj_to_submodl")
  # The intercept-only submodel, then the 2 best and x1 in each step:
  expect_identical(res$ok, 1L + 3L + 3L)
  expect_false("x1" %in% ranking(res$value)$fulldata)
})

test_that("the multilevel screen with a constant criterion gives another path", {
  skip_on_cran()
  refm <- refm_latml()
  vs_fwd <- vs_lat(refm, NULL, nterms_max = 3)
  testthat::local_mocked_bindings(
    latent_mlvl_rss = function(...) 1, .package = "projpred"
  )
  vs_const <- vs_lat(refm, 3, nterms_max = 3)
  expect_false(identical(ranking(vs_const)$fulldata, ranking(vs_fwd)$fulldata))
})

test_that("the multilevel screen drops a factor level of a zero-weight observation", {
  skip_on_cran()
  # Level "d" of `f` has only observation 1, which gets weight 0 (as in a fold
  # of the LOO CV):
  refm <- refm_latml(fml = y ~ x1 + x2 + x6 + x7 + f + (1 | g))
  refm$wobs[1] <- 0
  fetch_orig <- refm$fetch_data
  refm$fetch_data <- function(...) {
    d <- fetch_orig(...)
    d$f <- factor(as.character(d$f), levels = c("a", "b", "c", "d"))
    d$f[1] <- "d"
    d
  }
  set.seed(5521)
  p_ref <- get_refdist(refm, nclusters = 2)
  rss_na <- 0L
  rss_orig <- latent_mlvl_rss
  testthat::local_mocked_bindings(
    latent_mlvl_rss = function(terms_, dat, ...) {
      expect_false("d" %in% levels(dat$f))
      out <- rss_orig(terms_, dat, ...)
      rss_na <<- rss_na + is.na(out)
      out
    },
    .package = "projpred"
  )
  out <- muffle_lat(search_forward_latent_mlvl(
    p_ref = p_ref, refmodel = refm, nterms_max = 3, verbose = FALSE,
    search_terms = union("1", split_formula(refm$formula,
                                              data = refm$fetch_data())),
    screen_m = 2L
  ))
  expect_false(is.null(out))
  expect_length(out$predictor_ranking, 3)
  expect_identical(rss_na, 0L)
})

test_that("the multilevel screen uses the contrasts of the reference model", {
  skip_on_cran()
  # The screen fit uses the `contrasts` attribute: with a reduced contrast
  # matrix, the RSS of `x1 + f` is that of the regression on (1, x1, ctr[f, ]),
  # and it is larger than with the default contrasts.
  dat <- dat_lat(64, seed = 3)
  dat$.projpred_eta <- dat$x1 + c(a = 0, b = 1, c = -1)[as.character(dat$f)] +
    rnorm(64, sd = 0.3)
  dat$.projpred_w <- runif(64, 0.5, 2)
  fml_base <- .projpred_eta ~ 1
  ctr <- contr.sum(levels(dat$f))[, 1, drop = FALSE]
  rss_dflt <- latent_mlvl_rss(c("1", "x1", "f"), dat, fml_base)
  attr(dat$f, "contrasts") <- ctr
  rss_ctr <- latent_mlvl_rss(c("1", "x1", "f"), dat, fml_base)
  fit_ctr <- lm.wfit(cbind(1, dat$x1, ctr[as.integer(dat$f), ]),
                     dat$.projpred_eta, dat$.projpred_w)
  expect_equal(rss_ctr, sum(dat$.projpred_w * fit_ctr$residuals^2))
  expect_gt(rss_ctr - rss_dflt, 1e-6)

  # The screen passes the contrasts of the reference model to its fits. The
  # weight 0 of observation 1 (as in a fold of the LOO CV) checks the path with
  # removed rows; the other observations have all levels of `f`:
  refm <- refm_latml(fml = y ~ x1 + x2 + x6 + x7 + f + (1 | g))
  expect_true(all(table(refm$fetch_data()$f[-1]) > 0))
  ctr <- contr.sum(levels(refm$fetch_data()$f))[, 1, drop = FALSE]
  refm$contrasts <- list(f = ctr)
  refm$wobs[1] <- 0
  set.seed(5521)
  p_ref <- get_refdist(refm, nclusters = 2)
  ctrs_seen <- list()
  rss_orig <- latent_mlvl_rss
  testthat::local_mocked_bindings(
    latent_mlvl_rss = function(terms_, dat, ...) {
      ctrs_seen <<- c(ctrs_seen, list(attr(dat$f, "contrasts")))
      rss_orig(terms_, dat, ...)
    },
    .package = "projpred"
  )
  out <- muffle_lat(search_forward_latent_mlvl(
    p_ref = p_ref, refmodel = refm, nterms_max = 2, verbose = FALSE,
    search_terms = union("1", split_formula(refm$formula,
                                              data = refm$fetch_data())),
    screen_m = 2L
  ))
  expect_false(is.null(out))
  expect_gt(length(ctrs_seen), 0)
  for (ctr_seen in ctrs_seen) {
    expect_identical(ctr_seen, ctr)
  }
})

test_that("the multilevel screen returns `NULL` outside its conditions", {
  skip_on_cran()
  fun <- "search_forward_latent_mlvl"
  expect_lat_fallback(refm_latml(latent = FALSE), fun,
                      info = "traditional projection")
})

test_that("the multilevel screen prints its message and the progress", {
  skip_on_cran()
  refm <- refm_latml()
  msgs <- character()
  withCallingHandlers(
    vs_lat(refm, 2, nterms_max = 2, verbose = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("one multilevel fit to the mean latent predictor",
                        msgs)))
  expect_true(any(grepl("100% of terms selected", msgs, fixed = TRUE)))
})

test_that("cv_varsel() with validated search runs the multilevel screen", {
  skip_on_cran()
  refm <- refm_latml(n = 48)
  cvvs_lat <- function(search_screen) {
    muffle_lat(cv_varsel(refm, cv_method = "LOO", validate_search = TRUE,
                         search_screen = search_screen, nclusters = 2,
                         ndraws_pred = 10, nterms_max = 2, seed = 5521,
                         verbose = FALSE))
  }
  cv_fwd <- cvvs_lat(NULL)
  res <- count_calls_lat(cvvs_lat(10), fun = "search_forward_latent_mlvl")
  expect_identical(res$null, 0L)
  # The full-data search and one search per fold:
  expect_identical(res$ok, 49L)
  expect_identical(ranking(res$value), ranking(cv_fwd))
  expect_equal(performances(res$value, deltas = FALSE),
               performances(cv_fwd, deltas = FALSE))
})

test_that(paste(
  "the multilevel screen in the fold searches of LOO CV uses the perturbed",
  "observation weights and draws"
), {
  skip_on_cran()
  refm <- refm_latml(n = 48)
  rec <- local_fold_recorder("search_forward_latent_mlvl")
  rss_orig <- latent_mlvl_rss
  testthat::local_mocked_bindings(
    latent_mlvl_rss = function(terms_, dat, fml_base) {
      cur <- cur_search(rec)
      cur$n_rss <- (cur$n_rss %||% 0L) + 1L
      if (is.null(cur$dat)) cur$dat <- dat
      rss_orig(terms_, dat, fml_base)
    },
    .package = "projpred"
  )
  nloo <- 4L
  cv_lat_folds(refm, nloo)
  expect_fold_searches_perturbed(rec, refm, nloo, shift = "draws_shift")
  for (k in seq_len(nloo)) {
    cur <- rec$searches[[k + 1L]]
    prt <- rec$prt[[k]]
    info <- paste("fold", k)
    # The screen fits run, with the perturbed weights (without the left-out
    # observation):
    expect_gt(cur$n_rss %||% 0L, 0)
    keep <- prt$wobs > 0
    expect_equal(unname(cur$dat$.projpred_w), unname(prt$wobs[keep]),
                 info = info)
    # Their response is the mean latent predictor of `p_ref` without the
    # shift, plus the shift:
    p_unshifted <- rec$pclust[[k]]$unshifted
    wd <- p_unshifted$wdraws_prj
    eta_unshifted <- drop(unclass(p_unshifted$mu) %*% (wd / sum(wd)))
    expect_equal(unname(cur$dat$.projpred_eta - eta_unshifted[keep]),
                 unname(prt$draws_shift[keep]), tolerance = 1e-8,
                 info = info)
  }
})
