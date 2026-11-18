context("score-screened forward search (argument `search_screen`)")

# A binomial or Poisson reference model with simulated draws of the
# coefficients (no model fit is necessary): the draws are the maximum
# likelihood estimate plus N(0, 0.1^2) noise, so they are not a posterior.
# Predictors: 8 numeric columns and a factor `f` with 3 levels; `fml` selects
# the terms. `trials > 1` gives a binomial response with that many trials (the
# response is the number of successes, the observation weights are the trials).
# `offset = TRUE` adds an offset column `offs`. The fold fits of K-fold CV
# reuse the draws (no refit).
refm_screen <- function(family = binomial(),
                        fml = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + f,
                        n = 60, ndraws = 100, seed = 4127, trials = 1,
                        offset = FALSE, div_minimizer = NULL) {
  set.seed(seed)
  dat <- as.data.frame(matrix(rnorm(n * 8), n, 8,
                              dimnames = list(NULL, paste0("x", 1:8))))
  dat$f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  dat$offs <- if (offset) runif(n, -0.5, 0.5) else 0
  dat$trials <- trials
  eta <- 0.9 * dat$x1 - 0.7 * dat$x2 + 0.5 * dat$x3 +
    c(a = 0, b = 0.6, c = -0.6)[as.character(dat$f)]
  if (family$family == "binomial") {
    # The number of successes (projpred needs it, not the proportion):
    dat$y <- rbinom(n, trials, family$linkinv(eta + dat$offs))
    fml_glm <- update(fml, cbind(y, trials - y) ~ .)
  } else {
    dat$y <- rpois(n, family$linkinv(0.3 * eta + dat$offs))
    fml_glm <- fml
  }
  mm <- model.matrix(fml, dat)
  beta_hat <- coef(suppressWarnings(glm(fml_glm, family = family, data = dat,
                                        offset = offs)))
  beta_hat <- beta_hat[colnames(mm)]
  beta_hat[is.na(beta_hat)] <- 0
  beta <- beta_hat + matrix(rnorm(ncol(mm) * ndraws, sd = 0.1), ncol(mm))
  # The linear predictor without the offset:
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(delete.response(terms(fml)), nd) %*% fit$beta
  }
  extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                 extract_y = TRUE) {
    y_wobs_offs(newdata = newdata, wrhs = wrhs %||% ~ trials,
                orhs = orhs %||% if (offset) ~ offs else NULL,
                resp_form = if (extract_y) ~ y else NULL)
  }
  init_screen <- function(fit, called_from_cvrefbuilder = FALSE) {
    args <- list(fit, data = fit$data, formula = fml, family = family,
                 ref_predfun = ref_predfun,
                 extract_model_data = extract_model_data,
                 cvfun = function(folds) {
                   lapply(seq_len(max(folds)), function(k) {
                     list(beta = fit$beta, data = fit$data[folds != k, ])
                   })
                 },
                 cvrefbuilder = function(cvfit) {
                   init_screen(cvfit, called_from_cvrefbuilder = TRUE)
                 },
                 called_from_cvrefbuilder = called_from_cvrefbuilder)
    if (!is.null(div_minimizer)) args$div_minimizer <- div_minimizer
    do.call(init_refmodel, args)
  }
  init_screen(list(beta = beta, data = dat))
}

# Evaluates `expr` and counts the calls of `fun` (a projpred function) with a
# result (`ok`) and with `NULL` (`null`).
with_call_count <- function(expr, fun = "search_forward_screen") {
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

# Muffles the expected warnings of these small reference models (PSIS and
# subsampled LOO); other warnings stay visible.
muffle_known <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(paste0("Pareto k|PSIS|pareto-k|psis|subsampl|",
                     "non-integer #successes"), conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

vs_screen <- function(refm, search_screen, nterms_max = 5, ...) {
  muffle_known(varsel(refm, search_screen = search_screen, nclusters = 3,
                      ndraws_pred = 20, nterms_max = nterms_max, seed = 5521,
                      verbose = FALSE, ...))
}

# Expects that the screened search ran once and gave the plain forward result.
expect_screen_equal_forward <- function(refm, search_screen, info,
                                        nterms_max = 5) {
  vs_fwd <- vs_screen(refm, NULL, nterms_max = nterms_max)
  res_scr <- with_call_count(vs_screen(refm, search_screen,
                                       nterms_max = nterms_max))
  expect_identical(res_scr$ok, 1L, info = info)
  expect_identical(res_scr$null, 0L, info = info)
  expect_identical(ranking(res_scr$value)$fulldata, ranking(vs_fwd)$fulldata,
                   info = info)
  expect_equal(res_scr$value$ce, vs_fwd$ce, info = info)
  return(invisible(res_scr$value))
}

# Expects that the screened search returned `NULL` (plain forward search).
expect_screen_fallback <- function(refm, info, nterms_max = 4) {
  vs_fwd <- vs_screen(refm, NULL, nterms_max = nterms_max)
  res_scr <- with_call_count(vs_screen(refm, 2, nterms_max = nterms_max))
  expect_identical(res_scr$ok, 0L, info = info)
  expect_identical(res_scr$null, 1L, info = info)
  expect_identical(ranking(res_scr$value)$fulldata, ranking(vs_fwd)$fulldata,
                   info = info)
  expect_equal(res_scr$value$ce, vs_fwd$ce, info = info)
}

test_that("`search_screen` must be `NULL` or a single positive integer", {
  refm <- refm_screen()
  bad <- list(0, -1, 1.5, NA, NA_real_, Inf, "3", c(1, 2), TRUE)
  for (i in seq_along(bad)) {
    expect_error(vs_screen(refm, bad[[i]]),
                 "`search_screen` must be `NULL` or a single positive",
                 info = paste("bad value", i))
  }
})

test_that("`search_screen` with L1 search gives a message", {
  refm <- refm_screen()
  expect_message(vs_screen(refm, 3, method = "L1"),
                 "only takes effect if `method = \"forward\"`")
})

test_that(paste(
  "`search_screen` at least the number of candidates gives the plain forward",
  "search result through the screened search"
), {
  cases <- list(bernoulli = refm_screen(binomial()),
                binomial_trials = refm_screen(binomial(), trials = 4),
                poisson = refm_screen(poisson()))
  for (nm in names(cases)) {
    vs <- expect_screen_equal_forward(cases[[nm]], 20, info = nm)
    expect_identical(vs$args_search$search_screen, 20L, info = nm)
  }
  expect_true(all(cases$binomial_trials$wobs == 4))
  expect_null(vs_screen(cases$poisson, NULL)$args_search$search_screen)
})

test_that(paste(
  "`search_screen = 1` projects one candidate per step, and the approximate",
  "gain picks the forward choices for strong effects"
), {
  count_proj <- function(expr) {
    calls <- new.env()
    calls$n <- 0L
    ns <- asNamespace("projpred")
    suppressMessages(trace(
      "proj_to_submodl", where = ns, print = FALSE,
      tracer = bquote(assign("n", get("n", envir = .(calls)) + 1L,
                             envir = .(calls)))
    ))
    on.exit(suppressMessages(untrace("proj_to_submodl", where = ns)),
            add = TRUE)
    value <- expr
    return(list(value = value, n = calls$n))
  }
  for (fam in list(binomial(), poisson())) {
    refm <- refm_screen(fam, n = 200)
    res_scr <- count_proj(vs_screen(refm, 1, nterms_max = 4))
    # The evaluation of the final path alone:
    res_ev <- count_proj(vs_screen(
      refm, 1, nterms_max = 4,
      search_out = list(search_path = res_scr$value$search_path)
    ))
    # Size 0 (1 projection), then 1 projection for each of the sizes 1-4 (the
    # screening of size 1 starts from the intercept-only submodel):
    expect_identical(res_scr$n - res_ev$n, 1L + 4L, info = fam$family)
    # The first 3 forward choices are the strong terms; gains in a wrong order
    # would pick weak terms (x4-x8) instead:
    rk_fwd <- ranking(vs_screen(refm, NULL, nterms_max = 4))$fulldata
    rk_scr <- ranking(res_scr$value)$fulldata
    expect_identical(rk_scr[1:3], rk_fwd[1:3], info = fam$family)
    expect_false(any(rk_scr[1:3] %in% paste0("x", 4:8)), info = fam$family)
  }
})

test_that("`search_screen` applies to models with offsets", {
  for (fam in list(binomial(), poisson())) {
    refm <- refm_screen(fam, n = 200, offset = TRUE)
    expect_true(any(refm$offset != 0), info = fam$family)
    expect_screen_equal_forward(refm, 20, info = fam$family)
    # The submodels are fitted without the offsets, but `ce` includes them,
    # so the gains and `ce` can order the candidates differently. With 3
    # candidates per step, the first 3 terms are the forward choices here:
    rk_fwd <- ranking(vs_screen(refm, NULL, nterms_max = 4))$fulldata
    rk_scr <- ranking(vs_screen(refm, 3, nterms_max = 4))$fulldata
    expect_identical(rk_scr[1:3], rk_fwd[1:3], info = fam$family)
  }
})

test_that("`search_screen` gives verbose output", {
  refm <- refm_screen()
  expect_message(
    muffle_known(varsel(refm, search_screen = 2, nclusters = 3,
                        ndraws_pred = 20, nterms_max = 4, seed = 5521,
                        verbose = TRUE)),
    "100% of terms selected"
  )
  expect_message(
    muffle_known(varsel(refm_screen(binomial(link = "probit")),
                        search_screen = 2, nclusters = 3, ndraws_pred = 20,
                        nterms_max = 2, seed = 5521, verbose = TRUE)),
    "`search_screen` cannot be used for this model or data"
  )
})

test_that("`search_screen` applies to numeric interaction terms", {
  refm <- refm_screen(fml = y ~ x1 + x2 + x3 + x4 + x1:x2 + x3:x4)
  expect_screen_equal_forward(refm, 10, info = "numeric interactions",
                              nterms_max = 6)
})

test_that(paste(
  "`search_screen` falls back to the plain forward search for a",
  "non-canonical link, a factor in an interaction and a custom divergence",
  "minimizer"
), {
  expect_screen_fallback(refm_screen(binomial(link = "probit")),
                         info = "probit link")
  expect_screen_fallback(refm_screen(fml = y ~ x1 + x2 + f + x1:f),
                         info = "factor interaction")
  expect_screen_fallback(
    refm_screen(div_minimizer = function(...) divmin(...)),
    info = "custom divergence minimizer"
  )
})

test_that("`search_screen` does not change the search of Gaussian models", {
  set.seed(902)
  dat <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
  dat$y <- dat$x1 + rnorm(50)
  fml <- y ~ x1 + x2 + x3
  mm <- model.matrix(fml, dat)
  fit <- list(beta = qr.coef(qr(mm), dat$y) +
                matrix(rnorm(4 * 50, sd = 0.05), 4), data = dat)
  refm_gauss <- suppressWarnings(init_refmodel(
    fit, data = dat, formula = fml, family = gaussian(),
    ref_predfun = function(fit, newdata = NULL) {
      nd <- if (is.null(newdata)) fit$data else newdata
      model.matrix(delete.response(terms(fml)), nd) %*% fit$beta
    },
    dis = rep(1, 50)
  ))
  vs_fwd <- vs_screen(refm_gauss, NULL, nterms_max = 3)
  res_gauss <- with_call_count(
    with_call_count(vs_screen(refm_gauss, 1, nterms_max = 3)),
    fun = "search_forward_gauss"
  )
  # The Gaussian search ran and returned a result, so the screened search was
  # not called:
  expect_identical(res_gauss$ok, 1L)
  res_scr <- res_gauss$value
  expect_identical(res_scr$ok + res_scr$null, 0L)
  expect_identical(ranking(res_scr$value)$fulldata, ranking(vs_fwd)$fulldata)
  expect_equal(res_scr$value$ce, vs_fwd$ce)
})

test_that(paste(
  "`search_screen` reaches the fold searches of cv_varsel() (LOO and K-fold)",
  "and is reused by cv_varsel.vsel() and varsel.vsel()"
), {
  refm <- refm_screen()
  args_cv <- list(refm, validate_search = TRUE, search_screen = 3,
                  nclusters = 3, ndraws_pred = 20, nterms_max = 3,
                  seed = 5521, verbose = FALSE)
  nloo <- 4L
  res_loo <- with_call_count(muffle_known(do.call(
    cv_varsel, c(args_cv, list(cv_method = "LOO", nloo = nloo))
  )))
  # The full-data search and one search per fold:
  expect_identical(res_loo$ok, nloo + 1L)
  expect_identical(res_loo$null, 0L)
  expect_identical(res_loo$value$args_search$search_screen, 3L)

  K <- 3L
  res_kfold <- with_call_count(muffle_known(do.call(
    cv_varsel, c(args_cv, list(cv_method = "kfold", K = K))
  )))
  expect_identical(res_kfold$ok, K + 1L)
  expect_identical(res_kfold$null, 0L)

  vs <- vs_screen(refm, 3, nterms_max = 3)
  res_reuse <- with_call_count(muffle_known(cv_varsel(
    vs, cv_method = "LOO", nloo = nloo, validate_search = TRUE,
    ndraws_pred = 20, seed = 5521, verbose = FALSE
  )))
  # The full-data search is reused; the fold searches are screened:
  expect_identical(res_reuse$ok, nloo)
  expect_identical(res_reuse$value$args_search$search_screen, 3L)

  expect_error(varsel(vs, search_screen = 2),
               "`search_screen`.*cannot be specified")
  expect_message(
    muffle_known(cv_varsel(vs, search_screen = 2, cv_method = "LOO",
                           validate_search = FALSE, ndraws_pred = 20,
                           seed = 5521, verbose = FALSE)),
    "Argument `search_screen` ignored"
  )
  vs2 <- muffle_known(varsel(vs, ndraws_pred = 20, seed = 5521,
                             verbose = FALSE))
  expect_identical(vs2$args_search$search_screen, 3L)
  expect_identical(ranking(vs2)$fulldata, ranking(vs)$fulldata)
})

test_that(paste(
  "the screened search runs with the reference model of the parallel fold",
  "loops (without its fit) and gives the same ranking"
), {
  refm <- refm_screen()
  refm_slim <- refmodel_for_folds(refm)
  expect_null(refm_slim$fit)
  p_ref <- get_refdist(refm, nclusters = 3)
  search_terms <- union("1", split_formula(refm$formula,
                                           data = refm$fetch_data()))
  scr <- function(refmodel) {
    muffle_known(search_forward_screen(
      p_ref = p_ref, refmodel = refmodel, nterms_max = 4, verbose = FALSE,
      search_terms = search_terms, screen_m = 3L, search_control = NULL
    ))
  }
  out_full <- scr(refm)
  out_slim <- scr(refm_slim)
  expect_false(is.null(out_full))
  expect_false(is.null(out_slim))
  expect_identical(out_slim$predictor_ranking, out_full$predictor_ranking)
  expect_length(out_slim$predictor_ranking, 4)
})

test_that(paste(
  "the screened fold searches of LOO CV use the perturbed observation weights",
  "and draws"
), {
  refm <- refm_screen()
  # Records the results of perturb_refmodel() (one per fold) and, for each
  # screened search, its reference model, `p_ref` and the screen_gains() inputs:
  rec <- new.env()
  rec$prt <- list()
  rec$searches <- list()
  perturb_orig <- perturb_refmodel
  screen_orig <- search_forward_screen
  gains_orig <- screen_gains
  testthat::local_mocked_bindings(
    perturb_refmodel = function(...) {
      out <- perturb_orig(...)
      rec$prt <- c(rec$prt, list(out))
      out
    },
    search_forward_screen = function(p_ref, refmodel, ...) {
      cur <- new.env()
      cur$refmodel <- refmodel
      cur$p_ref <- p_ref
      cur$gains <- list()
      rec$searches <- c(rec$searches, list(cur))
      cur$out <- screen_orig(p_ref = p_ref, refmodel = refmodel, ...)
      cur$out
    },
    screen_gains = function(mu_ref, wobs, ...) {
      cur <- rec$searches[[length(rec$searches)]]
      cur$gains <- c(cur$gains, list(list(mu_ref = mu_ref, wobs = wobs)))
      gains_orig(mu_ref = mu_ref, wobs = wobs, ...)
    },
    .package = "projpred"
  )
  nloo <- 4L
  muffle_known(cv_varsel(
    refm, cv_method = "LOO", nloo = nloo, validate_search = TRUE,
    search_screen = 1, nclusters = 3, ndraws_pred = 20, nterms_max = 3,
    seed = 5521, parallel = FALSE, verbose = FALSE
  ))
  # The weight perturbation ran in each fold:
  expect_length(rec$prt, nloo)
  expect_false(any(vapply(rec$prt, is.null, logical(1))))
  # The full-data search and one search per fold, in this order, each with a
  # result of the screened search:
  expect_length(rec$searches, nloo + 1L)
  expect_false(any(vapply(rec$searches, function(cur) is.null(cur$out),
                          logical(1))))
  expect_equal(rec$searches[[1]]$refmodel$wobs, refm$wobs)
  for (k in seq_len(nloo)) {
    cur <- rec$searches[[k + 1L]]
    prt <- rec$prt[[k]]
    info <- paste("fold", k)
    # The search gets the perturbed reference model: weights with a zero for
    # the left-out observation, and shifted draws:
    expect_identical(cur$refmodel$wobs, prt$wobs, info = info)
    expect_identical(cur$refmodel$mu, prt$mu, info = info)
    expect_identical(sum(prt$wobs == 0), 1L, info = info)
    expect_false(isTRUE(all.equal(prt$wobs, refm$wobs)), info = info)
    expect_false(isTRUE(all.equal(prt$mu, refm$mu)), info = info)
    # The score statistic runs and uses these weights and the reference means
    # of the search:
    expect_gt(length(cur$gains), 0)
    for (g in cur$gains) {
      expect_identical(g$wobs, prt$wobs, info = info)
      expect_identical(g$mu_ref, unclass(cur$p_ref$mu), info = info)
    }
  }
})

test_that(paste(
  "screen_gains() equals the one-step reduction computed with solve()",
  "(binomial with trials, weighted draws, single and multi-column candidates)"
), {
  set.seed(731)
  n <- 40
  n_draws <- 3
  x <- matrix(rnorm(n * 6), n, 6)
  # Column 6 is in the span of the intercept and columns 1-2:
  x[, 6] <- 0.5 + x[, 1] - 2 * x[, 2]
  x_c <- scale(x, scale = FALSE)
  wobs <- sample(1:4, n, replace = TRUE)
  wdraws <- c(0.2, 0.5, 0.3)
  mu_ref <- matrix(runif(n * n_draws, 0.1, 0.9), n, n_draws)
  mu_a <- matrix(runif(n * n_draws, 0.2, 0.8), n, n_draws)
  var_fun <- binomial()$variance
  cols_chosen <- c(1L, 2L)
  cols_cands <- list(3L, c(4L, 5L), 6L)
  gains <- screen_gains(x_c = x_c, mu_ref = mu_ref, mu_a = mu_a, wobs = wobs,
                        wdraws = wdraws, var_fun = var_fun,
                        cols_cands = cols_cands, cols_chosen = cols_chosen)
  # The candidate in the span (column 6) has a singular information:
  gains_solve <- vapply(cols_cands[1:2], function(cj) {
    sum(vapply(seq_len(n_draws), function(s) {
      v <- diag(wobs * var_fun(mu_a[, s]))
      x_a <- cbind(1, x_c[, cols_chosen])
      x_j <- x_c[, cj, drop = FALSE]
      g <- -crossprod(x_j, wobs * (mu_ref[, s] - mu_a[, s]))
      info <- t(x_j) %*% v %*% x_j - t(x_j) %*% v %*% x_a %*%
        solve(t(x_a) %*% v %*% x_a, t(x_a) %*% v %*% x_j)
      wdraws[s] * drop(t(g) %*% solve(info, g)) / 2
    }, numeric(1))) / sum(wobs)
  }, numeric(1))
  expect_equal(gains[1:2], gains_solve, tolerance = 1e-10)
  expect_identical(gains[3], 0)
  expect_true(all(gains[1:2] > 0))
})

test_that(paste(
  "`search_screen` at least the number of candidates gives the plain forward",
  "result with user-defined `search_terms` and with a collinear column"
), {
  refm <- refm_screen()
  search_terms <- force_search_terms(
    forced_terms = "x1", optional_terms = c("x2", "x3", "x4", "f")
  )
  vs_fwd <- vs_screen(refm, NULL, search_terms = search_terms, nterms_max = 4)
  res_scr <- with_call_count(vs_screen(refm, 20, search_terms = search_terms,
                                       nterms_max = 4))
  expect_identical(res_scr$ok, 1L)
  expect_identical(ranking(res_scr$value)$fulldata, ranking(vs_fwd)$fulldata)
  expect_equal(res_scr$value$ce, vs_fwd$ce)

  # `I(2 * x1)` is collinear with `x1`:
  refm_col <- refm_screen(fml = y ~ x1 + x2 + x3 + I(2 * x1))
  expect_screen_equal_forward(refm_col, 10, info = "collinear column",
                              nterms_max = 4)
  res_col1 <- with_call_count(vs_screen(refm_col, 1, nterms_max = 4))
  expect_identical(res_col1$ok, 1L)
  expect_length(ranking(res_col1$value)$fulldata, 4)
})
