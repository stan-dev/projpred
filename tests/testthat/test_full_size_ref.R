# The submodel with all terms of the reference model gets the results of the
# reference model (global option `projpred.full_size_ref`, see
# full_size_ref_ok()).

# A Gaussian reference model with simulated draws of the coefficients (no model
# fit is necessary). With `mlvl = TRUE`, the formula has a group-level
# intercept, and the draws have group-level effects.
refm_full_ref <- function(dat, mlvl = FALSE, ndraws = 100, seed = 4127) {
  fml <- if (mlvl) y ~ x1 + x2 + (1 | g) else y ~ x1 + x2 + x3
  fml_fix <- if (mlvl) y ~ x1 + x2 else fml
  mm <- model.matrix(fml_fix, dat)
  beta_hat <- drop(solve(crossprod(mm), crossprod(mm, dat$y)))
  beta <- beta_hat + withr::with_seed(seed, {
    matrix(rnorm(ncol(mm) * ndraws, sd = 0.05), ncol(mm))
  })
  b <- NULL
  if (mlvl) {
    res <- dat$y - drop(mm %*% beta_hat)
    b_hat <- c(tapply(res, dat$g, mean))
    b <- b_hat + withr::with_seed(seed + 1L, {
      matrix(rnorm(length(b_hat) * ndraws, sd = 0.05), length(b_hat),
             dimnames = list(names(b_hat), NULL))
    })
  }
  # Dispersion draws that differ between the draws:
  dis <- withr::with_seed(seed + 2L, exp(rnorm(ndraws, sd = 0.2)))
  fit <- list(beta = beta, b = b, data = dat, fml_fix = fml_fix)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    eta <- model.matrix(delete.response(terms(fit$fml_fix)), nd) %*% fit$beta
    if (!is.null(fit$b)) {
      eta <- eta + fit$b[as.character(nd$g), , drop = FALSE]
    }
    unname(eta)
  }
  cvfun <- function(folds) {
    lapply(seq_len(max(folds)), function(k) {
      list(dat = dat[folds != k, , drop = FALSE])
    })
  }
  cvrefbuilder <- function(cvfit) {
    refm_full_ref(cvfit$dat, mlvl = mlvl, ndraws = ndraws, seed = seed)
  }
  suppressWarnings(init_refmodel(
    fit, data = dat, formula = fml, family = gaussian(),
    ref_predfun = ref_predfun, dis = dis, cvfun = cvfun,
    cvrefbuilder = cvrefbuilder
  ))
}

dat_full_ref <- function(n = 40, seed = 3517) {
  set.seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                    g = gl(5, n / 5))
  dat$y <- 0.8 * dat$x1 - 0.4 * dat$x2 +
    c(-0.6, -0.3, 0, 0.3, 0.6)[dat$g] + rnorm(n, sd = 0.7)
  return(dat)
}

# Evaluates `expr` and muffles the warnings whose message matches `pattern`.
muffle_warns <- function(expr, pattern) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(pattern, conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

# Evaluates `expr` and counts the calls of proj_to_submodl() for a submodel
# with at least `nterms` terms (without the intercept). Returns the value of
# `expr` and the count.
with_full_prj_count <- function(expr, nterms) {
  calls <- new.env()
  calls$n <- 0L
  ns <- asNamespace("projpred")
  suppressMessages(trace(
    "proj_to_submodl", where = ns, print = FALSE,
    tracer = bquote({
      trms <- setdiff(unique(unlist(predictor_terms)), "1")
      if (count_terms_chosen(c("1", trms)) - 1L >= .(nterms)) {
        assign("n", get("n", envir = .(calls)) + 1L, envir = .(calls))
      }
    })
  ))
  on.exit(suppressMessages(untrace("proj_to_submodl", where = ns)),
          add = TRUE)
  value <- expr
  return(list(value = value, n = calls$n))
}

# Evaluates `expr` and counts the calls of search_forward_gauss() (`search`) and
# perf_eval_gauss() (`eval`) with a result (not `NULL`), and the calls of
# perf_eval_gauss() with a result for the submodel with all search terms
# (`eval_full`). Returns the value of `expr` and the counts.
with_gauss_counts <- function(expr) {
  calls <- new.env()
  calls$search <- 0L
  calls$eval <- 0L
  calls$eval_full <- 0L
  add_one <- function(nm) {
    assign(nm, get(nm, envir = calls) + 1L, envir = calls)
  }
  ns <- asNamespace("projpred")
  suppressMessages(trace(
    "search_forward_gauss", where = ns, print = FALSE,
    exit = bquote(if (!is.null(returnValue(NULL))) .(add_one)("search"))
  ))
  suppressMessages(trace(
    "perf_eval_gauss", where = ns, print = FALSE,
    exit = bquote(if (!is.null(returnValue(NULL))) {
      .(add_one)("eval")
      if (length(search_path$predictor_ranking) %in% nterms) {
        .(add_one)("eval_full")
      }
    })
  ))
  on.exit({
    suppressMessages(untrace("search_forward_gauss", where = ns))
    suppressMessages(untrace("perf_eval_gauss", where = ns))
  }, add = TRUE)
  value <- expr
  return(list(value = value, search = calls$search, eval = calls$eval,
              eval_full = calls$eval_full))
}

# The checks for a `vsel` object `vs` whose search reached the submodel with all
# search terms.
expect_full_ref <- function(vs, info) {
  # The position of the submodel with all search terms in the ranking (with
  # custom `search_terms`, this can be smaller than its size):
  D <- length(vs$predictor_ranking)
  expect_true(ranking_is_full(vs$predictor_ranking, vs$nterms_all),
              info = info)
  expect_true(vs$full_ref, info = info)
  expect_identical(vs$summaries$sub[[D + 1]], vs$summaries$ref, info = info)
  if (!is.null(vs$summaries_fast)) {
    expect_identical(vs$summaries_fast$sub[[D + 1]], vs$summaries_fast$ref,
                     info = info)
  }
  # (Subsampled LOO CV can give a negative variance estimate for the RMSE.)
  smmry <- muffle_warns(
    summary(vs, stats = c("elpd", "rmse"), type = c("mean", "se", "upper"),
            deltas = TRUE),
    "negative variance estimate"
  )
  ps <- smmry$perf_sub
  expect_equal(ps$elpd[ps$size == D], 0, info = info)
  expect_equal(ps$elpd.se[ps$size == D], 0, info = info)
  expect_equal(ps$elpd.upper[ps$size == D], 0, info = info)
  expect_equal(ps$rmse[ps$size == D], 0, info = info)
  expect_identical(smmry$full_ref_size, D, info = info)
}

test_that(paste(
  "the submodel with all terms of the reference model gets the results of the",
  "reference model in varsel() and cv_varsel(), without a projection"
), {
  withr::local_options(projpred.mssg_ndraws = FALSE, projpred.warn_psis = FALSE,
                       projpred.warn_kfold_refits = FALSE,
                       projpred.full_size_ref = NULL)
  dat <- dat_full_ref()
  refm <- refm_full_ref(dat)
  for (gc in c(TRUE, FALSE)) {
    info_gc <- paste("gc:", gc)
    withr::with_options(list(projpred.search_gauss_crossprod = gc), {
      # Runs `expr` and checks that no projection onto the submodel with all
      # terms takes place (also not in the cross-product evaluation). Returns
      # the value of `expr` and the counts of with_gauss_counts():
      run <- function(expr) {
        # (The few draws of the evaluation in LOO CV give Pareto k warnings.)
        res <- muffle_warns(
          with_gauss_counts(with_full_prj_count(expr, nterms = 3L)),
          "Pareto k diagnostic"
        )
        expect_identical(res$value$n, 0L, info = info_gc)
        expect_identical(res$eval_full, 0L, info = info_gc)
        if (!gc) {
          expect_identical(res$search + res$eval, 0L, info = info_gc)
        }
        return(list(vs = res$value$value, search = res$search,
                    eval = res$eval))
      }
      res <- run(varsel(refm, nclusters = 3, ndraws_pred = 50))
      vs <- res$vs
      expect_full_ref(vs, info = paste("varsel,", info_gc))
      # The search did not project onto the submodel with all terms (with
      # `gc`, the search was the cross-product search):
      expect_identical(res$search, as.integer(gc), info = info_gc)
      expect_null(vs$search_path$outdmins[[4]])
      expect_false(is.null(vs$search_path$outdmins[[3]]))
      expect_true(all(diff(vs$ce) <= 1e-12), info = info_gc)
      # The smaller sizes are the same as without the option:
      vs_prj <- withr::with_options(
        list(projpred.full_size_ref = FALSE),
        varsel(refm, nclusters = 3, ndraws_pred = 50)
      )
      expect_false(vs_prj$full_ref)
      expect_identical(vs$summaries$sub[1:3], vs_prj$summaries$sub[1:3],
                       info = info_gc)
      expect_identical(vs$ce[1:3], vs_prj$ce[1:3], info = info_gc)
      # With all draws, the projection onto the submodel with all terms gives
      # (nearly) the results of the reference model:
      vs_all <- varsel(refm, nclusters = 3, ndraws_pred = 100)
      vs_all_prj <- withr::with_options(
        list(projpred.full_size_ref = FALSE),
        varsel(refm, nclusters = 3, ndraws_pred = 100)
      )
      expect_equal(vs_all$ce[4], vs_all_prj$ce[4], tolerance = 1e-9,
                   info = info_gc)
      expect_equal(unname(vs_all$summaries$sub[[4]]$lppd),
                   unname(vs_all_prj$summaries$sub[[4]]$lppd),
                   tolerance = 1e-4, info = info_gc)
      vs <- run(varsel(refm, nclusters = 3, refit_prj = FALSE))$vs
      expect_full_ref(vs, info = paste("refit_prj = FALSE,", info_gc))
      expect_null(vs$search_path$outdmins[[4]])
      vs <- run(cv_varsel(refm, nclusters = 3, ndraws_pred = 50,
                          validate_search = FALSE))$vs
      expect_full_ref(vs, info = paste("LOO, validate_search = FALSE,",
                                       info_gc))
      res <- run(cv_varsel(refm, nclusters = 2, nclusters_pred = 3))
      expect_full_ref(res$vs, info = paste("LOO,", info_gc))
      if (gc) {
        # The fold-wise searches and evaluations used the cross-products:
        expect_gt(res$search, 0L)
        expect_gt(res$eval, 0L)
      }
      vs <- run(cv_varsel(refm, nclusters = 2, nclusters_pred = 3,
                          nloo = 10))$vs
      expect_full_ref(vs, info = paste("subsampled LOO,", info_gc))
      vs <- run(cv_varsel(refm, cv_method = "kfold", K = 3, nclusters = 2,
                          nclusters_pred = 3))$vs
      expect_full_ref(vs, info = paste("K-fold,", info_gc))
      vs <- run(cv_varsel(refm, cv_method = "kfold", K = 3, nclusters = 2,
                          nclusters_pred = 3, validate_search = FALSE))$vs
      expect_full_ref(vs, info = paste("K-fold, validate_search = FALSE,",
                                       info_gc))
    })
  }
  # L1 search:
  res <- with_full_prj_count(
    varsel(refm, method = "L1", ndraws_pred = 50),
    nterms = 3L
  )
  expect_identical(res$n, 0L)
  expect_full_ref(res$value, info = "L1")
  # The parallel code path of the CV folds (without parallel workers, the folds
  # run in this R session):
  for (cv_meth in c("LOO", "kfold")) {
    args_cv <- list(refm, cv_method = cv_meth, nclusters = 2,
                    nclusters_pred = 3, parallel = TRUE)
    args_cv <- c(args_cv, if (cv_meth == "LOO") list(nloo = 5) else list(K = 3))
    res <- muffle_warns(
      with_full_prj_count(do.call(cv_varsel, args_cv), nterms = 3L),
      "Pareto k diagnostic|there are no parallel workers"
    )
    expect_identical(res$n, 0L, info = cv_meth)
    expect_full_ref(res$value, info = paste("parallel", cv_meth))
  }
  # Custom search terms that give all terms of the reference model:
  res <- with_full_prj_count(
    varsel(refm, nclusters = 3, ndraws_pred = 50,
           search_terms = c("x1", "x2 + x3")),
    nterms = 3L
  )
  expect_identical(res$n, 0L)
  expect_full_ref(res$value, info = "search_terms")
})

test_that(paste(
  "the submodel with all search terms is projected if it is not the reference",
  "model or if global option `projpred.full_size_ref` is `FALSE`"
), {
  withr::local_options(projpred.mssg_ndraws = FALSE, projpred.warn_psis = FALSE,
                       projpred.full_size_ref = NULL)
  dat <- dat_full_ref()
  refm <- refm_full_ref(dat)
  check_prj <- function(vs, n, info) {
    D <- length(vs$predictor_ranking)
    expect_false(vs$full_ref, info = info)
    expect_gt(n, 0L)
    expect_false(identical(vs$summaries$sub[[D + 1]], vs$summaries$ref),
                 info = info)
    expect_null(summary(vs)$full_ref_size)
  }
  res <- withr::with_options(
    list(projpred.full_size_ref = FALSE),
    with_full_prj_count(varsel(refm, nclusters = 3, ndraws_pred = 50), 3L)
  )
  check_prj(res$value, res$n, info = "option")
  res <- with_full_prj_count(
    varsel(refm, nclusters = 3, ndraws_pred = 50, search_terms = c("x1", "x2")),
    nterms = 2L
  )
  check_prj(res$value, res$n, info = "search_terms")
  vs <- varsel(refm, nclusters = 3, ndraws_pred = 50, nterms_max = 2)
  expect_false(vs$full_ref)
  # A search without the fit of the submodel with all terms, reused with
  # `refit_prj = FALSE` and the option `FALSE`: the evaluation projects onto
  # that submodel.
  vs <- varsel(refm, nclusters = 3, refit_prj = FALSE)
  expect_null(vs$search_path$outdmins[[4]])
  res <- withr::with_options(
    list(projpred.full_size_ref = FALSE),
    with_full_prj_count(varsel(vs, refit_prj = FALSE), 3L)
  )
  check_prj(res$value, res$n, info = "reused search")
  # project() with `predictor_terms` (all terms) projects, also with
  # `refit_prj = FALSE` (the search has no fit of this submodel, so the
  # evaluation projects with the clusters of the search):
  res <- with_full_prj_count(
    project(vs, predictor_terms = vs$predictor_ranking, ndraws = 20),
    3L
  )
  expect_identical(res$n, 1L)
  expect_false(is_ref_draws(res$value))
  expect_length(res$value$outdmin, 20)
  res <- with_full_prj_count(
    project(vs, predictor_terms = vs$predictor_ranking, refit_prj = FALSE),
    3L
  )
  expect_identical(res$n, 1L)
  expect_false(is_ref_draws(res$value))
  expect_length(res$value$outdmin, 3)
})

test_that(paste(
  "the submodel with all terms of a multilevel reference model is not",
  "projected and gives no warning about its projection"
), {
  skip_if_not_installed("lme4")
  withr::local_options(projpred.mssg_ndraws = FALSE, projpred.warn_psis = FALSE,
                       projpred.warn_instable_projections = TRUE,
                       projpred.full_size_ref = NULL)
  dat <- dat_full_ref()
  refm <- refm_full_ref(dat, mlvl = TRUE)
  res <- with_full_prj_count(
    expect_no_warning(
      vs <- varsel(refm, nclusters = 2, ndraws_pred = 20),
      message = "projection onto the full model"
    ),
    nterms = 3L
  )
  expect_identical(res$n, 0L)
  expect_full_ref(vs, info = "multilevel")
  expect_true("(1 | g)" %in% vs$predictor_ranking)
  # project() at the full size gives the reference model draws, without the
  # warning about the projection onto the full model:
  warns <- capture_warnings(prj <- project(vs, nterms = 3, ndraws = 20))
  expect_match(warns, "returns 20 draws of the reference model", all = FALSE)
  expect_false(any(grepl("projection onto the full model", warns)))
  expect_true(is_ref_draws(prj))
  # suggest_size() works, and with the lower bound, only the full size (delta 0
  # with SE 0) meets the criterion:
  expect_true(suggest_size(vs, warnings = FALSE) %in% 0:3)
  expect_identical(suggest_size(vs, type = "lower", warnings = FALSE), 3)
})

test_that("summary() and plot() mark the size of the reference model", {
  skip_if_not_installed("ggplot2")
  withr::local_options(projpred.mssg_ndraws = FALSE,
                       projpred.full_size_ref = NULL)
  dat <- dat_full_ref()
  refm <- refm_full_ref(dat)
  vs <- varsel(refm, nclusters = 3, ndraws_pred = 50)
  out <- capture.output(print(summary(vs)))
  expect_true(any(grepl("^Size 3 \\(all search terms\\) shows the results of",
                        out)))
  # Not if the table stops before the full size:
  out <- capture.output(print(summary(vs, nterms_max = 2)))
  expect_false(any(grepl("^Size 3", out)))
  shape_open <- function(pp) {
    vapply(pp$layers, function(l) {
      inherits(l$geom, "GeomPoint") && identical(l$aes_params$shape, 1)
    }, logical(1))
  }
  for (dlt in list(FALSE, TRUE, "mixed")) {
    pp <- plot(vs, deltas = dlt)
    is_open <- shape_open(pp)
    expect_identical(sum(is_open), 1L, info = paste("deltas:", dlt))
    expect_identical(pp$layers[[which(is_open)]]$data$size, 3,
                     info = paste("deltas:", dlt))
    expect_match(pp$labels$subtitle, "Open point: reference model",
                 info = paste("deltas:", dlt))
  }
  pp <- plot(vs, nterms_max = 2)
  expect_identical(sum(shape_open(pp)), 0L)
  # A `vsel` object of an older projpred version (without element `full_ref`):
  vs$full_ref <- NULL
  expect_null(summary(vs)$full_ref_size)
  expect_identical(sum(shape_open(plot(vs))), 0L)
})

test_that(paste(
  "project() returns the reference model draws for the submodel with all",
  "terms of the reference model, and the prediction functions use them"
), {
  withr::local_options(projpred.mssg_ndraws = FALSE,
                       projpred.full_size_ref = NULL)
  dat <- dat_full_ref()
  refm <- refm_full_ref(dat)
  vs <- varsel(refm, nclusters = 3, ndraws_pred = 50)
  warn_ref <- "returns 50 draws of the reference model, not a projection"
  res <- with_full_prj_count(
    expect_warning(prj <- project(vs, nterms = 3, ndraws = 50), warn_ref),
    nterms = 3L
  )
  expect_identical(res$n, 0L)
  expect_s3_class(prj, "projection")
  expect_true(is_ref_draws(prj))
  expect_equal(prj$nprjdraws, 50)
  expect_true(prj$const_wdraws_prj)
  draws <- prj$outdmin$draws
  expect_identical(sort(draws), which(!is.na(prj$cl_ref)))
  expect_identical(prj$cl_ref[draws], seq_len(50))
  # The thinned draws of the reference model:
  p_thin <- get_refdist(refm, ndraws = 50)
  expect_equal(prj$dis, p_thin$dis)
  expect_equal(prj$dis, refm$dis[draws])
  pl <- proj_linpred(prj)
  expect_equal(unname(pl$pred), unname(t(p_thin$mu)))
  eta_ref <- refm$ref_predfun(refm$fit)
  expect_equal(unname(pl$lpd), unname(t(dnorm(
    dat$y, eta_ref[, draws], rep(prj$dis, each = nrow(dat)), log = TRUE
  ))))
  # The projection onto the submodel with all terms reproduces these draws
  # (up to the convergence and the regularization of the fits):
  prj_full <- project(vs, predictor_terms = vs$predictor_ranking, ndraws = 50)
  expect_equal(prj$ce, prj_full$ce, tolerance = 1e-6)
  expect_equal(prj$dis, prj_full$dis, tolerance = 1e-4)
  expect_equal(unname(pl$pred), unname(proj_linpred(prj_full)$pred),
               tolerance = 1e-4)
  nd <- dat[1:5, ]
  expect_equal(
    unname(proj_linpred(prj, newdata = nd)$pred),
    unname(t(refm$ref_predfun(refm$fit, newdata = nd)[, draws]))
  )
  expect_identical(dim(proj_predict(prj, .seed = 1)), c(50L, nrow(dat)))
  expect_output(print(prj), "Number of reference model draws")
  expect_error(as.matrix(prj), "draws of the reference model, not projected")
  expect_error(mo_effects(prj), "draws of the reference model, not projected")
  # proj_linpred() of the `vsel` object:
  expect_warning(pl_vs <- proj_linpred(vs, nterms = 3, ndraws = 50), warn_ref)
  expect_equal(pl_vs$pred, pl$pred)
  # Several sizes (one warning):
  warns <- capture_warnings(prjs <- project(vs, nterms = 2:3, ndraws = 50))
  expect_length(warns, 1)
  expect_false(is_ref_draws(prjs[[1]]))
  expect_true(is_ref_draws(prjs[[2]]))
  expect_length(proj_linpred(prjs), 2)
  # `nclusters` and `refit_prj = FALSE`:
  expect_warning(prj7 <- project(vs, nterms = 3, nclusters = 7),
                 "returns 7 draws.*`nclusters` gives their number")
  expect_equal(prj7$nprjdraws, 7)
  expect_true(prj7$const_wdraws_prj)
  expect_warning(prj_f <- project(vs, nterms = 3, refit_prj = FALSE,
                                  ndraws = 50),
                 warn_ref)
  expect_identical(prj_f$outdmin, prj$outdmin)
  # A `vsel` object of cv_varsel():
  cvvs <- muffle_warns(
    cv_varsel(refm, validate_search = FALSE, nclusters = 3, ndraws_pred = 50),
    "Pareto k"
  )
  expect_warning(prj_cv <- project(cvvs, nterms = 3, ndraws = 50), warn_ref)
  expect_identical(prj_cv$outdmin, prj$outdmin)

  # The latent projection (response-scale predictions):
  dat_b <- dat
  dat_b$y <- withr::with_seed(2281, {
    rbinom(nrow(dat), 1, plogis(0.8 * dat$x1 - 0.4 * dat$x2))
  })
  fml <- y ~ x1 + x2 + x3
  mm <- model.matrix(fml, dat_b)
  beta_hat <- coef(glm(fml, family = binomial(), data = dat_b))
  beta <- beta_hat + withr::with_seed(2282, {
    matrix(rnorm(ncol(mm) * 60, sd = 0.1), ncol(mm))
  })
  refm_lat <- suppressWarnings(init_refmodel(
    list(beta = beta, data = dat_b), data = dat_b, formula = fml,
    family = binomial(),
    ref_predfun = function(fit, newdata = NULL) {
      nd <- if (is.null(newdata)) fit$data else newdata
      unname(model.matrix(delete.response(terms(fml)), nd) %*% fit$beta)
    },
    latent = TRUE, dis = rep(1, 60)
  ))
  vs_lat <- varsel(refm_lat, nclusters = 3, ndraws_pred = 30)
  expect_true(vs_lat$full_ref)
  expect_warning(prj_lat <- project(vs_lat, nterms = 3, ndraws = 30),
                 "returns 30 draws of the reference model")
  draws <- prj_lat$outdmin$draws
  eta_lat <- refm_lat$ref_predfun(refm_lat$fit)
  expect_equal(unname(proj_linpred(prj_lat)$pred),
               unname(t(eta_lat[, draws])))
  expect_equal(unname(proj_epred(prj_lat)),
               unname(t(plogis(eta_lat[, draws]))))
  expect_identical(dim(proj_predict(prj_lat, .seed = 1)),
                   c(30L, nrow(dat_b)))
})

test_that(paste(
  "the reference model draws of project() give the predictions of an",
  "rstanarm fit with trials and offsets, also for new data"
), {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  withr::local_options(projpred.mssg_ndraws = FALSE,
                       projpred.full_size_ref = NULL)
  dat <- dat_full_ref()
  n <- nrow(dat)
  dat$o <- withr::with_seed(5517, runif(n, -0.5, 0.5))
  dat$s <- withr::with_seed(5518, {
    rbinom(n, 5, plogis(0.8 * dat$x1 - 0.4 * dat$x2 + dat$o))
  })
  dat$f <- 5L - dat$s
  fit <- suppressWarnings(rstanarm::stan_glm(
    cbind(s, f) ~ x1 + x2 + x3, family = binomial(), data = dat,
    offset = dat$o, chains = 1, iter = 400, refresh = 0, seed = 7731
  ))
  vs <- suppressWarnings(varsel(fit, nclusters = 2, ndraws_pred = 20))
  expect_true(vs$full_ref)
  expect_warning(prj <- project(vs, nterms = 3, ndraws = 20),
                 "returns 20 draws of the reference model")
  draws <- prj$outdmin$draws
  # Training data:
  expect_equal(unname(proj_linpred(prj)$pred),
               unname(rstanarm::posterior_linpred(fit)[draws, ]))
  # New data with new offsets and trials:
  nd <- dat[1:6, ]
  offs_new <- seq(-1, 1, length.out = 6)
  pl_new <- proj_linpred(prj, newdata = nd, offsetnew = offs_new,
                         weightsnew = rep(3, 6))
  expect_equal(unname(pl_new$pred),
               unname(rstanarm::posterior_linpred(
                 fit, newdata = nd, offset = offs_new
               )[draws, ]))
  expect_identical(dim(pl_new$lpd), c(20L, 6L))
})
