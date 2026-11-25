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
    ref_predfun = ref_predfun, dis = rep(1, ndraws), cvfun = cvfun,
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
  # project() with `refit_prj = FALSE` at that size:
  res <- with_full_prj_count(project(vs, nterms = 3, refit_prj = FALSE), 3L)
  expect_identical(res$n, 1L)
  expect_s3_class(res$value, "projection")
  expect_length(res$value$outdmin, 3)
  expect_identical(predictor_terms(res$value), vs$predictor_ranking)
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
  # project() at the full size gives a projection (and the warning):
  expect_warning(prj <- project(vs, nterms = 3, ndraws = 20),
                 "projection onto the full model")
  expect_s3_class(prj, "projection")
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
