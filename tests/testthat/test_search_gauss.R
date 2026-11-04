context("forward search from cross-products for Gaussian submodels")

# A Gaussian reference model with simulated draws of the coefficients (no model
# fit is necessary). `fml` gives the predictor terms.
refm_search_gauss <- function(fml, n = 80, ndraws = 200, seed = 6283,
                              weighted = FALSE) {
  set.seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                    x4 = rnorm(n),
                    f = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  # A column that is collinear with `x1` and `x2`:
  dat$x5 <- dat$x1 + dat$x2
  dat$y <- 0.8 * dat$x1 - 0.4 * dat$x2 + 0.3 * dat$x1 * dat$x2 +
    c(a = 0, b = 0.5, c = -0.5)[as.character(dat$f)] + rnorm(n)
  extract_model_data <- NULL
  if (weighted) {
    # Observation weights from the column `w`:
    dat$w <- runif(n, 0.5, 2)
    extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                   extract_y = TRUE) {
      y_wobs_offs(newdata = newdata, wrhs = wrhs %||% ~ w, orhs = orhs,
                  resp_form = if (extract_y) ~ y else NULL)
    }
  }
  mm <- model.matrix(fml, dat)
  beta_hat <- drop(solve(crossprod(mm) + diag(1e-6, ncol(mm)),
                         crossprod(mm, dat$y)))
  beta <- beta_hat + matrix(rnorm(ncol(mm) * ndraws, sd = 0.05), ncol(mm))
  fit <- list(beta = beta, data = dat, fml = fml)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(delete.response(terms(fit$fml)), nd) %*% fit$beta
  }
  suppressWarnings(init_refmodel(fit, data = dat, formula = fml,
                                 family = gaussian(), ref_predfun = ref_predfun,
                                 extract_model_data = extract_model_data,
                                 dis = rep(1, ndraws)))
}

fml_search_gauss <- y ~ x1 + x2 + x3 + x4 + f + x1:x2 + poly(x3, 2)

# Evaluates `expr` and counts the calls of search_forward_gauss() that returned
# a search result (not `NULL`). Returns the value of `expr` and the count.
with_gauss_search_count <- function(expr) {
  calls <- new.env()
  calls$n <- 0L
  ns <- asNamespace("projpred")
  suppressMessages(trace(
    "search_forward_gauss", where = ns, print = FALSE,
    exit = bquote(if (!is.null(returnValue(NULL))) {
      assign("n", get("n", envir = .(calls)) + 1L, envir = .(calls))
    })
  ))
  on.exit(suppressMessages(untrace("search_forward_gauss", where = ns)),
          add = TRUE)
  value <- expr
  return(list(value = value, n = calls$n))
}

test_that(paste(
  "the cross-product search gives the same results as the search with",
  "draw-wise projections"
), {
  for (weighted in c(FALSE, TRUE)) {
    refm <- refm_search_gauss(fml_search_gauss, weighted = weighted)
    expect_identical(all(refm$wobs == 1), !weighted)
    for (nclusters_crr in c(1, 5)) {
      info_crr <- paste("weighted:", weighted, "nclusters:", nclusters_crr)
      res_crossprod <- with_gauss_search_count(suppressWarnings(varsel(
        refm, nclusters = nclusters_crr, nclusters_pred = 7, nterms_max = 6,
        seed = 3317, verbose = FALSE
      )))
      expect_identical(res_crossprod$n, 1L, info = info_crr)
      vs_crossprod <- res_crossprod$value
      opts_old <- options(projpred.search_gauss_crossprod = FALSE)
      on.exit(options(opts_old), add = TRUE)
      # The draw-wise reference fits each draw with glm_ridge(), not with the
      # joint Gaussian solve:
      vs_drawwise <- with_mocked_bindings(
        suppressWarnings(varsel(
          refm, nclusters = nclusters_crr, nclusters_pred = 7, nterms_max = 6,
          seed = 3317, verbose = FALSE
        )),
        glm_ridge_gauss_multi = function(...) NULL
      )
      options(opts_old)
      expect_identical(ranking(vs_crossprod)$fulldata,
                       ranking(vs_drawwise)$fulldata, info = info_crr)
      expect_equal(performances(vs_crossprod, deltas = FALSE),
                   performances(vs_drawwise, deltas = FALSE),
                   info = info_crr)
      expect_equal(vs_crossprod$search_path$outdmins,
                   vs_drawwise$search_path$outdmins, info = info_crr)
    }
  }
})

test_that(paste(
  "the cross-product search gives the same results as the search with",
  "draw-wise projections for custom search terms"
), {
  refm <- refm_search_gauss(fml_search_gauss)
  search_terms <- force_search_terms(forced_terms = c("x1", "x2"),
                                     optional_terms = c("x3", "x4", "f"))
  args_vs <- list(refm, search_terms = search_terms, nclusters = 5,
                  nclusters_pred = 7, nterms_max = 4, seed = 3317,
                  verbose = FALSE)
  res_crossprod <- with_gauss_search_count(
    suppressWarnings(do.call(varsel, args_vs))
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
})

test_that(paste(
  "the LOO CV with validated search gives the same results with and without",
  "the cross-product computations"
), {
  refm <- refm_search_gauss(fml_search_gauss, n = 40)
  args_cv <- list(refm, cv_method = "LOO", validate_search = TRUE, nloo = 10,
                  nclusters = 3, ndraws_pred = 20, nterms_max = 4, seed = 5821,
                  verbose = FALSE)
  cv_crossprod <- suppressWarnings(do.call(cv_varsel, args_cv))
  opts_old <- options(projpred.search_gauss_crossprod = FALSE)
  on.exit(options(opts_old), add = TRUE)
  cv_drawwise <- with_mocked_bindings(
    suppressWarnings(do.call(cv_varsel, args_cv)),
    glm_ridge_gauss_multi = function(...) NULL
  )
  options(opts_old)
  expect_identical(ranking(cv_crossprod), ranking(cv_drawwise))
  expect_equal(cv_crossprod$summaries$sub, cv_drawwise$summaries$sub,
               tolerance = 1e-10)
  expect_equal(performances(cv_crossprod, deltas = FALSE),
               performances(cv_drawwise, deltas = FALSE), tolerance = 1e-10)
})

test_that("perf_eval_gauss() applies in the LOO CV with validated search", {
  refm <- refm_search_gauss(fml_search_gauss, n = 40)
  n_used <- 0L
  trace("perf_eval_gauss", exit = function() {
    if (!is.null(returnValue())) n_used <<- n_used + 1L
  }, print = FALSE, where = asNamespace("projpred"))
  on.exit(untrace("perf_eval_gauss", where = asNamespace("projpred")))
  cv <- suppressWarnings(cv_varsel(
    refm, cv_method = "LOO", validate_search = TRUE, nloo = 5, nclusters = 3,
    ndraws_pred = 20, nterms_max = 3, seed = 5821, verbose = FALSE
  ))
  expect_identical(n_used, 5L)
})

test_that("search_forward_gauss() applies to Gaussian GLM-type submodels", {
  refm <- refm_search_gauss(fml_search_gauss)
  p_ref <- get_refdist(refm, nclusters = 3)
  search_terms <- union("1", split_formula(refm$formula,
                                           data = refm$fetch_data()))
  out <- search_forward_gauss(p_ref = p_ref, refmodel = refm, nterms_max = 4,
                              verbose = FALSE, search_terms = search_terms)
  expect_false(is.null(out))
  expect_length(out$predictor_ranking, 4)
  expect_length(out$outdmins, 5)
})

test_that(paste(
  "search_forward_gauss() with `return_outdmins = FALSE` gives the same",
  "ranking and no submodel fits"
), {
  refm <- refm_search_gauss(fml_search_gauss)
  p_ref <- get_refdist(refm, nclusters = 3)
  search_terms <- union("1", split_formula(refm$formula,
                                           data = refm$fetch_data()))
  out <- search_forward_gauss(p_ref = p_ref, refmodel = refm, nterms_max = 4,
                              verbose = FALSE, search_terms = search_terms)
  out_nofits <- search_forward_gauss(p_ref = p_ref, refmodel = refm,
                                     nterms_max = 4, verbose = FALSE,
                                     search_terms = search_terms,
                                     return_outdmins = FALSE)
  expect_identical(out_nofits$predictor_ranking, out$predictor_ranking)
  expect_null(out_nofits$outdmins)
})

test_that(paste(
  "the search of each left-out observation in the LOO CV with validated",
  "search does not compute the submodel fits along the ranking"
), {
  refm <- refm_search_gauss(y ~ x1 + x2 + x3 + x4 + f + x1:x2, n = 40)
  calls <- new.env()
  calls$n_nofits <- 0L
  ns <- asNamespace("projpred")
  suppressMessages(trace(
    "search_forward_gauss", where = ns, print = FALSE,
    tracer = bquote(if (!return_outdmins) {
      assign("n_nofits", get("n_nofits", envir = .(calls)) + 1L,
             envir = .(calls))
    })
  ))
  on.exit(suppressMessages(untrace("search_forward_gauss", where = ns)),
          add = TRUE)
  cv <- suppressWarnings(cv_varsel(
    refm, cv_method = "LOO", validate_search = TRUE, nloo = 5, nclusters = 3,
    ndraws_pred = 20, nterms_max = 3, seed = 5821, verbose = FALSE
  ))
  # One search per left-out observation; the search on the full data keeps
  # its submodel fits:
  expect_identical(calls$n_nofits, 5L)
})

test_that(paste(
  "the cross-product search and evaluation apply to the latent projection",
  "and give the same results as the draw-wise projections"
), {
  set.seed(7391)
  n <- 50
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                    x4 = rnorm(n))
  eta <- 0.9 * dat$x1 - 0.6 * dat$x2
  dat$y <- rbinom(n, 1, plogis(eta))
  ndraws <- 100
  mm <- model.matrix(~ x1 + x2 + x3 + x4, dat)
  beta <- c(0, 0.9, -0.6, 0, 0) +
    matrix(rnorm(ncol(mm) * ndraws, sd = 0.1), ncol(mm))
  fit <- list(beta = beta, data = dat)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(~ x1 + x2 + x3 + x4, nd) %*% fit$beta
  }
  refm <- suppressWarnings(init_refmodel(
    fit, data = dat, formula = y ~ x1 + x2 + x3 + x4, family = binomial(),
    ref_predfun = ref_predfun, dis = rep(1, ndraws), latent = TRUE
  ))
  expect_true(refm$family$for_latent)
  calls <- new.env()
  calls$n <- 0L
  suppressMessages(trace(
    "search_forward_gauss", where = asNamespace("projpred"), print = FALSE,
    exit = bquote(if (!is.null(returnValue(NULL))) {
      assign("n", get("n", envir = .(calls)) + 1L, envir = .(calls))
    })
  ))
  on.exit(suppressMessages(untrace("search_forward_gauss",
                                   where = asNamespace("projpred"))),
          add = TRUE)
  args_vs <- list(refm, nclusters = 3, nclusters_pred = 5, nterms_max = 3,
                  seed = 4410, verbose = FALSE)
  args_cv <- c(args_vs, list(cv_method = "LOO", validate_search = TRUE,
                             nloo = 5))
  vs_crossprod <- suppressWarnings(do.call(varsel, args_vs))
  cv_crossprod <- suppressWarnings(do.call(cv_varsel, args_cv))
  opts_old <- options(projpred.search_gauss_crossprod = FALSE)
  on.exit(options(opts_old), add = TRUE)
  vs_drawwise <- suppressWarnings(do.call(varsel, args_vs))
  cv_drawwise <- suppressWarnings(do.call(cv_varsel, args_cv))
  options(opts_old)
  expect_identical(ranking(vs_crossprod)$fulldata,
                   ranking(vs_drawwise)$fulldata)
  expect_equal(performances(vs_crossprod, deltas = FALSE),
               performances(vs_drawwise, deltas = FALSE))
  expect_identical(ranking(cv_crossprod), ranking(cv_drawwise))
  expect_equal(performances(cv_crossprod, deltas = FALSE),
               performances(cv_drawwise, deltas = FALSE), tolerance = 1e-10)
  # The cross-product search ran: once in varsel(), once for the full data in
  # cv_varsel() and once per left-out observation:
  expect_identical(calls$n, 7L)
})

test_that(paste(
  "search_forward_gauss() returns `NULL` for a factor in an interaction term"
), {
  refm <- refm_search_gauss(y ~ x1 + x2 + f + x1:f)
  p_ref <- get_refdist(refm, nclusters = 3)
  search_terms <- union("1", split_formula(refm$formula,
                                           data = refm$fetch_data()))
  out <- search_forward_gauss(p_ref = p_ref, refmodel = refm, nterms_max = 3,
                              verbose = FALSE, search_terms = search_terms)
  expect_null(out)
})

test_that("gauss_crossprod_regul() requires a positive `regul`", {
  refm <- refm_search_gauss(y ~ x1 + x2)
  expect_identical(gauss_crossprod_regul(refm, list()), 1e-4)
  expect_identical(gauss_crossprod_regul(refm, list(regul = 0.5)), 0.5)
  for (regul_crr in list(0, -1e-3, NA_real_)) {
    expect_null(gauss_crossprod_regul(refm, list(regul = regul_crr)),
                info = paste("regul:", regul_crr))
  }
})

test_that(paste(
  "a `regul` of zero or below zero uses the draw-wise path for collinear",
  "columns"
), {
  # With a `regul` of zero or below zero, both arms below use the same
  # draw-wise code. The test checks that the cross-product search does not run
  # and that no error occurs.
  refm <- refm_search_gauss(y ~ x1 + x2 + x3 + x5, n = 40)
  for (regul_crr in c(0, -1e-3)) {
    info_crr <- paste("regul:", regul_crr)
    args_vs <- list(refm, regul = regul_crr, nclusters = 3, nclusters_pred = 5,
                    nterms_max = 4, seed = 3317, verbose = FALSE)
    args_cv <- c(args_vs, list(cv_method = "LOO", validate_search = TRUE,
                               nloo = 5))
    res_vs <- with_gauss_search_count(
      suppressWarnings(do.call(varsel, args_vs))
    )
    expect_identical(res_vs$n, 0L, info = info_crr)
    vs_crossprod <- res_vs$value
    cv_crossprod <- suppressWarnings(do.call(cv_varsel, args_cv))
    opts_old <- options(projpred.search_gauss_crossprod = FALSE)
    on.exit(options(opts_old), add = TRUE)
    vs_drawwise <- suppressWarnings(do.call(varsel, args_vs))
    cv_drawwise <- suppressWarnings(do.call(cv_varsel, args_cv))
    options(opts_old)
    expect_identical(ranking(vs_crossprod)$fulldata,
                     ranking(vs_drawwise)$fulldata, info = info_crr)
    expect_equal(performances(vs_crossprod, deltas = FALSE),
                 performances(vs_drawwise, deltas = FALSE), info = info_crr)
    expect_identical(ranking(cv_crossprod), ranking(cv_drawwise),
                     info = info_crr)
    expect_equal(performances(cv_crossprod, deltas = FALSE),
                 performances(cv_drawwise, deltas = FALSE), info = info_crr)
  }
})

test_that(paste(
  "the cross-product search and evaluation return `NULL` for a matrix that",
  "is not positive definite"
), {
  # With `x5 = x1 + x2` and a negative `regul`, the matrix of the columns x1,
  # x2 and x5 is not positive definite:
  refm <- refm_search_gauss(y ~ x1 + x2 + x3 + x5, n = 40)
  p_ref <- get_refdist(refm, nclusters = 3)
  calls <- new.env()
  calls$n_all <- 0L
  calls$n_used <- 0L
  ns <- asNamespace("projpred")
  suppressMessages(trace(
    "perf_eval_gauss", where = ns, print = FALSE,
    exit = bquote({
      assign("n_all", get("n_all", envir = .(calls)) + 1L, envir = .(calls))
      if (!is.null(returnValue(NULL))) {
        assign("n_used", get("n_used", envir = .(calls)) + 1L,
               envir = .(calls))
      }
    })
  ))
  on.exit(suppressMessages(untrace("perf_eval_gauss", where = ns)),
          add = TRUE)
  local_mocked_bindings(gauss_crossprod_regul = function(...) -1)
  # Candidates with one column (crit_single()):
  search_terms <- union("1", split_formula(refm$formula,
                                           data = refm$fetch_data()))
  expect_null(search_forward_gauss(
    p_ref = p_ref, refmodel = refm, nterms_max = 4, verbose = FALSE,
    search_terms = search_terms
  ))
  # A candidate with several columns (crit_cand()):
  expect_null(search_forward_gauss(
    p_ref = p_ref, refmodel = refm, nterms_max = 4, verbose = FALSE,
    search_terms = c("1", "x3", "x1 + x2 + x5")
  ))
  cv <- suppressWarnings(cv_varsel(
    refm, cv_method = "LOO", validate_search = TRUE, nloo = 3, nclusters = 3,
    ndraws_pred = 20, nterms_max = 4, seed = 5821, verbose = FALSE
  ))
  expect_s3_class(cv, "vsel")
  # One evaluation per left-out observation, each on the draw-wise path:
  expect_identical(calls$n_all, 3L)
  expect_identical(calls$n_used, 0L)
})

test_that(paste(
  "terms with `+` inside a call work in the search and the evaluation, with",
  "and without the cross-product computations"
), {
  # A split of the term strings at each `+` gives the part "x3)" twice:
  refm <- refm_search_gauss(y ~ I(x1 + x3) + I(x2 + x3) + x4, n = 40)
  args_vs <- list(refm, nclusters = 3, nclusters_pred = 5, seed = 3317,
                  verbose = FALSE)
  args_cv <- c(args_vs, list(cv_method = "LOO", validate_search = TRUE,
                             nloo = 5))
  n_used <- 0L
  trace("perf_eval_gauss", exit = function() {
    if (!is.null(returnValue())) n_used <<- n_used + 1L
  }, print = FALSE, where = asNamespace("projpred"))
  on.exit(untrace("perf_eval_gauss", where = asNamespace("projpred")),
          add = TRUE)
  res_vs <- with_gauss_search_count(suppressWarnings(do.call(varsel, args_vs)))
  expect_identical(res_vs$n, 1L)
  expect_setequal(ranking(res_vs$value)$fulldata,
                  c("I(x1 + x3)", "I(x2 + x3)", "x4"))
  cv_crossprod <- suppressWarnings(do.call(cv_varsel, args_cv))
  expect_identical(n_used, 5L)
  opts_old <- options(projpred.search_gauss_crossprod = FALSE)
  on.exit(options(opts_old), add = TRUE)
  vs_drawwise <- suppressWarnings(do.call(varsel, args_vs))
  cv_drawwise <- suppressWarnings(do.call(cv_varsel, args_cv))
  options(opts_old)
  expect_identical(ranking(res_vs$value)$fulldata,
                   ranking(vs_drawwise)$fulldata)
  expect_equal(performances(res_vs$value, deltas = FALSE),
               performances(vs_drawwise, deltas = FALSE))
  expect_identical(ranking(cv_crossprod), ranking(cv_drawwise))
  expect_equal(performances(cv_crossprod, deltas = FALSE),
               performances(cv_drawwise, deltas = FALSE), tolerance = 1e-10)
})

test_that(paste(
  "the cross-product search and evaluation find functions from the",
  "environment of the reference model formula"
), {
  # `sq_local()` exists only in the environment of the formula:
  make_refm_local <- function() {
    sq_local <- function(x) x^2
    refm_search_gauss(y ~ x1 + sq_local(x2) + x3 + x4, n = 40)
  }
  refm <- make_refm_local()
  expect_false(exists("sq_local"))
  args_vs <- list(refm, nclusters = 3, nclusters_pred = 5, nterms_max = 4,
                  seed = 3317, verbose = FALSE)
  args_cv <- c(args_vs, list(cv_method = "LOO", validate_search = TRUE,
                             nloo = 5))
  n_used <- 0L
  trace("perf_eval_gauss", exit = function() {
    if (!is.null(returnValue())) n_used <<- n_used + 1L
  }, print = FALSE, where = asNamespace("projpred"))
  on.exit(untrace("perf_eval_gauss", where = asNamespace("projpred")),
          add = TRUE)
  res_vs <- with_gauss_search_count(suppressWarnings(do.call(varsel, args_vs)))
  expect_identical(res_vs$n, 1L)
  cv_crossprod <- suppressWarnings(do.call(cv_varsel, args_cv))
  expect_identical(n_used, 5L)
  opts_old <- options(projpred.search_gauss_crossprod = FALSE)
  on.exit(options(opts_old), add = TRUE)
  vs_drawwise <- suppressWarnings(do.call(varsel, args_vs))
  cv_drawwise <- suppressWarnings(do.call(cv_varsel, args_cv))
  options(opts_old)
  expect_identical(ranking(res_vs$value)$fulldata,
                   ranking(vs_drawwise)$fulldata)
  expect_equal(performances(res_vs$value, deltas = FALSE),
               performances(vs_drawwise, deltas = FALSE))
  expect_identical(ranking(cv_crossprod), ranking(cv_drawwise))
  expect_equal(performances(cv_crossprod, deltas = FALSE),
               performances(cv_drawwise, deltas = FALSE), tolerance = 1e-10)
})
