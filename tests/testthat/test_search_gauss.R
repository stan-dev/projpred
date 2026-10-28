context("forward search from cross-products for Gaussian submodels")

# A Gaussian reference model with simulated draws of the coefficients (no model
# fit is necessary). `fml` gives the predictor terms.
refm_search_gauss <- function(fml, n = 80, ndraws = 200, seed = 6283) {
  set.seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                    x4 = rnorm(n),
                    f = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  dat$y <- 0.8 * dat$x1 - 0.4 * dat$x2 + 0.3 * dat$x1 * dat$x2 +
    c(a = 0, b = 0.5, c = -0.5)[as.character(dat$f)] + rnorm(n)
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
                                 dis = rep(1, ndraws)))
}

fml_search_gauss <- y ~ x1 + x2 + x3 + x4 + f + x1:x2 + poly(x3, 2)

test_that(paste(
  "the cross-product search gives the same results as the search with",
  "draw-wise projections"
), {
  refm <- refm_search_gauss(fml_search_gauss)
  for (nclusters_crr in c(1, 5)) {
    vs_crossprod <- suppressWarnings(varsel(
      refm, nclusters = nclusters_crr, nclusters_pred = 7, nterms_max = 6,
      seed = 3317, verbose = FALSE
    ))
    opts_old <- options(projpred.search_gauss_crossprod = FALSE)
    vs_drawwise <- suppressWarnings(varsel(
      refm, nclusters = nclusters_crr, nclusters_pred = 7, nterms_max = 6,
      seed = 3317, verbose = FALSE
    ))
    options(opts_old)
    expect_identical(ranking(vs_crossprod)$fulldata,
                     ranking(vs_drawwise)$fulldata, info = nclusters_crr)
    expect_equal(performances(vs_crossprod, deltas = FALSE),
                 performances(vs_drawwise, deltas = FALSE),
                 info = nclusters_crr)
    expect_equal(vs_crossprod$search_path$outdmins,
                 vs_drawwise$search_path$outdmins, info = nclusters_crr)
  }
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
  cv_drawwise <- suppressWarnings(do.call(cv_varsel, args_cv))
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
