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
