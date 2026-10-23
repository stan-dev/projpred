context("submodel predictions for new data with data-dependent bases")

# A Gaussian reference model with a term that has a data-dependent basis (default
# poly()) and simulated draws of the coefficients. Its predictions for new data
# use the terms of the training model frame (with "predvars"), so the reference
# model itself has no problem with such terms.
refm_predvars <- function(fml = y ~ x1 + x2 + poly(x3, 2), n = 80,
                          ndraws = 200, seed = 6283) {
  set.seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  dat$y <- 0.8 * dat$x1 + 0.5 * dat$x3^2 + rnorm(n)
  mf <- model.frame(fml, dat)
  tt <- delete.response(attr(mf, "terms"))
  mm <- model.matrix(tt, mf)
  beta <- drop(solve(crossprod(mm), crossprod(mm, dat$y))) +
    matrix(rnorm(ncol(mm) * ndraws, sd = 0.05), ncol(mm))
  fit <- list(beta = beta, tt = tt, data = dat)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(fit$tt, model.frame(fit$tt, nd)) %*% fit$beta
  }
  suppressWarnings(init_refmodel(fit, data = dat, formula = fml,
                                 family = gaussian(), ref_predfun = ref_predfun,
                                 dis = rep(1, ndraws)))
}

test_that(paste(
  "projected predictions for new data with a poly() term equal the",
  "predictions for the same training observations"
), {
  refm <- refm_predvars()
  prj <- project(refm, predictor_terms = c("x1", "poly(x3, 2)"), ndraws = 20,
                 seed = 3906, verbose = FALSE)
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = refm$fetch_data()[1:10, ])$pred
  expect_equal(unname(pred_new), unname(pred_train[, 1:10]), tolerance = 1e-10)
})

test_that(paste(
  "projected predictions for new data with an ns() term equal the",
  "predictions for the same training observations"
), {
  skip_if_not_installed("splines")
  refm <- refm_predvars(y ~ x1 + x2 + splines::ns(x3, df = 3))
  prj <- project(refm, predictor_terms = c("x1", "splines::ns(x3, df = 3)"),
                 ndraws = 20, seed = 3906, verbose = FALSE)
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = refm$fetch_data()[1:10, ])$pred
  expect_equal(unname(pred_new), unname(pred_train[, 1:10]), tolerance = 1e-10)
})

test_that(paste(
  "predictions of L1-search submodels for new data with a poly() term equal",
  "the predictions for the same training observations"
), {
  refm <- refm_predvars()
  vs <- suppressWarnings(varsel(refm, method = "L1", refit_prj = FALSE,
                                nterms_max = 3, seed = 3906, verbose = FALSE))
  prj <- suppressWarnings(project(vs, nterms = 3, refit_prj = FALSE,
                                  verbose = FALSE))
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = refm$fetch_data()[1:10, ])$pred
  expect_equal(unname(drop(pred_new)), unname(drop(pred_train)[1:10]),
               tolerance = 1e-10)
})

test_that("the LOO CV with validated search works with a poly() term", {
  refm <- refm_predvars()
  cv <- suppressWarnings(cv_varsel(
    refm, cv_method = "LOO", validate_search = TRUE, nloo = 5, nclusters = 3,
    ndraws_pred = 20, nterms_max = 3, seed = 3906, verbose = FALSE
  ))
  expect_length(ranking(cv)$fulldata, 3)
})
