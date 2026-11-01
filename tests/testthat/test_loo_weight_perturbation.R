context("validated LOO: weight perturbation")

# A Gaussian reference model with fixed draws (fast, no MCMC).
cre_ref_wp <- function(n = 30, p = 6, S = 40) {
  set.seed(4711)
  x <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  dat <- data.frame(y = x[, 1] - 0.5 * x[, 2] + rnorm(n), x)
  beta <- matrix(rnorm((p + 1) * S, sd = 0.2), p + 1, S)
  beta[2, ] <- beta[2, ] + 1
  beta[3, ] <- beta[3, ] - 0.5
  ref_predfun <- function(fit, newdata = NULL) {
    if (is.null(newdata)) newdata <- dat
    cbind(1, as.matrix(newdata[, colnames(x), drop = FALSE])) %*% fit$beta
  }
  ref <- init_refmodel(
    object = list(beta = beta), data = dat,
    formula = reformulate(colnames(x), response = "y"), family = gaussian(),
    ref_predfun = ref_predfun, dis = exp(rnorm(S, sd = 0.1))
  )
  list(ref = ref, n = n, p = p, S = S)
}

run_cv_wp <- function(ref, alpha = 10, refit_prj = TRUE, ...) {
  suppressWarnings(cv_varsel(
    ref, cv_method = "LOO", validate_search = TRUE, method = "forward",
    nclusters = NULL, ndraws = 10, ndraws_pred = 40, refit_prj = refit_prj,
    loo_perturb_alpha = alpha, seed = 11, verbose = FALSE, ...
  ))
}

test_that("perturb_weights() sums to n and does not change the RNG state", {
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
  set.seed(1)
  rng_before <- .Random.seed
  w <- perturb_weights(50, 10, 123L)
  expect_identical(.Random.seed, rng_before)
  expect_equal(sum(w), 50)
  expect_true(all(w > 0))
  # The same weights for another RNG kind (as in the parallel workers):
  RNGkind("L'Ecuyer-CMRG")
  set.seed(2)
  rng_before <- .Random.seed
  expect_identical(perturb_weights(50, 10, 123L), w)
  expect_identical(.Random.seed, rng_before)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
})

test_that("perturb_refmodel() returns NULL for non-finite shifted draws", {
  r <- cre_ref_wp()
  refm <- r$ref
  ll <- t(refm$family$ll_fun(refm$mu_offs, refm$dis, refm$y, refm$wobs))
  lw <- weights(suppressWarnings(loo::psis(-ll, r_eff = NA)))
  i <- 5
  expect_false(is.null(perturb_refmodel(refm, i = i, ll = ll, lw_i = lw[, i],
                                        alpha = 10, seed = 99L)))
  for (val in c(NA, NaN, Inf, -Inf)) {
    refm_v <- refm
    refm_v$eta[3, 2] <- val
    expect_null(perturb_refmodel(refm_v, i = i, ll = ll, lw_i = lw[, i],
                                 alpha = 10, seed = 99L),
                info = format(val))
  }
})

test_that(paste(
  "perturb_refmodel() moves the posterior mean of the linear predictor to the",
  "posterior mean under the weights (first order)"
), {
  r <- cre_ref_wp()
  refm <- r$ref
  ll <- t(refm$family$ll_fun(refm$mu_offs, refm$dis, refm$y, refm$wobs))
  lw <- weights(suppressWarnings(loo::psis(-ll, r_eff = NA)))
  i <- 5
  om <- exp(lw[, i]) / sum(exp(lw[, i]))
  # With a large alpha the weights are close to 1 and the first-order shift
  # gives the exact weighted posterior mean up to O(1 / alpha).
  alpha <- 1e4
  rp <- perturb_refmodel(refm, i = i, ll = ll, lw_i = lw[, i], alpha = alpha,
                         seed = 99L)
  w <- perturb_weights(r$n, alpha, 99L)
  w[i] <- 0
  dw <- w - 1
  dw[i] <- 0
  # Exact importance weights of the weighted posterior relative to the fold
  # posterior:
  lv <- log(om) + drop(ll %*% dw)
  v <- exp(lv - max(lv))
  v <- v / sum(v)
  mean_exact <- drop(refm$eta %*% v)
  mean_fold <- drop(refm$eta %*% om)
  mean_shifted <- drop(rp$eta %*% om)
  expect_lt(max(abs(mean_shifted - mean_exact)),
            0.05 * max(abs(mean_exact - mean_fold)))
  expect_gt(max(abs(mean_exact - mean_fold)), 0)
  expect_equal(rp$mu, rp$eta)
  expect_equal(rp$wobs, refm$wobs * w)
  expect_identical(rp$wobs[i], 0)
})

test_that(paste(
  "weight perturbation changes the fold-wise searches, but not the",
  "performance evaluation"
), {
  skip_on_cran()
  r <- cre_ref_wp()
  cvvs <- run_cv_wp(r$ref)
  cvvs_off <- run_cv_wp(r$ref, alpha = Inf)
  rk <- ranking(cvvs)[["foldwise"]]
  rk_off <- ranking(cvvs_off)[["foldwise"]]
  expect_false(identical(rk, rk_off))
  # Reusing the fold-wise rankings runs only the performance evaluation. With
  # the weight perturbation turned off there, it gives the same results.
  cvvs_eval <- run_cv_wp(r$ref, alpha = Inf, search_out = list(
    search_path = cvvs$search_path, rk_foldwise = rk
  ))
  expect_identical(ranking(cvvs_eval)[["foldwise"]], rk)
  expect_equal(cvvs_eval$summaries$sub, cvvs$summaries$sub)
  # The same seed gives the same result:
  expect_equal(run_cv_wp(r$ref)$summaries$sub, cvvs$summaries$sub)
})

test_that("weight perturbation is off with `refit_prj = FALSE`", {
  skip_on_cran()
  r <- cre_ref_wp()
  expect_identical(
    ranking(run_cv_wp(r$ref, refit_prj = FALSE))[["foldwise"]],
    ranking(run_cv_wp(r$ref, alpha = Inf, refit_prj = FALSE))[["foldwise"]]
  )
})

test_that("weight perturbation works with `nloo < n`", {
  skip_on_cran()
  r <- cre_ref_wp()
  cvvs <- run_cv_wp(r$ref, nloo = 10)
  cvvs_off <- run_cv_wp(r$ref, alpha = Inf, nloo = 10)
  expect_identical(cvvs$loo_inds, cvvs_off$loo_inds)
  expect_false(identical(ranking(cvvs)[["foldwise"]],
                         ranking(cvvs_off)[["foldwise"]]))
})

test_that("weight perturbation gives the same result in parallel", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("doRNG")
  r <- cre_ref_wp()
  cvvs <- run_cv_wp(r$ref)
  doParallel::registerDoParallel(2)
  on.exit(doParallel::stopImplicitCluster())
  cvvs_prll <- run_cv_wp(r$ref, parallel = TRUE)
  expect_identical(ranking(cvvs_prll)[["foldwise"]],
                   ranking(cvvs)[["foldwise"]])
  expect_equal(cvvs_prll$summaries$sub, cvvs$summaries$sub)
})

test_that(paste(
  "argument `loo_perturb_alpha` is checked, stored, and reused by",
  "cv_varsel.vsel()"
), {
  skip_on_cran()
  r <- cre_ref_wp()
  expect_error(run_cv_wp(r$ref, alpha = -1), "must be a single positive")
  expect_error(run_cv_wp(r$ref, alpha = "a"), "must be a single positive")
  cvvs <- run_cv_wp(r$ref, alpha = 4)
  expect_identical(cvvs$loo_perturb_alpha, 4)
  cvvs2 <- suppressWarnings(cv_varsel(cvvs, validate_search = FALSE,
                                      verbose = FALSE))
  expect_identical(cvvs2$loo_perturb_alpha, 4)
  # Another `loo_perturb_alpha` runs the fold-wise searches again:
  expect_message(
    cvvs3 <- suppressWarnings(cv_varsel(
      cvvs, loo_perturb_alpha = 10, ndraws_pred = 40, seed = 11, verbose = FALSE
    )),
    "cannot be re-used"
  )
  expect_identical(cvvs3$loo_perturb_alpha, 10)
  expect_identical(ranking(cvvs3)[["foldwise"]],
                   ranking(run_cv_wp(r$ref, alpha = 10))[["foldwise"]])
})
