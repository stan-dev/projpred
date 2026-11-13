withr::local_seed(4710)

# A Gaussian reference model with fixed draws (fast, no MCMC).
cre_ref_wp <- function(n = 30, p = 6, S = 40) {
  withr::local_seed(4711)
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
  # The test sets the seed and the RNG kind; restore both at the end:
  withr::local_preserve_seed()
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
  "perturb_refmodel() returns NULL for an infinite linear predictor with a",
  "bounded inverse link and for infinite means with offsets"
), {
  withr::local_seed(4712)
  n <- 30
  S <- 40
  x <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("x1", "x2")))
  beta <- matrix(rnorm(3 * S, sd = 0.2), 3, S)
  cre_ref <- function(dat, family) {
    ref_predfun <- function(fit, newdata = NULL) {
      if (is.null(newdata)) newdata <- dat
      cbind(1, as.matrix(newdata[, colnames(x), drop = FALSE])) %*% fit$beta
    }
    suppressWarnings(init_refmodel(
      object = list(beta = beta), data = dat, formula = y ~ x1 + x2,
      family = family, ref_predfun = ref_predfun
    ))
  }
  i <- 5
  perturb <- function(refm) {
    ll <- t(refm$family$ll_fun(refm$mu_offs, refm$dis, refm$y, refm$wobs))
    lw <- weights(suppressWarnings(loo::psis(-ll, r_eff = NA)))
    list(ll = ll, lw_i = lw[, i])
  }
  # Logit link: an infinite value in `eta` (set here, not from an overflow of
  # the shift) gives `mu` 0 or 1. The shift of row k is eta[k, ] %*% cd. The
  # test puts the value into a draw s with cd_s > 0, so the shifted row is +Inf
  # and not NaN (a NaN would also make `mu` non-finite).
  refm <- cre_ref(data.frame(y = rbinom(n, 1, 0.5), x), binomial())
  pp <- perturb(refm)
  expect_false(is.null(perturb_refmodel(refm, i = i, ll = pp$ll,
                                        lw_i = pp$lw_i, alpha = 10,
                                        seed = 99L)))
  # `cd` as in perturb_refmodel():
  w <- perturb_weights(n, 10, 99L)
  w[i] <- 0
  dw <- w - 1
  dw[i] <- 0
  om <- exp(pp$lw_i) / sum(exp(pp$lw_i))
  cc <- drop(pp$ll %*% dw)
  cd <- om * (cc - sum(om * cc))
  s <- which.max(cd)
  expect_gt(cd[s], 0)
  refm_v <- refm
  refm_v$eta[3, s] <- Inf
  eta_shift <- refm_v$eta + drop(refm_v$eta %*% cd)
  expect_true(all(eta_shift[3, ] == Inf))
  expect_true(all(is.finite(refm_v$family$linkinv(eta_shift))))
  expect_null(perturb_refmodel(refm_v, i = i, ll = pp$ll, lw_i = pp$lw_i,
                               alpha = 10, seed = 99L))
  # Log link with an offset: `mu` is finite, but `mu_offs` is not.
  refm <- cre_ref(data.frame(y = rpois(n, 2), x), poisson())
  pp <- perturb(refm)
  refm_v <- refm
  refm_v$offset <- rep(0.1, n)
  expect_false(is.null(perturb_refmodel(refm_v, i = i, ll = pp$ll,
                                        lw_i = pp$lw_i, alpha = 10,
                                        seed = 99L)))
  refm_v$offset[3] <- 1000
  expect_null(perturb_refmodel(refm_v, i = i, ll = pp$ll, lw_i = pp$lw_i,
                               alpha = 10, seed = 99L))
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
  # With a large alpha the weights are close to 1, and the first-order shift
  # gives the weighted posterior mean with an error of O(1 / alpha).
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
  # Identity link without offsets: the draws are not shifted, and
  # `draws_shift` holds the shift.
  expect_identical(rp$eta, refm$eta)
  expect_false(is.null(rp$draws_shift))
  mean_shifted <- drop((rp$eta + rp$draws_shift) %*% om)
  expect_lt(max(abs(mean_shifted - mean_exact)),
            0.05 * max(abs(mean_exact - mean_fold)))
  expect_gt(max(abs(mean_exact - mean_fold)), 0)
  expect_identical(rp$mu, rp$eta)
  expect_identical(rp$mu_offs, rp$eta)
  expect_equal(rp$wobs, refm$wobs * w)
  expect_identical(rp$wobs[i], 0)
  # The shift is the first-order change of the draws:
  cc <- drop(ll %*% dw)
  cd <- om * (cc - sum(om * cc))
  expect_equal(rp$draws_shift, drop(refm$eta %*% cd), tolerance = 1e-14)
  # With offsets, the draws are shifted, and there is no `draws_shift`:
  refm_o <- refm
  refm_o$offset <- rep(0.1, r$n)
  rp_o <- perturb_refmodel(refm_o, i = i, ll = ll, lw_i = lw[, i],
                           alpha = alpha, seed = 99L)
  expect_null(rp_o$draws_shift)
  expect_identical(rp_o$eta, refm$eta + rp$draws_shift)
  expect_identical(rp_o$mu_offs, rp_o$eta + 0.1)
  expect_identical(rp_o$dis, rp$dis)
})

# Evaluates `expr` and counts the calls of search_forward_gauss() and
# perf_eval_gauss() that returned a result (not `NULL`).
with_gauss_counts <- function(expr) {
  calls <- new.env()
  calls$search <- 0L
  calls$eval <- 0L
  ns <- asNamespace("projpred")
  fns <- c(search = "search_forward_gauss", eval = "perf_eval_gauss")
  withr::defer(suppressMessages(for (fn in fns) untrace(fn, where = ns)))
  for (nm in names(fns)) {
    suppressMessages(trace(
      fns[[nm]], where = ns, print = FALSE,
      exit = bquote(if (!is.null(returnValue(NULL))) {
        assign(.(nm), get(.(nm), envir = .(calls)) + 1L, envir = .(calls))
      })
    ))
  }
  value <- expr
  return(list(value = value, search = calls$search, eval = calls$eval))
}

# Evaluates `expr` and counts the get_p_clust() calls that added `draws_shift`
# to the cluster means (and did not shift the draws).
with_shift_count <- function(expr) {
  calls <- new.env()
  calls$n <- 0L
  ns <- asNamespace("projpred")
  withr::defer(suppressMessages(untrace("get_p_clust", where = ns)))
  suppressMessages(trace(
    "get_p_clust", where = ns, print = FALSE,
    exit = bquote(if (!is.null(draws_shift)) {
      assign("n", get("n", envir = .(calls)) + 1L, envir = .(calls))
    })
  ))
  value <- expr
  return(list(value = value, n = calls$n))
}

test_that(paste(
  "get_p_clust() with `draws_shift` gives the results for the shifted",
  "draws"
), {
  r <- cre_ref_wp()
  refm <- r$ref
  S <- r$S
  withr::local_seed(3104)
  shift <- rnorm(r$n)
  w <- exp(rnorm(S))
  w_special <- w
  w_special[c(3, 8)] <- c(0, Inf)
  # A user-defined `predvar` turns off `family$predvar_mix`:
  fam_nomix <- refm$family
  fam_nomix$predvar <- function(mu, dis, wdraws = 1) {
    refm$family$predvar(mu, dis, wdraws)
  }
  cl_several <- rep(1:5, length.out = S)
  # The traced function is in the namespace (devtools::load_all() also puts an
  # untraced copy into the search path):
  get_p_clust_ns <- function(...) {
    get("get_p_clust", envir = asNamespace("projpred"))(...)
  }
  # `n_shift`: the number of calls that add the shift to the cluster means.
  cases <- list(
    several = list(cl = cl_several, w = w, n_shift = 1L),
    several_thinned = list(cl = replace(cl_several, 1:3, NA), w = w,
                           n_shift = 1L),
    student_t = list(cl = cl_several, w = w, n_shift = 1L,
                     family = extend_family(Student_t(nu = 4))),
    latent = list(cl = cl_several, w = w, n_shift = 1L,
                  family = extend_family(binomial(), latent = TRUE)),
    several_special = list(cl = cl_several, w = w_special, n_shift = 0L),
    nan_weight = list(cl = cl_several, w = replace(w, 4, NaN), n_shift = 0L),
    zero_weight_sum = list(cl = c(1, 1, rep(2, S - 2)),
                           w = replace(w, 1:2, 0), n_shift = 0L),
    empty_cluster = list(cl = c(1, rep(3, S - 1)), w = w, n_shift = 0L),
    one_draw = list(cl = seq_len(S), w = w, n_shift = 0L),
    thinned = list(cl = replace(rep(NA, S), c(2, 9, 30), 1:3), w = w,
                   n_shift = 0L),
    no_predvar_mix = list(cl = cl_several, w = w, n_shift = 0L,
                          family = fam_nomix)
  )
  for (case_name in names(cases)) {
    cs <- cases[[case_name]]
    args <- list(family = cs$family %||% refm$family, eta = refm$eta,
                 mu = refm$mu, mu_offs = refm$mu_offs, dis = refm$dis,
                 wdraws = cs$w, cl = cs$cl)
    res <- with_shift_count(do.call(get_p_clust_ns,
                                    c(args, list(draws_shift = shift))))
    expect_identical(res$n, cs$n_shift, info = case_name)
    args_ref <- modifyList(args, list(eta = refm$eta + shift,
                                      mu = refm$mu + shift,
                                      mu_offs = refm$mu_offs + shift))
    expect_equal(res$value, do.call(get_p_clust, args_ref),
                 tolerance = 1e-12, info = case_name)
  }
  # With the clustering in get_p_clust() (no `cl`):
  set.seed(77)
  res <- with_shift_count(get_p_clust_ns(
    refm$family, refm$eta, refm$mu, refm$mu_offs, refm$dis, nclusters = 4,
    draws_shift = shift
  ))
  expect_identical(res$n, 1L)
  set.seed(77)
  ref <- get_p_clust(refm$family, refm$eta + shift, refm$mu + shift,
                     refm$mu_offs + shift, refm$dis, nclusters = 4)
  expect_equal(res$value, ref, tolerance = 1e-12)
  # `as_moments = TRUE` returns the moments for one draw per cluster. With
  # `draws_shift`, it returns `NULL`: the callers of `as_moments` do not use
  # `draws_shift`.
  args_mom <- list(family = refm$family, eta = refm$eta, mu = refm$mu,
                   mu_offs = refm$mu_offs, dis = refm$dis, wdraws = w,
                   cl = seq_len(S), as_moments = TRUE)
  expect_identical(do.call(get_p_clust, args_mom)$moments_only, TRUE)
  expect_null(do.call(get_p_clust, c(args_mom, list(draws_shift = shift))))
})

test_that(paste(
  "LOO validated search with weight perturbation gives the same results with",
  "the shift in get_p_clust() and with shifted draws"
), {
  skip_on_cran()
  ns <- asNamespace("projpred")
  # A latent projection (binomial reference model) and a Gaussian one:
  withr::local_seed(7392)
  n <- 30
  S <- 40
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  dat$y <- rbinom(n, 1, plogis(0.9 * dat$x1 - 0.6 * dat$x2))
  fit <- list(beta = c(0, 0.9, -0.6, 0) +
                matrix(rnorm(4 * S, sd = 0.1), 4), data = dat)
  ref_lat <- suppressWarnings(init_refmodel(
    fit, data = dat, formula = y ~ x1 + x2 + x3, family = binomial(),
    ref_predfun = function(fit, newdata = NULL) {
      model.matrix(~ x1 + x2 + x3, newdata %||% fit$data) %*% fit$beta
    },
    dis = rep(1, S), latent = TRUE
  ))
  expect_true(ref_lat$family$for_latent)
  refs <- list(gaussian = cre_ref_wp()$ref, latent = ref_lat)
  # Several draws per search cluster, so that get_p_clust() uses the shift:
  run_cv <- function(ref) {
    suppressWarnings(cv_varsel(
      ref, cv_method = "LOO", validate_search = TRUE, method = "forward",
      nclusters = 5, ndraws_pred = 40, loo_perturb_alpha = 10, seed = 11,
      verbose = FALSE
    ))
  }
  res_shift <- lapply(refs, function(ref) {
    res <- with_gauss_counts(with_shift_count(run_cv(ref)))
    c(res$value, list(n_gauss = c(res$search, res$eval)))
  })
  # The reference: .select() gets the shifted draws.
  withr::defer(suppressMessages(untrace(".select", where = ns)))
  suppressMessages(trace(".select", where = ns, print = FALSE, tracer = quote({
    if (!is.null(refmodel$draws_shift)) {
      refmodel <- perturb_draws(refmodel, refmodel$draws_shift)
      refmodel$draws_shift <- NULL
    }
  })))
  res_full <- lapply(refs, function(ref) {
    res <- with_gauss_counts(with_shift_count(run_cv(ref)))
    c(res$value, list(n_gauss = c(res$search, res$eval)))
  })
  for (nm in names(refs)) {
    expect_identical(res_shift[[nm]]$n, as.integer(nrow(refs[[nm]]$mu)),
                     info = nm)
    expect_identical(res_full[[nm]]$n, 0L, info = nm)
    # The full-data search and each fold-wise search (with the shifted cluster
    # means) use the cross-product search. The evaluation (without the shift)
    # uses the cross-product evaluation in each fold. Both runs use the same
    # paths:
    nobs <- as.integer(nrow(refs[[nm]]$mu))
    expect_identical(res_shift[[nm]]$n_gauss, c(nobs + 1L, nobs), info = nm)
    expect_identical(res_full[[nm]]$n_gauss, res_shift[[nm]]$n_gauss,
                     info = nm)
    expect_identical(ranking(res_shift[[nm]]$value),
                     ranking(res_full[[nm]]$value), info = nm)
    expect_equal(res_shift[[nm]]$value$summaries$sub,
                 res_full[[nm]]$value$summaries$sub, tolerance = 1e-10,
                 info = nm)
  }
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
  # Re-using the fold-wise rankings runs only the performance evaluation. With
  # `alpha = Inf`, it gives the same summaries as `cvvs`.
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
  withr::defer(doParallel::stopImplicitCluster())
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
  # Without fold-wise searches, no search uses the weights:
  cvvs2 <- suppressWarnings(cv_varsel(cvvs, validate_search = FALSE,
                                      verbose = FALSE))
  expect_identical(cvvs2$loo_perturb_alpha, Inf)
  # The same value as an integer reuses the fold-wise searches:
  expect_no_message(
    cvvs4 <- suppressWarnings(cv_varsel(
      cvvs, loo_perturb_alpha = 4L, ndraws_pred = 40, seed = 11, verbose = FALSE
    )),
    message = "cannot be re-used"
  )
  expect_identical(cvvs4$loo_perturb_alpha, 4)
  expect_identical(ranking(cvvs4)[["foldwise"]], ranking(cvvs)[["foldwise"]])
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

test_that(paste(
  "the `vsel` object stores the `loo_perturb_alpha` that the fold-wise",
  "searches used"
), {
  skip_on_cran()
  r <- cre_ref_wp()
  rk_off <- ranking(run_cv_wp(r$ref, alpha = Inf))[["foldwise"]]
  # With `refit_prj = FALSE`, the searches do not use the weights:
  cvvs_f <- run_cv_wp(r$ref, refit_prj = FALSE)
  expect_identical(cvvs_f$loo_perturb_alpha, Inf)
  # cv_varsel.vsel() with `refit_prj = TRUE` reuses these searches and keeps
  # `Inf`:
  expect_no_message(
    cvvs_t <- suppressWarnings(cv_varsel(
      cvvs_f, refit_prj = TRUE, ndraws_pred = 40, seed = 11, verbose = FALSE
    )),
    message = "cannot be re-used"
  )
  expect_identical(cvvs_t$loo_perturb_alpha, Inf)
  expect_identical(ranking(cvvs_t)[["foldwise"]], rk_off)
  # An object without the element (older projpred versions) has searches
  # without the weights:
  cvvs_old <- cvvs_t
  cvvs_old$loo_perturb_alpha <- NULL
  expect_no_message(
    cvvs_o <- suppressWarnings(cv_varsel(
      cvvs_old, ndraws_pred = 40, seed = 11, verbose = FALSE
    )),
    message = "cannot be re-used"
  )
  expect_identical(cvvs_o$loo_perturb_alpha, Inf)
  # Without fold-wise searches to reuse, new searches use the default 10:
  cvvs_v <- suppressWarnings(cv_varsel(cvvs_f, validate_search = FALSE,
                                       verbose = FALSE))
  cvvs_n <- suppressWarnings(cv_varsel(
    cvvs_v, validate_search = TRUE, ndraws_pred = 40, seed = 11,
    verbose = FALSE
  ))
  expect_identical(cvvs_n$loo_perturb_alpha, 10)
  expect_identical(ranking(cvvs_n)[["foldwise"]],
                   ranking(run_cv_wp(r$ref))[["foldwise"]])
  # With `nloo < n`, cv_varsel.vsel() runs the fold-wise searches again and
  # keeps a finite value:
  cvvs_s <- run_cv_wp(r$ref, alpha = 4, nloo = 10)
  cvvs_s2 <- suppressMessages(suppressWarnings(cv_varsel(
    cvvs_s, ndraws_pred = 40, seed = 11, verbose = FALSE
  )))
  expect_identical(cvvs_s2$loo_perturb_alpha, 4)
})

test_that(paste(
  "weight perturbation gives the same results with and without the Gaussian",
  "cross-product search and evaluation"
), {
  skip_on_cran()
  r <- cre_ref_wp()
  res <- with_gauss_counts(run_cv_wp(r$ref))
  res_drawwise <- with_gauss_counts(withr::with_options(
    list(projpred.search_gauss_crossprod = FALSE),
    run_cv_wp(r$ref)
  ))
  # The default run uses the cross-product search in the full-data search and
  # in each fold, and the cross-product evaluation in each fold:
  expect_equal(res$search, r$n + 1)
  expect_equal(res$eval, r$n)
  expect_identical(res_drawwise$search, 0L)
  expect_identical(res_drawwise$eval, 0L)
  expect_identical(ranking(res$value)[["foldwise"]],
                   ranking(res_drawwise$value)[["foldwise"]])
  expect_equal(res$value$summaries$sub, res_drawwise$value$summaries$sub,
               tolerance = 1e-6)
})
