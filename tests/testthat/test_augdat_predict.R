context("augdat_predict")

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

if (exists(".Random.seed", envir = .GlobalEnv)) {
  rng_old <- get(".Random.seed", envir = .GlobalEnv)
}
set.seed(97531)
n_ap <- 120
dat_ap <- data.frame(a = rnorm(n_ap), b = rnorm(n_ap),
                     f = factor(sample(c("u", "v", "w"), n_ap, replace = TRUE)))
eta_ap <- with(dat_ap, 0.8 * a - 0.5 * b + 0.6 * (f == "v"))
dat_ap$yo <- cut(eta_ap + rlogis(n_ap), c(-Inf, -1, 0, 1, Inf),
                 labels = paste0("c", 1:4), ordered_result = TRUE)
dat_ap$yc <- factor(apply(cbind(0, eta_ap, -eta_ap, 0.3), 1, function(e) {
  sample(4, 1, prob = exp(e))
}))
newdat_ap <- data.frame(a = c(-1, 0.5, 2), b = c(0, 1, -1),
                        f = factor(c("w", "u", "v"), levels = c("u", "v", "w")))
rownames(newdat_ap) <- c("r1", "r2", "r3")

# subprd_augdat() without subprd_augdat_direct(), as on master:
lpreds_via_probs <- function(fits, newdata) {
  arrs <- lapply(fits, function(fit) {
    prbs <- predict(fit, newdata = newdata, type = "probs")
    if (inherits(fit, "polr")) {
      cprbs <- t(apply(prbs, 1, cumsum))[, -ncol(prbs), drop = FALSE]
      linkfun_raw(cprbs, link_nm = fit$method)
    } else {
      log(sweep(prbs[, -1, drop = FALSE], 1, prbs[, 1], "/"))
    }
  })
  do.call(abind::abind, list(arrs, rev.along = 0))
}

test_that("subprd_augdat_direct() gives the linear predictors of polr fits", {
  for (method in c("logistic", "probit", "cloglog")) {
    fits <- lapply(1:3, function(s) {
      wt <- runif(n_ap, 0.5, 1.5)
      suppressWarnings(MASS::polr(yo ~ a + b + f, data = cbind(dat_ap, wt = wt),
                                  weights = wt, method = method))
    })
    prd <- subprd_augdat_direct(fits, newdat_ap)
    expect_equal(prd, lpreds_via_probs(fits, newdat_ap), tolerance = 1e-10,
                 info = method)
    expect_identical(dimnames(prd), list(rownames(newdat_ap), c("c1", "c2", "c3"),
                                         NULL))
  }
  # Without predictors:
  fits_icpt <- lapply(1:2, function(s) MASS::polr(yo ~ 1, data = dat_ap))
  expect_equal(subprd_augdat_direct(fits_icpt, newdat_ap),
               lpreds_via_probs(fits_icpt, newdat_ap), tolerance = 1e-10)
})

test_that("subprd_augdat_direct() gives the linear predictors of multinom fits", {
  fits <- lapply(1:3, function(s) {
    wt <- runif(n_ap, 0.5, 1.5)
    nnet::multinom(yc ~ a + b + f, data = cbind(dat_ap, wt = wt), weights = wt,
                   trace = FALSE)
  })
  prd <- subprd_augdat_direct(fits, newdat_ap)
  expect_equal(prd, lpreds_via_probs(fits, newdat_ap), tolerance = 1e-10)
  expect_identical(dimnames(prd), list(rownames(newdat_ap), c("2", "3", "4"),
                                       NULL))
})

test_that("subprd_augdat_direct() returns `NULL` for other cases", {
  fit_p1 <- MASS::polr(yo ~ a, data = dat_ap)
  fit_p2 <- MASS::polr(yo ~ b, data = dat_ap)
  fit_m <- nnet::multinom(yc ~ a, data = dat_ap, trace = FALSE)
  expect_null(subprd_augdat_direct(list(fit_p1, fit_p2), newdat_ap))
  expect_null(subprd_augdat_direct(list(fit_p1, fit_m), newdat_ap))
  expect_null(subprd_augdat_direct(list(lm(a ~ b, data = dat_ap)), newdat_ap))
})

test_that("`ce_ptwise` of the augmented-data families sums over the categories", {
  fam <- extend_family(
    structure(list(family = "cumulative_rstanarm", link = "logit",
                   cats = paste0("c", 1:4)), class = "family"),
    augdat_y_unqs = paste0("c", 1:4), augdat_link = augdat_link_cumul,
    augdat_ilink = augdat_ilink_cumul, augdat_args_link = list(link = "logit"),
    augdat_args_ilink = list(link = "logit")
  )
  n_obs <- 6
  n_cat <- 4
  n_draws <- 3
  rand_prbs <- function() {
    arr <- array(runif(n_obs * n_cat * n_draws), c(n_obs, n_cat, n_draws))
    sweep(arr, c(1, 3), apply(arr, c(1, 3), sum), "/")
  }
  prbs_ref <- rand_prbs()
  prbs_sub <- rand_prbs()
  prbs_ref[1, 2, 1] <- 0
  prbs_sub[2, 3, 2] <- 0
  w_obs <- runif(n_obs)
  ce_expected <- matrix(NA_real_, n_obs, n_draws)
  for (i in seq_len(n_obs)) {
    for (s in seq_len(n_draws)) {
      pr <- prbs_ref[i, , s]
      ps <- prbs_sub[i, , s]
      pr[pr == 0] <- .Machine$double.eps
      ps[ps == 0] <- .Machine$double.eps
      ce_expected[i, s] <- w_obs[i] * (-sum(pr * log(ps)))
    }
  }
  expect_identical(
    fam$ce_ptwise(arr2augmat(prbs_ref), arr2augmat(prbs_sub), w_obs = w_obs),
    ce_expected
  )
})

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
