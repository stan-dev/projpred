context("cumul_native")

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

if (exists(".Random.seed", envir = .GlobalEnv)) {
  rng_old <- get(".Random.seed", envir = .GlobalEnv)
}
set.seed(8642)
n_cn <- 200
ncat_cn <- 4
x_cn <- matrix(rnorm(n_cn * 3), n_cn, 3, dimnames = list(NULL, c("a", "b", "c")))
eta_cn <- drop(x_cn %*% c(1, -0.5, 0.3))
links_cn <- c("logistic", "probit", "loglog", "cloglog", "cauchit")
pfun_cn <- list(logistic = plogis, probit = pnorm,
                loglog = function(q) exp(-exp(-q)),
                cloglog = function(q) 1 - exp(-exp(q)), cauchit = pcauchy)
# Category weights (n x C) as from the projected draws of a reference model:
wts_cn <- function(method, eta = eta_cn) {
  cp <- cbind(0, pfun_cn[[method]](outer(-eta, c(-1, 0, 1), "+")), 1)
  pmax(t(apply(cp, 1, diff)), 1e-10) * sample(1:3, length(eta), replace = TRUE)
}
# Augmented data (C copies of the observations) with the weights in column
# `wt`:
augdat_cn <- function(w) {
  data.frame(x_cn[rep(seq_len(n_cn), ncat_cn), ],
             y = factor(rep(seq_len(ncat_cn), each = n_cn)),
             wt = as.vector(w))
}

test_that("cumul_fit_multi_c() finds the optimum of the polr() objective", {
  for (method in setdiff(links_cn, "cloglog")) {
    w <- wts_cn(method)
    fit_polr <- suppressWarnings(MASS::polr(
      y ~ a + b + c, data = augdat_cn(w), weights = wt, method = method,
      start = c(1, -0.5, 0.3, -1, 0, 1),
      control = list(reltol = 1e-15, maxit = 5000)
    ))
    fit_nat <- cumul_fit_multi_c(x_cn, w, ncat_cn, rep(0, n_cn),
                                 match(method, links_cn), 1e-12, 100L)
    expect_identical(fit_nat$convergence, 0L, info = method)
    expect_equal(drop(fit_nat$beta), unname(coef(fit_polr)), tolerance = 1e-5,
                 info = method)
    expect_equal(drop(fit_nat$zeta), unname(fit_polr$zeta), tolerance = 1e-5,
                 info = method)
    expect_lt(2 * fit_nat$loss, fit_polr$deviance + 1e-6)
  }
  # cloglog (polr() often fails here because its probabilities of the highest
  # category round to zero): the gradient is zero at the optimum.
  w <- wts_cn("cloglog")
  fit_nat <- cumul_fit_multi_c(x_cn, w, ncat_cn, rep(0, n_cn), 4L, 1e-12, 100L)
  expect_identical(fit_nat$convergence, 0L)
  # Survival function of the cloglog latent variable (1 - pfun, without
  # cancellation):
  sfun <- function(q) exp(-exp(q))
  obj <- function(theta) {
    eta <- drop(x_cn %*% theta[1:3])
    g <- c(-Inf, theta[4:6], Inf)
    -sum(vapply(seq_len(ncat_cn), function(k) {
      z1 <- pmin(100, g[k + 1] - eta)
      z2 <- pmax(-100, g[k] - eta)
      pr <- ifelse(z2 > 0, sfun(z2) - sfun(z1),
                   pfun_cn$cloglog(z1) - pfun_cn$cloglog(z2))
      sum(w[, k] * log(pr))
    }, numeric(1)))
  }
  theta_hat <- c(fit_nat$beta, fit_nat$zeta)
  grad_num <- vapply(seq_along(theta_hat), function(j) {
    h <- 1e-6
    e <- replace(numeric(length(theta_hat)), j, h)
    (obj(theta_hat + e) - obj(theta_hat - e)) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(grad_num)), 1e-3)
  # Several weight matrices in one call give the same fits as one at a time:
  w2 <- wts_cn("logistic")
  fit_both <- cumul_fit_multi_c(x_cn, cbind(w, w2), ncat_cn, rep(0, n_cn), 1L,
                                1e-12, 100L)
  fit_2 <- cumul_fit_multi_c(x_cn, w2, ncat_cn, rep(0, n_cn), 1L, 1e-12, 100L)
  expect_identical(fit_both$beta[, 2], drop(fit_2$beta))
  expect_identical(fit_both$zeta[, 2], drop(fit_2$zeta))
})

test_that("cumul_native_multi() returns usable `polr` objects", {
  method <- "logistic"
  ws <- cbind(as.vector(wts_cn(method)), as.vector(wts_cn(method)))
  daug <- augdat_cn(ws[, 1])[c("a", "b", "c", "y")]
  fam <- structure(list(family = "cumulative", link = "logit"),
                   class = "family")
  fits <- cumul_native_multi(y ~ a + b, data = daug, family = fam,
                             ws_aug = ws, dot_args = list())
  expect_length(fits, 2)
  newdat <- data.frame(a = c(-1, 0, 2), b = c(0.5, 0, -1))
  for (s in 1:2) {
    daug$wt <- ws[, s]
    fit_polr <- suppressWarnings(MASS::polr(
      y ~ a + b, data = daug, weights = wt, method = method,
      control = list(reltol = 1e-15, maxit = 5000)
    ))
    expect_s3_class(fits[[s]], "polr")
    expect_identical(fits[[s]]$convergence, 0L)
    expect_identical(names(coef(fits[[s]])), names(coef(fit_polr)))
    expect_identical(names(fits[[s]]$zeta), names(fit_polr$zeta))
    expect_equal(coef(fits[[s]]), coef(fit_polr), tolerance = 1e-5)
    expect_equal(fits[[s]]$zeta, fit_polr$zeta, tolerance = 1e-5)
    expect_equal(predict(fits[[s]], newdata = newdat, type = "probs"),
                 predict(fit_polr, newdata = newdat, type = "probs"),
                 tolerance = 1e-5)
    expect_true(check_conv_s(fits[[s]]))
  }
  # Submodel without predictors: named empty coefficient vector, as from polr():
  fits_icpt <- cumul_native_multi(y ~ 1, data = daug, family = fam,
                                  ws_aug = ws, dot_args = list())
  expect_identical(coef(fits_icpt[[1]]), setNames(numeric(), character()))
  # Cases that use fit_cumul() per draw:
  expect_null(cumul_native_multi(y ~ a + b, data = daug, family = fam,
                                 ws_aug = ws, dot_args = list(Hess = TRUE)))
  opts_old <- options(projpred.augdat_native = FALSE)
  expect_null(cumul_native_multi(y ~ a + b, data = daug, family = fam,
                                 ws_aug = ws, dot_args = list()))
  options(opts_old)
  daug_perm <- daug[c(2:nrow(daug), 1), ]
  expect_null(cumul_native_multi(y ~ a + b, data = daug_perm, family = fam,
                                 ws_aug = ws, dot_args = list()))
})

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
