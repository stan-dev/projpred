withr::local_seed(8642)
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

test_that("cumul_fit_multi_c() returns the polr() objective, also at the bounds", {
  # Survival functions of the latent variables (1 - pfun, without
  # cancellation):
  sfun_cn <- list(logistic = function(q) plogis(q, lower.tail = FALSE),
                  probit = function(q) pnorm(q, lower.tail = FALSE),
                  loglog = function(q) -expm1(-exp(-q)),
                  cloglog = function(q) exp(-exp(q)),
                  cauchit = function(q) pcauchy(q, lower.tail = FALSE))
  # Offsets of +-150 for some observations, so that the bounds +-100 apply.
  # These observations have weight only in the highest or lowest category
  # (else the objective is infinite). No random numbers: later tests use the
  # random number state.
  offs <- rep_len(c(0, 0, 150, -150, 0.5), n_cn)
  for (method in links_cn) {
    cp <- cbind(0, pfun_cn[[method]](outer(-eta_cn, c(-1, 0, 1), "+")), 1)
    w <- pmax(t(apply(cp, 1, diff)), 1e-10)
    w[offs == 150, ] <- rep(c(0, 0, 0, 1), each = sum(offs == 150))
    w[offs == -150, ] <- rep(c(1, 0, 0, 0), each = sum(offs == -150))
    fit_nat <- cumul_fit_multi_c(x_cn, w, ncat_cn, offs,
                                 match(method, links_cn), 1e-12, 100L)
    eta <- offs + drop(x_cn %*% fit_nat$beta)
    g <- c(-Inf, fit_nat$zeta, Inf)
    loss <- -sum(vapply(seq_len(ncat_cn), function(k) {
      z1 <- pmin(100, g[k + 1] - eta)
      z2 <- pmax(-100, g[k] - eta)
      pr <- ifelse(z2 > 0, sfun_cn[[method]](z2) - sfun_cn[[method]](z1),
                   pfun_cn[[method]](z1) - pfun_cn[[method]](z2))
      # As in the C++ code, terms with weight 0 are left out:
      pos <- w[, k] > 0
      sum(w[pos, k] * log(pr[pos]))
    }, numeric(1)))
    expect_equal(drop(fit_nat$loss), loss, tolerance = 1e-10, info = method)
  }
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
  expect_null(withr::with_options(
    list(projpred.augdat_native = FALSE),
    cumul_native_multi(y ~ a + b, data = daug, family = fam,
                       ws_aug = ws, dot_args = list())
  ))
  daug_perm <- daug[c(2:nrow(daug), 1), ]
  expect_null(cumul_native_multi(y ~ a + b, data = daug_perm, family = fam,
                                 ws_aug = ws, dot_args = list()))
})

test_that("cumul_fit_multi_c() starts within the bounds for cauchit", {
  # One observation in the lowest category, with weight 0.4: its weight share
  # of 0.2 % gives a cauchit quantile below -100. No random numbers here.
  u <- eta_cn + 2 * sin(seq_len(n_cn) * 2.7)
  y <- 2L + as.integer(u > 0) + as.integer(u > 1.5)
  y[1] <- 1L
  wt <- replace(rep(1, n_cn), 1, 0.4)
  w <- matrix(0, n_cn, ncat_cn)
  w[cbind(seq_len(n_cn), y)] <- wt
  fit_nat <- cumul_fit_multi_c(x_cn, w, ncat_cn, rep(0, n_cn), 5L, 1e-12,
                               100L)
  expect_identical(fit_nat$convergence, 0L)
  expect_true(is.finite(drop(fit_nat$loss)))
  fit_polr <- suppressWarnings(MASS::polr(factor(y) ~ x_cn, weights = wt,
                                          method = "cauchit"))
  expect_lt(2 * drop(fit_nat$loss), fit_polr$deviance + 1e-6)
})

test_that("the native cumulative-link fits handle non-finite input", {
  # No random numbers here.
  cp <- cbind(0, plogis(outer(-eta_cn, c(-1, 0, 1), "+")), 1)
  w <- t(apply(cp, 1, diff))
  # An observation with weight in a middle category and both thresholds
  # beyond the bound: the objective is infinite at the start values.
  offs <- replace(rep(0, n_cn), 1, 150)
  fit_inf <- cumul_fit_multi_c(x_cn, w, ncat_cn, offs, 1L, 1e-12, 100L)
  expect_identical(fit_inf$convergence, 1L)
  expect_identical(drop(fit_inf$loss), Inf)
  # A non-finite predictor value (a NaN linear predictor at the start values):
  x_inf <- x_cn
  x_inf[1, 1] <- Inf
  fit_nan <- cumul_fit_multi_c(x_inf, w, ncat_cn, rep(0, n_cn), 1L, 1e-12,
                               100L)
  expect_identical(fit_nan$convergence, 1L)
  expect_identical(drop(fit_nan$loss), Inf)
  fit_nan_offs <- cumul_fit_multi_c(x_cn, w, ncat_cn,
                                    replace(rep(0, n_cn), 2, NaN), 1L, 1e-12,
                                    100L)
  expect_identical(fit_nan_offs$convergence, 1L)
  expect_identical(drop(fit_nan_offs$loss), Inf)
  # cumul_native_multi() leaves such data to polr():
  daug <- augdat_cn(w)[c("a", "b", "c", "y")]
  daug$a[1 + (seq_len(ncat_cn) - 1) * n_cn] <- Inf
  fam <- structure(list(family = "cumulative", link = "logit"),
                   class = "family")
  expect_null(cumul_native_multi(y ~ a + b, data = daug, family = fam,
                                 ws_aug = as.vector(w), dot_args = list()))
  daug$a <- augdat_cn(w)$a
  expect_null(cumul_native_multi(y ~ a + b, data = daug, family = fam,
                                 ws_aug = replace(as.vector(w), 3, Inf),
                                 dot_args = list()))
})
