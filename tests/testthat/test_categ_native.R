context("categ_native")

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

if (exists(".Random.seed", envir = .GlobalEnv)) {
  rng_old <- get(".Random.seed", envir = .GlobalEnv)
}
set.seed(24680)
n_ca <- 150
ncat_ca <- 4
dat_ca <- data.frame(a = rnorm(n_ca), b = rnorm(n_ca),
                     f = factor(sample(c("u", "v"), n_ca, replace = TRUE)))
eta_ca <- with(dat_ca, cbind(0, 0.8 * a, -0.6 * b + 0.5 * (f == "v"),
                             0.4 * a + 0.3 * b))
# Category weights (n x C) as from the projected draws of a reference model:
wts_ca <- function() {
  prb <- exp(eta_ca) / rowSums(exp(eta_ca))
  prb * sample(1:3, n_ca, replace = TRUE)
}
# Augmented data (C copies of the observations) with the weights in column
# `wt`:
augdat_ca <- function(w) {
  cbind(dat_ca[rep(seq_len(n_ca), ncat_ca), ],
        y = factor(rep(seq_len(ncat_ca), each = n_ca)),
        wt = as.vector(w))
}
fit_tight <- function(fml, daug) {
  nnet::multinom(fml, data = daug, weights = wt, trace = FALSE,
                 maxit = 10000, reltol = 1e-15, abstol = 0)
}

test_that("categ_fit_multi_c() finds the optimum of the multinom() objective", {
  w <- wts_ca()
  daug <- augdat_ca(w)
  fit_mn <- fit_tight(y ~ a + b + f, daug)
  x <- model.matrix(~ a + b + f, data = dat_ca)
  fit_nat <- categ_fit_multi_c(x, w, ncat_ca, 0L, 1e-12, 100L)
  expect_identical(fit_nat$convergence, 0L)
  expect_equal(t(fit_nat$coef[, , 1]), unname(coef(fit_mn)), tolerance = 1e-5)
  expect_lt(fit_nat$loss, fit_mn$value + 1e-8)
  # Several weight matrices in one call give the same fits as one at a time:
  w2 <- wts_ca()
  fit_both <- categ_fit_multi_c(x, cbind(w, w2), ncat_ca, 0L, 1e-12, 100L)
  fit_2 <- categ_fit_multi_c(x, w2, ncat_ca, 0L, 1e-12, 100L)
  expect_identical(fit_both$coef[, , 2], fit_2$coef[, , 1])
})

test_that("categ_native_multi() returns usable `multinom` objects", {
  ws <- cbind(as.vector(wts_ca()), as.vector(wts_ca()))
  daug <- augdat_ca(ws[, 1])[c("a", "b", "f", "y")]
  fam <- structure(list(family = "categorical", link = "logit"),
                   class = "family")
  fits <- categ_native_multi(y ~ a + f, data = daug, family = fam,
                             ws_aug = ws, dot_args = list())
  expect_length(fits, 2)
  newdat <- data.frame(a = c(-1, 0.5, 2),
                       f = factor(c("v", "u", "v"), levels = c("u", "v")))
  for (s in 1:2) {
    daug$wt <- ws[, s]
    fit_mn <- fit_tight(y ~ a + f, daug)
    expect_s3_class(fits[[s]], c("multinom", "nnet"))
    expect_identical(fits[[s]]$convergence, 0L)
    expect_identical(dimnames(coef(fits[[s]])), dimnames(coef(fit_mn)))
    expect_identical(fits[[s]]$lev, fit_mn$lev)
    expect_equal(coef(fits[[s]]), coef(fit_mn), tolerance = 1e-5)
    expect_equal(fits[[s]]$value, fit_mn$value, tolerance = 1e-8)
    expect_equal(predict(fits[[s]], newdata = newdat, type = "probs"),
                 predict(fit_mn, newdata = newdat, type = "probs"),
                 tolerance = 1e-5)
    expect_equal(unname(fits[[s]]$fitted.values), unname(fit_mn$fitted.values),
                 tolerance = 1e-5)
    expect_true(check_conv_s(fits[[s]]))
  }
  # Submodel without predictors:
  fits_icpt <- categ_native_multi(y ~ 1, data = daug, family = fam,
                                  ws_aug = ws, dot_args = list())
  daug$wt <- ws[, 1]
  expect_equal(coef(fits_icpt[[1]]), coef(fit_tight(y ~ 1, daug)),
               tolerance = 1e-5)
  # Cases that use fit_categ() per draw:
  expect_null(categ_native_multi(y ~ a, data = daug, family = fam,
                                 ws_aug = ws, dot_args = list(Hess = TRUE)))
  opts_old <- options(projpred.augdat_native = FALSE)
  expect_null(categ_native_multi(y ~ a, data = daug, family = fam,
                                 ws_aug = ws, dot_args = list()))
  options(opts_old)
  expect_null(categ_native_multi(y ~ a, data = daug[c(2:nrow(daug), 1), ],
                                 family = fam, ws_aug = ws, dot_args = list()))
})

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
