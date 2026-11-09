context("log-likelihood where the inverse link clamps the mean")

# The exact log-likelihood of the binomial (logit) and Poisson (log) models.
# The product with a zero count is 0, also for an infinite `eta`.
xlogp_tst <- function(k, logp) ifelse(k == 0, 0, k * logp)
ll_binom_exact <- function(eta, y, trials) {
  lchoose(trials, y) + xlogp_tst(y, plogis(eta, log.p = TRUE)) +
    xlogp_tst(trials - y, plogis(-eta, log.p = TRUE))
}
ll_poiss_exact <- function(eta, y, wobs) {
  # ifelse() keeps the dimensions of its first argument:
  wobs <- array(wobs, dim(eta))
  xlogp_tst(wobs, xlogp_tst(y, eta) - exp(eta) - lgamma(y + 1))
}

test_that("binomial ll_fun() with `eta` is exact where linkinv() clamps", {
  fam <- extend_family(binomial())
  eta <- cbind(c(-450, -40, -31, -30), c(450, 40, 31, 30), c(-5, 0, 5, 20),
               c(-Inf, Inf, -Inf, Inf))
  # Bernoulli, 4 trials, and a vector of trials with a zero:
  cases <- list(list(trials = 1, y = c(1, 0, 1, 0)),
                list(trials = 4, y = c(4, 0, 2, 3)),
                list(trials = c(1, 4, 2, 0), y = c(0, 4, 1, 0)))
  for (cs in cases) {
    trials <- cs$trials
    y <- cs$y
    y_mat <- matrix(y, nrow(eta), ncol(eta))
    t_mat <- matrix(trials, nrow(eta), ncol(eta))
    ll <- fam$ll_fun(fam$linkinv(eta), NULL, y, trials, eta = eta)
    info <- paste("trials:", paste(trials, collapse = ","))
    expect_equal(ll, ll_binom_exact(eta, y_mat, t_mat), tolerance = 1e-12,
                 info = info)
    # For |eta| <= 10, the values from `mu` are accurate and agree:
    ll_mu <- fam$ll_fun(fam$linkinv(eta), NULL, y, trials)
    is_small <- abs(eta) <= 10
    expect_equal(ll[is_small], ll_mu[is_small], tolerance = 1e-12,
                 info = info)
    # Without `eta`, the clamp gives values far from the exact ones (this is
    # the case that `eta` fixes):
    is_clamped <- abs(eta) > 30
    expect_gt(max(abs(ll_mu - ll_binom_exact(eta, y_mat, t_mat))[is_clamped]),
              100, label = info)
  }
  # Another link does not use `eta`:
  fam_probit <- extend_family(binomial(link = "probit"))
  eta_p <- cbind(c(-10, 10), c(-1, 1))
  expect_identical(
    fam_probit$ll_fun(fam_probit$linkinv(eta_p), NULL, c(1, 0), 1,
                      eta = eta_p),
    fam_probit$ll_fun(fam_probit$linkinv(eta_p), NULL, c(1, 0), 1)
  )
})

test_that("Poisson ll_fun() with `eta` is exact where linkinv() clamps", {
  fam <- extend_family(poisson())
  eta <- cbind(c(-100, -40, -37, -36), c(-Inf, -Inf, 0, 2))
  y <- c(0, 3, 1, 2)
  y_mat <- matrix(y, nrow(eta), ncol(eta))
  # The zero weight meets y = 3 and eta = -Inf (row 2, column 2):
  for (wobs in list(1, c(1, 2, 0.5, 1), c(1, 0, 0.5, 2))) {
    w_mat <- matrix(wobs, nrow(eta), ncol(eta))
    info <- paste("wobs:", paste(wobs, collapse = ","))
    ll <- fam$ll_fun(fam$linkinv(eta), NULL, y, wobs, eta = eta)
    expect_equal(ll, ll_poiss_exact(eta, y_mat, w_mat), tolerance = 1e-12,
                 info = info)
    ll_mu <- fam$ll_fun(fam$linkinv(eta), NULL, y, wobs)
    is_unclamped <- exp(eta) >= .Machine$double.eps
    expect_identical(ll[is_unclamped], ll_mu[is_unclamped], info = info)
    expect_gt(max(abs(ll_mu - ll_poiss_exact(eta, y_mat, w_mat))[!is_unclamped]),
              0.1, label = info)
  }
})

# A Bernoulli or Poisson reference model with simulated draws whose linear
# predictor reaches |eta| of about 100 to 300 for some draws (no model fit is
# necessary; the draws are not a posterior).
refm_extreme <- function(family, n = 40, ndraws = 60, seed = 5523,
                         offset = FALSE, n_extreme = 10) {
  withr::local_seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dat$offs <- if (offset) runif(n, -1, 1) else 0
  eta_true <- 1.2 * dat$x1 - 0.8 * dat$x2 + dat$offs
  dat$y <- if (family$family == "binomial") {
    rbinom(n, 1, plogis(eta_true))
  } else {
    rpois(n, exp(eta_true))
  }
  # `n_extreme` of the draws have a coefficient of x1 between 50 and 100:
  beta <- rbind(rnorm(ndraws, 0, 0.1), rnorm(ndraws, 1.2, 0.2),
                rnorm(ndraws, -0.8, 0.2))
  beta[2, seq_len(n_extreme)] <- sample(c(-1, 1), n_extreme, replace = TRUE) *
    runif(n_extreme, 50, 100)
  fml <- y ~ x1 + x2
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(delete.response(terms(fml)), nd) %*% fit$beta
  }
  extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                 extract_y = TRUE) {
    y_wobs_offs(newdata = newdata, wrhs = wrhs,
                orhs = orhs %||% if (offset) ~ offs else NULL,
                resp_form = if (extract_y) ~ y else NULL)
  }
  # The tests do not use K-fold CV, so the default `cvrefbuilder` is fine:
  withCallingHandlers(
    init_refmodel(list(beta = beta, data = dat), data = dat, formula = fml,
                  family = family, ref_predfun = ref_predfun,
                  extract_model_data = extract_model_data),
    warning = function(w) {
      if (grepl("`cvrefbuilder` was `NULL`", conditionMessage(w),
                fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

test_that("reference model LOO and predict() are exact for a large |eta|", {
  for (fam in list(binomial(), poisson())) for (offset in c(FALSE, TRUE)) {
    info <- paste(fam$family, "offset:", offset)
    # The exact log-likelihood (N x S) of the draws of `refm`:
    ll_exact_of <- function(refm) {
      eta <- refm$eta + refm$offset
      y_mat <- matrix(refm$y, nrow(eta), ncol(eta))
      if (fam$family == "binomial") {
        ll_binom_exact(eta, y_mat, 1)
      } else {
        ll_poiss_exact(eta, y_mat, 1)
      }
    }
    # predict.refmodel() with `ynew`, all draws with a large |eta| (with only
    # some such draws, the other draws dominate the mean over the draws):
    refm_all <- refm_extreme(fam, offset = offset, n_extreme = 60)
    expect_identical(any(refm_all$offset != 0), offset, info = info)
    expect_gt(min(apply(abs(refm_all$eta), 2, max)), 30)
    lpd <- predict(refm_all, newdata = refm_all$fetch_data(),
                   ynew = refm_all$y)
    lpd_exact <- apply(ll_exact_of(refm_all), 1, log_sum_exp) -
      log(ncol(refm_all$eta))
    expect_equal(lpd, lpd_exact, tolerance = 1e-10, info = info)
    # The PSIS-LOO values of the reference model in cv_varsel() (10 of the
    # draws with a large |eta|):
    refm <- refm_extreme(fam, offset = offset)
    ll_exact <- ll_exact_of(refm)
    expect_gt(max(abs(refm$eta)), 100)
    cvvs <- suppressWarnings(cv_varsel(
      refm, cv_method = "LOO", validate_search = FALSE, nterms_max = 1,
      method = "forward", verbose = FALSE
    ))
    psis_exact <- suppressWarnings(loo::psis(-t(ll_exact), cores = 1,
                                             r_eff = NA))
    lppd_exact <- apply(t(ll_exact) + weights(psis_exact), 2, log_sum_exp)
    expect_equal(cvvs$summaries$ref$lppd, lppd_exact, tolerance = 1e-10,
                 info = info)
  }
})
