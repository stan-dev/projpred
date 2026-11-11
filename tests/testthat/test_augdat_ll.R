context("augmented-data log-likelihood from the linear predictor")

# Augmented-data families of the tests (`ilink` is the inverse link that the
# family gets in extend_family()):
softmax_ref1 <- function(eta_arr) {
  eta_full <- abind::abind(array(0, dim = c(dim(eta_arr)[-3], 1)), eta_arr)
  prbs <- exp(sweep(eta_full, c(1, 2), apply(eta_full, c(1, 2), max)))
  sweep(prbs, c(1, 2), apply(prbs, c(1, 2), sum), "/")
}
fam_augdat <- function(fam_nm, link, ncat, ilink = NULL, args_ilink = NULL) {
  if (fam_nm == "binomial") {
    ilink <- ilink %||% augdat_ilink_binom
    lnk <- augdat_link_binom
  } else if (fam_nm == "categorical") {
    ilink <- ilink %||% softmax_ref1
    lnk <- function(prb_arr) {
      log(sweep(prb_arr[, , -1, drop = FALSE], c(1, 2), prb_arr[, , 1], "/"))
    }
  } else {
    ilink <- ilink %||% augdat_ilink_cumul
    lnk <- augdat_link_cumul
  }
  if (is.null(args_ilink) && fam_nm != "categorical") {
    args_ilink <- list(link = link)
  }
  extend_family(
    structure(list(family = fam_nm, link = link), class = "family"),
    augdat_y_unqs = paste0("c", seq_len(ncat)), augdat_link = lnk,
    augdat_ilink = ilink,
    augdat_args_link = if (fam_nm != "categorical") list(link = link),
    augdat_args_ilink = args_ilink %||% list()
  )
}
# The augmented-rows matrix of the linear predictor from an N x C_lat x S array:
eta_augmat <- function(eta_arr) arr2augmat(eta_arr)
# The log-likelihood with the observed category `y` (a factor) from `eta_arr`:
ll_of <- function(fam, eta_arr, y, wobs = 1, use_eta = TRUE) {
  eta <- eta_augmat(eta_arr)
  fam$ll_fun(fam$linkinv(eta), NULL, y, wobs, eta = if (use_eta) eta)
}
log1pexp_tst <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))

test_that("augdat ll_fun() with `eta` equals ll_cats() for moderate eta", {
  withr::local_seed(4411)
  n_obs <- 7L
  n_draws <- 5L
  cases <- c(
    lapply(c("logit", "logistic", "probit", "probit_approx", "cauchit",
             "cloglog"), function(l) list(fam = "cumulative", link = l)),
    lapply(c("logit", "probit", "cloglog"), function(l) {
      list(fam = "cumulative_rstanarm", link = l)
    }),
    lapply(c("logit", "probit", "cauchit", "cloglog"), function(l) {
      list(fam = "binomial", link = l)
    }),
    list(list(fam = "categorical", link = "logit"))
  )
  for (cs in cases) {
    info <- paste(cs$fam, cs$link)
    ncat <- if (cs$fam == "binomial") 2 else 4
    fam <- fam_augdat(cs$fam, cs$link, ncat)
    lp <- matrix(rnorm(n_obs * n_draws, sd = 1.5), n_obs, n_draws)
    if (cs$fam %in% c("cumulative", "cumulative_rstanarm")) {
      thr <- c(-1.2, 0.1, 1.3)
      eta_arr <- aperm(outer(thr, lp, function(t, l) t - l), c(2, 1, 3))
    } else {
      eta_arr <- array(rnorm(n_obs * (ncat - 1) * n_draws, sd = 2),
                       c(n_obs, ncat - 1, n_draws))
    }
    y <- factor(sample(paste0("c", seq_len(ncat)), n_obs, replace = TRUE),
                levels = paste0("c", seq_len(ncat)))
    ll_eta <- ll_of(fam, eta_arr, y)
    ll_mu <- ll_of(fam, eta_arr, y, use_eta = FALSE)
    expect_identical(dim(ll_eta), c(n_obs, n_draws), info = info)
    # The probabilities lose digits for small p (cloglog: already for an `eta`
    # above 3), so compare where p > exp(-10):
    sel <- ll_mu > -10
    expect_gt(mean(sel), 0.5, label = info)
    expect_equal(ll_eta[sel], ll_mu[sel], tolerance = 1e-10, info = info)
    # Observation weights, with a zero weight:
    wobs <- c(0, runif(n_obs - 1, 0.5, 2))
    expect_equal(ll_of(fam, eta_arr, y, wobs), wobs * ll_eta,
                 tolerance = 1e-14, info = info)
    expect_identical(ll_of(fam, eta_arr, y, wobs)[1, ], rep(0, n_draws),
                     info = info)
  }
})

test_that("augdat ll_fun() with `eta` is exact for small probabilities", {
  # One observation and one draw per call; `eta` holds the latent predictors of
  # the categories.
  ll_1 <- function(fam, eta_vec, cat, use_eta = TRUE) {
    eta_arr <- array(eta_vec, c(1, length(eta_vec), 1))
    y <- factor(cat, levels = fam$cats)
    ll_of(fam, eta_arr, y, use_eta = use_eta)[1, 1]
  }
  ll_bound <- log(.Machine$double.eps)
  # Cumulative logit: log(F(b) - F(a)) = b + log(1 - exp(a - b)) -
  # log(1 + exp(a)) - log(1 + exp(b)).
  fam <- fam_augdat("cumulative", "logit", 4)
  lp_diff_logit <- function(a, b) {
    b + log1p(-exp(a - b)) - log1pexp_tst(a) - log1pexp_tst(b)
  }
  for (shift in c(40, 60, 300)) {
    thr <- c(-1, 0, 1.5) + shift
    expect_equal(ll_1(fam, thr, "c2"), lp_diff_logit(thr[1], thr[2]),
                 tolerance = 1e-12, info = paste("shift:", shift))
    expect_equal(ll_1(fam, thr, "c4"), -log1pexp_tst(thr[3]),
                 tolerance = 1e-12, info = paste("shift:", shift))
    expect_equal(ll_1(fam, -rev(thr), "c1"), -log1pexp_tst(thr[3]),
                 tolerance = 1e-12, info = paste("shift:", shift))
    expect_equal(ll_1(fam, -rev(thr), "c3"),
                 lp_diff_logit(-thr[2], -thr[1]), tolerance = 1e-12,
                 info = paste("shift:", shift))
    # Without `eta`, the lower bound:
    expect_equal(ll_1(fam, thr, "c2", use_eta = FALSE), ll_bound)
  }
  # Cancellation: both cumulative probabilities close to 1 (p about 4e-14),
  # above the lower bound:
  thr <- c(30, 30.5, 35)
  expect_equal(ll_1(fam, thr, "c2"), lp_diff_logit(thr[1], thr[2]),
               tolerance = 1e-12)
  expect_gt(abs(ll_1(fam, thr, "c2", use_eta = FALSE) -
                  lp_diff_logit(thr[1], thr[2])), 1e-6)
  # Probit: with thresholds far apart, p is F(b) (or 1 - F(a)) up to a
  # relative error below 1e-100:
  fam <- fam_augdat("cumulative", "probit", 3)
  expect_equal(ll_1(fam, c(-50, -40), "c2"), pnorm(-40, log.p = TRUE),
               tolerance = 1e-12)
  expect_equal(ll_1(fam, c(40, 50), "c2"), pnorm(-40, log.p = TRUE),
               tolerance = 1e-12)
  expect_equal(ll_1(fam, c(-50, -40), "c3"),
               pnorm(-40, lower.tail = FALSE, log.p = TRUE))
  # Cloglog: F(x) = 1 - exp(-exp(x)) = -expm1(-exp(x)) (no cancellation in
  # the lower tail):
  fam <- fam_augdat("cumulative", "cloglog", 3)
  expect_equal(ll_1(fam, c(-60, -40), "c2"),
               log(-expm1(-exp(-40)) + expm1(-exp(-60))), tolerance = 1e-12)
  expect_equal(ll_1(fam, c(3, 5), "c2"),
               -exp(3) + log1p(-exp(exp(3) - exp(5))), tolerance = 1e-12)
  expect_equal(ll_1(fam, c(3, 5), "c3"), -exp(5), tolerance = 1e-12)
  # Binomial (augmented-data): the categories are 1 - F(eta) and F(eta).
  fam <- fam_augdat("binomial", "logit", 2)
  expect_equal(ll_1(fam, 50, "c1"), -log1pexp_tst(50), tolerance = 1e-12)
  expect_equal(ll_1(fam, -50, "c2"), -log1pexp_tst(50), tolerance = 1e-12)
  # Categorical: log-softmax with the first category as the reference.
  fam <- fam_augdat("categorical", "logit", 3)
  lse <- function(x) max(x) + log(sum(exp(x - max(x))))
  expect_equal(ll_1(fam, c(800, -3), "c1"), -lse(c(0, 800, -3)),
               tolerance = 1e-12)
  expect_equal(ll_1(fam, c(800, -3), "c3"), -3 - lse(c(0, 800, -3)),
               tolerance = 1e-12)
  expect_equal(ll_1(fam, c(800, -3), "c1", use_eta = FALSE), ll_bound)
})

test_that("augdat ll_fun() keeps ll_cats() for `NA` and other cases", {
  withr::local_seed(2207)
  fam <- fam_augdat("cumulative", "logit", 3)
  eta_arr <- array(c(-0.5, 0.2, 1, 1.4, 2, 2.1), c(2, 2, 3))
  eta_arr[, 2, ] <- eta_arr[, 1, ] + 1
  eta_arr[2, 1, 3] <- -60
  eta_arr[2, 2, 3] <- -59
  y <- factor(c("c2", "c3"), levels = fam$cats)
  eta_na <- eta_arr
  eta_na[1, 2, 2] <- NA
  ll_na <- ll_of(fam, eta_na, y)
  ll_mu <- ll_of(fam, eta_na, y, use_eta = FALSE)
  expect_identical(ll_na[1, 2], ll_mu[1, 2])
  expect_true(is.na(ll_mu[1, 2]))
  expect_equal(ll_na[-3], ll_of(fam, eta_arr, y)[-3], tolerance = 1e-14)
  # No formula (family `sratio`): ll_cats().
  fam_s <- extend_family(
    structure(list(family = "sratio", link = "logit"), class = "family"),
    augdat_y_unqs = c("c1", "c2", "c3"), augdat_link = augdat_link_cumul,
    augdat_ilink = augdat_ilink_cumul, augdat_args_link = list(link = "logit"),
    augdat_args_ilink = list(link = "logit")
  )
  expect_identical(ll_of(fam_s, eta_arr, y),
                   ll_of(fam_s, eta_arr, y, use_eta = FALSE))
  # An `augdat_ilink` that differs from the standard one for the family and the
  # link (here: probit instead of logit): a warning and ll_cats().
  fam_w <- fam_augdat("cumulative", "logit", 3,
                      args_ilink = list(link = "probit"))
  expect_warning(ll_w <- ll_of(fam_w, eta_arr, y),
                 "differ from those of the standard inverse link")
  expect_identical(ll_w, ll_of(fam_w, eta_arr, y, use_eta = FALSE))
  # Wrong dimensions of `eta`:
  eta <- eta_augmat(eta_arr)
  expect_error(fam$ll_fun(fam$linkinv(eta), NULL, y, 1, eta = eta[, 1:2]),
               "Unexpected dimensions of `eta`")
})

# A cumulative reference model (augmented-data projection, no Stan fit) with
# draws that have a large linear predictor:
refm_augdat_extreme <- function(link = "logit", offset = FALSE, n = 40,
                                ndraws = 60, n_extreme = 10, seed = 3381) {
  withr::local_seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dat$offs <- if (offset) runif(n, -1, 1) else 0
  thr <- c(-1, 0.3, 1.2)
  lp_true <- 1.2 * dat$x1 - 0.8 * dat$x2 + dat$offs
  dat$y <- cut(lp_true + rlogis(n), c(-Inf, thr, Inf),
               labels = paste0("c", 1:4), ordered_result = TRUE)
  beta <- rbind(rnorm(ndraws, 1.2, 0.2), rnorm(ndraws, -0.8, 0.2))
  beta[1, seq_len(n_extreme)] <- sample(c(-1, 1), n_extreme, replace = TRUE) *
    runif(n_extreme, 50, 100)
  thr_draws <- t(sapply(thr, function(t) t + rnorm(ndraws, 0, 0.05)))
  fml <- y ~ x1 + x2
  # The latent predictors `thr_k - lp` (S x N x C_lat), offsets included:
  eta_of <- function(fit, nd) {
    lp <- t(as.matrix(nd[, c("x1", "x2")]) %*% fit$beta) +
      matrix(nd$offs, nrow = ncol(fit$beta), ncol = nrow(nd), byrow = TRUE)
    simplify2array(lapply(seq_len(nrow(fit$thr)), function(k) {
      fit$thr[k, ] - lp
    }))
  }
  ref_predfun <- function(fit, newdata = NULL) {
    eta_of(fit, if (is.null(newdata)) fit$data else newdata)
  }
  extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                 extract_y = TRUE) {
    y_wobs_offs(newdata = newdata, wrhs = wrhs,
                orhs = orhs %||% if (offset) ~ offs else NULL,
                resp_form = if (extract_y) ~ y else NULL)
  }
  fit <- list(beta = beta, thr = thr_draws, data = dat)
  refm <- withCallingHandlers(
    init_refmodel(fit, data = dat, formula = fml,
                  family = structure(list(family = "cumulative_rstanarm",
                                          link = link), class = "family"),
                  ref_predfun = ref_predfun,
                  extract_model_data = extract_model_data,
                  augdat_y_unqs = paste0("c", 1:4),
                  augdat_link = augdat_link_cumul,
                  augdat_ilink = augdat_ilink_cumul,
                  augdat_args_link = list(link = link),
                  augdat_args_ilink = list(link = link)),
    warning = function(w) {
      if (grepl("`cvrefbuilder` was `NULL`", conditionMessage(w),
                fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  # The exact N x S log-likelihood (logit link), from the draws (not from
  # `refm`), with the thresholds -Inf and Inf for the first and last category:
  stopifnot(link == "logit")
  eta_arr <- eta_of(fit, dat)
  y_int <- as.integer(dat$y)
  ll_exact <- t(sapply(seq_len(n), function(i) {
    e <- cbind(-Inf, eta_arr[, i, ], Inf)
    a <- e[, y_int[i]]
    b <- e[, y_int[i] + 1]
    if (y_int[i] == 1) {
      -log1pexp_tst(-b)
    } else if (y_int[i] == 4) {
      -log1pexp_tst(a)
    } else {
      b + log1p(-exp(a - b)) - log1pexp_tst(a) - log1pexp_tst(b)
    }
  }))
  list(refm = refm, ll_exact = ll_exact)
}

test_that("augmented-data reference model LOO and predict() are exact", {
  skip_if_not_installed("MASS")
  for (offset in c(FALSE, TRUE)) {
    info <- paste("offset:", offset)
    # predict.refmodel() with `ynew`, all draws with a large |eta|:
    obj_all <- refm_augdat_extreme(offset = offset, n_extreme = 60)
    refm_all <- obj_all$refm
    expect_identical(any(refm_all$offset != 0), offset, info = info)
    lpd <- predict(refm_all, newdata = refm_all$fetch_data(),
                   ynew = refm_all$y)
    lpd_exact <- apply(obj_all$ll_exact, 1, log_sum_exp) -
      log(ncol(obj_all$ll_exact))
    expect_lt(min(lpd_exact), log(.Machine$double.eps))
    expect_equal(lpd, lpd_exact, tolerance = 1e-10, info = info)
    # The PSIS-LOO values of the reference model in cv_varsel():
    obj <- refm_augdat_extreme(offset = offset)
    expect_lt(min(obj$ll_exact), log(.Machine$double.eps))
    cvvs <- suppressWarnings(cv_varsel(
      obj$refm, cv_method = "LOO", validate_search = FALSE, nterms_max = 1,
      method = "forward", verbose = FALSE
    ))
    psis_exact <- suppressWarnings(loo::psis(-t(obj$ll_exact), cores = 1,
                                             r_eff = NA))
    lppd_exact <- apply(t(obj$ll_exact) + weights(psis_exact), 2,
                        log_sum_exp)
    expect_equal(cvvs$summaries$ref$lppd, lppd_exact, tolerance = 1e-10,
                 info = info)
  }
})
