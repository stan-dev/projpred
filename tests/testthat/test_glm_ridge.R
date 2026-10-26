context("ridge")

# tests for ridge regression, currently untested combinations
# - gaussian with inverse-link
# - binomial with log or cloglog-link
# - poisson with sqrt or id-link
# - Gamma with inverse or id-link
# - everything except gaussian with id-link for ridge penalty

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

if (exists(".Random.seed", envir = .GlobalEnv)) {
  rng_old <- get(".Random.seed", envir = .GlobalEnv)
}
set.seed(1235)
n <- 40
nterms <- 10
nterms_fit <- nterms - 5
x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
b <- runif(nterms) - 0.5
dis <- runif(1, 1, 2)
x_tr <- x[, 1:nterms_fit]
weights <- sample(1:4, n, replace = TRUE)
offset <- rnorm(n, 0, 1)

tol <- 3e-03
# some link-functions seem to need higher thresh-argument for glm_ridge
# (gaussian-log, binomial-cauchit, Gamma-log)
extra_thresh <- 1e-10

test_that(paste(
  "glmfun: gradients should give the same results as finite",
  "differences"
), {
  fdiffu <- function(f, x, h = 1e-3, order = 1) {
    ## function for computing derivative of univariate function f at x using
    ## finite difference.
    if (order != 1 && order != 2) {
      stop("Order must be either 1 or 2.")
    }
    n <- length(x)
    df <- rep(0, n)
    for (i in 1:n) {
      if (order == 1) {
        df[i] <- (f(x[i] + h) - f(x[i])) / h
      } else if (order == 2) {
        df[i] <- (f(x[i] + h) - 2 * f(x[i]) + f(x[i] - h)) / h^2
      }
    }
    return(df)
  }

  fams <- list(
    extend_family(gaussian(link = "identity")),
    extend_family(gaussian(link = "log")),
    extend_family(binomial(link = "logit")),
    extend_family(binomial(link = "probit")),
    extend_family(binomial(link = "cauchit")),
    extend_family(poisson(link = "log")),
    extend_family(Student_t(nu = 3, link = "identity")),
    extend_family(Student_t(nu = 4, link = "log")),
    extend_family(Student_t(nu = 7, link = "inverse"))
  )

  n <- 10
  weights <- sample(1:4, n, replace = TRUE)
  offset <- rnorm(n, 0, 1)

  for (i in seq_along(fams)) {
    fam <- fams[[i]]
    if (fam$family == "gaussian" || fam$family == "Student_t") {
      y <- rnorm(n)
    } else if (fam$family == "binomial") {
      y <- rbinom(n, 1, 0.6)
    } else if (fam$family == "poisson") {
      y <- rpois(n, 1)
    }

    devfun <- function(f)
      pseudo_data(f, y, fam, weights = weights, offset = offset)$loss
    zfun <- function(f)
      pseudo_data(f, y, fam, weights = weights, offset = offset)$z
    wfun <- function(f)
      pseudo_data(f, y, fam, weights = weights, offset = offset)$w
    gradan <- function(f)
      sum(pseudo_data(f, y, fam, weights = weights, offset = offset)$grad)
    gradfd <- function(f)
      fdiffu(devfun, f, h = 1e-5) # finite difference

    # compare analytic and finite difference gradients
    fval <- seq(-5, 5, len = 100)
    gan <- sapply(fval, gradan)
    gfd <- sapply(fval, gradfd)
    expect_equal(gan, gfd, tol = 1e-4 * max(abs(gan), abs(gfd)))
  }
})

test_that("glm_ridge: gaussian, id-link, intercept, lambda = 0", {
  fam <- extend_family(gaussian(link = "identity"))
  y <- rnorm(n, x %*% b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol
  )
})

test_that("glm_ridge: gaussian, id-link, no intercept, lambda = 0", {
  fam <- extend_family(gaussian(link = "identity"))
  y <- rnorm(n, x %*% b, dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = FALSE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: gaussian, id-link, intercept, lambda = 0.5", {
  fam <- extend_family(gaussian(link = "identity"))
  y <- rnorm(n, x %*% b, dis)
  lambda <- 0.5

  ridgefit <- glm_ridge(x_tr, y,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE
  )
  # analytic solution, no penalty on the intercept term
  penalty <- 0.5 * diag(c(0, rep(lambda, nterms_fit)))
  exp_beta <- c(solve(
    crossprod(cbind(1, x_tr) * sqrt(weights)) + penalty,
    crossprod(cbind(1, x_tr) * weights, y - offset)
  ))

  expect_equal(exp_beta, c(ridgefit$beta0, ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: gaussian, log-link, intercept, lambda = 0", {
  fam <- extend_family(gaussian(link = "log"))
  # intercept of 4 to ensure that y are positive
  y <- rnorm(n, fam$linkinv(x %*% b + 4), dis)
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
                        family = fam, lambda = lambda, weights = weights,
                        offset = offset, intercept = TRUE,
                        thresh = extra_thresh
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol
  )
})

test_that("glm_ridge: binomial, logit-link, intercept, lambda = 0", {
  fam <- extend_family(binomial(link = "logit"))
  y <- rbinom(n, weights, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights - y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y / weights,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol
  )
})

test_that("glm_ridge: binomial, logit-link, no intercept, lambda = 0", {
  fam <- extend_family(binomial(link = "logit"))
  y <- rbinom(n, weights, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights - y) ~ x_tr - 1, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y / weights,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = FALSE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

test_that("glm_ridge: binomial, probit-link, intercept, lambda = 0", {
  fam <- extend_family(binomial(link = "probit"))
  y <- rbinom(n, weights, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights - y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y / weights,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE,
                        thresh = extra_thresh
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol
  )
})

test_that("glm_ridge: binomial, cauchit-link, intercept, lambda = 0", {
  fam <- extend_family(binomial(link = "cauchit"))
  y <- rbinom(n, weights, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(cbind(y, weights - y) ~ x_tr, family = fam, offset = offset)
  ridgefit <- glm_ridge(x_tr, y / weights,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE,
                        thresh = extra_thresh
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol
  )
})

test_that("glm_ridge: poisson, log-link, intercept, lambda = 0", {
  fam <- extend_family(poisson(link = "log"))
  y <- rpois(n, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = TRUE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
               tolerance = tol
  )
})

test_that("glm_ridge: poisson, log-link, no intercept, lambda = 0", {
  fam <- extend_family(poisson(link = "log"))
  y <- rpois(n, fam$linkinv(x %*% b))
  lambda <- 0

  glmfit <- glm(y ~ x_tr - 1, family = fam, weights = weights, offset = offset)
  ridgefit <- glm_ridge(x_tr, y,
                        family = fam, lambda = lambda,
                        weights = weights, offset = offset, intercept = FALSE
  )

  expect_equal(unname(coef(glmfit)), c(ridgefit$beta), tolerance = tol)
})

# test_that("glm_ridge: Gamma, log-link, intercept, lambda = 0", {
#   fam <- extend_family(Gamma(link = 'log'))
#   y <- rgamma(n, fam$linkinv(x%*%b + 1))
#   lambda <- 0
#
#   glmfit <- glm(y ~ x_tr, family = fam, weights = weights, offset = offset)
#   ridgefit <- glm_ridge(x_tr, y, family = fam, lambda = lambda,
#                         weights = weights, offset = offset, intercept = TRUE,
#                         thresh = extra_thresh)
#
#   expect_equal(unname(coef(glmfit)), c(ridgefit$beta0, ridgefit$beta),
#                tolerance = tol)
# })

test_that("glm_ridge: warning messages are returned and printed", {
  fam <- extend_family(binomial())
  x_w <- x[, 1:2]
  y_w <- rbinom(n, 1, fam$linkinv(4 * x_w[, 1] - 3 * x_w[, 2]))
  # Too few quadratic approximation updates:
  expect_output(
    fit_print <- glm_ridge(x_w, y_w, family = fam, qa_updates_max = 2),
    "glm_ridge warning: maximum number of quadratic approximation updates"
  )
  expect_silent(
    fit_quiet <- glm_ridge(x_w, y_w, family = fam, qa_updates_max = 2,
                           print_msgs = FALSE)
  )
  expect_length(fit_quiet$msgs, 1)
  expect_match(fit_quiet$msgs,
               "^glm_ridge warning: maximum number of quadratic approximation")
  expect_identical(fit_print[setdiff(names(fit_print), "msgs")],
                   fit_quiet[setdiff(names(fit_quiet), "msgs")])
  # Converged fit:
  fit_conv <- glm_ridge(x_w, y_w, family = fam, print_msgs = FALSE)
  expect_identical(fit_conv$msgs, character())
  # The submodel fitter turns the messages into one R warning:
  dat_w <- data.frame(y = y_w, x1 = x_w[, 1], x2 = x_w[, 2])
  expect_warning(
    fit_glm_ridge_callback(y ~ x1 + x2, data = dat_w, family = fam,
                           qa_updates_max = 2),
    "maximum number of quadratic approximation updates"
  )
})

test_that("glm_ridge: native pseudo-observations equal the R callback", {
  with_native <- function(native, expr) {
    opts_old <- options(projpred.glm_native = native)
    on.exit(options(opts_old))
    expr
  }
  eta_n <- drop(x[, 1:4] %*% c(0.8, -0.6, 0.4, 0.2))
  fams <- list(
    list(binomial("logit"), rbinom(n, weights, plogis(eta_n)) / weights),
    list(binomial("probit"), rbinom(n, weights, pnorm(eta_n)) / weights),
    list(binomial("cloglog"),
         rbinom(n, weights, 1 - exp(-exp(eta_n - 0.5))) / weights),
    list(binomial("cauchit"), rbinom(n, weights, pcauchy(eta_n)) / weights),
    list(poisson("log"), rpois(n, exp(eta_n / 2))),
    list(gaussian("log"), exp(eta_n / 3) + rnorm(n, 0, 0.1)),
    list(gaussian("identity"), eta_n + rnorm(n))
  )
  for (fam_y in fams) {
    fam <- extend_family(fam_y[[1]])
    y_n <- fam_y[[2]]
    expect_false(is.null(native_pseudo_codes(fam)))
    for (intcpt in c(TRUE, FALSE)) {
      for (x_n in list(x_tr, NULL)) {
        if (is.null(x_n) && !intcpt) next
        fit_args <- list(x = x_n, y = y_n, family = fam, lambda = 1e-4,
                         weights = weights, offset = offset,
                         intercept = intcpt)
        fit_native <- with_native(TRUE, do.call(glm_ridge, fit_args))
        fit_r <- with_native(FALSE, do.call(glm_ridge, fit_args))
        info_str <- paste(fam$family, fam$link, intcpt, is.null(x_n))
        expect_equal(fit_native, fit_r, tolerance = 1e-12, info = info_str)
      }
    }
  }
  # Families without native pseudo-observations:
  expect_null(native_pseudo_codes(extend_family(Gamma())))
  expect_null(native_pseudo_codes(extend_family(poisson("sqrt"))))
  fam_mod <- extend_family(binomial())
  fam_mod$linkinv <- function(eta) 1 / (1 + exp(-eta))
  expect_null(native_pseudo_codes(fam_mod))
  expect_null(with_native(FALSE,
                          native_pseudo_codes(extend_family(binomial()))))
})

test_that(paste(
  "glm_ridge_native_multi() gives the same submodel fits as the draw-wise",
  "fitter"
), {
  n_draws <- 4
  eta_d <- x[, 1:3] %*% matrix(rnorm(3 * n_draws, 0, 0.6), 3, n_draws)
  fams_d <- list(binomial("logit"), binomial("probit"), poisson("log"),
                 gaussian("log"))
  for (fam_d in fams_d) {
    fam <- extend_family(fam_d)
    mu_d <- fam$linkinv(eta_d + if (fam$family == "gaussian") 1 else 0)
    dat_d <- data.frame(x[, 1:4], xf = factor(rep(c("a", "b"), length.out = n)))
    resp_nms <- paste0(".y", seq_len(n_draws))
    dat_d[resp_nms] <- mu_d
    fmls <- lapply(resp_nms, function(resp_nm) {
      as.formula(paste(resp_nm, "~ X1 + X2 + X3 + X4 + xf"))
    })
    shared_x <- shared_model_matrix(fmls, data = dat_d)
    for (qa_max in list(NULL, 2)) {
      dot_args <- list(data = dat_d, family = fam, weights = weights)
      if (!is.null(qa_max)) dot_args$qa_updates_max <- qa_max
      out_multi <- glm_ridge_native_multi(fmls, shared_x = shared_x,
                                          dot_args = dot_args)
      sdivmin_shared <- glm_ridge_callback_shared(shared_x)
      capts_drawwise <- list()
      out_drawwise <- lapply(seq_along(fmls), function(s) {
        capts_drawwise[[s]] <<- capt_conds(
          fit_s <- do.call(sdivmin_shared, c(list(formula = fmls[[s]]),
                                             dot_args))
        )
        fit_s
      })
      info_str <- paste(fam$family, fam$link, is.null(qa_max))
      expect_identical(attr(out_multi, "mssgs_warns_capts"), capts_drawwise,
                       info = info_str)
      attr(out_multi, "mssgs_warns_capts") <- NULL
      expect_equal(out_multi, out_drawwise, tolerance = 1e-12, info = info_str)
      if (!is.null(qa_max) && fam$family != "gaussian") {
        expect_true(any(lengths(capts_drawwise) > 0), info = info_str)
      }
    }
  }
  # Cases that use the draw-wise fitter:
  dot_args <- list(data = dat_d, family = extend_family(Gamma()))
  expect_null(glm_ridge_native_multi(fmls, shared_x, dot_args))
  dot_args <- list(data = dat_d, family = extend_family(binomial()),
                   offset = rep(0, n))
  expect_null(glm_ridge_native_multi(fmls, shared_x, dot_args))
})

test_that("capt_conds() collects messages and warnings", {
  expect_silent(capt <- capt_conds({
    message("a message")
    warning("line 1\nline 2")
    warning("line 1\nline 2")
    invisible(1)
  }))
  expect_identical(capt, c("a message", "line 1", "line 2"))
  expect_identical(capt_conds(invisible(1)), character())
  expect_error(capt_conds(stop("an error")), "an error")
})

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
