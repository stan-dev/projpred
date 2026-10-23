context("formula and model-matrix shortcuts of the projections")

test_that(paste(
  "validate_response_formula() gives the same formulas as one update() per",
  "response"
), {
  rhs_list <- c(
    "x1 + x2", "x1 * f", "x1 + offset(o)", "x1 + (1 | g)",
    "x1 + (x2 | g) + (1 | h)", "s(x1) + x2", "1"
  )
  resps <- paste0(".y_draw.", 1:4)
  for (rhs in rhs_list) {
    fml <- as.formula(paste0("cbind(", paste(resps, collapse = ", "), ") ~ ",
                             rhs))
    fmls_expected <- lapply(resps, function(r) update(fml, paste0(r, " ~ .")))
    expect_identical(validate_response_formula(fml), fmls_expected,
                     info = rhs)
  }
})

test_that(paste(
  "subprd_fitted() gives the training-data predictions of subprd() and",
  "returns NULL if a model matrix has fewer rows than the data"
), {
  set.seed(8413)
  n <- 30
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                    f = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  dat$y1 <- dat$x1 - 0.5 * dat$x2 + (dat$f == "b") + rnorm(n, sd = 0.3)
  dat$y2 <- dat$y1 + rnorm(n, sd = 0.1)
  fmls <- list(y1 ~ x1 * f + x2, y2 ~ x1 * f + x2)
  fits <- lapply(fmls, fit_glm_ridge_callback, data = dat,
                 family = extend_family(gaussian()), weights = rep(1, n))
  expect_equal(subprd_fitted(fits, nobs = n), subprd(fits, newdata = dat),
               tolerance = 1e-12)
  # Intercept-only submodels:
  fits_icpt <- lapply(list(y1 ~ 1, y2 ~ 1), fit_glm_ridge_callback,
                      data = dat, family = extend_family(gaussian()),
                      weights = rep(1, n))
  expect_equal(subprd_fitted(fits_icpt, nobs = n),
               subprd(fits_icpt, newdata = dat), tolerance = 1e-12)
  # If model.frame() removed rows with missing values, the stored model
  # matrices have fewer rows than the data:
  fits_short <- lapply(fits, function(fit) {
    fit$x <- fit$x[-5, , drop = FALSE]
    fit
  })
  expect_null(subprd_fitted(fits_short, nobs = n))
})
