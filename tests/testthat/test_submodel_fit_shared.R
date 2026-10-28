context("shared model matrix and joint Gaussian solve for submodel fits")

# Data with one response column per projected draw, as divmin() gets them, and
# the draw-wise formulas. The formulas include a factor and an interaction.
shared_fit_data <- function(n = 40, S = 5, seed = 4721) {
  set.seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                    f = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  for (s in seq_len(S)) {
    dat[[paste0("y", s)]] <- 0.5 * dat$x1 - 0.3 * dat$x2 + (dat$f == "b") +
      rnorm(n, sd = 0.3)
  }
  fmls <- lapply(seq_len(S), function(s) {
    as.formula(paste0("y", s, " ~ x1 + x2 + f + x1:x2"))
  })
  list(data = dat, formulas = fmls, w = runif(n, 0.5, 2))
}

test_that(paste(
  "the joint Gaussian solve gives the same fits as the draw-wise",
  "fit_glm_ridge_callback()"
), {
  d <- shared_fit_data()
  shared_x <- shared_model_matrix(d$formulas, data = d$data)
  expect_false(is.null(shared_x))
  fits_multi <- glm_ridge_gauss_multi(
    d$formulas, shared_x = shared_x,
    dot_args = list(data = d$data, family = gaussian(), weights = d$w)
  )
  expect_length(fits_multi, length(d$formulas))
  fits_loop <- lapply(d$formulas, fit_glm_ridge_callback, data = d$data,
                      family = gaussian(), weights = d$w)
  for (s in seq_along(d$formulas)) {
    expect_s3_class(fits_multi[[s]], "subfit")
    expect_identical(names(fits_multi[[s]]), names(fits_loop[[s]]))
    expect_equal(unname(fits_multi[[s]]$alpha), unname(fits_loop[[s]]$alpha),
                 tolerance = 1e-10)
    expect_equal(fits_multi[[s]]$beta, fits_loop[[s]]$beta, tolerance = 1e-10)
    expect_identical(fits_multi[[s]]$x, fits_loop[[s]]$x)
    expect_identical(fits_multi[[s]]$y, fits_loop[[s]]$y)
    expect_identical(fits_multi[[s]]$xlvls, fits_loop[[s]]$xlvls)
    expect_equal(fits_multi[[s]]$w, fits_loop[[s]]$w)
  }
})

test_that("the joint Gaussian solve falls back for other families", {
  d <- shared_fit_data()
  shared_x <- shared_model_matrix(d$formulas, data = d$data)
  expect_null(glm_ridge_gauss_multi(
    d$formulas, shared_x = shared_x,
    dot_args = list(data = d$data, family = gaussian(link = "log"))
  ))
  expect_null(glm_ridge_gauss_multi(
    d$formulas, shared_x = shared_x,
    dot_args = list(data = d$data, family = poisson())
  ))
})

test_that(paste(
  "the joint Gaussian solve falls back for a singular system or non-finite",
  "coefficients"
), {
  d <- shared_fit_data()
  # Without penalty, a duplicated column makes the system singular:
  dat_dup <- d$data
  dat_dup$x3 <- dat_dup$x1
  fmls_dup <- lapply(d$formulas, update, . ~ . + x3)
  shared_x_dup <- shared_model_matrix(fmls_dup, data = dat_dup)
  expect_false(is.null(shared_x_dup))
  expect_null(glm_ridge_gauss_multi(
    fmls_dup, shared_x = shared_x_dup,
    dot_args = list(data = dat_dup, family = gaussian(), weights = d$w,
                    regul = 0)
  ))
  # With the default penalty, the same system has a solution:
  expect_false(is.null(glm_ridge_gauss_multi(
    fmls_dup, shared_x = shared_x_dup,
    dot_args = list(data = dat_dup, family = gaussian(), weights = d$w)
  )))
  # An infinite response gives non-finite coefficients:
  dat_inf <- d$data
  dat_inf$y2[4] <- Inf
  shared_x_inf <- shared_model_matrix(d$formulas, data = dat_inf)
  expect_false(is.null(shared_x_inf))
  expect_null(glm_ridge_gauss_multi(
    d$formulas, shared_x = shared_x_inf,
    dot_args = list(data = dat_inf, family = gaussian(), weights = d$w)
  ))
})

test_that("shared_model_matrix() returns NULL if the model matrix differs", {
  d <- shared_fit_data()
  # Different right-hand sides:
  fmls_rhs <- d$formulas
  fmls_rhs[[2]] <- y2 ~ x1 + x2
  expect_null(shared_model_matrix(fmls_rhs, data = d$data))
  # A response that is not a plain column of the data:
  fmls_resp <- d$formulas
  fmls_resp[[1]] <- log(abs(y1)) ~ x1 + x2 + f + x1:x2
  expect_null(shared_model_matrix(fmls_resp, data = d$data))
  # A missing predictor value:
  dat_na <- d$data
  dat_na$x2[3] <- NA
  expect_null(shared_model_matrix(d$formulas, data = dat_na))
})

test_that(paste(
  "the joint predictions for new data equal the draw-wise predictions of",
  "predict.subfit()"
), {
  d <- shared_fit_data()
  shared_x <- shared_model_matrix(d$formulas, data = d$data)
  fits <- glm_ridge_gauss_multi(
    d$formulas, shared_x = shared_x,
    dot_args = list(data = d$data, family = gaussian(), weights = d$w)
  )
  nd <- d$data[c(2, 5, 7, 11), ]
  prd_loop <- do.call(cbind, lapply(fits, predict, newdata = nd))
  prd_shared <- subprd_subfit_shared(fits, nd)
  expect_equal(unname(prd_shared), unname(prd_loop), tolerance = 1e-12)
  expect_equal(unname(subprd(fits, nd)), unname(prd_loop), tolerance = 1e-12)
})

test_that("project() uses the joint Gaussian solve for a Gaussian submodel", {
  d <- shared_fit_data(n = 60)
  dat <- d$data[, c("x1", "x2", "f")]
  dat$y <- d$data$y1
  fml <- y ~ x1 + x2 + f
  mm <- model.matrix(fml, dat)
  ndraws <- 100
  beta <- drop(solve(crossprod(mm), crossprod(mm, dat$y))) +
    matrix(rnorm(ncol(mm) * ndraws, sd = 0.05), ncol(mm))
  fit <- list(beta = beta, data = dat)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(~ x1 + x2 + f, nd) %*% fit$beta
  }
  refm <- suppressWarnings(init_refmodel(
    fit, data = dat, formula = fml, family = gaussian(),
    ref_predfun = ref_predfun, dis = rep(0.3, ndraws)
  ))
  calls <- new.env()
  calls$n <- 0L
  ns <- environment(glm_ridge_gauss_multi)
  suppressMessages(trace(
    "glm_ridge_gauss_multi", where = ns, print = FALSE,
    exit = bquote(if (!is.null(returnValue(NULL))) {
      assign("n", get("n", envir = .(calls)) + 1L, envir = .(calls))
    })
  ))
  on.exit(suppressMessages(untrace("glm_ridge_gauss_multi", where = ns)),
          add = TRUE)
  prj <- project(refm, predictor_terms = c("x1", "f"), nclusters = 5,
                 seed = 3906, verbose = FALSE)
  expect_gt(calls$n, 0L)
  expect_length(prj$outdmin, 5)
})
