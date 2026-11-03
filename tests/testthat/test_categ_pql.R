# projpred's PQL fit of multilevel submodels for the brms::categorical() family
# (categ_pql_fit()) against mclogit::mblogit(). These tests need no reference
# model fits, so they also run on CRAN.

# Data with one row per observation and a response matrix of category
# probabilities (the data form of the augmented-data projection).
categ_pql_test_data <- function(seed = 5, n = 200, n_g = 20, sd_g = 0.6) {
  set.seed(seed)
  g <- factor(sample(n_g, n, replace = TRUE))
  h <- factor(paste0("h", sample(8, n, replace = TRUE)))
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  b <- matrix(rnorm(n_g * 2, sd = sd_g), n_g, 2)
  bh <- matrix(rnorm(8 * 2, sd = 0.4), 8, 2)
  eta <- cbind(0,
               0.3 + 0.5 * x1 - 0.3 * x2 + b[g, 1] + bh[h, 1],
               -0.2 - 0.4 * x1 + 0.2 * x2 + b[g, 2] + 0.5 * x1 * b[g, 1] +
                 bh[h, 2])
  P <- exp(eta) / rowSums(exp(eta))
  colnames(P) <- c("a", "b", "c")
  dat <- data.frame(x1 = x1, x2 = x2, g = g, h = h)
  dat$P <- P
  dat
}

# mblogit() with the corrected sum of the diagonal blocks (see the header of
# R/categ_pql.R), so that both fits use the same gradient:
mblogit_fixed <- function(...) {
  fixed_sum_blockDiag <- function(x, n) {
    nblks <- nrow(x) %/% n
    offs <- rep((seq_len(nblks) - 1L) * n, each = n * n)
    i <- rep(rep(seq_len(n), n), nblks) + offs
    j <- rep(rep(seq_len(n), each = n), nblks) + offs
    Matrix::Matrix(matrix(rowSums(matrix(x[cbind(i, j)], nrow = n * n)), n, n))
  }
  local_mocked_bindings(sum_blockDiag = fixed_sum_blockDiag,
                        .package = "mclogit")
  out <- utils::capture.output(fit <- suppressWarnings(mclogit::mblogit(...)))
  fit
}

# `tol_vc`: tolerance for the group-level parameters, which are less precise
# where the criterion is flat (for example a nearly singular covariance).
expect_pql_like_mblogit <- function(fit, fit_ref, tol = 1e-4, tol_vc = tol) {
  expect_s3_class(fit, "mmblogit")
  expect_equal(fit$coefficients, fit_ref$coefficients, tolerance = tol)
  expect_equal(fit$coefmat, fit_ref$coefmat, tolerance = tol)
  expect_equal(unname(fit$VarCov), unname(fit_ref$VarCov), tolerance = tol_vc)
  expect_identical(lapply(fit$VarCov, dimnames),
                   lapply(fit_ref$VarCov, dimnames))
  expect_identical(names(fit$groups), names(fit_ref$groups))
  expect_identical(lapply(fit$groups, levels), lapply(fit_ref$groups, levels))
  expect_equal(lapply(fit$random.effects, as.vector),
               lapply(fit_ref$random.effects, as.vector), tolerance = tol_vc)
  expect_equal(fit$deviance, fit_ref$deviance, tolerance = tol)
  nd <- fit_ref$data[1:7, ]
  expect_equal(predict(fit, newdata = nd, conditional = FALSE),
               predict(fit_ref, newdata = nd, conditional = FALSE),
               tolerance = tol)
}

test_that("categ_pql_fit() agrees with mblogit() for a group-level intercept", {
  skip_if_not_installed("mclogit")
  dat <- categ_pql_test_data()
  fit_ref <- mblogit_fixed(P ~ x1 + x2, random = ~ 1 | g, data = dat)
  fit <- suppressWarnings(
    categ_pql_fit(P ~ x1 + x2, random = list(~ 1 | g), data = dat)
  )
  expect_pql_like_mblogit(fit, fit_ref)
})

test_that("categ_pql_fit() agrees with mblogit() for a slope and two groups", {
  skip_if_not_installed("mclogit")
  dat <- categ_pql_test_data(seed = 13)
  fit_ref <- mblogit_fixed(P ~ x1 + x2, random = list(~ 1 + x1 | g, ~ 1 | h),
                           data = dat)
  fit <- suppressWarnings(
    categ_pql_fit(P ~ x1 + x2, random = list(~ 1 + x1 | g, ~ 1 | h),
                  data = dat)
  )
  # The 4 x 4 covariance matrix of `g` is nearly singular:
  expect_pql_like_mblogit(fit, fit_ref, tol = 1e-3, tol_vc = 0.1)
})

test_that("the gradient of the inner PQL objective is correct", {
  set.seed(3)
  d <- c(2L, 1L)
  m <- c(3L, 2L)
  p <- 2L
  q <- sum(d * m)
  A <- matrix(rnorm(40 * (p + q + 1)), 40)
  G <- crossprod(A)
  ix <- seq_len(p)
  iz <- p + seq_len(q)
  iy <- p + q + 1L
  z_idx <- split(seq_len(q), rep(seq_along(d), d * m))
  pll <- categ_pql_objective(
    XWX = G[ix, ix], XWy = G[ix, iy], ZWX = G[iz, ix], ZWZ = G[iz, iz],
    ZWy = G[iz, iy], yWy = G[iy, iy], z_idx = z_idx, d = d, m = m
  )
  lambda <- c(1.3, 0.4, 0.8, 1.1)
  grad <- pll(lambda, gradient = TRUE)$gradient
  h <- 1e-6
  grad_num <- vapply(seq_along(lambda), function(i) {
    e <- replace(numeric(length(lambda)), i, h)
    (pll(lambda + e)$logLik - pll(lambda - e)$logLik) / (2 * h)
  }, numeric(1))
  expect_equal(grad, grad_num, tolerance = 1e-6)
})

test_that("fit_categ_mlvl() uses categ_pql_fit() for augmented data", {
  skip_if_not_installed("mclogit")
  dat <- categ_pql_test_data()
  P <- dat$P
  dat$P <- NULL
  # Augmented data: the observations in one block per response category.
  dat_aug <- dat[rep(seq_len(nrow(dat)), 3), ]
  dat_aug$y <- factor(rep(colnames(P), each = nrow(dat)))
  args <- list(formula = y ~ x1 + x2 + (1 | g),
               projpred_formula_no_random = y ~ x1 + x2,
               projpred_random = ~ (1 | g), data = dat_aug,
               family = list(link = "logit"), weights = as.vector(P))
  fit <- suppressWarnings(do.call(fit_categ_mlvl, args))
  expect_s3_class(fit, "mmblogit")
  expect_identical(fit$response.type, "matrix")
  # The global option selects mblogit():
  opts_old <- options(projpred.categ_mlvl_mblogit = TRUE)
  on.exit(options(opts_old), add = TRUE)
  fit_mb <- suppressWarnings(do.call(fit_categ_mlvl, args))
  expect_identical(fit_mb$response.type, "factor")
  # The parts of the fit that projpred uses have the same structure:
  expect_identical(names(get_subparams(fit, nm_scheme = "brms")),
                   names(get_subparams(fit_mb, nm_scheme = "brms")))
  expect_identical(dim(repair_re(fit, newdata = dat[1:4, ])),
                   dim(repair_re(fit_mb, newdata = dat[1:4, ])))
  # Other settings for mblogit() also select mblogit():
  options(opts_old)
  fit_mb2 <- suppressWarnings(do.call(fit_categ_mlvl,
                                      c(args, list(inner.optimizer = "BFGS"))))
  expect_identical(fit_mb2$response.type, "factor")
})

test_that("categ_mat_resp() needs the augmented-data layout", {
  dat <- data.frame(x = rep(1:3, 2), y = factor(rep(c("a", "b"), each = 3)))
  res <- categ_mat_resp(y ~ x, dat, weights = c(0.2, 0.5, 0.9, 0.8, 0.5, 0.1))
  expect_identical(deparse(res$formula), "projpred_internal_resp_mat ~ x")
  expect_identical(dim(res$data$projpred_internal_resp_mat), c(3L, 2L))
  expect_identical(colnames(res$data$projpred_internal_resp_mat), c("a", "b"))
  expect_null(categ_mat_resp(y ~ x, dat[c(1, 4, 2, 5, 3, 6), ],
                             weights = rep(0.5, 6)))
  dat_diff <- dat
  dat_diff$x[4] <- 7L
  expect_null(categ_mat_resp(y ~ x, dat_diff, weights = rep(0.5, 6)))
})
