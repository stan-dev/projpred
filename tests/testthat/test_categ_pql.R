# Tests of categ_pql_fit(), projpred's PQL fit of multilevel submodels for the
# brms::categorical() family, against mclogit::mblogit(). These tests need no
# reference model fits, so they also run on CRAN.

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
  skip_if_not(exists("sum_blockDiag", envir = asNamespace("mclogit"),
                     inherits = FALSE))
  fixed_sum_blockDiag <- function(x, n) {
    nblks <- nrow(x) %/% n
    offs <- rep((seq_len(nblks) - 1L) * n, each = n * n)
    i <- rep(rep(seq_len(n), n), nblks) + offs
    j <- rep(rep(seq_len(n), each = n), nblks) + offs
    Matrix::Matrix(matrix(rowSums(matrix(x[cbind(i, j)], nrow = n * n)), n, n))
  }
  local_mocked_bindings(sum_blockDiag = fixed_sum_blockDiag,
                        .package = "mclogit")
  utils::capture.output(fit <- suppressWarnings(mclogit::mblogit(...)))
  fit
}

# `tol_vc`: tolerance for the group-level parameters, which are less precise
# where the criterion is flat (for example a nearly singular covariance).
expect_pql_like_mblogit <- function(fit, fit_ref, tol = 1e-4, tol_vc = tol) {
  expect_s3_class(fit, "mmblogit")
  expect_equal(fit$coefficients, fit_ref$coefficients, tolerance = tol)
  expect_equal(fit$coefmat, fit_ref$coefmat, tolerance = tol)
  expect_equal(fit$VarCov, fit_ref$VarCov, tolerance = tol_vc)
  expect_identical(names(fit$groups), names(fit_ref$groups))
  expect_identical(lapply(fit$groups, levels), lapply(fit_ref$groups, levels))
  ranef_ref <- lapply(seq_along(fit_ref$random.effects), function(k) {
    fit_ref$random.effects[[k]]
  })
  expect_equal(fit$random.effects, ranef_ref, tolerance = tol_vc)
  expect_equal(fit$deviance, fit_ref$deviance, tolerance = tol)
  expect_identical(fit$random, fit_ref$random)
  expect_identical(fit$D, fit_ref$D)
  # predict.mmblogit() fails for factor levels that the fit does not know:
  nd <- droplevels(fit_ref$data[1:7, ])
  expect_equal(predict(fit, newdata = nd, conditional = FALSE),
               predict(fit_ref, newdata = nd, conditional = FALSE),
               tolerance = tol)
}

# Where a group-level covariance matrix is nearly singular, the criterion is
# flat and nlminb() can report a singular convergence.  Allow only this warning:
with_inner_warnings_only <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    expect_match(conditionMessage(w), "^Inner iterations did not converge")
    invokeRestart("muffleWarning")
  })
}

test_that("categ_pql_fit() agrees with mblogit() for a group-level intercept", {
  dat <- categ_pql_test_data()
  fit_ref <- mblogit_fixed(P ~ x1 + x2, random = ~ 1 | g, data = dat)
  fit <- with_inner_warnings_only(
    categ_pql_fit(P ~ x1 + x2, random = list(~ 1 | g), data = dat)
  )
  expect_true(fit$converged)
  expect_pql_like_mblogit(fit, fit_ref)
})

test_that("categ_pql_fit() agrees with mblogit() for a slope and two groups", {
  dat <- categ_pql_test_data(seed = 13)
  fit_ref <- mblogit_fixed(P ~ x1 + x2, random = list(~ 1 + x1 | g, ~ 1 | h),
                           data = dat)
  fit <- with_inner_warnings_only(
    categ_pql_fit(P ~ x1 + x2, random = list(~ 1 + x1 | g, ~ 1 | h),
                  data = dat)
  )
  expect_true(fit$converged)
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

test_that("categ_pql_fit() handles zero weights, factors and character groups", {
  dat <- categ_pql_test_data()
  dat$P[1:3, ] <- 0
  dat$f <- factor(rep(c("u", "v"), length.out = nrow(dat)),
                  levels = c("u", "v", "unused"))
  dat$g_chr <- as.character(dat$g)
  fit <- with_inner_warnings_only(
    categ_pql_fit(P ~ x1 + f, random = list(~ 1 | g_chr), data = dat)
  )
  expect_identical(fit$xlevels, list(f = c("u", "v")))
  expect_identical(levels(fit$groups$g_chr), sort(unique(dat$g_chr)))
  fit_ref <- mblogit_fixed(P ~ x1 + f, random = ~ 1 | g_chr, data = dat)
  expect_pql_like_mblogit(fit, fit_ref)
})

test_that("categ_pql_fit() agrees with mblogit() for two categories", {
  # One logit, so the group-level covariance blocks are 1 x 1:
  dat <- categ_pql_test_data()
  dat$P <- cbind(a = dat$P[, "a"], b = dat$P[, "b"] + dat$P[, "c"])
  fit_ref <- mblogit_fixed(P ~ x1 + x2, random = ~ 1 | g, data = dat)
  fit <- with_inner_warnings_only(
    categ_pql_fit(P ~ x1 + x2, random = list(~ 1 | g), data = dat)
  )
  expect_identical(rownames(fit$random.effects[[1]]), levels(dat$g))
  expect_pql_like_mblogit(fit, fit_ref)
})

# categ_pql_inner() or categ_pql_eval() fails at call `k`.
fail_at <- function(fun, k, fail) {
  # Evaluate `fun` now, before the mock replaces the binding:
  force(fun)
  n_calls <- 0L
  function(...) {
    n_calls <<- n_calls + 1L
    if (n_calls == k) return(fail(...))
    fun(...)
  }
}

test_that("categ_pql_fit() steps back after numeric problems", {
  dat <- categ_pql_test_data()
  args <- list(formula = P ~ x1 + x2, random = list(~ 1 | g), data = dat)
  fit_1 <- suppressWarnings(do.call(categ_pql_fit, c(args, list(maxit = 1))))
  inner_fail <- function(...) {
    categ_pql_numeric(stop("singular test matrix"))
  }
  # In the first iteration, the caller falls back to mblogit():
  local({
    local_mocked_bindings(categ_pql_inner = fail_at(categ_pql_inner, 1,
                                                    inner_fail))
    expect_null(do.call(categ_pql_fit, args))
  })
  # Later, the fit steps back to the parameters of the previous iteration:
  local({
    local_mocked_bindings(categ_pql_inner = fail_at(categ_pql_inner, 2,
                                                    inner_fail))
    msgs <- character()
    fit <- withCallingHandlers(do.call(categ_pql_fit, args),
                               warning = function(w) {
                                 msgs <<- c(msgs, conditionMessage(w))
                                 invokeRestart("muffleWarning")
                               })
    expect_true(any(grepl(paste0("^Numeric problems in inner iteration ",
                                 "\\(singular test matrix\\), stepping back"),
                          msgs)))
    expect_true(any(msgs == "Algorithm stopped without convergence"))
    expect_false(fit$converged)
    expect_equal(fit$coefficients, fit_1$coefficients)
    expect_equal(fit$linear.predictors, fit_1$linear.predictors)
  })
  # Other errors are not caught:
  local({
    local_mocked_bindings(categ_pql_inner = fail_at(
      categ_pql_inner, 2, function(...) stop("programming error")
    ))
    expect_error(suppressWarnings(do.call(categ_pql_fit, args)),
                 "programming error")
  })
})

test_that("categ_pql_fit() backs up after a non-finite deviance", {
  dat <- categ_pql_test_data()
  args <- list(formula = P ~ x1 + x2, random = list(~ 1 | g), data = dat)
  fit_1 <- suppressWarnings(do.call(categ_pql_fit, c(args, list(maxit = 1))))
  eval_orig <- categ_pql_eval
  eval_inf <- function(...) {
    res <- eval_orig(...)
    res$deviance <- Inf
    res
  }
  local({
    local_mocked_bindings(categ_pql_eval = fail_at(categ_pql_eval, 1,
                                                   eval_inf))
    expect_null(with_inner_warnings_only(do.call(categ_pql_fit, args)))
  })
  local({
    local_mocked_bindings(categ_pql_eval = fail_at(categ_pql_eval, 2,
                                                   eval_inf))
    msgs <- character()
    fit <- withCallingHandlers(do.call(categ_pql_fit, args),
                               warning = function(w) {
                                 msgs <<- c(msgs, conditionMessage(w))
                                 invokeRestart("muffleWarning")
                               })
    expect_true(any(msgs == "Non-finite deviance, backing up"))
    expect_false(fit$converged)
    expect_equal(fit$coefficients, fit_1$coefficients)
    expect_equal(fit$linear.predictors, fit_1$linear.predictors)
  })
})

test_that("categ_pql_fit() returns NULL for cases left to mblogit()", {
  dat <- categ_pql_test_data()
  dat$x1_copy <- dat$x1
  # Rank-deficient model matrix:
  expect_null(categ_pql_fit(P ~ x1 + x1_copy, random = list(~ 1 | g),
                            data = dat))
  # A group variable in two terms:
  expect_null(categ_pql_fit(P ~ x1, random = list(~ 1 | g, ~ 0 + x1 | g),
                            data = dat))
  # A nested group term:
  expect_null(categ_pql_fit(P ~ x1, random = list(~ 1 | h / g), data = dat))
  # Too many cells in the group-level design matrix:
  expect_null(categ_pql_fit(P ~ x1, random = list(~ 1 | g), data = dat,
                            max_cells = 100))
})

# Augmented data: the observations in one block per response category.
categ_pql_aug_args <- function(dat) {
  P <- dat$P
  dat$P <- NULL
  dat_aug <- dat[rep(seq_len(nrow(dat)), ncol(P)), ]
  dat_aug$y <- factor(rep(colnames(P), each = nrow(dat)))
  list(formula = y ~ x1 + x2 + (1 | g),
       projpred_formula_no_random = y ~ x1 + x2,
       projpred_random = ~ (1 | g), data = dat_aug,
       family = list(link = "logit"), weights = as.vector(P))
}

test_that("fit_categ_mlvl() uses categ_pql_fit() for augmented data", {
  dat <- categ_pql_test_data()
  args <- categ_pql_aug_args(dat)
  fit <- suppressWarnings(do.call(fit_categ_mlvl, args))
  expect_identical(fit$response.type, "matrix")
  # `avoid.increase` has no effect in the PQL algorithm of mclogit:
  fit_ai <- suppressWarnings(
    do.call(fit_categ_mlvl, c(args, list(avoid.increase = TRUE)))
  )
  expect_identical(fit_ai$response.type, "matrix")
  # Other mblogit() settings select mblogit():
  fit_mb <- suppressWarnings(
    do.call(fit_categ_mlvl, c(args, list(inner.optimizer = "BFGS")))
  )
  expect_identical(fit_mb$response.type, "factor")
  local({
    local_mocked_bindings(mblogit = function(...) {
      cat("converged\n")
      "mblogit() called"
    }, .package = "mclogit")
    expect_identical(
      do.call(fit_categ_mlvl, c(args, list(break.on.increase = TRUE))),
      "mblogit() called"
    )
  })
  # The parts of the fit that projpred uses have the same structure:
  expect_identical(names(get_subparams(fit, nm_scheme = "brms")),
                   names(get_subparams(fit_mb, nm_scheme = "brms")))
  expect_identical(dim(repair_re(fit, newdata = dat[1:4, ])),
                   dim(repair_re(fit_mb, newdata = dat[1:4, ])))
  # Data in another layout go to mblogit():
  args_shuffled <- args
  args_shuffled$data <- args$data[rev(seq_len(nrow(args$data))), ]
  args_shuffled$weights <- rev(args$weights)
  fit_shuffled <- suppressWarnings(do.call(fit_categ_mlvl, args_shuffled))
  expect_identical(fit_shuffled$response.type, "factor")
})

test_that("global option `projpred.categ_mlvl_mblogit` selects mblogit()", {
  opts_old <- options(projpred.categ_mlvl_mblogit = TRUE)
  on.exit(options(opts_old), add = TRUE)
  fit <- suppressWarnings(
    do.call(fit_categ_mlvl, categ_pql_aug_args(categ_pql_test_data()))
  )
  expect_identical(fit$response.type, "factor")
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
