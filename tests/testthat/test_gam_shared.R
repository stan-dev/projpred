context("gam_shared")

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

n_gs <- 120L
idx_gs <- seq_len(n_gs)
dat_gs <- data.frame(x1 = sin(idx_gs), x2 = cos(1.3 * idx_gs),
                     x3 = sin(0.7 * idx_gs + 1),
                     f = factor(c("a", "b", "c")[idx_gs %% 3L + 1L]))
n_draws_gs <- 4L
for (s in seq_len(n_draws_gs)) {
  dat_gs[[paste0(".y.", s)]] <- sin(1.5 * dat_gs$x3) + 0.5 * dat_gs$x1 +
    0.3 * cos(s * idx_gs)
}
wts_gs <- 1 + 0.4 * sin(2 * idx_gs)
fml_gs <- as.formula(paste0(
  "cbind(", paste0(".y.", seq_len(n_draws_gs), collapse = ", "), ") ~ ",
  "s(x3) + s(x2) + x1 + f"
))
divmin_gs <- function(shared) {
  opts_old <- options(projpred.gam_shared = shared)
  on.exit(options(opts_old))
  suppressWarnings(divmin(
    formula = fml_gs, data = dat_gs, family = extend_family(gaussian()),
    weights = wts_gs, projpred_var = matrix(NA, n_gs, n_draws_gs),
    projpred_ws_aug = matrix(NA, n_gs, n_draws_gs), verbose_divmin = FALSE
  ))
}

test_that("divmin() gives the same gam() fits with and without sharing", {
  fits_shared <- divmin_gs(TRUE)
  fits_drawwise <- divmin_gs(FALSE)
  expect_length(fits_shared, n_draws_gs)
  for (s in seq_len(n_draws_gs)) {
    expect_s3_class(fits_shared[[s]], "gam")
    expect_equal(fits_shared[[s]], fits_drawwise[[s]], tolerance = 0,
                 info = s)
    expect_identical(formula(fits_shared[[s]]),
                     formula(fits_drawwise[[s]]), info = s)
  }
})

test_that("subprd() gives the predictions of predict.gam() for gam fits", {
  fits <- divmin_gs(FALSE)
  for (nd in list(dat_gs, dat_gs[c(7, 3, 50), ])) {
    nd_icpt <- cbind(`(Intercept)` = rep(1, nrow(nd)), nd)
    prd_expected <- do.call(cbind, lapply(fits, function(fit) {
      predict(fit, newdata = nd_icpt)
    }))
    expect_identical(subprd(fits, nd), prd_expected)
  }
  # Fits of different submodels:
  fit_other <- mgcv::gam(.y.1 ~ s(x3) + x1, data = dat_gs)
  expect_null(subprd_gam_direct(list(fits[[1]], fit_other), dat_gs))
})

test_that("gam_shared_callback() returns `NULL` for other cases", {
  fmls <- validate_response_formula(fml_gs)
  dot_args <- list(data = dat_gs, family = extend_family(gaussian()),
                   weights = wts_gs)
  expect_true(is.function(gam_shared_callback(fmls, dot_args)))
  opts_old <- options(projpred.gam_shared = FALSE)
  expect_null(gam_shared_callback(fmls, dot_args))
  options(opts_old)
  dat_na <- dat_gs
  dat_na$.y.3[5] <- NA
  expect_null(gam_shared_callback(fmls, c(list(data = dat_na),
                                          dot_args[-1])))
})

test_that("shared_setup_callback() builds the setup once, at its first call", {
  n_built <- 0L
  make_shared <- function(formulas, dot_args) {
    n_built <<- n_built + 1L
    n_setup <- nrow(dot_args[["data"]])
    function(formula, ...) c("shared", n_setup)
  }
  make_null <- function(formulas, dot_args) {
    n_built <<- n_built + 1L
    NULL
  }
  plain <- function(formula, ...) c("plain", nrow(list(...)[["data"]]))
  fmls <- validate_response_formula(fml_gs)
  dat_red <- dat_gs[-1, ]
  cb <- shared_setup_callback(make_shared, fmls, plain)
  expect_identical(n_built, 0L)
  expect_identical(cb(fmls[[1]], data = dat_red), c("shared", n_gs - 1L))
  expect_identical(cb(fmls[[2]], data = dat_red), c("shared", n_gs - 1L))
  # Other data than at the first call:
  expect_identical(cb(fmls[[3]], data = dat_gs), c("plain", n_gs))
  expect_identical(n_built, 1L)
  n_built <- 0L
  cb_null <- shared_setup_callback(make_null, fmls, plain)
  expect_identical(cb_null(fmls[[1]], data = dat_red), c("plain", n_gs - 1L))
  expect_identical(cb_null(fmls[[2]], data = dat_red), c("plain", n_gs - 1L))
  expect_identical(n_built, 1L)
})

test_that(paste(
  "divmin() gives the same gam() fits with and without sharing if an",
  "observation has weight 0"
), {
  opts_old <- options(projpred.gam_shared = TRUE)
  on.exit(options(opts_old))
  wts_0 <- replace(wts_gs, 1, 0)
  divmin_w0 <- function(shared) {
    options(projpred.gam_shared = shared)
    suppressWarnings(divmin(
      formula = fml_gs, data = dat_gs, family = extend_family(gaussian()),
      weights = wts_0, projpred_var = matrix(NA, n_gs, n_draws_gs),
      projpred_ws_aug = matrix(NA, n_gs, n_draws_gs), verbose_divmin = FALSE
    ))
  }
  fits_shared <- divmin_w0(TRUE)
  fits_drawwise <- divmin_w0(FALSE)
  for (s in seq_len(n_draws_gs)) {
    expect_equal(fits_shared[[s]], fits_drawwise[[s]], tolerance = 0,
                 info = s)
  }
})

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
