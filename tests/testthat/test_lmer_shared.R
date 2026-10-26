context("lmer_shared")

# These tests run on CRAN on purpose: lmer_shared_fit() repeats the steps of
# lme4::lmer() after a shared lme4::lFormula(), so a change in lme4 has to make
# them fail (also in the reverse-dependency checks of lme4).

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

test_that("lmer_shared_fit() gives the fits of lme4::lmer()", {
  expect_true(lmer_shared_check())
  lmer_shared_state$ok <- NULL
  expect_true(lmer_shared_ok())
})

n_ls <- 80L
idx_ls <- seq_len(n_ls)
dat_ls <- data.frame(x1 = sin(idx_ls), x2 = cos(1.3 * idx_ls),
                     f = factor(c("a", "b", "c")[idx_ls %% 3L + 1L]),
                     g = factor(paste0("g", idx_ls %% 8L)),
                     h = factor(paste0("h", idx_ls %% 5L)))
u_ls <- sin(1:8)
n_draws_ls <- 4L
for (s in seq_len(n_draws_ls)) {
  dat_ls[[paste0(".y.", s)]] <- dat_ls$x1 - 0.5 * dat_ls$x2 +
    u_ls[as.integer(dat_ls$g)] * (1 + 0.2 * s) + 0.3 * cos(s * idx_ls) +
    0.3 * dat_ls$x1 * u_ls[as.integer(dat_ls$g)]
}
wts_ls <- 1 + 0.4 * sin(2 * idx_ls)
fml_ls <- as.formula(paste0(
  "cbind(", paste0(".y.", seq_len(n_draws_ls), collapse = ", "), ") ~ ",
  "x1 + poly(x2, 2, raw = TRUE) + f + (1 + x1 | g) + (1 | h)"
))
divmin_ls <- function(shared) {
  opts_old <- options(projpred.lmer_shared = shared)
  on.exit(options(opts_old))
  suppressWarnings(divmin(
    formula = fml_ls, data = dat_ls, family = extend_family(gaussian()),
    weights = wts_ls, projpred_var = matrix(NA, n_ls, n_draws_ls),
    projpred_ws_aug = matrix(NA, n_ls, n_draws_ls), verbose_divmin = FALSE
  ))
}

test_that("divmin() gives the same lmer() fits with and without sharing", {
  fits_shared <- divmin_ls(TRUE)
  fits_drawwise <- divmin_ls(FALSE)
  expect_length(fits_shared, n_draws_ls)
  for (s in seq_len(n_draws_ls)) {
    expect_s4_class(fits_shared[[s]], "lmerMod")
    expect_equal(fits_shared[[s]], fits_drawwise[[s]], tolerance = 0,
                 info = s)
    expect_identical(lme4::getME(fits_shared[[s]], "theta"),
                     lme4::getME(fits_drawwise[[s]], "theta"), info = s)
    expect_identical(formula(fits_shared[[s]]),
                     formula(fits_drawwise[[s]]), info = s)
  }
})

test_that(paste(
  "divmin() gives the same lmer() fits with and without sharing if an",
  "observation has weight 0"
), {
  wts_bu <- wts_ls
  on.exit(wts_ls <<- wts_bu)
  wts_ls <<- replace(wts_ls, 1, 0)
  fits_shared <- divmin_ls(TRUE)
  fits_drawwise <- divmin_ls(FALSE)
  for (s in seq_len(n_draws_ls)) {
    expect_equal(fits_shared[[s]], fits_drawwise[[s]], tolerance = 0,
                 info = s)
    expect_identical(lme4::getME(fits_shared[[s]], "theta"),
                     lme4::getME(fits_drawwise[[s]], "theta"), info = s)
  }
})

test_that("lmer_shared_callback() returns `NULL` for other cases", {
  fmls <- validate_response_formula(fml_ls)
  dot_args <- list(data = dat_ls, family = extend_family(gaussian()),
                   weights = wts_ls)
  expect_true(is.function(lmer_shared_callback(fmls, dot_args)))
  opts_old <- options(projpred.lmer_shared = FALSE)
  expect_null(lmer_shared_callback(fmls, dot_args))
  options(opts_old)
  dat_na <- dat_ls
  dat_na$.y.2[3] <- NA
  expect_null(lmer_shared_callback(fmls, c(list(data = dat_na),
                                           dot_args[-1])))
})

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
