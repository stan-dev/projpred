context("lmer_shared")

# These tests run on CRAN on purpose: with lme4 >= 2.0-1, lmer_shared_fit()
# repeats the steps of lme4::lmer() after a shared lme4::lFormula(), so a change
# in lme4 has to make them fail (also in the reverse-dependency checks of lme4).
# With older lme4, divmin() fits each draw with lmer().

lme4_new_ls <- utils::packageVersion("lme4") >= lmer_shared_lme4_min

test_that("lmer_shared_fit() gives the fits of lme4::lmer()", {
  skip_if_not(lme4_new_ls, paste("lme4 is older than", lmer_shared_lme4_min))
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
  divmin_call <- function() {
    suppressWarnings(divmin(
      formula = fml_ls, data = dat_ls, family = extend_family(gaussian()),
      weights = wts_ls, projpred_var = matrix(NA, n_ls, n_draws_ls),
      projpred_ws_aug = matrix(NA, n_ls, n_draws_ls), verbose_divmin = FALSE
    ))
  }
  if (!shared) {
    return(divmin_call())
  }
  # The shared fits must not fall back to draw-wise fits:
  with_mocked_bindings(
    divmin_call(),
    fit_glmer_callback = function(...) stop("draw-wise fallback")
  )
}

# New data with new levels of `g` and `h`, so that the predictions of the
# lmer() fits use predict.merMod() and the model frame of the fit:
idx_new_ls <- seq_len(6)
dat_new_ls <- data.frame(x1 = cos(idx_new_ls), x2 = sin(2 * idx_new_ls),
                         f = factor(c("a", "b", "c", "a", "b", "c")),
                         g = factor(c("g1", "g2", "g3", "g_new", "g5", "g6")),
                         h = factor(c("h0", "h1", "h_new", "h3", "h4", "h0")))

test_that("divmin() gives the same lmer() fits with and without sharing", {
  skip_if_not(lme4_new_ls, paste("lme4 is older than", lmer_shared_lme4_min))
  fits_shared <- divmin_ls(TRUE)
  fits_drawwise <- divmin_ls(FALSE)
  expect_length(fits_shared, n_draws_ls)
  for (s in seq_len(n_draws_ls)) {
    info_s <- paste("draw", s)
    expect_s4_class(fits_shared[[s]], "lmerMod")
    expect_equal(fits_shared[[s]], fits_drawwise[[s]], tolerance = 0,
                 info = info_s)
    expect_identical(lme4::getME(fits_shared[[s]], "theta"),
                     lme4::getME(fits_drawwise[[s]], "theta"), info = info_s)
    expect_identical(formula(fits_shared[[s]]),
                     formula(fits_drawwise[[s]]), info = info_s)
    # all.equal() does not compare these attributes:
    expect_identical(lmer_frame_attrs(fits_shared[[s]]),
                     lmer_frame_attrs(fits_drawwise[[s]]), info = info_s)
  }
})

test_that(paste(
  "divmin() gives the same lmer() fits with and without sharing if an",
  "observation has weight 0"
), {
  skip_if_not(lme4_new_ls, paste("lme4 is older than", lmer_shared_lme4_min))
  wts_bu <- wts_ls
  on.exit(wts_ls <<- wts_bu)
  wts_ls <<- replace(wts_ls, 1, 0)
  fits_shared <- divmin_ls(TRUE)
  fits_drawwise <- divmin_ls(FALSE)
  for (s in seq_len(n_draws_ls)) {
    info_s <- paste("draw", s)
    expect_equal(fits_shared[[s]], fits_drawwise[[s]], tolerance = 0,
                 info = info_s)
    expect_identical(lme4::getME(fits_shared[[s]], "theta"),
                     lme4::getME(fits_drawwise[[s]], "theta"), info = info_s)
    expect_identical(lmer_frame_attrs(fits_shared[[s]]),
                     lmer_frame_attrs(fits_drawwise[[s]]), info = info_s)
  }
})

test_that(paste(
  "lmer() fits with and without sharing give the same predictions for new",
  "data"
), {
  skip_if_not(lme4_new_ls, paste("lme4 is older than", lmer_shared_lme4_min))
  fits_shared <- divmin_ls(TRUE)
  fits_drawwise <- divmin_ls(FALSE)
  for (s in seq_len(n_draws_ls)) {
    info_s <- paste("draw", s)
    expect_identical(names(model.frame(fits_shared[[s]], fixed.only = TRUE)),
                     names(model.frame(fits_drawwise[[s]], fixed.only = TRUE)),
                     info = info_s)
    expect_equal(
      predict(fits_shared[[s]], newdata = dat_new_ls,
              allow.new.levels = TRUE),
      predict(fits_drawwise[[s]], newdata = dat_new_ls,
              allow.new.levels = TRUE),
      tolerance = 0, info = info_s
    )
  }
})

test_that(paste(
  "proj_linpred() gives the same predictions for new data with and without",
  "sharing the lmer() setup"
), {
  skip_if_not(lme4_new_ls, paste("lme4 is older than", lmer_shared_lme4_min))
  dat_ref <- dat_ls[c("x1", "x2", "f", "g", "h")]
  dat_ref$y <- dat_ls$.y.1
  nobs_ref <- nrow(dat_ref)
  # A custom reference model with 5 draws and a fixed dispersion parameter:
  mu_ref <- sapply(seq_len(5), function(s) {
    dat_ref$y + 0.1 * s * sin(s * seq_len(nobs_ref))
  })
  fit_ref <- list(mu = mu_ref)
  refm <- suppressWarnings(init_refmodel(
    fit_ref, data = dat_ref,
    formula = y ~ x1 + x2 + f + (1 | g) + (1 | h), family = gaussian(),
    ref_predfun = function(fit, newdata = NULL) {
      if (!is.null(newdata)) stop("unexpected call with newdata")
      fit$mu
    },
    dis = rep(0.5, 5)
  ))
  pl <- lapply(c(shared = TRUE, drawwise = FALSE), function(shared) {
    opts_old <- options(projpred.lmer_shared = shared)
    on.exit(options(opts_old))
    project_call <- function() {
      suppressWarnings(project(
        refm, predictor_terms = c("x1", "x2", "f", "(1 | g)", "(1 | h)"),
        ndraws = 5, verbose = 0
      ))
    }
    if (shared) {
      # The shared fits must not fall back to draw-wise fits:
      prj <- with_mocked_bindings(
        project_call(),
        fit_glmer_callback = function(...) stop("draw-wise fallback")
      )
    } else {
      prj <- project_call()
    }
    set.seed(5867)
    proj_linpred(prj, newdata = dat_new_ls)$pred
  })
  expect_equal(pl$shared, pl$drawwise, tolerance = 0)
})

test_that("lmer_shared_callback() returns `NULL` for other cases", {
  fmls <- validate_response_formula(fml_ls)
  dot_args <- list(data = dat_ls, family = extend_family(gaussian()),
                   weights = wts_ls)
  if (lme4_new_ls) {
    expect_true(is.function(lmer_shared_callback(fmls, dot_args)))
  } else {
    lmer_shared_state$ok <- NULL
    expect_false(lmer_shared_ok())
    expect_null(lmer_shared_callback(fmls, dot_args))
  }
  opts_old <- options(projpred.lmer_shared = FALSE)
  expect_null(lmer_shared_callback(fmls, dot_args))
  options(opts_old)
  dat_na <- dat_ls
  dat_na$.y.2[3] <- NA
  expect_null(lmer_shared_callback(fmls, c(list(data = dat_na),
                                           dot_args[-1])))
})
