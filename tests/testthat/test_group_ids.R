# Group terms with a brms correlation ID (`(1 |p| g)`), group terms with `||`,
# and several group terms of one grouping factor. These tests run on CRAN: the
# reference model is a custom one (no Stan fit).

# Data with 8 groups of 10 observations, and draws of the linear predictor:
# the fitted values of an lmer() fit plus a small noise.
ids_inputs <- withr::with_seed(4, {
  g <- factor(paste0("g", rep(1:8, each = 10)))
  x <- rnorm(80)
  b0 <- rnorm(8, sd = 0.8)
  b1 <- rnorm(8, sd = 0.5)
  dat <- data.frame(y = 0.5 + 0.7 * x + b0[g] + b1[g] * x +
                      rnorm(80, sd = 0.5),
                    x = x, z = rnorm(80), g = g)
  fit_lmer <- lme4::lmer(y ~ x + (x | g), data = dat)
  eta <- fitted(fit_lmer) + matrix(rnorm(80 * 20, sd = 0.05), 80, 20)
  list(dat = dat, eta = eta)
})

refm_ids <- function(formula, ...) {
  # Without the test setup, a custom reference model without `cvrefbuilder`
  # gives a warning:
  withr::local_options(projpred.warn_cvrefbuilder_NULL = FALSE)
  init_refmodel(list(eta = ids_inputs$eta), data = ids_inputs$dat,
                formula = formula, family = gaussian(),
                ref_predfun = function(fit, newdata = NULL) {
                  stopifnot(is.null(newdata))
                  fit$eta
                },
                dis = rep(0.5, ncol(ids_inputs$eta)), ...)
}

test_that("rm_group_ids() removes the correlation IDs of group terms", {
  fml <- y ~ x + (1 |p| g) + (0 + x |p| g) + (1 | h) + I(m[, 1])
  expect_message(out <- rm_group_ids(fml), "`p`")
  expect_identical(
    out$formula,
    `environment<-`(y ~ x + (1 | g) + (0 + x | g) + (1 | h) + I(m[, 1]),
                    environment(fml))
  )
  expect_identical(out$grp_fctrs_shared_id, "g")
  expect_silent(out_q <- rm_group_ids(fml, verbose = FALSE))
  expect_identical(out_q, out)

  # No ID: no change and no message.
  fml_noid <- y ~ x + (1 + x | g)
  expect_silent(out_noid <- rm_group_ids(fml_noid))
  expect_identical(out_noid$formula, fml_noid)
  expect_length(out_noid$grp_fctrs_shared_id, 0)

  # One-sided formula, a two-column response and an offset:
  expect_identical(rm_group_ids(~ (1 |p| g), verbose = FALSE)$formula,
                   `environment<-`(~ (1 | g), environment(fml)))
  fml_cb <- cbind(s, f) ~ x + (x |q| g) + offset(o)
  expect_identical(
    rm_group_ids(fml_cb, verbose = FALSE)$formula,
    `environment<-`(cbind(s, f) ~ x + (x | g) + offset(o), environment(fml))
  )

  # A `|` inside another call does not change:
  fml_or <- y ~ I(a | b | c) + (1 | g)
  expect_silent(out_or <- rm_group_ids(fml_or))
  expect_identical(out_or$formula, fml_or)

  # Terms of one grouping factor with other IDs or without an ID:
  fml_mix <- y ~ (1 |p| g) + (0 + x | g) + (1 |p| h) + (0 + x |q| h)
  expect_length(rm_group_ids(fml_mix, verbose = FALSE)$grp_fctrs_shared_id, 0)
})

test_that("extract_terms_response() gives informative errors", {
  expect_error(extract_terms_response(y ~ x + (x || g)),
               "does not support group terms with `||`", fixed = TRUE)
  expect_error(extract_terms_response(y ~ x + (1 |p| g)),
               "more than one `|`", fixed = TRUE)
})

test_that(paste(
  "init_refmodel() removes the correlation ID with a message, and project()",
  "works"
), {
  expect_message(refm <- refm_ids(y ~ x + (1 + x |p| g)),
                 "correlation ID(s) `p`", fixed = TRUE)
  fml_txt <- paste(deparse(refm$formula), collapse = " ")
  expect_true(grepl("(1 + x | g)", fml_txt, fixed = TRUE), info = fml_txt)
  # No message for the fold-wise reference models:
  expect_no_message(refm_ids(y ~ x + (1 + x |p| g),
                             called_from_cvrefbuilder = TRUE))
  prj <- project(refm, predictor_terms = c("x", "(x | g)"), nclusters = 2,
                 verbose = FALSE)
  expect_s4_class(prj$outdmin[[1]], "lmerMod")
  expect_identical(names(lme4::VarCorr(prj$outdmin[[1]])), "g")
})

test_that("init_refmodel() stops for a group term with `||`", {
  expect_error(refm_ids(y ~ x + (x || g)),
               "does not support group terms with `||`", fixed = TRUE)
})

test_that(paste(
  "init_refmodel() warns if a grouping factor has more than one group term,",
  "except if they all have the same ID"
), {
  wrn <- "more than one group term for the grouping factor(s) `g`"
  expect_warning(refm_ids(y ~ x + (1 | g) + (0 + x | g)), wrn, fixed = TRUE)
  expect_warning(
    expect_message(refm_ids(y ~ x + (1 |p| g) + (0 + x | g)), "`p`"),
    wrn, fixed = TRUE
  )
  expect_no_warning(refm_ids(y ~ x + (1 | g) + (0 + x | g),
                             called_from_cvrefbuilder = TRUE))
  expect_no_warning(refm_ids(y ~ x + (1 + x | g)))
  expect_no_warning(
    expect_message(refm_ids(y ~ x + (1 |p| g) + (0 + x |p| g)), "`p`")
  )
})
