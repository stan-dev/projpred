# Group terms with a brms correlation ID (`(1 |p| g)`), group terms with `||`,
# and several group terms of one grouping factor. These tests run on CRAN: the
# reference model is a custom one (no Stan fit).

# Data with 8 groups of 10 observations. The state of the random number
# generator is restored at the end.
dat_ids <- function() {
  withr::local_seed(4)
  g <- factor(paste0("g", rep(1:8, each = 10)))
  x <- rnorm(80)
  b0 <- rnorm(8, sd = 0.8)
  b1 <- rnorm(8, sd = 0.5)
  data.frame(y = 0.5 + 0.7 * x + b0[g] + b1[g] * x + rnorm(80, sd = 0.5),
             x = x, g = g)
}

# A custom reference model: the draws of the linear predictor are the fitted
# values of an lmer() fit plus a small noise.
refm_ids <- function(formula, ...) {
  dat <- dat_ids()
  fit_lmer <- lme4::lmer(y ~ x + (x | g), data = dat)
  withr::local_seed(5)
  ndraws <- 20
  eta <- fitted(fit_lmer) + matrix(rnorm(nrow(dat) * ndraws, sd = 0.05),
                                   nrow(dat), ndraws)
  fit <- list(eta = eta)
  init_refmodel(fit, data = dat, formula = formula, family = gaussian(),
                ref_predfun = function(fit, newdata = NULL) {
                  stopifnot(is.null(newdata))
                  fit$eta
                },
                dis = rep(0.5, ndraws), ...)
}

test_that("rm_group_ids() removes the correlation IDs of group terms", {
  fml <- y ~ x + (1 |p| g) + (0 + x |p| g) + (1 | h) + I(m[, 1])
  expect_message(fml_out <- rm_group_ids(fml), "`p`")
  expect_identical(
    fml_out,
    `environment<-`(y ~ x + (1 | g) + (0 + x | g) + (1 | h) + I(m[, 1]),
                    environment(fml))
  )
  expect_silent(fml_out_q <- rm_group_ids(fml, verbose = FALSE))
  expect_identical(fml_out_q, fml_out)
  # No ID: no change and no message.
  fml_noid <- y ~ x + (1 + x | g)
  expect_silent(expect_identical(rm_group_ids(fml_noid), fml_noid))
})

test_that("extract_terms_response() gives informative errors for `||`", {
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
  # Without the message in the fold-wise reference models:
  expect_no_message(refm_ids(y ~ x + (1 + x |p| g),
                             called_from_cvrefbuilder = TRUE))
  prj <- suppressWarnings(project(refm, predictor_terms = c("x", "(x | g)"),
                                  nclusters = 2, verbose = FALSE))
  expect_s4_class(prj$outdmin[[1]], "lmerMod")
  expect_identical(names(lme4::VarCorr(prj$outdmin[[1]])), "g")
})

test_that("init_refmodel() stops for a group term with `||`", {
  expect_error(refm_ids(y ~ x + (x || g)),
               "does not support group terms with `||`", fixed = TRUE)
})

test_that(paste(
  "init_refmodel() warns if a grouping factor has more than one group term"
), {
  expect_warning(refm_ids(y ~ x + (1 | g) + (0 + x | g)),
                 "more than one group term for the grouping factor(s) `g`",
                 fixed = TRUE)
  expect_no_warning(refm_ids(y ~ x + (1 | g) + (0 + x | g),
                             called_from_cvrefbuilder = TRUE))
  expect_no_warning(refm_ids(y ~ x + (1 + x | g)))
})
