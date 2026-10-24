context("validated LOO: fold projections")

# A Gaussian reference model with fixed draws and more predictors than
# observations, so that submodels with at least n - 1 terms can interpolate
# the reference predictions at the observations used in the projection.
cre_ref_interp <- function(n = 8, p = 9, S = 20) {
  set.seed(5309)
  x <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  dat <- data.frame(y = x[, 1] + rnorm(n), x)
  beta <- matrix(rnorm((p + 1) * S, sd = 0.3), p + 1, S)
  beta[2, ] <- beta[2, ] + 1
  ref_predfun <- function(fit, newdata = NULL) {
    if (is.null(newdata)) newdata <- dat
    cbind(1, as.matrix(newdata[, colnames(x), drop = FALSE])) %*% fit$beta
  }
  ref <- init_refmodel(
    object = list(beta = beta), data = dat,
    formula = reformulate(colnames(x), response = "y"), family = gaussian(),
    ref_predfun = ref_predfun, dis = rep(1, S)
  )
  list(ref = ref, n = n, p = p, S = S)
}

test_that(paste(
  "cv_varsel() with LOO and `validate_search = TRUE` does not use the",
  "left-out observation in the projections of its fold"
), {
  skip_on_cran()
  r <- cre_ref_interp()
  cvvs <- suppressWarnings(cv_varsel(
    r$ref, cv_method = "LOO", validate_search = TRUE, method = "forward",
    nterms_max = r$p, nclusters = NULL, ndraws = 5, ndraws_pred = r$S,
    refit_prj = TRUE, seed = 1, verbose = FALSE
  ))
  lppd <- sapply(cvvs$summaries$sub, `[[`, "lppd")
  expect_equal(dim(lppd), c(r$n, r$p + 1))
  # Columns `r$n` to `r$p + 1` are the sizes n - 1 to p. If the left-out
  # observation were used in the fold projections, these submodels would
  # interpolate the reference predictions at it, and their fold lppd would be
  # the same for all these sizes.
  sizes_interp <- seq(r$n, r$p + 1L)
  spread <- apply(lppd[, sizes_interp, drop = FALSE], 1, function(z) diff(range(z)))
  expect_gt(max(spread), 1e-2)
})
