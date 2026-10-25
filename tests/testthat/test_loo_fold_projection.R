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

test_that(paste(
  "an observation with weight 0 does not affect the projection (all model",
  "types, families and projection types of the test setups)"
), {
  skip_if_not(run_prj)
  # Changes the reference target of observation `i` in a way that changes the
  # projection if the observation has a positive weight:
  change_target <- function(p_ref, refm, i) {
    fam <- refm$family
    if (fam$for_augdat) {
      # Reverse the order of the category probabilities of observation `i`:
      n <- refm$nobs
      idx <- i + n * (seq_len(nrow(p_ref$mu) / n) - 1L)
      p_ref$mu[idx, ] <- p_ref$mu[rev(idx), ]
    } else if (fam$for_latent || fam$family %in% c("gaussian", "Student_t")) {
      p_ref$mu[i, ] <- p_ref$mu[i, ] + 5
    } else if (fam$family == "binomial") {
      p_ref$mu[i, ] <- 1 - p_ref$mu[i, ]
    } else {
      p_ref$mu[i, ] <- 3 * p_ref$mu[i, ] + 1
    }
    p_ref
  }
  # One linear term and, depending on the model type, one group-level term
  # and one smooth term, so that each submodel fitter is used:
  terms_case <- function(refm) {
    trms <- setdiff(split_formula(refm$formula, data = refm$fetch_data()), "1")
    trms <- trms[!grepl("\\+", trms)]
    is_grp <- grepl("\\|", trms)
    is_smooth <- grepl("^s\\(", trms)
    c(trms[!is_grp & !is_smooth][1], trms[is_grp][1], trms[is_smooth][1])[
      c(TRUE, any(is_grp), any(is_smooth))
    ]
  }
  # One test setup per combination of model type, family, projection type and
  # observation weights:
  case_key <- function(tstsetup) {
    refm <- refmods[[tstsetup]]
    mod_nm <- sub("^[^.]+\\.([^.]+)\\..*$", "\\1", tstsetup)
    paste(mod_nm, refm$family$family, refm$family$for_latent,
          refm$family$for_augdat, any(refm$wobs != 1))
  }
  keys <- vapply(names(refmods), case_key, character(1))
  tstsetups <- names(refmods)[!duplicated(keys)]
  expect_gt(length(tstsetups), 5)
  for (tstsetup in tstsetups) {
    refm <- refmods[[tstsetup]]
    trms <- terms_case(refm)
    p_ref <- get_refdist(refm, ndraws = 3)
    p_ref_chg <- change_target(p_ref, refm, i = 1)
    refm_w0 <- refm
    refm_w0$wobs[1] <- 0
    # The predictions of the submodel fits for the training data:
    prd <- function(p, r) {
      outdmin <- suppressWarnings(suppressMessages(proj_to_submodl(
        predictor_terms = trms, p_ref = p, refmodel = r, verbose = FALSE
      )$outdmin))
      unclass(refm$mu_fun(outdmin, offset = refm$offset))
    }
    # With weight 0, the changed target of observation 1 does not matter:
    expect_equal(prd(p_ref_chg, refm_w0), prd(p_ref, refm_w0),
                 tolerance = 1e-5, info = tstsetup)
    # With its original weight, it does:
    expect_gt(max(abs(prd(p_ref_chg, refm) - prd(p_ref, refm))), 1e-4,
              label = tstsetup)
  }
})
