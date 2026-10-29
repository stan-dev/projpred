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
  # The augmented-data projection supports only these families:
  is_projectable <- function(tstsetup) {
    fam <- refmods[[tstsetup]]$family
    !fam$for_augdat ||
      fam$family %in% c("cumulative", "cumulative_rstanarm", "categorical")
  }
  keys <- vapply(names(refmods), case_key, character(1))
  tstsetups <- names(refmods)[!duplicated(keys)]
  tstsetups <- Filter(is_projectable, tstsetups)
  expect_gt(length(tstsetups), 5)
  for (tstsetup in tstsetups) {
    refm <- refmods[[tstsetup]]
    trms <- terms_case(refm)
    p_ref <- get_refdist(refm, ndraws = 3)
    p_ref_chg <- change_target(p_ref, refm, i = 1)
    refm_w0 <- refm
    refm_w0$wobs[1] <- 0
    # The same with other predictor values of observation 1 (in the training
    # data of the projection only):
    refm_w0_x <- refm_w0
    refm_w0_x$fetch_data <- function(obs = NULL) {
      dat <- refm_w0$fetch_data(obs = obs)
      if (is.null(obs)) {
        trms_ref <- delete.response(terms(refm$formula))
        offs_vars <- unlist(lapply(attr(trms_ref, "offset"), function(i) {
          all.vars(attr(trms_ref, "variables")[[i + 1L]])
        }))
        x_nms <- setdiff(intersect(all.vars(trms_ref), names(dat)), offs_vars)
        x_nms <- x_nms[vapply(dat[x_nms], is.numeric, logical(1))]
        dat[1, x_nms] <- dat[1, x_nms] + 3
      }
      dat
    }
    is_mlvl_add <- any(grepl("\\||^s\\(", trms))
    # The predictions of the submodel fits for the training data:
    prd <- function(p, r) {
      outdmin <- suppressWarnings(suppressMessages(proj_to_submodl(
        predictor_terms = trms, p_ref = p, refmodel = r, verbose = FALSE
      )$outdmin))
      unclass(refm$mu_fun(outdmin, offset = refm$offset))
    }
    # An error in one setup is a failure for that setup only:
    res <- tryCatch(
      list(w0_chg = prd(p_ref_chg, refm_w0), w0 = prd(p_ref, refm_w0),
           w1_chg = prd(p_ref_chg, refm), w1 = prd(p_ref, refm),
           w0_x = if (is_mlvl_add) prd(p_ref, refm_w0_x)),
      error = function(e) conditionMessage(e)
    )
    if (is.character(res)) {
      expect_true(FALSE, label = paste(tstsetup, "error:", res))
      next
    }
    # With weight 0, the changed target of observation 1 does not matter:
    expect_equal(res$w0_chg, res$w0, tolerance = 1e-5, info = tstsetup)
    # The multilevel and additive fitters get the data without the rows with
    # weight 0, so the predictor values of observation 1 must not matter
    # either:
    if (is_mlvl_add) {
      expect_equal(res$w0_x, res$w0, tolerance = 1e-5, info = tstsetup)
    }
    # With its original weight, it does:
    expect_gt(max(abs(res$w1_chg - res$w1)), 1e-4, label = tstsetup)
  }
})

test_that(paste(
  "the submodel fits predict a factor level that only an observation with",
  "weight 0 has as the reference level"
), {
  skip_on_cran()
  set.seed(2)
  n <- 40
  dat <- data.frame(x = rnorm(n), g = factor(rep(paste0("g", 1:5), 8)),
                    xc = factor(c(rep(c("a", "b"), length.out = n - 1), "r")))
  for (s in 1:2) dat[[paste0(".y.", s)]] <- dat$x + rnorm(n)
  # Observation n is the only one with level "r", and it has weight 0:
  wobs <- c(rep(1, n - 1), 0)
  dat_ref <- dat
  dat_ref$xc[n] <- "a"
  dat_ref$xc <- droplevels(dat_ref$xc)
  divmin_w0 <- function(fml) {
    suppressWarnings(divmin(
      formula = fml, data = dat, family = extend_family(gaussian()),
      weights = wobs, projpred_var = matrix(NA, n, 2),
      projpred_ws_aug = matrix(NA, n, 2), verbose_divmin = FALSE
    ))
  }
  fmls <- list(glmm = cbind(.y.1, .y.2) ~ x + xc + (1 | g),
               gam = cbind(.y.1, .y.2) ~ s(x) + xc)
  for (nm in names(fmls)) {
    fits <- divmin_w0(fmls[[nm]])
    expect_identical(attr(fits, "projpred_lost_lvls"),
                     list(xc = list(lost = "r", ref = "a")), info = nm)
    prd <- subprd(fits, recode_lost_levels(
      dat, attr(fits, "projpred_lost_lvls")
    ))
    expect_identical(prd, subprd(fits, dat_ref), info = nm)
  }
  # Augmented-data projection (rows: the observations within the response
  # categories):
  ncat <- 3
  lv <- paste0("c", seq_len(ncat))
  dat_aug <- dat[rep(seq_len(n), ncat), ]
  ws_aug <- matrix(runif(n * ncat * 2), ncol = 2)
  for (fam_nm in c("cumulative", "categorical")) {
    if (fam_nm == "cumulative") {
      skip_if_not_installed("ordinal")
    } else {
      skip_if_not_installed("mclogit")
    }
    dat_aug$y <- factor(rep(lv, each = n), levels = lv,
                        ordered = fam_nm == "cumulative")
    fam <- structure(list(family = fam_nm, link = "logit"), class = "family")
    # With and without group-level terms (clmm() or polr(), mblogit() or
    # multinom()):
    for (fml in list(y ~ x + xc + (1 | g), y ~ x + xc)) {
      info_txt <- paste(fam_nm, deparse(fml))
      fits <- suppressWarnings(divmin_augdat(
        formula = fml, data = dat_aug, family = fam,
        weights = wobs, projpred_var = matrix(NA, n * ncat, 2),
        projpred_ws_aug = ws_aug, verbose_divmin = FALSE
      ))
      expect_identical(attr(fits, "projpred_lost_lvls"),
                       list(xc = list(lost = "r", ref = "a")), info = info_txt)
      prd <- subprd_augdat(fits, recode_lost_levels(
        dat, attr(fits, "projpred_lost_lvls")
      ))
      expect_identical(prd, subprd_augdat(fits, dat_ref), info = info_txt)
    }
  }
  # Other factor-like predictors: a lost first level, a character variable and
  # a numeric variable inside factor():
  dat2 <- dat
  dat2$xf <- factor(c(rep(c("a", "b"), length.out = n - 1), "r"),
                    levels = c("r", "a", "b"))
  dat2$xs <- c(rep(c("u", "v"), length.out = n - 1), "w")
  dat2$z <- c(rep(1:3, length.out = n - 1), 9)
  fits2 <- suppressWarnings(divmin(
    formula = cbind(.y.1, .y.2) ~ x + xf + xs + factor(z) + (1 | g),
    data = dat2, family = extend_family(gaussian()), weights = wobs,
    projpred_var = matrix(NA, n, 2), projpred_ws_aug = matrix(NA, n, 2),
    verbose_divmin = FALSE
  ))
  expect_identical(attr(fits2, "projpred_lost_lvls"), list(
    xf = list(lost = "r", ref = "a"), xs = list(lost = "w", ref = "u"),
    z = list(lost = 9, ref = 1)
  ))
  dat2_ref <- dat2
  dat2_ref$xf[n] <- "a"
  dat2_ref$xf <- droplevels(dat2_ref$xf)
  dat2_ref$xs[n] <- "u"
  dat2_ref$z[n] <- 1
  expect_identical(
    subprd(fits2, recode_lost_levels(dat2, attr(fits2, "projpred_lost_lvls"))),
    subprd(fits2, dat2_ref)
  )
  # A grouping factor is not recoded, also if it is a fixed-effect predictor:
  # its lost level is a new group.
  dat3 <- dat
  dat3$g <- factor(c(as.character(dat$g[-n]), "g9"))
  expect_null(lost_levels(dat3, keep = wobs != 0,
                          formula = .y.1 ~ x + g + (1 | g)))
  # Through proj_to_submodl() and the `mu_fun()` of a reference model:
  n_drws <- 4
  dat_r <- dat[, c("x", "g", "xc")]
  dat_r$y <- dat$.y.1
  refm <- suppressWarnings(init_refmodel(
    object = list(beta = rnorm(n_drws, 1, 0.1)), data = dat_r,
    formula = y ~ x + xc + (1 | g), family = gaussian(),
    ref_predfun = function(fit, newdata = NULL) {
      outer((newdata %||% dat_r)$x, fit$beta)
    },
    dis = rep(1, n_drws)
  ))
  refm$wobs[n] <- 0
  submodl <- suppressWarnings(proj_to_submodl(
    predictor_terms = c("x", "xc", "(1 | g)"),
    p_ref = get_refdist(refm, ndraws = 2), refmodel = refm, verbose = FALSE
  ))
  expect_identical(attr(submodl$outdmin, "projpred_lost_lvls"),
                   list(xc = list(lost = "r", ref = "a")))
  dat_r_ref <- dat_r
  dat_r_ref$xc[n] <- "a"
  dat_r_ref$xc <- droplevels(dat_r_ref$xc)
  expect_equal(refm$mu_fun(submodl$outdmin),
               refm$mu_fun(submodl$outdmin, newdata = dat_r_ref))
  # No level is lost if observation n has a positive weight:
  wobs[n] <- 1
  expect_null(attr(divmin_w0(fmls$glmm), "projpred_lost_lvls"))
})
