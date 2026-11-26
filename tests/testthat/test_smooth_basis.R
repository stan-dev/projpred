# Projected draws of smooth terms in the basis of the reference model fit (see
# `R/smooth_basis.R`). The reference models are custom reference models without
# Stan: an empty `brmsfit` (brms::brm(empty = TRUE)) and an object of class
# `stanreg` with the elements `jam` and `x` that rstanarm::stan_gamm4() makes
# (from mgcv::jagam(diagonalize = TRUE), as in stan_gamm4()).

n_sb <- 150L
idx_sb <- seq_len(n_sb)
dat_sb <- data.frame(x = (sin(1.7 * idx_sb) + 1) / 2,
                     z = (cos(2.3 * idx_sb) + 1) / 2,
                     w = (sin(0.9 * idx_sb + 1) + 1) / 2,
                     v = cos(0.4 * idx_sb),
                     wt = 1 + 0.5 * sin(3.1 * idx_sb),
                     tr = 5 + idx_sb %% 4)
eta_sb <- sin(2 * pi * dat_sb$x) + dat_sb$z * dat_sb$w + 0.5 * dat_sb$v
dat_sb$y <- eta_sb + 0.3 * sin(5.1 * idx_sb)
dat_sb$yb <- round(dat_sb$tr * plogis(eta_sb))
dat_sb$g <- factor(paste0("g", idx_sb %% 6))
ndr_sb <- 25L
eta_draws_sb <- sapply(seq_len(ndr_sb), function(s) {
  eta_sb + 0.05 * sin(s * idx_sb + s)
})
# With group effects (for the GAMM submodels):
eff_g_sb <- c(-0.6, -0.3, 0, 0.2, 0.4, 0.5)
eta_draws_g_sb <- eta_draws_sb + eff_g_sb[as.integer(dat_sb$g)]

# The warnings of init_refmodel() for these custom reference models:
muffle_refmodel_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("cvrefbuilder|experimental", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

# A projection with the custom reference model of `fit`. `wcol`: the column of
# the observation weights (or `NULL`). `data`: the reference model data.
prj_sb <- function(fit, fml, predictor_terms, family = gaussian(),
                   eta_draws = eta_draws_sb, wcol = NULL, data = dat_sb,
                   ...) {
  resp <- as.character(fml[[2]])
  refm <- muffle_refmodel_warnings(init_refmodel(
    fit, data = data, formula = fml, family = family,
    ref_predfun = function(fit, newdata = NULL) {
      if (!is.null(newdata)) stop("Unexpected `newdata`.")
      eta_draws
    },
    extract_model_data = function(object, newdata, wrhs = NULL, orhs = NULL,
                                  extract_y = TRUE) {
      list(y = newdata[[resp]],
           weights = if (is.null(wcol)) rep(1, nrow(newdata)) else
             newdata[[wcol]],
           offset = rep(0, nrow(newdata)))
    },
    dis = if (family$family == "gaussian") rep(0.3, ncol(eta_draws))
  ))
  project(refm, predictor_terms = predictor_terms, verbose = FALSE, ...)
}

brmsfit_sb <- function(fml) {
  suppressMessages(brms::brm(fml, data = dat_sb, empty = TRUE))
}

stanreg_sb <- function(fml) {
  dat_j <- dat_sb
  dat_j$temp_y <- 1
  jd <- mgcv::jagam(update(fml, temp_y ~ .), family = gaussian(),
                    data = dat_j, file = tempfile(fileext = ".jags"),
                    diagonalize = TRUE)
  X <- jd$jags.data$X
  colnames(X) <- jd$pregam$term.names
  structure(list(jam = jd$pregam, x = X), class = "stanreg")
}

# The brms basis (`Xs`, `Zs_*`) with the brms parameter names as column names.
brms_basis_sb <- function(fit) {
  sdat <- brms::standata(fit)
  X <- sdat$Xs
  colnames(X) <- paste0("bs_", colnames(X))
  for (nm in grep("^Zs_", names(sdat), value = TRUE)) {
    jk <- as.integer(strsplit(sub("^Zs_", "", nm), "_")[[1]])
    lab <- sub("_1$", "",
               colnames(sdat$Xs)[attr(sdat$Xs, "smcols")[[jk[1]]][1]])
    Z <- sdat[[nm]]
    colnames(Z) <- paste0("s_", lab, "_", jk[2], "[", seq_len(ncol(Z)), "]")
    X <- cbind(X, Z)
  }
  return(X)
}

# The SD blocks of a brms basis from brms_basis_sb(): for each SD column, the
# coefficient columns.
brms_blocks_sb <- function(X_sm) {
  s_nms <- grep("^s_", colnames(X_sm), value = TRUE)
  sd_of <- sub("^s_", "sds_", sub("\\[[0-9]+\\]$", "", s_nms))
  split(s_nms, sd_of)
}

# Checks that the columns of the smooth terms in `m` (from as.matrix()) and the
# basis `X_sm` (columns named as the parameters) give the linear predictor of
# the projection, with `X_par` (the intercept and parametric columns):
expect_lp_rebuilt <- function(m, prj, X_par, X_sm, info, ...) {
  X <- cbind(X_par, X_sm[, intersect(colnames(X_sm), colnames(m)),
                         drop = FALSE])
  lp <- m[, colnames(X), drop = FALSE] %*% t(X)
  lp_prj <- proj_linpred(prj, transform = FALSE, integrated = FALSE,
                         ...)$pred
  # (A clustered projection gives the weights of the draws as an attribute.)
  expect_equal(matrix(lp, nrow(lp)), matrix(lp_prj, nrow(lp_prj)),
               tolerance = 1e-8, info = info)
}

# Checks the SD columns: with the penalty sig2 / SD^2 on the coefficients of
# each SD block (as in the mixed-model form of the smooth terms), the penalized
# weighted least-squares fit to the working response of each submodel fit (at
# the working weights of the fit) gives the projected coefficients. `blocks`:
# for each SD column, the coefficient columns.
expect_sds_reproduce_fit <- function(m, prj, X_par, X_sm, blocks, info) {
  X <- cbind(X_par, X_sm)
  for (s in seq_along(prj$outdmin)) {
    fit_s <- prj$outdmin[[s]]
    pen <- setNames(rep(0, ncol(X)), colnames(X))
    for (sd_nm in names(blocks)) {
      pen[blocks[[sd_nm]]] <- fit_s$sig2 / m[s, sd_nm]^2
    }
    eta_s <- fit_s$linear.predictors
    z_s <- eta_s + (fit_s$y - fit_s$fitted.values) /
      fit_s$family$mu.eta(eta_s)
    wX <- X * fit_s$weights
    theta <- solve(crossprod(wX, X) + diag(pen), crossprod(wX, z_s))
    expect_equal(drop(theta), m[s, colnames(X)], tolerance = 1e-6,
                 info = paste(info, "draw", s))
  }
}

test_that("as.matrix.projection() gives smooth terms in the brms basis", {
  skip_if_not_installed("brms")
  fml <- y ~ s(x) + t2(z, w) + v
  fit <- brmsfit_sb(fml)
  # With observation weights:
  prj <- prj_sb(fit, fml, c("s(x)", "t2(z, w)", "v"), wcol = "wt",
                ndraws = ndr_sb)
  expect_warning(m <- as.matrix(prj), NA)
  X_sm <- brms_basis_sb(fit)
  sd_nms <- c("sds_sx_1", paste0("sds_t2zw_", 1:3))
  expect_setequal(colnames(m), c("b_Intercept", "b_v", colnames(X_sm), sd_nms,
                                 "sigma"))
  X_par <- cbind(b_Intercept = 1, b_v = dat_sb$v)
  expect_lp_rebuilt(m, prj, X_par, X_sm, info = "brms")
  blocks <- brms_blocks_sb(X_sm)
  expect_setequal(names(blocks), sd_nms)
  expect_sds_reproduce_fit(m, prj, X_par, X_sm, blocks, info = "brms")

  # A submodel with a subset of the smooth terms:
  prj_x <- prj_sb(fit, fml, "s(x)", ndraws = ndr_sb)
  m_x <- as.matrix(prj_x)
  expect_setequal(colnames(m_x),
                  c("b_Intercept", "bs_sx_1", "sds_sx_1",
                    grep("^s_sx_1\\[", colnames(X_sm), value = TRUE),
                    "sigma"))
  expect_lp_rebuilt(m_x, prj_x, cbind(b_Intercept = rep(1, n_sb)), X_sm,
                    info = "brms, s(x)")

  # Clustered projection:
  prj_cl <- prj_sb(fit, fml, c("s(x)", "t2(z, w)", "v"), nclusters = 3)
  m_cl <- as.matrix(prj_cl, allow_nonconst_wdraws_prj = TRUE)
  expect_identical(colnames(m_cl), colnames(m))
  expect_lp_rebuilt(m_cl, prj_cl, X_par, X_sm, info = "brms, clusters",
                    allow_nonconst_wdraws_prj = TRUE)

  # With `nm_scheme = "rstanarm"`, a brmsfit gives no basis:
  expect_warning(m_r <- as.matrix(prj, nm_scheme = "rstanarm"),
                 "does not give the bases of the smooth terms")
  expect_all_true_tst(c("s(x).1", "t2(z,w).1") %in% colnames(m_r))
})

test_that(paste(
  "as.matrix.projection() gives smooth terms in the rstanarm basis, without",
  "the SD of the null space"
), {
  fml <- y ~ s(x) + t2(z, w) + v
  fit <- stanreg_sb(fml)
  sd_nms <- paste0("smooth_sd[", c("s(x)1", "t2(z,w)rr", "t2(z,w)nr",
                                   "t2(z,w)rn"), "]")
  X <- fit$x
  X_par <- X[, c("(Intercept)", "v")]
  X_sm <- X[, setdiff(colnames(X), colnames(X_par))]
  # The blocks of rstanarm::stan_gamm4(): by the ranks of the penalties.
  blocks <- list()
  sp_nms <- names(fit$jam$sp)
  sp_i <- 0L
  for (sm in fit$jam$smooth) {
    col_i <- sm$first.para - 1L
    for (r in sm$rank) {
      sp_i <- sp_i + 1L
      blocks[[paste0("smooth_sd[", sp_nms[sp_i], "]")]] <-
        colnames(X)[col_i + seq_len(r)]
      col_i <- col_i + r
    }
  }
  blocks <- blocks[sd_nms]
  for (fam_nm in c("gaussian", "binomial")) {
    if (fam_nm == "gaussian") {
      prj <- prj_sb(fit, fml, c("s(x)", "t2(z, w)", "v"), ndraws = ndr_sb)
    } else {
      fml_b <- yb ~ s(x) + t2(z, w) + v
      prj <- prj_sb(fit, fml_b, c("s(x)", "t2(z, w)", "v"),
                    family = binomial(), wcol = "tr", ndraws = ndr_sb)
    }
    expect_warning(m <- as.matrix(prj), NA, info = fam_nm)
    expect_setequal(colnames(m),
                    c(colnames(X), sd_nms,
                      if (fam_nm == "gaussian") "sigma"))
    expect_lp_rebuilt(m, prj, X_par, X_sm, info = fam_nm)
    expect_sds_reproduce_fit(m, prj, X_par, X_sm, blocks, info = fam_nm)
  }
})

test_that(paste(
  "as.matrix.projection() leaves out the SD columns of a smooth term that",
  "overlaps another one"
), {
  skip_if_not_installed("brms")
  fml <- y ~ s(x) + t2(x, z)
  fit <- brmsfit_sb(fml)
  prj <- prj_sb(fit, fml, c("s(x)", "t2(x, z)"), ndraws = ndr_sb)
  expect_warning(m <- as.matrix(prj), "fewer basis functions")
  expect_true("sds_sx_1" %in% colnames(m))
  expect_false(any(grepl("^sds_t2xz", colnames(m))))
  expect_lp_rebuilt(m, prj, cbind(b_Intercept = rep(1, n_sb)),
                    brms_basis_sb(fit), info = "overlap")
})

test_that(paste(
  "as.matrix.projection() gives `NA` SDs if the submodel reproduces the",
  "reference model exactly"
), {
  skip_if_not_installed("brms")
  fml <- y ~ s(x)
  fit <- brmsfit_sb(fml)
  X_sm <- brms_basis_sb(fit)
  # Reference draws in the span of the submodel basis:
  eta_span <- sapply(seq_len(ndr_sb), function(s) {
    drop(cbind(1, X_sm) %*% (0.3 * sin(s * seq_len(ncol(X_sm) + 1L))))
  })
  prj <- prj_sb(fit, fml, "s(x)", eta_draws = eta_span, ndraws = ndr_sb)
  expect_warning(m <- as.matrix(prj),
                 paste("not identified in", ndr_sb, "of", ndr_sb))
  expect_all_true_tst(is.na(m[, "sds_sx_1"]))
  expect_lp_rebuilt(m, prj, cbind(b_Intercept = rep(1, n_sb)), X_sm,
                    info = "span")
})

test_that(paste(
  "as.matrix.projection() keeps the mgcv basis with a warning if the bases do",
  "not match"
), {
  skip_if_not_installed("brms")
  fml <- y ~ s(x)
  fit <- brmsfit_sb(fml)
  # The reference model data in another order than the data of the fit:
  dat_rev <- dat_sb[rev(idx_sb), ]
  prj <- prj_sb(fit, fml, "s(x)", data = dat_rev,
                eta_draws = eta_draws_sb[rev(idx_sb), ], ndraws = ndr_sb)
  expect_warning(m <- as.matrix(prj), "does not match a smooth term")
  expect_true("b_s(x).1" %in% colnames(m))
})

test_that(paste(
  "as.matrix.projection() keeps the mgcv basis without a warning for a",
  "custom reference model fit and the default `nm_scheme`"
), {
  fml <- y ~ s(x) + v
  prj <- prj_sb(structure(list(), class = "custom_fit_sb"), fml,
                c("s(x)", "v"), ndraws = ndr_sb)
  expect_warning(m <- as.matrix(prj), NA)
  expect_identical(colnames(m),
                   c(names(coef(prj$outdmin[[1]])), "sigma"))
  expect_warning(as.matrix(prj, nm_scheme = "brms"),
                 "does not give the bases of the smooth terms")
})

test_that(paste(
  "as.matrix.projection() gives the smooth and the group-level terms of GAMM",
  "submodels in the brms basis"
), {
  skip_if_not_installed("brms")
  skip_if_not_installed("gamm4")
  # The known warnings of the multilevel projections and of lme4:
  prj_gamm <- function(...) {
    withCallingHandlers(prj_sb(...), warning = function(w) {
      if (grepl(paste("issue #323|can be near zero|might not have converged",
                      "failed to converge", sep = "|"),
                conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    })
  }
  fml <- y ~ s(x) + v + (1 | g)
  fit <- brmsfit_sb(fml)
  X_sm <- brms_basis_sb(fit)
  r_nms <- paste0("r_g[", levels(dat_sb$g), ",Intercept]")
  Z_g <- outer(as.integer(dat_sb$g), seq_len(nlevels(dat_sb$g)), "==") * 1
  colnames(Z_g) <- r_nms
  X_par <- cbind(b_Intercept = 1, b_v = dat_sb$v, Z_g)
  for (fam_nm in c("gaussian", "binomial")) {
    if (fam_nm == "gaussian") {
      prj <- prj_gamm(fit, fml, c("s(x)", "v", "(1 | g)"),
                    eta_draws = eta_draws_g_sb, ndraws = ndr_sb)
    } else {
      prj <- prj_gamm(fit, update(fml, yb ~ .), c("s(x)", "v", "(1 | g)"),
                    family = binomial(), wcol = "tr",
                    eta_draws = eta_draws_g_sb, ndraws = ndr_sb)
    }
    expect_s3_class(prj$outdmin[[1]], "gamm4")
    expect_warning(m <- as.matrix(prj), NA, info = fam_nm)
    expect_setequal(colnames(m),
                    c("b_Intercept", "b_v", colnames(X_sm), "sds_sx_1",
                      "sd_g__Intercept", r_nms,
                      if (fam_nm == "gaussian") "sigma"))
    expect_lp_rebuilt(m, prj, X_par, X_sm, info = fam_nm)
    # The group SD and the smooth SD are those of the lme4 part of the fits
    # (the smooth term is the variance component `Xr` there):
    sds_mer <- t(vapply(prj$outdmin, function(fit_s) {
      vc <- lme4::VarCorr(fit_s$mer)
      c(attr(vc$g, "stddev"), attr(vc$Xr, "stddev"))
    }, numeric(2)))
    expect_equal(unname(m[, "sd_g__Intercept"]), unname(sds_mer[, 1]),
                 tolerance = 1e-10, info = fam_nm)
    expect_equal(unname(m[, "sds_sx_1"]), unname(sds_mer[, 2]),
                 tolerance = 1e-6, info = fam_nm)
    # as_draws_matrix() gives the same draws:
    expect_identical(unclass(posterior::as_draws_matrix(prj))[, colnames(m)],
                     unclass(posterior::as_draws_matrix(m))[, colnames(m)],
                     info = fam_nm)
  }

  # A t2() term (several variance components in the lme4 part):
  fml_t2 <- y ~ s(x) + t2(z, w) + v + (1 | g)
  fit_t2 <- brmsfit_sb(fml_t2)
  prj_t2 <- prj_gamm(fit_t2, fml_t2, c("s(x)", "t2(z, w)", "v", "(1 | g)"),
                   eta_draws = eta_draws_g_sb, ndraws = ndr_sb)
  expect_warning(m_t2 <- as.matrix(prj_t2), NA)
  X_sm_t2 <- brms_basis_sb(fit_t2)
  expect_setequal(colnames(m_t2),
                  c("b_Intercept", "b_v", colnames(X_sm_t2), "sds_sx_1",
                    paste0("sds_t2zw_", 1:3), "sd_g__Intercept", r_nms,
                    "sigma"))
  expect_lp_rebuilt(m_t2, prj_t2, X_par, X_sm_t2, info = "t2")
  # The SDs of each smooth term are those of its variance components in the
  # lme4 part (`lmer.name`; compared as sets, because the order of the
  # penalties of a t2() term can differ):
  for (s in seq_along(prj_t2$outdmin)) {
    fit_s <- prj_t2$outdmin[[s]]
    vc <- lme4::VarCorr(fit_s$mer)
    sm <- fit_s$gam$smooth
    sd_lme4 <- lapply(sm, function(sm_j) {
      vapply(sm_j$lmer.name, function(nm) attr(vc[[nm]], "stddev"),
             numeric(1))
    })
    expect_equal(unname(m_t2[s, "sds_sx_1"]), unname(sd_lme4[[1]]),
                 tolerance = 1e-6, info = paste("t2, draw", s))
    expect_equal(sort(unname(m_t2[s, paste0("sds_t2zw_", 1:3)])),
                 sort(unname(sd_lme4[[2]])), tolerance = 1e-6,
                 info = paste("t2, draw", s))
  }

  # rstanarm naming (a stan_gamm4() fit without Stan; the group-level term
  # is only in the reference model formula of the projection):
  fit_r <- stanreg_sb(y ~ s(x) + v)
  prj_r <- prj_gamm(fit_r, fml, c("s(x)", "v", "(1 | g)"),
                    eta_draws = eta_draws_g_sb, ndraws = ndr_sb)
  expect_warning(m_r <- as.matrix(prj_r), NA)
  b_nms <- paste0("b[(Intercept) g:", levels(dat_sb$g), "]")
  expect_setequal(colnames(m_r),
                  c(colnames(fit_r$x), "smooth_sd[s(x)1]",
                    "Sigma[g:(Intercept),(Intercept)]", b_nms, "sigma"))
  vc_r <- t(vapply(prj_r$outdmin, function(fit_s) {
    vc <- lme4::VarCorr(fit_s$mer)
    c(vc$g[1, 1], attr(vc$Xr, "stddev"))
  }, numeric(2)))
  # The variance (not the SD) of the group-level intercepts:
  expect_equal(unname(m_r[, "Sigma[g:(Intercept),(Intercept)]"]),
               unname(vc_r[, 1]), tolerance = 1e-10)
  expect_equal(unname(m_r[, "smooth_sd[s(x)1]"]), unname(vc_r[, 2]),
               tolerance = 1e-6)
  Z_g_r <- Z_g
  colnames(Z_g_r) <- b_nms
  X_r <- fit_r$x
  expect_lp_rebuilt(m_r, prj_r, cbind(X_r[, c("(Intercept)", "v")], Z_g_r),
                    X_r[, setdiff(colnames(X_r), c("(Intercept)", "v"))],
                    info = "rstanarm")

  # Without group-level effects in the predictions (only the group SDs):
  opts_old <- options(projpred.mlvl_pred_new = TRUE)
  on.exit(options(opts_old), add = TRUE)
  m_new <- as.matrix(prj)
  expect_false(any(grepl("^r_", colnames(m_new))))
  expect_true("sd_g__Intercept" %in% colnames(m_new))
  m_new_r <- as.matrix(prj_r)
  expect_false(any(grepl("^b\\[", colnames(m_new_r))))
  expect_true("Sigma[g:(Intercept),(Intercept)]" %in% colnames(m_new_r))
})

test_that(paste(
  "as.matrix.projection() maps the smooth terms if the reference model data",
  "is a model frame"
), {
  # (As the data of a rstanarm::stan_gamm4() fit with argument `random`.)
  dat_mf <- model.frame(y ~ x + v, data = dat_sb)
  fml <- y ~ s(x) + v
  fit <- stanreg_sb(fml)
  prj <- prj_sb(fit, fml, c("s(x)", "v"), data = dat_mf, ndraws = ndr_sb)
  expect_false(is.null(attr(prj$refmodel$fetch_data(), "terms")))
  expect_warning(m <- as.matrix(prj), NA)
  expect_true("smooth_sd[s(x)1]" %in% colnames(m))
  X <- fit$x
  expect_lp_rebuilt(m, prj, X[, c("(Intercept)", "v")],
                    X[, setdiff(colnames(X), c("(Intercept)", "v"))],
                    info = "model frame")
})

test_that("brms_rename() renames as brms does", {
  expect_identical(brms_rename(c("s(x)", "t2(x,z)", "s(x.1)", "I(x^2)")),
                   c("sx", "t2xz", "sx.1", "IxE2"))
})
