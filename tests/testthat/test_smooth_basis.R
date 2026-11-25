context("smooth_basis")

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
ndr_sb <- 25L
eta_draws_sb <- sapply(seq_len(ndr_sb), function(s) {
  eta_sb + 0.05 * sin(s * idx_sb + s)
})

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
  expect_true(all(c("s(x).1", "t2(z,w).1") %in% colnames(m_r)))
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
  expect_true(all(is.na(m[, "sds_sx_1"])))
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

test_that("brms_rename() renames as brms does", {
  expect_identical(brms_rename(c("s(x)", "t2(x,z)", "s(x.1)", "I(x^2)")),
                   c("sx", "t2xz", "sx.1", "IxE2"))
})
