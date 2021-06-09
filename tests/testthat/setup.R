#_______________________________________________________________________________
# Setup for the unit tests
#_______________________________________________________________________________

# General setup -----------------------------------------------------------

seed_tst <- 1235
set.seed(seed_tst)
source(testthat::test_path("helpers", "SW.R"))
fam_nms <- setNames(nm = c("gauss", "binom"))

# rstanarm setup ----------------------------------------------------------

chains_tst <- 2L
iter_tst <- 500L

# projpred setup ----------------------------------------------------------

nclusters_tst <- 2L
nclusters_pred_tst <- 3L
nresample_clusters_tst <- 100L
nresample_clusters_default <- 1000L # Adopt this if the default is changed.
seed2_tst <- 866028
solterms_tst <- c("x.2", "x.4")

# GLMs --------------------------------------------------------------------

## Data -------------------------------------------------------------------

n_obs <- 40L
nterms_glm <- 5L
x_pop <- matrix(rnorm(n_obs * nterms_glm, 0, 1), n_obs, nterms_glm)
b_pop <- runif(nterms_glm) - 0.5
icpt <- -0.42
offs <- rnorm(n_obs)
eta_glm <- icpt +
  x_pop %*% b_pop +
  offs
disp <- runif(1L, 1, 2)
w_obs <- sample(1:4, n_obs, replace = TRUE)
f_gauss <- gaussian()
f_binom <- binomial()
df_glm <- data.frame(y_gauss = rnorm(n_obs, f_gauss$linkinv(eta_glm), disp),
                     y_binom = rbinom(n_obs, w_obs, f_binom$linkinv(eta_glm)),
                     x = x_pop, w_obs_col = w_obs, offs_col = offs)
ys_glm <- lapply(fam_nms, function(fam_nm) {
  df_glm[[paste0("y_", fam_nm)]]
})

## Fit --------------------------------------------------------------------

# Notes:
#   * Argument `weights` is not needed when using the cbind() syntax (for the
#     binomial family with > 1 trials).
SW({
  fit_gauss_glm <- rstanarm::stan_glm(
    y_gauss ~ x.1 + x.2 + x.3 + x.4 + x.5,
    family = f_gauss, data = df_glm,
    weights = w_obs, offset = offs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_binom_glm <- rstanarm::stan_glm(
    cbind(y_binom, w_obs_col - y_binom) ~ x.1 + x.2 + x.3 + x.4 + x.5,
    family = f_binom, data = df_glm,
    offset = offs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst
  )
})
fits_glm <- lapply(fam_nms, function(fam_nm) {
  get(paste0("fit_", fam_nm, "_glm"))
})

## projpred ---------------------------------------------------------------

# For the binomial family with > 1 trials, we currently expect the warning
# "Using formula(x) is deprecated when x is a character vector of length > 1"
# (see GitHub issue #136), so temporarily wrap the following call in SW():
SW(refmods_glm <- lapply(fits_glm, get_refmodel))
vss_glm <- lapply(refmods_glm, varsel,
                  nclusters = nclusters_tst,
                  nclusters_pred = nclusters_pred_tst,
                  nterms_max = nterms_glm, verbose = FALSE)

# GLMMs -------------------------------------------------------------------

## Data -------------------------------------------------------------------

ngr <- 8L
xgr <- gl(n = ngr, k = floor(n_obs / ngr), length = n_obs,
          labels = paste0("gr", seq_len(ngr)))
bgr_icpts <- rnorm(ngr, sd = 0.8)
bgr_x.1 <- rnorm(ngr, sd = 0.8)
eta_glmm <- icpt +
  x_pop %*% b_pop +
  bgr_icpts[xgr] +
  bgr_x.1[xgr] * x_pop[, 1] +
  offs
df_glmm <- data.frame(y_gauss = rnorm(n_obs, f_gauss$linkinv(eta_glmm), disp),
                      y_binom = rbinom(n_obs, w_obs, f_binom$linkinv(eta_glmm)),
                      x = x_pop, x.gr = xgr, w_obs_col = w_obs, offs_col = offs)
ys_glmm <- lapply(fam_nms, function(fam_nm) {
  df_glmm[[paste0("y_", fam_nm)]]
})

# Add the number of multilevel terms to the number of population-level terms to
# obtain the total number of terms:
nterms_glmm <- nterms_glm + length(c("(1 | x.gr)", "(x.1 | x.gr)"))

## Fit --------------------------------------------------------------------

# Notes:
#   * Argument `weights` is not needed when using the cbind() syntax (for the
#     binomial family with > 1 trials).
SW({
  fit_gauss_glmm <- rstanarm::stan_glmer(
    y_gauss ~ x.1 + x.2 + x.3 + x.4 + x.5 + (x.1 | x.gr),
    family = f_gauss, data = df_glmm,
    weights = w_obs, offset = offs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_binom_glmm <- rstanarm::stan_glmer(
    cbind(y_binom, w_obs_col - y_binom) ~
      x.1 + x.2 + x.3 + x.4 + x.5 + (x.1 | x.gr),
    family = f_binom, data = df_glmm,
    offset = offs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst
  )
})
fits_glmm <- lapply(fam_nms, function(fam_nm) {
  get(paste0("fit_", fam_nm, "_glmm"))
})

## projpred ---------------------------------------------------------------

# For the binomial family with > 1 trials, we currently expect the warning
# "Using formula(x) is deprecated when x is a character vector of length > 1"
# (see GitHub issue #136), so temporarily wrap the following call in SW():
SW(refmods_glmm <- lapply(fits_glmm, get_refmodel))
vss_glmm <- lapply(refmods_glmm, varsel,
                   nclusters = nclusters_tst,
                   nclusters_pred = nclusters_pred_tst,
                   nterms_max = nterms_glmm, verbose = FALSE)

# GAMs --------------------------------------------------------------------

## Data -------------------------------------------------------------------
## Note: An alternative to mgcv::gamSim() might be the example from
## `?mgcv::concurvity` or deriving an own dataset based on the dataset for the
## GLMs above.

.Random.seed_gauss <- .Random.seed
df_gam_gauss <- mgcv::gamSim(eg = 5, n = n_obs, dist = "normal", scale = disp,
                             verbose = FALSE)
.Random.seed_bu <- .Random.seed
.Random.seed <- .Random.seed_gauss
df_gam_binom <- mgcv::gamSim(eg = 5, n = n_obs, dist = "normal", scale = 0,
                             verbose = FALSE)
.Random.seed <- .Random.seed_bu
rm(.Random.seed_gauss)
rm(.Random.seed_bu)
stopifnot(identical(
  df_gam_gauss[, setdiff(names(df_gam_gauss), "y")],
  df_gam_binom[, setdiff(names(df_gam_binom), "y")]
))
### Somehow mgcv::gamSim() always simulates 200 observations, not `n`:
df_gam_gauss <- head(df_gam_gauss, n_obs)
df_gam_binom <- head(df_gam_binom, n_obs)
###
df_gam_binom$y <- df_gam_binom$y - 4 * as.numeric(df_gam_binom$x0)
df_gam_binom$y <- rbinom(n_obs, w_obs, f_binom$linkinv(df_gam_binom$y))
df_gam <- data.frame(y_gauss = df_gam_gauss$y,
                     y_binom = df_gam_binom$y,
                     df_gam_gauss[, setdiff(names(df_gam_gauss), "y")],
                     w_obs_col = w_obs, offs_col = offs)
names(df_gam) <- sub("^x", "x.", names(df_gam))
### For shifting the enumeration:
# names(df_gam)[grep("^x", names(df_gam))] <- paste0(
#   "x.",
#   as.numeric(sub("^x", "", grep("^x", names(df_gam), value = TRUE))) + 1
# )
###
ys_gam <- lapply(fam_nms, function(fam_nm) {
  df_gam[[paste0("y_", fam_nm)]]
})
rm(df_gam_gauss)
rm(df_gam_binom)

nterms_gam <- length("x.0") + 2L * length(c("s(x.1)", "s(x.2)", "s(x.3)"))

## Fit --------------------------------------------------------------------

# Notes:
#   * Argument `weights` is not needed when using the cbind() syntax (for the
#     binomial family with > 1 trials).
#   * Argument `offset` is not supported by rstanarm::stan_gamm4(). Instead, use
#     offset() in the formula.
SW({
  fit_gauss_gam <- rstanarm::stan_gamm4(
    y_gauss ~ x.0 + s(x.1) + s(x.2) + s(x.3) + offset(offs_col),
    random = NULL,
    family = f_gauss, data = df_gam,
    weights = w_obs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_binom_gam <- rstanarm::stan_gamm4(
    cbind(y_binom, w_obs_col - y_binom) ~
      x.0 + s(x.1) + s(x.2) + s(x.3) + offset(offs_col),
    random = NULL,
    family = f_binom, data = df_gam,
    chains = chains_tst, seed = seed_tst, iter = iter_tst
  )
})
fits_gam <- lapply(fam_nms, function(fam_nm) {
  get(paste0("fit_", fam_nm, "_gam"))
})

## projpred ---------------------------------------------------------------

# For the binomial family with > 1 trials, we currently expect the warning
# "Using formula(x) is deprecated when x is a character vector of length > 1"
# (see GitHub issue #136), so temporarily wrap the following call in SW():
SW(refmods_gam <- lapply(fits_gam, get_refmodel))
vss_gam <- lapply(refmods_gam, varsel,
                  nclusters = nclusters_tst,
                  nclusters_pred = nclusters_pred_tst,
                  nterms_max = nterms_gam, verbose = FALSE)
