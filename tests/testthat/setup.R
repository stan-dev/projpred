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

# Notes:
#   * Have a look at the source code of mgcv::gamSim() for the underlying truth.
#   * An alternative to mgcv::gamSim() might be the example from
#     `?mgcv::concurvity` or deriving an own dataset based on the dataset for
#     the GLMs above.
df_gam <- mgcv::gamSim(eg = 5, n = n_obs, dist = "normal", scale = 0,
                       verbose = FALSE)
### A bug in mgcv::gamSim() causes `eg = 5` to always simulate 200 observations,
### not `n = n_obs`:
df_gam <- head(df_gam, n_obs)
###
eta_gam <- df_gam$y - 4 * as.numeric(df_gam$x0)
df_gam <- data.frame(y_gauss = rnorm(n_obs, mean = eta_gam, sd = disp),
                     y_binom = rbinom(n_obs, w_obs, f_binom$linkinv(eta_gam)),
                     df_gam[, setdiff(names(df_gam), "y")],
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

nterms_gam <- length("x.0") + 2L * length(c("s(x.1)", "s(x.2)", "s(x.3)"))

## Fit --------------------------------------------------------------------

# Notes:
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

# GAMMs -------------------------------------------------------------------

## Data -------------------------------------------------------------------

# Notes:
#   * Have a look at the source code of mgcv::gamSim() for the underlying truth.
#   * An alternative to mgcv::gamSim() might be the example from
#     `?mgcv::concurvity` or deriving an own dataset based on the dataset for
#     the GLMMs above.
### For `eg = 6`, mgcv::gamSim() requires `n = n_obs` to be divisible by 4
### (which is probably a bug):
stopifnot(identical(n_obs %% 4L, 0L))
###
### A bug in mgcv::gamSim() causes `eg = 6` to not respect `verbose = FALSE`
### properly, so we use `invisible(capture.output(<...>))`:
invisible(capture.output(
  df_gamm <- mgcv::gamSim(eg = 6, n = n_obs, dist = "normal", scale = 0,
                          verbose = FALSE)
))
###
eta_gamm <- df_gamm$y - 7 * as.numeric(df_gamm$fac)
df_gamm <- data.frame(y_gauss = rnorm(n_obs, mean = eta_gamm, sd = disp),
                      y_binom = rbinom(n_obs, w_obs, f_binom$linkinv(eta_gamm)),
                      df_gamm[, setdiff(
                        names(df_gamm),
                        c("y", "f", grep("^x", names(df_gamm), value = TRUE),
                          "f3")
                      )],
                      w_obs_col = w_obs, offs_col = offs)
names(df_gamm)[names(df_gamm) == "fac"] <- "x.gr"
names(df_gamm)[grep("^f", names(df_gamm))] <- paste0(
  "x.",
  as.numeric(sub("^f", "", grep("^f", names(df_gamm), value = TRUE))) + 1
)
df_gamm_colsBegin <- c(grep("^y", names(df_gamm), value = TRUE),
                       sort(grep("^x", names(df_gamm), value = TRUE)))
df_gamm <- df_gamm[, c(df_gamm_colsBegin,
                       setdiff(names(df_gamm), df_gamm_colsBegin))]
rm(df_gamm_colsBegin)
### For shifting the enumeration:
# names(df_gamm)[grep("^x", names(df_gamm))] <- paste0(
#   "x.",
#   as.numeric(sub("^x", "", grep("^x", names(df_gamm), value = TRUE))) + 1
# )
###
ys_gamm <- lapply(fam_nms, function(fam_nm) {
  df_gamm[[paste0("y_", fam_nm)]]
})

nterms_gamm <- 2L * length(c("s(x.1)", "s(x.2)", "s(x.3)")) +
  length(c("(1 | x.gr)", "(x.1 | x.gr)"))

## Fit --------------------------------------------------------------------

# Notes:
#   * In the presence of multilevel terms (argument `random`),
#     rstanarm::stan_gamm4() seems to be unable to support an offset() in the
#     formula. Therefore, omit the offset here.
#   * In the presence of multilevel terms (argument `random`),
#     rstanarm::stan_gamm4() seems to be unable to support `QR = TRUE`.
#     Therefore, omit `QR = TRUE` here.
#   * In the presence of multilevel terms (argument `random`),
#     rstanarm::stan_gamm4() seems to be unable to support the cbind() syntax
#     (for the binomial family with > 1 trials). Therefore, use argument
#     `weights` instead, together with `y_binom` transformed to a proportion.
df_gamm$y_binom <- df_gamm$y_binom / w_obs
SW({
  fit_gauss_gamm <- rstanarm::stan_gamm4(
    y_gauss ~ s(x.1) + s(x.2) + s(x.3), # + offset(offs_col)
    random = ~ (x.1 | x.gr),
    family = f_gauss, data = df_gamm,
    weights = w_obs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst # , QR = TRUE
  )
  fit_binom_gamm <- rstanarm::stan_gamm4(
    y_binom ~ s(x.1) + s(x.2) + s(x.3), # + offset(offs_col)
    random = ~ (x.1 | x.gr),
    family = f_binom, data = df_gamm,
    weights = w_obs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst
  )
})
fits_gamm <- lapply(fam_nms, function(fam_nm) {
  get(paste0("fit_", fam_nm, "_gamm"))
})

## projpred ---------------------------------------------------------------

# For the binomial family with > 1 trials, we currently expect the warning
# "Using formula(x) is deprecated when x is a character vector of length > 1"
# (see GitHub issue #136), so temporarily wrap the following call in SW():
SW(refmods_gamm <- lapply(fits_gamm, get_refmodel))
### To avoid the error
### `Error in ranef(fit$mer) : could not find function "ranef"`:
library(lme4)
###
vss_gamm <- lapply(refmods_gamm, varsel,
                   nclusters = nclusters_tst,
                   nclusters_pred = nclusters_pred_tst,
                   nterms_max = nterms_gamm, verbose = FALSE)
### Clean up (belongs to the code "To avoid the error [...]" above):
detach("package:lme4")
###
