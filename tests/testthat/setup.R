#_______________________________________________________________________________
# Setup for the unit tests
#_______________________________________________________________________________

# General setup -----------------------------------------------------------

seed_tst <- 1235
set.seed(seed_tst)
source(testthat::test_path("helpers", "SW.R"))
mod_nms <- setNames(nm = c("glm", "glmm")) # , "gam", "gamm"
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

# Data --------------------------------------------------------------------

n_obs <- 45L

## Nonpooled ("fixed") effects --------------------------------------------

nterms_cont <- 3L
x_cont <- matrix(rnorm(n_obs * nterms_cont), n_obs, nterms_cont)
b_cont <- runif(nterms_cont, min = -0.5, max = 0.5)

nlvl_fix <- c(3L, 2L)
nlvl_fix <- setNames(nlvl_fix, seq_along(nlvl_fix))
if (length(nlvl_fix) <= 1) {
  names(nlvl_fix) <- paste0("z.", names(nlvl_fix))
}
nterms_cate <- length(nlvl_fix)
x_cate_list <- lapply(nlvl_fix, function(nlvl_fix_i) {
  x_cate <- gl(n = nlvl_fix_i, k = floor(n_obs / nlvl_fix_i), length = n_obs,
               labels = paste0("lvl", seq_len(nlvl_fix_i)))
  b_cate <- runif(nlvl_fix_i, min = -0.5, max = 0.5)
  ### Using a model.matrix() approach:
  # x_cate_mat <- model.matrix(~ 0 + x_cate)
  # eta_cate <- x_cate_mat %*% b_cate
  ###
  ### Using an indexing approach:
  eta_cate <- b_cate[x_cate]
  ###
  return(nlist(x_cate, eta_cate, b_cate))
})

icpt <- -0.42
offs <- rnorm(n_obs)
eta_glm <- icpt +
  x_cont %*% b_cont +
  do.call("+", lapply(x_cate_list, "[[", "eta_cate")) +
  offs

nterms_glm <- nterms_cont + nterms_cate

## Partially pooled ("random") effects ------------------------------------

nlvl_ran <- c(8L)
nlvl_ran <- setNames(nlvl_ran, seq_along(nlvl_ran))
if (length(nlvl_ran) <= 1) {
  names(nlvl_ran) <- paste0("z.", names(nlvl_ran))
}
# Multiply by 2 because of random intercepts and random slopes for
# `x_cont[, 1]`:
nterms_z <- length(nlvl_ran) * 2L
z_list <- lapply(nlvl_ran, function(nlvl_ran_i) {
  z <- gl(n = nlvl_ran_i, k = floor(n_obs / nlvl_ran_i), length = n_obs,
          labels = paste0("lvl", seq_len(nlvl_ran_i)))
  r_icpts <- rnorm(nlvl_ran_i, sd = 0.8)
  r_xco1 <- rnorm(nlvl_ran_i, sd = 0.8)
  eta_z <- r_icpts[z] + r_xco1[z] * x_cont[, 1]
  return(nlist(z, eta_z, r_icpts, r_xco1))
})
eta_glmm <- eta_glm +
  do.call("+", lapply(z_list, "[[", "eta_z"))

## Responses and dataset --------------------------------------------------

f_gauss <- gaussian()
f_binom <- binomial()
disp <- runif(1L, 1, 2)
wobs_tst <- sample(1:4, n_obs, replace = TRUE)
data_tst <- data.frame(
  y_gauss_glm = rnorm(n_obs, f_gauss$linkinv(eta_glm), disp),
  y_binom_glm = rbinom(n_obs, wobs_tst, f_binom$linkinv(eta_glm)),
  y_gauss_glmm = rnorm(n_obs, f_gauss$linkinv(eta_glmm), disp),
  y_binom_glmm = rbinom(n_obs, wobs_tst, f_binom$linkinv(eta_glmm)),
  xco = x_cont, xca = lapply(x_cate_list, "[[", "x_cate"),
  z = lapply(z_list, "[[", "z"),
  wobs_col = wobs_tst, offs_col = offs,
  check.names = FALSE
)
ys <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    data_tst[[paste("y", fam_nm, mod_nm, sep = "_")]]
  })
})

# Add the number of multilevel terms to the number of population-level terms to
# obtain the total number of terms:
nterms_glmm <- nterms_glm + nterms_z

# GLMs --------------------------------------------------------------------

## Fit --------------------------------------------------------------------

# Notes:
#   * Argument `weights` is not needed when using the cbind() syntax (for the
#     binomial family with > 1 trials).
SW({
  fit_gauss_glm <- rstanarm::stan_glm(
    y_gauss_glm ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2,
    family = f_gauss, data = data_tst,
    weights = wobs_tst, offset = offs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_binom_glm <- rstanarm::stan_glm(
    cbind(y_binom_glm, wobs_col - y_binom_glm) ~
      xco.1 + xco.2 + xco.3 + xca.1 + xca.2,
    family = f_binom, data = data_tst,
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

## Fit --------------------------------------------------------------------

SW({
  fit_gauss_glmm <- rstanarm::stan_glmer(
    y_gauss_glmm ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 + (xco.1 | z.1),
    family = f_gauss, data = data_tst,
    weights = wobs_tst, offset = offs,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_binom_glmm <- rstanarm::stan_glmer(
    cbind(y_binom_glmm, wobs_col - y_binom_glmm) ~
      xco.1 + xco.2 + xco.3 + xca.1 + xca.2 + (xco.1 | z.1),
    family = f_binom, data = data_tst,
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

### TODO:
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
                     y_binom = rbinom(n_obs, wobs_tst, f_binom$linkinv(eta_gam)),
                     df_gam[, setdiff(names(df_gam), "y")],
                     wobs_col = wobs_tst, offs_col = offs)
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
solterms_tst_gam <- c("s(x.1)")
###

## Fit --------------------------------------------------------------------

# Notes:
#   * Argument `offset` is not supported by rstanarm::stan_gamm4(). Instead, use
#     offset() in the formula.
SW({
  fit_gauss_gam <- rstanarm::stan_gamm4(
    y_gauss ~ x.0 + s(x.1) + s(x.2) + s(x.3) + offset(offs_col),
    random = NULL,
    family = f_gauss, data = df_gam,
    weights = wobs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_binom_gam <- rstanarm::stan_gamm4(
    cbind(y_binom, wobs_col - y_binom) ~
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

###_____________________________________________________________________________
### Currently deactivated because of issue #146:
# # GAMMs -------------------------------------------------------------------
#
# ## Data -------------------------------------------------------------------
#
# # Notes:
# #   * Have a look at the source code of mgcv::gamSim() for the underlying truth.
# #   * An alternative to mgcv::gamSim() might be the example from
# #     `?mgcv::concurvity` or deriving an own dataset based on the dataset for
# #     the GLMMs above.
# ### For `eg = 6`, mgcv::gamSim() requires `n = n_obs` to be divisible by 4
# ### (which is probably a bug):
# stopifnot(identical(n_obs %% 4L, 0L))
# ###
# ### A bug in mgcv::gamSim() causes `eg = 6` to not respect `verbose = FALSE`
# ### properly, so we use `invisible(capture.output(<...>))`:
# invisible(capture.output(
#   df_gamm <- mgcv::gamSim(eg = 6, n = n_obs, dist = "normal", scale = 0,
#                           verbose = FALSE)
# ))
# ###
# eta_gamm <- df_gamm$y - 7 * as.numeric(df_gamm$fac)
# df_gamm <- data.frame(y_gauss = rnorm(n_obs, mean = eta_gamm, sd = disp),
#                       y_binom = rbinom(n_obs, wobs_tst, f_binom$linkinv(eta_gamm)),
#                       df_gamm[, setdiff(
#                         names(df_gamm),
#                         c("y", "f", grep("^x", names(df_gamm), value = TRUE),
#                           "f3")
#                       )],
#                       wobs_col = wobs_tst, offs_col = offs)
# names(df_gamm)[names(df_gamm) == "fac"] <- "x.gr"
# names(df_gamm)[grep("^f", names(df_gamm))] <- paste0(
#   "x.",
#   as.numeric(sub("^f", "", grep("^f", names(df_gamm), value = TRUE))) + 1
# )
# df_gamm_colsBegin <- c(grep("^y", names(df_gamm), value = TRUE),
#                        sort(grep("^x", names(df_gamm), value = TRUE)))
# df_gamm <- df_gamm[, c(df_gamm_colsBegin,
#                        setdiff(names(df_gamm), df_gamm_colsBegin))]
# rm(df_gamm_colsBegin)
# ### For shifting the enumeration:
# # names(df_gamm)[grep("^x", names(df_gamm))] <- paste0(
# #   "x.",
# #   as.numeric(sub("^x", "", grep("^x", names(df_gamm), value = TRUE))) + 1
# # )
# ###
# ys_gamm <- lapply(fam_nms, function(fam_nm) {
#   df_gamm[[paste0("y_", fam_nm)]]
# })
#
# nterms_gamm <- 2L * length(c("s(x.1)", "s(x.2)", "s(x.3)")) +
#   length(c("(1 | x.gr)", "(x.1 | x.gr)"))
#
# ## Fit --------------------------------------------------------------------
#
# # Notes:
# #   * In the presence of multilevel terms (argument `random`),
# #     rstanarm::stan_gamm4() seems to be unable to support an offset() in the
# #     formula. Therefore, omit the offset here.
# #   * In the presence of multilevel terms (argument `random`),
# #     rstanarm::stan_gamm4() seems to be unable to support `QR = TRUE`.
# #     Therefore, omit `QR = TRUE` here.
# #   * In the presence of multilevel terms (argument `random`),
# #     rstanarm::stan_gamm4() seems to be unable to support the cbind() syntax
# #     (for the binomial family with > 1 trials). Therefore, use argument
# #     `weights` instead, together with `y_binom` transformed to a proportion.
# df_gamm$y_binom <- df_gamm$y_binom / wobs_tst
# SW({
#   fit_gauss_gamm <- rstanarm::stan_gamm4(
#     y_gauss ~ s(x.1) + s(x.2) + s(x.3), # + offset(offs_col)
#     random = ~ (x.1 | x.gr),
#     family = f_gauss, data = df_gamm,
#     weights = wobs_tst,
#     chains = chains_tst, seed = seed_tst, iter = iter_tst # , QR = TRUE
#   )
#   fit_binom_gamm <- rstanarm::stan_gamm4(
#     y_binom ~ s(x.1) + s(x.2) + s(x.3), # + offset(offs_col)
#     random = ~ (x.1 | x.gr),
#     family = f_binom, data = df_gamm,
#     weights = wobs_tst,
#     chains = chains_tst, seed = seed_tst, iter = iter_tst
#   )
# })
# fits_gamm <- lapply(fam_nms, function(fam_nm) {
#   get(paste0("fit_", fam_nm, "_gamm"))
# })
#
# ## projpred ---------------------------------------------------------------
#
# # For the binomial family with > 1 trials, we currently expect the warning
# # "Using formula(x) is deprecated when x is a character vector of length > 1"
# # (see GitHub issue #136), so temporarily wrap the following call in SW():
# SW(refmods_gamm <- lapply(fits_gamm, get_refmodel))
# ### To avoid issue #144:
# library(lme4)
# ###
# vss_gamm <- lapply(refmods_gamm, varsel,
#                    nclusters = nclusters_tst,
#                    nclusters_pred = nclusters_pred_tst,
#                    nterms_max = nterms_gamm, verbose = FALSE)
# ### Clean up (belongs to the fix for issue #144 above):
# detach("package:lme4")
# ###
###_____________________________________________________________________________
