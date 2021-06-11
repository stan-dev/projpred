#_______________________________________________________________________________
# Setup for the unit tests
#_______________________________________________________________________________

# General setup -----------------------------------------------------------

seed_tst <- 1235
set.seed(seed_tst)
source(testthat::test_path("helpers", "SW.R"))
mod_nms <- setNames(nm = c("glm", "glmm", "gam")) # , "gamm"
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
solterms <- c("x.2", "x.4")

# Data --------------------------------------------------------------------

n_tst <- 45L

## Nonpooled ("fixed") effects --------------------------------------------

nterms_cont <- 3L
x_cont <- matrix(rnorm(n_tst * nterms_cont), n_tst, nterms_cont)
b_cont <- runif(nterms_cont, min = -0.5, max = 0.5)

nlvl_fix <- c(3L, 2L)
nlvl_fix <- setNames(nlvl_fix, seq_along(nlvl_fix))
if (length(nlvl_fix) <= 1) {
  names(nlvl_fix) <- paste0("z.", names(nlvl_fix))
}
nterms_cate <- length(nlvl_fix)
x_cate_list <- lapply(nlvl_fix, function(nlvl_fix_i) {
  x_cate <- gl(n = nlvl_fix_i, k = floor(n_tst / nlvl_fix_i), length = n_tst,
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
offs_tst <- rnorm(n_tst)
eta_glm <- icpt +
  x_cont %*% b_cont +
  do.call("+", lapply(x_cate_list, "[[", "eta_cate")) +
  offs_tst

nterms_glm <- nterms_cont + nterms_cate

## Partially pooled ("random") effects ------------------------------------

nlvl_ran <- c(8L)
nlvl_ran <- setNames(nlvl_ran, seq_along(nlvl_ran))
if (length(nlvl_ran) <= 1) {
  names(nlvl_ran) <- paste0("z.", names(nlvl_ran))
}
# Multiply by 2 because of the random intercept as well as the random slope (the
# latter for `x_cont[, 1]`):
nterms_z <- length(nlvl_ran) * 2L
z_list <- lapply(nlvl_ran, function(nlvl_ran_i) {
  z <- gl(n = nlvl_ran_i, k = floor(n_tst / nlvl_ran_i), length = n_tst,
          labels = paste0("lvl", seq_len(nlvl_ran_i)))
  r_icpts <- rnorm(nlvl_ran_i, sd = 0.8)
  r_xco1 <- rnorm(nlvl_ran_i, sd = 0.8)
  eta_z <- r_icpts[z] + r_xco1[z] * x_cont[, 1]
  return(nlist(z, eta_z, r_icpts, r_xco1))
})
eta_glmm <- eta_glm +
  do.call("+", lapply(z_list, "[[", "eta_z"))

nterms_glmm <- nterms_glm + nterms_z

## Nonlinear (smoothed) effects -------------------------------------------

# For simplicity, always use the same nonlinear function (could be extended in
# the future):
s_mat <- apply(x_cont, 2, function(x, a = -3, b = 0.75, c = 0.5) {
  b * (x - c)^2 + a
})
nterms_s <- ncol(s_mat)
eta_gam <- eta_glm +
  rowSums(s_mat)

# Multiply by 2 because of the baseline linear term as well as the standard
# deviation for the wiggliness around it:
nterms_gam <- nterms_glm + 2L * nterms_s

## Responses and dataset --------------------------------------------------

f_gauss <- gaussian()
f_binom <- binomial()
dis_tst <- runif(1L, 1, 2)
wobs_tst <- sample(1:4, n_tst, replace = TRUE)
dat <- data.frame(
  y_glm_gauss = rnorm(n_tst, f_gauss$linkinv(eta_glm), dis_tst),
  y_glm_binom = rbinom(n_tst, wobs_tst, f_binom$linkinv(eta_glm)),
  y_glmm_gauss = rnorm(n_tst, f_gauss$linkinv(eta_glmm), dis_tst),
  y_glmm_binom = rbinom(n_tst, wobs_tst, f_binom$linkinv(eta_glmm)),
  y_gam_gauss = rnorm(n_tst, f_gauss$linkinv(eta_gam), dis_tst),
  y_gam_binom = rbinom(n_tst, wobs_tst, f_binom$linkinv(eta_gam)),
  xco = x_cont, xca = lapply(x_cate_list, "[[", "x_cate"),
  z = lapply(z_list, "[[", "z"),
  s = s_mat,
  wobs_col = wobs_tst, offs_col = offs_tst,
  check.names = FALSE
)
ys <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    dat[[paste("y", mod_nm, fam_nm, sep = "_")]]
  })
})

# Fits --------------------------------------------------------------------

# Notes:
#   * In principle, one could also generalize the construction of the formulas.
#     However, rstanarm seems to have problems with that (see
#     <https://github.com/stan-dev/projpred/issues/65#issuecomment-765345522>).
#   * Argument `weights` is not needed when using the cbind() syntax (for the
#     binomial family with > 1 trials).
#   * Argument `offset` is not supported by rstanarm::stan_gamm4(). Instead, use
#     offset() in the formula.

## GLMs -------------------------------------------------------------------

SW({
  fit_glm_gauss <- rstanarm::stan_glm(
    y_glm_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2,
    family = f_gauss, data = dat,
    weights = wobs_tst, offset = offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_glm_binom <- rstanarm::stan_glm(
    cbind(y_glm_binom, wobs_col - y_glm_binom) ~
      xco.1 + xco.2 + xco.3 + xca.1 + xca.2,
    family = f_binom, data = dat,
    offset = offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst
  )
})

## GLMMs ------------------------------------------------------------------

SW({
  fit_glmm_gauss <- rstanarm::stan_glmer(
    y_glmm_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 + (xco.1 | z.1),
    family = f_gauss, data = dat,
    weights = wobs_tst, offset = offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_glmm_binom <- rstanarm::stan_glmer(
    cbind(y_glmm_binom, wobs_col - y_glmm_binom) ~
      xco.1 + xco.2 + xco.3 + xca.1 + xca.2 + (xco.1 | z.1),
    family = f_binom, data = dat,
    offset = offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst
  )
})

## GAMs -------------------------------------------------------------------

SW({
  fit_gam_gauss <- rstanarm::stan_gamm4(
    y_gam_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 +
      s(s.1) + s(s.2) + s(s.3) + offset(offs_col),
    family = f_gauss, data = dat,
    weights = wobs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_gam_binom <- rstanarm::stan_gamm4(
    cbind(y_gam_binom, wobs_col - y_gam_binom) ~
      xco.1 + xco.2 + xco.3 + xca.1 + xca.2 +
      s(s.1) + s(s.2) + s(s.3) + offset(offs_col),
    family = f_binom, data = dat,
    chains = chains_tst, seed = seed_tst, iter = iter_tst
  )
})

## List -------------------------------------------------------------------

fits <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    get(paste("fit", mod_nm, fam_nm, sep = "_"))
  })
})
rm(list = grep("^fit_", ls(), value = TRUE))

# projpred ----------------------------------------------------------------

# For the binomial family with > 1 trials, we currently expect the warning
# "Using formula(x) is deprecated when x is a character vector of length > 1"
# (see GitHub issue #136), so temporarily wrap the following call in SW():
SW(refmods <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    get_refmodel(fits[[mod_nm]][[fam_nm]])
  })
}))
SW(vss <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    varsel(refmods[[mod_nm]][[fam_nm]],
           nclusters = nclusters_tst,
           nclusters_pred = nclusters_pred_tst,
           nterms_max = nterms_glm, verbose = FALSE)
  })
}))

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
# ### For `eg = 6`, mgcv::gamSim() requires `n = n_tst` to be divisible by 4
# ### (which is probably a bug):
# stopifnot(identical(n_tst %% 4L, 0L))
# ###
# ### A bug in mgcv::gamSim() causes `eg = 6` to not respect `verbose = FALSE`
# ### properly, so we use `invisible(capture.output(<...>))`:
# invisible(capture.output(
#   df_gamm <- mgcv::gamSim(eg = 6, n = n_tst, dist = "normal", scale = 0,
#                           verbose = FALSE)
# ))
# ###
# eta_gamm <- df_gamm$y - 7 * as.numeric(df_gamm$fac)
# df_gamm <- data.frame(y_gauss = rnorm(n_tst, mean = eta_gamm, sd = dis_tst),
#                       y_binom = rbinom(n_tst, wobs_tst, f_binom$linkinv(eta_gamm)),
#                       df_gamm[, setdiff(
#                         names(df_gamm),
#                         c("y", "f", grep("^x", names(df_gamm), value = TRUE),
#                           "f3")
#                       )],
#                       wobs_col = wobs_tst, offs_col = offs_tst)
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
