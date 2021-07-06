#_______________________________________________________________________________
# Setup for the unit tests
#_______________________________________________________________________________

# General setup -----------------------------------------------------------

# When debugging interactively without needing the "vsel" objects, these
# switches may be set to `FALSE` to source() this script faster:
run_vs <- identical(Sys.getenv("NOT_CRAN"), "true")
run_cvvs <- run_vs

set.seed(8541351)

source(testthat::test_path("helpers", "SW.R"))
source(testthat::test_path("helpers", "unlist_cust.R"))
source(testthat::test_path("helpers", "testers.R"))

mod_nms <- setNames(nm = c("glm", "glmm", "gam", "gamm"))
### Exclude GAMs because of issue #150:
mod_nms <- setNames(nm = setdiff(mod_nms, "gam"))
###
### Exclude GAMMs because of issue #148:
mod_nms <- setNames(nm = setdiff(mod_nms, "gamm"))
###
fam_nms <- setNames(nm = c("gauss", "binom", "poiss"))

seed_tst <- 74345

# rstanarm setup ----------------------------------------------------------

chains_tst <- 2L
iter_tst <- 500L

# projpred setup ----------------------------------------------------------

## Defaults ---------------------------------------------------------------

ndraws_default <- 400L # Adopt this if the default is changed.
nresample_clusters_default <- 1000L # Adopt this if the default is changed.
projection_nms <- c(
  "dis", "kl", "weights", "solution_terms", "sub_fit", "family",
  "p_type", "intercept", "extract_model_data", "refmodel"
)
sub_fit_nms <- c("alpha", "beta", "w", "formula", "x", "y")

## Customized -------------------------------------------------------------

seed2_tst <- 866028
nclusters_tst <- 2L
nclusters_pred_tst <- 3L
ndr_ncl_pred_tst <- list(
  default = list(),
  noclust = list(ndraws = 25L),
  clust = list(nclusters = nclusters_pred_tst),
  clust_draws = list(ndraws = 3L),
  clust1 = list(nclusters = 1L)
)
### Because of issue #149:
# solterms_x <- c("xco.2", "xca.1")
solterms_x <- c("xco.2", "xco.1")
###
solterms_z <- c("(1 | z.1)", "(xco.1 | z.1)")
solterms_s <- c("s(s.1)", "s(s.2)")

# Data --------------------------------------------------------------------

n_tst <- 45L

## GLMs --------------------------------------------------------------------
## Add nonpooled ("fixed") effects to the intercept-(and-offset-)only model

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

## GLMMs ------------------------------------------------------------------
## Add partially pooled ("random") effects to the GLMs

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

## GAMs -------------------------------------------------------------------
## Add nonlinear (smoothed) effects to the GLMs

# For simplicity, always use the same nonlinear function (could be extended in
# the future):
s_mat <- apply(x_cont, 2, function(x, a = -3, b = 0.75, c = 0.5) {
  b * (x - c)^2 + a
})
s_sum <- rowSums(s_mat)
nterms_s <- ncol(s_mat)
eta_gam <- eta_glm +
  s_sum

# Multiply by 2 because of the baseline linear term as well as the standard
# deviation for the wiggliness around it:
nterms_gam <- nterms_glm + 2L * nterms_s

## GAMMs ------------------------------------------------------------------
## Add nonlinear (smoothed) effects to the GLMMs

eta_gamm <- eta_glmm +
  s_sum

nterms_gamm <- nterms_glmm + 2L * nterms_s

## Combined dataset -------------------------------------------------------

f_gauss <- gaussian()
f_binom <- binomial()
f_poiss <- poisson()
dis_tst <- runif(1L, 1, 2)
wobs_tst <- sample(1:4, n_tst, replace = TRUE)
dat <- data.frame(
  y_glm_gauss = rnorm(n_tst, f_gauss$linkinv(eta_glm), dis_tst),
  y_glm_binom = rbinom(n_tst, wobs_tst, f_binom$linkinv(eta_glm)),
  y_glm_poiss = rpois(n_tst, f_poiss$linkinv(eta_glm)),
  y_glmm_gauss = rnorm(n_tst, f_gauss$linkinv(eta_glmm), dis_tst),
  y_glmm_binom = rbinom(n_tst, wobs_tst, f_binom$linkinv(eta_glmm)),
  y_glmm_poiss = rpois(n_tst, f_poiss$linkinv(eta_glmm)),
  y_gam_gauss = rnorm(n_tst, f_gauss$linkinv(eta_gam), dis_tst),
  y_gam_binom = rbinom(n_tst, wobs_tst, f_binom$linkinv(eta_gam)),
  y_gam_poiss = rpois(n_tst, f_poiss$linkinv(eta_gam)),
  y_gamm_gauss = rnorm(n_tst, f_gauss$linkinv(eta_gamm), dis_tst),
  y_gamm_binom = rbinom(n_tst, wobs_tst, f_binom$linkinv(eta_gamm)),
  y_gamm_poiss = rpois(n_tst, f_poiss$linkinv(eta_gamm)),
  xco = x_cont, xca = lapply(x_cate_list, "[[", "x_cate"),
  z = lapply(z_list, "[[", "z"),
  s = s_mat,
  wobs_col = wobs_tst, offs_col = offs_tst,
  check.names = FALSE
)

## nterms -----------------------------------------------------------------

ntermss <- sapply(mod_nms, function(mod_nm) {
  get(paste("nterms", mod_nm, sep = "_"))
})
nterms_max_tst <- min(ntermss)

nterms_unavail <- list(
  single = nterms_max_tst + 130L,
  vec = c(nterms_max_tst + 130L, nterms_max_tst + 290L)
)
nterms_avail <- list(
  default = NULL,
  empty = 0L,
  single = nterms_max_tst %/% 2L,
  subvec = as.integer(round(seq(0, nterms_max_tst, length.out = 3))),
  full = 0:nterms_max_tst
)

# Fits --------------------------------------------------------------------

# Notes:
#   * In principle, one could also generalize the construction of the formulas.
#     However, rstanarm seems to have problems with that (see
#     <https://github.com/stan-dev/projpred/issues/65#issuecomment-765345522>).
#   * Argument `weights` is not needed when using the cbind() syntax (for the
#     binomial family with > 1 trials).

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
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_glm_poiss <- rstanarm::stan_glm(
    y_glm_poiss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2,
    family = f_poiss, data = dat,
    weights = wobs_tst, offset = offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
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
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
  fit_glmm_poiss <- rstanarm::stan_glmer(
    y_glmm_poiss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 + (xco.1 | z.1),
    family = f_poiss, data = dat,
    weights = wobs_tst, offset = offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  )
})

## GAMs -------------------------------------------------------------------

# Notes:
#   * Argument `offset` is not supported by rstanarm::stan_gamm4(). Instead, use
#     offset() in the formula.

### Exclude GAMs because of issue #150:
# SW({
#   fit_gam_gauss <- rstanarm::stan_gamm4(
#     y_gam_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 +
#       s(s.1) + s(s.2) + s(s.3) + offset(offs_col),
#     family = f_gauss, data = dat,
#     weights = wobs_tst,
#     chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
#   )
#   fit_gam_binom <- rstanarm::stan_gamm4(
#     cbind(y_gam_binom, wobs_col - y_gam_binom) ~
#       xco.1 + xco.2 + xco.3 + xca.1 + xca.2 +
#       s(s.1) + s(s.2) + s(s.3) + offset(offs_col),
#     family = f_binom, data = dat,
#     chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
#   )
#   fit_gam_poiss <- rstanarm::stan_gamm4(
#     y_gam_poiss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 +
#       s(s.1) + s(s.2) + s(s.3) + offset(offs_col),
#     family = f_poiss, data = dat,
#     weights = wobs_tst,
#     chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
#   )
# })
###

## GAMMs ------------------------------------------------------------------

# Notes:
#   * In the presence of multilevel terms (argument `random`),
#     rstanarm::stan_gamm4() seems to be unable to support an offset() in the
#     formula. Therefore, omit the offset here.

### Exclude GAMMs because of issue #148:
# SW({
#   fit_gamm_gauss <- rstanarm::stan_gamm4(
#     y_gamm_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 +
#       s(s.1) + s(s.2) + s(s.3), # + offset(offs_col)
#     random = ~ (xco.1 | z.1),
#     family = f_gauss, data = dat,
#     weights = wobs_tst,
#     chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
#   )
#   fit_gamm_binom <- rstanarm::stan_gamm4(
#     cbind(y_gamm_binom, wobs_col - y_gamm_binom) ~
#       xco.1 + xco.2 + xco.3 + xca.1 + xca.2 +
#       s(s.1) + s(s.2) + s(s.3), # + offset(offs_col)
#     random = ~ (xco.1 | z.1),
#     family = f_binom, data = dat,
#     chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
#   )
#   fit_gamm_poiss <- rstanarm::stan_gamm4(
#     y_gamm_poiss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 +
#       s(s.1) + s(s.2) + s(s.3), # + offset(offs_col)
#     random = ~ (xco.1 | z.1),
#     family = f_poiss, data = dat,
#     weights = wobs_tst,
#     chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
#   )
# })
###

## List -------------------------------------------------------------------

fits <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    get(paste("fit", mod_nm, fam_nm, sep = "_"))
  })
})
rm(list = grep("^fit_", ls(), value = TRUE))

# projpred ----------------------------------------------------------------

## Reference model --------------------------------------------------------

# For the binomial family with > 1 trials, we currently expect the warning
# "Using formula(x) is deprecated when x is a character vector of length > 1"
# (see GitHub issue #136), so temporarily wrap the following call in SW():
SW(refmods <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    get_refmodel(fits[[mod_nm]][[fam_nm]])
  })
}))

## Variable selection -----------------------------------------------------

### Exclude GAMMs because of issue #148:
# ### To avoid issue #144 (for GAMMs):
# if (run_vs || run_cvvs) {
#   library(lme4)
# }
# ###
###
if (run_vs) {
  args_vs <- lapply(mod_nms, function(mod_nm) {
    lapply(fam_nms, function(fam_nm) {
      return(nlist(
        mod_nm, fam_nm, nclusters = nclusters_tst,
        nclusters_pred = nclusters_pred_tst, nterms_max = nterms_max_tst,
        verbose = FALSE, seed = seed_tst
      ))
    })
  })
  args_vs <- unlist_cust(args_vs)

  vss <- lapply(args_vs, function(args_vs_i) {
    do.call(varsel, c(
      list(object = refmods[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]]),
      args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
    ))
  })
}
if (run_cvvs) {
  args_cvvs <- lapply(mod_nms, function(mod_nm) {
    lapply(fam_nms, function(fam_nm) {
      return(nlist(
        mod_nm, fam_nm, nclusters = nclusters_tst,
        nclusters_pred = nclusters_pred_tst, nterms_max = nterms_max_tst,
        verbose = FALSE, seed = seed_tst
      ))
    })
  })
  args_cvvs <- unlist_cust(args_cvvs)

  # Use SW() because of occasional warnings concerning Pareto k diagnostics:
  SW(cvvss <- lapply(args_cvvs, function(args_cvvs_i) {
    do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$mod_nm]][[args_cvvs_i$fam_nm]]),
      args_cvvs_i[setdiff(names(args_cvvs_i), c("mod_nm", "fam_nm"))]
    ))
  }))
}
### Exclude GAMMs because of issue #148:
# ### Clean up (belongs to the fix for issue #144 above):
# if (run_vs || run_cvvs) {
#   detach("package:lme4")
# }
# ###
###

## Projection -------------------------------------------------------------

### From "refmodel" -------------------------------------------------------

args_prj <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    solterms <- nlist(empty = character(), solterms_x)
    if (mod_nm %in% c("glmm", "gamm")) {
      solterms <- c(solterms,
                    nlist(solterms_z, solterms_xz = c(solterms_x, solterms_z)))
    }
    if (mod_nm %in% c("gam", "gamm")) {
      solterms <- c(solterms,
                    nlist(solterms_s, solterms_xs = c(solterms_x, solterms_s)))
    }
    if (mod_nm == "gamm") {
      solterms <- c(solterms,
                    nlist(solterms_sz = c(solterms_s, solterms_z),
                          solterms_xsz = c(solterms_x, solterms_s, solterms_z)))
    }
    if (fam_nm != "gauss") {
      solterms <- tail(solterms, 1)
    }
    lapply(solterms, function(solterms_i) {
      if (mod_nm == "glm" && fam_nm == "gauss") {
        ndr_ncl_pred <- ndr_ncl_pred_tst
      } else {
        ndr_ncl_pred <- ndr_ncl_pred_tst["clust"]
      }
      lapply(ndr_ncl_pred, function(ndr_ncl_pred_i) {
        return(c(
          nlist(mod_nm, fam_nm, solution_terms = solterms_i, seed = seed_tst),
          ndr_ncl_pred_i
        ))
      })
    })
  })
})
args_prj <- unlist_cust(args_prj)

prjs <- lapply(args_prj, function(args_prj_i) {
  do.call(project, c(
    list(object = refmods[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]]),
    args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
  ))
})

### From "vsel" -----------------------------------------------------------

args_prj_vs <- lapply(mod_nms["glm"], function(mod_nm) {
  lapply(fam_nms["gauss"], function(fam_nm) {
    lapply(nterms_avail, function(nterms_crr) {
      args_out <- nlist(mod_nm, fam_nm,
                        nclusters = nclusters_pred_tst, seed = seed_tst)
      if (!is.null(nterms_crr)) {
        args_out <- c(args_out, list(nterms = nterms_crr))
      }
      return(args_out)
    })
  })
})
args_prj_vs <- unlist_cust(args_prj_vs)

if (run_vs) {
  prjs_vs <- lapply(args_prj_vs, function(args_prj_vs_i) {
    tstsetup <- grep(
      paste0("^", args_prj_vs_i$mod_nm, "\\.", args_prj_vs_i$fam_nm),
      names(vss),
      value = TRUE
    )
    stopifnot(length(tstsetup) == 1)
    do.call(project, c(
      list(object = vss[[tstsetup]]),
      args_prj_vs_i[setdiff(names(args_prj_vs_i), c("mod_nm", "fam_nm"))]
    ))
  })
}

args_prj_cvvs <- args_prj_vs

if (run_cvvs) {
  # Use SW() because of occasional pwrssUpdate() warnings:
  SW(prjs_cvvs <- lapply(args_prj_cvvs, function(args_prj_cvvs_i) {
    tstsetup <- grep(
      paste0("^", args_prj_cvvs_i$mod_nm, "\\.", args_prj_cvvs_i$fam_nm),
      names(cvvss),
      value = TRUE
    )
    stopifnot(length(tstsetup) == 1)
    do.call(project, c(
      list(object = cvvss[[tstsetup]]),
      args_prj_cvvs_i[setdiff(names(args_prj_cvvs_i), c("mod_nm", "fam_nm"))]
    ))
  }))
}
