context("latent projection")

# Setup -------------------------------------------------------------------

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

if (exists(".Random.seed", envir = .GlobalEnv)) {
  rng_old <- get(".Random.seed", envir = .GlobalEnv)
}

# latent.R ----------------------------------------------------------------

set.seed(seed2_tst)

S <- 100L
P <- 3L
draws <- matrix(rnorm(S * P), nrow = S, ncol = P)

S_cl <- 10L
cl_draws <- sample.int(S_cl, size = S, replace = TRUE)
w_draws <- rgamma(S, shape = 4)

S_th <- 50L
idxs_thin <- round(seq(1, S, length.out = S_th))
th_draws <- rep(NA, S)
th_draws[idxs_thin] <- seq_len(S_th)

## cl_agg -----------------------------------------------------------------

test_that("cl_agg() with default arguments is the identity", {
  expect_equal(draws, cl_agg(draws), tolerance = .Machine$double.eps)
})

test_that(paste(
  "cl_agg() gives the correct output structure for clustering with constant",
  "(default) weights"
), {
  draws_cl <- cl_agg(draws, cl = cl_draws)
  expect_identical(dim(draws_cl), c(S_cl, P))
})

test_that(paste(
  "cl_agg() gives the correct output structure for clustering with nonconstant",
  "weights"
), {
  draws_cl <- cl_agg(draws, cl = cl_draws, wdraws = w_draws)
  expect_identical(dim(draws_cl), c(S_cl, P))
})

test_that(paste(
  "cl_agg() gives the correct output structure for thinning (with constant",
  "(default) weights)"
), {
  draws_th <- cl_agg(draws, cl = th_draws)
  expect_identical(dim(draws_th), c(S_th, P))
})

# Linear predictors -------------------------------------------------------

test_that(paste(
  "for the latent projection, `<refmodel_object>$mu` is the same as",
  "`<refmodel_object>$eta`"
), {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.latent\\.", names(prjs), value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    refmod_crr <- refmods[[args_prj_i$tstsetup_ref]]
    expect_identical(refmod_crr$mu, refmod_crr$eta, info = tstsetup)
  }
})

test_that(paste(
  "`<refmodel_object>$eta` in case of the `brnll` family is the same",
  "no matter whether `<refmodel_object>` was set up for the latent or the",
  "traditional projection"
), {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.brnll\\..*\\.latent\\.", names(prjs), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.latent\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(prjs)) next

    args_prj_i <- args_prj[[tstsetup]]
    args_prj_i_trad <- args_prj[[tstsetup_trad]]
    refmod_crr <- refmods[[args_prj_i$tstsetup_ref]]
    refmod_crr_trad <- refmods[[args_prj_i_trad$tstsetup_ref]]

    expect_identical(refmod_crr$eta, refmod_crr_trad$eta, info = tstsetup)
  }
})

# Comparison with traditional projection ----------------------------------

## Clustering -------------------------------------------------------------

test_that(paste(
  "clustering `<refmodel_object>$mu` in case of the `brnll` family is the same",
  "no matter whether `<refmodel_object>` was set up for the latent or the",
  "traditional projection"
), {
  skip_if_not(run_prj)
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  tstsetups <- grep("\\.brnll\\..*\\.latent\\.", names(prjs), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_trad <- sub("\\.latent\\.", ".trad_compare.", tstsetup)
    if (!tstsetup_trad %in% names(prjs)) next

    args_prj_i <- args_prj[[tstsetup]]
    args_prj_i_trad <- args_prj[[tstsetup_trad]]
    refmod_crr <- refmods[[args_prj_i$tstsetup_ref]]
    refmod_crr_trad <- refmods[[args_prj_i_trad$tstsetup_ref]]
    set.seed(args_prj_i$seed)
    pref_lat <- get_refdist(refmod_crr,
                            ndraws = args_prj_i$ndraws,
                            nclusters = args_prj_i$nclusters)
    set.seed(args_prj_i_trad$seed)
    pref_trad <- get_refdist(refmod_crr_trad,
                             ndraws = args_prj_i_trad$ndraws,
                             nclusters = args_prj_i_trad$nclusters)

    expect_identical(
      pref_lat[setdiff(names(pref_lat), c("mu", "mu_offs", "var", "dis"))],
      pref_trad[setdiff(names(pref_trad), c("mu", "mu_offs", "var", "dis"))],
      info = tstsetup
    )

    mu_lat_oscale <- refmod_crr$family$latent_ilink(t(refmod_crr$mu))
    mu_lat_oscale_cl <- sapply(
      seq_len(max(pref_lat$cl, na.rm = TRUE)),
      function(cl_idx) {
        # We don't use `eps` here, so there are minor differences compared to
        # .get_pclust():
        colMeans(mu_lat_oscale[which(pref_lat$cl == cl_idx), , drop = FALSE])
      }
    )
    expect_equal(mu_lat_oscale_cl, pref_trad$mu, tolerance = 1e-10,
                 info = tstsetup)

    mu_offs_lat_oscale <- refmod_crr$family$latent_ilink(t(refmod_crr$mu_offs))
    mu_offs_lat_oscale_cl <- sapply(
      seq_len(max(pref_lat$cl, na.rm = TRUE)),
      function(cl_idx) {
        # We don't use `eps` here, so there are minor differences compared to
        # .get_pclust():
        colMeans(
          mu_offs_lat_oscale[which(pref_lat$cl == cl_idx), , drop = FALSE]
        )
      }
    )
    expect_equal(mu_offs_lat_oscale_cl, pref_trad$mu_offs, tolerance = 1e-10,
                 info = tstsetup)
  }
  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})

## Gaussian family --------------------------------------------------------

test_that(paste(
  "for the gaussian() family, the latent projection is the same as the",
  "traditional projection (when setting `dis` appropriately)"
), {
  skip_if_not(run_prj)
  tstsetups <- grep("\\.gauss\\.", names(prjs), value = TRUE)
  if (any(grepl("\\.without_wobs\\.", tstsetups))) {
    message("The test \"for the gaussian() family, the latent projection is ",
            "the same as the traditional projection (when setting `dis` ",
            "appropriately)\" could be simplified.")
    tstsetups <- grep("\\.without_wobs\\.", tstsetups, value = TRUE,
                      invert = TRUE)
  }
  if (any(grepl("^brms\\.", tstsetups))) {
    message("The test \"for the gaussian() family, the latent projection is ",
            "the same as the traditional projection (when setting `dis` ",
            "appropriately)\" could be extended (to `brmsfit`s).")
    tstsetups <- grep("^brms\\.", tstsetups, value = TRUE, invert = TRUE)
  }
  fits_no_wobs <- NULL
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])
    if (!args_prj_i$tstsetup_fit %in% names(fits_no_wobs)) {
      args_fit_i <- args_fit[[args_prj_i$tstsetup_fit]]
      args_fit_i$weights <- NULL
      fit_fun_crr <- switch(args_fit_i$mod_nm,
                            "glm" = "stan_glm",
                            "glmm" = "stan_glmer",
                            "stan_gamm4")
      fits_no_wobs[[args_prj_i$tstsetup_fit]] <- suppressWarnings(
        do.call(get(fit_fun_crr, asNamespace(args_fit_i$pkg_nm)),
                excl_nonargs(args_fit_i))
      )
      rm(args_fit_i)
      rm(fit_fun_crr)
    }
    fit_crr <- fits_no_wobs[[args_prj_i$tstsetup_fit]]
    dis_crr <- as.matrix(fit_crr)[, "sigma"]
    refmod_crr_trad <- suppressMessages(do.call(get_refmodel, c(
      list(object = fit_crr),
      excl_nonargs(args_ref[[args_prj_i$tstsetup_ref]])
    )))
    prj_crr_trad <- do.call(project, c(
      list(object = refmod_crr_trad),
      excl_nonargs(args_prj_i)
    ))
    refmod_crr_lat <- suppressMessages(do.call(get_refmodel, c(
      list(object = fit_crr),
      excl_nonargs(args_ref[[args_prj_i$tstsetup_ref]]),
      list(latent = TRUE, dis = dis_crr)
    )))
    prj_crr_lat <- do.call(project, c(
      list(object = refmod_crr_lat),
      excl_nonargs(args_prj_i)
    ))
    tstsetup_mod <- sub("\\.with_wobs\\.", ".without_wobs.", tstsetup)
    if (ndr_ncl$clust_used) {
      warn_prjmat_expect <- "The projected draws have different .*weights"
    } else {
      warn_prjmat_expect <- NA
    }
    expect_warning(prjmat_crr_lat <- as.matrix(prj_crr_lat),
                   warn_prjmat_expect, info = tstsetup_mod)
    expect_warning(prjmat_crr_trad <- as.matrix(prj_crr_trad),
                   warn_prjmat_expect, info = tstsetup_mod)
    expect_identical(prjmat_crr_lat, prjmat_crr_trad, info = tstsetup_mod)
    prj_el_excl <- c("outdmin", "refmodel")
    expect_identical(prj_crr_lat[setdiff(names(prj_crr_lat), prj_el_excl)],
                     prj_crr_trad[setdiff(names(prj_crr_trad), prj_el_excl)],
                     info = tstsetup_mod)
  }
})

# Teardown ----------------------------------------------------------------

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
