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
    pref_lat <- .get_refdist(refmod_crr,
                             ndraws = args_prj_i$ndraws,
                             nclusters = args_prj_i$nclusters)
    set.seed(args_prj_i_trad$seed)
    pref_trad <- .get_refdist(refmod_crr_trad,
                              ndraws = args_prj_i_trad$ndraws,
                              nclusters = args_prj_i_trad$nclusters)

    expect_identical(
      pref_lat[setdiff(names(pref_lat), c("mu", "var", "dis"))],
      pref_trad[setdiff(names(pref_trad), c("mu", "var", "dis"))],
      info = tstsetup
    )
    mu_lat_Orig <- refmod_crr$family$latent_ilink(t(refmod_crr$mu))
    mu_lat_Orig_cl <- sapply(
      seq_len(max(pref_lat$cl, na.rm = TRUE)),
      function(cl_idx) {
        # We don't use `eps` here, so there are minor differences compared to
        # .get_pclust():
        colMeans(mu_lat_Orig[which(pref_lat$cl == cl_idx), , drop = FALSE])
      }
    )
    expect_equal(mu_lat_Orig_cl, pref_trad$mu, tolerance = 1e-10,
                 info = tstsetup)
  }
  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})

# Teardown ----------------------------------------------------------------

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
