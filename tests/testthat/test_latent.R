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

# Teardown ----------------------------------------------------------------

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
