# Switches of the tests ------------------------------------------------------

# setup.R sources this file first. These switches may be set to `FALSE` to
# save time (e.g., when debugging interactively). setup.R turns off some of
# them if a package or more than one core is not available:
# Run more tests, at the downside of increased runtime?:
run_more <- identical(Sys.getenv("RUN_MORE"), "true")
# Run project()?:
run_prj <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run varsel()?:
run_vs <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run cv_varsel()?:
run_cvvs <- run_vs
# Run cv_varsel() with `validate_search = TRUE` always (`TRUE`) or just for L1
# search (`FALSE`)?:
run_valsearch_always <- FALSE
# Run cv_varsel() with `validate_search = TRUE` also for the latent and the
# augmented-data projection (and the corresponding traditional projection
# settings which are used for comparison)? (Only relevant if
# `run_valsearch_always = FALSE`.):
run_valsearch_aug_lat <- FALSE
# Run the `cvfits` test for all possible test setups (`TRUE`) or just for the
# first one among the GLMMs (`FALSE`; note that if there is no GLMM available in
# that test, the first test setup among those for K-fold CV is used)?:
run_cvfits_all <- run_more
# Run tests for `brmsfit`s?:
run_brms <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run snapshot tests?:
# Notes:
#   * For general information about snapshot tests, see, e.g.,
#   `?expect_snapshot` and `vignette("snapshotting", package = "testthat")`.
#   * The snapshot tests are at least OS-dependent (perhaps even
#   machine-dependent), so they only make sense locally. Therefore, we don't run
#   the snapshot tests on CRAN or continuous integration (CI) systems. The
#   detection of a CI system by the help of environment variable `CI` needs
#   special care, see <https://github.com/r-lib/testthat/issues/825> and the
#   source code of `testthat:::on_ci()`.
#   * The last of the following conditions avoids that the snapshot tests are
#   run by a local `R CMD check` (at least in RStudio). The reason for avoiding
#   this is that in `R CMD check`, the previous snapshots are not available (at
#   least as long as they are listed in the `.Rbuildignore` file), so they would
#   be re-created, which would throw a lot of test warnings (which could obscure
#   potentially important warnings).
run_snaps <- Sys.getenv("RUN_SNAPS")
if (identical(run_snaps, "")) {
  run_snaps <- identical(Sys.getenv("NOT_CRAN"), "true") &&
    !identical(toupper(Sys.getenv("CI")), "TRUE") &&
    identical(Sys.getenv("_R_CHECK_FORCE_SUGGESTS_"), "")
} else {
  run_snaps <- as.logical(toupper(run_snaps))
  stopifnot(isTRUE(run_snaps) || isFALSE(run_snaps))
}
# Run the tests that need parallel workers (see test_parallel.R)?:
# Notes:
#   * We don't run these tests on CRAN or continuous integration (CI) systems
#   because parallelization might require special care there.
#   * The parallelization is not tested on Windows.
run_prll <- identical(Sys.getenv("NOT_CRAN"), "true") &&
  !identical(toupper(Sys.getenv("CI")), "TRUE") &&
  !identical(.Platform$OS.type, "windows")
# Run all test scripts (following the setup script) in a completely random RNG
# state? (The tests should still pass then, because in all situations where RNG
# is used, a specific seed is supposed to be set.):
run_randRNG <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run tests for additive models (GAMs and GAMMs)?:
run_additive <- TRUE

# Use a factor or an integer response for ordinal and categorical families?:
use_fac <- TRUE

# Use polym() instead of poly()?:
use_polym <- FALSE
