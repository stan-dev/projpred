context("predictive variance of the Gaussian and Student-t families")

# Per-observation weighted variance across the draws, computed row by row (the
# reference for rowvar_wdraws()):
rowvar_loop <- function(mu, wdraws) {
  vapply(seq_len(nrow(mu)), function(i) {
    sum(wdraws * (mu[i, ] - weighted.mean(mu[i, ], wdraws))^2)
  }, numeric(1))
}

predvar_inputs <- function(n = 30, S = 7, seed = 8152) {
  set.seed(seed)
  mu <- matrix(rnorm(n * S, mean = 2), n, S)
  wdraws <- runif(S, 0.2, 1)
  list(mu = mu, dis = runif(S, 0.5, 1.5), wdraws = wdraws / sum(wdraws))
}

test_that("rowvar_wdraws() equals the row-wise weighted variance", {
  inp <- predvar_inputs()
  expect_equal(rowvar_wdraws(inp$mu, inp$wdraws),
               rowvar_loop(inp$mu, inp$wdraws), tolerance = 1e-12)
  # Equal weights:
  w_eq <- rep(1 / ncol(inp$mu), ncol(inp$mu))
  expect_equal(rowvar_wdraws(inp$mu, w_eq), rowvar_loop(inp$mu, w_eq),
               tolerance = 1e-12)
  # One draw: the variance is exactly zero.
  expect_identical(rowvar_wdraws(inp$mu[, 1, drop = FALSE], 1),
                   rep(0, nrow(inp$mu)))
  # A missing value gives a missing variance for its observation only.
  mu_na <- inp$mu
  mu_na[4, 2] <- NA
  v_na <- rowvar_wdraws(mu_na, inp$wdraws)
  expect_identical(which(is.na(v_na)), 4L)
  expect_equal(v_na[-4], rowvar_loop(mu_na, inp$wdraws)[-4], tolerance = 1e-12)
})

test_that(paste(
  "predvar() of the Gaussian and Student-t families adds the weighted",
  "variance of mu to the mean of the squared dispersion"
), {
  inp <- predvar_inputs()
  v_mu <- rowvar_loop(inp$mu, inp$wdraws)
  fam_gauss <- extend_family(gaussian())
  expect_equal(fam_gauss$predvar(inp$mu, inp$dis, inp$wdraws),
               sum(inp$wdraws * inp$dis^2) + v_mu, tolerance = 1e-12)
  fam_t <- extend_family(Student_t(nu = 5))
  expect_equal(fam_t$predvar(inp$mu, inp$dis, inp$wdraws),
               5 / 3 * sum(inp$wdraws * inp$dis^2) + v_mu, tolerance = 1e-12)
})
