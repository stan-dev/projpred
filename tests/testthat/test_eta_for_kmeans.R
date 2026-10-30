context("eta_for_kmeans()")

test_that(paste(
  "eta_for_kmeans() reduces low-rank draws and keeps their distances and the",
  "k-means clustering"
), {
  set.seed(8571)
  x <- cbind(1, matrix(rnorm(500 * 5), 500, 5))
  beta <- matrix(rnorm(6 * 300), 6, 300)
  eta <- x %*% beta
  coords <- eta_for_kmeans(eta)
  expect_identical(nrow(coords), 300L)
  expect_lte(ncol(coords), 6L)
  expect_equal(as.matrix(dist(coords)), as.matrix(dist(t(eta))),
               tolerance = 1e-10)
  set.seed(2318)
  cl_full <- kmeans(t(eta), 5, iter.max = 50)$cluster
  set.seed(2318)
  cl_coords <- kmeans(coords, 5, iter.max = 50)$cluster
  expect_identical(cl_coords, cl_full)
})

test_that("eta_for_kmeans() returns `t(eta)` for full-rank draws", {
  set.seed(8571)
  eta <- matrix(rnorm(300 * 200), 300, 200)
  expect_identical(eta_for_kmeans(eta), t(eta))
})

test_that("eta_for_kmeans() does not change the state of the RNG", {
  set.seed(8571)
  x <- cbind(1, matrix(rnorm(500 * 5), 500, 5))
  eta <- x %*% matrix(rnorm(6 * 300), 6, 300)
  set.seed(5544)
  invisible(eta_for_kmeans(eta))
  u_after <- runif(1)
  set.seed(5544)
  expect_identical(runif(1), u_after)
})
