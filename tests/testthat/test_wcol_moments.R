context("weighted column moments in compiled code")

test_that(paste(
  "wcol_moments_c() gives the weighted means and sums of squares of scaled",
  "columns"
), {
  set.seed(4417)
  n <- 31
  S <- 9
  mu <- matrix(rnorm(n * S, mean = 50), n, S)
  wobs <- runif(n, 0.5, 2)
  sum_w <- sum(wobs)
  ref_moments <- function(m) {
    ybar <- colSums(m * wobs) / sum_w
    list(mean = ybar, ss = colSums(sweep(m, 2, ybar)^2 * wobs))
  }
  # All columns, scale 1:
  mom <- wcol_moments_c(mu, seq_len(S) - 1L, rep(1, S), wobs, sum_w)
  expect_equal(mom, ref_moments(mu), tolerance = 1e-12)
  # Some columns in another order, scaled:
  idx <- c(4L, 2L, 7L)
  sc <- c(0.5, 1 - 1e-10, 2)
  mom_idx <- wcol_moments_c(mu, idx - 1L, sc, wobs, sum_w)
  expect_equal(mom_idx,
               ref_moments(mu[, idx] * rep(sc, each = n)), tolerance = 1e-12)
  # Non-finite values propagate to their column only:
  mu[3, 2] <- NaN
  mom_na <- wcol_moments_c(mu, seq_len(S) - 1L, rep(1, S), wobs, sum_w)
  expect_identical(which(is.na(mom_na$ss)), 2L)
  expect_identical(which(is.na(mom_na$mean)), 2L)
  # As in colSums(mu * wobs), a weight of 0 and an infinite value give NaN:
  wobs_0 <- wobs
  wobs_0[7] <- 0
  mu_inf <- mu
  mu_inf[7, 6] <- Inf
  mom_inf <- wcol_moments_c(mu_inf, seq_len(S) - 1L, rep(1, S), wobs_0,
                            sum(wobs_0))
  expect_identical(is.na(mom_inf$mean), is.na(colSums(mu_inf * wobs_0)))
  expect_error(wcol_moments_c(mu, 0:1, 1, wobs, sum_w), "dimensions")
  expect_error(wcol_moments_c(mu, c(0L, S), c(1, 1), wobs, sum_w),
               "out of range")
})

test_that(paste(
  "get_p_clust(as_moments = TRUE) describes the centers and predictive",
  "variances of get_p_clust()"
), {
  set.seed(9310)
  n <- 25
  S <- 30
  mu <- matrix(rnorm(n * S, mean = 3), n, S)
  dis <- runif(S, 0.5, 1.5)
  wdraws <- exp(rnorm(S))
  thin <- rep(NA, S)
  thin[c(2, 9, 17, 30)] <- 1:4
  for (fam in list(extend_family(gaussian()),
                   extend_family(Student_t(nu = 5)))) {
    for (cl in list(seq_len(S), sample(S), thin)) {
      args <- list(family = fam, eta = mu, mu = mu, mu_offs = mu + 1,
                   dis = dis, wdraws = wdraws, cl = cl)
      full <- do.call(get_p_clust, args)
      mom <- do.call(get_p_clust, c(args, as_moments = TRUE))
      expect_true(isTRUE(mom$moments_only))
      ws <- rep_len(mom$ws, length(mom$idx))
      expect_equal(unclass(full$mu_offs),
                   (mu + 1)[, mom$idx] * rep(ws, each = n), tolerance = 1e-15)
      expect_equal(unclass(full$var),
                   matrix(mom$draws_predvar, n, length(mom$idx), byrow = TRUE),
                   tolerance = 1e-15)
      for (el in c("dis", "wdraws_prj", "const_wdraws_prj", "nprjdraws",
                   "cl", "wdraws_orig", "clust_used")) {
        expect_identical(mom[[el]], full[[el]], info = el)
      }
    }
  }
  # Cases without the light version:
  args <- list(family = extend_family(gaussian()), eta = mu, mu = mu,
               mu_offs = mu, dis = dis, wdraws = wdraws, as_moments = TRUE)
  expect_null(do.call(get_p_clust,
                      c(args, list(cl = rep(1:5, length.out = S)))))
  w_0 <- wdraws
  w_0[3] <- 0
  expect_null(do.call(get_p_clust,
                      modifyList(args, list(cl = seq_len(S), wdraws = w_0))))
  # A non-finite `mu_offs` is checked in perf_eval_gauss() (through the sums of
  # squares), not here:
  args$mu_offs[2, 2] <- Inf
  expect_false(is.null(do.call(get_p_clust, c(args, list(cl = seq_len(S))))))
})
