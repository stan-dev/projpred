# Tests of the negative binomial (NB) family helpers (R/negbinomial.R). This
# file builds its own inputs and does not need the reference models of
# setup.R.

# The sum over y of p(y) (lgamma(y + phi) - lgamma(phi) - lgamma(y + 1)) and
# its first two derivatives in `phi`, for y ~ NB(mu, r), by the long sum:
nb_geval_exact <- function(mu, r, phi) {
  k <- qnbinom(1e-17, size = r, mu = mu, lower.tail = FALSE)
  y <- 0:k
  p <- dnbinom(y, size = r, mu = mu)
  c(sum(p * (lgamma(y + phi) - lgamma(phi) - lgamma(y + 1))),
    sum(p * (digamma(y + phi) - digamma(phi))),
    sum(p * (trigamma(y + phi) - trigamma(phi))))
}

test_that("nb_geval() agrees with the long sum", {
  cases <- expand.grid(mu = c(0.01, 3, 100, 1e4), r = c(0.3, 1, 10, 1000))
  for (i in seq_len(nrow(cases))) {
    mu <- cases$mu[i]
    r <- cases$r[i]
    gs <- nb_gsum(mu, r, 1)
    for (phi in r * c(0.3, 1, 3)) {
      ex <- nb_geval_exact(mu, r, phi)
      ap <- nb_geval(gs, phi)
      info <- paste("mu:", mu, "r:", r, "phi:", phi)
      expect_lt(abs(ap[1] - ex[1]), 1e-4, label = info)
      expect_lt(abs(ap[2] - ex[2]), 1e-6, label = info)
      expect_lt(abs(ap[3] - ex[3]), 1e-8, label = info)
    }
  }
})

test_that("nb_gsum() with the quadrature agrees with the long sum", {
  cases <- expand.grid(mu = c(0.01, 3, 60, 100, 1e3, 1e4),
                       r = c(0.1, 0.3, 1, 10, 1000, 1e5, 1e8, Inf))
  for (i in seq_len(nrow(cases))) {
    mu <- cases$mu[i]
    r <- cases$r[i]
    if (!is.finite(r) && mu > 1e3) {
      next
    }
    gs_q <- nb_gsum(mu, r, 1, method = "quadrature")
    gs_e <- nb_gsum(mu, r, 1, method = "exact")
    phis <- if (is.finite(r)) r * c(0.3, 1, 3) else c(1, 100)
    for (phi in phis) {
      ex <- if (is.finite(r)) nb_geval_exact(mu, r, phi) else nb_geval(gs_e, phi)
      ap <- nb_geval(gs_q, phi)
      info <- paste("mu:", mu, "r:", r, "phi:", phi)
      expect_lt(abs(ap[1] - ex[1]), 1e-4, label = info)
      expect_lt(abs(ap[2] - ex[2]), 1e-6, label = info)
      expect_lt(abs(ap[3] - ex[3]), 1e-8, label = info)
    }
    # The total mass is 1:
    expect_equal(sum(gs_q$p0) + sum(gs_q$bins[, "m"]), 1, tolerance = 1e-9,
                 label = paste("mass, mu:", mu, "r:", r))
  }
  # A mixture of draws with weights:
  mu_draws <- matrix(c(1, 300, 2500, 40, 120, 9), nrow = 2, byrow = TRUE)
  shp <- rep(c(0.7, 3, 25), each = 2)
  gs_q <- nb_gsum(mu_draws, shp, c(0.4, 1.2), v = c(0.5, 0.3, 0.2),
                  method = "quadrature")
  gs_e <- nb_gsum(mu_draws, shp, c(0.4, 1.2), v = c(0.5, 0.3, 0.2),
                  method = "exact")
  for (phi in c(0.5, 4, 40)) {
    expect_lt(max(abs(nb_geval(gs_q, phi) - nb_geval(gs_e, phi)) /
                    c(1e-4, 1e-6, 1e-8)), 1)
  }
  # The option selects the method (the default is the quadrature):
  expect_identical(nb_gsum(mu_draws, shp, 1),
                   nb_gsum(mu_draws, shp, 1, method = "quadrature"))
  opts_old <- options(projpred.nb_gsum = "exact")
  on.exit(options(opts_old), add = TRUE)
  expect_identical(nb_gsum(mu_draws, shp, 1),
                   nb_gsum(mu_draws, shp, 1, method = "exact"))
  expect_false(identical(nb_gsum(mu_draws, shp, 1),
                         nb_gsum(mu_draws, shp, 1, method = "quadrature")))
  options(projpred.nb_gsum = "other")
  expect_error(nb_gsum(1, 1, 1), "must be")
})

test_that("nb_gsum() stops for invalid means and shapes", {
  for (bad in list(list(NA, 1), list(Inf, 1), list(-1, 1), list(1, 0),
                   list(1, NA), list(NaN, 2))) {
    expect_error(nb_gsum(bad[[1]], bad[[2]], 1), "finite reference means",
                 info = paste(unlist(bad), collapse = " "))
  }
})

test_that("nb_gsum() handles the Poisson limit and weights", {
  # Shape Inf (Poisson), against the long Poisson sum:
  gs <- nb_gsum(c(2, 150), Inf, 1)
  for (phi in c(1, 50)) {
    ex <- sum(vapply(c(2, 150), function(mu) {
      y <- 0:qpois(1e-17, mu, lower.tail = FALSE)
      sum(dpois(y, mu) * (digamma(y + phi) - digamma(phi)))
    }, numeric(1)))
    expect_equal(nb_geval(gs, phi)[2], ex, tolerance = 1e-8)
  }
  # Weights: an observation with weight 0 does not enter, and the result is
  # linear in the weights:
  mu <- c(1, 20, 300)
  r <- c(2, 0.5, 5)
  w <- c(0.5, 0, 2)
  gs_w <- nb_gsum(mu, r, w)
  gs_sub <- nb_gsum(mu[-2], r[-2], w[-2])
  expect_identical(gs_w, gs_sub)
  sep <- 0.5 * nb_geval(nb_gsum(mu[1], r[1], 1), 3) +
    2 * nb_geval(nb_gsum(mu[3], r[3], 1), 3)
  expect_equal(nb_geval(gs_w, 3), sep, tolerance = 1e-6)
})

test_that("nb_gsum() gives the mixture of the NB distributions of draws", {
  mu_draws <- matrix(c(1, 3, 10, 2, 50, 400), nrow = 2, byrow = TRUE)
  shape_draws <- c(0.5, 2, 4)
  v <- c(0.2, 0.3, 0.5)
  w <- c(1.5, 0.7)
  gs <- nb_gsum(mu_draws, rep(shape_draws, each = 2), w, v = v)
  for (phi in c(0.3, 2, 20)) {
    ex <- Reduce(`+`, lapply(1:2, function(i) {
      w[i] * Reduce(`+`, lapply(1:3, function(s) {
        v[s] * nb_geval_exact(mu_draws[i, s], shape_draws[s], phi)
      }))
    }))
    ap <- nb_geval(gs, phi)
    info <- paste("phi:", phi)
    expect_lt(abs(ap[1] - ex[1]), 1e-4, label = info)
    expect_lt(abs(ap[2] - ex[2]), 1e-6, label = info)
    expect_lt(abs(ap[3] - ex[3]), 1e-8, label = info)
  }
})

test_that("nb_phi_step() maximizes the expected log density in the shape", {
  set.seed(5)
  n <- 40
  mu_r <- exp(rnorm(n, 1.5, 1))
  w <- runif(n, 0.5, 1.5)
  # The submodel mean equals the reference mean (identity case):
  gs <- nb_gsum(mu_r, 2.5, w)
  expect_equal(nb_phi_step(mu_r, mu_r, w, gs, phi0 = 0.7), 2.5,
               tolerance = 1e-8)
  # Another submodel mean, against optimize() on the long sum:
  mu <- mu_r * exp(rnorm(n, 0, 0.3))
  shape_r <- runif(n, 1, 4)
  gs <- nb_gsum(mu_r, shape_r, w)
  obj_exact <- function(phi) {
    sum(w * vapply(seq_len(n), function(i) {
      y <- 0:qnbinom(1e-17, size = shape_r[i], mu = mu_r[i],
                     lower.tail = FALSE)
      sum(dnbinom(y, size = shape_r[i], mu = mu_r[i]) *
            dnbinom(y, size = phi, mu = mu[i], log = TRUE))
    }, numeric(1)))
  }
  opt <- optimize(function(u) obj_exact(exp(u)), c(-5, 5), maximum = TRUE,
                  tol = 1e-10)
  expect_equal(nb_phi_step(mu, mu_r, w, gs), exp(opt$maximum),
               tolerance = 1e-6)
  # Without overdispersion in the reference, the shape becomes large (the
  # objective is flat there, so the search can stop before the bound):
  gs_pois <- nb_gsum(mu_r, Inf, w)
  phi_pois <- nb_phi_step(mu_r, mu_r, w, gs_pois)
  expect_gt(phi_pois, 1e5)
  expect_lte(phi_pois, nb_phi_max())
})

test_that("extend_family() maps the NB family names to `negbinomial`", {
  fam <- extend_family(MASS::negative.binomial(theta = 2))
  expect_identical(fam$family, "negbinomial")
  expect_identical(fam$link, "log")
  expect_true(has_dispersion(fam))
  fam_sqrt <- extend_family(MASS::negative.binomial(theta = 2, link = "sqrt"))
  expect_equal(fam_sqrt$linkinv(3), 9)
  fam_rstan <- structure(list(family = "neg_binomial_2", link = "log"),
                         class = "family")
  expect_identical(extend_family(fam_rstan)$family, "negbinomial")
  expect_error(
    extend_family(structure(list(family = "negbinomial", link = "softplus"),
                            class = "family")),
    "supports only the links"
  )
})

test_that("the NB family functions agree with dnbinom()", {
  fam <- extend_family(MASS::negative.binomial(theta = 2))
  mu <- matrix(c(0.5, 3, 40, 1, 2, 30), ncol = 2)
  dis <- c(1.5, 8)
  y <- c(0, 4, 35)
  ll <- fam$ll_fun(mu, dis, y, weights = c(1, 2, 1))
  expect_equal(ll[, 2], c(1, 2, 1) * dnbinom(y, size = 8, mu = mu[, 2],
                                             log = TRUE))
  # The reduced deviance differs from -2 * the log density by a term that does
  # not depend on `mu`:
  d1 <- fam$deviance(mu[, 1, drop = FALSE], y, dis = 1.5)
  d2 <- fam$deviance(mu[, 2, drop = FALSE], y, dis = 1.5)
  expect_equal(as.vector(d2 - d1),
               -2 * (dnbinom(y, size = 1.5, mu = mu[, 2], log = TRUE) -
                       dnbinom(y, size = 1.5, mu = mu[, 1], log = TRUE)))
  set.seed(3)
  ppd <- fam$ppd(mu[, 1], 1.5)
  expect_true(all(ppd == round(ppd)))
  # ce() and dis_fun() for one reference draw with the submodel equal to the
  # reference:
  pref <- list(mu = mu[, 1, drop = FALSE], dis = 1.5)
  # (With 3 observations, the error of the bins in nb_gsum() moves the shape
  # by about 1e-8.)
  phi <- fam$dis_fun(pref, list(mu = mu[, 1, drop = FALSE]), wobs = 1)
  expect_equal(phi, 1.5, tolerance = 1e-6)
  ce_ref <- fam$ce(pref, list(weights = rep(1, 3)),
                   list(mu = mu[, 1, drop = FALSE], dis = phi))
  ce_other <- fam$ce(pref, list(weights = rep(1, 3)),
                     list(mu = mu[, 1, drop = FALSE], dis = 3))
  expect_lt(ce_ref, ce_other)
})

test_that("get_p_clust() keeps the draws of the NB clusters", {
  fam <- extend_family(MASS::negative.binomial(theta = 2))
  set.seed(11)
  n <- 15
  S <- 12
  eta <- matrix(rnorm(n * S, 1, 0.5), n, S)
  mu <- exp(eta)
  dis <- runif(S, 1, 5)
  wdraws <- runif(S)
  cl <- rep(1:3, each = 4)
  p_ref <- get_p_clust(family = fam, eta = eta, mu = mu, mu_offs = mu,
                       dis = dis, wdraws = wdraws, cl = cl)
  wobs <- runif(n)
  wobs[3] <- 0
  for (j in 1:3) {
    ind <- which(cl == j)
    expect_identical(p_ref$nb_mix$ind[[j]], ind)
    gs_mix <- nb_gsum(mu[, ind], rep(dis[ind], each = n), wobs,
                      v = wdraws[ind] / sum(wdraws[ind]))
    expect_equal(nb_pref_gsum(p_ref, j, wobs), gs_mix)
  }
  # The cache returns the same result and is renewed for other weights:
  gs_2 <- nb_pref_gsum(p_ref, 2, wobs)
  expect_identical(gs_2, p_ref$nb_cache[["2"]]$gs)
  gs_other <- nb_pref_gsum(p_ref, 2, rep(1, n))
  expect_identical(p_ref$nb_cache[["2"]]$wobs, rep(1, n))
  expect_false(isTRUE(all.equal(gs_other, gs_2)))
  # One draw per cluster: no `nb_mix`; the shape is the one of the draw:
  p_ref1 <- get_p_clust(family = fam, eta = eta, mu = mu, mu_offs = mu,
                        dis = dis, cl = seq_len(S))
  expect_null(p_ref1$nb_mix)
  expect_equal(nb_pref_gsum(p_ref1, 4, wobs),
               nb_gsum(p_ref1$mu[, 4], dis[4], wobs), tolerance = 1e-9)
  # Other families get no NB elements:
  p_ref_pois <- get_p_clust(family = extend_family(poisson()), eta = eta,
                            mu = mu, mu_offs = mu, dis = dis, cl = cl)
  expect_false(any(c("nb_mix", "nb_cache") %in% names(p_ref_pois)))
})

test_that("glm_ridge() with a fixed NB shape agrees with glm()", {
  set.seed(7)
  n <- 60
  x <- matrix(rnorm(n * 2), n, 2)
  mu_r <- exp(0.5 + drop(x %*% c(0.4, -0.3)) + rnorm(n, 0, 0.2))
  w <- runif(n, 0.5, 1.5)
  for (link in c("log", "sqrt")) {
    fam <- nb_fixed_family(2.5, link = link)
    fit_r <- glm_ridge(x, mu_r, family = fam, weights = w, thresh = 1e-12)
    fit_g <- glm(mu_r ~ x, family = fam, weights = w,
                 control = glm.control(epsilon = 1e-14, maxit = 100))
    expect_equal(c(fit_r$beta0, fit_r$beta), unname(coef(fit_g)),
                 tolerance = 1e-6, info = link)
  }
  # The aic of the fixed family is finite for a non-integer response:
  expect_true(is.finite(fam$aic(mu_r, 1, mu_r, w, 0)))
})

test_that("nb_gsum_c() checks its arguments", {
  # The quadrature needs y0 >= 100; the exact method does not:
  expect_error(nb_gsum(5, 2, 1, y0 = 50L, method = "quadrature"),
               "invalid `y0`")
  expect_no_error(nb_gsum(5, 2, 1, y0 = 50L, method = "exact"))
  m1 <- matrix(1, 1, 1)
  expect_error(nb_gsum_c(matrix(1, 2, 1), m1, c(1, 1), 1, 100L, 0.25, 1e-17,
                         1L),
               "dimensions")
  expect_error(nb_gsum_c(m1, m1, 1, 1, 100L, 0.25, 1e-17, 2L),
               "unknown method")
  expect_error(nb_gsum_c(m1, m1, 1, 1, 0L, 0.25, 1e-17, 0L), "invalid `y0`")
  expect_error(nb_gsum_c(m1, m1, 1, 1, 100L, 0, 1e-17, 1L), "invalid `y0`")
  expect_error(nb_gsum_c(-m1, m1, 1, 1, 100L, 0.25, 1e-17, 1L),
               "invalid mean or shape")
})

test_that("nb_gsum() gives accurate probabilities for large shapes", {
  # The log pmf from sums that do not cancel. (dnbinom() has a relative error
  # of about 2e-9 at r = 1e8 and mu = 3.)
  y <- 0:20
  for (r in c(99, 1e5, 1e8)) {
    lp <- vapply(y, function(k) {
      sum(log1p((seq_len(k) - 1) / r)) - lgamma(k + 1) + k * log(3) -
        (r + k) * log1p(3 / r)
    }, numeric(1))
    for (method in c("exact", "quadrature")) {
      gs <- nb_gsum(3, r, 1, method = method)
      expect_lt(max(abs(log(gs$p0[y + 1]) - lp)), 1e-11,
                label = paste("r:", r, method))
    }
  }
})
