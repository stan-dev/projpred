withr::local_seed(6130)

# The cluster loop of get_p_clust() for given cluster indices `cl` (the
# reference for the computation of all clusters at once):
p_clust_loop <- function(family, mu, mu_offs, dis, wdraws, cl) {
  nclusters <- max(cl, na.rm = TRUE)
  centers <- matrix(0, nrow = dim(mu)[1], ncol = nclusters)
  centers_offs <- matrix(0, nrow = dim(mu_offs)[1], ncol = nclusters)
  wcluster <- rep(0, nclusters)
  dis_agg <- rep(NA_real_, nclusters)
  predvar <- matrix(nrow = dim(mu)[1], ncol = nclusters)
  eps <- 1e-10
  for (j in seq_len(nclusters)) {
    ind <- which(cl == j)
    ws <- wdraws[ind] / sum(wdraws[ind]) * (1 - eps)
    centers[, j] <- mu[, ind, drop = FALSE] %*% ws
    centers_offs[, j] <- mu_offs[, ind, drop = FALSE] %*% ws
    wcluster[j] <- sum(wdraws[ind])
    dis_agg[j] <- crossprod(dis[ind], ws)
    predvar[, j] <- family$predvar(mu_offs[, ind, drop = FALSE], dis[ind], ws)
  }
  list(mu = unname(centers), mu_offs = unname(centers_offs), var = predvar,
       dis = dis_agg, wdraws_prj = wcluster / sum(wcluster))
}

p_clust_inputs <- function(family, n = 23, S = 40, seed = 6129) {
  withr::local_seed(seed)
  eta <- matrix(rnorm(n * S), n, S)
  mu <- family$linkinv(eta)
  dis <- if (family$family %in% c("gaussian", "Student_t")) {
    runif(S, 0.5, 1.5)
  } else {
    rep(NA, S)
  }
  list(eta = eta, mu = mu, mu_offs = mu + 0.25, dis = dis,
       wdraws = exp(rnorm(S)))
}

# With one draw per cluster, get_p_clust() gives identical() results; with
# several draws in a cluster, the results differ by rounding.
expect_p_clust_loop <- function(family, inp, cl, wdraws = inp$wdraws,
                                info = NULL) {
  out <- get_p_clust(family = family, eta = inp$eta, mu = inp$mu,
                     mu_offs = inp$mu_offs, dis = inp$dis, wdraws = wdraws,
                     cl = cl)
  ref <- p_clust_loop(family, mu = inp$mu, mu_offs = inp$mu_offs,
                      dis = inp$dis, wdraws = wdraws, cl = cl)
  one_draw <- all(tabulate(cl[!is.na(cl)], nbins = max(cl, na.rm = TRUE)) ==
                    1L)
  for (el in names(ref)) {
    if (one_draw) {
      expect_identical(unclass(out[[el]]), ref[[el]], info = paste(info, el))
    } else {
      # Arithmetic on NA can give NaN on some platforms, so compare the
      # positions of NA and NaN first:
      out_el <- unclass(out[[el]])
      ref_el <- ref[[el]]
      expect_equal(is.na(out_el), is.na(ref_el), info = paste(info, el, "NA"))
      expect_equal(out_el[!is.na(out_el)], ref_el[!is.na(ref_el)],
                   tolerance = 1e-12, info = paste(info, el))
    }
  }
}

p_clust_cases <- function(inp) {
  S <- ncol(inp$mu)
  thin <- rep(NA, S)
  thin[c(1, 5, 9, 20, 33)] <- 1:5
  w_special <- inp$wdraws
  w_special[c(2, 7, 11)] <- c(0, NA, Inf)
  mu_nonfinite <- inp
  mu_nonfinite$mu_offs[3, 4] <- NaN
  mu_nonfinite$mu_offs[6, 8] <- Inf
  w_zero_cluster <- inp$wdraws
  w_zero_cluster[seq(2, S, by = 8)] <- 0
  dis_na <- inp
  dis_na$dis[5] <- NA
  list(
    one_cluster = list(cl = rep(1, S)),
    several_draws_dis_na = list(cl = rep(1:8, length.out = S), inp = dis_na),
    all_draws = list(cl = seq_len(S)),
    permuted = list(cl = sample(S)),
    thinned = list(cl = thin),
    special_weights = list(cl = seq_len(S), w = w_special),
    nonfinite_mu = list(cl = seq_len(S), inp = mu_nonfinite),
    several_draws = list(cl = rep(1:8, length.out = S)),
    several_draws_thinned = list(cl = ifelse(seq_len(S) %% 3 == 0, NA,
                                             rep(1:6, length.out = S))),
    several_draws_special_weights = list(cl = rep(1:8, length.out = S),
                                         w = w_special),
    several_draws_zero_cluster = list(cl = rep(1:8, length.out = S),
                                      w = w_zero_cluster),
    several_draws_nonfinite_mu = list(cl = rep(1:8, length.out = S),
                                      inp = mu_nonfinite),
    empty_cluster = list(cl = c(1, rep(3, S - 1)))
  )
}

test_that("get_p_clust() gives the results of the cluster loop", {
  fams <- list(extend_family(gaussian()), extend_family(Student_t(nu = 4)),
               extend_family(binomial()), extend_family(poisson()))
  for (fam in fams) {
    inp <- p_clust_inputs(fam)
    cases <- p_clust_cases(inp)
    for (case_name in names(cases)) {
      cs <- cases[[case_name]]
      expect_p_clust_loop(fam, cs$inp %||% inp, cl = cs$cl,
                          wdraws = cs$w %||% inp$wdraws,
                          info = paste(fam$family, case_name))
    }
    # One draw, one observation:
    inp_1 <- p_clust_inputs(fam, n = 1, S = 1)
    expect_p_clust_loop(fam, inp_1, cl = 1L, info = paste(fam$family, "S 1"))
  }
})

test_that(paste(
  "get_p_clust() uses a replaced predvar() or a family without",
  "predvar_mix()"
), {
  inp <- p_clust_inputs(extend_family(gaussian()))
  S <- ncol(inp$mu)
  fam_repl <- extend_family(gaussian())
  fam_repl$predvar <- function(mu, dis, wdraws = 1) {
    rep(sum(wdraws * dis), NROW(mu))
  }
  fam_none <- extend_family(gaussian())
  fam_none$predvar_mix <- NULL
  w_special <- inp$wdraws
  w_special[c(2, 7, 11)] <- c(0, NA, Inf)
  for (cl in list(seq_len(S), sample(S), rep(1:8, length.out = S))) {
    for (w in list(inp$wdraws, w_special)) {
      expect_p_clust_loop(fam_repl, inp, cl = cl, wdraws = w,
                          info = "replaced")
      expect_p_clust_loop(fam_none, inp, cl = cl, wdraws = w, info = "none")
    }
  }
})

test_that("the built-in families give get_p_clust() its predvar_mix() path", {
  fams <- list(gaussian = extend_family(gaussian()),
               Student_t = extend_family(Student_t(nu = 4)),
               binomial = extend_family(binomial()),
               poisson = extend_family(poisson()),
               latent = extend_family(binomial(), latent = TRUE))
  for (fam in fams) {
    expect_true(is.function(fam$predvar_mix), info = fam$family)
    expect_identical(fam$predvar, attr(fam$predvar_mix, "predvar"),
                     info = fam$family)
  }
  uses_moments <- vapply(fams, function(fam) {
    attr(fam$predvar_mix, "uses_moments")
  }, logical(1))
  expect_identical(uses_moments, c(gaussian = TRUE, Student_t = TRUE,
                                   binomial = FALSE, poisson = FALSE,
                                   latent = TRUE))
  expect_null(extend_family(Gamma())$predvar_mix)
})

test_that("clust_moments_c() checks its inputs", {
  mu <- matrix(rnorm(12), 3, 4)
  w <- rep(0.5, 4)
  cl <- c(0L, 0L, 1L, 1L)
  expect_error(clust_moments_c(mu, mu[, -1], rep(1, 4), w, cl, 2L, TRUE),
               "dimensions")
  expect_error(clust_moments_c(mu, mu, rep(1, 4), w, cl[-1], 2L, TRUE),
               "dimensions")
  expect_error(clust_moments_c(mu, mu, rep(1, 4), w, c(0L, 0L, 1L, 2L), 2L,
                               TRUE),
               "out of range")
  expect_error(clust_moments_c(mu, mu, rep(1, 4), w, c(0L, NA, 1L, 1L), 2L,
                               TRUE),
               "out of range")
})

test_that(paste(
  "clust_moments_c() without the variances gives the same means and NA",
  "variances"
), {
  mu <- matrix(rnorm(15), 3, 5)
  mu_offs <- mu + 0.5
  dis <- runif(5, 0.5, 1.5)
  w <- c(0.2, 0.8, 0.5, 0.5, 0)
  cl <- c(0L, 0L, 1L, 1L, -1L)
  with_var <- clust_moments_c(mu, mu_offs, dis, w, cl, 2L, TRUE)
  without_var <- clust_moments_c(mu, mu_offs, dis, w, cl, 2L, FALSE)
  for (el in c("mu", "mu_offs", "dis")) {
    expect_identical(without_var[[el]], with_var[[el]], info = el)
  }
  expect_true(all(is.na(without_var$var_offs)))
  expect_identical(dim(without_var$var_offs), dim(with_var$var_offs))
  expect_true(all(is.na(without_var$dis2)))
  expect_length(without_var$dis2, 2L)
  expect_false(anyNA(with_var$var_offs))
  expect_false(anyNA(with_var$dis2))
})

test_that(paste(
  "get_p_clust() computes the variances in clust_moments_c() only if",
  "predvar_mix() uses them"
), {
  need_var_args <- logical()
  clust_moments_orig <- clust_moments_c
  fams <- list(gaussian = extend_family(gaussian()),
               Student_t = extend_family(Student_t(nu = 4)),
               latent = extend_family(binomial(), latent = TRUE),
               binomial = extend_family(binomial()),
               poisson = extend_family(poisson()))
  # A `predvar_mix()` without attribute "uses_moments" gets the variances:
  fam_noattr <- extend_family(gaussian())
  attr(fam_noattr$predvar_mix, "uses_moments") <- NULL
  fams$noattr <- fam_noattr
  # Without the predvar_mix() path (a replaced predvar(), or no predvar_mix()),
  # the cluster loop computes the predictive variances:
  fam_repl <- extend_family(gaussian())
  fam_repl$predvar <- function(mu, dis, wdraws = 1) {
    rep(sum(wdraws * dis), NROW(mu))
  }
  fams$replaced <- fam_repl
  fam_none <- extend_family(gaussian())
  fam_none$predvar_mix <- NULL
  fams$none <- fam_none
  need_var_exp <- c(gaussian = TRUE, Student_t = TRUE, latent = TRUE,
                    binomial = FALSE, poisson = FALSE, noattr = TRUE,
                    replaced = FALSE, none = FALSE)
  with_mocked_bindings(
    clust_moments_c = function(mu, mu_offs, dis, wn, cl, nclusters,
                               need_var) {
      need_var_args <<- c(need_var_args, need_var)
      clust_moments_orig(mu, mu_offs, dis, wn, cl, nclusters, need_var)
    },
    {
      for (fam_name in names(fams)) {
        fam <- fams[[fam_name]]
        inp <- p_clust_inputs(fam)
        S <- ncol(inp$mu)
        need_var_args <- logical()
        expect_p_clust_loop(fam, inp, cl = rep(1:8, length.out = S),
                            info = fam_name)
        expect_identical(need_var_args, unname(need_var_exp[fam_name]),
                         info = fam_name)
      }
    }
  )
})
