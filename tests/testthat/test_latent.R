# Setup -------------------------------------------------------------------

withr::local_preserve_seed()

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

test_that("cl_agg() gives the results of a loop over the clusters", {
  # The loop over the clusters (the reference for the computation in one pass):
  cl_agg_loop <- function(draws, cl = seq_len(nrow(draws)),
                          wdraws = rep(1, nrow(draws)), eps_wdraws = 0) {
    n_cl <- max(cl, na.rm = TRUE)
    do.call(rbind, lapply(seq_len(n_cl), function(i_cl) {
      idxs_draws_i <- which(cl == i_cl)
      wdraws_i <- wdraws[idxs_draws_i] / sum(wdraws[idxs_draws_i]) *
        (1 - eps_wdraws)
      wdraws_i %*% draws[idxs_draws_i, , drop = FALSE]
    }))
  }
  draws_nm <- draws
  dimnames(draws_nm) <- list(paste0("d", seq_len(S)), paste0("p", seq_len(P)))
  # No random numbers here (later tests use the random number state):
  draws_int <- matrix(rep_len(c(3L, 7L, 1L, 9L, 4L), S * P), nrow = S)
  draws_lgl <- draws_int > 4L
  draws_onecol <- draws[, 1, drop = FALSE]
  colnames(draws_onecol) <- "prec"
  perm <- order(sin(seq_len(S)))
  cl_list <- list(singleton = seq_len(S), perm = perm,
                  perm_dbl = as.numeric(perm), cl = cl_draws,
                  thin = th_draws,
                  empty = replace(cl_draws, cl_draws == 3L, 4L),
                  near_int = replace(as.numeric(seq_len(S)), 2, 2 + 1e-14))
  w_list <- list(const = rep(1, S), gamma = w_draws,
                 zero = replace(w_draws, 5, 0),
                 inf = replace(w_draws, 5, Inf), na = replace(w_draws, 5, NA))
  draws_list <- list(dbl = draws, named = draws_nm, int = draws_int,
                     lgl = draws_lgl, onecol = draws_onecol,
                     draws_matrix = posterior::as_draws_matrix(draws))
  for (cl_nm in names(cl_list)) {
    for (w_nm in names(w_list)) {
      for (d_nm in names(draws_list)) {
        for (eps in c(0, 0.1)) {
          args <- list(draws_list[[d_nm]], cl = cl_list[[cl_nm]],
                       wdraws = w_list[[w_nm]], eps_wdraws = eps)
          expect_identical(do.call(cl_agg, args), do.call(cl_agg_loop, args),
                           info = paste(cl_nm, w_nm, d_nm, eps))
        }
      }
    }
  }
  # Large cluster indices as doubles (match() compares numbers):
  draws_3 <- draws[1:3, ]
  expect_identical(cl_agg(draws_3, cl = c(1, 2, 1e5)),
                   cl_agg_loop(draws_3, cl = c(1, 2, 1e5)))
  expect_identical(cl_agg(draws_3, cl = c(0, 0, 0)),
                   cl_agg_loop(draws_3, cl = c(0, 0, 0)))
  expect_identical(cl_agg(draws_3, cl = c(0.5, 0.5, 0.2)),
                   cl_agg_loop(draws_3, cl = c(0.5, 0.5, 0.2)))
  expect_identical(cl_agg(draws_3, eps_wdraws = NA),
                   cl_agg_loop(draws_3, eps_wdraws = NA))
  expect_error(cl_agg(matrix(as.character(draws_3), nrow = 3)),
               "requires numeric/complex matrix/vector arguments")
  # Without weights (the old loop gives a matrix with zero rows):
  expect_identical(cl_agg(draws_3, wdraws = NULL),
                   cl_agg_loop(draws_3, wdraws = NULL))
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
  withr::local_preserve_seed()
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
    prjmat_crr_lat <- as.matrix(
      prj_crr_lat,
      allow_nonconst_wdraws_prj = ndr_ncl$clust_used_gt1
    )
    prjmat_crr_trad <- as.matrix(
      prj_crr_trad,
      allow_nonconst_wdraws_prj = ndr_ncl$clust_used_gt1
    )
    expect_identical(prjmat_crr_lat, prjmat_crr_trad, info = tstsetup_mod)
    prj_el_excl <- c("outdmin", "refmodel")
    expect_identical(prj_crr_lat[setdiff(names(prj_crr_lat), prj_el_excl)],
                     prj_crr_trad[setdiff(names(prj_crr_trad), prj_el_excl)],
                     info = tstsetup_mod)
  }
})

# Censored observations ---------------------------------------------------

test_that(paste0(
  "it is possible to handle censored observations in `latent_ll_oscale`"
), {
  tstsetups <- names(refmods[vapply(refmods, function(refmod_i) {
    refmod_i$family$for_latent
  }, FUN.VALUE = logical(1))])
  tstsetups <- head(tstsetups, 1)
  for (tstsetup in tstsetups) {
    refmod_i <- refmods[[tstsetup]]
    refmod_i$family$latent_ll_oscale <- structure(
      function(ilpreds,
               dis = rep(NA, nrow(ilpreds)),
               y_oscale,
               wobs = rep(1, ncol(ilpreds)),
               cens,
               cl_ref,
               wdraws_ref = rep(1, length(cl_ref))) {
        return(matrix(cens,
                      nrow = nrow(as.matrix(refmod_i$fit)),
                      ncol = length(cens),
                      byrow = TRUE))
      },
      cens_var = ~ dummy_cens_var
    )
    dat_cens_dummy <- head(refmod_i$fetch_data(), 3)
    dat_cens_dummy$dummy_cens_var <- rep(log(0.5), nrow(dat_cens_dummy))
    refpred_cens_dummy <- predict(refmod_i,
                                  newdata = dat_cens_dummy,
                                  ynew = rep(1, nrow(dat_cens_dummy)))
    expect_equal(refpred_cens_dummy, rep(log(0.5), nrow(dat_cens_dummy)))
  }
})

## Observation weights ----------------------------------------------------

test_that(paste(
  "the latent and the augmented-data projection stop for binomial trials > 1",
  "and name the workarounds"
), {
  set.seed(5301)
  n <- 30
  ndraws <- 20
  dat <- data.frame(x1 = rnorm(n), trials = sample(2:10, n, replace = TRUE))
  dat$y <- rbinom(n, dat$trials, plogis(0.5 * dat$x1))
  fit <- list(beta = rbind(rnorm(ndraws, 0, 0.1), rnorm(ndraws, 0.5, 0.1)),
              data = dat)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(~ x1, nd) %*% fit$beta
  }
  extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                 extract_y = TRUE) {
    y_wobs_offs(newdata = newdata, wrhs = wrhs %||% ~ trials, orhs = orhs,
                resp_form = if (extract_y) ~ y else NULL)
  }
  withr::local_options(projpred.warn_cvrefbuilder_NULL = FALSE)
  # Trials as observation weights (as in a brms model with `trials()`):
  expect_error(
    init_refmodel(
      fit, data = dat, formula = y ~ x1, family = binomial(),
      ref_predfun = ref_predfun, extract_model_data = extract_model_data,
      dis = rep(1, ndraws), latent = TRUE
    ),
    paste("latent projection may not be combined with observation weights",
          "\\(other than 1\\)\\. For the binomial family,",
          ".*\\(`latent = FALSE`\\).*de-aggregate")
  )
  # Trials in a 2-column response. The message does not name the traditional
  # projection, which reads only the first column with the default
  # `extract_model_data`:
  err_2col <- tryCatch(
    init_refmodel(
      fit, data = dat, formula = cbind(y, trials - y) ~ x1,
      family = binomial(), ref_predfun = ref_predfun, dis = rep(1, ndraws),
      latent = TRUE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err_2col, paste(
    "latent projection may not be combined with a 2-column response\\.",
    "As a workaround, de-aggregate"
  ))
  expect_no_match(err_2col, "traditional", fixed = TRUE)
  # The same for the augmented-data projection:
  err_2col_aug <- tryCatch(
    init_refmodel(
      fit, data = dat, formula = cbind(y, trials - y) ~ x1,
      family = binomial(), ref_predfun = ref_predfun,
      augdat_y_unqs = c("0", "1"), augdat_link = augdat_link_binom,
      augdat_ilink = augdat_ilink_binom
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err_2col_aug, paste(
    "augmented-data projection may not be combined with a 2-column",
    "response\\. As a workaround, de-aggregate"
  ))
  expect_no_match(err_2col_aug, "traditional", fixed = TRUE)
})

test_that(paste(
  "check_wobs_aug_lat() names the workarounds only for the binomial family"
), {
  wobs <- c(1, 3)
  fam_lat_binom <- extend_family(binomial(), latent = TRUE)
  fam_lat_poiss <- extend_family(poisson(), latent = TRUE)
  fam_aug_binom <- extend_family(binomial(),
                                 augdat_y_unqs = c("0", "1"),
                                 augdat_link = augdat_link_binom,
                                 augdat_ilink = augdat_ilink_binom)
  fam_trad_binom <- extend_family(binomial())
  expect_error(check_wobs_aug_lat(fam_lat_binom, wobs),
               "latent projection .*`latent = FALSE`.*de-aggregate")
  err_aug <- tryCatch(check_wobs_aug_lat(fam_aug_binom, wobs),
                      error = function(e) conditionMessage(e))
  expect_match(err_aug, "augmented-data projection .*traditional .*de-aggregate")
  expect_no_match(err_aug, "latent = FALSE", fixed = TRUE)
  err_poiss <- tryCatch(check_wobs_aug_lat(fam_lat_poiss, wobs),
                        error = function(e) conditionMessage(e))
  expect_match(err_poiss, "latent projection may not be combined",
               fixed = TRUE)
  expect_no_match(err_poiss, "de-aggregate", fixed = TRUE)
  expect_true(check_wobs_aug_lat(fam_lat_binom, c(1, 1)))
  expect_true(check_wobs_aug_lat(fam_trad_binom, wobs))
})

test_that(paste(
  "extend_family() stops for the brms families with `trials()` other than",
  "binomial(), in the latent and in the traditional projection"
), {
  for (fam_nm in c("beta_binomial", "zero_inflated_binomial",
                   "zero_inflated_beta_binomial")) {
    # A family object with the brms family name is sufficient, because the
    # check reads only `family$family`:
    fam <- binomial()
    fam$family <- fam_nm
    for (latent in c(FALSE, TRUE)) {
      expect_error(
        extend_family(fam, latent = latent),
        paste0("Currently, projpred does not support the brms family `",
               fam_nm, "()`."),
        fixed = TRUE, info = paste(fam_nm, "latent:", latent)
      )
    }
  }
})

test_that(paste(
  "get_refmodel() with `latent = TRUE` and without `dis` gives no message",
  "about the `NA` dispersion, and as.matrix.projection() has no `sigma` column"
), {
  skip_if_not_installed("rstanarm")
  set.seed(4521)
  dat <- data.frame(y = rpois(30, 2), x1 = rnorm(30), x2 = rnorm(30))
  fit <- suppressWarnings(rstanarm::stan_glm(
    y ~ x1 + x2, family = poisson(), data = dat, chains = 1, iter = 200,
    refresh = 0, seed = 7412
  ))
  expect_no_message(
    refm <- get_refmodel(fit, latent = TRUE),
    message = "will consist of only `NA`s"
  )
  expect_true(all(is.na(refm$dis)))
  # Without a latent dispersion, as.matrix() has no `sigma` column:
  prj <- project(refm, predictor_terms = "x1", ndraws = 20, verbose = 0)
  expect_false("sigma" %in% colnames(as.matrix(prj)))
  expect_false(
    "sigma" %in% posterior::variables(posterior::as_draws_matrix(prj))
  )
  # With a latent dispersion, it has one:
  refm_dis <- suppressMessages(
    get_refmodel(fit, latent = TRUE, dis = rep(1, length(refm$dis)))
  )
  prj_dis <- project(refm_dis, predictor_terms = "x1", ndraws = 20,
                     verbose = 0)
  expect_true("sigma" %in% colnames(as.matrix(prj_dis)))
})

## Expected values (`latent_epred`) ----------------------------------------

# A log-normal latent reference model with fixed draws. The output of
# `latent_ilink` is `meanlog`, and `latent_epred` (if used) gives the mean.
cre_ref_lnorm <- function(use_epred = TRUE, epred = NULL) {
  withr::local_seed(5302)
  n <- 40
  S <- 30
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dat$y <- rlnorm(n, meanlog = 0.5 * dat$x1, sdlog = 0.5)
  fit <- list(beta = rbind(rnorm(S, 0, 0.05), rnorm(S, 0.5, 0.05),
                           rnorm(S, 0, 0.05)),
              data = dat)
  sigma <- rgamma(S, shape = 50, rate = 100)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(~ x1 + x2, nd) %*% fit$beta
  }
  ilink <- function(lpreds, cl_ref, wdraws_ref = rep(1, length(cl_ref))) {
    lpreds
  }
  ll <- function(ilpreds, dis, y_oscale, wobs = rep(1, ncol(ilpreds)), cens,
                 cl_ref, wdraws_ref = rep(1, length(cl_ref))) {
    y_mat <- matrix(y_oscale, nrow(ilpreds), ncol(ilpreds), byrow = TRUE)
    dlnorm(y_mat, meanlog = ilpreds, sdlog = dis, log = TRUE)
  }
  ppd <- function(ilpreds_resamp, dis_resamp,
                  wobs = rep(1, ncol(ilpreds_resamp)), cl_ref,
                  wdraws_ref = rep(1, length(cl_ref)), idxs_prjdraws) {
    matrix(rlnorm(length(ilpreds_resamp), ilpreds_resamp, dis_resamp),
           nrow(ilpreds_resamp))
  }
  if (is.null(epred)) {
    epred <- function(ilpreds, dis, cl_ref,
                      wdraws_ref = rep(1, length(cl_ref))) {
      exp(ilpreds + dis^2 / 2)
    }
  }
  withr::local_options(projpred.warn_cvrefbuilder_NULL = FALSE)
  ref <- init_refmodel(
    fit, data = dat, formula = y ~ x1 + x2, family = gaussian(),
    ref_predfun = ref_predfun, dis = sigma, latent = TRUE,
    latent_ilink = ilink, latent_ll_oscale = ll, latent_ppd_oscale = ppd,
    latent_epred = if (use_epred) epred
  )
  list(ref = ref, dat = dat, sigma = sigma, fit = fit,
       ref_predfun = ref_predfun, eta = ref_predfun(fit))
}

# Muffles the Pareto k warnings of the small LOO CV runs below.
muffle_pareto <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Pareto", conditionMessage(w), ignore.case = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that(paste(
  "`latent_epred` gives the response-scale predictions of predict.refmodel(),",
  "proj_linpred() and proj_epred(), and does not change the log-likelihood"
), {
  r_e <- cre_ref_lnorm()
  r_p <- cre_ref_lnorm(use_epred = FALSE)
  # predict.refmodel():
  epred_ref <- exp(t(r_e$eta) + r_e$sigma^2 / 2)
  expect_equal(unname(predict(r_e$ref, type = "response")),
               unname(colMeans(epred_ref)))
  expect_equal(unname(predict(r_p$ref, type = "response")),
               unname(colMeans(t(r_p$eta))))
  expect_equal(predict(r_e$ref, ynew = r_e$dat$y),
               predict(r_p$ref, ynew = r_p$dat$y))
  # proj_linpred() and proj_epred(), with all draws:
  prj_e <- project(r_e$ref, predictor_terms = "x1", ndraws = 30,
                   verbose = FALSE)
  prj_p <- project(r_p$ref, predictor_terms = "x1", ndraws = 30,
                   verbose = FALSE)
  pl_lat <- proj_linpred(prj_e, transform = FALSE)$pred
  pl_e <- proj_linpred(prj_e, transform = TRUE)
  pl_p <- proj_linpred(prj_p, transform = TRUE)
  expect_equal(pl_e$pred, exp(pl_lat + prj_e$dis^2 / 2))
  expect_equal(pl_p$pred, pl_lat)
  expect_equal(pl_e$lpd, pl_p$lpd)
  expect_equal(proj_epred(prj_e), pl_e$pred)
  pl_e_int <- proj_linpred(prj_e, transform = TRUE, integrated = TRUE)$pred
  expect_equal(unname(c(pl_e_int)), unname(colMeans(pl_e$pred)))
  # With clusters, `latent_epred` gets the projected `dis` of the clusters:
  prj_c <- project(r_e$ref, predictor_terms = "x1", nclusters = 10, seed = 1,
                   verbose = FALSE)
  pl_lat_c <- proj_linpred(prj_c, transform = FALSE,
                           allow_nonconst_wdraws_prj = TRUE)$pred
  pl_e_c <- proj_linpred(prj_c, transform = TRUE,
                         allow_nonconst_wdraws_prj = TRUE)$pred
  expect_length(prj_c$dis, 10)
  expect_equal(unclass(pl_e_c), unclass(exp(pl_lat_c + prj_c$dis^2 / 2)))
})

test_that(paste(
  "`latent_epred` gives the expected values for `\"mse\"`, `\"rmse\"` and",
  "`\"R2\"` in varsel() and cv_varsel(), and does not change the lppd"
), {
  skip_on_cran()
  r_e <- cre_ref_lnorm()
  r_p <- cre_ref_lnorm(use_epred = FALSE)
  # With all draws (`ndraws_pred = 30`), the submodel projections are the
  # same as those of project():
  vs_e <- varsel(r_e$ref, nterms_max = 2, ndraws_pred = 30, seed = 1,
                 verbose = FALSE)
  vs_p <- varsel(r_p$ref, nterms_max = 2, ndraws_pred = 30, seed = 1,
                 verbose = FALSE)
  expect_equal(unname(vs_e$summaries$ref$oscale$mu),
               unname(colMeans(exp(t(r_e$eta) + r_e$sigma^2 / 2))))
  expect_equal(unname(vs_p$summaries$ref$oscale$mu),
               unname(colMeans(t(r_p$eta))))
  for (k in seq_along(vs_e$summaries$sub)) {
    info_k <- paste("size:", k - 1)
    sub_e <- vs_e$summaries$sub[[k]]$oscale
    sub_p <- vs_p$summaries$sub[[k]]$oscale
    expect_equal(sub_e$lppd, sub_p$lppd, info = info_k)
    prj_k <- project(vs_e, nterms = k - 1, ndraws = 30, verbose = FALSE)
    pl_lat_k <- proj_linpred(prj_k, transform = FALSE)$pred
    expect_equal(unname(sub_e$mu),
                 unname(colMeans(exp(pl_lat_k + prj_k$dis^2 / 2))),
                 info = info_k)
  }
  smmry_e <- summary(vs_e, stats = c("elpd", "mse"))$perf_sub
  smmry_p <- summary(vs_p, stats = c("elpd", "mse"))$perf_sub
  expect_equal(smmry_e$elpd, smmry_p$elpd)
  expect_false(isTRUE(all.equal(smmry_e$mse, smmry_p$mse)))
  # Test data (`d_test`):
  withr::local_seed(5303)
  dat_new <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  dat_new$y <- rlnorm(20, meanlog = 0.5 * dat_new$x1, sdlog = 0.5)
  eta_new <- r_e$ref_predfun(r_e$fit, newdata = dat_new)
  d_test <- list(data = dat_new, offset = rep(0, 20), weights = rep(1, 20),
                 y = rowMeans(eta_new), y_oscale = dat_new$y)
  vs_t <- varsel(r_e$ref, d_test = d_test, nterms_max = 1, ndraws_pred = 30,
                 seed = 1, verbose = FALSE)
  expect_equal(unname(vs_t$summaries$ref$oscale$mu),
               unname(colMeans(exp(t(eta_new) + r_e$sigma^2 / 2))))
  # LOO CV:
  for (vsel in c(FALSE, TRUE)) {
    cvvs_e <- muffle_pareto(cv_varsel(
      r_e$ref, cv_method = "LOO", validate_search = vsel, nterms_max = 2,
      nclusters_pred = 10, seed = 1, verbose = FALSE
    ))
    cvvs_p <- muffle_pareto(cv_varsel(
      r_p$ref, cv_method = "LOO", validate_search = vsel, nterms_max = 2,
      nclusters_pred = 10, seed = 1, verbose = FALSE
    ))
    info_vsel <- paste("validate_search:", vsel)
    expect_equal(cvvs_e$summaries$ref$oscale$lppd,
                 cvvs_p$summaries$ref$oscale$lppd, info = info_vsel)
    expect_true(all(cvvs_e$summaries$ref$oscale$mu >
                      cvvs_p$summaries$ref$oscale$mu), info = info_vsel)
    for (k in seq_along(cvvs_e$summaries$sub)) {
      sub_e <- cvvs_e$summaries$sub[[k]]$oscale
      sub_p <- cvvs_p$summaries$sub[[k]]$oscale
      info_k <- paste(info_vsel, "size:", k - 1)
      expect_equal(sub_e$lppd, sub_p$lppd, info = info_k)
      expect_true(all(sub_e$mu > sub_p$mu), info = info_k)
    }
  }
})

test_that("extend_family() and the predictions check `latent_epred`", {
  expect_error(
    extend_family(binomial(), latent = TRUE, latent_y_unqs = c("0", "1"),
                  latent_epred = function(ilpreds, dis, cl_ref, wdraws_ref) {
                    ilpreds
                  }),
    "Argument `latent_epred` cannot be used if `family$cats` is not `NULL`.",
    fixed = TRUE
  )
  expect_error(
    extend_family(poisson(), latent = TRUE, latent_epred = TRUE),
    "Argument `latent_epred` needs to be a function or `NULL`.", fixed = TRUE
  )
  expect_error(
    extend_family(poisson(), latent = TRUE,
                  latent_epred = function(ilpreds, cl_ref) ilpreds),
    "Please add argument(s) `dis`, `wdraws_ref` to the `latent_epred`",
    fixed = TRUE
  )
  expect_s3_class(
    extend_family(poisson(), latent = TRUE,
                  latent_epred = function(ilpreds, ...) ilpreds),
    "family"
  )
  r_bad <- cre_ref_lnorm(epred = function(ilpreds, dis, cl_ref, wdraws_ref) {
    ilpreds[, -1, drop = FALSE]
  })
  expect_error(predict(r_bad$ref, type = "response"),
               "Unexpected structure for the output of `latent_epred`.",
               fixed = TRUE)
  # Only `NA`s from `latent_epred` (with a finite `latent_ilink` output) give
  # the message of proj_linpred():
  r_na <- cre_ref_lnorm(epred = function(ilpreds, dis, cl_ref, wdraws_ref) {
    ilpreds * NA
  })
  prj_na <- suppressMessages(
    project(r_na$ref, predictor_terms = "x1", ndraws = 30, verbose = FALSE)
  )
  expect_message(
    pl_na <- proj_linpred(prj_na, transform = TRUE),
    "`latent_ilink` (or `latent_epred`) returned only `NA`s", fixed = TRUE
  )
  expect_true(all(is.na(pl_na$pred)))
  expect_false(anyNA(pl_na$lpd))
})
