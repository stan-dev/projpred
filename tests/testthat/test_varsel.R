# varsel() ----------------------------------------------------------------

context("varsel()")

test_that(paste(
  "`object` of class \"refmodel\", `method`, `nterms_max`, `nclusters`, and",
  "`nclusters_pred` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(vss)) {
    tstsetup_ref <- args_vs[[tstsetup]]$tstsetup_ref
    mod_crr <- args_vs[[tstsetup]]$mod_nm
    fam_crr <- args_vs[[tstsetup]]$fam_nm
    prj_crr <- args_vs[[tstsetup]]$prj_nm
    meth_exp_crr <- args_vs[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    vsel_tester(
      vss[[tstsetup]],
      refmod_expected = refmods[[tstsetup_ref]],
      solterms_len_expected = args_vs[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      search_trms_empty_size =
        length(args_vs[[tstsetup]]$search_terms) &&
        all(grepl("\\+", args_vs[[tstsetup]]$search_terms)),
      info_str = tstsetup
    )
  }
})

test_that("invalid `object` fails", {
  expect_error(varsel(rnorm(5), verbose = FALSE),
               "no applicable method")
})

test_that("invalid `method` fails", {
  for (tstsetup in names(refmods)) {
    expect_error(varsel(refmods[[tstsetup]], method = "k-fold"),
                 "Unknown search method",
                 info = tstsetup)
    if (args_ref[[tstsetup]]$mod_nm != "glm") {
      expect_error(varsel(refmods[[tstsetup]], method = "L1"),
                   paste("^L1 search is only supported for reference models",
                         "without multilevel and without additive",
                         "\\(\"smoothing\"\\) terms\\.$"),
                   info = tstsetup)
    }
    if (args_ref[[tstsetup]]$mod_nm == "glm" &&
        args_ref[[tstsetup]]$prj_nm == "augdat") {
      expect_error(varsel(refmods[[tstsetup]], method = "L1"),
                   paste("^Currently, the augmented-data projection may not be",
                         "combined with an L1 search\\.$"),
                   info = tstsetup)
    }
  }
})

test_that("`seed` works (and restores the RNG state afterwards)", {
  skip_if_not(run_vs)
  # To save time:
  tstsetups <- grep("\\.glm\\.gauss\\.", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    vs_orig <- vss[[tstsetup]]
    rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_repr1 <- .Random.seed
    vs_repr <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$tstsetup_ref]]),
      excl_nonargs(args_vs_i)
    ))
    .Random.seed_repr2 <- .Random.seed
    rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
    # Expected equality:
    expect_equal(vs_repr, vs_orig, info = tstsetup)
    expect_equal(.Random.seed_repr2, .Random.seed_repr1, info = tstsetup)
    # Expected inequality:
    expect_false(isTRUE(all.equal(rand_new, rand_orig)), info = tstsetup)
  }
})

## d_test -----------------------------------------------------------------

test_that(paste(
  "`d_test` set to the training data gives the same results as its default"
), {
  skip_if_not(run_vs)
  tstsetups <- names(vss)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    tstsetup_ref <- args_vs_i$tstsetup_ref
    pkg_crr <- args_vs_i$pkg_nm
    mod_crr <- args_vs_i$mod_nm
    fam_crr <- args_vs_i$fam_nm
    prj_crr <- args_vs_i$prj_nm
    if (!all(refmods[[tstsetup_ref]]$offset == 0)) {
      offs_crr <- offs_tst
    } else {
      offs_crr <- rep(0, nobsv)
    }
    if (!all(refmods[[tstsetup_ref]]$wobs == 1)) {
      wobs_crr <- wobs_tst
    } else {
      wobs_crr <- rep(1, nobsv)
    }
    formul_fit_crr <- args_fit[[args_vs_i$tstsetup_fit]]$formula
    dat_crr <- get_dat_formul(formul_crr = formul_fit_crr,
                              needs_adj = grepl("\\.spclformul", tstsetup))
    d_test_crr <- list(
      data = dat,
      offset = offs_crr,
      weights = wobs_crr,
      y = dat_crr[[stdize_lhs(formul_fit_crr)$y_nm]]
    )
    y_oscale_crr <- d_test_crr$y
    if (prj_crr %in% c("latent", "augdat")) {
      if (use_fac) {
        yunqs <- yunq_chr
      } else {
        yunqs <- as.character(yunq_num)
      }
      if (!(prj_crr == "latent" && fam_crr == "brnll")) {
        lvls_crr <- args_ref[[args_vs_i$tstsetup_ref]]$augdat_y_unqs
        lvls_crr <- lvls_crr %||%
          args_ref[[args_vs_i$tstsetup_ref]]$latent_y_unqs
        lvls_crr <- lvls_crr %||% yunqs
        y_oscale_crr <- factor(as.character(y_oscale_crr), levels = lvls_crr,
                               ordered = is.ordered(y_oscale_crr))
      }
      if (prj_crr == "augdat") {
        d_test_crr$y <- y_oscale_crr
      } else if (prj_crr == "latent") {
        d_test_crr$y <- colMeans(
          unname(posterior_linpred(fits[[args_vs_i$tstsetup_fit]]))
        )
      }
    }
    d_test_crr$y_oscale <- y_oscale_crr
    if (prj_crr == "augdat" && fam_crr == "cumul") {
      warn_expected <- "non-integer #successes in a binomial glm!"
    } else if (!is.null(args_vs_i$avoid.increase)) {
      warn_expected <- warn_mclogit
    } else {
      warn_expected <- NA
    }
    expect_warning(
      vs_repr <- do.call(varsel, c(
        list(object = refmods[[tstsetup_ref]], d_test = d_test_crr),
        excl_nonargs(args_vs_i)
      )),
      warn_expected
    )
    meth_exp_crr <- args_vs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    vsel_tester(
      vs_repr,
      refmod_expected = refmods[[tstsetup_ref]],
      ywtest_expected = setNames(
        as.data.frame(d_test_crr[nms_y_wobs_test(wobs_nm = "weights")]),
        nms_y_wobs_test()
      ),
      solterms_len_expected = args_vs_i$nterms_max,
      method_expected = meth_exp_crr,
      search_trms_empty_size =
        length(args_vs_i$search_terms) &&
        all(grepl("\\+", args_vs_i$search_terms)),
      info_str = tstsetup
    )
    expect_equal(vs_repr[setdiff(names(vs_repr),
                                 c("type_test", "y_wobs_test"))],
                 vss[[tstsetup]][setdiff(names(vss[[tstsetup]]),
                                         c("type_test", "y_wobs_test"))],
                 info = tstsetup)
    y_wobs_test_orig <- vss[[tstsetup]]$y_wobs_test
    if (pkg_crr == "brms") {
      # brms seems to set argument `contrasts`, but this is not important for
      # projpred, so ignore it in the comparison:
      attr(y_wobs_test_orig$y, "contrasts") <- NULL
      attr(y_wobs_test_orig$y_oscale, "contrasts") <- NULL
    }
    expect_equal(vs_repr$y_wobs_test, y_wobs_test_orig, info = tstsetup)
  }
})

test_that(paste(
  "`d_test` set to actual test data gives a `<vsel_object>$summaries$sub`",
  "object that can be reproduced by proj_linpred() and a",
  "`<vsel_object>$summaries$ref` object that can be reproduced by",
  "posterior_epred() and log_lik()"
), {
  skip_if_not(run_vs)
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  tstsetups <- names(vss)
  ### TODO (GAMMs): Currently, the following test setups (can) lead to the error
  ### ```
  ### Error in t(as.matrix(b$reTrms$Zt[ii, ])) %*%
  ### as.matrix(c(as.matrix(ranef[[i]]))) :
  ###   non-conformable arguments
  ### ```
  ### thrown by predict.gamm4(). This needs to be fixed. For now, exclude these
  ### test setups:
  tstsetups <- grep("\\.gamm\\.", tstsetups, value = TRUE, invert = TRUE)
  ###
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    tstsetup_ref <- args_vs_i$tstsetup_ref
    pkg_crr <- args_vs_i$pkg_nm
    mod_crr <- args_vs_i$mod_nm
    fam_crr <- args_vs_i$fam_nm
    prj_crr <- args_vs_i$prj_nm
    if (!all(refmods[[tstsetup_ref]]$offset == 0)) {
      offs_crr <- offs_indep
    } else {
      offs_crr <- rep(0, nobsv_indep)
    }
    if (!all(refmods[[tstsetup_ref]]$wobs == 1)) {
      wobs_crr <- wobs_indep
    } else {
      wobs_crr <- rep(1, nobsv_indep)
    }
    formul_fit_crr <- args_fit[[args_vs_i$tstsetup_fit]]$formula
    y_nm_crr <- stdize_lhs(formul_fit_crr)$y_nm
    dat_indep_crr <- get_dat_formul(
      formul_crr = formul_fit_crr,
      needs_adj = grepl("\\.spclformul", tstsetup),
      dat_crr = dat_indep
    )
    d_test_crr <- list(
      data = dat_indep,
      offset = offs_crr,
      weights = wobs_crr,
      y = dat_indep_crr[[y_nm_crr]]
    )
    y_oscale_crr <- d_test_crr$y
    if (prj_crr %in% c("latent", "augdat")) {
      if (use_fac) {
        yunqs <- yunq_chr
      } else {
        yunqs <- as.character(yunq_num)
      }
      if (!(prj_crr == "latent" && fam_crr == "brnll")) {
        lvls_crr <- args_ref[[args_vs_i$tstsetup_ref]]$augdat_y_unqs
        lvls_crr <- lvls_crr %||%
          args_ref[[args_vs_i$tstsetup_ref]]$latent_y_unqs
        lvls_crr <- lvls_crr %||% yunqs
        y_oscale_crr <- factor(as.character(y_oscale_crr), levels = lvls_crr,
                               ordered = is.ordered(y_oscale_crr))
      }
      if (prj_crr == "augdat") {
        d_test_crr$y <- y_oscale_crr
      } else if (prj_crr == "latent") {
        if (pkg_crr == "rstanarm") {
          post_linpred <- posterior_linpred(fits[[args_vs_i$tstsetup_fit]],
                                            newdata = dat_indep,
                                            offset = d_test_crr$offset)
        } else {
          post_linpred <- posterior_linpred(fits[[args_vs_i$tstsetup_fit]],
                                            newdata = dat_indep)
        }
        d_test_crr$y <- colMeans(unname(post_linpred))
      }
    }
    d_test_crr$y_oscale <- y_oscale_crr
    if (prj_crr == "augdat" && fam_crr == "cumul") {
      warn_expected <- "non-integer #successes in a binomial glm!"
    } else if (!is.null(args_vs_i$avoid.increase)) {
      warn_expected <- warn_mclogit
    } else {
      warn_expected <- NA
    }
    expect_warning(
      vs_indep <- do.call(varsel, c(
        list(object = refmods[[tstsetup_ref]], d_test = d_test_crr),
        excl_nonargs(args_vs_i)
      )),
      warn_expected
    )
    meth_exp_crr <- args_vs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    vsel_tester(
      vs_indep,
      refmod_expected = refmods[[tstsetup_ref]],
      ywtest_expected = setNames(
        as.data.frame(d_test_crr[nms_y_wobs_test(wobs_nm = "weights")]),
        nms_y_wobs_test()
      ),
      solterms_len_expected = args_vs_i$nterms_max,
      method_expected = meth_exp_crr,
      search_trms_empty_size =
        length(args_vs_i$search_terms) &&
        all(grepl("\\+", args_vs_i$search_terms)),
      info_str = tstsetup
    )

    ### Summaries for the submodels -------------------------------------------

    if (!(getOption("projpred.mlvl_pred_new", FALSE) &&
          mod_crr %in% c("glmm", "gamm") &&
          any(grepl("\\|", vs_indep$solution_terms)))) {
      # In the negation of this case (i.e., multilevel models with option
      # `projpred.mlvl_pred_new` being set to `TRUE`), proj_linpred() can't be
      # used to calculate the reference model's performance statistics because
      # proj_linpred()'s argument `.seed` cannot be set such that the
      # .Random.seed from inside proj_linpred() at the place where the new
      # group-level effects are drawn coincides with .Random.seed from inside
      # varsel() at the place where the new group-level effects are drawn (not
      # even `.seed = NA` with an appropriate preparation is possible).

      if (!is.null(args_vs_i$avoid.increase)) {
        warn_expected <- NA
      }
      # For getting the correct seed in proj_linpred():
      set.seed(args_vs_i$seed)
      p_sel_dummy <- get_refdist(refmods[[tstsetup_ref]],
                                 nclusters = vs_indep$nprjdraws_search)
      expect_warning(
        pl_indep <- proj_linpred(
          vs_indep,
          newdata = dat_indep_crr,
          offsetnew = d_test_crr$offset,
          weightsnew = d_test_crr$weights,
          transform = TRUE,
          integrated = TRUE,
          .seed = NA,
          nterms = c(0L, seq_along(vs_indep$solution_terms)),
          nclusters = args_vs_i$nclusters_pred,
          seed = NA
        ),
        warn_expected
      )
      summ_sub_ch <- lapply(pl_indep, function(pl_indep_k) {
        names(pl_indep_k)[names(pl_indep_k) == "pred"] <- "mu"
        names(pl_indep_k)[names(pl_indep_k) == "lpd"] <- "lppd"
        pl_indep_k$mu <- unname(drop(pl_indep_k$mu))
        pl_indep_k$lppd <- drop(pl_indep_k$lppd)
        if (!is.null(refmods[[tstsetup_ref]]$family$cats)) {
          pl_indep_k$mu <- structure(as.vector(pl_indep_k$mu),
                                     class = "augvec",
                                     nobs_orig = nrow(pl_indep_k$mu))
        }
        return(pl_indep_k)
      })
      if (prj_crr == "latent") {
        # For getting the correct seed in proj_linpred():
        set.seed(args_vs_i$seed)
        p_sel_dummy <- get_refdist(refmods[[tstsetup_ref]],
                                   nclusters = vs_indep$nprjdraws_search)
        dat_indep_crr[[paste0(".", y_nm_crr)]] <- d_test_crr$y
        pl_indep_lat <- proj_linpred(
          vs_indep,
          newdata = dat_indep_crr,
          offsetnew = d_test_crr$offset,
          weightsnew = d_test_crr$weights,
          transform = FALSE,
          integrated = TRUE,
          .seed = NA,
          nterms = c(0L, seq_along(vs_indep$solution_terms)),
          nclusters = args_vs_i$nclusters_pred,
          seed = NA
        )
        y_lat_mat <- matrix(d_test_crr$y, nrow = args_vs_i$nclusters_pred,
                            ncol = nobsv_indep, byrow = TRUE)
        summ_sub_ch_lat <- lapply(seq_along(pl_indep_lat), function(k_idx) {
          pl_indep_k <- pl_indep_lat[[k_idx]]
          names(pl_indep_k)[names(pl_indep_k) == "pred"] <- "mu"
          names(pl_indep_k)[names(pl_indep_k) == "lpd"] <- "lppd"
          pl_indep_k$mu <- unname(drop(pl_indep_k$mu))
          pl_indep_k$lppd <- drop(pl_indep_k$lppd)
          return(pl_indep_k)
        })
        summ_sub_ch <- lapply(seq_along(summ_sub_ch), function(k_idx) {
          c(summ_sub_ch_lat[[k_idx]], list("oscale" = summ_sub_ch[[k_idx]]))
        })
      }
      names(summ_sub_ch) <- NULL
      expect_equal(vs_indep$summaries$sub, summ_sub_ch,
                   tolerance = .Machine$double.eps, info = tstsetup)
    }

    ### Summaries for the reference model -------------------------------------

    if (getOption("projpred.mlvl_pred_new", FALSE)) {
      dat_indep_crr$z.1 <- as.factor(paste0("NEW_", dat_indep_crr$z.1))
    }
    if (pkg_crr == "rstanarm") {
      mu_new <- rstantools::posterior_epred(refmods[[tstsetup_ref]]$fit,
                                            newdata = dat_indep_crr,
                                            offset = d_test_crr$offset)
      if (fam_crr == "cumul") {
        eta_new <- rstantools::posterior_linpred(refmods[[tstsetup_ref]]$fit,
                                                 newdata = dat_indep_crr,
                                                 offset = d_test_crr$offset)
        # The following shows that in case of an rstanarm::stan_polr() fit,
        # rstantools::posterior_epred() returns the linear predictors with a
        # threshold of zero, transformed to response scale (which is not really
        # helpful):
        mu_new_ch <- augdat_ilink_cumul(
          array(eta_new, dim = c(dim(eta_new), 1L)),
          link = link_str
        )
        stopifnot(isTRUE(all.equal(unname(mu_new), mu_new_ch[, , 1],
                                   tolerance = .Machine$double.eps)))
        # Therefore, `mu_new` has to be adapted to incorporate the correct
        # thresholds:
        drws <- as.matrix(refmods[[tstsetup_ref]]$fit)
        drws_thres <- drws[, grep("\\|", colnames(drws))]
        mu_new <- apply(drws_thres, 2, function(thres_vec) {
          thres_vec - eta_new
        }, simplify = FALSE)
        mu_new <- abind::abind(mu_new, rev.along = 0)
        mu_new <- augdat_ilink_cumul(mu_new, link = link_str)
      }
      if (grepl("\\.without_wobs", tstsetup)) {
        lppd_new <- rstantools::log_lik(refmods[[tstsetup_ref]]$fit,
                                        newdata = dat_indep_crr,
                                        offset = d_test_crr$offset)
      } else {
        # Currently, rstanarm issue #567 causes an error to be thrown when
        # calling log_lik(). Therefore, use the following dummy which guarantees
        # test success:
        lppd_new <- matrix(vs_indep$summaries$ref$lppd,
                           nrow = nrefdraws, ncol = nobsv_indep, byrow = TRUE)
      }
      if (prj_crr == "latent") {
        mu_new_lat <- rstantools::posterior_linpred(refmods[[tstsetup_ref]]$fit,
                                                    newdata = dat_indep_crr,
                                                    offset = d_test_crr$offset)
      }
    } else if (pkg_crr == "brms") {
      expr_seed <- expression({
        set.seed(seed2_tst)
        kfold_seed_dummy <- sample.int(.Machine$integer.max, 1)
        refprd_seed_dummy <- sample.int(.Machine$integer.max, 1)
        set.seed(refprd_seed_dummy)
      })
      eval(expr_seed)
      mu_new <- rstantools::posterior_epred(refmods[[tstsetup_ref]]$fit,
                                            newdata = dat_indep_crr,
                                            allow_new_levels = TRUE,
                                            sample_new_levels = "gaussian")
      if (fam_crr == "binom") {
        # Compared to rstanarm, brms uses a different convention for the
        # binomial family: The values returned by posterior_epred() are not
        # probabilities, but the expected values on the scale of the response
        # (so the probabilities multiplied by the number of trials). Thus, we
        # have to revert this here:
        mu_new <- mu_new / matrix(wobs_indep, nrow = nrow(mu_new),
                                  ncol = ncol(mu_new), byrow = TRUE)
      }
      eval(expr_seed)
      lppd_new <- rstantools::log_lik(refmods[[tstsetup_ref]]$fit,
                                      newdata = dat_indep_crr,
                                      allow_new_levels = TRUE,
                                      sample_new_levels = "gaussian")
      if (prj_crr == "latent") {
        eval(expr_seed)
        mu_new_lat <- rstantools::posterior_linpred(
          refmods[[tstsetup_ref]]$fit,
          newdata = dat_indep_crr,
          allow_new_levels = TRUE,
          sample_new_levels = "gaussian"
        )
      }
    }
    if (length(dim(mu_new)) == 2) {
      mu_new <- colMeans(mu_new)
    } else if (length(dim(mu_new)) == 3) {
      # In fact, we have `identical(colMeans(mu_new), apply(mu_new, c(2, 3),
      # mean))` giving `TRUE`, but it's better to be explicit:
      mu_new <- apply(mu_new, c(2, 3), mean)
      mu_new <- structure(as.vector(mu_new), class = "augvec",
                          nobs_orig = nobsv_indep)
    } else {
      stop("Unexpected number of margins for `mu_new`.")
    }
    summ_ref_ch <- list(
      mu = unname(mu_new),
      lppd = unname(apply(lppd_new, 2, log_sum_exp) - log(nrefdraws))
    )
    if (prj_crr == "augdat" && fam_crr %in% c("brnll", "binom")) {
      summ_ref_ch$mu <- structure(c(1 - summ_ref_ch$mu, summ_ref_ch$mu),
                                  class = "augvec",
                                  nobs_orig = length(summ_ref_ch$mu))
    }
    if (prj_crr == "latent") {
      y_lat_mat <- matrix(d_test_crr$y, nrow = nrefdraws, ncol = nobsv_indep,
                          byrow = TRUE)
      lppd_new_lat <- dnorm(y_lat_mat, mean = mu_new_lat,
                            sd = refmods[[tstsetup_ref]]$dis, log = TRUE)
      summ_ref_ch_lat <- list(
        mu = unname(colMeans(mu_new_lat)),
        lppd = unname(apply(lppd_new_lat, 2, log_sum_exp) - log(nrefdraws))
      )
      summ_ref_ch <- c(summ_ref_ch_lat, list("oscale" = summ_ref_ch))
    }
    expect_equal(vs_indep$summaries$ref, summ_ref_ch,
                 tolerance = 1e3 * .Machine$double.eps, info = tstsetup)
    lppd_ref_ch2 <- unname(loo::elpd(lppd_new)$pointwise[, "elpd"])
    if (prj_crr == "latent") {
      lppd_ref_ch2_oscale <- lppd_ref_ch2
      expect_equal(vs_indep$summaries$ref$oscale$lppd, lppd_ref_ch2_oscale,
                   tolerance = 1e1 * .Machine$double.eps, info = tstsetup)
      lppd_ref_ch2 <- loo::elpd(lppd_new_lat)$pointwise[, "elpd"]
    }
    expect_equal(vs_indep$summaries$ref$lppd, lppd_ref_ch2,
                 tolerance = 1e2 * .Machine$double.eps, info = tstsetup)
  }
  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})

## refit_prj --------------------------------------------------------------

test_that("`refit_prj` works", {
  skip_if_not(run_vs)
  if (run_more) {
    tstsetups <- names(vss)
  } else {
    tstsetups <- head(grep("\\.glm\\.", names(vss), value = TRUE), 1)
  }
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    args_vs_i$refit_prj <- FALSE
    if (args_vs_i$prj_nm == "augdat" && args_vs_i$fam_nm == "cumul") {
      warn_expected <- "non-integer #successes in a binomial glm!"
    } else if (!is.null(args_vs_i$avoid.increase)) {
      warn_expected <- warn_mclogit
    } else {
      warn_expected <- NA
    }
    expect_warning(
      vs_reuse <- do.call(varsel, c(
        list(object = refmods[[args_vs_i$tstsetup_ref]]),
        excl_nonargs(args_vs_i)
      )),
      warn_expected,
      info = tstsetup
    )
    mod_crr <- args_vs_i$mod_nm
    fam_crr <- args_vs_i$fam_nm
    prj_crr <- args_vs_i$prj_nm
    meth_exp_crr <- args_vs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    extra_tol_crr <- 1.1
    if (meth_exp_crr == "L1" &&
        any(grepl(":", ranking(vs_reuse)[["fulldata"]]))) {
      ### Testing for non-increasing element `ce` (for increasing model size)
      ### doesn't make sense if the ranking of predictors involved in
      ### interactions has been changed, so we choose a higher `extra_tol`:
      extra_tol_crr <- 1.2
      ###
    }
    vsel_tester(
      vs_reuse,
      refmod_expected = refmods[[args_vs_i$tstsetup_ref]],
      solterms_len_expected = args_vs_i$nterms_max,
      method_expected = meth_exp_crr,
      refit_prj_expected = FALSE,
      search_trms_empty_size =
        length(args_vs_i$search_terms) &&
        all(grepl("\\+", args_vs_i$search_terms)),
      extra_tol = extra_tol_crr,
      info_str = tstsetup
    )
  }
})

## Regularization ---------------------------------------------------------

# In fact, `regul` is already checked in `test_project.R`, so the `regul` tests
# could be omitted here since varsel() and cv_varsel() also pass `regul` to
# get_submodl_prj() (usually via get_submodls(), just like project()). This
# doesn't hold for L1 search, though. So for L1 search, the `regul` tests are
# still needed.

test_that(paste(
  "for GLMs with L1 search, `regul` only has an effect on prediction, not on",
  "selection"
), {
  skip_if_not(run_vs)
  regul_tst <- c(regul_default, 1e-1, 1e2)
  stopifnot(regul_tst[1] == regul_default)
  stopifnot(all(diff(regul_tst) > 0))
  tstsetups <- setdiff(
    setdiff(grep("\\.glm\\.", names(vss), value = TRUE),
            grep("\\.glm\\..*\\.forward", names(vss), value = TRUE)),
    grep("\\.glm\\..*\\.augdat\\.", names(vss), value = TRUE)
  )
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    m_max <- args_vs_i$nterms_max + 1L
    ssq_regul_prd <- array(dim = c(length(regul_tst), m_max))
    for (j in seq_along(regul_tst)) {
      if (regul_tst[j] == regul_default) {
        vs_regul <- vss[[tstsetup]]
      } else {
        vs_regul <- do.call(varsel, c(
          list(object = refmods[[args_vs_i$tstsetup_ref]],
               regul = regul_tst[j]),
          excl_nonargs(args_vs_i)
        ))
        vsel_tester(
          vs_regul,
          refmod_expected = refmods[[args_vs_i$tstsetup_ref]],
          solterms_len_expected = args_vs_i$nterms_max,
          method_expected = "L1",
          info_str = tstsetup
        )
        # Expect equality for all components not related to prediction:
        expect_equal(vs_regul[setdiff(vsel_nms, vsel_nms_pred)],
                     vss[[tstsetup]][setdiff(vsel_nms, vsel_nms_pred)],
                     info = paste(tstsetup, j, sep = "__"))
        # Expect inequality for the components related to prediction (but note
        # that the components from `vsel_nms_pred_opt` can be, but don't need to
        # be differing):
        for (vsel_nm in setdiff(vsel_nms_pred, vsel_nms_pred_opt)) {
          expect_false(isTRUE(all.equal(vs_regul[[vsel_nm]],
                                        vss[[tstsetup]][[vsel_nm]])),
                       info = paste(tstsetup, j, vsel_nm, sep = "__"))
        }
      }
      # Check the inequality of the prediction components in detail: Expect a
      # reduction of the sum of the squared coefficients (excluding the
      # intercept) for increasing `regul`:
      for (m in seq_len(m_max)) {
        # Since varsel() doesn't output object `p_sub`, use the linear predictor
        # here (instead of the coefficients themselves, which would only be
        # accessible from `p_sub`):
        mu_jm_regul <- vs_regul$refmodel$family$linkfun(
          vs_regul$summaries$sub[[m]]$mu
        )
        if (grepl("\\.with_offs", tstsetup)) {
          mu_jm_regul <- mu_jm_regul - offs_tst
        }
        # In fact, `sum((mu - offset - intercept)^2)` would make more sense than
        # `var(mu - offset) = sum((mu - offset - mean(mu - offset))^2)` but
        # since varsel() doesn't output object `p_sub`, the intercept from the
        # prediction is not accessible here.
        ssq_regul_prd[j, m] <- var(mu_jm_regul)
      }
    }
    # For the intercept-only model, the linear predictor consists only
    # of the intercept, so we expect no variation in `mu_jm_regul`:
    expect_true(all(ssq_regul_prd[, 1] <= 1e-5), info = tstsetup)
    # All other (i.e., not intercept-only) models:
    for (j in seq_len(dim(ssq_regul_prd)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_prd)[2])[-1]) {
        expect_lt(ssq_regul_prd[!!j, !!m], ssq_regul_prd[j - 1, m])
      }
    }
  }
})

test_that(paste(
  "for GLMs with forward search, `regul` has an expected effect on selection",
  "as well as on prediction"
), {
  skip_if_not(run_vs)
  regul_tst <- c(regul_default, 1e-1, 1e2)
  stopifnot(regul_tst[1] == regul_default)
  stopifnot(all(diff(regul_tst) > 0))
  tstsetups <- union(grep("\\.glm\\..*\\.forward", names(vss), value = TRUE),
                     grep("\\.glm\\..*\\.augdat\\.", names(vss), value = TRUE))
  tstsetups <- grep(fam_nms_aug_regex, tstsetups, value = TRUE, invert = TRUE)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    m_max <- args_vs_i$nterms_max + 1L
    if (length(args_vs_i$search_terms) &&
        all(grepl("\\+", args_vs_i$search_terms))) {
      # This is the "empty_size" setting, so we have to subtract the skipped
      # model size (see issue #307):
      m_max <- m_max - 1L
    }
    ncl_crr <- args_vs_i$nclusters
    ssq_regul_sel_alpha <- array(dim = c(length(regul_tst), m_max, ncl_crr))
    ssq_regul_sel_beta <- array(dim = c(length(regul_tst), m_max, ncl_crr))
    ssq_regul_prd <- array(dim = c(length(regul_tst), m_max))
    for (j in seq_along(regul_tst)) {
      if (regul_tst[j] == regul_default) {
        vs_regul <- vss[[tstsetup]]
      } else {
        vs_regul <- do.call(varsel, c(
          list(object = refmods[[args_vs_i$tstsetup_ref]],
               regul = regul_tst[j]),
          excl_nonargs(args_vs_i)
        ))
        vsel_tester(
          vs_regul,
          refmod_expected = refmods[[args_vs_i$tstsetup_ref]],
          solterms_len_expected = args_vs_i$nterms_max,
          method_expected = "forward",
          search_trms_empty_size =
            length(args_vs_i$search_terms) &&
            all(grepl("\\+", args_vs_i$search_terms)),
          info_str = tstsetup
        )
      }
      for (m in seq_len(m_max)) {
        # Selection:
        outdmin_jm_regul <- vs_regul$search_path$outdmins[[m]]
        if (ncl_crr == 1) {
          outdmin_jm_regul <- list(outdmin_jm_regul)
        } else {
          stopifnot(identical(ncl_crr, length(outdmin_jm_regul)))
        }
        for (nn in seq_len(ncl_crr)) {
          stopifnot(length(outdmin_jm_regul[[nn]]$alpha) == 1)
          ssq_regul_sel_alpha[j, m, nn] <- outdmin_jm_regul[[nn]]$alpha^2
          if (length(outdmin_jm_regul[[nn]]$beta) > 0) {
            ssq_regul_sel_beta[j, m, nn] <- sum(outdmin_jm_regul[[nn]]$beta^2)
          }
        }
        # Prediction:
        # Since varsel() doesn't output object `p_sub`, use the linear predictor
        # here (instead of the coefficients themselves, which would only be
        # accessible from `p_sub`):
        mu_jm_regul <- vs_regul$summaries$sub[[m]]$mu
        if (args_vs_i$prj_nm == "augdat") {
          mu_jm_regul <- augvec2augmat(mu_jm_regul)
        }
        mu_jm_regul <- vs_regul$refmodel$family$linkfun(mu_jm_regul)
        if (grepl("\\.with_offs", tstsetup)) {
          mu_jm_regul <- mu_jm_regul - offs_tst
        }
        # In fact, `sum((mu - offset - intercept)^2)` would make more sense than
        # `var(mu - offset) = sum((mu - offset - mean(mu - offset))^2)` but
        # since varsel() doesn't output object `p_sub`, the intercept from the
        # prediction is not accessible here.
        if (args_vs_i$prj_nm == "augdat") {
          # Take the maximum variance across the response categories (i.e., the
          # worst-case scenario):
          mu_jm_regul <- augmat2arr(mu_jm_regul)
          mu_jm_regul <- matrix(mu_jm_regul,
                                nrow = dim(mu_jm_regul)[1],
                                ncol = dim(mu_jm_regul)[2])
          var_jm_regul <- max(apply(mu_jm_regul, 2, var))
        } else {
          var_jm_regul <- var(mu_jm_regul)
        }
        ssq_regul_prd[j, m] <- var_jm_regul
      }
    }
    # Selection:
    # For the intercept-only model:
    for (nn in seq_len(dim(ssq_regul_sel_alpha)[3])) {
      expect_length(unique(ssq_regul_sel_alpha[, 1, !!nn]), 1)
    }
    expect_true(all(is.na(ssq_regul_sel_beta[, 1, ])), info = tstsetup)
    # All other (i.e., not intercept-only) models (note: as discussed at issue
    # #169, the intercept is not tested here to stay the same):
    ssq_regul_sel_beta_cond <- array(
      dim = dim(ssq_regul_sel_beta) + c(-1L, -1L, 0L)
    )
    for (j in seq_len(dim(ssq_regul_sel_beta)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_sel_beta)[2])[-1]) {
        for (nn in seq_len(dim(ssq_regul_sel_beta)[3])) {
          ssq_regul_sel_beta_cond[j - 1, m - 1, nn] <-
            ssq_regul_sel_beta[j, m, nn] < ssq_regul_sel_beta[j - 1, m, nn]
        }
      }
    }
    sum_as_unexpected <- 0L
    expect_true(sum(!ssq_regul_sel_beta_cond) <= sum_as_unexpected,
                info = tstsetup)
    # Prediction:
    # For the intercept-only model, the linear predictor consists only
    # of the intercept, so we expect no variation in `mu_jm_regul`:
    expect_true(all(ssq_regul_prd[, 1] <= 1e-12), info = tstsetup)
    # All other (i.e., not intercept-only) models:
    for (j in seq_len(dim(ssq_regul_prd)[1])[-1]) {
      for (m in seq_len(dim(ssq_regul_prd)[2])[-1]) {
        expect_lt(ssq_regul_prd[!!j, !!m], ssq_regul_prd[j - 1, m])
      }
    }
  }
})

## Penalty ----------------------------------------------------------------

test_that("`penalty` of invalid length fails", {
  skip_if_not(run_vs)
  tstsetups <- setdiff(
    setdiff(grep("\\.glm\\.", names(args_vs), value = TRUE),
            grep("\\.glm\\..*\\.forward", names(args_vs), value = TRUE)),
    grep("\\.glm\\..*\\.augdat\\.", names(args_vs), value = TRUE)
  )
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    formul_crr <- get_formul_from_fit(fits[[args_vs_i$tstsetup_fit]])
    formul_crr <- rm_addresp(formul_crr)
    penal_possbl <- get_penal_possbl(formul_crr)
    len_penal <- length(penal_possbl)
    # The `penalty` objects to be tested:
    penal_tst <- list(rep(1, len_penal + 1), rep(1, len_penal - 1))
    for (penal_crr in penal_tst) {
      expect_error(
        do.call(varsel, c(
          list(object = refmods[[args_vs_i$tstsetup_ref]],
               penalty = penal_crr),
          excl_nonargs(args_vs_i)
        )),
        paste0("^Incorrect length of penalty vector \\(should be ",
               len_penal, "\\)\\.$"),
        info = paste(tstsetup, which(sapply(penal_tst, identical, penal_crr)),
                     sep = "__")
      )
    }
  }
})

test_that("for forward search, `penalty` has no effect", {
  skip_if_not(run_vs)
  penal_tst <- 2
  tstsetups <- union(
    union(grep("\\.forward", names(vss), value = TRUE),
          grep("\\.glm\\.", names(vss), value = TRUE, invert = TRUE)),
    grep("\\.augdat\\.", names(vss), value = TRUE)
  )
  # To save time:
  if (!run_more) {
    tstsetups <- head(tstsetups, 1)
  }
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    if (args_vs_i$prj_nm == "augdat" && args_vs_i$fam_nm == "cumul") {
      warn_expected <- "non-integer #successes in a binomial glm!"
    } else if (!is.null(args_vs_i$avoid.increase)) {
      warn_expected <- warn_mclogit
    } else {
      warn_expected <- NA
    }
    expect_warning(
      vs_penal <- do.call(varsel, c(
        list(object = refmods[[args_vs_i$tstsetup_ref]],
             penalty = penal_tst),
        excl_nonargs(args_vs_i)
      )),
      warn_expected
    )
    expect_equal(vs_penal, vss[[tstsetup]], info = tstsetup)
  }
})

test_that("for L1 search, `penalty` has an expected effect", {
  skip_if_not(run_vs)
  tstsetups <- setdiff(
    setdiff(grep("\\.glm\\.", names(vss), value = TRUE),
            grep("\\.glm\\..*\\.forward", names(vss), value = TRUE)),
    grep("\\.glm\\..*\\.augdat\\.", names(vss), value = TRUE)
  )
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]

    formul_crr <- get_formul_from_fit(fits[[args_vs_i$tstsetup_fit]])
    formul_crr <- rm_addresp(formul_crr)
    penal_possbl <- get_penal_possbl(formul_crr)
    len_penal <- length(penal_possbl)
    penal_crr <- rep(1, len_penal)
    stopifnot(len_penal >= 3)
    # TODO: This test should be extended to also test the case where a
    # categorical predictor (more precisely, one of its dummy variables) or a
    # poly() term (more precisely, one of its lower-order terms resulting from
    # the expansion of the poly() term) gets zero or infinite penalty. For now,
    # the following code ensures that no categorical predictors and no poly()
    # terms get zero or infinite penalty.
    idx_cat <- grep("xca\\.", penal_possbl)
    idx_poly <- grep("poly[m]*\\(", penal_possbl)
    # Two predictors without cost:
    idx_penal_0 <- head(setdiff(seq_along(penal_crr),
                                c(idx_cat, idx_poly)),
                        2)
    stopifnot(length(idx_penal_0) == 2)
    # One predictor with infinite penalty:
    idx_penal_Inf <- head(setdiff(seq_along(penal_crr),
                                  c(idx_penal_0, idx_cat, idx_poly)),
                          1)
    stopifnot(length(idx_penal_Inf) == 1)
    penal_crr[idx_penal_0] <- 0
    penal_crr[idx_penal_Inf] <- Inf

    vs_penal <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$tstsetup_ref]],
           penalty = penal_crr),
      excl_nonargs(args_vs_i, nms_excl_add = "nterms_max")
    ))
    nterms_max_crr <- count_terms_in_formula(formul_crr) - 1L
    vsel_tester(
      vs_penal,
      refmod_expected = refmods[[args_vs_i$tstsetup_ref]],
      solterms_len_expected = nterms_max_crr,
      method_expected = "L1",
      info_str = tstsetup
    )
    # Check that the variables with no cost are selected first and the ones
    # with infinite penalty last:
    solterms_penal <- vs_penal$solution_terms
    solterms_penal <- sub("(I\\(.*as\\.logical\\(.*\\)\\))", "\\1TRUE",
                          solterms_penal)
    expect_identical(solterms_penal[seq_along(idx_penal_0)],
                     penal_possbl[idx_penal_0],
                     info = tstsetup)
    expect_identical(rev(solterms_penal)[seq_along(idx_penal_Inf)],
                     rev(penal_possbl[idx_penal_Inf]),
                     info = tstsetup)
  }
})

## L1 search and interactions ---------------------------------------------

test_that("L1 search handles three-way (second-order) interactions correctly", {
  skip_if_not(run_vs)
  skip_if_not_installed("rstanarm")
  warn_L1_ia_orig <- options(projpred.warn_L1_interactions = TRUE)
  main_terms_in_ia <- c("xca.2", "xco.3", "xco.1")
  all_ias_split <- lapply(seq_along(main_terms_in_ia), combn,
                          x = main_terms_in_ia, simplify = FALSE)
  all_ias <- unlist(lapply(all_ias_split, function(ia_split) {
    lapply(ia_split, all_ia_perms, is_split = TRUE)
  }))
  trms_universe_split_bu <- trms_universe_split
  trms_universe_split <<- union(trms_universe_split, all_ias)
  tstsetup <- head(grep("^rstanarm\\.glm", names(fits), value = TRUE), 1)
  args_fit_i <- args_fit[[tstsetup]]
  stopifnot(!(args_fit_i$pkg_nm == "rstanarm" && args_fit_i$fam_nm == "cumul"))
  fit_fun_nm <- get_fit_fun_nm(args_fit_i)
  args_fit_i$formula <- update(args_fit_i$formula,
                               . ~ . + xca.2 * xco.3 * xco.1)
  fit <- suppressWarnings(do.call(
    get(fit_fun_nm, asNamespace(args_fit_i$pkg_nm)),
    excl_nonargs(args_fit_i)
  ))
  args_ref_i <- args_ref[[paste0(tstsetup, ".trad")]]
  refmod <- do.call(get_refmodel, c(
    list(object = fit),
    excl_nonargs(args_ref_i)
  ))
  args_vs_i <- args_vs[[paste0(tstsetup,
                               ".trad.default_meth.default_search_trms")]]
  args_vs_i$refit_prj <- FALSE
  args_vs_i$nterms_max <- NULL
  expect_warning(
    vs <- do.call(varsel, c(
      list(object = refmod),
      excl_nonargs(args_vs_i)
    )),
    "was selected before all.+lower-order interaction terms have been selected",
    info = tstsetup
  )
  vsel_tester(
    vs,
    refmod_expected = refmod,
    solterms_len_expected = count_terms_in_formula(refmod$formula) - 1L,
    method_expected = "L1",
    refit_prj_expected = FALSE,
    ### Testing for non-increasing element `ce` (for increasing model size)
    ### doesn't make sense if the ranking of predictors involved in interactions
    ### has been changed, so we choose a higher `extra_tol` than by default:
    extra_tol = 1.2,
    ###
    info_str = tstsetup
  )
  rk <- ranking(vs)[["fulldata"]]
  expect_true(
    all(sapply(grep(":", rk), function(ia_idx) {
      main_terms_in_ia <- strsplit(rk[ia_idx], ":")[[1]]
      all_ias_split <- lapply(seq_len(length(main_terms_in_ia) - 1L), combn,
                              x = main_terms_in_ia, simplify = FALSE)
      ias_lower <- unlist(lapply(all_ias_split, function(ia_split) {
        lapply(ia_split, all_ia_perms, is_split = TRUE)
      }))
      return(all(which(rk %in% ias_lower) < ia_idx))
    })),
    info = tstsetup
  )
  trms_universe_split <<- trms_universe_split_bu
  options(warn_L1_ia_orig)
})

## search_terms -----------------------------------------------------------

test_that(paste(
  "including all terms in `search_terms` gives the same results as the default",
  "`search_terms`"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.alltrms", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    tstsetup_default <- sub("\\.alltrms", "\\.default_search_trms", tstsetup)
    if (!tstsetup_default %in% names(vss)) next
    expect_identical(vss[[tstsetup]], vss[[tstsetup_default]], info = tstsetup)
  }
})

test_that(paste(
  "forcing the inclusion of a term in the candidate models via `search_terms`",
  "works as expected"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.fixed", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    # In principle, `search_trms_tst$fixed$search_terms[1]` could be used
    # instead of `"xco.1"`, but that would seem like the forced term always has
    # to come first in `search_terms` (which is not the case):
    expect_identical(vss[[tstsetup]]$solution_terms[1], "xco.1",
                     info = tstsetup)
  }
})

test_that(paste(
  "forcing the exclusion of a term in the candidate models via `search_terms`",
  "works as expected"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.excluded", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    expect_false("xco.1" %in% vss[[tstsetup]]$solution_terms, info = tstsetup)
  }
})

test_that(paste(
  "forcing the skipping of a model size via `search_terms` works as expected"
), {
  skip_if_not(run_vs)
  tstsetups <- grep("\\.empty_size", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    soltrms_out <- vss[[tstsetup]]$solution_terms
    expect_true(
      grepl("\\+", soltrms_out[1]) && !any(grepl("\\+", soltrms_out[-1])),
      info = tstsetup
    )
  }
})

# cv_varsel() -------------------------------------------------------------

context("cv_varsel()")

test_that(paste(
  "`object` of class \"refmodel\", `method`, `cv_method`, `nterms_max`,",
  "`nclusters`, and `nclusters_pred` work"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(cvvss)) {
    mod_crr <- args_cvvs[[tstsetup]]$mod_nm
    fam_crr <- args_cvvs[[tstsetup]]$fam_nm
    prj_crr <- args_cvvs[[tstsetup]]$prj_nm
    meth_exp_crr <- args_cvvs[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    vsel_tester(
      cvvss[[tstsetup]],
      with_cv = TRUE,
      refmod_expected = refmods[[args_cvvs[[tstsetup]]$tstsetup_ref]],
      solterms_len_expected = args_cvvs[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = args_cvvs[[tstsetup]]$cv_method,
      valsearch_expected = args_cvvs[[tstsetup]]$validate_search,
      search_trms_empty_size =
        length(args_cvvs[[tstsetup]]$search_terms) &&
        all(grepl("\\+", args_cvvs[[tstsetup]]$search_terms)),
      info_str = tstsetup
    )
  }
})

test_that("invalid `object` fails", {
  expect_error(cv_varsel(rnorm(5)),
               "^no applicable method for")
})

test_that("invalid `method` fails", {
  for (tstsetup in names(refmods)) {
    expect_error(cv_varsel(refmods[[tstsetup]], method = "k-fold"),
                 "^Unknown search method$",
                 info = tstsetup)
    if (args_ref[[tstsetup]]$mod_nm != "glm") {
      expect_error(cv_varsel(refmods[[tstsetup]], method = "L1"),
                   paste("^L1 search is only supported for reference models",
                         "without multilevel and without additive",
                         "\\(\"smoothing\"\\) terms\\.$"),
                   info = tstsetup)
    }
    if (args_ref[[tstsetup]]$mod_nm == "glm" &&
        args_ref[[tstsetup]]$prj_nm == "augdat") {
      expect_error(cv_varsel(refmods[[tstsetup]], method = "L1"),
                   paste("^Currently, the augmented-data projection may not be",
                         "combined with an L1 search\\.$"),
                   info = tstsetup)
    }
  }
})

test_that("invalid `cv_method` fails", {
  for (tstsetup in names(refmods)) {
    expect_error(
      suppressWarnings(cv_varsel(refmods[[tstsetup]], cv_method = "k-fold")),
      "^Unknown `cv_method`\\.$",
      info = tstsetup
    )
  }
})

test_that("`seed` works (and restores the RNG state afterwards)", {
  skip_if_not(run_cvvs)
  # To save time:
  tstsetups <- union(
    grep("\\.glm\\.gauss", names(cvvss), value = TRUE),
    # Important for testing get_refmodel.brmsfit()'s internal `kfold_seed` (and
    # also `refprd_seed` if we are lucky and get a fold which separates out at
    # least one group):
    grep("^brms\\.(glmm|gamm)\\..*\\.kfold", names(cvvss), value = TRUE)
  )
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    cvvs_orig <- cvvss[[tstsetup]]
    rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_repr1 <- .Random.seed
    cvvs_repr <- suppressWarnings(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]]),
      excl_nonargs(args_cvvs_i)
    )))
    .Random.seed_repr2 <- .Random.seed
    rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
    # Expected equality:
    expect_equal(cvvs_repr, cvvs_orig, info = tstsetup)
    expect_equal(.Random.seed_repr2, .Random.seed_repr1, info = tstsetup)
    # Expected inequality:
    expect_false(isTRUE(all.equal(rand_new, rand_orig)), info = tstsetup)
  }
})

## refit_prj --------------------------------------------------------------

test_that("`refit_prj` works", {
  skip_if_not(run_cvvs)
  if (run_more) {
    tstsetups <- names(cvvss)
  } else {
    tstsetups <- head(grep("\\.glm\\.", names(cvvss), value = TRUE), 1)
  }
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    args_cvvs_i$refit_prj <- FALSE
    cvvs_reuse <- suppressWarnings(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]]),
      excl_nonargs(args_cvvs_i)
    )))
    mod_crr <- args_cvvs_i$mod_nm
    fam_crr <- args_cvvs_i$fam_nm
    prj_crr <- args_cvvs_i$prj_nm
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    vsel_tester(
      cvvs_reuse,
      with_cv = TRUE,
      refmod_expected = refmods[[args_cvvs_i$tstsetup_ref]],
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      refit_prj_expected = FALSE,
      cv_method_expected = args_cvvs_i$cv_method,
      valsearch_expected = args_cvvs_i$validate_search,
      search_trms_empty_size =
        length(args_cvvs_i$search_terms) &&
        all(grepl("\\+", args_cvvs_i$search_terms)),
      info_str = tstsetup
    )
  }
})

## nloo -------------------------------------------------------------------

test_that("invalid `nloo` fails", {
  for (tstsetup in names(refmods)) {
    # Use suppressWarnings() because of occasional warnings concerning Pareto k
    # diagnostics:
    expect_error(suppressWarnings(cv_varsel(refmods[[tstsetup]], nloo = -1)),
                 "^nloo must be at least 1$",
                 info = tstsetup)
  }
})

test_that(paste(
  "setting `nloo` at least as large as the number of observations doesn't",
  "change results"
), {
  skip_if_not(run_cvvs)
  nloo_tst <- nobsv + 1L
  tstsetups <- grep("\\.glm\\.gauss\\..*\\.default_cvmeth", names(cvvss),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    # Use suppressWarnings() because of occasional warnings concerning Pareto k
    # diagnostics:
    cvvs_nloo <- suppressWarnings(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]],
           nloo = nloo_tst),
      excl_nonargs(args_cvvs_i)
    )))
    expect_equal(cvvs_nloo, cvvss[[tstsetup]], info = tstsetup)
  }
})

test_that("setting `nloo` smaller than the number of observations works", {
  skip_if_not(run_cvvs)
  nloo_tst <- nobsv %/% 5L
  tstsetups <- grep("\\.glm\\.gauss\\..*\\.default_cvmeth", names(cvvss),
                    value = TRUE)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    tstsetup_ref <- args_cvvs_i$tstsetup_ref
    mod_crr <- args_cvvs_i$mod_nm
    fam_crr <- args_cvvs_i$fam_nm
    prj_crr <- args_cvvs_i$prj_nm
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    # Use suppressWarnings() because of occasional warnings concerning Pareto k
    # diagnostics and also because of the warning concerning subsampled LOO CV
    # (see issue #94):
    cvvs_nloo <- suppressWarnings(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]],
           nloo = nloo_tst),
      excl_nonargs(args_cvvs_i)
    )))
    vsel_tester(
      cvvs_nloo,
      with_cv = TRUE,
      refmod_expected = refmods[[tstsetup_ref]],
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "LOO",
      valsearch_expected = args_cvvs_i$validate_search,
      nloo_expected = nloo_tst,
      search_trms_empty_size =
        length(args_cvvs_i$search_terms) &&
        all(grepl("\\+", args_cvvs_i$search_terms)),
      info_str = tstsetup
    )
    # Expected equality for most components with a few exceptions:
    expect_equal(cvvs_nloo[setdiff(vsel_nms, vsel_nms_nloo)],
                 cvvss[[tstsetup]][setdiff(vsel_nms, vsel_nms_nloo)],
                 info = tstsetup)
    # Expected inequality for the exceptions (but note that the components from
    # `vsel_nms_nloo_opt` can be, but don't need to be differing):
    for (vsel_nm in setdiff(vsel_nms_nloo, vsel_nms_nloo_opt)) {
      expect_false(isTRUE(all.equal(cvvs_nloo[[vsel_nm]],
                                    cvvss[[tstsetup]][[vsel_nm]])),
                   info = paste(tstsetup, vsel_nm, sep = "__"))
    }
  }
})

## validate_search --------------------------------------------------------

test_that("`validate_search` works", {
  skip_if_not(run_cvvs)
  tstsetups <- grep("\\.default_cvmeth", names(cvvss), value = TRUE)
  if (!run_valsearch_always) {
    has_valsearch_true <- sapply(tstsetups, function(tstsetup_cvvs) {
      !isFALSE(args_cvvs[[tstsetup_cvvs]]$validate_search)
    })
    tstsetups <- tstsetups[has_valsearch_true]
  }
  suggsize_cond <- setNames(rep(NA, length(tstsetups)), nm = tstsetups)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    stopifnot(is.null(args_cvvs_i$validate_search) ||
                isTRUE(args_cvvs_i$validate_search))
    tstsetup_ref <- args_cvvs_i$tstsetup_ref
    mod_crr <- args_cvvs_i$mod_nm
    fam_crr <- args_cvvs_i$fam_nm
    prj_crr <- args_cvvs_i$prj_nm
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    # Use suppressWarnings() because of occasional warnings concerning Pareto k
    # diagnostics:
    cvvs_valsearch <- suppressWarnings(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]],
           validate_search = FALSE),
      excl_nonargs(args_cvvs_i)
    )))
    vsel_tester(
      cvvs_valsearch,
      with_cv = TRUE,
      refmod_expected = refmods[[tstsetup_ref]],
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "LOO",
      valsearch_expected = FALSE,
      search_trms_empty_size =
        length(args_cvvs_i$search_terms) &&
        all(grepl("\\+", args_cvvs_i$search_terms)),
      info_str = tstsetup
    )
    # Expected equality for most components with a few exceptions:
    expect_equal(cvvs_valsearch[setdiff(vsel_nms, vsel_nms_valsearch)],
                 cvvss[[tstsetup]][setdiff(vsel_nms, vsel_nms_valsearch)],
                 info = tstsetup)
    expect_identical(cvvs_valsearch$summaries$ref,
                     cvvss[[tstsetup]]$summaries$ref,
                     info = tstsetup)
    # Expected inequality for the exceptions (but note that the components from
    # `vsel_nms_valsearch_opt` can be, but don't need to be differing):
    for (vsel_nm in setdiff(vsel_nms_valsearch, vsel_nms_valsearch_opt)) {
      expect_false(isTRUE(all.equal(cvvs_valsearch[[vsel_nm]],
                                    cvvss[[tstsetup]][[vsel_nm]])),
                   info = paste(tstsetup, vsel_nm, sep = "__"))
    }
    # Check the expected inequalities more specifically:
    # Without a validated search, we expect increased LPPDs (and consequently
    # also an increased ELPD) in the submodels (since the hold-out fold was
    # included in the dataset for fitting the submodels):
    tol_crr <- 2e-1
    # Allow for just a small proportion of extreme differences:
    prop_as_expected <- 0.9
    for (j in seq_along(cvvs_valsearch$summaries$sub)) {
      expect_true(mean(cvvs_valsearch$summaries$sub[[j]]$lppd >=
                         cvvss[[tstsetup]]$summaries$sub[[j]]$lppd - tol_crr) >=
                    prop_as_expected,
                  info = paste(tstsetup, j, sep = "__"))
    }
    expect_true(all(summary(cvvs_valsearch)$selection$elpd.loo >=
                      summary(cvvss[[tstsetup]])$selection$elpd.loo),
                info = tstsetup)
    # Without a validated search, we expect overfitting in the suggested model
    # size:
    sgg_size_valsearch <- suggest_size(cvvs_valsearch, warnings = FALSE)
    sgg_size <- suggest_size(cvvss[[tstsetup]], warnings = FALSE)
    if (!is.na(sgg_size_valsearch) & !is.na(sgg_size)) {
      suggsize_cond[tstsetup] <- sgg_size_valsearch >= sgg_size
    }
  }
  sum_as_unexpected <- 2L
  expect_true(sum(!suggsize_cond, na.rm = TRUE) <= sum_as_unexpected)
})

## Arguments specific to K-fold CV ----------------------------------------

test_that("invalid `K` fails", {
  skip_if_not(length(fits) > 0)
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = 1),
               "^`K` must be at least 2\\.$")
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = 1000),
               "^`K` cannot exceed the number of observations\\.$")
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = c(4, 9)),
               "^`K` must be a single integer value\\.$")
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = "a"),
               "^`K` must be a single integer value\\.$")
  expect_error(cv_varsel(refmods[[1]], cv_method = "kfold", K = dat),
               "^`K` must be a single integer value\\.$")
})

test_that(paste(
  "`cvfits` (actually passed to init_refmodel()) works for rstanarm reference",
  "models"
), {
  skip_if_not(run_cvvs)
  tstsetups <- grep("^rstanarm\\..*\\.kfold", names(cvvss), value = TRUE)
  if (!run_cvfits_all) {
    tstsetups_tmp <- head(grep("\\.glmm\\.", tstsetups, value = TRUE), 1)
    if (length(tstsetups_tmp) == 0) {
      tstsetups_tmp <- head(tstsetups, 1)
    }
    tstsetups <- tstsetups_tmp
  }
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    tstsetup_fit <- args_cvvs_i$tstsetup_fit
    mod_crr <- args_cvvs_i$mod_nm
    fam_crr <- args_cvvs_i$fam_nm
    prj_crr <- args_cvvs_i$prj_nm
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    fit_crr <- fits[[tstsetup_fit]]
    K_crr <- args_cvvs_i$K

    # Refit `K_crr` times (note: below, the seed for constructing `folds_vec`
    # had to be changed in some cases to avoid unfavorable PRNG situations,
    # leading to technical issues such as nonconvergence of the submodel fitter;
    # this is also tied to the value of `seed_tst`):
    if (grepl("\\.glmm\\.", tstsetup)) {
      # Perform a grouped K-fold CV to test an edge case where all observations
      # belonging to the same level of a variable with group-level effects are
      # in the same fold, so prediction is performed for new levels (see, e.g.,
      # brms's GitHub issue #1286):
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rng_old <- get(".Random.seed", envir = .GlobalEnv)
      }
      # Make the construction of the CV folds reproducible:
      set.seed(seed2_tst * 3L)
      folds_vec <- loo::kfold_split_grouped(K = K_crr, x = dat$z.1)
      if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
    } else {
      folds_vec <- cv_folds(nobsv, K = K_crr, seed = seed2_tst)
    }
    # Additionally to suppressWarnings(), suppressMessages() could be used here
    # (but is not necessary since messages seem to be suppressed within
    # test_that()'s `code`); furthermore, try() is used because rstanarm
    # sometimes fails to refit:
    kfold_obj <- try(suppressWarnings(kfold(fit_crr,
                                            K = K_crr,
                                            folds = folds_vec,
                                            save_fits = TRUE,
                                            cores = 1)),
                     silent = TRUE)
    if (inherits(kfold_obj, "try-error")) {
      cat("Could not test `tstsetup = \"", tstsetup, "\"` in the rstanarm ",
          "`cvfits` test. Error message: \"",
          attr(kfold_obj, "condition")$message, "\"\n", sep = "")
      next
    }
    kfold_obj <- structure(list(fits = kfold_obj$fits[, "fit"]),
                           K = K_crr,
                           folds = folds_vec)

    # Create `"refmodel"` object with `cvfits`:
    refmod_crr <- do.call(get_refmodel, c(
      list(object = fit_crr, cvfits = kfold_obj),
      excl_nonargs(args_ref[[args_cvvs_i$tstsetup_ref]])
    ))

    # Run cv_varsel():
    cvvs_cvfits <- do.call(cv_varsel, c(
      list(object = refmod_crr),
      excl_nonargs(args_cvvs_i, nms_excl_add = "K")
    ))

    # Checks:
    vsel_tester(
      cvvs_cvfits,
      with_cv = TRUE,
      refmod_expected = refmod_crr,
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "kfold",
      valsearch_expected = args_cvvs_i$validate_search,
      search_trms_empty_size =
        length(args_cvvs_i$search_terms) &&
        all(grepl("\\+", args_cvvs_i$search_terms)),
      info_str = tstsetup
    )
    # Expected equality for some components:
    # TODO: Currently, `check.environment = FALSE` is needed. The reason is
    # probably that in the divergence minimizers, the projpred-extended family
    # is passed to argument `family` of the external model fitting functions
    # like lme4::glmer(). This should be fixed and then `check.environment =
    # FALSE` should be removed.
    expect_equal(cvvs_cvfits[setdiff(vsel_nms, vsel_nms_cvfits)],
                 cvvss[[tstsetup]][setdiff(vsel_nms, vsel_nms_cvfits)],
                 check.environment = FALSE,
                 info = tstsetup)
    # Expected inequality for the remaining components (but note that the
    # components from `vsel_nms_cvfits_opt` can be, but don't need to be
    # differing):
    for (vsel_nm in setdiff(vsel_nms_cvfits, vsel_nms_cvfits_opt)) {
      expect_false(isTRUE(all.equal(cvvs_cvfits[[vsel_nm]],
                                    cvvss[[tstsetup]][[vsel_nm]])),
                   info = paste(tstsetup, vsel_nm, sep = "__"))
    }
  }
})

test_that(paste(
  "`cvfits` (actually passed to init_refmodel()) works for brms reference",
  "models"
), {
  skip_if_not(run_cvvs)
  skip_if_not(packageVersion("brms") >= "2.16.4")
  tstsetups <- grep("^brms\\..*\\.kfold", names(cvvss), value = TRUE)
  if (!run_cvfits_all) {
    tstsetups_tmp <- head(grep("\\.glmm\\.", tstsetups, value = TRUE), 1)
    if (length(tstsetups_tmp) == 0) {
      tstsetups_tmp <- head(tstsetups, 1)
    }
    tstsetups <- tstsetups_tmp
  }
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    tstsetup_fit <- args_cvvs_i$tstsetup_fit
    mod_crr <- args_cvvs_i$mod_nm
    fam_crr <- args_cvvs_i$fam_nm
    prj_crr <- args_cvvs_i$prj_nm
    meth_exp_crr <- args_cvvs_i$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm" && prj_crr != "augdat",
                             "L1", "forward")
    }
    fit_crr <- fits[[tstsetup_fit]]
    K_crr <- args_cvvs_i$K

    # Refit `K_crr` times (note: below, the seed for constructing `folds_vec`
    # had to be changed in some cases to avoid unfavorable PRNG situations,
    # leading to technical issues such as nonconvergence of the submodel fitter;
    # this is also tied to the value of `seed_tst`):
    if (grepl("\\.glmm\\.", tstsetup)) {
      # Perform a grouped K-fold CV to test an edge case where all observations
      # belonging to the same level of a variable with group-level effects are
      # in the same fold, so prediction is performed for new levels (see, e.g.,
      # brms's GitHub issue #1286):
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rng_old <- get(".Random.seed", envir = .GlobalEnv)
      }
      # Make the construction of the CV folds reproducible:
      set.seed(seed2_tst + 10L)
      folds_vec <- loo::kfold_split_grouped(K = K_crr, x = dat$z.1)
      if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
    } else if (grepl("\\.gam\\.", tstsetup)) {
      folds_vec <- cv_folds(nobsv, K = K_crr, seed = seed2_tst + 10L)
    } else {
      folds_vec <- cv_folds(nobsv, K = K_crr, seed = seed2_tst)
    }
    kfold_obj <- kfold(fit_crr,
                       K = K_crr,
                       folds = folds_vec,
                       save_fits = TRUE,
                       seed = seed_fit)
    kfold_obj <- structure(list(fits = kfold_obj$fits[, "fit"]),
                           K = K_crr,
                           folds = folds_vec)

    # Create `"refmodel"` object with `cvfits`:
    refmod_crr <- do.call(get_refmodel, c(
      list(object = fit_crr, cvfits = kfold_obj),
      excl_nonargs(args_ref[[args_cvvs_i$tstsetup_ref]])
    ))

    # Run cv_varsel():
    cvvs_cvfits <- try(
      do.call(cv_varsel, c(
        list(object = refmod_crr),
        excl_nonargs(args_cvvs_i, nms_excl_add = "K")
      )),
      silent = TRUE
    )
    if (inherits(cvvs_cvfits, "try-error")) {
      cat("Failure for `tstsetup = \"", tstsetup, "\"` in the brms ",
          "`cvfits` test. Error message: \"",
          attr(cvvs_cvfits, "condition")$message, "\"\n", sep = "")
      # Check that this is a "pwrssUpdate" failure in lme4, so for solving this,
      # we would either need to tweak the lme4 tuning parameters manually (via
      # `...`) or change the data-generating mechanism here in the tests (to
      # obtain less extreme or more data):
      expect_true(grepl("pwrssUpdate", attr(cvvs_cvfits, "condition")$message),
                  info = tstsetup)
      # Furthermore, this should only occur in the `run_more = TRUE` case, so it
      # can be skipped (because there are enough other `tstsetups` for which
      # this works):
      expect_true(run_more, info = tstsetup)
      next
    }
    # Test the reproducibility of ref_predfun() when applied to new observations
    # (should be ensured by get_refmodel.brmsfit()'s internal `refprd_seed`):
    runif(1)
    cvvs_cvfits_repr <- do.call(cv_varsel, c(
      list(object = refmod_crr),
      excl_nonargs(args_cvvs_i, nms_excl_add = "K")
    ))

    # Checks:
    expect_equal(cvvs_cvfits, cvvs_cvfits_repr, info = tstsetup)
    vsel_tester(
      cvvs_cvfits,
      with_cv = TRUE,
      refmod_expected = refmod_crr,
      solterms_len_expected = args_cvvs_i$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "kfold",
      valsearch_expected = args_cvvs_i$validate_search,
      search_trms_empty_size =
        length(args_cvvs_i$search_terms) &&
        all(grepl("\\+", args_cvvs_i$search_terms)),
      info_str = tstsetup
    )
    # Expected equality for some components:
    # TODO: Currently, `check.environment = FALSE` is needed. The reason is
    # probably that in the divergence minimizers, the projpred-extended family
    # is passed to argument `family` of the external model fitting functions
    # like lme4::glmer(). This should be fixed and then `check.environment =
    # FALSE` should be removed.
    expect_equal(cvvs_cvfits[setdiff(vsel_nms, vsel_nms_cvfits)],
                 cvvss[[tstsetup]][setdiff(vsel_nms, vsel_nms_cvfits)],
                 check.environment = FALSE,
                 info = tstsetup)
    # Expected inequality for the remaining components (but note that the
    # components from `vsel_nms_cvfits_opt` can be, but don't need to be
    # differing):
    for (vsel_nm in setdiff(vsel_nms_cvfits, vsel_nms_cvfits_opt)) {
      expect_false(isTRUE(all.equal(cvvs_cvfits[[vsel_nm]],
                                    cvvss[[tstsetup]][[vsel_nm]])),
                   info = paste(tstsetup, vsel_nm, sep = "__"))
    }
  }
})
