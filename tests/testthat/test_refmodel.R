# get_refmodel() ----------------------------------------------------------

test_that("`object` of class `stanreg` or `brmsfit` works", {
  skip_if_no_setups(refmods)
  for (tstsetup in names(refmods)) {
    tstsetup_fit <- args_ref[[tstsetup]]$tstsetup_fit
    with_spclformul_crr <- grepl("\\.spclformul", tstsetup)
    if (args_ref[[tstsetup]]$fam_nm == "binom" ||
        grepl("\\.with_wobs", tstsetup)) {
      wobs_expected_crr <- wobs_tst
    } else {
      wobs_expected_crr <- rep(1, nobsv)
    }
    if (grepl("\\.with_offs", tstsetup)) {
      offs_expected_crr <- offs_tst
    } else {
      offs_expected_crr <- rep(0, nobsv)
    }
    if (args_ref[[tstsetup]]$prj_nm == "latent") {
      fam_orig_expected <- f_gauss
    } else if (args_ref[[tstsetup]]$pkg_nm == "brms" &&
               args_ref[[tstsetup]]$fam_nm %in% fam_nms_aug) {
      fam_orig_expected <- eval(args_fit[[tstsetup_fit]]$family)
    } else {
      if (args_ref[[tstsetup]]$pkg_nm == "rstanarm" &&
          args_ref[[tstsetup]]$fam_nm == "cumul") {
        f_cumul <- get_f_cumul()
      }
      fam_orig_expected <- get(paste0("f_", args_fit[[tstsetup_fit]]$fam_nm))
    }
    refmodel_tester(
      refmods[[tstsetup]],
      pkg_nm = args_ref[[tstsetup]]$pkg_nm,
      fit_expected = fits[[tstsetup_fit]],
      with_spclformul = with_spclformul_crr,
      wobs_expected = wobs_expected_crr,
      offs_expected = offs_expected_crr,
      fam_orig = fam_orig_expected,
      mod_nm = args_ref[[tstsetup]]$mod_nm,
      fam_nm = args_ref[[tstsetup]]$fam_nm,
      augdat_expected = args_ref[[tstsetup]]$prj_nm == "augdat",
      latent_expected = args_ref[[tstsetup]]$prj_nm == "latent",
      info_str = tstsetup
    )
  }
})

test_that("missing `data` fails", {
  skip_if_not_installed("rstanarm")
  fit_nodata <- suppressWarnings(rstanarm::stan_glm(
    dat$y_glm_gauss ~ dat$xco.1 + dat$xco.2 + dat$xco.3 +
      dat$xca.1 + dat$xca.2 + offset(dat$offs_col),
    family = f_gauss,
    weights = dat$wobs_col,
    chains = chains_tst, seed = seed_fit, iter = iter_tst, QR = TRUE,
    refresh = 0
  ))
  expect_error(
    get_refmodel(fit_nodata),
    paste("^is\\.data\\.frame\\(data\\) is not TRUE$")
  )
})

test_that("`formula` as a character string fails", {
  skip_if_not_installed("rstanarm")
  # If `formula` is a character string, rstanarm::stan_glm() is not able to find
  # objects supplied to arguments `weights` or `offset`, at least when using
  # devtools::test():
  fit_str <- suppressWarnings(rstanarm::stan_glm(
    "y_glm_gauss ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2",
    family = f_gauss, data = dat,
    chains = chains_tst, seed = seed_fit, iter = iter_tst, QR = TRUE,
    refresh = 0
  ))
  expect_error(get_refmodel(fit_str),
               "^inherits\\(formula, \"formula\"\\) is not TRUE$")
})

test_that("reference models lacking an intercept work", {
  args_fit_i <- args_fit$rstanarm.glm.gauss.stdformul.with_wobs.with_offs
  skip_if_not(!is.null(args_fit_i))
  fit_fun_nm <- get_fit_fun_nm(args_fit_i)
  fit_no_icpt <- suppressWarnings(do.call(
    get(fit_fun_nm, asNamespace(args_fit_i$pkg_nm)),
    c(list(formula = update(args_fit_i$formula, . ~ . - 1)),
      excl_nonargs(args_fit_i, nms_excl_add = "formula"))
  ))
  expect_message(
    refmod_no_icpt <- get_refmodel(fit_no_icpt),
    "Adding an intercept to `formula`",
    info = "rstanarm.glm.gauss.stdformul.with_wobs.with_offs"
  )
  nms_compare <- c("formula", "div_minimizer", "y", "wobs", "wdraws_ref",
                   "offset", "y_oscale")
  expect_equal(
    refmod_no_icpt[nms_compare],
    refmods$rstanarm.glm.gauss.stdformul.with_wobs.with_offs[nms_compare],
    tolerance = .Machine$double.eps,
    info = "rstanarm.glm.gauss.stdformul.with_wobs.with_offs"
  )
})

test_that("offsets specified via argument `offset` work", {
  args_fit_i <- args_fit$rstanarm.glm.gauss.stdformul.with_wobs.with_offs
  skip_if_not(!is.null(args_fit_i))
  fit_fun_nm <- get_fit_fun_nm(args_fit_i)
  upd_no_offs <- paste(". ~", sub(" \\+ offset\\(offs_col\\)", "",
                                  as.character(args_fit_i$formula[3])))
  fit_offs_arg <- suppressWarnings(do.call(
    get(fit_fun_nm, asNamespace(args_fit_i$pkg_nm)),
    c(list(formula = update(args_fit_i$formula, upd_no_offs),
           offset = offs_tst),
      excl_nonargs(args_fit_i, nms_excl_add = "formula"))
  ))
  refmod_offs_arg <- get_refmodel(fit_offs_arg)
  expect_equal(
    as.matrix(fit_offs_arg),
    as.matrix(fits$rstanarm.glm.gauss.stdformul.with_wobs.with_offs),
    tolerance = 1e-12,
    info = "rstanarm.glm.gauss.stdformul.with_wobs.with_offs"
  )
  nms_compare <- c("div_minimizer", "eta", "mu", "mu_offs", "dis", "y", "wobs",
                   "wdraws_ref", "offset", "y_oscale")
  expect_equal(
    refmod_offs_arg[nms_compare],
    refmods$rstanarm.glm.gauss.stdformul.with_wobs.with_offs[nms_compare],
    tolerance = .Machine$double.eps,
    info = "rstanarm.glm.gauss.stdformul.with_wobs.with_offs"
  )
})

test_that(paste(
  "binomial family with 1-column response and weights which are not all ones",
  "errors"
), {
  skip_if_not_installed("rstanarm")
  dat_prop <- within(dat, {
    ybinprop_glm <- y_glm_binom / wobs_col
  })
  fit_binom_1col_wobs <- suppressWarnings(rstanarm::stan_glm(
    ybinprop_glm ~ xco.1 + xco.2 + xco.3 + xca.1 + xca.2 + offset(offs_col),
    family = f_binom, data = dat_prop,
    weights = wobs_tst,
    chains = chains_tst, seed = seed_fit, iter = iter_tst, QR = TRUE,
    refresh = 0
  ))
  if ("rstanarm.glm.binom.stdformul.without_wobs.with_offs" %in% names(fits)) {
    expect_equal(
      as.matrix(fit_binom_1col_wobs),
      as.matrix(fits$rstanarm.glm.binom.stdformul.without_wobs.with_offs)
    )
  }
  expect_error(get_refmodel(fit_binom_1col_wobs),
               "response values must be numbers of successes")
})

test_that("function calls in group terms fail", {
  tstsetup <- "brms.glmm.brnll.stdformul.without_wobs.without_offs"
  args_fit_i <- args_fit[[tstsetup]]
  skip_if_not(!is.null(args_fit_i))
  fit_gr <- fits[[tstsetup]]
  fit_gr$formula <- update(fit_gr$formula,
                           . ~ . - (xco.1 | z.1) + (xco.1 | gr(z.1)))
  expect_error(
    refmod_gr <- get_refmodel(fit_gr),
    paste("Function calls on the right-hand side of a group-term `|` character",
          "are not allowed\\."),
    info = tstsetup
  )
})

test_that("extra arguments in s() or t2() terms fail", {
  args_fit_i <- args_fit$rstanarm.gam.gauss.spclformul.with_wobs.without_offs
  skip_if_not(!is.null(args_fit_i))
  fit_fun_nm <- get_fit_fun_nm(args_fit_i)
  fit_s <- suppressWarnings(do.call(
    get(fit_fun_nm, asNamespace(args_fit_i$pkg_nm)),
    c(list(formula = update(args_fit_i$formula,
                            . ~ . - s(s.1) + s(s.1, bs = "cr"))),
      excl_nonargs(args_fit_i, nms_excl_add = "formula"))
  ))
  expect_error(
    refmod_s <- get_refmodel(fit_s),
    "arguments other than predictors are not allowed",
    info = paste0("rstanarm.gam.gauss.stdformul.with_wobs.without_offs", "__s")
  )
  fit_t2 <- suppressWarnings(do.call(
    get(fit_fun_nm, asNamespace(args_fit_i$pkg_nm)),
    c(list(formula = update(args_fit_i$formula,
                            . ~ . - s(s.1) + t2(s.1, bs = "tp"))),
      excl_nonargs(args_fit_i, nms_excl_add = "formula"))
  ))
  expect_error(
    refmod_t2 <- get_refmodel(fit_t2),
    "arguments other than predictors are not allowed",
    info = paste0("rstanarm.gam.gauss.stdformul.with_wobs.without_offs", "__t2")
  )
})

test_that(paste(
  "the formula of an rstanarm GAMM reference model has no linear terms for",
  "the variables of the smooth terms, as for an rstanarm GAM"
), {
  skip_if_no_setups(refmods)
  skip_if_not(run_additive)
  tstsetups <- grep("^rstanarm\\.gam{1,2}\\.", names(refmods), value = TRUE)
  expect_true(any(grepl("^rstanarm\\.gamm\\.", tstsetups)))
  for (tstsetup in tstsetups) {
    fml <- refmods[[tstsetup]]$formula
    expect_identical(smooth_covered_terms(fml), character(0), info = tstsetup)
    expect_true(formula_contains_additive_terms(fml), info = tstsetup)
    if (grepl("^rstanarm\\.gamm\\.", tstsetup)) {
      expect_true(formula_contains_group_terms(fml), info = tstsetup)
    }
  }
})

test_that("get_refmodel() is idempotent", {
  skip_if_no_setups(refmods)
  for (tstsetup in names(refmods)) {
    expect_identical(get_refmodel(refmods[[tstsetup]]),
                     refmods[[tstsetup]],
                     info = tstsetup)
  }
})

test_that(paste(
  "the `cvrefbuilder` of get_refmodel.stanreg() takes the default `dis` from",
  "the fold-wise fit and passes a user-supplied `dis` unchanged"
), {
  tstsetups <- grep("^rstanarm\\..*\\.gauss\\.", names(fits), value = TRUE)
  skip_if_not(length(tstsetups) >= 2)
  # The second fit stands in for a fold refit:
  fit_full <- fits[[tstsetups[1]]]
  fit_fold <- fits[[tstsetups[2]]]
  sigma_full <- data.frame(fit_full)[, "sigma"]
  sigma_fold <- data.frame(fit_fold)[, "sigma"]
  skip_if(isTRUE(all.equal(sigma_full, sigma_fold)))

  refmod <- get_refmodel(fit_full)
  expect_identical(refmod$dis, sigma_full)
  expect_identical(refmod$cvrefbuilder(fit_fold)$dis, sigma_fold)

  dis_user <- rep(0.7, length(sigma_full))
  refmod_user <- get_refmodel(fit_full, dis = dis_user)
  expect_identical(refmod_user$cvrefbuilder(fit_fold)$dis, dis_user)
})

test_that(paste(
  "get_refmodel.stanreg() uses user-supplied `cvfun` and `cvrefbuilder`",
  "functions in K-fold CV"
), {
  tstsetups <- grep("^rstanarm\\..*\\.gauss\\.", names(fits), value = TRUE)
  skip_if_not(length(tstsetups) >= 1)
  fit <- fits[[tstsetups[1]]]

  # The user-supplied `cvfun` returns the full-data fit for each fold, so no
  # refit is needed:
  cvfun_folds <- list()
  cvfun_user <- function(folds) {
    cvfun_folds[[length(cvfun_folds) + 1L]] <<- folds
    rep(list(fit), max(folds))
  }
  cvrefbuilder_ks <- integer()
  cvrefbuilder_user <- function(cvfit) {
    cvrefbuilder_ks <<- c(cvrefbuilder_ks, cvfit$projpred_k)
    get_refmodel(cvfit, called_from_cvrefbuilder = TRUE)
  }
  refmod <- get_refmodel(fit, cvfun = cvfun_user,
                         cvrefbuilder = cvrefbuilder_user)
  expect_identical(refmod$cvfun, cvfun_user)
  expect_identical(refmod$cvrefbuilder, cvrefbuilder_user)

  folds <- rep_len(seq_len(K_tst), refmod$nobs)
  cvfits <- run_cvfun(refmod, folds = folds)
  expect_length(cvfun_folds, 1)
  expect_identical(cvfun_folds[[1]], folds)
  expect_identical(attr(cvfits, "folds"), folds)
  expect_length(cvfits, K_tst)

  cvvs <- cv_varsel(refmod, cv_method = "kfold", K = K_tst, cvfits = cvfits,
                    nclusters = nclusters_tst,
                    nclusters_pred = nclusters_pred_tst, nterms_max = 1,
                    verbose = FALSE, seed = seed_tst)
  expect_s3_class(cvvs, "vsel")
  expect_identical(cvrefbuilder_ks, seq_len(K_tst))
})

test_that(paste(
  "get_refmodel.stanreg() with only one of `cvfun` and `cvrefbuilder` keeps",
  "the default of the other, and the default `cvrefbuilder` does not pass a",
  "user-supplied `cvfun` to the folds"
), {
  tstsetups <- grep("^rstanarm\\..*\\.gauss\\.", names(fits), value = TRUE)
  skip_if_not(length(tstsetups) >= 2)
  # The second fit stands in for a fold refit:
  fit_full <- fits[[tstsetups[1]]]
  fit_fold <- fits[[tstsetups[2]]]
  sigma_fold <- data.frame(fit_fold)[, "sigma"]
  refmod_dflt <- get_refmodel(fit_full)

  cvfun_user <- function(folds) stop("Not used.")
  refmod <- get_refmodel(fit_full, cvfun = cvfun_user)
  expect_identical(refmod$cvfun, cvfun_user)
  refmod_fold <- refmod$cvrefbuilder(fit_fold)
  expect_false(identical(refmod_fold$cvfun, cvfun_user))
  expect_identical(refmod_fold$dis, sigma_fold)

  cvrefbuilder_user <- function(cvfit) stop("Not used.")
  refmod <- get_refmodel(fit_full, cvrefbuilder = cvrefbuilder_user)
  expect_identical(refmod$cvrefbuilder, cvrefbuilder_user)
  # The default `cvfun` is kept:
  expect_identical(body(refmod$cvfun), body(refmod_dflt$cvfun))
})

test_that(paste(
  "get_refmodel.stanreg() rejects a `cvfun` or `cvrefbuilder` that is not a",
  "function"
), {
  tstsetups <- grep("^rstanarm\\.", names(fits), value = TRUE)
  skip_if_not(length(tstsetups) >= 1)
  fit <- fits[[tstsetups[1]]]
  expect_error(get_refmodel(fit, cvfun = "a"),
               "^Argument `cvfun` must be a function or `NULL`\\.$")
  expect_error(get_refmodel(fit, cvrefbuilder = 1),
               "^Argument `cvrefbuilder` must be a function or `NULL`\\.$")
})

# predict.refmodel() ------------------------------------------------------

test_that("invalid `type` fails", {
  skip_if_not(length(fits) > 0)
  expect_error(predict(refmods[[1]], dat, type = "zzz"),
               "^type should be one of")
})

test_that("invalid `ynew` fails", {
  skip_if_not(length(fits) > 0)
  expect_error(predict(refmods[[1]], dat, ynew = dat),
               "^Argument `ynew` must be a numeric vector\\.$")
})

test_that(paste(
  "`object` of class `refmodel` and arguments `newdata`, `ynew`, and `type`",
  "work"
), {
  skip_if_no_setups(refmods)
  for (tstsetup in names(refmods)) {
    pkg_crr <- args_ref[[tstsetup]]$pkg_nm
    mod_crr <- args_ref[[tstsetup]]$mod_nm
    fam_crr <- args_ref[[tstsetup]]$fam_nm
    prj_crr <- args_ref[[tstsetup]]$prj_nm

    if (grepl("\\.with_wobs|\\.binom", tstsetup)) {
      wobs_crr <- wobs_tst
    } else {
      wobs_crr <- NULL
    }
    if (grepl("\\.with_offs", tstsetup)) {
      offs_crr <- offs_tst
    } else {
      offs_crr <- NULL
    }

    y_crr <- dat[, paste("y", mod_crr, fam_crr, sep = "_")]
    if (prj_crr == "latent") {
      dat_crr <- dat
      if (pkg_crr == "rstanarm" && grepl("\\.with_offs\\.", tstsetup)) {
        dat_crr$projpred_internal_offs_stanreg <- 0
      }
      y_nm <- stdize_lhs(refmods[[tstsetup]]$formula)$y_nm
      y_crr_link <- rowMeans(refmods[[tstsetup]]$ref_predfun(
        fit = refmods[[tstsetup]]$fit, newdata = dat_crr, excl_offs = FALSE,
        mlvl_allrandom = getOption("projpred.mlvl_proj_ref_new", FALSE)
      ))
    } else {
      y_crr_link <- y_crr
    }

    # Without `ynew`:
    expect_warning(
      withCallingHandlers(
        predref_resp <- predict(refmods[[tstsetup]], dat, weightsnew = wobs_crr,
                                offsetnew = offs_crr, type = "response"),
        warning = function(w) {
          if (grepl("(nobars|subbars|findbars|mkReTrms).*reformulas", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      get_warn_wrhs_orhs(tstsetup, weightsnew = wobs_crr,
                         offsetnew = offs_crr),
      info = tstsetup
    )

    expect_warning(
      predref_link <- predict(refmods[[tstsetup]], dat, weightsnew = wobs_crr,
                              offsetnew = offs_crr, type = "link"),
      get_warn_wrhs_orhs(tstsetup, weightsnew = wobs_crr,
                         offsetnew = offs_crr),
      info = tstsetup
    )

    # With `ynew`:
    expect_warning(
      predref_ynew_resp <- predict(refmods[[tstsetup]], dat,
                                   weightsnew = wobs_crr, offsetnew = offs_crr,
                                   ynew = y_crr, type = "response"),
      get_warn_wrhs_orhs(tstsetup, weightsnew = wobs_crr,
                         offsetnew = offs_crr),
      info = tstsetup
    )
    expect_warning(
      predref_ynew_link <- predict(refmods[[tstsetup]], dat,
                                   weightsnew = wobs_crr, offsetnew = offs_crr,
                                   ynew = y_crr_link, type = "link"),
      get_warn_wrhs_orhs(tstsetup, weightsnew = wobs_crr,
                         offsetnew = offs_crr),
      info = tstsetup
    )

    # Checks without `ynew`:
    if (prj_crr %in% c("latent", "augdat")) {
      if (prj_crr == "augdat" || !is.null(refmods[[tstsetup]]$family$cats)) {
        expect_identical(dim(predref_resp),
                         c(nobsv, length(refmods[[tstsetup]]$family$cats)),
                         info = tstsetup)
        expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                    info = tstsetup)
      } else {
        expect_true(is.vector(predref_resp, "double"), info = tstsetup)
        expect_length(predref_resp, nobsv)
        if (fam_crr %in% c("brnll", "binom")) {
          expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                      info = tstsetup)
        }
      }
      if (prj_crr == "augdat") {
        expect_identical(dim(predref_link),
                         c(nobsv, length(refmods[[tstsetup]]$family$cats) - 1L),
                         info = tstsetup)
      } else if (prj_crr == "latent") {
        expect_true(is.vector(predref_link, "double"), info = tstsetup)
        expect_length(predref_link, nobsv)
      }
    } else {
      expect_true(is.vector(predref_resp, "double"), info = tstsetup)
      expect_length(predref_resp, nobsv)
      if (fam_crr %in% c("brnll", "binom")) {
        expect_true(all(predref_resp >= 0 & predref_resp <= 1),
                    info = tstsetup)
      }
      expect_true(is.vector(predref_link, "double"), info = tstsetup)
      expect_length(predref_link, nobsv)
      if (fam_crr == "gauss") {
        expect_equal(predref_resp, predref_link, info = tstsetup)
      }
    }

    # Checks with `ynew`:
    if (prj_crr != "latent") {
      expect_equal(predref_ynew_resp, predref_ynew_link, info = tstsetup)
    } else {
      expect_false(isTRUE(all.equal(predref_ynew_resp, predref_ynew_link)),
                   info = tstsetup)
    }
    expect_true(is.vector(predref_ynew_resp, "double"), info = tstsetup)
    expect_length(predref_ynew_resp, nobsv)
    expect_false(isTRUE(all.equal(predref_ynew_resp, predref_resp)),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(predref_ynew_resp, predref_link)),
                 info = tstsetup)

    # Snapshots:
    if (run_snaps) {
      if (testthat_ed_max2) local_edition(3)
      withr::with_options(list(width = 145), {
        expect_snapshot({
          print(tstsetup)
          print(rlang::hash(predref_resp))
          print(rlang::hash(predref_link))
          print(rlang::hash(predref_ynew_resp))
          print(rlang::hash(predref_ynew_link))
        })
      })
      if (testthat_ed_max2) local_edition(2)
    }
  }
})

# 2-column binomial response ----------------------------------------------

test_that(paste(
  "init_refmodel() with a 2-column binomial response and the default",
  "`extract_model_data` uses the numbers of trials as observation weights"
), {
  withr::local_options(projpred.warn_cvrefbuilder_NULL = FALSE)
  set.seed(5302)
  n <- 40
  ndraws <- 20
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                    trials = sample(2:10, n, replace = TRUE))
  dat$y <- rbinom(n, dat$trials, plogis(0.5 * dat$x1))
  fit <- list(beta = rbind(rnorm(ndraws, 0, 0.1), rnorm(ndraws, 0.5, 0.1),
                           rnorm(ndraws, 0, 0.1)),
              data = dat)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(~ x1 + x2, nd) %*% fit$beta
  }
  refm_2col <- init_refmodel(
    fit, data = dat, formula = cbind(y, trials - y) ~ x1 + x2,
    family = binomial(), ref_predfun = ref_predfun
  )
  # The same model with a 1-column response and the trials as weights:
  emd_1col <- function(object, newdata, wrhs = NULL, orhs = NULL,
                       extract_y = TRUE) {
    y_wobs_offs(newdata = newdata, wrhs = wrhs %||% ~ trials, orhs = orhs,
                resp_form = if (extract_y) ~ y else NULL)
  }
  refm_1col <- init_refmodel(
    fit, data = dat, formula = y ~ x1 + x2, family = binomial(),
    ref_predfun = ref_predfun, extract_model_data = emd_1col
  )
  expect_equal(refm_2col$wobs, dat$trials)
  expect_equal(refm_2col$y, refm_1col$y)
  prj_2col <- suppressMessages(
    project(refm_2col, predictor_terms = "x1", ndraws = ndraws, seed = 5302)
  )
  prj_1col <- suppressMessages(
    project(refm_1col, predictor_terms = "x1", ndraws = ndraws, seed = 5302)
  )
  expect_equal(as.matrix(prj_2col), as.matrix(prj_1col))
  expect_equal(proj_linpred(prj_2col, newdata = dat)$lpd,
               proj_linpred(prj_1col, newdata = dat)$lpd)
  # Without the response columns in `newdata`, the weights need `weightsnew`:
  nd <- dat[1:5, c("x1", "x2")]
  expect_error(proj_linpred(prj_2col, newdata = nd),
               "or supply argument `weightsnew`", fixed = TRUE)
  # `pred` does not use the weights; this checks that the call with
  # `weightsnew` does not stop:
  expect_equal(
    proj_linpred(prj_2col, newdata = nd, weightsnew = dat$trials[1:5])$pred,
    proj_linpred(prj_1col, newdata = nd, weightsnew = dat$trials[1:5])$pred
  )
  expect_error(predict(refm_2col, newdata = nd),
               "or supply argument `weightsnew`", fixed = TRUE)
  expect_equal(
    predict(refm_2col, newdata = nd, weightsnew = dat$trials[1:5]),
    predict(refm_1col, newdata = nd, weightsnew = dat$trials[1:5])
  )
  expect_error(proj_predict(prj_2col, newdata = nd),
               "or supply argument `weightsnew`", fixed = TRUE)
  expect_equal(
    proj_predict(prj_2col, newdata = nd, weightsnew = dat$trials[1:5],
                 .seed = 5302),
    proj_predict(prj_1col, newdata = nd, weightsnew = dat$trials[1:5],
                 .seed = 5302)
  )

  # Two response columns that are both columns of `data`: the second column
  # must stay unchanged in the stored data.
  dat$fail <- dat$trials - dat$y
  fit$data <- dat
  refm_sf <- init_refmodel(
    fit, data = dat, formula = cbind(y, fail) ~ x1 + x2,
    family = binomial(), ref_predfun = ref_predfun
  )
  expect_equal(refm_sf$wobs, dat$trials)
  expect_equal(refm_sf$fetch_data()$fail, dat$fail)
  prj_sf <- suppressMessages(
    project(refm_sf, predictor_terms = "x1", ndraws = ndraws, seed = 5302)
  )
  expect_equal(proj_linpred(prj_sf)$lpd, proj_linpred(prj_1col)$lpd)

  # A constant number of trials from the formula environment: only the
  # response variables that are columns of `data` must be in `newdata`.
  nt <- 10L
  dat_nt <- dat
  dat_nt$y <- rbinom(n, nt, plogis(0.5 * dat_nt$x1))
  fit_nt <- fit
  fit_nt$data <- dat_nt
  refm_nt <- init_refmodel(
    fit_nt, data = dat_nt, formula = cbind(y, nt - y) ~ x1 + x2,
    family = binomial(), ref_predfun = ref_predfun
  )
  expect_equal(refm_nt$wobs, rep(nt, n))
  prj_nt <- suppressMessages(
    project(refm_nt, predictor_terms = "x1", ndraws = ndraws, seed = 5302)
  )
  expect_no_error(
    proj_linpred(prj_nt, newdata = dat_nt[1:5, c("x1", "x2", "y")])
  )
})
