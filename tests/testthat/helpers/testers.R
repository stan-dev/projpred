# A helper function for testing the structure of an expected extended `"family"`
# object.
#
# @param extfam An object of class `"family"` (at least expected so) which has
#   been extended by projpred.
# @param fam_orig The original object of class `"family"` which has been used as
#   input to extend_family().
# @param extfam_nms_add2 A character vector of additional element names which
#   do not exist in the output of extend_family() (needed for testing
#   `get_refmodel([...])$family`).
# @param from_brms A single logical value indicating whether the extended family
#   was created for a reference model based on a `"brmsfit"` (`TRUE`) or not
#   (`FALSE`).
# @param augdat_expected A single logical value indicating whether the extended
#   family is expected to be for augmented-data projection (`TRUE`) or not
#   (`FALSE`).
# @param latent_expected A single logical value indicating whether the reference
#   model is expected to be for latent projection (`TRUE`) or not (`FALSE`).
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
extfam_tester <- function(extfam,
                          fam_orig,
                          extfam_nms_add2 = character(),
                          from_brms = FALSE,
                          augdat_expected = FALSE,
                          latent_expected = FALSE,
                          info_str) {
  # General structure tests -------------------------------------------------

  ## For `fam_orig` ---------------------------------------------------------

  expect_s3_class(fam_orig, "family")
  expect_type(fam_orig, "list")
  # Basic names of `fam_orig` (will be adapted later if necessary):
  fam_orig_nms <- c("family", "link")
  # Check for element `family`:
  expect_true("family" %in% names(fam_orig), info = info_str)
  expect_true(is.vector(fam_orig$family, "character"), info = info_str)
  expect_length(fam_orig$family, 1)
  # Family "cumulative_rstanarm" is an artificial family, so it lacks some
  # elements present in actual families:
  if (fam_orig$family != "cumulative_rstanarm") {
    fam_orig_nms <- c(fam_orig_nms, "linkfun", "linkinv")
  }
  # Families from brms which require the augmented-data projection:
  bfam_nms <- c("categorical", "cumulative", "cratio", "sratio", "acat")
  # Check the names of `fam_orig`:
  if (fam_orig$family %in% bfam_nms) {
    # Just check that the basic names exist (to avoid depending too strongly on
    # brms internals):
    expect_true(all(fam_orig_nms %in% names(fam_orig)), info = info_str)
  } else {
    if (fam_orig$family != "cumulative_rstanarm") {
      fam_orig_nms <- c(fam_orig_nms, "variance", "dev.resids", "aic", "mu.eta",
                        "initialize", "validmu", "valideta")
      if (fam_orig$family %in% c("binomial", "poisson")) {
        fam_orig_nms <- c(fam_orig_nms, "simulate")
      }
    }
    expect_named(fam_orig, fam_orig_nms, info = info_str)
  }

  ## For `extfam` -----------------------------------------------------------

  expect_s3_class(extfam, "family")
  expect_type(extfam, "list")
  expect_true("for_augdat" %in% names(extfam), info = info_str)
  expect_true(isTRUE(extfam$for_augdat) || isFALSE(extfam$for_augdat),
              info = info_str)
  expect_identical(extfam$for_augdat, augdat_expected, info = info_str)
  expect_true("for_latent" %in% names(extfam), info = info_str)
  expect_true(isTRUE(extfam$for_latent) || isFALSE(extfam$for_latent),
              info = info_str)
  expect_identical(extfam$for_latent, latent_expected, info = info_str)
  extfam_nms_add <- c("kl", "dis_fun", "predvar", "ll_fun", "deviance", "ppd",
                      "for_latent", "for_augdat", "is_extended")
  if (extfam$for_augdat) {
    extfam_nms_add <- setdiff(extfam_nms_add, "deviance")
    extfam_nms_add <- c(extfam_nms_add, "cats", "kl_ptwise")
    if (extfam$family == "categorical") {
      extfam_nms_add <- c(extfam_nms_add, "refcat")
    }
    if (extfam$family == "cumulative_rstanarm") {
      extfam_nms_add <- c(extfam_nms_add, "linkfun", "linkinv")
    }
  } else if (extfam$for_latent) {
    extfam_nms_add <- c(extfam_nms_add, "familyOrig", "linkOrig",
                        "respOrig_possible", "ppdOrig_possible", "latent_ilink",
                        "latent_llOrig", "latent_ppdOrig")
    if (extfam$familyOrig != "binomial") {
      extfam_nms_add <- c(extfam_nms_add, "cats")
    }
  }
  extfam_nms_add <- c(extfam_nms_add, extfam_nms_add2)
  extfam_nms <- c(fam_orig_nms, extfam_nms_add)
  if (fam_orig$family %in% bfam_nms) {
    expect_true(all(extfam_nms %in% names(extfam)), info = info_str)
  } else {
    expect_named(extfam, extfam_nms,
                 ignore.order = extfam$for_augdat || extfam$for_latent,
                 info = info_str)
  }

  # Detailed tests ----------------------------------------------------------

  ## For `fam_orig` ---------------------------------------------------------

  fam_orig_ch <- structure(
    extfam[setdiff(names(fam_orig), c("dpars", "multi_dpars"))],
    class = if (fam_orig$family %in% bfam_nms) {
      c("brmsfamily", "family")
    } else {
      "family"
    }
  )
  if (extfam$family == "binomial") {
    fam_orig_ch$initialize <- fam_orig$initialize
    if (extfam$for_augdat) {
      expect_identical(
        get("augdat_link", envir = environment(fam_orig_ch$linkfun)),
        augdat_link_binom,
        info = info_str
      )
      expect_identical(
        get("augdat_ilink", envir = environment(fam_orig_ch$linkinv)),
        augdat_ilink_binom,
        info = info_str
      )
      fam_orig_ch$linkfun <- fam_orig$linkfun
      fam_orig_ch$linkinv <- fam_orig$linkinv
    }
  } else if (extfam$for_augdat) {
    fam_orig_ch$linkfun <- fam_orig$linkfun
    fam_orig_ch$linkinv <- fam_orig$linkinv
  }
  if (!from_brms) {
    expect_identical(fam_orig_ch, fam_orig,
                     ignore.environment = extfam$for_latent, info = info_str)
  } else if (extfam$family %in% bfam_nms) {
    expect_identical(
      fam_orig_ch,
      structure(fam_orig[setdiff(names(fam_orig), c("dpars", "multi_dpars"))],
                class = class(fam_orig)),
      ignore.environment = TRUE,
      info = info_str
    )
  } else {
    expect_identical(fam_orig_ch, fam_orig, ignore.environment = TRUE,
                     info = info_str)
  }

  ## For `extfam` -----------------------------------------------------------

  if ("familyOrig" %in% names(extfam)) {
    expect_true(extfam$familyOrig %in% fam_nms_long, info = info_str)
  }
  if ("linkOrig" %in% names(extfam)) {
    expect_true(extfam$linkOrig %in% link_str, info = info_str)
  }
  if ("respOrig_possible" %in% names(extfam)) {
    expect_true(extfam$respOrig_possible, info = info_str)
  }
  if ("ppdOrig_possible" %in% names(extfam)) {
    expect_true(extfam$ppdOrig_possible, info = info_str)
  }
  if ("latent_llOrig" %in% names(extfam)) {
    if (extfam$familyOrig == "binomial") {
      expect_identical(extfam$latent_llOrig, latent_llOrig_binom_nocats,
                       info = info_str)
    } else {
      expect_identical(extfam$latent_llOrig, latent_llOrig_cats,
                       info = info_str)
    }
  }
  if ("latent_ppdOrig" %in% names(extfam)) {
    if (extfam$familyOrig == "binomial") {
      expect_identical(extfam$latent_ppdOrig, latent_ppdOrig_binom_nocats,
                       info = info_str)
    } else {
      expect_identical(extfam$latent_ppdOrig, latent_ppdOrig_cats,
                       info = info_str)
    }
  }
  if ("refcat" %in% names(extfam)) {
    expect_true(is.vector(extfam$refcat, "character"), info = info_str)
    expect_length(extfam$refcat, 1)
  }
  if ("cats" %in% names(extfam)) {
    expect_true(is.vector(extfam$cats, "character"), info = info_str)
  }
  expect_true(extfam$is_extended, info = info_str)
  el_nms_clos <- setdiff(
    extfam_nms_add,
    c("familyOrig", "linkOrig", "respOrig_possible", "ppdOrig_possible",
      "refcat", "cats", "for_latent", "for_augdat", "is_extended")
  )
  for (el_nm in el_nms_clos) {
    expect_type(extfam[[el_nm]], "closure")
  }

  if (extfam$for_augdat) {
    arr_pr <- abind::abind(array(c(0.7, 0.6, 0.3, 0.4), dim = c(2, 1, 2)),
                           array(c(0.1, 0.2, 0.9, 0.8), dim = c(2, 1, 2)),
                           along = 2)
    augm_pr <- arr2augmat(arr_pr, margin_draws = 1)
    expect_equal(extfam$linkinv(extfam$linkfun(augm_pr)), augm_pr,
                 info = info_str)
    # We expect an N x S matrix:
    expect_equal(extfam$kl_ptwise(mu_ref = augm_pr, mu_sub = augm_pr),
                 matrix(0, nrow = 2, ncol = 2),
                 info = info_str)
    # We expect a vector of length S:
    expect_equal(extfam$kl(pref = list(mu = augm_pr),
                           data = list(weights = rep(1, 2)),
                           psub = list(mu = augm_pr)),
                 numeric(2),
                 info = info_str)
    # We expect a vector of length S:
    expect_equal(extfam$dis_fun(pref = list(mu = augm_pr),
                                psub = list(mu = augm_pr)),
                 rep(NA, 2),
                 info = info_str)
    # We expect a vector of length N_augcat = nrow(augm_pr):
    expect_equal(extfam$predvar(mu = augm_pr, dis = NA),
                 rep(NA, 4),
                 info = info_str)
    # We expect an N x S matrix:
    expect_equal(extfam$ll_fun(mu = augm_pr,
                               y = factor(head(letters, 2)[c(2, 1)])),
                 matrix(log(c(0.3, 0.1, 0.4, 0.2)), ncol = 2),
                 info = info_str)
    # We expect a vector of length N:
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(seed2_tst)
    ppd_draws <- extfam$ppd(mu = augm_pr[, 1])
    expect_true(is.vector(ppd_draws, "integer"), info = info_str)
    expect_length(ppd_draws, 2)
  }
  # TODO: For the traditional (and latent) projection, add some mathematical
  # checks (i.e., check that the calculations for the objects listed in
  # `extfam_nms_add` are mathematically correct).

  # Output ------------------------------------------------------------------

  return(invisible(TRUE))
}

# A helper function for testing the structure of an expected `"refmodel"`
# object.
#
# @param refmod An object of class `"refmodel"` (at least expected so).
# @param is_datafit A single logical value indicating whether the reference
#   model is expected to be a `"datafit"` (`TRUE`) or not (`FALSE`).
# @param pkg_nm A single character string specifying the name of the package
#   upon whose fit the reference model (or `"datafit"`) is based.
# @param fit_expected The expected `refmod$fit` object.
# @param formul_expected The expected `refmod$formula` object. A cbind()
#   expression on the left-hand side of the formula is handled automatically.
# @param data_expected The original dataset used for the reference model fit or
#   as input to get_refmodel() or init_refmodel(). Internal changes (i.e.,
#   inside of projpred) of this dataset are performed automatically.
# @param with_spclformul A single logical value indicating whether the reference
#   model has a special formula (`TRUE`) or not (`FALSE`).
# @param nobsv_expected A single integer value giving the expected number of
#   observations.
# @param wobs_expected The expected numeric vector of observation weights.
# @param offs_expected The expected numeric vector of offsets.
# @param nrefdraws_expected A single integer value giving the expected number of
#   posterior draws in the reference model.
# @param fam_orig The original object of class `"family"` which has been used as
#   input to extend_family().
# @param mod_nm A single character string specifying the type of model (see
#   object `mod_nms`.
# @param fam_nm A single character string specifying the family (see object
#   `fam_nms`.
# @param augdat_expected A single logical value indicating whether the reference
#   model is expected to be for augmented-data projection (`TRUE`) or not
#   (`FALSE`).
# @param latent_expected A single logical value indicating whether the reference
#   model is expected to be for latent projection (`TRUE`) or not (`FALSE`).
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
refmodel_tester <- function(
    refmod,
    is_datafit = FALSE,
    pkg_nm,
    fit_expected,
    formul_expected = get_formul_from_fit(fit_expected),
    data_expected = dat,
    with_spclformul = FALSE,
    nobsv_expected = nobsv,
    wobs_expected = wobs_tst,
    offs_expected = offs_tst,
    nrefdraws_expected = chains_tst * (iter_tst %/% 2L),
    fam_orig,
    mod_nm,
    fam_nm,
    augdat_expected = FALSE,
    latent_expected = FALSE,
    info_str
) {
  # Preparations:
  needs_wobs_added <- !is_datafit && pkg_nm == "rstanarm" &&
    length(refmod$fit$weights) > 0
  if (needs_wobs_added) {
    data_expected$projpred_internal_wobs_stanreg <- refmod$fit$weights
  }
  needs_offs_added <- !is_datafit && pkg_nm == "rstanarm" &&
    length(refmod$fit$offset) > 0
  if (needs_offs_added) {
    data_expected$projpred_internal_offs_stanreg <- refmod$fit$offset
  }
  if (refmod$family$for_latent) {
    formul_expected[[2]] <- str2lang(
      paste0(".", as.character(formul_expected[[2]]))
    )
  }
  if (!is.null(attr(terms(formul_expected), "offset"))) {
    # In the reference model, the offset() term is placed last:
    formul_expected <- update(formul_expected,
                              . ~ . - offset(offs_col) + offset(offs_col))
  }
  stdized_lhs <- stdize_lhs(formul_expected)
  formul_expected <- stdized_lhs$fml
  if (with_spclformul) {
    # Reference models take arithmetic expressions on the left-hand side of
    # the formula into account:
    y_spclformul <- stdized_lhs$y_nm_orig
    y_spclformul_new <- stdized_lhs$y_nm
    formul_expected <- update(
      formul_expected,
      as.formula(paste(y_spclformul_new, "~ ."))
    )
    data_expected <- within(data_expected, {
      assign(y_spclformul_new, eval(str2lang(y_spclformul)))
    })
  }
  if (!is_datafit && pkg_nm == "rstanarm" &&
      refmod$fit$stan_function == "stan_gamm4" &&
      refmod$family$familyOrig %||% refmod$family$family == "binomial") {
    # A column added internally by rstanarm which is not relevant for projpred:
    data_expected$temp_y <- 1
  }
  has_grp <- mod_nm %in% c("glmm", "gamm")
  has_add <- mod_nm %in% c("gam", "gamm")
  is_gamm <- mod_nm == "gamm"

  # Test the general structure of the object:
  refmod_nms <- c(
    "fit", "formula", "div_minimizer", "family", "mu", "eta", "dis", "y",
    "loglik", "intercept", "proj_predfun", "fetch_data", "wobs", "wsample",
    "offset", "cvfun", "cvfits", "extract_model_data", "ref_predfun",
    "cvrefbuilder", "yOrig"
  )
  refmod_class_expected <- "refmodel"
  if (is_datafit) {
    refmod_class_expected <- c("datafit", refmod_class_expected)
  }
  expect_s3_class(refmod, refmod_class_expected, exact = TRUE)
  expect_type(refmod, "list")
  expect_named(refmod, refmod_nms, info = info_str)

  # fit
  expect_identical(refmod$fit, fit_expected, info = info_str)

  # formula
  if (!is_gamm) {
    # TODO (GAMMs): Adapt the expected formula to GAMMs.
    if (is_datafit && pkg_nm == "brms") {
      expect_equal(refmod$formula, formul_expected, info = info_str)
    } else {
      expect_identical(refmod$formula, formul_expected, info = info_str)
    }
  }

  # div_minimizer
  expect_type(refmod$div_minimizer, "closure")

  # family
  extfam_tester(refmod$family,
                fam_orig = fam_orig,
                extfam_nms_add2 = "mu_fun",
                from_brms = (pkg_nm == "brms"),
                augdat_expected = augdat_expected,
                latent_expected = latent_expected,
                info_str = info_str)

  # mu
  ### Not needed because of the more precise test below:
  # expect_true(is.matrix(refmod$mu), info = info_str)
  # expect_type(refmod$mu, "double")
  # expect_identical(dim(refmod$mu), c(nobsv_expected, nrefdraws_expected),
  #                  info = info_str)
  ###
  if (!is_datafit) {
    ### Helpful for debugging:
    # mu_expected_ch <- unname(t(posterior_linpred(refmod$fit)))
    ###
    if (pkg_nm == "rstanarm") {
      # In this case, the linear predictors are calculated manually because of
      # the offset issues in rstanarm.
      drws <- as.matrix(refmod$fit)
      if ("(Intercept)" %in% colnames(drws)) {
        drws_icpt <- drws[, "(Intercept)"]
      } else {
        drws_icpt <- numeric(nrow(drws))
        drws_thres <- drws[, grep("\\|", colnames(drws))]
      }
      drws_beta_cont <- drws[
        ,
        setdiff(grep("xco\\.", colnames(drws), value = TRUE),
                grep("z\\.", colnames(drws), value = TRUE)),
        drop = FALSE
      ]
      mm_cont <- model.matrix(
        as.formula(paste("~",
                         paste(colnames(drws_beta_cont), collapse = " + "))),
        data = data_expected
      )
      stopifnot(identical(c("(Intercept)", colnames(drws_beta_cont)),
                          colnames(mm_cont)))
      mu_expected <- cbind(drws_icpt, drws_beta_cont) %*% t(mm_cont)
      cate_post <- lapply(names(x_cate_list), function(x_cate_idx) {
        sapply(x_cate_list[[x_cate_idx]]$x_cate, function(lvl_obs_i) {
          if (lvl_obs_i != "lvl1") {
            return(drws[, paste0("xca.", x_cate_idx, lvl_obs_i)])
          } else {
            return(matrix(0, nrow = nrow(drws)))
          }
        })
      })
      mu_expected <- mu_expected + do.call("+", cate_post)
      if (has_grp) {
        r_post <- lapply(names(z_list), function(z_nm) {
          unname(
            drws[, paste0("b[(Intercept) ", z_nm, ":", z_list[[z_nm]]$z, "]")] +
              drws[, paste0("b[xco.1 ", z_nm, ":", z_list[[z_nm]]$z, "]")] %*%
              diag(x_cont[, 1])
          )
        })
        mu_expected <- mu_expected + do.call("+", r_post)
      }
      if (has_add) {
        drws_beta_s <- drws[
          ,
          grep("^s\\(", colnames(drws), value = TRUE),
          drop = FALSE
        ]
        ### TODO (GAMs and GAMMs): Do this manually:
        mm_s <- rstanarm:::pp_data(refmod$fit)$x
        mm_s <- mm_s[, grep("^s\\(", colnames(mm_s), value = TRUE),
                     drop = FALSE]
        ###
        mu_expected <- mu_expected + drws_beta_s %*% t(mm_s)
      }
      mu_expected <- unname(mu_expected)
    } else if (pkg_nm == "brms") {
      mu_expected <- sweep(posterior_linpred(refmod$fit), 2L, offs_expected)
      if (fam_nm %in% fam_nms_ordin) {
        drws <- as.matrix(refmod$fit)
        drws_thres <- drws[, grep("b_Intercept\\[", colnames(drws))]
      }
    }
    if (refmod$family$family != "gaussian") {
      if (refmod$family$family %in% fam_nms_aug_long) {
        if (refmod$family$family %in% fam_nms_ordin_long) {
          if (refmod$family$family %in% c("cumulative", "cumulative_rstanarm",
                                          "sratio")) {
            mu_expected <- apply(drws_thres, 2, function(thres_vec) {
              thres_vec - mu_expected
            }, simplify = FALSE)
          } else if (refmod$family$family %in% c("cratio", "acat")) {
            mu_expected <- apply(drws_thres, 2, function(thres_vec) {
              mu_expected - thres_vec
            }, simplify = FALSE)
          }
          mu_expected <- do.call(abind::abind, c(mu_expected, rev.along = 0))
        }
        mu_expected <- arr2augmat(mu_expected, margin_draws = 1)
        mu_expected <- refmod$family$linkinv(mu_expected)
      } else {
        mu_expected <- fam_orig$linkinv(mu_expected)
      }
    }
    if (refmod$family$for_augdat && refmod$family$family == "binomial") {
      mu_expected <- cbind(1 - mu_expected,
                           mu_expected)
    }
    if (!refmod$family$family %in% fam_nms_aug_long) {
      mu_expected <- t(mu_expected)
    }
    if (refmod$family$for_augdat && refmod$family$family == "binomial") {
      mu_expected <- structure(mu_expected,
                               nobs_orig = nobsv,
                               class = "augmat")
    }
    expect_equal(refmod$mu, mu_expected, info = info_str)
  } else {
    if (refmod$family$family != "binomial") {
      expect_identical(refmod$mu, as.matrix(refmod$y), info = info_str)
    } else {
      expect_identical(refmod$mu, as.matrix(refmod$y / refmod$wobs),
                       info = info_str)
    }
  }

  # eta
  # In principle, it would be desirable to compare `refmod$eta` to
  # `refmod$family$linkfun(refmod$mu)`, but numerical underflow and overflow can
  # make this problematic. (Here in the unit tests, we generate rather extreme
  # linear predictors, which should be avoided in the first place, but doesn't
  # seem to be that simple.)
  if (refmod$family$family %in% c("binomial")) {
    eta_cut <- refmod$eta
    mu_cut <- refmod$mu
    tol_ex <- 1e-12
    eta_cut[eta_cut < f_binom$linkfun(tol_ex)] <- f_binom$linkfun(tol_ex)
    eta_cut[eta_cut > f_binom$linkfun(1 - tol_ex)] <-
      f_binom$linkfun(1 - tol_ex)
    mu_cut[mu_cut < tol_ex] <- tol_ex
    mu_cut[mu_cut > 1 - tol_ex] <- 1 - tol_ex
    expect_equal(eta_cut, refmod$family$linkfun(mu_cut), info = info_str)
  } else if (refmod$family$family %in% fam_nms_aug_long &&
             (any(abs(refmod$mu - 0) <= .Machine$double.eps) ||
              any(abs(refmod$mu - 1) <= .Machine$double.eps))) {
    # The degenerate probabilities in `refmod$mu` are probably due to numerical
    # underflow and overflow (for zeros and ones, respectively), so applying the
    # link function would lead to infinite values. Thus, the only sensible (and
    # quickly feasible) check is:
    expect_equal(refmod$mu, refmod$family$linkinv(refmod$eta), info = info_str)
  } else {
    expect_equal(refmod$eta, refmod$family$linkfun(refmod$mu), info = info_str)
  }

  # dis
  if (refmod$family$family == "gaussian") {
    if (is_datafit) {
      expect_identical(refmod$dis, 0, info = info_str)
    } else if (latent_expected) {
      expect_identical(refmod$dis, rep(1.6, nrefdraws_expected),
                       info = info_str)
    } else {
      expect_true(is.vector(refmod$dis, "double"), info = info_str)
      expect_length(refmod$dis, nrefdraws_expected)
      expect_true(all(refmod$dis > 0), info = info_str)
    }
  } else {
    expect_identical(refmod$dis, rep(NA, nrefdraws_expected), info = info_str)
  }

  # y
  ### Not needed because of the more precise test below:
  # expect_true(is.vector(refmod$y, "numeric"), info = info_str)
  # expect_length(refmod$y, nobsv_expected)
  ###
  if (!with_spclformul) {
    y_expected <- dat[[paste("y", mod_nm, fam_nm, sep = "_")]]
    if (!is_datafit && pkg_nm == "brms" && packageVersion("brms") < "2.16.11" &&
        fam_nm == "brnll") {
      # Fixed (as a side effect) by brms PR #1314:
      y_expected <- as.numeric(y_expected)
    }
    if (latent_expected) {
      y_expected <- unname(colMeans(posterior_linpred(fit_expected)))
    }
  } else {
    y_expected <- data_expected[[y_spclformul_new]]
  }
  if (refmod$family$for_augdat) {
    y_expected <- as.factor(y_expected)
    if (fam_nm %in% fam_nms_aug && pkg_nm == "brms") {
      # brms seems to set argument `contrasts`, but this is not important for
      # projpred, so ignore it in the comparison:
      attr(y_expected, "contrasts") <- attr(refmod$y, "contrasts")
    }
  }
  expect_identical(refmod$y, y_expected, info = info_str)

  # loglik
  if (!is_datafit) {
    expect_true(is.matrix(refmod$loglik), info = info_str)
    expect_type(refmod$loglik, "double")
    expect_identical(dim(refmod$loglik), c(nrefdraws_expected, nobsv_expected),
                     info = info_str)
  } else {
    expect_null(refmod$loglik, info = info_str)
  }

  # intercept
  expect_type(refmod$intercept, "logical")
  expect_length(refmod$intercept, 1)
  expect_false(is.na(refmod$intercept), info = info_str)
  # As long as models without an intercept are not supported by projpred:
  expect_true(refmod$intercept, info = info_str)

  # proj_predfun
  expect_type(refmod$proj_predfun, "closure")

  # fetch_data
  expect_type(refmod$fetch_data, "closure")
  if (!is_gamm) {
    # TODO (GAMMs): Adapt this to GAMMs.
    if (latent_expected) {
      data_expected[[stdized_lhs$y_nm]] <- y_expected
    }
    if ((!is_datafit && pkg_nm != "brms") ||
        (is_datafit && (pkg_nm == "brms" || fam_nm != "binom"))) {
      expect_identical(refmod$fetch_data(), data_expected, info = info_str)
    } else if (!is_datafit && pkg_nm == "brms") {
      if (!with_spclformul) {
        refdat_colnms <- as.character(
          attr(terms(formul_expected), "variables")
        )[-1]
        refdat_colnms <- sub(".*\\|[[:blank:]]*", "", refdat_colnms)
        refdat_colnms <- sub("s\\((.*)\\)", "\\1", refdat_colnms)
        if (!all(wobs_expected == 1)) {
          refdat_colnms <- c(head(refdat_colnms, 1),
                             "wobs_col",
                             tail(refdat_colnms, -1))
        }
        refdat_colnms <- sub("^offset\\((.*)\\)$", "\\1", refdat_colnms)
        if (latent_expected) {
          # Re-order:
          refdat_colnms <- c(sub("^\\.", "", stdized_lhs$y_nm),
                             setdiff(refdat_colnms, stdized_lhs$y_nm),
                             stdized_lhs$y_nm)
        }
        refdat_ch <- data_expected[, refdat_colnms, drop = FALSE]
        expect_equal(refmod$fetch_data(), refdat_ch, check.attributes = FALSE,
                     info = info_str)
      } else {
        # TODO: The check in this case is not optimal yet because it subsets
        # `refmod$fetch_data()` to only those columns which also exist in the
        # rebuilt dataset. It would be more desirable to check *all* columns of
        # `refmod$fetch_data()`.
        refdat_ch <- model.frame(formul_expected, data = data_expected)
        if (!all(wobs_expected == 1)) {
          refdat_ch <- cbind(refdat_ch, wobs_col = data_expected$wobs_col)
        }
        names(refdat_ch) <- sub("^offset\\((.*)\\)$", "\\1", names(refdat_ch))
        expect_equal(refmod$fetch_data()[, names(refdat_ch), drop = FALSE],
                     refdat_ch, check.attributes = FALSE, info = info_str)
      }
    } else if (is_datafit && pkg_nm != "brms" && fam_nm == "binom") {
      refdat_ch <- data_expected
      y_nm <- paste("y", mod_nm, fam_nm, sep = "_")
      refdat_ch$dummy_nm <- refdat_ch$wobs_col - refdat_ch[, y_nm]
      names(refdat_ch)[names(refdat_ch) == "dummy_nm"] <- paste("wobs_col -",
                                                                y_nm)
      expect_identical(refmod$fetch_data(), refdat_ch, info = info_str)
    } else {
      stop("This case should not occur.")
    }
  }

  # wobs
  ### Not needed because of the more precise test below:
  # expect_true(is.vector(refmod$wobs, "numeric"), info = info_str)
  # expect_length(refmod$wobs, nobsv_expected)
  # expect_true(all(refmod$wobs > 0), info = info_str)
  ###
  if (!is_gamm) {
    # TODO (GAMMs): Adapt the expected observation weights to GAMMs.
    expect_identical(refmod$wobs, wobs_expected, info = info_str)
  }

  # wsample
  expect_true(is.vector(refmod$wsample, "double"), info = info_str)
  expect_length(refmod$wsample, nrefdraws_expected)
  expect_true(all(refmod$wsample > 0), info = info_str)

  # offset
  ### Not needed because of the more precise test below:
  # expect_true(is.vector(refmod$offset, "double"), info = info_str)
  # expect_length(refmod$offset, nobsv_expected)
  ###
  expect_identical(refmod$offset, offs_expected, info = info_str)

  # cvfun
  expect_type(refmod$cvfun, "closure")

  # cvfits
  expect_null(refmod$cvfits, info = info_str)
  expect_false(is.null(refmod$cvfun) && is.null(refmod$cvfits), info = info_str)

  # extract_model_data
  expect_type(refmod$extract_model_data, "closure")

  # ref_predfun
  expect_type(refmod$ref_predfun, "closure")

  # cvrefbuilder
  expect_type(refmod$cvrefbuilder, "closure")

  # yOrig
  if (latent_expected) {
    yOrig_expected <- data_expected[[sub("^\\.", "", stdized_lhs$y_nm)]]
    if (!is.null(refmod$family$cats) && !is.factor(yOrig_expected)) {
      yOrig_expected <- as.factor(yOrig_expected)
    }
    if (pkg_nm == "brms") {
      # brms seems to set argument `contrasts`, but this is not important for
      # projpred, so ignore it in the comparison:
      attr(yOrig_expected, "contrasts") <- attr(refmod$yOrig, "contrasts")
    }
  } else {
    yOrig_expected <- y_expected
  }
  expect_identical(refmod$yOrig, yOrig_expected, info = info_str)

  return(invisible(TRUE))
}

# A helper function for testing the structure of a list of fits for the same
# single submodel, with that submodel coming from a traditional projection.
#
# @inheritParams submodl_tester
# @param seq_extensive_tests A numeric vector of indexes from
#   `seq_along(submodl_totest)` indicating those elements of `submodl_totest`
#   for which more extensive tests should be conducted.
# @param nobsv The number of (possibly augmented) observations.
#
# @return `TRUE` (invisible).
submodl_tester_trad <- function(
    submodl_totest,
    nprjdraws_expected,
    sub_formul,
    sub_data,
    sub_fam,
    has_grp = formula_contains_group_terms(sub_formul[[1]]),
    has_add = formula_contains_additive_terms(sub_formul[[1]]),
    wobs_expected = wobs_tst,
    solterms_vsel_L1_search = NULL,
    with_offs = FALSE,
    augdat_cats = NULL,
    allow_w_zero = FALSE,
    check_y_from_resp = TRUE,
    seq_extensive_tests = seq_extensive_tests,
    nobsv = nobsv,
    info_str
) {
  from_vsel_L1_search <- !is.null(solterms_vsel_L1_search)

  if (!has_grp && !has_add) {
    sub_x_expected <- model.matrix(sub_formul[[1]], data = sub_data)
    sub_x_expected <- sub_x_expected[
      , colnames(sub_x_expected) != "(Intercept)", drop = FALSE
    ]
    ncoefs <- ncol(sub_x_expected)
    if (from_vsel_L1_search) {
      # Unfortunately, model.matrix() uses terms() and there seems to be no way
      # to set `keep.order = TRUE` in that internal terms() call. Thus, we have
      # to reorder the columns manually:
      if (length(solterms_vsel_L1_search) > 0) {
        terms_contr_expd <- lapply(solterms_vsel_L1_search, function(term_crr) {
          if (!is.factor(sub_data[[term_crr]])) {
            return(term_crr)
          } else {
            return(paste0(term_crr, levels(sub_data[[term_crr]])[-1]))
          }
        })
        terms_contr_expd <- unlist(terms_contr_expd)
        colnms_x <- colnames(sub_x_expected)
        colnms_x <- unlist(lapply(colnms_x, function(colnm_x) {
          if (!colnm_x %in% terms_contr_expd) {
            colon_found <- gregexpr(":", colnm_x)
            if (length(colon_found) != 1 || colon_found == -1) {
              stop("The following code is not general enough. It needs to be ",
                   "adapted.")
            }
            return(paste(rev(strsplit(colnm_x, ":")[[1]]), collapse = ":"))
          } else {
            return(colnm_x)
          }
        }))
        colnames(sub_x_expected) <- colnms_x
        sub_x_expected <- sub_x_expected[
          ,
          colnames(sub_x_expected)[order(match(colnames(sub_x_expected),
                                               terms_contr_expd))],
          drop = FALSE
        ]
      }
    }
    subfit_nms <- c("alpha", "beta", "w", "formula", "x", "y")
    if (from_vsel_L1_search) {
      subfit_nms <- setdiff(subfit_nms, "y")
    }
    for (j in seq_along(submodl_totest)) {
      expect_s3_class(submodl_totest[[!!j]], "subfit")
      expect_type(submodl_totest[[!!j]], "list")
      expect_named(submodl_totest[[!!j]], subfit_nms, info = info_str)

      if (j %in% seq_extensive_tests) {
        expect_true(is.vector(submodl_totest[[!!j]]$alpha, "double"),
                    info = info_str)
        expect_length(submodl_totest[[!!j]]$alpha, 1)

        if (ncoefs > 0 || !from_vsel_L1_search) {
          expect_true(is.matrix(submodl_totest[[!!j]]$beta), info = info_str)
          expect_true(is.numeric(submodl_totest[[!!j]]$beta), info = info_str)
          expect_identical(dim(submodl_totest[[!!j]]$beta), c(ncoefs, 1L),
                           info = info_str)
        } else if (ncoefs == 0) {
          expect_null(submodl_totest[[!!j]]$beta, info = info_str)
        }

        if (!from_vsel_L1_search) {
          expect_true(is.matrix(submodl_totest[[!!j]]$w), info = info_str)
          expect_type(submodl_totest[[!!j]]$w, "double")
          expect_identical(dim(submodl_totest[[!!j]]$w), c(nobsv, 1L),
                           info = info_str)
        } else {
          expect_true(is.vector(submodl_totest[[!!j]]$w, "double"),
                      info = info_str)
          expect_length(submodl_totest[[!!j]]$w, nobsv)
        }
        if (allow_w_zero) {
          expect_true(all(submodl_totest[[!!j]]$w >= 0), info = info_str)
        } else {
          expect_true(all(submodl_totest[[!!j]]$w > 0), info = info_str)
        }
        if (sub_fam == "gaussian") {
          # Note: For non-Gaussian families, a comparison of
          # `submodl_totest[[j]]$w` with `wobs_expected` doesn't make sense
          # since `glm_ridge(<...>)$w` contains the weights of the
          # pseudo-Gaussian observations as calculated in pseudo_data().
          expect_equal(as.vector(submodl_totest[[!!j]]$w),
                       wobs_expected %||% rep(1, nobsv),
                       info = info_str)
        }

        expect_s3_class(submodl_totest[[!!j]]$formula, "formula")
        if (!grepl(":", as.character(submodl_totest[[j]]$formula)[3])) {
          expect_equal(submodl_totest[[!!j]]$formula, sub_formul[[!!j]],
                       info = info_str)
        } else {
          # The order of terms as well as the order of individual terms within
          # ":" interaction terms might be changed in the reference model:
          expect_equal(submodl_totest[[!!j]]$formula[[2]],
                       sub_formul[[!!j]][[2]],
                       info = info_str)
          trms_to_test <- labels(terms(submodl_totest[[j]]$formula))
          trms_ch <- labels(terms(sub_formul[[j]]))
          expect_true(setequal(c(trms_to_test, revIA(trms_to_test)),
                               c(trms_ch, revIA(trms_ch))),
                      info = info_str)
        }

        x_to_test <- submodl_totest[[j]]$x
        x_ch <- sub_x_expected
        if (!identical(dimnames(x_to_test)[[2]],
                       dimnames(x_ch)[[2]])) {
          # Try reversing the order of individual terms within ":" interaction
          # terms:
          dimnames(x_ch)[[2]][grep(":", dimnames(x_ch)[[2]])] <- revIA(
            dimnames(x_ch)[[2]]
          )
        }
        expect_identical(x_to_test, x_ch,
                         info = info_str)

        if (!from_vsel_L1_search) {
          y_ch <- setNames(eval(str2lang(as.character(sub_formul[[j]])[2]),
                                sub_data),
                           seq_len(nobsv))
          if (!is.null(augdat_cats)) {
            stopifnot(is.factor(y_ch), identical(levels(y_ch), c("0", "1")))
            y_ch <- setNames(as.integer(y_ch) - 1L, names(y_ch))
          }
          expect_identical(submodl_totest[[!!j]]$y, y_ch, info = info_str)
        }
      }
    }
  } else if (has_grp && !has_add) {
    if (sub_fam == "gaussian") {
      for (j in seq_along(submodl_totest)) {
        expect_s4_class(submodl_totest[[!!j]], "lmerMod")
      }
    } else {
      for (j in seq_along(submodl_totest)) {
        expect_s4_class(submodl_totest[[!!j]], "glmerMod")
      }
    }

    sub_trms_for_mf <- labels(terms(sub_formul[[1]]))
    sub_trms_for_mf <- sub(
      "^[[:blank:]]*(.*)[[:blank:]]*\\|[[:blank:]]*(.*)[[:blank:]]*$",
      "\\1 + \\2",
      sub_trms_for_mf
    )
    sub_trms_for_mf <- unique(sub_trms_for_mf)
    sub_formul_for_mf <- as.formula(paste(
      "~", paste(sub_trms_for_mf, collapse = " + ")
    ))
    sub_mf_expected <- model.frame(sub_formul_for_mf, data = sub_data)

    sub_formul_for_mm <- grep("\\|", labels(terms(sub_formul[[1]])),
                              value = TRUE, invert = TRUE)
    if (length(sub_formul_for_mm) == 0) {
      sub_formul_for_mm <- "1"
    }
    sub_formul_for_mm <- as.formula(paste(
      "~", paste(sub_formul_for_mm, collapse = " + ")
    ))
    mm_expected <- model.matrix(sub_formul_for_mm, data = sub_mf_expected)
    attr(mm_expected, "msgScaleX") <- character()

    for (j in seq_extensive_tests) {
      # formula
      if (!with_offs) {
        sub_formul_expected <- flatten_formula(sub_formul[[j]])
      } else {
        sub_formul_expected <- sub_formul[[j]]
      }
      expect_equal(submodl_totest[[!!j]]@call[["formula"]],
                   sub_formul_expected,
                   info = info_str)

      # resp
      if (!with_offs) {
        offs_expected <- numeric(nobsv)
      } else {
        offs_expected <- offs_tst
      }
      expect_identical(submodl_totest[[!!j]]@resp$offset,
                       offs_expected,
                       info = info_str)
      if (!is.null(wobs_expected)) {
        if (is.matrix(wobs_expected)) {
          wobs_expected_j <- wobs_expected[, j]
        } else {
          wobs_expected_j <- wobs_expected
        }
        expect_equal(submodl_totest[[!!j]]@resp$weights,
                     wobs_expected_j,
                     info = info_str)
      } else {
        expect_equal(submodl_totest[[!!j]]@resp$weights,
                     rep(1, nobsv),
                     info = info_str)
      }
      y_from_resp <- submodl_totest[[j]]@resp$y
      if (!is.null(augdat_cats)) {
        y_from_resp <- as.factor(y_from_resp)
      }
      if (check_y_from_resp) {
        expect_equal(y_from_resp,
                     eval(str2lang(as.character(sub_formul[[!!j]])[2]),
                          sub_data),
                     info = info_str)
      }

      # frame
      expect_identical(submodl_totest[[!!j]]@frame,
                       model.frame(submodl_totest[[!!j]]),
                       info = info_str)
      if (check_y_from_resp) {
        expect_equal(
          submodl_totest[[!!j]]@frame[[
            grep("y_|ybinprop", names(submodl_totest[[!!j]]@frame),
                 value = TRUE)
          ]],
          y_from_resp,
          info = info_str
        )
      }
      if (!is.null(wobs_expected)) {
        expect_equal(structure(submodl_totest[[!!j]]@frame$`(weights)`,
                               nobs_orig = NULL,
                               class = NULL),
                     submodl_totest[[!!j]]@resp$weights,
                     info = info_str)
      } else {
        expect_null(submodl_totest[[!!j]]@frame$`(weights)`,
                    info = info_str)
      }
      if (with_offs) {
        expect_equal(submodl_totest[[!!j]]@frame$`offset(offs_col)`,
                     offs_expected,
                     info = info_str)
      }
      frame_nms <- grep("y_|ybinprop|^\\(weights\\)$|^offset\\(.*\\)$",
                        names(submodl_totest[[j]]@frame),
                        value = TRUE,
                        invert = TRUE)
      expect_setequal(frame_nms, names(sub_mf_expected))
      expect_equal(
        submodl_totest[[!!j]]@frame[frame_nms],
        sub_mf_expected[frame_nms],
        info = info_str
      )

      # model.matrix()
      expect_identical(model.matrix(submodl_totest[[!!j]]), mm_expected,
                       info = info_str)

      # flist
      expect_type(submodl_totest[[!!j]]@flist, "list")
      expect_length(submodl_totest[[!!j]]@flist, length(nlvl_ran))
      z_nms <- intersect(names(submodl_totest[[j]]@flist),
                         names(submodl_totest[[j]]@frame))
      expect_identical(submodl_totest[[!!j]]@flist[z_nms],
                       as.list(submodl_totest[[!!j]]@frame[z_nms]),
                       info = info_str)

      # coef()
      coefs_crr <- coef(submodl_totest[[j]])
      expect_type(coefs_crr, "list")
      expect_length(coefs_crr, length(nlvl_ran))
      for (zz in seq_len(length(nlvl_ran))) {
        expect_true(is.data.frame(coefs_crr[[zz]]),
                    info = paste(info_str, j, zz, sep = "__"))
        expect_identical(nrow(coefs_crr[[zz]]),
                         unname(nlvl_ran[zz]),
                         info = paste(info_str, j, zz, sep = "__"))
        expect_true(all(sapply(coefs_crr[[zz]], is.numeric)),
                    info = paste(info_str, j, zz, sep = "__"))
      }
    }
  } else if (!has_grp && has_add) {
    for (j in seq_along(submodl_totest)) {
      expect_s3_class(submodl_totest[[!!j]], "gam")
    }
    # TODO (GAMs): Add more expectations for GAMs.
  } else if (has_grp && has_add) {
    for (j in seq_along(submodl_totest)) {
      expect_s3_class(submodl_totest[[!!j]], "gamm4")
    }
    # TODO (GAMMs): Add more expectations for GAMMs.
  }

  return(invisible(TRUE))
}

# A helper function for testing the structure of a list of fits for the same
# single submodel, with that submodel coming from an augmented-data projection.
#
# @inheritParams submodl_tester_trad
#
# @return `TRUE` (invisible).
submodl_tester_aug <- function(
    submodl_totest,
    nprjdraws_expected,
    sub_formul,
    sub_data,
    sub_fam,
    has_grp = formula_contains_group_terms(sub_formul[[1]]),
    has_add = formula_contains_additive_terms(sub_formul[[1]]),
    wobs_expected = wobs_tst,
    solterms_vsel_L1_search = NULL,
    with_offs = FALSE,
    augdat_cats = NULL,
    allow_w_zero = FALSE,
    check_y_from_resp = TRUE,
    seq_extensive_tests = seq_extensive_tests,
    nobsv = nobsv,
    info_str
) {
  sub_formul <- sub_formul[[1]]
  if (has_add) {
    stop("This case should not occur (yet).")
  } else if (!has_grp) {
    if (sub_fam %in% c("cumulative", "cumulative_rstanarm")) {
      for (j in seq_along(submodl_totest)) {
        expect_s3_class(submodl_totest[[!!j]], "polr")
      }
      for (j in seq_extensive_tests) {
        # coef()
        coefs_crr <- coef(submodl_totest[[j]])
        expect_true(is.vector(coefs_crr, "numeric"), info = info_str)
        expect_named(
          coefs_crr,
          grep("Intercept", colnames(model.matrix(sub_formul, data = sub_data)),
               value = TRUE, invert = TRUE),
          info = info_str
        )

        # zeta
        zeta_crr <- submodl_totest[[j]]$zeta
        expect_true(is.vector(zeta_crr, "numeric"), info = info_str)
        expect_named(
          zeta_crr,
          paste(head(augdat_cats, -1), tail(augdat_cats, -1), sep = "|"),
          info = info_str
        )

        # lev
        expect_identical(submodl_totest[[!!j]]$lev, augdat_cats,
                         info = info_str)

        # method
        expect_identical(submodl_totest[[!!j]]$method, "logistic",
                         info = info_str)
      }
    } else if (sub_fam == "categorical") {
      for (j in seq_along(submodl_totest)) {
        expect_s3_class(submodl_totest[[!!j]], c("multinom", "nnet"))
      }
      for (j in seq_extensive_tests) {
        # coef()
        coefs_crr <- coef(submodl_totest[[j]])
        expect_true(is.matrix(coefs_crr), info = info_str)
        expect_true(is.numeric(coefs_crr), info = info_str)
        expect_identical(
          dimnames(coefs_crr),
          list(tail(augdat_cats, -1),
               colnames(model.matrix(sub_formul, data = sub_data))),
          info = info_str
        )

        # lev
        expect_identical(submodl_totest[[!!j]]$lev, augdat_cats,
                         info = info_str)
      }
    } else {
      stop("Unexpected `sub_fam` value of `", sub_fam, "`. Info: ", info_str)
    }
  } else if (has_grp) {
    grp_trms_for_coef <- extract_terms_response(sub_formul)$group_terms
    grp_trms_for_coef <- sub("[[:blank:]]*\\|.*$", "", grp_trms_for_coef)
    coef_nms <- strsplit(grp_trms_for_coef, "[[:blank:]]*\\+[[:blank:]]*")
    coef_nms <- union("1", coef_nms)
    coef_nms <- sub("^1$", "(Intercept)", unlist(coef_nms))
    if (sub_fam %in% c("cumulative", "cumulative_rstanarm")) {
      for (j in seq_along(submodl_totest)) {
        expect_s3_class(submodl_totest[[!!j]], "clmm")
      }
      for (j in seq_extensive_tests) {
        # alpha
        alpha_crr <- submodl_totest[[j]]$alpha
        expect_true(is.vector(alpha_crr, "numeric"), info = info_str)
        expect_named(
          alpha_crr,
          paste(head(augdat_cats, -1), tail(augdat_cats, -1), sep = "|"),
          info = info_str
        )

        # beta
        coefs_crr <- submodl_totest[[j]]$beta
        expect_true(is.vector(coefs_crr, "numeric"), info = info_str)
        ### A quick-and-dirty workaround to get rid of group-level terms:
        stopifnot(identical(trms_grp, c("(xco.1 | z.1)")))
        sub_formul_no_grp <- update(sub_formul,
                                    . ~ . - (1 | z.1) - (xco.1 | z.1))
        ###
        coefs_nms_expected <- grep(
          "Intercept",
          colnames(model.matrix(sub_formul_no_grp, data = sub_data)),
          value = TRUE, invert = TRUE
        )
        if (length(coefs_nms_expected)) {
          expect_named(coefs_crr, coefs_nms_expected, info = info_str)
        } else {
          expect_length(coefs_crr, 0)
        }

        # ordinal::ranef()
        ranef_crr <- ordinal::ranef(submodl_totest[[j]])
        expect_type(ranef_crr, "list")
        expect_named(ranef_crr, "z.1", info = info_str)
        expect_true(is.data.frame(ranef_crr[["z.1"]]), info = info_str)
        expect_named(ranef_crr[["z.1"]], coef_nms, info = info_str)
        expect_true(is.vector(ranef_crr[["z.1"]][["(Intercept)"]], "numeric"),
                    info = info_str)
        expect_length(ranef_crr[["z.1"]][["(Intercept)"]], nlevels(dat$z.1))
        expect_true(is.vector(ranef_crr[["z.1"]][["xco.1"]], "numeric"),
                    info = info_str)
        expect_length(ranef_crr[["z.1"]][["xco.1"]], nlevels(dat$z.1))

        # ordinal::VarCorr()
        VarCorr_crr <- ordinal::VarCorr(submodl_totest[[j]])
        expect_type(VarCorr_crr, "list")
        expect_named(VarCorr_crr, "z.1", info = info_str)
        expect_true(is.matrix(VarCorr_crr[["z.1"]]), info = info_str)
        expect_true(is.numeric(VarCorr_crr[["z.1"]]), info = info_str)
        expect_identical(dimnames(VarCorr_crr[["z.1"]]),
                         replicate(2, coef_nms, simplify = FALSE),
                         info = info_str)
        expect_true(is.vector(attr(VarCorr_crr[["z.1"]], "stddev"), "numeric"),
                    info = info_str)
        expect_named(attr(VarCorr_crr[["z.1"]], "stddev"), coef_nms,
                     info = info_str)
        expect_true(is.matrix(attr(VarCorr_crr[["z.1"]], "correlation")),
                    info = info_str)
        expect_true(is.numeric(attr(VarCorr_crr[["z.1"]], "correlation")),
                    info = info_str)
        expect_identical(dimnames(attr(VarCorr_crr[["z.1"]], "correlation")),
                         replicate(2, coef_nms, simplify = FALSE),
                         info = info_str)

        # formula()
        formula_crr <- formula(submodl_totest[[j]])
        # Unfortunately, there are some minor caveats to take care of when
        # comparing `formula_crr` with `sub_formul`: (i) the intercept (`1`)
        # needs to included explicitly, (ii) group-level terms need to be
        # "flattened" (but flatten_formula() would omit offset terms; adding an
        # argument `incl_offs` to flatten_formula() could solve this, but the
        # internal update() call would then still move offset terms to the end):
        sub_trms <- labels(terms(sub_formul))
        nongrp_trms <- grep("\\|", sub_trms, value = TRUE, invert = TRUE)
        grp_trms <- grep("\\|", sub_trms, value = TRUE)
        grp_trms_flat <- flatten_group_terms(grp_trms)
        offs_trms <- if (with_offs) "offset(offs_col)" else NULL
        expect_identical(
          formula_crr,
          # Add the intercept explicitly and use the same environment:
          as.formula(paste(as.character(sub_formul[[2]]), "~ 1 +",
                           paste(c(nongrp_trms, offs_trms, grp_trms_flat),
                                 collapse = " + ")),
                     env = environment(formula_crr)),
          info = info_str
        )

        # xlevels
        xlevels_crr <- submodl_totest[[j]]$xlevels
        xca_nms <- grep("^xca\\.", labels(terms(sub_formul)), value = TRUE)
        if (length(xca_nms)) {
          expect_identical(xlevels_crr, lapply(dat[xca_nms], levels),
                           info = info_str)
        } else {
          expect_null(xlevels_crr, info = info_str)
        }
      }
    } else if (sub_fam == "categorical") {
      for (j in seq_along(submodl_totest)) {
        expect_s3_class(submodl_totest[[!!j]],
                        c("mmblogit", "mblogit", "mmclogit", "mclogit", "lm"))
      }
      coef_nms <- sub("^\\(Intercept\\)$", "1", coef_nms)
      for (j in seq_extensive_tests) {
        # coefmat
        coefs_crr <- submodl_totest[[j]]$coefmat
        expect_true(is.matrix(coefs_crr), info = info_str)
        expect_true(is.numeric(coefs_crr), info = info_str)
        ### A quick-and-dirty workaround to get rid of group-level terms:
        stopifnot(identical(trms_grp, c("(xco.1 | z.1)")))
        sub_formul_no_grp <- update(sub_formul,
                                    . ~ . - (1 | z.1) - (xco.1 | z.1))
        ###
        expect_identical(
          dimnames(coefs_crr),
          list("Response categories" = tail(augdat_cats, -1),
               "Predictors" = colnames(model.matrix(sub_formul_no_grp,
                                                    data = sub_data))),
          info = info_str
        )

        # random.effects
        ranef_crr <- submodl_totest[[j]]$random.effects
        expect_type(ranef_crr, "list")
        expect_length(ranef_crr, 1)
        expect_named(ranef_crr, NULL, info = info_str)
        expect_true(is.matrix(ranef_crr[[1]]), info = info_str)
        expect_true(is.numeric(ranef_crr[[1]]), info = info_str)
        if (packageVersion("Matrix") >= "1.5-0") {
          expect_null(dimnames(ranef_crr[[1]]), info = info_str)
        } else {
          expect_identical(dimnames(ranef_crr[[1]]),
                           replicate(2, NULL, simplify = FALSE),
                           info = info_str)
        }
        expect_identical(dim(ranef_crr[[1]]),
                         c(nthres * length(coef_nms) * nlevels(dat$z.1), 1L),
                         info = info_str)

        # VarCov
        VarCorr_crr <- submodl_totest[[j]]$VarCov
        expect_type(VarCorr_crr, "list")
        expect_named(VarCorr_crr, "z.1", info = info_str)
        expect_true(is.matrix(VarCorr_crr[["z.1"]]), info = info_str)
        expect_true(is.numeric(VarCorr_crr[["z.1"]]), info = info_str)
        coef_nms_y <- unlist(lapply(coef_nms, function(coef_nm) {
          paste(tail(augdat_cats, -1), coef_nm, sep = "~")
        }))
        expect_identical(dimnames(VarCorr_crr[["z.1"]]),
                         replicate(2, coef_nms_y, simplify = FALSE),
                         info = info_str)

        # D
        D_crr <- submodl_totest[[j]]$D
        expect_identical(D_crr, contrasts(as.factor(augdat_cats)),
                         info = info_str)

        # random
        random_crr <- submodl_totest[[j]]$random
        expect_type(random_crr, "list")
        expect_length(random_crr, 1)
        expect_named(random_crr, NULL, info = info_str)
        expect_named(random_crr[[1]], c("formula", "groups"), info = info_str)
        coef_nms_no_icpt <- setdiff(coef_nms, "1")
        if (length(coef_nms_no_icpt)) {
          expect_equal(
            random_crr[[1]]$formula,
            as.formula(paste("~", paste(coef_nms_no_icpt, collapse = "+"))),
            info = info_str
          )
        } else {
          expect_equal(random_crr[[1]]$formula, ~ 1, info = info_str)
        }
        expect_identical(random_crr[[1]]$groups, "z.1", info = info_str)

        # groups
        groups_crr <- submodl_totest[[j]]$groups
        expect_type(groups_crr, "list")
        expect_named(groups_crr, "z.1", info = info_str)
        expect_identical(levels(groups_crr[["z.1"]]), levels(dat$z.1),
                         info = info_str)

        # xlevels
        xlevels_crr <- submodl_totest[[j]]$xlevels
        xca_nms <- grep("^xca\\.", labels(terms(sub_formul)), value = TRUE)
        if (length(xca_nms)) {
          expect_identical(xlevels_crr, lapply(dat[xca_nms], levels),
                           info = info_str)
        } else {
          expect_null(xlevels_crr, info = info_str)
        }
      }
    } else {
      stop("Unexpected `sub_fam` value of `", sub_fam, "`. Info: ", info_str)
    }
  }

  return(invisible(TRUE))
}

# A helper function for testing the structure of a list of fits for the same
# single submodel.
#
# @param submodl_totest The `submodl` object (a list of fits for a single
#   submodel, with one fit per projected draw) to test.
# @param nprjdraws_expected A single numeric value giving the expected number of
#   projected draws.
# @param sub_formul A list of formulas for the submodel (with one element per
#   projected draw).
# @param sub_data The dataset used for fitting the submodel.
# @param sub_fam A single character string giving the submodel's family.
# @param has_grp A single logical value indicating whether the fits in
#   `submodl_totest` are expected to be of class `"lmerMod"` or `"glmerMod"`
#   (if, at the same time, `has_add` is `FALSE`).
# @param has_add A single logical value indicating whether the fits in
#   `submodl_totest` are expected to be of class `"gam"` or `"gamm4"` (depending
#   on whether the submodel is non-multilevel or multilevel, respectively).
# @param wobs_expected The expected numeric vector of observation weights.
# @param solterms_vsel_L1_search If `submodl_totest` comes from the L1
#   `search_path` of an object of class `"vsel"`, provide here the solution
#   terms. Otherwise, use `NULL`.
# @param with_offs A single logical value indicating whether `submodl_totest` is
#   expected to include offsets (`TRUE`) or not (`FALSE`).
# @param augdat_cats A character vector of response levels in case of the
#   augmented-data projection. Needs to be `NULL` for the traditional and the
#   latent projection.
# @param allow_w_zero A single logical value indicating whether observation
#   weights are allowed to have a value of zero (`TRUE`) or not (`FALSE`).
# @param check_y_from_resp A single logical value indicating whether to check
#   elements `submodl_totest[[j]]@resp$y` for GLMMs (`TRUE`) or not (`FALSE`).
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
submodl_tester <- function(
    submodl_totest,
    nprjdraws_expected,
    sub_formul,
    sub_data,
    sub_fam,
    has_grp = formula_contains_group_terms(sub_formul[[1]]),
    has_add = formula_contains_additive_terms(sub_formul[[1]]),
    wobs_expected = wobs_tst,
    solterms_vsel_L1_search = NULL,
    with_offs = FALSE,
    augdat_cats = NULL,
    allow_w_zero = FALSE,
    check_y_from_resp = TRUE,
    info_str
) {
  expect_type(submodl_totest, "list")
  expect_length(submodl_totest, nprjdraws_expected)

  seq_extensive_tests <- unique(round(
    seq(1, length(submodl_totest),
        length.out = min(length(submodl_totest), nclusters_pred_tst))
  ))

  if (!is.null(augdat_cats)) {
    nobsv <- nobsv * length(augdat_cats)
    expect_length(sub_formul, 1)
    sub_formul <- replicate(nprjdraws_expected, sub_formul[[1]],
                            simplify = FALSE)
  }

  if (!is.null(augdat_cats) && sub_fam %in% fam_nms_aug_long) {
    submodl_tester_aug(
      submodl_totest = submodl_totest,
      nprjdraws_expected = nprjdraws_expected,
      sub_formul = sub_formul,
      sub_data = sub_data,
      sub_fam = sub_fam,
      has_grp = has_grp,
      has_add = has_add,
      wobs_expected = wobs_expected,
      solterms_vsel_L1_search = solterms_vsel_L1_search,
      with_offs = with_offs,
      augdat_cats = augdat_cats,
      allow_w_zero = allow_w_zero,
      check_y_from_resp = check_y_from_resp,
      seq_extensive_tests = seq_extensive_tests,
      nobsv = nobsv,
      info_str = info_str
    )
  } else {
    submodl_tester_trad(
      submodl_totest = submodl_totest,
      nprjdraws_expected = nprjdraws_expected,
      sub_formul = sub_formul,
      sub_data = sub_data,
      sub_fam = sub_fam,
      has_grp = has_grp,
      has_add = has_add,
      wobs_expected = wobs_expected,
      solterms_vsel_L1_search = solterms_vsel_L1_search,
      with_offs = with_offs,
      augdat_cats = augdat_cats,
      allow_w_zero = allow_w_zero,
      check_y_from_resp = check_y_from_resp,
      seq_extensive_tests = seq_extensive_tests,
      nobsv = nobsv,
      info_str = info_str
    )
  }

  return(invisible(TRUE))
}

# A helper function for testing the structure of .get_refdist()'s output.
#
# @param refd Output of .get_refdist().
# @param nprjdraws_expected A single numeric value giving the expected number of
#   projected draws.
# @param clust_expected A single logical value giving the expected value for
#   `refd$clust_used`.
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
refdist_tester <- function(refd,
                           nobsv_expected = nobsv,
                           nprjdraws_expected = nclusters_pred_tst,
                           clust_expected = TRUE,
                           info_str) {
  expect_named(
    refd, c("mu", "var", "dis", "weights", "cl", "wsample_orig", "clust_used"),
    info = info_str
  )
  expect_identical(dim(refd$mu), c(nobsv_expected, nprjdraws_expected),
                   info = info_str)
  expect_identical(dim(refd$var), c(nobsv_expected, nprjdraws_expected),
                   info = info_str)
  expect_true(is.vector(refd$dis) && is.atomic(refd$dis),
              info = info_str)
  expect_length(refd$dis, nprjdraws_expected)
  expect_true(is.vector(refd$weights) && is.atomic(refd$weights),
              info = info_str)
  expect_length(refd$weights, nprjdraws_expected)
  expect_true(is.vector(refd$cl, "integer"), info = info_str)
  expect_length(refd$cl, nrefdraws)
  expect_identical(refd$wsample_orig, rep(1, nrefdraws), info = info_str)
  expect_identical(refd$clust_used, clust_expected, info = info_str)
  return(invisible(TRUE))
}

# A helper function for testing the structure of an expected `"projection"`
# object.
#
# @param p An object of class `"projection"` (at least expected so).
# @param solterms_expected Either a single numeric value giving the expected
#   number of solution terms (not counting the intercept, even for the
#   intercept-only model), a character vector giving the expected solution
#   terms, or `NULL` for not testing the solution terms at all.
# @param nprjdraws_expected A single numeric value giving the expected number of
#   projected draws.
# @param p_type_expected A single logical value giving the expected value for
#   `p$p_type`.
# @param seed_expected The seed which was used for clustering the posterior
#   draws of the reference model.
# @param prjdraw_weights_expected The expected weights for the projected draws
#   or `NULL` for not testing these weights at all.
# @param from_vsel_L1_search A single logical value indicating whether `p` uses
#   the L1 `search_path` of an object of class `"vsel"` for extracting the
#   subfit(s) (`TRUE`) or not (`FALSE`).
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
projection_tester <- function(p,
                              refmod_expected,
                              solterms_expected,
                              nprjdraws_expected,
                              p_type_expected,
                              seed_expected = seed_tst,
                              prjdraw_weights_expected = NULL,
                              from_vsel_L1_search = FALSE,
                              info_str = "") {
  expect_s3_class(p, "projection")
  expect_type(p, "list")
  # Check the names using `ignore.order = FALSE` because an incorrect
  # order would mean that the documentation of project()'s return value
  # would have to be updated:
  expect_named(
    p,
    c("dis", "kl", "weights", "solution_terms", "submodl", "cl_ref",
      "wdraws_ref", "p_type", "refmodel"),
    info = info_str
  )

  # refmodel
  # Note: Extensive tests for `"refmodel"`s and `"datafit"`s may be run via
  # refmodel_tester().
  expect_identical(p$refmodel, refmod_expected, info = info_str)

  # solution_terms
  if (is.numeric(solterms_expected)) {
    expect_length(p$solution_terms, solterms_expected)
    # Same check, but using count_terms_chosen():
    expect_equal(count_terms_chosen(p$solution_terms, add_icpt = TRUE),
                 solterms_expected + 1, info = info_str)
  } else if (is.character(solterms_expected)) {
    expect_identical(p$solution_terms, solterms_expected, info = info_str)
  }

  # submodl
  sub_trms_crr <- p$solution_terms
  if (length(sub_trms_crr) == 0) {
    sub_trms_crr <- as.character(as.numeric(p$refmodel$intercept))
  }
  if (!from_vsel_L1_search) {
    y_nm <- as.character(p$refmodel$formula)[2]
    solterms_vsel_L1_search_crr <- NULL
  } else {
    y_nm <- ""
    solterms_vsel_L1_search_crr <- p$solution_terms
  }
  y_nms <- y_nm
  # For checking for the augmented-data projection case, we use the "unsafer"
  # `p$refmodel$family$for_augdat` here instead of an extra argument such as
  # `augdat_expected` because such an argument already exists in
  # extfam_tester():
  if (!p$refmodel$family$for_augdat) {
    y_nms <- paste0(".", y_nms)
  }
  # A preliminary check for `nprjdraws_expected` (doesn't work for "datafit"s
  # and, because of issue #131, for submodels which are GAMMs):
  sub_formul_crr_rhs <- as.formula(paste(
    "~", paste(sub_trms_crr, collapse = " + ")
  ))
  if (all(grepl("\\+", sub_trms_crr))) {
    # Avoid duplicated terms in the "empty_size" `search_terms` setting:
    sub_formul_crr_rhs <- update(sub_formul_crr_rhs, . ~ .)
  }
  if (!inherits(p$refmodel, "datafit") &&
      !(formula_contains_additive_terms(sub_formul_crr_rhs) &&
        formula_contains_group_terms(sub_formul_crr_rhs))) {
    # Number of projected draws in as.matrix.projection() (note that more
    # extensive tests for as.matrix.projection() may be found in
    # "test_as_matrix.R"):
    expect_identical(suppressWarnings(NROW(as.matrix(p))), nprjdraws_expected,
                     info = info_str)
  }
  if (!p$refmodel$family$for_augdat && nprjdraws_expected > 1) {
    y_nms <- paste0(y_nms, ".", seq_len(nprjdraws_expected))
  }
  sub_formul_crr <- lapply(y_nms, function(y_nm_i) {
    fml_tmp <- as.formula(paste(
      y_nm_i, "~", paste(sub_trms_crr, collapse = " + ")
    ))
    if (all(grepl("\\+", sub_trms_crr))) {
      # Avoid duplicated terms in the "empty_size" `search_terms` setting:
      fml_tmp <- update(fml_tmp, . ~ .)
    }
    return(fml_tmp)
  })
  sub_data_crr <- p$refmodel$fetch_data()
  if (p_type_expected) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
    }
    set.seed(seed_expected)
    clust_ref <- .get_refdist(p$refmodel, nclusters = nprjdraws_expected)
  } else {
    clust_ref <- .get_refdist(p$refmodel, ndraws = nprjdraws_expected)
  }
  if (p$refmodel$family$for_augdat) {
    # Create the augmented dataset:
    y_unqs <- p$refmodel$family$cats
    sub_data_crr <- do.call(rbind, lapply(y_unqs, function(y_unq) {
      sub_data_crr_j <- sub_data_crr
      sub_data_crr_j[[y_nms]] <- y_unq
      return(sub_data_crr_j)
    }))
    sub_data_crr[[y_nms]] <- factor(sub_data_crr[[y_nms]], levels = y_unqs)
  } else {
    for (i in seq_len(nprjdraws_expected)) {
      sub_data_crr[[y_nms[i]]] <- clust_ref$mu[, i]
    }
  }
  if (p$refmodel$family$for_augdat) {
    wobs_expected_crr <- unclass(clust_ref$mu)
  } else {
    wobs_expected_crr <- p$refmodel$wobs
  }
  if (p$refmodel$family$for_augdat) {
    augdat_cats_crr <- p$refmodel$family$cats
  } else {
    augdat_cats_crr <- NULL
  }
  submodl_tester(p$submodl,
                 nprjdraws_expected = nprjdraws_expected,
                 sub_formul = sub_formul_crr,
                 sub_data = sub_data_crr,
                 sub_fam = p$refmodel$family$family,
                 wobs_expected = wobs_expected_crr,
                 solterms_vsel_L1_search = solterms_vsel_L1_search_crr,
                 augdat_cats = augdat_cats_crr,
                 info_str = info_str)

  # dis
  expect_length(p$dis, nprjdraws_expected)

  # kl
  expect_type(p$kl, "double")
  expect_length(p$kl, 1)
  expect_true(!is.na(p$kl), info = info_str)
  expect_gte(p$kl, 0)

  # weights
  expect_length(p$weights, nprjdraws_expected)
  if (!is.null(prjdraw_weights_expected)) {
    expect_identical(p$weights, prjdraw_weights_expected, info = info_str)
  }
  if (nprjdraws_expected == 1) {
    expect_identical(p$weights, 1, info = info_str)
  }

  # cl_ref
  expect_true(is.vector(p$cl_ref, "numeric"), info = info_str)
  expect_length(p$cl_ref, length(p$refmodel$wsample))

  # wdraws_ref
  expect_identical(p$wdraws_ref, rep(1, length(p$refmodel$wsample)),
                   info = info_str)

  # p_type
  expect_identical(p$p_type, p_type_expected, info = info_str)

  return(invisible(TRUE))
}

# A helper function for testing the structure of an expected `"proj_list"`
# object.
#
# @param p An object of (informal) class `"proj_list"` (at least expected so).
# @param len_expected The expected length of `p`.
# @param is_seq A single logical value indicating whether `p` is expected to be
#   sequential (i.e., the number of solution terms increases by 1 from one
#   element of `p` to the next).
# @param extra_tol A single numeric value giving the relative tolerance when
#   checking the monotonicity of the KL divergences. Because this is a
#   *relative* tolerance, 1 is the neutral value.
# @param info_str A single character string giving information to be printed in
#   case of failure.
# @param ... Arguments passed to projection_tester(), apart from
#   projection_tester()'s arguments `p`, `solterms_expected`, and `info_str`.
#
# @return `TRUE` (invisible).
proj_list_tester <- function(p,
                             len_expected = nterms_max_tst + 1L,
                             is_seq = TRUE,
                             extra_tol = 1.1,
                             info_str = "",
                             ...) {
  expect_type(p, "list")
  expect_length(p, len_expected)
  expect_true(.is_proj_list(p), info = info_str)

  for (j in seq_along(p)) {
    if (is_seq) {
      # The j-th element should have j solution terms (not counting the
      # intercept, even for the intercept-only model):
      solterms_expected_crr <- j - 1
    } else {
      solterms_expected_crr <- NULL
    }
    projection_tester(p[[j]],
                      solterms_expected = solterms_expected_crr,
                      info_str = paste(info_str, j, sep = "__"),
                      ...)
  }
  if (is_seq) {
    # For a sequential `"proj_list"` object and training data, `kl` should be
    # non-increasing for increasing model size:
    klseq <- sapply(p, function(x) sum(x$kl))
    expect_true(all(tail(klseq, -1) <= extra_tol * head(klseq, -1)),
                info = info_str)
    ### Too unsafe because `length(klseq)` is usually small:
    # prop_as_expected <- 0.8
    # expect_true(
    #   mean(tail(klseq, -1) <= extra_tol * head(klseq, -1)) >=
    #     prop_as_expected,
    #   info = info_str
    # )
    ###
  }
  return(invisible(TRUE))
}

# A helper function for testing the structure of an object returned by
# proj_linpred().
#
# @param pl An object resulting from a call to proj_linpred().
# @param len_expected The number of `"projection"` objects used for `pl`.
# @param nprjdraws_expected The expected number of projected draws in `pl`.
# @param nobsv_expected The expected number of observations in `pl`.
# @param lpd_null_expected A single logical value indicating whether output
#   element `lpd` is expected to be `NULL` (`TRUE`) or not (`FALSE`).
# @param ncats_nlats_expected A list of length `len_expected`. If the
#   augmented-data projection is expected to have been applied in case of
#   element `j`, then element `j` has to be a single integer value giving the
#   number of response categories or latent response categories (depending on
#   whether the linear predictor was transformed to response scale or not).
#   Else, element `j` has to be `integer()`.
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
pl_tester <- function(pl,
                      len_expected = 1,
                      nprjdraws_expected = nclusters_pred_tst,
                      nobsv_expected = nobsv,
                      lpd_null_expected = FALSE,
                      ncats_nlats_expected = replicate(len_expected,
                                                       integer(),
                                                       simplify = FALSE),
                      info_str) {
  if (len_expected == 1) {
    pl <- list(pl)
  } else {
    expect_type(pl, "list")
    expect_length(pl, len_expected)
  }
  for (j in seq_along(pl)) {
    expect_named(pl[[!!j]], c("pred", "lpd"), info = info_str)
    expect_identical(
      dim(pl[[!!j]]$pred),
      c(nprjdraws_expected, nobsv_expected, ncats_nlats_expected[[!!j]]),
      info = info_str
    )
    if (!lpd_null_expected) {
      expect_identical(dim(pl[[!!j]]$lpd),
                       c(nprjdraws_expected, nobsv_expected),
                       info = info_str)
    } else {
      expect_null(pl[[!!j]]$lpd, info = info_str)
    }
  }
  return(invisible(TRUE))
}

# A helper function for testing the structure of an object returned by
# proj_predict().
#
# @param pp An object resulting from a call to proj_predict().
# @param len_expected The number of `"projection"` objects used for `pp`.
# @param nprjdraws_out_expected The expected number of projected draws in `pp`.
#   In contrast to argument `nprjdraws_expected` of pl_tester(), this also needs
#   to take proj_predict()'s argument `nresample_clusters` into account.
# @param nobsv_expected The expected number of observations in `pp`.
# @param lpd_null_expected A single logical value indicating whether output
#   element `lpd` is expected to be `NULL` (`TRUE`) or not (`FALSE`).
# @param cats_expected A list of length `len_expected`. If the
#   augmented-data projection is expected to have been applied in case of
#   element `j`, then element `j` has to be a character vector giving the
#   response categories. Else, element `j` has to be `NULL`.
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
pp_tester <- function(pp,
                      len_expected = 1,
                      nprjdraws_out_expected = nresample_clusters_default,
                      nobsv_expected = nobsv,
                      cats_expected = replicate(len_expected,
                                                NULL,
                                                simplify = FALSE),
                      info_str) {
  if (len_expected == 1) {
    pp <- list(pp)
  } else {
    expect_type(pp, "list")
    expect_length(pp, len_expected)
  }
  for (j in seq_along(pp)) {
    expect_identical(dim(pp[[!!j]]),
                     c(nprjdraws_out_expected, nobsv_expected),
                     info = info_str)
    expect_identical(attr(pp[[!!j]], "cats"), cats_expected[[j]],
                     info = info_str)
  }
  return(invisible(TRUE))
}

# A helper function for testing the structure of an expected `"vsel"` object.
#
# @param vs An object of class `"vsel"` (at least expected so).
# @param with_cv A single logical value indicating whether `vs` was created by
#   cv_varsel() (`TRUE`) or not (`FALSE`).
# @param refmod_expected The expected `vs$refmodel` object.
# @param dtest_expected If `vs` was created with a non-`NULL` argument `d_test`
#   (which is only possible for varsel()), then this needs to be the expected
#   `vs$d_test` object. Otherwise, this needs to be `NULL`.
# @param solterms_len_expected A single numeric value giving the expected number
#   of solution terms (not counting the intercept, even for the intercept-only
#   model).
# @param method_expected The expected `vs$method` object.
# @param cv_method_expected The expected `vs$cv_method` object.
# @param valsearch_expected The expected `vs$validate_search` object.
# @param cl_search_expected The expected `vs$clust_used_search` object.
# @param cl_eval_expected The expected `vs$clust_used_eval` object.
# @param nprjdraws_search_expected The expected `vs$nprjdraws_search` object.
# @param nprjdraws_eval_expected The expected `vs$nprjdraws_eval` object.
# @param seed_expected The seed which was used for clustering the posterior
#   draws of the reference model.
# @param nloo_expected Only relevant if `with_cv` is `TRUE`. The value which was
#   used for argument `nloo` of cv_varsel().
# @param search_trms_empty_size A single logical value indicating whether
#   `search_terms` was constructed in a way that causes a model size to be
#   without candidate models.
# @param extra_tol A single numeric value giving the relative tolerance when
#   checking the monotonicity of the KL divergences. Because this is a
#   *relative* tolerance, 1 is the neutral value.
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
vsel_tester <- function(
    vs,
    with_cv = FALSE,
    from_datafit = FALSE,
    refmod_expected,
    dtest_expected = NULL,
    solterms_len_expected,
    method_expected,
    cv_method_expected = NULL,
    valsearch_expected = NULL,
    cl_search_expected = !from_datafit,
    cl_eval_expected = !from_datafit,
    nprjdraws_search_expected = if (!from_datafit) nclusters_tst else 1L,
    nprjdraws_eval_expected = if (!from_datafit) nclusters_pred_tst else 1L,
    seed_expected = seed_tst,
    nloo_expected = NULL,
    search_trms_empty_size = FALSE,
    extra_tol = 1.1,
    info_str = ""
) {
  # Preparations:
  if (with_cv) {
    vsel_nms <- vsel_nms_cv
    vsel_smmrs_sub_nms <- c(vsel_smmrs_sub_nms, "wcv")

    if (is.null(cv_method_expected)) {
      cv_method_expected <- "LOO"
    }
    if (is.null(valsearch_expected)) {
      valsearch_expected <- TRUE
    }

    if (cv_method_expected == "LOO") {
      # Re-order:
      vsel_smmrs_sub_nms[1:2] <- vsel_smmrs_sub_nms[2:1]
      vsel_smmrs_ref_nms[1:2] <- vsel_smmrs_ref_nms[2:1]
    }
  }
  method_expected <- tolower(method_expected)
  if (method_expected == "l1") {
    cl_search_expected <- !from_datafit
    nprjdraws_search_expected <- 1
  }
  if (search_trms_empty_size) {
    # This is the "empty_size" setting, so we have to subtract the skipped model
    # size (see issue #307):
    solterms_len_expected <- solterms_len_expected - 1L
  }

  # Test the general structure of the object:
  expect_s3_class(vs, "vsel")
  expect_type(vs, "list")
  expect_named(vs, vsel_nms, info = info_str)

  # refmodel
  expect_identical(vs$refmodel, refmod_expected, info = info_str)

  # search_path
  expect_type(vs$search_path, "list")
  expect_named(vs$search_path, c("solution_terms", "submodls", "p_sel"),
               info = info_str)
  expect_identical(vs$search_path$solution_terms, vs$solution_terms,
                   info = info_str)
  expect_type(vs$search_path$submodls, "list")
  expect_length(vs$search_path$submodls, solterms_len_expected + 1)
  from_vsel_L1_search <- method_expected == "l1"
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state_old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed_expected)
  if (cl_search_expected) {
    clust_ref <- .get_refdist(vs$refmodel,
                              nclusters = nprjdraws_search_expected)
  } else {
    clust_ref <- .get_refdist(vs$refmodel,
                              ndraws = nprjdraws_search_expected)
  }
  if (!from_vsel_L1_search) {
    y_nm <- as.character(vs$refmodel$formula)[2]
    solterms_vsel_L1_search_crr <- NULL
  } else {
    y_nm <- ""
    solterms_vsel_L1_search_crr <- vs$solution_terms
  }
  y_nms <- y_nm
  # For checking for the augmented-data projection case, we use the "unsafer"
  # `vs$refmodel$family$for_augdat` here instead of an extra argument such as
  # `augdat_expected` because such an argument already exists in
  # extfam_tester():
  if (!vs$refmodel$family$for_augdat) {
    y_nms <- paste0(".", y_nms)
  }
  if (!vs$refmodel$family$for_augdat && nprjdraws_search_expected > 1) {
    y_nms <- paste0(y_nms, ".", seq_len(nprjdraws_search_expected))
  }
  sub_data_crr <- vs$refmodel$fetch_data()
  if (vs$refmodel$family$for_augdat) {
    # Create the augmented dataset:
    y_unqs <- vs$refmodel$family$cats
    sub_data_crr <- do.call(rbind, lapply(y_unqs, function(y_unq) {
      sub_data_crr_j <- sub_data_crr
      sub_data_crr_j[[y_nms]] <- y_unq
      return(sub_data_crr_j)
    }))
    sub_data_crr[[y_nms]] <- factor(sub_data_crr[[y_nms]], levels = y_unqs)
  } else {
    for (i in seq_len(nprjdraws_search_expected)) {
      sub_data_crr[[y_nms[i]]] <- clust_ref$mu[, i]
    }
  }
  if (vs$refmodel$family$for_augdat) {
    wobs_expected_crr <- unclass(clust_ref$mu)
  } else {
    wobs_expected_crr <- vs$refmodel$wobs
  }
  solterms_for_sub <- c(as.character(as.numeric(vs$refmodel$intercept)),
                        vs$solution_terms)
  for (i in seq_along(vs$search_path$submodls)) {
    sub_trms_crr <- head(solterms_for_sub, i)
    if (length(sub_trms_crr) > 1) {
      sub_trms_crr <- setdiff(sub_trms_crr, "1")
    }
    sub_formul_crr <- lapply(y_nms, function(y_nm_i) {
      fml_tmp <- as.formula(paste(
        y_nm_i, "~", paste(sub_trms_crr, collapse = " + ")
      ))
      if (all(grepl("\\+", sub_trms_crr))) {
        # Avoid duplicated terms in the "empty_size" `search_terms` setting:
        fml_tmp <- update(fml_tmp, . ~ .)
      }
      return(fml_tmp)
    })
    if (vs$refmodel$family$for_augdat) {
      augdat_cats_crr <- vs$refmodel$family$cats
    } else {
      augdat_cats_crr <- NULL
    }
    submodl_tester(
      vs$search_path$submodls[[i]],
      nprjdraws_expected = nprjdraws_search_expected,
      sub_formul = sub_formul_crr,
      sub_data = sub_data_crr,
      sub_fam = vs$refmodel$family$family,
      wobs_expected = wobs_expected_crr,
      solterms_vsel_L1_search = solterms_vsel_L1_search_crr,
      augdat_cats = augdat_cats_crr,
      info_str = paste(info_str, i, sep = "__")
    )
  }
  if (!vs$refmodel$family$for_augdat) {
    ncats <- 1L
  } else {
    ncats <- length(vs$refmodel$family$cats)
  }
  nobsv_aug <- nobsv * ncats
  expect_type(vs$search_path$p_sel, "list")
  expect_named(vs$search_path$p_sel,
               c("mu", "var", "dis", "weights", "cl", "wsample_orig",
                 "clust_used"),
               info = info_str)
  expect_true(is.matrix(vs$search_path$p_sel$mu), info = info_str)
  expect_true(is.numeric(vs$search_path$p_sel$mu), info = info_str)
  expect_equal(dim(vs$search_path$p_sel$mu),
               c(nobsv_aug, nprjdraws_search_expected),
               info = info_str)
  if (vs$refmodel$family$for_augdat) {
    expect_s3_class(vs$search_path$p_sel$mu, "augmat")
  }
  expect_true(is.matrix(vs$search_path$p_sel$var), info = info_str)
  if (vs$refmodel$family$family == "gaussian") {
    expect_type(vs$search_path$p_sel$var, "double")
  } else {
    expect_true(all(is.na(vs$search_path$p_sel$var)), info = info_str)
  }
  expect_equal(dim(vs$search_path$p_sel$var),
               c(nobsv_aug, nprjdraws_search_expected),
               info = info_str)
  if (vs$refmodel$family$for_augdat) {
    expect_s3_class(vs$search_path$p_sel$var, "augmat")
  }
  expect_true(is.vector(vs$search_path$p_sel$dis) &&
                is.atomic(vs$search_path$p_sel$dis),
              info = info_str)
  expect_length(vs$search_path$p_sel$dis, nprjdraws_search_expected)
  expect_type(vs$search_path$p_sel$weights, "double")
  expect_length(vs$search_path$p_sel$weights, nprjdraws_search_expected)
  expect_true(is.numeric(vs$search_path$p_sel$cl), info = info_str)
  expect_length(vs$search_path$p_sel$cl, ncol(vs$refmodel$mu))
  expect_identical(vs$search_path$p_sel$wsample_orig,
                   rep(1, ncol(vs$refmodel$mu)), info = info_str)
  expect_identical(vs$search_path$p_sel$clust_used, cl_search_expected,
                   info = info_str)

  # d_test
  if (is.null(dtest_expected)) {
    expect_type(vs$d_test, "list")
    expect_named(vs$d_test, nms_d_test(), info = info_str)
    dtest_type <- cv_method_expected
    if (length(dtest_type) == 0) {
      dtest_type <- "train"
    }
    expect_identical(vs$d_test$type, dtest_type, info = info_str)
    expect_null(vs$d_test$data, info = info_str)
    expect_identical(vs$d_test$offset, vs$refmodel$offset, info = info_str)
    expect_identical(vs$d_test$weights, vs$refmodel$wobs, info = info_str)
    expect_identical(vs$d_test$y, vs$refmodel$y, info = info_str)
    expect_identical(vs$d_test$yOrig, vs$refmodel$yOrig, info = info_str)
  } else {
    expect_identical(vs$d_test, dtest_expected, info = info_str)
  }

  # summaries
  expect_type(vs$summaries, "list")
  expect_named(vs$summaries, c("sub", "ref"), info = info_str)
  expect_type(vs$summaries$sub, "list")
  expect_length(vs$summaries$sub, solterms_len_expected + 1)
  if (with_cv) {
    if (is.null(nloo_expected) || nloo_expected > nobsv) {
      nloo_expected <- nobsv
    }
  }
  nobsv_summ <- nobsv
  nobsv_summ_aug <- nobsv_aug
  if (!is.null(dtest_expected)) {
    nobsv_summ <- nrow(dtest_expected$data)
    nobsv_summ_aug <- nobsv_summ * ncats
  }
  if (vs$refmodel$family$for_latent) {
    vsel_smmrs_sub_nms <- c(vsel_smmrs_sub_nms, "Orig")
    if ("wcv" %in% vsel_smmrs_sub_nms &&
        identical(cv_method_expected, "kfold")) {
      vsel_smmrs_sub_nms[vsel_smmrs_sub_nms %in% c("wcv", "Orig")] <- c("Orig",
                                                                        "wcv")
    }
    vsel_smmrs_ref_nms <- c(vsel_smmrs_ref_nms, "Orig")
  }
  smmrs_sub_j_tester <- function(smmrs_sub_j, tests_Orig = FALSE) {
    if (tests_Orig) {
      vsel_smmrs_sub_nms <- setdiff(vsel_smmrs_sub_nms, "Orig")
      if (!is.null(vs$refmodel$family$cats)) {
        ncats <- length(vs$refmodel$family$cats)
      } else {
        ncats <- 1L
      }
      nobsv_summ_aug <- nobsv_summ * ncats
    }
    expect_type(smmrs_sub_j, "list")
    expect_named(smmrs_sub_j, vsel_smmrs_sub_nms, info = info_str)
    expect_type(smmrs_sub_j$mu, "double")
    expect_length(smmrs_sub_j$mu, nobsv_summ_aug)
    if (vs$refmodel$family$for_augdat ||
        (vs$refmodel$family$for_latent && tests_Orig &&
         !is.null(vs$refmodel$family$cats))) {
      expect_s3_class(smmrs_sub_j$mu, "augvec")
    }
    if (with_cv) {
      expect_identical(sum(!is.na(smmrs_sub_j$mu)), nloo_expected * ncats,
                       info = info_str)
    } else {
      expect_true(all(!is.na(smmrs_sub_j$mu)), info = info_str)
    }
    expect_type(smmrs_sub_j$lppd, "double")
    expect_length(smmrs_sub_j$lppd, nobsv_summ)
    if (with_cv) {
      expect_identical(sum(!is.na(smmrs_sub_j$lppd)), nloo_expected,
                       info = info_str)
    } else {
      expect_true(all(!is.na(smmrs_sub_j$lppd)), info = info_str)
    }
    if (with_cv) {
      expect_type(smmrs_sub_j$wcv, "double")
      expect_length(smmrs_sub_j$wcv, nobsv)
      expect_true(all(!is.na(smmrs_sub_j$wcv)), info = info_str)
      if (nloo_expected == nobsv) {
        expect_equal(smmrs_sub_j$wcv, rep(1 / nobsv, nobsv), info = info_str)
      } else {
        expect_true(any(smmrs_sub_j$wcv != rep(1 / nobsv, nobsv)),
                    info = info_str)
      }
    }
    return(invisible(TRUE))
  }
  for (j in seq_along(vs$summaries$sub)) {
    smmrs_sub_j_tester(vs$summaries$sub[[j]])
    if (vs$refmodel$family$for_latent) {
      smmrs_sub_j_tester(vs$summaries$sub[[j]]$Orig, tests_Orig = TRUE)
    }
  }
  smmrs_ref_tester <- function(smmrs_ref, tests_Orig = FALSE) {
    if (tests_Orig) {
      vsel_smmrs_ref_nms <- setdiff(vsel_smmrs_ref_nms, "Orig")
      if (!is.null(vs$refmodel$family$cats)) {
        ncats <- length(vs$refmodel$family$cats)
      } else {
        ncats <- 1L
      }
      nobsv_summ_aug <- nobsv_summ * ncats
    }
    expect_type(smmrs_ref, "list")
    expect_named(smmrs_ref, vsel_smmrs_ref_nms, info = info_str)
    if (!from_datafit) {
      expect_type(smmrs_ref$mu, "double")
    }
    expect_length(smmrs_ref$mu, nobsv_summ_aug)
    if (vs$refmodel$family$for_augdat ||
        (vs$refmodel$family$for_latent && tests_Orig &&
         !is.null(vs$refmodel$family$cats))) {
      expect_s3_class(smmrs_ref$mu, "augvec")
    }
    if (!from_datafit) {
      expect_true(all(!is.na(smmrs_ref$mu)), info = info_str)
    } else {
      expect_true(all(is.na(smmrs_ref$mu)), info = info_str)
    }
    if (!from_datafit) {
      expect_type(smmrs_ref$lppd, "double")
    }
    expect_length(smmrs_ref$lppd, nobsv_summ)
    if (!from_datafit) {
      expect_true(all(!is.na(smmrs_ref$lppd)), info = info_str)
    } else {
      expect_true(all(is.na(smmrs_ref$lppd)), info = info_str)
    }
    return(invisible(TRUE))
  }
  smmrs_ref_tester(vs$summaries$ref)
  if (vs$refmodel$family$for_latent) {
    smmrs_ref_tester(vs$summaries$ref$Orig, tests_Orig = TRUE)
  }

  # solution_terms
  expect_type(vs$solution_terms, "character")
  expect_length(vs$solution_terms, solterms_len_expected)
  soltrms <- vs$solution_terms
  for (soltrms_plus in grep("\\+", soltrms, value = TRUE)) {
    soltrms <- setdiff(soltrms, soltrms_plus)
    soltrms <- c(soltrms, labels(terms(as.formula(paste(". ~", soltrms_plus)))))
  }
  expect_true(
    all(soltrms %in% split_formula(vs$refmodel$formula,
                                   add_main_effects = FALSE)),
    info = info_str
  )

  # kl
  expect_type(vs$kl, "double")
  expect_length(vs$kl, solterms_len_expected + 1)
  expect_true(all(vs$kl >= 0), info = info_str)
  # Expected to be non-increasing for increasing model size:
  expect_true(all(tail(vs$kl, -1) <= extra_tol * head(vs$kl, -1)),
              info = info_str)
  ### Too unsafe because `length(vs$kl)` is usually small:
  # prop_as_expected <- 0.8
  # expect_true(
  #   mean(tail(vs$kl, -1) <= extra_tol * head(vs$kl, -1)) >=
  #     prop_as_expected,
  #   info = info_str
  # )
  ###

  # pct_solution_terms_cv
  if (with_cv) {
    expect_true(is.matrix(vs$pct_solution_terms_cv), info = info_str)
    expect_type(vs$pct_solution_terms_cv, "double")
    expect_identical(dim(vs$pct_solution_terms_cv),
                     c(solterms_len_expected, 1L + solterms_len_expected),
                     info = info_str)
    expect_identical(colnames(vs$pct_solution_terms_cv),
                     c("size", vs$solution_terms),
                     info = info_str)
    expect_identical(vs$pct_solution_terms_cv[, "size"],
                     as.numeric(seq_len(solterms_len_expected)),
                     info = info_str)
    pct_nonsize_nms <- setdiff(colnames(vs$pct_solution_terms_cv), "size")
    pct_solterms <- vs$pct_solution_terms_cv[, pct_nonsize_nms, drop = FALSE]
    expect_false(anyNA(pct_solterms), info = info_str)
    expect_true(all(pct_solterms >= 0 & pct_solterms <= 1), info = info_str)
    if (isFALSE(vs$validate_search) &&
        !identical(cv_method_expected, "kfold")) {
      expect_true(all(pct_solterms %in% c(0, 1)), info = info_str)
      # More specifically:
      pct_solterms_ch <- matrix(0, nrow = nrow(pct_solterms),
                                ncol = ncol(pct_solterms))
      diag(pct_solterms_ch) <- 1
      colnames(pct_solterms_ch) <- pct_nonsize_nms
      expect_identical(pct_solterms_ch, pct_solterms, info = info_str)
    }
  }

  # nterms_max
  nterms_max_expected <- solterms_len_expected + 1
  if (search_trms_empty_size) {
    nterms_max_expected <- nterms_max_expected + 1
  }
  expect_equal(vs$nterms_max, nterms_max_expected, info = info_str)

  # nterms_all
  expect_identical(vs$nterms_all, count_terms_in_formula(vs$refmodel$formula),
                   info = info_str)

  # method
  expect_identical(vs$method, method_expected, info = info_str)

  # cv_method
  expect_identical(vs$cv_method, cv_method_expected, info = info_str)

  # validate_search
  expect_identical(vs$validate_search, valsearch_expected, info = info_str)

  # clust_used_search
  expect_equal(vs$clust_used_search, cl_search_expected, info = info_str)

  # clust_used_eval
  expect_equal(vs$clust_used_eval, cl_eval_expected, info = info_str)

  # nprjdraws_search
  expect_equal(vs$nprjdraws_search, nprjdraws_search_expected, info = info_str)

  # nprjdraws_eval
  expect_equal(vs$nprjdraws_eval, nprjdraws_eval_expected, info = info_str)

  return(invisible(TRUE))
}

# A helper function for testing the structure of an object as returned by
# summary.vsel().
#
# @param smmry An object as returned by summary.vsel().
# @param vsel_expected The `"vsel"` object which was used in the summary.vsel()
#   call.
# @param nterms_max_expected A single numeric value as supplied to
#   summary.vsel()'s argument `nterms_max`.
# @param respOrig_expected A single logical value indicating whether
#   element `respOrig` is expected to be `TRUE` or `FALSE`.
# @param search_trms_empty_size A single logical value indicating whether
#   `search_terms` was constructed in a way that causes a model size to be
#   without candidate models.
# @param info_str A single character string giving information to be printed in
#   case of failure.
# @param ... Arguments passed to smmry_sel_tester(), apart from
#   smmry_sel_tester()'s arguments `smmry_sel`, `nterms_max_expected`, and
#   `info_str`.
#
# @return `TRUE` (invisible).
smmry_tester <- function(smmry, vsel_expected, nterms_max_expected = NULL,
                         respOrig_expected = TRUE,
                         search_trms_empty_size = FALSE, info_str, ...) {
  expect_s3_class(smmry, "vselsummary")
  expect_type(smmry, "list")
  pct_solterms_nm <- if ("pct_solution_terms_cv" %in% names(vsel_expected)) {
    "pct_solution_terms_cv"
  } else {
    character()
  }
  expect_named(
    smmry,
    c("formula", "family", "nobs_train", "nobs_test", "method", "cv_method",
      "validate_search", "clust_used_search", "clust_used_eval",
      "nprjdraws_search", "nprjdraws_eval", "search_included", "nterms",
      pct_solterms_nm, "selection", "respOrig"),
    info = info_str
  )

  for (nm in c(
    "method", "cv_method", "validate_search", "clust_used_search",
    "clust_used_eval", "nprjdraws_search", "nprjdraws_eval", pct_solterms_nm
  )) {
    expect_identical(smmry[[nm]], vsel_expected[[nm]],
                     info = paste(info_str, nm, sep = "__"))
  }
  expect_identical(smmry$family, vsel_expected$refmodel$family,
                   info = info_str)
  expect_identical(smmry$formula, vsel_expected$refmodel$formula,
                   info = info_str)
  expect_null(smmry$fit, info = info_str)
  expect_identical(smmry$nobs_train, length(vsel_expected$refmodel$y),
                   info = info_str)
  expect_identical(smmry$nobs_test, nrow(vsel_expected$d_test$data),
                   info = info_str)
  # In summary.vsel(), `nterms_max` and output element `nterms` do not count the
  # intercept (whereas `vsel_expected$nterms_max` does):
  if (is.null(nterms_max_expected)) {
    nterms_ch <- vsel_expected$nterms_max - 1
  } else {
    nterms_ch <- nterms_max_expected
  }
  if (search_trms_empty_size) {
    # This is the "empty_size" setting, so we have to subtract the skipped model
    # size (see issue #307):
    nterms_ch <- nterms_ch - 1
  }
  expect_identical(smmry$nterms, nterms_ch,
                   info = info_str)
  expect_true(smmry$search_included %in% c("search included",
                                           "search not included"),
              info = info_str)
  expect_identical(smmry$search_included == "search included",
                   isTRUE(vsel_expected$validate_search),
                   info = info_str)
  smmry_sel_tester(smmry$selection,
                   summaries_ref = vsel_expected$summaries$ref,
                   nterms_max_expected = nterms_max_expected,
                   info_str = info_str, ...)
  expect_identical(smmry$respOrig, respOrig_expected, info = info_str)

  return(invisible(TRUE))
}

# A helper function for testing the structure of a `data.frame` as returned by
# summary.vsel() in its output element `selection`.
#
# @param smmry_sel A `data.frame` as returned by summary.vsel() in its output
#   element `selection`.
# @param summaries_ref The reference model's summaries, as stored in
#   `<vsel_object>$summaries$ref`.
# @param stats_expected A character vector of expected `stats` (see the
#   corresponding argument of summary.vsel()). Use `NULL` for the default.
# @param type_expected A character vector of expected `type`s (see the
#   corresponding argument of summary.vsel()). Use `NULL` for the default.
# @param nterms_max_expected A single numeric value as supplied to
#   summary.vsel()'s argument `nterms_max`.
# @param cv_method_expected Either `character()` for the no-CV case or a single
#   character string giving the CV method (see argument `cv_method` of
#   cv_varsel()).
# @param solterms_expected A character vector giving the expected solution terms
#   (not counting the intercept, even for the intercept-only model).
# @param from_datafit A single logical value indicating whether an object of
#   class `"datafit"` was used for creating the `"vsel"` object (from which
#   `smmry_sel` was created) (`TRUE`) or not (`FALSE`).
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
smmry_sel_tester <- function(
    smmry_sel,
    summaries_ref,
    stats_expected = NULL,
    type_expected = NULL,
    nterms_max_expected = NULL,
    cv_method_expected = character(),
    solterms_expected,
    from_datafit = FALSE,
    info_str
) {
  if (is.null(stats_expected)) {
    stats_expected <- "elpd"
  }
  if (is.null(type_expected)) {
    type_expected <- c("mean", "se", "diff", "diff.se")
  }
  if (is.null(nterms_max_expected)) {
    nterms_max_expected <- length(solterms_expected)
  } else {
    solterms_expected <- head(solterms_expected, nterms_max_expected)
  }

  expect_s3_class(smmry_sel, "data.frame")

  # Rows:
  expect_identical(nrow(smmry_sel), nterms_max_expected + 1L,
                   info = info_str)

  # Columns:
  smmry_nms <- c("size", "solution_terms")
  ### Requires R >= 4.0.1:
  # stats_mean_name <- paste0(
  #   stats_expected,
  #   paste0(".", tolower(cv_method_expected), recycle0 = TRUE)
  # )
  ###
  ### Without relying on R >= 4.0.1:
  if (length(cv_method_expected) == 0) {
    stats_mean_name <- stats_expected
  } else {
    stopifnot(length(cv_method_expected) == 1)
    stats_mean_name <- paste(stats_expected, tolower(cv_method_expected),
                             sep = ".")
  }
  ###
  if (length(stats_expected) == 1) {
    smmry_nms <- c(smmry_nms, stats_mean_name, setdiff(type_expected, "mean"))
  } else {
    smmry_nms <- c(
      smmry_nms,
      sapply(seq_along(stats_expected), function(stat_idx) {
        c(stats_mean_name[stat_idx],
          paste(stats_expected[stat_idx], setdiff(type_expected, "mean"),
                sep = "."))
      })
    )
  }
  expect_named(smmry_sel, smmry_nms, info = info_str)

  # Columns in detail:
  expect_identical(smmry_sel$size, seq_len(nrow(smmry_sel)) - 1,
                   info = info_str)
  expect_identical(smmry_sel$solution_terms,
                   c(NA_character_, solterms_expected),
                   info = info_str)
  if ("diff" %in% type_expected) {
    if (length(stats_expected) == 1) {
      diff_nm <- "diff"
    } else {
      diff_nm <- paste(stats_expected, "diff", sep = ".")
    }
    for (stat_idx in seq_along(stats_expected)) {
      if (!from_datafit) {
        expect_equal(
          diff(smmry_sel[, stats_mean_name[stat_idx]]),
          diff(smmry_sel[, diff_nm[stat_idx]]),
          info = info_str
        )
        if (stats_expected[stat_idx] == "elpd") {
          stat_ref <- sum(summaries_ref$lppd)
        } else if (stats_expected[stat_idx] == "mlpd") {
          stat_ref <- mean(summaries_ref$lppd)
        } else {
          # TODO: Implement `stat_ref` for the remaining `stats`.
          stat_ref <- NULL
        }
        if (!is.null(stat_ref)) {
          expect_equal(
            smmry_sel[, stats_mean_name[stat_idx]] - stat_ref,
            smmry_sel[, diff_nm[stat_idx]],
            tolerance = 1e-12, info = info_str
          )
        }
      } else {
        expect_true(all(is.na(smmry_sel[, diff_nm[stat_idx]])), info = info_str)
      }
    }
  }
  if ("lower" %in% type_expected) {
    if (length(stats_expected) == 1) {
      lower_nm <- "lower"
    } else {
      lower_nm <- paste(stats_expected, "lower", sep = ".")
    }
    for (stat_idx in seq_along(stats_expected)) {
      if (!stats_expected[stat_idx] %in% c("rmse", "auc")) {
        # RMSE and AUC are excluded here because of PR #347.
        expect_true(all(smmry_sel[, stats_mean_name[stat_idx]] >=
                          smmry_sel[, lower_nm[stat_idx]]),
                    info = info_str)
      }
    }
  }
  if ("upper" %in% type_expected) {
    if (length(stats_expected) == 1) {
      upper_nm <- "upper"
    } else {
      upper_nm <- paste(stats_expected, "upper", sep = ".")
    }
    for (stat_idx in seq_along(stats_expected)) {
      if (!stats_expected[stat_idx] %in% c("rmse", "auc")) {
        # RMSE and AUC are excluded here because of PR #347.
        expect_true(all(smmry_sel[, stats_mean_name[stat_idx]] <=
                          smmry_sel[, upper_nm[stat_idx]]),
                    info = info_str)
      }
    }
  }

  return(invisible(TRUE))
}
