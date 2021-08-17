# A helper function for testing the structure of an expected extended `"family"`
# object
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
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
extfam_tester <- function(extfam,
                          fam_orig,
                          extfam_nms_add2 = character(),
                          from_brms = FALSE,
                          info_str) {
  # Some minimal checks for `fam_orig`:
  expect_s3_class(fam_orig, "family")
  expect_type(fam_orig, "list")
  fam_orig_nms <- fam_orig_nms_common <- c(
    "family", "link", "linkfun", "linkinv", "variance", "dev.resids", "aic",
    "mu.eta", "initialize", "validmu", "valideta"
  )
  expect_true(all(fam_orig_nms_common %in% names(fam_orig)), info = info_str)
  if (fam_orig$family %in% c("binomial", "poisson")) {
    fam_orig_nms <- c(fam_orig_nms, "simulate")
  }
  expect_named(fam_orig, fam_orig_nms, info = info_str)

  # Now the checks for `extfam` (first starting with the general structure):
  extfam_nms_add <- c("kl", "dis_fun", "predvar", "ll_fun", "deviance", "ppd",
                      extfam_nms_add2)
  extfam_nms <- c(names(fam_orig), extfam_nms_add)
  expect_s3_class(extfam, "family")
  expect_type(extfam, "list")
  expect_named(extfam, extfam_nms, info = info_str)

  fam_orig_ch <- structure(extfam[names(fam_orig)], class = "family")
  if (extfam$family == "binomial") {
    fam_orig_ch$initialize <- fam_orig$initialize
  }
  if (!from_brms) {
    expect_identical(fam_orig_ch, fam_orig, info = info_str)
  } else {
    expect_identical(fam_orig_ch, fam_orig, ignore.environment = TRUE,
                     info = info_str)
  }

  for (el_nm in extfam_nms_add) {
    expect_type(extfam[[el_nm]], "closure")
  }

  # TODO: Add some mathematical checks (i.e., check that the calculations for
  # the objects listed in `extfam_nms_add` are mathematically correct).

  return(invisible(TRUE))
}

# A helper function for testing the structure of an expected `"refmodel"` object
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
  if (!is.null(attr(terms(formul_expected), "offset"))) {
    # In the reference model, the offset() term is placed last:
    formul_expected <- update(formul_expected,
                              . ~ . - offset(offs_col) + offset(offs_col))
  }
  formul_expected <- rm_cbind(formul_expected)
  formul_expected <- rm_addresp(formul_expected)
  if (with_spclformul) {
    # Reference models take arithmetic expressions on the left-hand side of
    # the formula into account:
    y_spclformul <- as.character(formul_expected)[2]
    y_spclformul_new <- gsub("\\(|\\)", "", y_spclformul)
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
      refmod$family$family == "binomial") {
    data_expected$temp_y <- 1
  }
  has_grp <- mod_nm %in% c("glmm", "gamm")
  has_add <- mod_nm %in% c("gam", "gamm")
  is_gamm <- mod_nm == "gamm"

  # Test the general structure of the object:
  refmod_nms <- c(
    "fit", "formula", "div_minimizer", "family", "mu", "dis", "y", "loglik",
    "intercept", "proj_predfun", "fetch_data", "wobs", "wsample", "offset",
    "cvfun", "cvfits", "extract_model_data", "ref_predfun"
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
  extfam_tester(refmod$family, fam_orig = fam_orig,
                extfam_nms_add2 = "mu_fun",
                from_brms = (pkg_nm == "brms"),
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
      drws <- as.matrix(refmod$fit)
      drws_icpt <- drws[, "(Intercept)"]
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
        mm_s <- `%:::%`("rstanarm", "pp_data")(refmod$fit)$x
        mm_s <- mm_s[, grep("^s\\(", colnames(mm_s), value = TRUE),
                     drop = FALSE]
        ###
        mu_expected <- mu_expected + drws_beta_s %*% t(mm_s)
      }
      mu_expected <- unname(mu_expected)
    } else if (pkg_nm == "brms") {
      # TODO (brms): Do this manually (but posterior_linpred.brmsfit() should
      # already be tested extensively in brms, so perhaps we don't need a manual
      # calculation here (the reason for the manual calculation for "stanreg"s
      # were the offset issues in rstanarm)):
      mu_expected <- posterior_linpred(refmod$fit) - matrix(
        offs_expected,
        nrow = nrefdraws_expected,
        ncol = nobsv_expected,
        byrow = TRUE
      )
    }
    if (refmod$family$family != "gaussian") {
      mu_expected <- fam_orig$linkinv(mu_expected)
    }
    expect_equal(refmod$mu, t(mu_expected), info = info_str)
  } else {
    if (refmod$family$family != "binomial") {
      expect_identical(refmod$mu, as.matrix(refmod$y), info = info_str)
    } else {
      expect_identical(refmod$mu, as.matrix(refmod$y / refmod$wobs),
                       info = info_str)
    }
  }

  # dis
  if (refmod$family$family == "gaussian") {
    expect_true(is.vector(refmod$dis, "double"), info = info_str)
    expect_length(refmod$dis, nrefdraws_expected)
    if (!is_datafit) {
      expect_true(all(refmod$dis > 0), info = info_str)
    } else {
      expect_identical(refmod$dis, 0, info = info_str)
    }
  } else {
    expect_identical(refmod$dis, rep(0, nrefdraws_expected), info = info_str)
  }

  # y
  ### Not needed because of the more precise test below:
  # expect_true(is.vector(refmod$y, "numeric"), info = info_str)
  # expect_length(refmod$y, nobsv_expected)
  ###
  if (!with_spclformul) {
    y_expected <- dat[[paste("y", mod_nm, fam_nm, sep = "_")]]
    if (!is_datafit && pkg_nm == "brms" && fam_nm == "brnll") {
      y_expected <- as.numeric(y_expected)
    }
  } else {
    y_expected <- data_expected[[y_spclformul_new]]
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
        refdat_ch <- data_expected[, refdat_colnms, drop = FALSE]
        expect_equal(refmod$fetch_data(), refdat_ch, check.attributes = FALSE,
                     info = info_str)
      } else {
        # TODO: The check in this case is not perfect yet.
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

  return(invisible(TRUE))
}

# A helper function for testing the structure of a list of subfits (whose
# elements must not necessarily be of class `"subfit"`) for the same single
# submodel
#
# @param sub_fit_obj The list of subfits to test.
# @param nprjdraws_expected A single numeric value giving the expected number of
#   projected draws.
# @param sub_formul A list of formulas for the submodel (with one element per
#   projected draw).
# @param sub_data The dataset used for fitting the submodel.
# @param sub_fam A single character string giving the submodel's family.
# @param has_grp A single logical value indicating whether `sub_fit_obj` is
#   expected to be of class `"lmerMod"` or `"glmerMod"` (if, at the same time,
#   `has_add` is `FALSE`).
# @param has_add A single logical value indicating whether `sub_fit_obj` is
#   expected to be of class `"gam"` or `"gamm4"` (depending on whether the
#   submodel is non-multilevel or multilevel, respectively).
# @param wobs_expected The expected numeric vector of observation weights.
# @param ref_formul The formula of the reference model . Should only be needed
#   if `sub_fit_obj` comes from the L1 `search_path` of an object of class
#   `"vsel"`. Otherwise, use `NULL`.
# @param with_offs A single logical value indicating whether `sub_fit_obj` is
#   expected to include offsets (`TRUE`) or not (`FALSE`).
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
sub_fit_tester <- function(
  sub_fit_obj,
  nprjdraws_expected,
  sub_formul,
  sub_data,
  sub_fam,
  has_grp = formula_contains_group_terms(sub_formul[[1]]),
  has_add = formula_contains_additive_terms(sub_formul[[1]]),
  wobs_expected = wobs_tst,
  ref_formul = NULL,
  with_offs = FALSE,
  info_str
) {
  if (nprjdraws_expected > 1) {
    expect_type(sub_fit_obj, "list")
    expect_length(sub_fit_obj, nprjdraws_expected)
    sub_fit_totest <- sub_fit_obj
  } else {
    sub_fit_totest <- list(sub_fit_obj)
  }

  from_vsel_L1_search <- !is.null(ref_formul)

  seq_extensive_tests <- unique(round(
    seq(1, length(sub_fit_totest),
        length.out = min(length(sub_fit_totest), nclusters_pred_tst))
  ))

  if (!has_grp && !has_add) {
    sub_trms <- labels(terms(sub_formul[[1]]))
    if (length(sub_trms) > 0) {
      ncoefs <- sum(sapply(sub_trms, function(trm_i) {
        ncol(model.matrix(
          as.formula(paste("~ 0 +", trm_i)),
          data = sub_data
        ))
      }))
      ### As discussed in issue #149, the following might be more appropriate:
      # ncoefs <- ncol(model.matrix(
      #   as.formula(paste("~", paste(sub_trms, collapse = " + "))),
      #   data = sub_data
      # )) - 1L
      ###
    } else {
      ncoefs <- 0L
    }
    if (!from_vsel_L1_search) {
      formul_for_mm <- sub_formul[[1]]
    } else {
      formul_for_mm <- ref_formul
    }
    sub_trms_for_mm <- labels(terms(formul_for_mm))
    if (length(sub_trms_for_mm) > 0) {
      if (any(grepl("xca\\.", sub_trms_for_mm))) {
        sub_contr <- lapply(
          setNames(nm = grep("xca\\.", sub_trms_for_mm, value = TRUE)),
          function(x_nm) {
            contrasts(get(x_nm, envir = as.environment(sub_data)),
                      contrasts = FALSE)
          }
        )
      } else {
        sub_contr <- NULL
      }
    } else {
      sub_contr <- NULL
    }
    sub_x_expected <- model.matrix(update(formul_for_mm, NULL ~ . + 0),
                                   data = sub_data,
                                   contrasts.arg = sub_contr)
    if (from_vsel_L1_search) {
      subfit_nms <- setdiff(subfit_nms, "y")
    }
    for (j in seq_along(sub_fit_totest)) {
      expect_s3_class(sub_fit_totest[[!!j]], "subfit")
      expect_type(sub_fit_totest[[!!j]], "list")
      expect_named(sub_fit_totest[[!!j]], subfit_nms, info = info_str)

      if (j %in% seq_extensive_tests) {
        expect_true(is.vector(sub_fit_totest[[!!j]]$alpha, "double"),
                    info = info_str)
        expect_length(sub_fit_totest[[!!j]]$alpha, 1)

        if (length(sub_trms) > 0 || !from_vsel_L1_search) {
          expect_true(is.matrix(sub_fit_totest[[!!j]]$beta), info = info_str)
          expect_true(is.numeric(sub_fit_totest[[!!j]]$beta), info = info_str)
          expect_identical(dim(sub_fit_totest[[!!j]]$beta), c(ncoefs, 1L),
                           info = info_str)
        } else if (length(sub_trms) == 0) {
          expect_null(sub_fit_totest[[!!j]]$beta, info = info_str)
        }

        if (!from_vsel_L1_search) {
          expect_true(is.matrix(sub_fit_totest[[!!j]]$w), info = info_str)
          expect_type(sub_fit_totest[[!!j]]$w, "double")
          expect_identical(dim(sub_fit_totest[[!!j]]$w), c(nobsv, 1L),
                           info = info_str)
        } else {
          expect_true(is.vector(sub_fit_totest[[!!j]]$w, "double"),
                      info = info_str)
          expect_length(sub_fit_totest[[!!j]]$w, nobsv)
        }
        expect_true(all(sub_fit_totest[[!!j]]$w > 0), info = info_str)

        expect_s3_class(sub_fit_totest[[!!j]]$formula, "formula")
        if (!grepl(":", as.character(sub_fit_totest[[j]]$formula)[3])) {
          expect_equal(sub_fit_totest[[!!j]]$formula, sub_formul[[!!j]],
                       info = info_str)
        } else {
          # The order of interactions might be changed in the reference model:
          expect_equal(sub_fit_totest[[!!j]]$formula[[2]],
                       sub_formul[[!!j]][[2]],
                       info = info_str)
          expect_equal(labels(terms(sub_fit_totest[[!!j]]$formula)),
                       labels(terms(sub_formul[[!!j]])),
                       info = info_str)
        }

        expect_identical(sub_fit_totest[[!!j]]$x, sub_x_expected,
                         info = info_str)

        if (!from_vsel_L1_search) {
          y_ch <- setNames(eval(str2lang(as.character(sub_formul[[j]])[2]),
                                sub_data),
                           seq_len(nobsv))
          expect_identical(sub_fit_totest[[!!j]]$y, y_ch, info = info_str)
        }
      }
    }
  } else if (has_grp && !has_add) {
    if (sub_fam == "gaussian") {
      for (j in seq_along(sub_fit_totest)) {
        expect_s4_class(sub_fit_totest[[!!j]], "lmerMod")
      }
    } else {
      for (j in seq_along(sub_fit_totest)) {
        expect_s4_class(sub_fit_totest[[!!j]], "glmerMod")
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
      expect_equal(sub_fit_totest[[!!j]]@call[["formula"]],
                   sub_formul_expected,
                   info = info_str)

      # resp
      if (!with_offs) {
        offs_expected <- numeric(nobsv)
      } else {
        offs_expected <- offs_tst
      }
      expect_identical(sub_fit_totest[[!!j]]@resp$offset,
                       offs_expected,
                       info = info_str)
      if (!is.null(wobs_expected)) {
        expect_equal(sub_fit_totest[[!!j]]@resp$weights,
                     wobs_expected,
                     info = info_str)
      } else {
        expect_equal(sub_fit_totest[[!!j]]@resp$weights,
                     rep(1, nobsv),
                     info = info_str)
      }
      expect_equal(sub_fit_totest[[!!j]]@resp$y,
                   eval(str2lang(as.character(sub_formul[[!!j]])[2]),
                        sub_data),
                   info = info_str)

      # frame
      expect_identical(sub_fit_totest[[!!j]]@frame,
                       model.frame(sub_fit_totest[[!!j]]),
                       info = info_str)
      expect_equal(
        sub_fit_totest[[!!j]]@frame[[
          grep("y_|ybinprop", names(sub_fit_totest[[!!j]]@frame), value = TRUE)
        ]],
        sub_fit_totest[[!!j]]@resp$y,
        info = info_str
      )
      if (!is.null(wobs_expected)) {
        expect_equal(sub_fit_totest[[!!j]]@frame$`(weights)`,
                     sub_fit_totest[[!!j]]@resp$weights,
                     info = info_str)
      } else {
        expect_null(sub_fit_totest[[!!j]]@frame$`(weights)`,
                    info = info_str)
      }
      if (with_offs) {
        expect_equal(sub_fit_totest[[!!j]]@frame$`offset(offs_col)`,
                     offs_expected,
                     info = info_str)
      }
      frame_nms <- grep("y_|ybinprop|^\\(weights\\)$|^offset\\(.*\\)$",
                        names(sub_fit_totest[[j]]@frame),
                        value = TRUE,
                        invert = TRUE)
      expect_setequal(frame_nms, names(sub_mf_expected))
      expect_equal(
        sub_fit_totest[[!!j]]@frame[frame_nms],
        sub_mf_expected[frame_nms],
        info = info_str
      )

      # model.matrix()
      expect_identical(model.matrix(sub_fit_totest[[!!j]]), mm_expected,
                       info = info_str)

      # flist
      expect_type(sub_fit_totest[[!!j]]@flist, "list")
      expect_length(sub_fit_totest[[!!j]]@flist, length(nlvl_ran))
      z_nms <- intersect(names(sub_fit_totest[[j]]@flist),
                         names(sub_fit_totest[[j]]@frame))
      expect_identical(sub_fit_totest[[!!j]]@flist[z_nms],
                       as.list(sub_fit_totest[[!!j]]@frame[z_nms]),
                       info = info_str)

      # coef()
      coefs_crr <- coef(sub_fit_totest[[!!j]])
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
    for (j in seq_along(sub_fit_totest)) {
      expect_s3_class(sub_fit_totest[[!!j]], "gam")
    }
    # TODO (GAMs): Add more expectations for GAMs.
  } else if (has_grp && has_add) {
    for (j in seq_along(sub_fit_totest)) {
      expect_s3_class(sub_fit_totest[[!!j]], "gamm4")
    }
    # TODO (GAMMs): Add more expectations for GAMMs.
  }

  return(invisible(TRUE))
}

# A helper function for testing the structure of an expected `"projection"`
# object
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
# @param fam_expected The expected `"family"` object or `NULL` for not testing
#   the family object at all.
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
                              fam_expected = NULL,
                              prjdraw_weights_expected = NULL,
                              from_vsel_L1_search = FALSE,
                              info_str = "") {
  expect_s3_class(p, "projection")
  expect_type(p, "list")
  # Check the names using `ignore.order = FALSE` because an incorrect
  # order would mean that the documentation of project()'s return value
  # would have to be updated:
  expect_named(p, projection_nms, info = info_str)

  # refmodel
  # Note: Extensive tests for `"refmodel"`s and `"datafit"`s may be run via
  # refmodel_tester().
  expect_identical(p$refmodel, refmod_expected, info = info_str)

  # extract_model_data
  expect_identical(p$extract_model_data, p$refmodel$extract_model_data,
                   info = info_str)

  # intercept
  expect_identical(p$intercept, p$refmodel$intercept, info = info_str)

  # family
  expect_identical(p$family, p$refmodel$family, info = info_str)
  if (!is.null(fam_expected)) {
    expect_identical(p$family, fam_expected, info = info_str)
  }

  # solution_terms
  if (is.numeric(solterms_expected)) {
    expect_length(p$solution_terms, solterms_expected)
    # Same check, but using count_terms_chosen():
    expect_equal(count_terms_chosen(p$solution_terms, add_icpt = TRUE),
                 solterms_expected + 1, info = info_str)
  } else if (is.character(solterms_expected)) {
    expect_identical(p$solution_terms, solterms_expected, info = info_str)
  }

  # sub_fit
  sub_trms_crr <- p$solution_terms
  if (length(sub_trms_crr) == 0) {
    sub_trms_crr <- as.character(as.numeric(p$intercept))
  }
  if (!from_vsel_L1_search) {
    y_nm <- as.character(p$refmodel$formula)[2]
    ref_formul_crr <- NULL
  } else {
    y_nm <- ""
    ref_formul_crr <- p$refmodel$formula
  }
  y_nms <- paste0(".", y_nm)
  # A preliminary check for `nprjdraws_expected` (doesn't work for "datafit"s
  # and, because of issue #131, for submodels with multilevel terms which belong
  # to an additive reference model):
  sub_formul_crr_rhs <- as.formula(paste(
    "~", paste(sub_trms_crr, collapse = " + ")
  ))
  if (!inherits(p$refmodel, "datafit") &&
      !(formula_contains_additive_terms(p$refmodel$formula) &&
        formula_contains_group_terms(sub_formul_crr_rhs))) {
    # Number of projected draws in as.matrix.projection() (note that more
    # extensive tests for as.matrix.projection() may be found in
    # "test_as_matrix.R"):
    expect_identical(suppressWarnings(NROW(as.matrix(p))), nprjdraws_expected,
                     info = info_str)
  }
  if (nprjdraws_expected > 1) {
    y_nms <- paste0(y_nms, ".", seq_len(nprjdraws_expected))
  }
  sub_formul_crr <- lapply(y_nms, function(y_nm_i) {
    as.formula(paste(
      y_nm_i, "~", paste(sub_trms_crr, collapse = " + ")
    ))
  })
  sub_data_crr <- p$refmodel$fetch_data()
  if (p_type_expected) {
    clust_ref <- .get_refdist(p$refmodel,
                              nclusters = nprjdraws_expected,
                              seed = seed_expected)
  } else {
    clust_ref <- .get_refdist(p$refmodel,
                              ndraws = nprjdraws_expected,
                              seed = seed_expected)
  }
  for (i in seq_len(nprjdraws_expected)) {
    sub_data_crr[[y_nms[i]]] <- clust_ref$mu[, i]
  }
  sub_fit_tester(p$sub_fit,
                 nprjdraws_expected = nprjdraws_expected,
                 sub_formul = sub_formul_crr,
                 sub_data = sub_data_crr,
                 sub_fam = p$family$family,
                 has_add = formula_contains_additive_terms(p$refmodel$formula),
                 wobs_expected = p$refmodel$wobs,
                 ref_formul = ref_formul_crr,
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

  # p_type
  expect_identical(p$p_type, p_type_expected, info = info_str)

  return(invisible(TRUE))
}

# A helper function for testing the structure of an expected `"proj_list"`
# object
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
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
pl_tester <- function(pl,
                      len_expected = 1,
                      nprjdraws_expected = nclusters_pred_tst,
                      nobsv_expected = nobsv,
                      lpd_null_expected = FALSE,
                      info_str) {
  if (len_expected == 1) {
    pl <- list(pl)
  } else {
    expect_type(pl, "list")
    expect_length(pl, len_expected)
  }
  for (j in seq_along(pl)) {
    expect_named(pl[[!!j]], c("pred", "lpd"), info = info_str)
    expect_identical(dim(pl[[!!j]]$pred),
                     c(nprjdraws_expected, nobsv_expected),
                     info = info_str)
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
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
pp_tester <- function(pp,
                      len_expected = 1,
                      nprjdraws_out_expected = nresample_clusters_default,
                      nobsv_expected = nobsv,
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
  }
  return(invisible(TRUE))
}

# A helper function for testing the structure of an expected `"vsel"` object
#
# @param vs An object of class `"vsel"` (at least expected so).
# @param with_cv A single logical value indicating whether `vs` was created by
#   cv_varsel() (`TRUE`) or not (`FALSE`).
# @param refmod_expected The expected `vs$refmodel` object.
# @param dtest_expected The expected `vs$d_test` object.
# @param solterms_len_expected A single numeric value giving the expected number
#   of solution terms (not counting the intercept, even for the intercept-only
#   model).
# @param method_expected The expected `vs$method` object.
# @param cv_method_expected The expected `vs$cv_method` object.
# @param valsearch_expected The expected `vs$validate_search` object.
# @param ndraws_expected The expected `vs$ndraws` object.
# @param ndraws_pred_expected The expected `vs$ndraws_pred` object.
# @param nclusters_expected The expected `vs$nclusters` object (not adapted for
#   L1 search).
# @param nclusters_pred_expected The expected `vs$nclusters_pred` object.
# @param seed_expected The seed which was used for clustering the posterior
#   draws of the reference model.
# @param nloo_expected Only relevant if `with_cv` is `TRUE`. The value which was
#   used for argument `nloo` of cv_varsel().
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
  ndraws_expected = if (!from_datafit) ndraws_default else 1L,
  ndraws_pred_expected = if (!from_datafit) ndraws_pred_default else 1L,
  nclusters_expected = NULL,
  nclusters_pred_expected = NULL,
  seed_expected = seed_tst,
  nloo_expected = NULL,
  extra_tol = 1.1,
  info_str = ""
) {
  # Preparations:
  dtest_type <- "train"
  if (with_cv) {
    vsel_nms <- vsel_nms_cv
    vsel_smmrs_sub_nms <- c("lppd", "mu", "w")
    vsel_smmrs_ref_nms <- c("lppd", "mu")

    if (is.null(cv_method_expected)) {
      cv_method_expected <- "LOO"
    }
    if (is.null(valsearch_expected)) {
      valsearch_expected <- TRUE
    }

    dtest_type <- cv_method_expected
    if (cv_method_expected == "LOO") {
      # Re-order:
      dtest_nms <- dtest_nms[c(1, 5, 2, 4, 3, 6)]
    } else if (cv_method_expected == "kfold") {
      # Re-order and remove `"data"`:
      dtest_nms <- dtest_nms[c(1, 4, 2, 6, 5)]
      # Re-order:
      vsel_smmrs_sub_nms[1:2] <- vsel_smmrs_sub_nms[2:1]
      vsel_smmrs_ref_nms[1:2] <- vsel_smmrs_ref_nms[2:1]
    }
  }
  method_expected <- tolower(method_expected)
  if (method_expected == "l1") {
    nclusters_expected <- 1
  }

  # Test the general structure of the object:
  expect_s3_class(vs, "vsel")
  expect_type(vs, "list")
  expect_named(vs, vsel_nms, info = info_str)

  # refmodel
  expect_identical(vs$refmodel, refmod_expected, info = info_str)

  # search_path
  expect_type(vs$search_path, "list")
  expect_named(vs$search_path, searchpth_nms, info = info_str)
  expect_identical(vs$search_path$solution_terms, vs$solution_terms,
                   info = info_str)
  expect_type(vs$search_path$sub_fits, "list")
  expect_length(vs$search_path$sub_fits, solterms_len_expected + 1)
  from_vsel_L1_search <- method_expected == "l1"
  clust_ref <- .get_refdist(vs$refmodel,
                            ndraws = ndraws_expected,
                            nclusters = nclusters_expected,
                            seed = seed_expected)
  nprjdraws_expected <- ncol(clust_ref$mu)
  if (!from_vsel_L1_search) {
    y_nm <- as.character(vs$refmodel$formula)[2]
    ref_formul_crr <- NULL
  } else {
    y_nm <- ""
    ref_formul_crr <- vs$refmodel$formula
  }
  y_nms <- paste0(".", y_nm)
  if (nprjdraws_expected > 1) {
    y_nms <- paste0(y_nms, ".", seq_len(nprjdraws_expected))
  }
  sub_data_crr <- vs$refmodel$fetch_data()
  for (i in seq_len(nprjdraws_expected)) {
    sub_data_crr[[y_nms[i]]] <- clust_ref$mu[, i]
  }
  solterms_for_subfits <- c(as.character(as.numeric(vs$refmodel$intercept)),
                            vs$solution_terms)
  for (i in seq_along(vs$search_path$sub_fits)) {
    sub_trms_crr <- head(solterms_for_subfits, i)
    if (length(sub_trms_crr) > 1) {
      sub_trms_crr <- setdiff(sub_trms_crr, "1")
    }
    sub_formul_crr <- lapply(y_nms, function(y_nm_i) {
      as.formula(paste(
        y_nm_i, "~", paste(sub_trms_crr, collapse = " + ")
      ))
    })
    sub_fit_tester(
      vs$search_path$sub_fits[[i]],
      nprjdraws_expected = nprjdraws_expected,
      sub_formul = sub_formul_crr,
      sub_data = sub_data_crr,
      sub_fam = vs$family$family,
      has_add = formula_contains_additive_terms(vs$refmodel$formula),
      wobs_expected = vs$refmodel$wobs,
      ref_formul = ref_formul_crr,
      info_str = paste(info_str, i, sep = "__")
    )
  }
  expect_type(vs$search_path$p_sel, "list")
  expect_named(vs$search_path$p_sel, psel_nms, info = info_str)
  expect_true(is.matrix(vs$search_path$p_sel$mu), info = info_str)
  expect_type(vs$search_path$p_sel$mu, "double")
  expect_equal(dim(vs$search_path$p_sel$mu), c(nobsv, nclusters_expected),
               info = info_str)
  if (vs$family$family == "gaussian") {
    expect_true(is.matrix(vs$search_path$p_sel$var), info = info_str)
    expect_type(vs$search_path$p_sel$var, "double")
    expect_equal(dim(vs$search_path$p_sel$var), c(nobsv, nclusters_expected),
                 info = info_str)
  } else {
    expect_type(vs$search_path$p_sel$var, "double")
    expect_length(vs$search_path$p_sel$var, nclusters_expected)
  }
  expect_type(vs$search_path$p_sel$weights, "double")
  expect_length(vs$search_path$p_sel$weights, nclusters_expected)
  expect_true(is.numeric(vs$search_path$p_sel$cl), info = info_str)
  expect_length(vs$search_path$p_sel$cl, ncol(vs$refmodel$mu))

  # d_test
  if (is.null(dtest_expected)) {
    expect_type(vs$d_test, "list")
    expect_named(vs$d_test, dtest_nms, info = info_str)
    if (identical(cv_method_expected, "kfold")) {
      expect_identical(vs$d_test$y[order(vs$d_test$test_points)],
                       vs$refmodel$y, info = info_str)
      expect_identical(vs$d_test$test_points[order(vs$d_test$test_points)],
                       seq_len(nobsv), info = info_str)
      expect_identical(vs$d_test$weights[order(vs$d_test$test_points)],
                       vs$refmodel$wobs, info = info_str)
    } else {
      expect_identical(vs$d_test$y, vs$refmodel$y, info = info_str)
      expect_identical(vs$d_test$test_points, seq_len(nobsv), info = info_str)
      expect_identical(vs$d_test$weights, vs$refmodel$wobs, info = info_str)
    }
    expect_null(vs$d_test$data, info = info_str)
    expect_identical(vs$d_test$type, dtest_type, info = info_str)
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
  for (j in seq_along(vs$summaries$sub)) {
    expect_named(vs$summaries$sub[[!!j]], vsel_smmrs_sub_nms, info = info_str)
    expect_type(vs$summaries$sub[[!!j]]$mu, "double")
    expect_length(vs$summaries$sub[[!!j]]$mu, nobsv)
    if (with_cv) {
      expect_identical(sum(!is.na(vs$summaries$sub[[!!j]]$mu)),
                       nloo_expected, info = info_str)
    } else {
      expect_true(all(!is.na(vs$summaries$sub[[!!j]]$mu)), info = info_str)
    }
    expect_type(vs$summaries$sub[[!!j]]$lppd, "double")
    expect_length(vs$summaries$sub[[!!j]]$lppd, nobsv)
    if (with_cv) {
      expect_identical(sum(!is.na(vs$summaries$sub[[!!j]]$lppd)),
                       nloo_expected, info = info_str)
    } else {
      expect_true(all(!is.na(vs$summaries$sub[[!!j]]$lppd)), info = info_str)
    }
    if (with_cv) {
      expect_type(vs$summaries$sub[[!!j]]$w, "double")
      expect_length(vs$summaries$sub[[!!j]]$w, nobsv)
      expect_true(all(!is.na(vs$summaries$sub[[!!j]]$w)), info = info_str)
      if (nloo_expected == nobsv) {
        expect_equal(vs$summaries$sub[[!!j]]$w, rep(1 / nobsv, nobsv),
                     info = info_str)
      } else {
        expect_true(any(vs$summaries$sub[[!!j]]$w != rep(1 / nobsv, nobsv)),
                    info = info_str)
      }
    }
  }
  expect_type(vs$summaries$ref, "list")
  expect_named(vs$summaries$ref, vsel_smmrs_ref_nms, info = info_str)
  expect_length(vs$summaries$ref$mu, nobsv)
  if (!from_datafit) {
    expect_true(all(!is.na(vs$summaries$ref$mu)), info = info_str)
  } else {
    expect_true(all(is.na(vs$summaries$ref$mu)), info = info_str)
  }
  expect_length(vs$summaries$ref$lppd, nobsv)
  if (!from_datafit) {
    expect_true(all(!is.na(vs$summaries$ref$lppd)), info = info_str)
  } else {
    expect_true(all(is.na(vs$summaries$ref$lppd)), info = info_str)
  }

  # family
  expect_s3_class(vs$family, "family")
  expect_identical(vs$family, refmod_expected$family, info = info_str)

  # solution_terms
  expect_type(vs$solution_terms, "character")
  expect_length(vs$solution_terms, solterms_len_expected)
  expect_true(
    all(vs$solution_terms %in% split_formula(vs$refmodel$formula,
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
  expect_identical(vs$nterms_max, solterms_len_expected + 1, info = info_str)

  # nterms_all
  expect_identical(vs$nterms_all, count_terms_in_formula(vs$refmodel$formula),
                   info = info_str)

  # method
  expect_identical(vs$method, method_expected, info = info_str)

  # cv_method
  expect_identical(vs$cv_method, cv_method_expected, info = info_str)

  # validate_search
  expect_identical(vs$validate_search, valsearch_expected, info = info_str)

  # ndraws
  expect_equal(vs$ndraws, ndraws_expected, info = info_str)

  # ndraws_pred
  expect_equal(vs$ndraws_pred, ndraws_pred_expected, info = info_str)

  # nclusters
  expect_equal(vs$nclusters, nclusters_expected, info = info_str)

  # nclusters_pred
  expect_equal(vs$nclusters_pred, nclusters_pred_expected, info = info_str)

  # suggested_size
  expect_type(vs$suggested_size, "double")
  expect_length(vs$suggested_size, 1)

  # summary
  smmry_sel_tester(
    vs$summary,
    summaries_ref = vs$summaries$ref,
    cv_method_expected = if (with_cv) cv_method_expected else character(),
    solterms_expected = vs$solution_terms,
    from_datafit = from_datafit,
    info_str = info_str
  )

  return(invisible(TRUE))
}

# A helper function for testing the structure of an object as returned by
# summary.vsel()
#
# @param smmry An object as returned by summary.vsel().
# @param vsel_expected The `"vsel"` object which was used in the summary.vsel()
#   call.
# @param nterms_max_expected A single numeric value as supplied to
#   summary.vsel()'s argument `nterms_max`.
# @param info_str A single character string giving information to be printed in
#   case of failure.
# @param ... Arguments passed to smmry_sel_tester(), apart from
#   smmry_sel_tester()'s arguments `smmry_sel`, `nterms_max_expected`, and
#   `info_str`.
#
# @return `TRUE` (invisible).
smmry_tester <- function(smmry, vsel_expected, nterms_max_expected = NULL,
                         info_str, ...) {
  expect_s3_class(smmry, "vselsummary")
  expect_type(smmry, "list")
  pct_solterms_nm <- if ("pct_solution_terms_cv" %in% names(vsel_expected)) {
    "pct_solution_terms_cv"
  } else {
    character()
  }
  expect_named(
    smmry,
    c("formula", "fit", "family", "nobs", "method", "cv_method",
      "validate_search", "ndraws", "ndraws_pred", "nclusters", "nclusters_pred",
      "search_included", "nterms", pct_solterms_nm, "suggested_size",
      "selection"),
    info = info_str
  )

  for (nm in c(
    "family", "method", "cv_method", "validate_search", "ndraws", "ndraws_pred",
    "nclusters", "nclusters_pred", pct_solterms_nm, "suggested_size"
  )) {
    expect_identical(smmry[[nm]], vsel_expected[[nm]],
                     info = paste(info_str, nm, sep = "__"))
  }
  expect_identical(smmry$formula, vsel_expected$refmodel$formula,
                   info = info_str)
  expect_null(smmry$fit, info = info_str)
  expect_identical(smmry$nobs, length(vsel_expected$refmodel$y),
                   info = info_str)
  # In summary.vsel(), `nterms_max` and output element `nterms` do not count the
  # intercept (whereas `vsel_expected$nterms_max` does):
  if (is.null(nterms_max_expected)) {
    nterms_ch <- vsel_expected$nterms_max - 1
  } else {
    nterms_ch <- nterms_max_expected
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

  return(invisible(TRUE))
}

# A helper function for testing the structure of a `data.frame` as returned by
# summary.vsel() in its output element `selection`
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
        expect_equal(smmry_sel[, diff_nm[stat_idx]], numeric(nrow(smmry_sel)),
                     info = info_str)
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
      expect_true(all(smmry_sel[, stats_mean_name[stat_idx]] >=
                        smmry_sel[, lower_nm[stat_idx]]),
                  info = info_str)
    }
  }
  if ("upper" %in% type_expected) {
    if (length(stats_expected) == 1) {
      upper_nm <- "upper"
    } else {
      upper_nm <- paste(stats_expected, "upper", sep = ".")
    }
    for (stat_idx in seq_along(stats_expected)) {
      expect_true(all(smmry_sel[, stats_mean_name[stat_idx]] <=
                        smmry_sel[, upper_nm[stat_idx]]),
                  info = info_str)
    }
  }

  return(invisible(TRUE))
}
