# A helper function for testing the structure of an expected "projection" object
#
# @param p An object of class "projection" (at least expected so).
# @param solterms_expected Either a single numeric value giving the expected
#   number of solution terms (not counting the intercept, even for the
#   intercept-only model), a character vector giving the expected solution
#   terms, or `NULL` for not testing the solution terms at all.
# @param nprjdraws_expected A single numeric value giving the expected number of
#   projected draws.
# @param p_type_expected A single logical value giving the expected value for
#   `p$p_type`.
# @param fam_expected The expected "family" object or `NULL` for not testing the
#   family object at all.
# @param prjdraw_weights_expected The expected weights for the projected draws
#   or `NULL` for not testing these weights at all.
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
projection_tester <- function(p,
                              solterms_expected,
                              nprjdraws_expected,
                              p_type_expected,
                              fam_expected = NULL,
                              prjdraw_weights_expected = NULL,
                              from_datafit = FALSE,
                              info_str = "") {
  expect_s3_class(p, "projection")
  # Check the names using `ignore.order = FALSE` because an incorrect
  # order would mean that the documentation of project()'s return value
  # would have to be updated:
  expect_named(p, projection_nms, info = info_str)

  if (from_datafit) {
    subfit_nms <- setdiff(subfit_nms, "y")
  }

  if (nprjdraws_expected > 1) {
    expect_length(p$sub_fit, nprjdraws_expected)
    sub_fit_totest <- p$sub_fit
  } else {
    sub_fit_totest <- list(p$sub_fit)
  }
  has_grp <- formula_contains_group_terms(p$refmodel$formula)
  has_add <- formula_contains_additive_terms(p$refmodel$formula)
  for (j in seq_along(sub_fit_totest)) {
    if (!has_grp && !has_add) {
      expect_s3_class(sub_fit_totest[[!!j]], "subfit")
      expect_type(sub_fit_totest[[!!j]], "list")
      expect_named(sub_fit_totest[[!!j]], subfit_nms, info = info_str)
      expect_length(sub_fit_totest[[!!j]]$alpha, 1)
      if (length(p$solution_terms) > 0) {
        expect_identical(ncol(sub_fit_totest[[!!j]]$beta), 1L, info = info_str)
        solterms_len_expected <- sum(sapply(p$solution_terms, function(trm_i) {
          ncol(model.matrix(
            as.formula(paste("~ 0 +", trm_i)),
            data = p$refmodel$fetch_data()
          ))
        }))
        ### As discussed in issue #149, the following might be more appropriate:
        # solterms_len_expected <- ncol(model.matrix(
        #   as.formula(paste("~", paste(p$solution_terms, collapse = " + "))),
        #   data = p$refmodel$fetch_data()
        # )) - 1L
        ###
        expect_equal(nrow(sub_fit_totest[[!!j]]$beta), solterms_len_expected,
                     info = info_str)
      } else {
        if (!from_datafit) {
          expect_identical(dim(sub_fit_totest[[!!j]]$beta), c(0L, 1L),
                           info = info_str)
        } else {
          expect_null(sub_fit_totest[[!!j]]$beta, info = info_str)
        }
      }
    } else if (has_grp && !has_add) {
      inherits(sub_fit_totest[[!!j]], c("lmerMod", "glmerMod"))
    } else if (has_add) {
      stop("Still to-do.")
    }
  }
  expect_length(p$weights, nprjdraws_expected)
  expect_length(p$dis, nprjdraws_expected)
  if (!from_datafit) {
    # Number of projected draws in as.matrix.projection() (note that more
    # extensive tests for as.matrix.projection() may be found in
    # "test_as_matrix.R"):
    SW(nprjdraws <- NROW(as.matrix(p)))
    expect_identical(nprjdraws, nprjdraws_expected, info = info_str)
  }
  expect_identical(p$p_type, p_type_expected, info = info_str)
  if (!is.null(fam_expected)) {
    expect_identical(p$family, fam_expected, info = info_str)
  }
  if (!is.null(prjdraw_weights_expected)) {
    expect_identical(p$weights, prjdraw_weights_expected, info = info_str)
  }
  if (is.numeric(solterms_expected)) {
    expect_length(p$solution_terms, solterms_expected)
    # Same check, but using count_terms_chosen():
    expect_equal(count_terms_chosen(p$solution_terms, add_icpt = TRUE),
                 solterms_expected + 1, info = info_str)
  } else if (is.character(solterms_expected)) {
    expect_identical(p$solution_terms, solterms_expected, info = info_str)
  }
  if (nprjdraws_expected == 1) {
    expect_identical(p$weights, 1, info = info_str)
  }
  return(invisible(TRUE))
}
# A helper function for testing the structure of an expected "proj_list" object
#
# @param p An object of (informal) class "proj_list" (at least expected so).
# @param len_expected The expected length of `p`.
# @param is_seq A single logical value indicating whether `p` is expected to be
#   sequential (i.e., the number of solution terms increases by 1 from one
#   element of `p` to the next).
# @param info_str A single character string giving information to be printed in
#   case of failure.
# @param ... Arguments passed to projection_tester(), apart from
#   projection_tester()'s arguments `p`, `solterms_expected`, `from_datafit`,
#   and `info_str`.
#
# @return `TRUE` (invisible).
proj_list_tester <- function(p,
                             len_expected = nterms_max_tst + 1L,
                             is_seq = TRUE,
                             from_datafit = FALSE,
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
                      from_datafit = from_datafit,
                      info_str = paste(info_str, j, sep = "__"),
                      ...)
  }
  if (is_seq) {
    # kl should be non-increasing on training data
    klseq <- sapply(p, function(x) sum(x$kl))
    if (from_datafit) {
      # For some "datafit"s, we need to allow for a certain tolerance:
      expect_true(all(diff(klseq) < 3e-1), info = info_str)
    } else {
      expect_identical(klseq, cummin(klseq), info = info_str)
    }
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

# A helper function for testing the structure of an expected "vsel" object
#
# @param vs An object of class "vsel" (at least expected so).
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
# @param nclusters_expected The expected `vs$nclusters` object (not adopted for
#   L1 search).
# @param nclusters_pred_expected The expected `vs$nclusters_pred` object.
# @param nloo_expected Only relevant if `with_cv` is `TRUE`. The value which was
#   used for argument `nloo` of cv_varsel().
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
  nloo_expected = NULL,
  info_str = ""
) {
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

  expect_s3_class(vs, "vsel")
  expect_named(vs, vsel_nms, info = info_str)

  # refmodel
  refmod_class_expected <- "refmodel"
  if (from_datafit) {
    refmod_class_expected <- c("datafit", refmod_class_expected)
  }
  expect_s3_class(vs$refmodel, refmod_class_expected, exact = TRUE)
  expect_identical(vs$refmodel, refmod_expected, info = info_str)

  # search_path
  expect_type(vs$search_path, "list")
  expect_named(vs$search_path, searchpth_nms, info = info_str)
  expect_identical(vs$search_path$solution_terms, vs$solution_terms,
                   info = info_str)
  expect_type(vs$search_path$sub_fits, "list")
  expect_length(vs$search_path$sub_fits, solterms_len_expected + 1)
  for (i in seq_along(vs$search_path$sub_fits)) {
    if (!is.null(nclusters_expected) && nclusters_expected == 1) {
      sub_fits_totest <- list(vs$search_path$sub_fits[[i]])
    } else {
      expect_length(vs$search_path$sub_fits[[!!i]], nclusters_expected)
      sub_fits_totest <- vs$search_path$sub_fits[[i]]
    }
    for (j in seq_along(sub_fits_totest)) {
      expect_true(inherits(sub_fits_totest[[!!j]],
                           get_as.matrix_cls_projpred()),
                  info = paste(info_str, i, sep = "__"))
    }
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
  # Expected to be decreasing:
  if (!from_datafit) {
    expect_identical(vs$kl, cummin(vs$kl), info = info_str)
  } else {
    # For some "datafit"s, we need to allow for a certain tolerance:
    expect_true(all(diff(vs$kl) < 3e-1), info = info_str)
  }

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
    expect_true(all(pct_solterms >= 0 & pct_solterms <= 1), info = info_str)
    ### Excluded because of issue #173:
    # if (isFALSE(vs$validate_search)) {
    #   expect_true(all(pct_solterms %in% c(0, 1)), info = info_str)
    #   # More specifically:
    #   pct_solterms_ch <- matrix(0, nrow = nrow(pct_solterms),
    #                             ncol = ncol(pct_solterms))
    #   pct_solterms_ch[lower.tri(pct_solterms_ch)] <- 1
    #   diag(pct_solterms_ch) <- 1
    #   colnames(pct_solterms_ch) <- pct_nonsize_nms
    #   expect_identical(pct_solterms_ch, pct_solterms, info = info_str)
    # }
    ###
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
  smmry_sel_tester(smmry$selection, nterms_max_expected = nterms_max_expected,
                   info_str = info_str, ...)

  return(invisible(TRUE))
}

# A helper function for testing the structure of a `data.frame` as returned by
# summary.vsel() in its output element `selection`
#
# @param smmry_sel A `data.frame` as returned by summary.vsel() in its output
#   element `selection`.
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
