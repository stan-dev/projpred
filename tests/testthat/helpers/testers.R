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
#   projection_tester()'s arguments `p`, `solterms_expected`, and `info_str`.
#
# @return `TRUE` (invisible).
#
proj_list_tester <- function(p,
                             len_expected = nterms_max_tst + 1L,
                             is_seq = TRUE,
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
    # kl should be non-increasing on training data
    klseq <- sapply(p, function(x) sum(x$kl))
    expect_identical(klseq, cummin(klseq), info = info_str)
    ### Check with tolerance:
    # expect_true(all(diff(klseq) < 1e-1), info = info_str)
    ###
  }
  return(invisible(TRUE))
}

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
#
projection_tester <- function(p,
                              solterms_expected,
                              nprjdraws_expected,
                              p_type_expected,
                              fam_expected = NULL,
                              prjdraw_weights_expected = NULL,
                              info_str = "") {
  expect_s3_class(p, "projection")
  # Check the names using `ignore.order = FALSE` because an incorrect
  # order would mean that the documentation of project()'s return value
  # would have to be updated:
  expect_named(p, projection_nms, info = info_str)
  if (nprjdraws_expected > 1) {
    expect_length(p$sub_fit, nprjdraws_expected)
  }
  expect_length(p$weights, nprjdraws_expected)
  expect_length(p$dis, nprjdraws_expected)
  # Number of projected draws in as.matrix.projection() (note that more
  # extensive tests for as.matrix.projection() may be found in
  # "test_as_matrix.R"):
  SW(nprjdraws <- NROW(as.matrix(p)))
  expect_identical(nprjdraws, nprjdraws_expected, info = info_str)
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
# @param nloo_expected The value which was used for argument `nloo` of
#   cv_varsel(). Leave empty (missing) to use the number of observations in the
#   original dataset (which is also the way to handle this argument for
#   `!with_cv`).
# @param info_str A single character string giving information to be printed in
#   case of failure.
#
# @return `TRUE` (invisible).
#
vsel_tester <- function(vs,
                        with_cv = FALSE,
                        refmod_expected,
                        dtest_expected = NULL,
                        solterms_len_expected,
                        method_expected,
                        cv_method_expected = NULL,
                        valsearch_expected = NULL,
                        ndraws_expected = ndraws_default,
                        ndraws_pred_expected = ndraws_pred_default,
                        nclusters_expected = NULL,
                        nclusters_pred_expected = NULL,
                        nloo_expected,
                        info_str = "") {
  dtest_type <- "train"
  if (with_cv) {
    vsel_nms <- vsel_nms_cv

    # Note: As mentioned in issue #167, component `"offset"` should in fact be
    # added here, too:
    dtest_nms <- setdiff(dtest_nms, "offset")
    # Need to re-order:
    dtest_nms <- dtest_nms[c(1, 5, 2, 4, 3)]

    vsel_smmrs_sub_nms <- c("lppd", "mu", "w")
    vsel_smmrs_ref_nms <- c("lppd", "mu")

    if (is.null(valsearch_expected)) {
      valsearch_expected <- TRUE
    }

    if (identical(cv_method_expected, "LOO")) {
      dtest_type <- "LOO"
      vsel_smmry_nms <- sub("^elpd$", "elpd.loo", vsel_smmry_nms)
    } else {
      stop("Probably need to adopt this.")
    }
  }
  method_expected <- tolower(method_expected)
  if (method_expected == "l1") {
    nclusters_expected <- 1
  }

  expect_s3_class(vs, "vsel")
  expect_named(vs, vsel_nms, info = info_str)

  # refmodel
  expect_s3_class(vs$refmodel, "refmodel")
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
  expect_equal(dim(vs$search_path$p_sel$mu), c(n_tst, nclusters_expected),
               info = info_str)
  if (vs$family$family == "gaussian") {
    expect_true(is.matrix(vs$search_path$p_sel$var), info = info_str)
    expect_type(vs$search_path$p_sel$var, "double")
    expect_equal(dim(vs$search_path$p_sel$var), c(n_tst, nclusters_expected),
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
    expect_identical(vs$d_test$y, vs$refmodel$y, info = info_str)
    expect_identical(vs$d_test$test_points, seq_len(n_tst), info = info_str)
    expect_null(vs$d_test$data, info = info_str)
    expect_identical(vs$d_test$weights, vs$refmodel$wobs, info = info_str)
    expect_identical(vs$d_test$type, dtest_type, info = info_str)
  } else {
    expect_identical(vs$d_test, dtest_expected, info = info_str)
  }

  # summaries
  expect_type(vs$summaries, "list")
  expect_named(vs$summaries, c("sub", "ref"), info = info_str)
  expect_type(vs$summaries$sub, "list")
  expect_length(vs$summaries$sub, solterms_len_expected + 1)
  nobsv <- nrow(vs$refmodel$fetch_data())
  if (missing(nloo_expected) ||
      is.null(nloo_expected) ||
      nloo_expected > nobsv) {
    nloo_expected <- nobsv
  }
  for (j in seq_along(vs$summaries$sub)) {
    expect_named(vs$summaries$sub[[!!j]], vsel_smmrs_sub_nms, info = info_str)
    expect_type(vs$summaries$sub[[!!j]]$mu, "double")
    expect_length(vs$summaries$sub[[!!j]]$mu, n_tst)
    expect_identical(sum(!is.na(vs$summaries$sub[[!!j]]$mu)),
                     nloo_expected, info = info_str)
    expect_type(vs$summaries$sub[[!!j]]$lppd, "double")
    expect_length(vs$summaries$sub[[!!j]]$lppd, n_tst)
    expect_identical(sum(!is.na(vs$summaries$sub[[!!j]]$lppd)),
                     nloo_expected, info = info_str)
    if (with_cv) {
      expect_type(vs$summaries$sub[[!!j]]$w, "double")
      expect_length(vs$summaries$sub[[!!j]]$w, n_tst)
      expect_true(all(!is.na(vs$summaries$sub[[!!j]]$w)), info = info_str)
      if (nloo_expected == nobsv) {
        expect_equal(vs$summaries$sub[[!!j]]$w, rep(1 / n_tst, n_tst),
                     info = info_str)
      } else {
        expect_true(any(vs$summaries$sub[[!!j]]$w != rep(1 / n_tst, n_tst)),
                    info = info_str)
      }
    }
  }
  expect_type(vs$summaries$ref, "list")
  expect_named(vs$summaries$ref, vsel_smmrs_ref_nms, info = info_str)
  expect_length(vs$summaries$ref$mu, n_tst)
  expect_true(all(!is.na(vs$summaries$ref$mu)), info = info_str)
  expect_length(vs$summaries$ref$lppd, n_tst)
  expect_true(all(!is.na(vs$summaries$ref$lppd)), info = info_str)

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
  expect_identical(vs$kl, cummin(vs$kl), info = info_str)
  ### Check with tolerance:
  # expect_true(all(diff(vs$kl) < 1e-1), info = info_str)
  ###

  # pct_solution_terms_cv
  if (with_cv) {
    expect_true(is.matrix(vs$pct_solution_terms_cv), info = info_str)
    expect_type(vs$pct_solution_terms_cv, "double")
    expect_identical(dim(vs$pct_solution_terms_cv),
                     c(solterms_len_expected, 1L + solterms_len_expected),
                     info = info_str)
    expect_identical(vs$pct_solution_terms_cv[, "size"],
                     as.numeric(seq_len(solterms_len_expected)),
                     info = info_str)
    if (isFALSE(vs$validate_search)) {
      expect_true(
        all(vs$pct_solution_terms_cv[
          ,
          -grep("size", colnames(vs$pct_solution_terms_cv))
        ] %in% c(0, 1)),
        info = info_str
      )
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
  expect_s3_class(vs$summary, "data.frame")
  expect_named(vs$summary, vsel_smmry_nms, info = info_str)
  expect_identical(nrow(vs$summary), solterms_len_expected + 1L,
                   info = info_str)
  expect_identical(vs$summary$size, seq_len(nrow(vs$summary)) - 1,
                   info = info_str)
  expect_identical(vs$summary$solution_terms,
                   c(NA_character_, vs$solution_terms),
                   info = info_str)
  expect_equal(diff(vs$summary[, grep("^elpd", vsel_smmry_nms, value = TRUE)]),
               diff(vs$summary$diff), info = info_str)

  return(invisible(TRUE))
}
