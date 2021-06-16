# A helper function for testing the output structure from project() applied to
# a "vsel" object
#
# @param p An object returned by project().
# @param fam_expected The expected "family" object.
# @param info_str A single character string giving info to be printed in case of
#   failure.
#
# @return `TRUE` (invisible).
#
prj_vsel_str_tester <- function(p, fam_expected, info_str = "") {
  expect_type(p, "list")
  expect_length(p, nterms_max_tst + 1)
  expect_true(.is_proj_list(p), info = info_str)

  prjdraw_weights <- p[[1]]$weights
  for (j in seq_along(p)) {
    expect_s3_class(p[[!!j]], "projection")
    # Check the names using `ignore.order = FALSE` because an incorrect
    # order would mean that the documentation of project()'s return value
    # would have to be updated:
    expect_named(p[[!!j]], projection_nms, info = info_str)
    # Number of projected draws should be equal to the default of `ndraws`
    # (note that more extensive tests for as.matrix.projection() may be
    # found in "test_as_matrix.R"):
    expect_length(p[[!!j]]$sub_fit, nclusters_pred_tst)
    expect_length(p[[!!j]]$weights, nclusters_pred_tst)
    expect_length(p[[!!j]]$dis, nclusters_pred_tst)
    SW(nprjdraws <- NROW(as.matrix(p[[!!j]])))
    expect_identical(nprjdraws, nclusters_pred_tst, info = info_str)
    # The j-th element should have j solution terms (usually excluding the
    # intercept, but counting it for `j == 1`):
    expect_length(p[[!!j]]$solution_terms, max(j - 1, 1))
    # Same check, but using count_terms_chosen():
    expect_equal(count_terms_chosen(p[[!!j]]$solution_terms), !!j,
                 info = info_str)
    expect_identical(p[[!!j]]$family, fam_expected, info = info_str)
    # All submodels should use the same clustering:
    expect_identical(p[[!!j]]$weights, prjdraw_weights, info = info_str)
  }
  # kl should be non-increasing on training data
  klseq <- sapply(p, function(x) sum(x$kl))
  # Remove intercept from the comparison:
  klseq <- klseq[-1]
  expect_identical(klseq, cummin(klseq), info = info_str)
  ### Check with tolerance:
  # expect_true(all(diff(klseq) - 1e-1 < 0), info = info_str)
  ###
  return(invisible(TRUE))
}

# A helper function for testing the output structure from project() applied to
# a "refmodel" object
#
# @param p An object returned by project().
# @param info_str A single character string giving info to be printed in case of
#   failure.
#
# @return `TRUE` (invisible).
#
prj_refmodel_str_tester <- function(p, solterms_expected, nprjdraws_expected,
                                    info_str = "") {
  expect_s3_class(p, "projection")
  expect_named(p, projection_nms, info = info_str)
  if (nprjdraws_expected == 1) {
    nprjdraws_sub_fit <- length(sub_fit_nms)
  } else {
    nprjdraws_sub_fit <- nprjdraws_expected
  }
  expect_length(p$sub_fit, nprjdraws_sub_fit)
  expect_length(p$weights, nprjdraws_expected)
  expect_length(p$dis, nprjdraws_expected)
  SW(nprjdraws <- NROW(as.matrix(p)))
  expect_identical(nprjdraws, nprjdraws_expected, info = info_str)
  solterms_out <- if (length(solterms_expected) == 0) "1" else solterms_expected
  expect_identical(p$solution_terms, solterms_out, info = info_str)
  return(invisible(TRUE))
}
