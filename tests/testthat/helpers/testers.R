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
                      info_str = paste(info_str, j, sep = "_"),
                      ...)
  }
  if (is_seq) {
    # kl should be non-increasing on training data
    klseq <- sapply(p, function(x) sum(x$kl))
    # Remove intercept from the comparison:
    klseq <- klseq[-1]
    expect_identical(klseq, cummin(klseq), info = info_str)
    ### Check with tolerance:
    # expect_true(all(diff(klseq) - 1e-1 < 0), info = info_str)
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
  if (nprjdraws_expected == 1) {
    nprjdraws_sub_fit <- length(sub_fit_nms)
  } else {
    nprjdraws_sub_fit <- nprjdraws_expected
  }
  expect_length(p$sub_fit, nprjdraws_sub_fit)
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
    # Do not count the intercept, except for the intercept-only model:
    expect_length(p$solution_terms, max(solterms_expected, 1))
    # Same check, but using count_terms_chosen() which counts the intercept:
    expect_equal(count_terms_chosen(p$solution_terms), solterms_expected + 1,
                 info = info_str)
  } else if (is.character(solterms_expected)) {
    if (length(solterms_expected) == 0) {
      solterms_out <- "1"
    } else {
      solterms_out <- solterms_expected
    }
    expect_identical(p$solution_terms, solterms_out, info = info_str)
  }
  if (nprjdraws_expected == 1) {
    expect_identical(p$weights, 1, info = info_str)
  }
  return(invisible(TRUE))
}
