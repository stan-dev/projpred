prj_vsel_str_tester <- function(p) {
  expect_type(p, "list")
  expect_length(p, nterms + 1)
  expect_true(.is_proj_list(p), info = tstsetup)

  prjdraw_weights <- p[[1]]$weights
  for (j in seq_along(p)) {
    expect_s3_class(p[[!!j]], "projection")
    # Check the names using `ignore.order = FALSE` because an incorrect
    # order would mean that the documentation of project()'s return value
    # would have to be updated:
    expect_named(p[[!!j]], projection_nms, info = tstsetup)
    # Number of projected draws should be equal to the default of `ndraws`
    # (note that more extensive tests for as.matrix.projection() may be
    # found in "test_as_matrix.R"):
    expect_length(p[[!!j]]$sub_fit, nclusters_pred_tst)
    expect_length(p[[!!j]]$weights, nclusters_pred_tst)
    expect_length(p[[!!j]]$dis, nclusters_pred_tst)
    SW(nprjdraws <- NROW(as.matrix(p[[!!j]])))
    expect_identical(nprjdraws, nclusters_pred_tst, info = tstsetup)
    # The j-th element should have j solution terms (usually excluding the
    # intercept, but counting it for `j == 1`):
    expect_length(p[[!!j]]$solution_terms, max(j - 1, 1))
    # Same check, but using count_terms_chosen():
    expect_equal(count_terms_chosen(p[[!!j]]$solution_terms), !!j, info = tstsetup)
    expect_identical(p[[!!j]]$family, vs_list[[!!i]]$family, info = tstsetup)
    # All submodels should use the same clustering:
    expect_identical(p[[!!j]]$weights, prjdraw_weights, info = tstsetup)
  }
  # kl should be non-increasing on training data
  klseq <- sapply(p, function(x) sum(x$kl))
  # Remove intercept from the comparison:
  klseq <- klseq[-1]
  expect_identical(klseq, cummin(klseq), info = tstsetup)
  ### Check with tolerance:
  # expect_true(all(diff(klseq) - 1e-1 < 0), info = tstsetup)
  ###
  return(invisible(TRUE))
}
