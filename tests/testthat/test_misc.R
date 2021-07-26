context("miscellaneous")

test_that(paste(
  ".get_refdist(): `seed` works (and restores the RNG state afterwards)"
), {
  for (tstsetup in names(refmods)) {
    .Random.seed_orig1 <- .Random.seed
    refdist_orig <- .get_refdist(refmods[[tstsetup]], nclusters = 10,
                                 seed = seed_tst)
    .Random.seed_orig2 <- .Random.seed
    rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_new1 <- .Random.seed
    refdist_new <- .get_refdist(refmods[[tstsetup]], nclusters = 10,
                                seed = seed_tst + 1L)
    .Random.seed_new2 <- .Random.seed
    rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
    .Random.seed_repr1 <- .Random.seed
    refdist_repr <- .get_refdist(refmods[[tstsetup]], nclusters = 10,
                                 seed = seed_tst)
    .Random.seed_repr2 <- .Random.seed
    # Expected equality:
    expect_equal(refdist_repr, refdist_orig, info = tstsetup)
    expect_equal(.Random.seed_orig2, .Random.seed_orig1, info = tstsetup)
    expect_equal(.Random.seed_new2, .Random.seed_new1, info = tstsetup)
    expect_equal(.Random.seed_repr2, .Random.seed_repr1, info = tstsetup)
    # Expected inequality:
    expect_false(isTRUE(all.equal(refdist_new, refdist_orig)),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(rand_new, rand_orig)), info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_new2, .Random.seed_orig2)),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_repr2, .Random.seed_orig2)),
                 info = tstsetup)
    expect_false(isTRUE(all.equal(.Random.seed_repr2, .Random.seed_new2)),
                 info = tstsetup)
  }
})

test_that("rstanarm: special formulas work", {
  tstsetups <- grep("\\.spclformul", names(fits), value = TRUE)
  # Compare the "special formula" fit with the corresponding "standard formula"
  # fit (which does not have the special formula elements):
  for (tstsetup in tstsetups) {
    mf_spclformul <- fits[[tstsetup]]$model
    nms_spclformul <- setdiff(
      grep("y_|xco", names(mf_spclformul), value = TRUE),
      "xco.1"
    )

    tstsetup_stdformul <- sub("\\.spclformul", ".stdformul", tstsetup)
    stopifnot(tstsetup_stdformul != tstsetup)
    mf_stdformul <- fits[[tstsetup_stdformul]]$model
    nms_stdformul <- setdiff(
      grep("y_|xco", names(mf_stdformul), value = TRUE),
      "xco.1"
    )

    expect_equal(mf_spclformul[, setdiff(names(mf_spclformul), nms_spclformul)],
                 mf_stdformul[, setdiff(names(mf_stdformul), nms_stdformul)],
                 info = tstsetup)
    # Check arithmetic expressions:
    for (nm_spclformul in nms_spclformul) {
      expect_equal(mf_spclformul[, nm_spclformul],
                   eval(str2lang(nm_spclformul), envir = dat),
                   info = paste(tstsetup, nm_spclformul, sep = "__"))
    }
  }
})
