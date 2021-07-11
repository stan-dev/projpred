context("miscellaneous")

test_that(paste(
  "specifying `seed` for .get_refdist() correctly leads to reproducible",
  "results (and restores the RNG state afterwards)"
), {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      .Random.seed_orig1 <- .Random.seed
      refdist_orig <- .get_refdist(refmods[[mod_nm]][[fam_nm]], nclusters = 10,
                                   seed = seed_tst)
      .Random.seed_orig2 <- .Random.seed
      rand_orig <- runif(1) # Just to advance `.Random.seed[2]`.
      .Random.seed_new1 <- .Random.seed
      refdist_new <- .get_refdist(refmods[[mod_nm]][[fam_nm]], nclusters = 10,
                                  seed = seed_tst + 1L)
      .Random.seed_new2 <- .Random.seed
      rand_new <- runif(1) # Just to advance `.Random.seed[2]`.
      .Random.seed_repr1 <- .Random.seed
      refdist_repr <- .get_refdist(refmods[[mod_nm]][[fam_nm]], nclusters = 10,
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
  }
})
