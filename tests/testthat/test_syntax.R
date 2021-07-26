context("special formulas")

test_that("special formulas work", {
  tstsetup <- "glm__gauss__spclformul"
  SW(fit_glm_gauss_spclformul <- rstanarm::stan_glm(
    log(abs(y_glm_gauss)) ~
      xco.1 + I(xco.1^2) + log(abs(xco.2)) * poly(xco.3, 3) + xca.1 + xca.2,
    family = f_gauss, data = dat,
    weights = wobs_tst, offset = offs_tst,
    chains = chains_tst, seed = seed_tst, iter = iter_tst, QR = TRUE
  ))
  mf_spclformul <- fit_glm_gauss_spclformul$model

  ### TODO: Move this to `test_misc.R`:
  # Compare with standard fit (excluding the special formula elements):
  nms_spclformul <- setdiff(grep("y_|xco", names(mf_spclformul), value = TRUE),
                            "xco.1")
  mf <- fits$glm$gauss$model
  nms <- setdiff(grep("y_|xco", names(mf), value = TRUE), "xco.1")
  expect_equal(mf_spclformul[, setdiff(names(mf_spclformul), nms_spclformul)],
               mf[, setdiff(names(mf), nms)], info = tstsetup)
  # Check arithmetic expressions:
  expect_equal(mf_spclformul$`log(abs(y_glm_gauss))`, log(abs(dat$y_glm_gauss)),
               info = tstsetup)
  for (nm_spclformul in nms_spclformul) {
    expect_equal(mf_spclformul[, nm_spclformul],
                 eval(str2lang(nm_spclformul), envir = dat),
                 info = paste(tstsetup, nm_spclformul, sep = "__"))
  }
  ###

  # TODO: Add tests for refmod_spclformul$div_minimizer() to check that the
  # arithmetic expressions on the right-hand side of the formula are taken there
  # into account.
})
