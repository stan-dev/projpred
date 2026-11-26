test_that("is_known_setup_err() accepts an error only for its model types", {
  msg_bs <- "singular matrix in 'backsolve'. First zero in diagonal [1]"
  expect_true(is_known_setup_err(msg_bs, mod_nm = "gamm"))
  expect_false(is_known_setup_err(msg_bs, mod_nm = "glmm"))
  expect_false(is_known_setup_err(msg_bs, mod_nm = "gam"))
  expect_false(is_known_setup_err(msg_bs, mod_nm = "gamm", kfold = TRUE))
  msg_pw <- "pwrssUpdate did not converge in (maxit) iterations"
  expect_true(is_known_setup_err(msg_pw, mod_nm = "glmm", kfold = TRUE))
  expect_true(is_known_setup_err(msg_pw, mod_nm = "gamm", kfold = TRUE))
  expect_false(is_known_setup_err(msg_pw, mod_nm = "glm", kfold = TRUE))
  expect_false(is_known_setup_err(msg_pw, mod_nm = "gam", kfold = TRUE))
  expect_false(is_known_setup_err(msg_pw, mod_nm = "glmm"))
  expect_false(is_known_setup_err("other error", mod_nm = "gamm"))
  expect_false(is_known_setup_err(NULL, mod_nm = "gamm"))
  expect_false(is_known_setup_err(msg_bs, mod_nm = NULL))
})

test_that("the getters report the dropped test setups", {
  # A getter (see helpers/fixtures.R) drops a test setup only after a known
  # upstream error (see is_known_setup_err()). This test builds all objects of
  # the getters, one object type after the other. The skip message lists the
  # dropped test setups. (The fits are never dropped; a failed fit stops this
  # test.)
  stages <- list(fit = names(args_fit), ref = names(args_ref))
  if (run_vs) {
    stages <- c(stages, list(vs = names(args_vs)))
  }
  if (run_cvvs) {
    stages <- c(stages, list(cvfits = tstsetups_ref_kfold,
                             cvvs = names(args_cvvs)))
  }
  if (run_prj) {
    stages <- c(stages, list(prj = names(args_prj), pl = names(args_prj),
                             pe = names(args_prj), pp = names(args_prj)))
  }
  for (vsel_type in c("vs", "cvvs")) {
    if (!get(paste0("run_", vsel_type))) {
      next
    }
    nms_prj <- names(get(paste0("args_prj_", vsel_type)))
    stages_vsel <- list(
      prj = nms_prj, pl = nms_prj, pe = nms_prj, pp = nms_prj,
      smmry = names(get(paste0("args_smmry_", vsel_type))),
      perf = names(get(paste0("args_perf_", vsel_type))),
      plot = names(get(paste0("args_plot_", vsel_type))),
      rk = names(get(paste0("args_rk_", vsel_type)))
    )
    names(stages_vsel) <- paste0(names(stages_vsel), "_", vsel_type)
    stages <- c(stages, stages_vsel)
  }
  if (run_cvvs) {
    stages <- c(stages, list(pr_cvvs = names(args_pr_cvvs),
                             plotpr = names(args_plotpr)))
  }
  drops <- lapply(setNames(nm = names(stages)), function(type) {
    dropped_tst(type, stages[[type]])
  })
  drops <- drops[lengths(drops) > 0]
  if (length(drops) > 0) {
    skip(paste0(
      "The getters dropped test setups after a known upstream error: ",
      paste(names(drops), vapply(drops, paste, character(1), collapse = ", "),
            sep = ": ", collapse = "; ")
    ))
  }
  expect_length(drops, 0)
})
