context("parallel")

# Setup -------------------------------------------------------------------

if (run_prll) {
  cv_prev <- options(projpred.prll_cv = FALSE)
  trigger_default <- options(projpred.prll_prj_trigger = 0L)
}

# project() ---------------------------------------------------------------

test_that("project() in parallel gives the same results as sequentially", {
  skip_if_not(run_prll)
  skip_if_not(run_prj)
  tstsetups <- grep("\\.glm\\.", names(prjs), value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    p_repr <- do.call(project, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]]),
      excl_nonargs(args_prj_i)
    ))
    expect_equal(p_repr, prjs[[tstsetup]], info = tstsetup)
  }
})

# varsel() ----------------------------------------------------------------

test_that("varsel() in parallel gives the same results as sequentially", {
  skip_if_not(run_prll)
  skip_if_not(run_vs)
  tstsetups <- grep("\\.glm\\.", names(vss), value = TRUE)
  for (tstsetup in tstsetups) {
    args_vs_i <- args_vs[[tstsetup]]
    vs_repr <- do.call(varsel, c(
      list(object = refmods[[args_vs_i$tstsetup_ref]]),
      excl_nonargs(args_vs_i)
    ))
    expect_equal(vs_repr, vss[[tstsetup]], info = tstsetup)
  }
})

# cv_varsel() -------------------------------------------------------------

test_that("cv_varsel() in parallel gives the same results as sequentially", {
  skip_if_not(run_prll)
  skip_if_not(run_cvvs)
  tstsetups <- grep("\\.glm\\.", names(cvvss), value = TRUE)
  for (tstsetup in tstsetups) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    # Use suppressWarnings() because of occasional warnings concerning Pareto k
    # diagnostics:
    cvvs_expr <- expression(suppressWarnings(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]]),
      excl_nonargs(args_cvvs_i)
    ))))
    cvvs_repr <- eval(cvvs_expr)
    ### For the original (expected) cv_varsel() output, we cannot simply use
    ### `cvvss[[tstsetup]]` because `cvvss[[tstsetup]]` made use of the CV
    ### parallelization. So we need to run cv_varsel() again, but sequentially:
    trigger_prev <- options(trigger_default)
    cvvs_orig <- eval(cvvs_expr)
    options(trigger_prev)
    ###
    expect_equal(cvvs_repr, cvvs_orig, info = tstsetup)
  }
})

# Teardown ----------------------------------------------------------------

if (run_prll) {
  if (dopar_backend == "doParallel") {
    doParallel::stopImplicitCluster()
  } else if (dopar_backend == "doFuture") {
    future::plan(future::sequential)
    options(export_default)
    rm(export_default)
  } else {
    stop("Unrecognized `dopar_backend`.")
  }

  options(trigger_default)
  rm(trigger_default)
  options(cv_prev)
  rm(cv_prev)
}
