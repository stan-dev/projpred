context("parallel")

# Setup -------------------------------------------------------------------

if (run_prll) {
  trigger_default <- options(projpred.prll_prj_trigger = 0L)

  if (dopar_backend == "doParallel") {
    doParallel::registerDoParallel(ncores)
  } else if (dopar_backend == "doFuture") {
    doFuture::registerDoFuture()
    export_default <- options(doFuture.foreach.export = ".export")
    if (future_plan == "callr") {
      future::plan(future.callr::callr, workers = ncores)
    } else if (future_plan == "multisession") {
      future::plan(future::multisession, workers = ncores)
    } else if (future_plan == "multicore") {
      future::plan(future::multicore, workers = ncores)
    } else {
      stop("Unrecognized `future_plan`.")
    }
  } else {
    stop("Unrecognized `dopar_backend`.")
  }
  stopifnot(identical(foreach::getDoParWorkers(), ncores))
}

# project() ---------------------------------------------------------------

test_that("project() in parallel gives the same results as sequentially", {
  skip_if_not(run_prll)
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
    cvvs_repr <- suppressWarnings(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]]),
      excl_nonargs(args_cvvs_i)
    )))
    expect_equal(cvvs_repr, cvvss[[tstsetup]], info = tstsetup)
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
}
