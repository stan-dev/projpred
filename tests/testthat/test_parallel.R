# Setup -------------------------------------------------------------------

# withr resets the options and the parallel backend at the end of this file,
# also if the file stops with an error.
if (run_prll) {
  withr::local_options(list(projpred.parallel_proj_trigger = 0L))

  if (dopar_backend == "doParallel") {
    doParallel::registerDoParallel(ncores)
    withr::defer(doParallel::stopImplicitCluster())
  } else if (dopar_backend == "doFuture") {
    doFuture::registerDoFuture()
    withr::local_options(list(doFuture.foreach.export = ".export"))
    withr::defer(future::plan(future::sequential))
    if (future_plan == "multicore") {
      future::plan(future::multicore, workers = ncores)
    } else if (future_plan == "multisession") {
      future::plan(future::multisession, workers = ncores)
    } else if (future_plan == "callr") {
      future::plan(future.callr::callr, workers = ncores)
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
    # Use suppressWarnings() because test_that() somehow redirects stderr() and
    # so throws warnings that projpred wants to capture internally:
    cvvs_repr <- suppressWarnings(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]],
           cvfits = if (identical(args_cvvs_i$cv_method, "kfold")) {
             cvfitss[[args_cvvs_i$tstsetup_ref]]
           } else {
             refmods[[args_cvvs_i$tstsetup_ref]]$cvfits # should be `NULL`
           }),
      excl_nonargs(args_cvvs_i)
    )))
    expect_equal(cvvs_repr, cvvss[[tstsetup]], info = tstsetup)
  }
})

# The two tests below do not need `run_prll`: the first one registers its own
# (sequential) foreach backend, and the second one only reads code.

test_that(paste(
  "the fold function of the parallel cv_varsel() holds only the objects that",
  "it uses"
), {
  skip_if_not(run_cvvs)
  skip_if_not_installed("foreach")
  skip_if_not_installed("doRNG")
  skip_if_not_installed("iterators")
  # A GLM with K-fold fits (in the setup, these are the reference models
  # without observation weights), preferably with the traditional projection:
  tstsetups_ref <- names(cvfitss)[!sapply(cvfitss, inherits, "try-error")]
  tstsetup_ref <- c(grep("\\.glm\\..*\\.trad", tstsetups_ref, value = TRUE),
                    grep("\\.glm\\.", tstsetups_ref, value = TRUE))[1]
  skip_if(is.na(tstsetup_ref))

  # A foreach backend that records the objects in the environment of the fold
  # function (these objects go to the workers with backends that keep the
  # environment, e.g., doFuture) and then runs the iterations sequentially.
  # It also runs the foreach loops of the projection parallelization (option
  # `projpred.parallel_proj_trigger`), which have no fold function.
  env_nms <- list()
  do_record <- function(obj, expr, envir, data) {
    fun_nm <- intersect(c("one_obs", "one_fold"), obj$export)
    if (length(fun_nm) == 1L) {
      env_nms[[fun_nm]] <<- ls(environment(get(fun_nm, envir)),
                               all.names = TRUE)
    }
    it <- iterators::iter(obj)
    accumulator <- foreach::makeAccum(it)
    results <- lapply(as.list(it), function(args_i) {
      eval(expr, list2env(args_i, parent = envir))
    })
    accumulator(results, seq_along(results))
    foreach::getResult(it)
  }
  # Restore the registered backend at the end. foreach has no exported getter
  # for the registered backend, so read its internal store (the setup of this
  # file can have registered doParallel or doFuture):
  fe_globals <- foreach:::.foreachGlobals
  nms_reg <- c("fun", "data", "info")
  reg_old <- mget(nms_reg, envir = fe_globals, ifnotfound = list(NULL))
  withr::defer({
    if (all(nms_reg %in% ls(fe_globals))) rm(list = nms_reg, envir = fe_globals)
    if (!is.null(reg_old$fun)) {
      foreach::setDoPar(reg_old$fun, data = reg_old$data, info = reg_old$info)
    }
  })
  foreach::setDoPar(do_record, info = function(data, item) {
    switch(item, name = "record", version = "0", workers = 1L, NULL)
  })

  refmod <- refmods[[tstsetup_ref]]
  for (cv_meth in c("LOO", "kfold")) {
    # Use suppressWarnings() because of Pareto-k warnings of the small data:
    cvvs <- suppressWarnings(cv_varsel(
      refmod, cv_method = cv_meth,
      cvfits = if (cv_meth == "kfold") cvfitss[[tstsetup_ref]],
      validate_search = TRUE, nclusters = nclusters_tst,
      nclusters_pred = nclusters_pred_tst, nterms_max = nterms_max_tst,
      parallel = TRUE, seed = seed_tst, verbose = FALSE
    ))
    expect_true(inherits(cvvs, "vsel"), info = cv_meth)
  }
  expect_true("refmodel" %in% env_nms$one_obs)
  expect_false(any(c("psisloo", "loglik_forPSIS", "refdist_pred") %in%
                     env_nms$one_obs))
  expect_true("refmodel" %in% env_nms$one_fold)
  expect_false("list_cv" %in% env_nms$one_fold)
})

test_that(paste(
  "the object lists of the fold functions of cv_varsel() match their free",
  "variables"
), {
  skip_if_not_installed("codetools")
  # The right-hand sides of all assignments `nm <- ...` in `expr`:
  rhs_of <- function(expr, nm) {
    out <- list()
    walk <- function(x) {
      if (!is.call(x)) return(invisible())
      if (identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(nm))) {
        out[[length(out) + 1L]] <<- x[[3]]
      }
      for (j in seq_along(x)) {
        # Skip empty arguments, as in `x[, 1]`:
        if (is.symbol(x[[j]]) && identical(deparse(x[[j]]), "")) next
        walk(x[[j]])
      }
    }
    walk(expr)
    out
  }
  for (outer_nm in c("loo_varsel", "kfold_varsel")) {
    inner_nm <- c(loo_varsel = "one_obs", kfold_varsel = "one_fold")[[outer_nm]]
    outer_fun <- get(outer_nm, envir = asNamespace("projpred"))
    rhss <- rhs_of(body(outer_fun), inner_nm)
    is_def <- vapply(rhss, function(r) {
      is.call(r) && identical(r[[1]], as.name("function"))
    }, logical(1))
    is_slim <- vapply(rhss, function(r) {
      is.call(r) && identical(r[[1]], as.name("fun_with_objs"))
    }, logical(1))
    expect_identical(c(sum(is_def), sum(is_slim)), c(1L, 1L), info = outer_nm)
    inner_fun <- eval(rhss[[which(is_def)]])
    globs <- codetools::findGlobals(inner_fun, merge = FALSE)
    outer_locals <- unique(c(
      names(formals(outer_fun)),
      codetools::findFuncLocals(formals(outer_fun), body(outer_fun))
    ))
    listed <- eval(rhss[[which(is_slim)]][[3]])
    expect_setequal(intersect(globs$variables, outer_locals), listed)
    expect_length(intersect(globs$functions, outer_locals), 0)
  }
})
