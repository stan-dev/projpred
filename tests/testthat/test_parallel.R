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

# The four tests below do not need `run_prll`: the first three register their
# own (sequential) foreach backend, and the fourth one only reads code.

# A foreach backend for the tests below. It records the objects in the
# enclosing environment of the fold function (these objects go to the workers
# with backends that keep the environment, e.g., doFuture) and the option
# `doFuture.foreach.export` during the fold loops. For the loops of the
# projection parallelization (option `projpred.parallel_proj_trigger`), it
# records the same for each call, in `draw_env_nms` and `draw_export_opt`. Then
# it runs the iterations sequentially. With `export_only = TRUE`, the loop body
# sees only the `.export` objects, the iteration arguments and the search path,
# as with doFuture and `doFuture.foreach.export = ".export"`. With
# `fail = TRUE`, it stops in the first fold loop. When `frame` exits, the
# previous backend is registered again. foreach has no exported getter for the
# registered backend, so this reads its internal store (the setup of this file
# can have registered doParallel or doFuture).
local_record_backend <- function(fail = FALSE, export_only = FALSE,
                                 frame = parent.frame()) {
  rec <- new.env()
  rec$env_nms <- list()
  rec$export_opt <- list()
  rec$draw_env_nms <- list()
  rec$draw_export_opt <- list()
  do_record <- function(obj, expr, envir, data) {
    fun_nm <- intersect(c("one_obs", "one_fold"), obj$export)
    if (length(fun_nm) == 1L) {
      rec$env_nms[[fun_nm]] <- ls(environment(get(fun_nm, envir)),
                                  all.names = TRUE)
      rec$export_opt[fun_nm] <- list(getOption("doFuture.foreach.export"))
      if (fail) stop("Stop of the recording foreach backend.")
    }
    if ("one_draw" %in% obj$export) {
      rec$draw_env_nms <- c(rec$draw_env_nms, list(
        ls(environment(get("one_draw", envir)), all.names = TRUE)
      ))
      rec$draw_export_opt <- c(rec$draw_export_opt,
                               list(getOption("doFuture.foreach.export")))
    }
    if (export_only) {
      envir <- list2env(mget(obj$export, envir = envir, inherits = TRUE),
                        parent = globalenv())
    }
    it <- iterators::iter(obj)
    accumulator <- foreach::makeAccum(it)
    results <- lapply(as.list(it), function(args_i) {
      eval(expr, list2env(args_i, parent = envir))
    })
    accumulator(results, seq_along(results))
    foreach::getResult(it)
  }
  fe_globals <- foreach:::.foreachGlobals
  nms_reg <- c("fun", "data", "info")
  reg_old <- mget(nms_reg, envir = fe_globals, ifnotfound = list(NULL))
  withr::defer({
    if (all(nms_reg %in% ls(fe_globals))) rm(list = nms_reg, envir = fe_globals)
    if (!is.null(reg_old$fun)) {
      foreach::setDoPar(reg_old$fun, data = reg_old$data, info = reg_old$info)
    }
  }, envir = frame)
  foreach::setDoPar(do_record, info = function(data, item) {
    switch(item, name = "record", version = "0", workers = 1L, NULL)
  })
  rec
}

# A GLM with K-fold fits (in the setup, these are the reference models without
# observation weights), preferably with the traditional projection (`NA` if
# there is none):
tstsetup_ref_par <- function() {
  tstsetups_ref <- names(cvfitss)[!sapply(cvfitss, inherits, "try-error")]
  c(grep("\\.glm\\..*\\.trad", tstsetups_ref, value = TRUE),
    grep("\\.glm\\.", tstsetups_ref, value = TRUE))[1]
}

cvvs_par <- function(tstsetup_ref, cv_meth) {
  # Use suppressWarnings() because of Pareto-k warnings of the small data:
  suppressWarnings(cv_varsel(
    refmods[[tstsetup_ref]], cv_method = cv_meth,
    cvfits = if (cv_meth == "kfold") cvfitss[[tstsetup_ref]],
    validate_search = TRUE, nclusters = nclusters_tst,
    nclusters_pred = nclusters_pred_tst, nterms_max = nterms_max_tst,
    parallel = TRUE, seed = seed_tst, verbose = FALSE
  ))
}

test_that(paste(
  "the fold function of the parallel cv_varsel() holds only the objects that",
  "it uses"
), {
  skip_if_not(run_cvvs)
  skip_if_not_installed("foreach")
  skip_if_not_installed("doRNG")
  skip_if_not_installed("iterators")
  tstsetup_ref <- tstsetup_ref_par()
  skip_if(is.na(tstsetup_ref))
  rec <- local_record_backend()
  for (cv_meth in c("LOO", "kfold")) {
    cvvs <- cvvs_par(tstsetup_ref, cv_meth)
    expect_true(inherits(cvvs, "vsel"), info = cv_meth)
  }
  expect_true("refmodel" %in% rec$env_nms$one_obs)
  expect_false(any(c("psisloo", "loglik_forPSIS", "refdist_pred") %in%
                     rec$env_nms$one_obs))
  expect_true("refmodel" %in% rec$env_nms$one_fold)
  expect_false("list_cv" %in% rec$env_nms$one_fold)
})

test_that(paste(
  "the parallel cv_varsel() sets `doFuture.foreach.export = \".export\"` only",
  "for its fold loops"
), {
  skip_if_not(run_cvvs)
  skip_if_not_installed("foreach")
  skip_if_not_installed("doRNG")
  skip_if_not_installed("iterators")
  tstsetup_ref <- tstsetup_ref_par()
  skip_if(is.na(tstsetup_ref))
  rec <- local_record_backend()
  # A value set by the user, and the option not set (the default):
  for (opt_user in list(".export-and-automatic", NULL)) {
    withr::local_options(list(doFuture.foreach.export = opt_user))
    for (cv_meth in c("LOO", "kfold")) {
      rec$export_opt <- list()
      cvvs_par(tstsetup_ref, cv_meth)
      fun_nm <- c(LOO = "one_obs", kfold = "one_fold")[[cv_meth]]
      info_i <- paste(cv_meth, if (is.null(opt_user)) "NULL" else opt_user)
      expect_identical(rec$export_opt[[fun_nm]], ".export", info = info_i)
      expect_identical(getOption("doFuture.foreach.export"), opt_user,
                       info = info_i)
    }
  }
  # An error in the fold loop:
  withr::local_options(list(doFuture.foreach.export = ".export-and-automatic"))
  local({
    local_record_backend(fail = TRUE)
    for (cv_meth in c("LOO", "kfold")) {
      expect_error(cvvs_par(tstsetup_ref, cv_meth),
                   "Stop of the recording foreach backend", info = cv_meth)
      expect_identical(getOption("doFuture.foreach.export"),
                       ".export-and-automatic", info = cv_meth)
    }
  })
})

test_that(paste(
  "the parallel projection of single draws works with only the `.export`",
  "objects and gives the same results as the sequential one"
), {
  skip_if_not(run_prj)
  skip_if_not_installed("foreach")
  skip_if_not_installed("iterators")
  # The GLM projections fit all draws in one call, so use a multilevel
  # submodel (loop in divmin()) and an augmented-data projection (loop in
  # divmin_augdat()):
  first_prj <- function(cond) {
    tstsetups_cond <- names(Filter(cond, args_prj[names(prjs)]))
    c(tstsetups_cond, NA_character_)[1]
  }
  tstsetups <- c(
    divmin = first_prj(function(args_prj_i) {
      args_prj_i$mod_nm == "glmm" && args_prj_i$prj_nm == "trad" &&
        any(grepl("|", args_prj_i$predictor_terms, fixed = TRUE))
    }),
    divmin_augdat = first_prj(function(args_prj_i) {
      args_prj_i$prj_nm == "augdat" && args_prj_i$fam_nm == "brnll"
    })
  )
  skip_if(all(is.na(tstsetups)))
  nms_expected <- list(
    divmin = c("capt_fun", "dot_args", "projpred_random", "sdivmin"),
    divmin_augdat = c("data", "dot_args", "family", "formula",
                      "projpred_formula_no_random", "projpred_random",
                      "sdivmin")
  )
  withr::local_options(list(
    projpred.parallel_proj_trigger = 0L,
    doFuture.foreach.export = ".export-and-automatic"
  ))
  rec <- local_record_backend(export_only = TRUE)
  for (fun_nm in names(tstsetups)[!is.na(tstsetups)]) {
    tstsetup <- tstsetups[[fun_nm]]
    args_prj_i <- args_prj[[tstsetup]]
    rec$draw_env_nms <- list()
    rec$draw_export_opt <- list()
    p_repr <- suppressWarnings(do.call(project, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]]),
      excl_nonargs(args_prj_i)
    )))
    expect_equal(p_repr, prjs[[tstsetup]], info = tstsetup)
    expect_gt(length(rec$draw_env_nms), 0)
    for (env_nms in rec$draw_env_nms) {
      expect_setequal(env_nms, nms_expected[[fun_nm]])
    }
    for (export_opt in rec$draw_export_opt) {
      expect_identical(export_opt, ".export", info = tstsetup)
    }
    expect_identical(getOption("doFuture.foreach.export"),
                     ".export-and-automatic", info = tstsetup)
  }
})

test_that(paste(
  "the object lists of the loop functions of cv_varsel() and divmin() match",
  "their free variables"
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
  inner_nms <- c(loo_varsel = "one_obs", kfold_varsel = "one_fold",
                 divmin = "one_draw", divmin_augdat = "one_draw")
  for (outer_nm in names(inner_nms)) {
    inner_nm <- inner_nms[[outer_nm]]
    outer_fun <- get(outer_nm, envir = asNamespace("projpred"))
    rhss <- rhs_of(body(outer_fun), inner_nm)
    is_def <- vapply(rhss, function(r) {
      is.call(r) && identical(r[[1]], as.name("function"))
    }, logical(1))
    is_slim <- vapply(rhss, function(r) {
      is.call(r) && identical(r[[1]], as.name("fun_with_objs"))
    }, logical(1))
    expect_identical(sum(is_slim), 1L, info = outer_nm)
    slim_call <- rhss[[which(is_slim)]]
    # The function is defined in a separate assignment or inside the
    # fun_with_objs() call:
    if (is.symbol(slim_call[[2]])) {
      expect_identical(sum(is_def), 1L, info = outer_nm)
      inner_fun <- eval(rhss[[which(is_def)]])
    } else {
      expect_identical(sum(is_def), 0L, info = outer_nm)
      inner_fun <- eval(slim_call[[2]])
    }
    globs <- codetools::findGlobals(inner_fun, merge = FALSE)
    outer_locals <- unique(c(
      names(formals(outer_fun)),
      codetools::findFuncLocals(formals(outer_fun), body(outer_fun))
    ))
    listed <- eval(slim_call[[3]])
    expect_setequal(intersect(c(globs$variables, globs$functions),
                              outer_locals),
                    listed)
    # With `doFuture.foreach.export = ".export"`, the workers find all other
    # free variables only if they are in the projpred namespace, its imports
    # or base (not in the global environment):
    others <- setdiff(c(globs$variables, globs$functions), outer_locals)
    ns <- asNamespace("projpred")
    in_ns <- vapply(others, function(nm) {
      exists(nm, envir = ns, inherits = FALSE) ||
        exists(nm, envir = parent.env(ns), inherits = FALSE) ||
        exists(nm, envir = baseenv(), inherits = FALSE)
    }, logical(1))
    expect_true(all(in_ns), info = paste(outer_nm, ":",
                                         paste(others[!in_ns], collapse = ", ")))
  }
})
