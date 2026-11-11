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
# sees only the `.export` objects, the iteration arguments and base R. This is
# stricter than doFuture with `doFuture.foreach.export = ".export"`, whose
# workers also attach the `.packages`. With `fail = TRUE`, it stops in the
# first fold loop or projection loop. When `frame` exits, the
# previous backend is registered again (see local_dopar_restore()).
local_record_backend <- function(fail = FALSE, export_only = FALSE,
                                 frame = parent.frame()) {
  rec <- new.env()
  rec$env_nms <- list()
  rec$fit_null <- logical()
  rec$export_opt <- list()
  rec$draw_env_nms <- list()
  rec$draw_export_opt <- list()
  do_record <- function(obj, expr, envir, data) {
    fun_nm <- intersect(c("one_obs", "one_fold"), obj$export)
    if (length(fun_nm) == 1L) {
      rec$env_nms[[fun_nm]] <- ls(environment(get(fun_nm, envir)),
                                  all.names = TRUE)
      rec$fit_null[fun_nm] <- is.null(
        environment(get(fun_nm, envir))$refmodel$fit
      )
      if (fun_nm == "one_fold") {
        # The reference models of the folds (iteration argument `list_cv_k`):
        folds <- get("list_cv", envir)
        rec$fold_fits_null <- all(vapply(folds, function(fold) {
          is.null(fold$refmodel$fit) && !is.null(fold$eta_test)
        }, logical(1)))
      }
      rec$export_opt[fun_nm] <- list(getOption("doFuture.foreach.export"))
      if (fail) stop("Stop of the recording foreach backend.")
    }
    if ("one_draw" %in% obj$export) {
      rec$draw_env_nms <- c(rec$draw_env_nms, list(
        ls(environment(get("one_draw", envir)), all.names = TRUE)
      ))
      rec$draw_export_opt <- c(rec$draw_export_opt,
                               list(getOption("doFuture.foreach.export")))
      if (fail) stop("Stop of the recording foreach backend.")
    }
    if (export_only) {
      envir <- list2env(mget(obj$export, envir = envir, inherits = TRUE),
                        parent = baseenv())
    }
    it <- iterators::iter(obj)
    accumulator <- foreach::makeAccum(it)
    results <- lapply(as.list(it), function(args_i) {
      eval(expr, list2env(args_i, parent = envir))
    })
    accumulator(results, seq_along(results))
    foreach::getResult(it)
  }
  local_dopar_restore(frame)
  foreach::setDoPar(do_record, info = function(data, item) {
    switch(item, name = "record", version = "0", workers = 1L, NULL)
  })
  rec
}

# When `frame` exits, the foreach backend that is registered now is registered
# again. foreach has no exported getter for the registered backend, so this
# reads its internal store (the setup of this file can have registered
# doParallel or doFuture).
local_dopar_restore <- function(frame = parent.frame()) {
  fe_globals <- foreach:::.foreachGlobals
  nms_reg <- c("fun", "data", "info")
  reg_old <- mget(nms_reg, envir = fe_globals, ifnotfound = list(NULL))
  withr::defer({
    if (all(nms_reg %in% ls(fe_globals))) rm(list = nms_reg, envir = fe_globals)
    if (!is.null(reg_old$fun)) {
      foreach::setDoPar(reg_old$fun, data = reg_old$data, info = reg_old$info)
    }
  }, envir = frame)
  invisible()
}

# A GLM with K-fold fits (in the setup, these are the reference models without
# observation weights), preferably with the traditional projection (`NA` if
# there is none). With `pkg`, a GLM of that package (e.g., "rstanarm"):
tstsetup_ref_par <- function(pkg = "[^.]*") {
  tstsetups_ref <- names(cvfitss)[!sapply(cvfitss, inherits, "try-error")]
  tstsetups_ref <- grep(paste0("^", pkg, "\\.glm\\."), tstsetups_ref,
                        value = TRUE)
  c(grep("\\.trad", tstsetups_ref, value = TRUE), tstsetups_ref,
    NA_character_)[1]
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
  # The folds get the reference model without its fit (see
  # refmodel_for_folds()):
  expect_false(is.null(refmods[[tstsetup_ref]]$fit))
  expect_true(isTRUE(unname(rec$fit_null["one_obs"])))
  expect_true(isTRUE(unname(rec$fit_null["one_fold"])))
  expect_true(isTRUE(rec$fold_fits_null))
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
  withr::local_options(list(projpred.parallel_proj_trigger = 0L))
  rec <- local_record_backend(export_only = TRUE)
  for (fun_nm in names(tstsetups)[!is.na(tstsetups)]) {
    tstsetup <- tstsetups[[fun_nm]]
    args_prj_i <- args_prj[[tstsetup]]
    prj_i <- function() {
      suppressWarnings(do.call(project, c(
        list(object = refmods[[args_prj_i$tstsetup_ref]]),
        excl_nonargs(args_prj_i)
      )))
    }
    # A value set by the user, and the option not set (the default):
    for (opt_user in list(".export-and-automatic", NULL)) {
      withr::local_options(list(doFuture.foreach.export = opt_user))
      info_i <- paste(tstsetup, if (is.null(opt_user)) "NULL" else opt_user)
      rec$draw_env_nms <- list()
      rec$draw_export_opt <- list()
      expect_equal(prj_i(), prjs[[tstsetup]], info = info_i)
      expect_gt(length(rec$draw_env_nms), 0)
      for (env_nms in rec$draw_env_nms) {
        expect_setequal(env_nms, nms_expected[[fun_nm]])
      }
      for (export_opt in rec$draw_export_opt) {
        expect_identical(export_opt, ".export", info = info_i)
      }
      expect_identical(getOption("doFuture.foreach.export"), opt_user,
                       info = info_i)
    }
    # An error in the projection loop:
    withr::local_options(list(
      doFuture.foreach.export = ".export-and-automatic"
    ))
    local({
      local_record_backend(fail = TRUE)
      expect_error(prj_i(), "Stop of the recording foreach backend",
                   info = tstsetup)
      expect_identical(getOption("doFuture.foreach.export"),
                       ".export-and-automatic", info = tstsetup)
    })
  }
})

test_that(paste(
  "the object lists of the loop functions of cv_varsel(), divmin() and",
  "divmin_augdat() match their free variables"
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

test_that(paste(
  "refmodel_for_folds() gives a reference model without the fit, with the",
  "objects that mu_fun() and fetch_data() use"
), {
  skip_if_not_installed("codetools")
  skip_if(length(refmods) == 0)
  # The objects of the environment of `fun` that `fun` uses:
  used_objs <- function(fun) {
    globs <- codetools::findGlobals(fun, merge = FALSE)
    intersect(c(globs$variables, globs$functions),
              ls(environment(fun), all.names = TRUE))
  }
  topenv_nm <- function(fun) {
    env <- environment(fun)
    if (is.null(env)) "base" else environmentName(topenv(env))
  }
  # The path from `x` to the first object of a class in `cls` that `x` reaches
  # through list elements, attributes, and environments (except namespaces
  # and the global environment), or `NULL`:
  path_to_class <- function(x, cls) {
    seen <- new.env(parent = emptyenv())
    walk <- function(x, path) {
      if (inherits(x, cls)) return(path)
      if (is.environment(x)) {
        if (isNamespace(x) || identical(x, globalenv()) ||
            identical(x, baseenv()) || identical(x, emptyenv())) {
          return(NULL)
        }
        key <- utils::capture.output(print(x))[1]
        if (exists(key, envir = seen, inherits = FALSE)) return(NULL)
        assign(key, TRUE, envir = seen)
        for (nm in ls(x, all.names = TRUE)) {
          p <- walk(get0(nm, envir = x, inherits = FALSE),
                    paste0(path, " > env$", nm))
          if (!is.null(p)) return(p)
        }
        return(walk(parent.env(x), paste0(path, " > parent.env")))
      }
      if (is.function(x)) {
        p <- walk(environment(x), paste0(path, " > environment"))
        if (!is.null(p)) return(p)
      }
      if (is.list(x)) {
        nms <- if (is.null(names(x))) rep("", length(x)) else names(x)
        for (j in seq_along(x)) {
          p <- walk(x[[j]], paste0(path, " > [[", j, ":", nms[j], "]]"))
          if (!is.null(p)) return(p)
        }
      }
      at <- attributes(x)
      for (anm in names(at)) {
        p <- walk(at[[anm]], paste0(path, " > attr ", anm))
        if (!is.null(p)) return(p)
      }
      NULL
    }
    walk(x, "x")
  }
  fit_cls <- c("stanreg", "brmsfit", "stanfit")
  pkgs_copied <- character()
  for (tstsetup in names(refmods)) {
    refm <- refmods[[tstsetup]]
    refm_f <- projpred:::refmodel_for_folds(refm)
    # For the latent and the augmented-data projection and for family functions
    # from other packages, the reference model stays as it is:
    fam_pkgs <- vapply(Filter(is.function, refm$family), topenv_nm,
                       character(1))
    if (isTRUE(refm$family$for_latent) || isTRUE(refm$family$for_augdat) ||
        !all(fam_pkgs %in% c("base", "stats", "projpred"))) {
      expect_identical(refm_f, refm, info = tstsetup)
      next
    }
    pkgs_copied <- c(pkgs_copied, sub("\\..*", "", tstsetup))
    expect_true("fit" %in% names(refm_f), info = tstsetup)
    expect_null(refm_f$fit, info = tstsetup)
    expect_null(refm_f$cvfits, info = tstsetup)
    expect_error(refm_f$ref_predfun(NULL),
                 "not available in the parallel CV folds", info = tstsetup)
    # The copy changes only the environments of the formula attributes of the
    # data (e.g., "terms" if the data is a model frame):
    no_fmls <- function(d) {
      for (anm in names(attributes(d))) {
        if (inherits(attr(d, anm), "formula")) attr(d, anm) <- NULL
      }
      d
    }
    expect_identical(no_fmls(refm_f$fetch_data()), no_fmls(refm$fetch_data()),
                     info = tstsetup)
    expect_identical(no_fmls(refm_f$fetch_data(obs = 1:3)),
                     no_fmls(refm$fetch_data(obs = 1:3)), info = tstsetup)
    env_f <- environment(refm_f$mu_fun)
    expect_identical(environment(refm_f$fetch_data), env_f, info = tstsetup)
    expect_false(identical(env_f, environment(refm$mu_fun)), info = tstsetup)
    # The new environment holds exactly the objects of the init_refmodel()
    # environment that mu_fun() and fetch_data() use:
    expect_setequal(ls(env_f, all.names = TRUE),
                    union(used_objs(refm$mu_fun), used_objs(refm$fetch_data)))
    # The original reaches the fit, the copy does not:
    if (!is.null(refm$fit)) {
      expect_false(is.null(path_to_class(refm, fit_cls)), info = tstsetup)
      path_f <- path_to_class(refm_f, fit_cls)
      expect_null(path_f, info = paste(tstsetup, ":", path_f))
      expect_lt(length(serialize(refm_f, NULL)),
                length(serialize(refm, NULL)), label = tstsetup)
    }
  }
  # At least one reference model of each package gets the copy:
  pkgs_all <- intersect(sub("\\..*", "", names(refmods)),
                        c("rstanarm", "brms"))
  expect_gt(length(pkgs_all), 0)
  expect_setequal(intersect(pkgs_copied, pkgs_all), pkgs_all)
  # formula_for_folds(): a formula from a function that also holds a fit (as
  # the refit call of rstanarm::kfold()). The new environment holds only the
  # objects that the formula uses and that are not data columns:
  fit_dummy <- structure(list(), class = "stanreg")
  make_fml <- function() {
    x <- fit_dummy
    object <- fit_dummy
    knots <- c(-1, 1)
    y ~ x + splines::ns(z, knots = knots)
  }
  fml <- make_fml()
  fml_f <- projpred:::formula_for_folds(fml, data_nms = c("y", "x", "z"),
                                        env_parent = globalenv())
  expect_false(is.null(path_to_class(fml, fit_cls)))
  expect_null(path_to_class(fml_f, fit_cls))
  expect_identical(ls(environment(fml_f), all.names = TRUE), "knots")
  expect_identical(environment(fml_f)$knots, c(-1, 1))
  expect_identical(parent.env(environment(fml_f)), globalenv())
  expect_identical(deparse(fml_f), deparse(fml))
  expect_identical(class(fml_f), class(fml))
  # Objects of the ancestors of the formula environment (up to `env_parent`)
  # are also kept, the nearest one first. A `.` and an argument without a
  # value are not objects of the formula:
  make_fml_outer <- function(arg) {
    knots <- 3
    offs <- 4
    make_fml_inner <- function() {
      knots <- c(-1, 1)
      y ~ . + splines::ns(z, knots = knots) + offset(offs) + arg
    }
    make_fml_inner()
  }
  fml_o <- make_fml_outer()
  fml_o_f <- projpred:::formula_for_folds(fml_o, data_nms = c("y", "z"),
                                          env_parent = globalenv())
  expect_identical(as.list(environment(fml_o_f), all.names = TRUE)[
    c("knots", "offs")
  ], list(knots = c(-1, 1), offs = 4))
  expect_setequal(ls(environment(fml_o_f), all.names = TRUE),
                  c("knots", "offs"))
  # With the environment of make_fml_outer() as `env_parent`, the new
  # environment holds only the object of make_fml_inner():
  env_o <- parent.env(environment(fml_o))
  fml_o_f <- projpred:::formula_for_folds(fml_o, data_nms = c("y", "z"),
                                          env_parent = env_o)
  expect_identical(ls(environment(fml_o_f), all.names = TRUE), "knots")
  expect_identical(parent.env(environment(fml_o_f)), env_o)
  # A formula with the global environment, a namespace or `env_parent` as its
  # environment stays as it is:
  for (env in list(globalenv(), asNamespace("stats"))) {
    environment(fml) <- env
    expect_identical(
      projpred:::formula_for_folds(fml, data_nms = "y",
                                   env_parent = environment(make_fml())),
      fml
    )
  }
  expect_identical(
    projpred:::formula_for_folds(fml_f, data_nms = "y",
                                 env_parent = environment(fml_f)),
    fml_f
  )
})

test_that(paste(
  "the workers of the parallel CV folds do not load the package of the",
  "reference model fit"
), {
  skip_if_not(run_prll)
  skip_if_not(run_cvvs)
  skip_if_not_installed("doParallel")
  skip_if_not_installed("doRNG")
  tstsetups_ref <- vapply(pkg_nms, tstsetup_ref_par, character(1))
  skip_if(all(is.na(tstsetups_ref)))
  # Workers that stay alive between the loops and that do not get the loaded
  # namespaces of this R session (as forked workers do). With
  # pkgload::load_all() (e.g., in devtools::test()), each worker loads the
  # same projpred source, not an installed projpred:
  cl <- parallel::makePSOCKcluster(2L)
  withr::defer(parallel::stopCluster(cl))
  init_wrk <- function(libs, path) {
    .libPaths(libs)
    if (file.exists(file.path(path, "R", "projpred-package.R"))) {
      # (pkgload is installed if this R session uses it, so it is not in
      # `Suggests:`.)
      load_all <- get("load_all", envir = asNamespace("pkgload"))
      load_all(path, compile = FALSE, helpers = FALSE, quiet = TRUE)
    }
    invisible()
  }
  environment(init_wrk) <- baseenv()
  parallel::clusterCall(cl, init_wrk, libs = .libPaths(),
                        path = getNamespaceInfo("projpred", "path"))
  local_dopar_restore()
  doParallel::registerDoParallel(cl)
  expect_identical(foreach::getDoParWorkers(), 2L)
  # Before each run, the workers get the default RNG kind. doRNG sets
  # L'Ecuyer-CMRG on each worker that runs an iteration, so this shows that
  # both workers ran folds:
  rng_reset <- function() RNGkind("Mersenne-Twister")
  environment(rng_reset) <- baseenv()
  rng_kind <- function() RNGkind()[1]
  environment(rng_kind) <- baseenv()
  check_wrk <- function(label) {
    expect_identical(unlist(parallel::clusterCall(cl, rng_kind)),
                     rep("L'Ecuyer-CMRG", 2), label = label)
    ns <- parallel::clusterCall(cl, loadedNamespaces)
    for (pkg_fit in pkg_nms) {
      expect_false(any(vapply(ns, function(ns_i) pkg_fit %in% ns_i,
                              logical(1))),
                   label = paste(label, pkg_fit))
    }
  }
  # The check comes after each run, so that it names the first CV method
  # that makes the workers load a package of a fit:
  for (tstsetup_ref in tstsetups_ref[!is.na(tstsetups_ref)]) {
    for (cv_meth in c("LOO", "kfold")) {
      parallel::clusterCall(cl, rng_reset)
      expect_s3_class(cvvs_par(tstsetup_ref, cv_meth), "vsel")
      check_wrk(paste(tstsetup_ref, cv_meth))
    }
  }
  # The fits of setup.R have a formula object in their calls. Usually, the
  # call has the formula as an expression. rstanarm::kfold() then evaluates
  # it in the refit call, so the formulas of the fold reference models get an
  # environment that holds the fits:
  tstsetup_ref <- tstsetups_ref[["rstanarm"]]
  skip_if(is.null(tstsetup_ref) || is.na(tstsetup_ref))
  fit_expr <- refmods[[tstsetup_ref]]$fit
  fit_expr$call$formula <- str2lang(paste(deparse(formula(fit_expr)),
                                          collapse = " "))
  refm_expr <- suppressWarnings(get_refmodel(fit_expr))
  cvfits_expr <- try(suppressWarnings(run_cvfun(refm_expr, K = K_tst,
                                                seed = seed3_tst)),
                     silent = TRUE)
  skip_if(inherits(cvfits_expr, "try-error"))
  fml_fold <- projpred:::get_kfold(refm_expr, K = K_tst, cvfits = cvfits_expr,
                                   verbose = FALSE)[[1]]$refmodel$formula
  expect_false(identical(environment(fml_fold),
                         environment(refm_expr$formula)))
  parallel::clusterCall(cl, rng_reset)
  cvvs_expr <- suppressWarnings(cv_varsel(
    refm_expr, cv_method = "kfold", cvfits = cvfits_expr,
    validate_search = TRUE, nclusters = nclusters_tst,
    nclusters_pred = nclusters_pred_tst, nterms_max = nterms_max_tst,
    parallel = TRUE, seed = seed_tst, verbose = FALSE
  ))
  expect_s3_class(cvvs_expr, "vsel")
  check_wrk("rstanarm fit with formula expression, kfold")
})
