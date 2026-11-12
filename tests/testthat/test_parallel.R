# Setup -------------------------------------------------------------------

# The plans of the future package for the parallel tests: a PSOCK cluster
# ("psock") and, if installed, the workers of the future.mirai package:
par_backends_tst <- if (run_prll &&
                        requireNamespace("future", quietly = TRUE)) {
  c("psock",
    if (requireNamespace("future.mirai", quietly = TRUE)) "future.mirai")
} else {
  character()
}

# The numbers of workers of the parallel tests. R CMD check permits only 2
# processes (see `_R_CHECK_LIMIT_CORES_`):
n_wrks_tst <- if (tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_")) %in%
                  c("", "false")) {
  2:3
} else {
  2L
}

# Sets up `n` workers of the plan `backend` (see `par_backends_tst`) until
# `frame` exits. The workers stay alive between the parallel loops and do not
# get the loaded namespaces of this R session (as forked workers do). With
# pkgload::load_all() (e.g., in devtools::test()), each worker loads the same
# projpred source, not an installed projpred. Returns a function that gives
# the loaded namespaces of each worker.
local_par_workers <- function(backend, n = 2L, frame = parent.frame()) {
  init_wrk <- function(libs, path) {
    .libPaths(libs)
    if (file.exists(file.path(path, "R", "projpred-package.R")) &&
        !"projpred" %in% loadedNamespaces()) {
      # (pkgload is installed if this R session uses it, so it is not in
      # `Suggests:`.)
      load_all <- get("load_all", envir = asNamespace("pkgload"))
      load_all(path, compile = FALSE, helpers = FALSE, quiet = TRUE)
    }
    invisible()
  }
  environment(init_wrk) <- baseenv()
  args_init <- list(libs = .libPaths(),
                    path = getNamespaceInfo("projpred", "path"))
  if (backend == "psock") {
    cl <- parallel::makePSOCKcluster(n)
    withr::defer(parallel::stopCluster(cl), envir = frame)
    parallel::clusterCall(cl, init_wrk, libs = args_init$libs,
                          path = args_init$path)
    plan_old <- future::plan(future::cluster, workers = cl)
    withr::defer(future::plan(plan_old), envir = frame)
    ns_fun <- function() parallel::clusterCall(cl, loadedNamespaces)
  } else if (backend == "future.mirai") {
    plan_old <- future::plan(future.mirai::mirai_multisession, workers = n)
    withr::defer(future::plan(plan_old), envir = frame)
    # future cannot address a single worker, so run `fun` in more futures
    # than there are workers, with a pause that makes all workers take
    # one. Gives one result per worker:
    on_workers <- function(fun, args = list()) {
      one <- function(fun, args) {
        Sys.sleep(0.3)
        list(pid = Sys.getpid(), out = do.call(fun, args))
      }
      environment(one) <- baseenv()
      res <- future::value(lapply(seq_len(4L * n), function(i) {
        future::future(one(fun, args),
                       globals = list(one = one, fun = fun, args = args),
                       seed = NULL, envir = baseenv())
      }))
      res <- res[!duplicated(vapply(res, "[[", integer(1), "pid"))]
      skip_if_not(length(res) == n, "The futures did not reach all workers.")
      lapply(res, "[[", "out")
    }
    on_workers(init_wrk, args_init)
    ns_fun <- function() on_workers(loadedNamespaces)
  } else {
    stop("Unrecognized `backend`.")
  }
  invisible(ns_fun)
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

cvvs_par <- function(tstsetup_ref, cv_meth, parallel = TRUE) {
  # Use suppressWarnings() because of Pareto-k warnings of the small data:
  suppressWarnings(cv_varsel(
    refmods[[tstsetup_ref]], cv_method = cv_meth,
    cvfits = if (cv_meth == "kfold") cvfitss[[tstsetup_ref]],
    validate_search = TRUE, nclusters = nclusters_tst,
    nclusters_pred = nclusters_pred_tst, nterms_max = nterms_max_tst,
    parallel = parallel, seed = seed_tst, verbose = FALSE
  ))
}

# par_lapply() ------------------------------------------------------------

test_that(paste(
  "par_lapply() gives the results of lapply(), with RNG streams that do not",
  "depend on the backend or the number of workers"
), {
  skip_if_not(run_prll)
  skip_if(length(par_backends_tst) == 0)
  # The functions go to the workers with their environments, so do not use
  # the environment of this test:
  fun_rng <- function(i, a) c(i + a, stats::runif(1))
  environment(fun_rng) <- baseenv()
  fun_err <- function(i) if (i == 3L) stop("fold error") else i
  environment(fun_err) <- baseenv()
  fun_opt <- function(i) {
    list(getOption("projpred.tst_opt"),
         getOption("projpred.parallel_proj_trigger"),
         getOption("projpred.parallel_cv"))
  }
  environment(fun_opt) <- baseenv()
  fun_glob <- function(i) i + projpred_tst_glob
  environment(fun_glob) <- globalenv()
  # A function whose environment has the global environment as its parent
  # (as a function that a script makes in local() or in another function):
  fun_fact <- local({
    a <- 1
    function(i) i + a + projpred_tst_glob - 1
  }, envir = new.env(parent = globalenv()))
  fun_get <- function(i) i + get(paste0("projpred_tst_", "glob"))
  environment(fun_get) <- globalenv()
  fun_has_glob <- function(i) {
    exists("projpred_tst_glob", envir = globalenv(), inherits = FALSE)
  }
  environment(fun_has_glob) <- baseenv()
  assign("projpred_tst_glob", 5, envir = globalenv())
  withr::defer(rm(list = "projpred_tst_glob", envir = globalenv()))
  rng_kind_old <- RNGkind()
  # The global options that par_lapply() reads:
  opts_tst <- function() {
    opts <- options()
    opts[grepl("^projpred\\.", names(opts)) |
           names(opts) %in% c("contrasts", "na.action")]
  }
  res_first <- NULL
  for (backend in par_backends_tst) {
    for (n_wrk in n_wrks_tst) {
      local({
        local_par_workers(backend, n_wrk)
        info_i <- paste(backend, n_wrk)
        expect_identical(projpred:::par_backend(), "future", info = info_i)
        res <- withr::with_seed(seed_tst, projpred:::par_lapply(
          1:7, fun_rng, a = 10, .backend = "future", .rng = TRUE
        ))
        expect_identical(vapply(res, "[", numeric(1), 1), 11:17 + 0,
                         info = info_i)
        expect_identical(length(unique(vapply(res, "[", numeric(1), 2))), 7L,
                         info = info_i)
        if (is.null(res_first)) {
          res_first <<- res
        } else {
          expect_identical(res, res_first, info = info_i)
        }
        expect_identical(RNGkind(), rng_kind_old, info = info_i)
        # More chunks than workers:
        withr::with_options(
          list(projpred.parallel_chunks_per_worker = 2L),
          res_chunks <- withr::with_seed(seed_tst, projpred:::par_lapply(
            1:7, fun_rng, a = 10, .backend = "future", .rng = TRUE
          ))
        )
        expect_identical(res_chunks, res_first, info = info_i)
        # No elements, and fewer elements than workers:
        expect_identical(
          projpred:::par_lapply(list(), fun_err, .backend = "future"), list(),
          info = info_i
        )
        expect_identical(
          projpred:::par_lapply(1L, fun_err, .backend = "future"), list(1L),
          info = info_i
        )
        # An error in a task stops the call, and the next call works:
        opts_before <- opts_tst()
        expect_error(projpred:::par_lapply(1:4, fun_err, .backend = "future"),
                     "fold error", info = info_i)
        expect_identical(opts_tst(), opts_before, info = info_i)
        expect_identical(
          projpred:::par_lapply(1:2, fun_err, .backend = "future"),
          list(1L, 2L), info = info_i
        )
        # The workers get the projpred options, except those of the
        # parallelization, only for the call:
        withr::with_options(
          list(projpred.tst_opt = "a", projpred.parallel_proj_trigger = 0L,
               projpred.parallel_cv = TRUE),
          res_opt <- projpred:::par_lapply(1:2, fun_opt, .backend = "future")
        )
        expect_identical(res_opt, rep(list(list("a", NULL, NULL)), 2),
                         info = info_i)
        expect_identical(
          projpred:::par_lapply(1:2, fun_opt, .backend = "future"),
          rep(list(list(NULL, NULL, NULL)), 2), info = info_i
        )
        # The workers get the objects of the global environment that a
        # function of the global environment uses (found from its code), but
        # only for the call:
        if (requireNamespace("globals", quietly = TRUE)) {
          expect_identical(
            projpred:::par_lapply(1:2, fun_glob, .backend = "future"),
            list(6, 7), info = info_i
          )
          expect_identical(
            projpred:::par_lapply(1:2, fun_fact, .backend = "future"),
            list(6, 7), info = info_i
          )
        }
        expect_identical(
          projpred:::par_lapply(1:2, fun_has_glob, .backend = "future"),
          list(FALSE, FALSE), info = info_i
        )
        # A name that the code creates at run time needs global option
        # `projpred.export_to_workers`:
        expect_error(projpred:::par_lapply(1:2, fun_get, .backend = "future"),
                     "projpred_tst_glob", info = info_i)
        withr::with_options(
          list(projpred.export_to_workers = "projpred_tst_glob"),
          res_get <- projpred:::par_lapply(1:2, fun_get, .backend = "future")
        )
        expect_identical(res_get, list(6, 7), info = info_i)
        # With the option set to `FALSE`, the workers get no such objects:
        withr::with_options(
          list(projpred.export_to_workers = FALSE),
          expect_error(
            projpred:::par_lapply(1:2, fun_glob, .backend = "future"),
            "projpred_tst_glob", info = info_i
          )
        )
      })
    }
  }
  # Without a backend, the elements run in this R session, with the same
  # streams, and the RNG kind and the RNG state after the call are those of
  # a call with workers (one draw for the streams):
  res_local <- withr::with_seed(seed_tst, {
    res <- projpred:::par_lapply(1:7, fun_rng, a = 10, .backend = NULL,
                                 .rng = TRUE)
    list(res = res, next_draw = stats::runif(1))
  })
  expect_identical(res_local$res, res_first)
  expect_identical(RNGkind(), rng_kind_old)
  next_draw <- withr::with_seed(seed_tst, {
    sample.int(.Machine$integer.max, 6L)
    stats::runif(1)
  })
  expect_identical(res_local$next_draw, next_draw)
  # Without a backend, a function can change an object of the global
  # environment (the objects are not exported and restored then):
  assign("projpred_tst_count", 0, envir = globalenv())
  withr::defer(rm(list = "projpred_tst_count", envir = globalenv()))
  fun_count <- function(i) {
    projpred_tst_count <<- projpred_tst_count + 1
    i
  }
  environment(fun_count) <- globalenv()
  expect_identical(projpred:::par_lapply(1:3, fun_count, .backend = NULL),
                   list(1L, 2L, 3L))
  expect_identical(get("projpred_tst_count", envir = globalenv()), 3)
  # The sequential loop with streams gives the same results:
  res_seq <- withr::with_seed(seed_tst, {
    res <- projpred:::seq_lapply_rng(7L, function(i) c(i + 10, stats::runif(1)))
    list(res = res, next_draw = stats::runif(1))
  })
  expect_identical(res_seq$res, res_first)
  expect_identical(res_seq$next_draw, next_draw)
})

test_that(paste(
  "workers that do not load the projpred of this R session stop the parallel",
  "run"
), {
  skip_if_not(run_prll)
  skip_if_not("psock" %in% par_backends_tst)
  # Only with pkgload::load_all(): a new worker then loads the installed
  # projpred, or none.
  skip_if_not(file.exists(file.path(getNamespaceInfo("projpred", "path"), "R",
                                    "projpred-package.R")))
  cl <- parallel::makePSOCKcluster(1L)
  withr::defer(parallel::stopCluster(cl))
  plan_old <- future::plan(future::cluster, workers = cl)
  withr::defer(future::plan(plan_old))
  # No result of an earlier check of other workers:
  par_state <- get("par_state", envir = asNamespace("projpred"))
  par_state$key_checked <- NULL
  expect_error(projpred:::par_lapply(1:2, identity, .backend = "future"),
               "The parallel workers (load|cannot load) projpred")
})

test_that("a parallel run without workers warns and runs in this R session", {
  skip_if_not(run_cvvs)
  skip_if_not_installed("future")
  skip_if_not(inherits(future::plan(), "sequential"))
  tstsetup_ref <- tstsetup_ref_par()
  skip_if(is.na(tstsetup_ref))
  expect_null(suppressWarnings(projpred:::par_backend()))
  warns <- character()
  cvvs <- withCallingHandlers(
    cv_varsel(
      refmods[[tstsetup_ref]], cv_method = "LOO", validate_search = TRUE,
      nclusters = nclusters_tst, nclusters_pred = nclusters_pred_tst,
      nterms_max = nterms_max_tst, parallel = TRUE, seed = seed_tst,
      verbose = FALSE
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("there are no parallel workers", warns)))
  cvvs_seq <- cvvs_par(tstsetup_ref, "LOO", parallel = FALSE)
  expect_identical(ranking(cvvs)[["foldwise"]],
                   ranking(cvvs_seq)[["foldwise"]])
  expect_equal(cvvs$summaries, cvvs_seq$summaries)
  # The loops over the projected draws warn only once per R session:
  par_state <- get("par_state", envir = asNamespace("projpred"))
  warned_old <- par_state$warned
  withr::defer(par_state$warned <- warned_old)
  par_state$warned <- NULL
  expect_warning(projpred:::par_backend(warn_once = TRUE),
                 "there are no parallel workers")
  expect_no_warning(projpred:::par_backend(warn_once = TRUE))
})

test_that(paste(
  "the workers get the objects of the global environment that a function of",
  "a custom reference model uses"
), {
  skip_if_not(run_prll)
  skip_if_not(run_cvvs)
  skip_if(length(par_backends_tst) == 0)
  skip_if_not_installed("globals")
  tstsetup_ref <- tstsetup_ref_par()
  skip_if(is.na(tstsetup_ref))
  # A divergence minimizer as a user defines it in a script: a function of
  # the global environment that uses another function and (through it) an
  # object of the global environment:
  assign("projpred_tst_flag", TRUE, envir = globalenv())
  check_tst <- function() stopifnot(projpred_tst_flag)
  environment(check_tst) <- globalenv()
  assign("projpred_tst_check", check_tst, envir = globalenv())
  withr::defer(rm(list = c("projpred_tst_flag", "projpred_tst_check"),
                  envir = globalenv()))
  divmin_tst <- function(...) {
    projpred_tst_check()
    get("divmin", envir = asNamespace("projpred"))(...)
  }
  environment(divmin_tst) <- globalenv()
  refm <- refmods[[tstsetup_ref]]
  refm$div_minimizer <- divmin_tst
  expect_true(any(vapply(projpred:::global_closures(list(refm)), identical,
                         logical(1), divmin_tst)))
  expect_true(all(c("projpred_tst_flag", "projpred_tst_check") %in%
                    names(projpred:::par_globals(list(refm)))))
  cvvs_custom <- function(parallel) {
    suppressWarnings(cv_varsel(
      refm, cv_method = "LOO", validate_search = TRUE,
      nclusters = nclusters_tst, nclusters_pred = nclusters_pred_tst,
      nterms_max = nterms_max_tst, parallel = parallel, seed = seed_tst,
      verbose = FALSE
    ))
  }
  cvvs_seq <- cvvs_custom(FALSE)
  for (backend in par_backends_tst) {
    local({
      local_par_workers(backend)
      cvvs <- cvvs_custom(TRUE)
      expect_identical(ranking(cvvs)[["foldwise"]],
                       ranking(cvvs_seq)[["foldwise"]], info = backend)
      expect_equal(cvvs$summaries, cvvs_seq$summaries, info = backend)
      # Without the export, the function fails on the workers:
      withr::with_options(
        list(projpred.export_to_workers = FALSE),
        expect_error(cvvs_custom(TRUE), "projpred_tst_check", info = backend)
      )
    })
  }
})

# project() ---------------------------------------------------------------

test_that("project() in parallel gives the same results as sequentially", {
  skip_if_not(run_prll)
  skip_if_not(run_prj)
  skip_if(length(par_backends_tst) == 0)
  withr::local_options(list(projpred.parallel_proj_trigger = 0L))
  tstsetups <- grep("\\.glm\\.", names(prjs), value = TRUE)
  # Most GLM projections fit all draws in one call, so add a multilevel
  # submodel (loop in divmin()) and an augmented-data projection (loop in
  # divmin_augdat()):
  first_prj <- function(cond) {
    names(Filter(cond, args_prj[names(prjs)]))[1]
  }
  tstsetups_loop <- c(
    first_prj(function(args_prj_i) {
      args_prj_i$mod_nm == "glmm" && args_prj_i$prj_nm == "trad" &&
        any(grepl("|", args_prj_i$predictor_terms, fixed = TRUE))
    }),
    first_prj(function(args_prj_i) {
      args_prj_i$prj_nm == "augdat" && args_prj_i$fam_nm == "brnll"
    })
  )
  tstsetups_loop <- tstsetups_loop[!is.na(tstsetups_loop)]
  tstsetups <- union(tstsetups, tstsetups_loop)
  # Counts the calls of par_lapply() with workers:
  n_par <- 0L
  par_lapply_orig <- get("par_lapply", envir = asNamespace("projpred"))
  testthat::local_mocked_bindings(
    par_lapply = function(..., .backend) {
      if (!is.null(.backend)) n_par <<- n_par + 1L
      par_lapply_orig(..., .backend = .backend)
    },
    .package = "projpred"
  )
  for (backend in par_backends_tst) {
    local({
      local_par_workers(backend)
      for (tstsetup in tstsetups) {
        args_prj_i <- args_prj[[tstsetup]]
        n_par <<- 0L
        # Use suppressWarnings() because of convergence warnings of the
        # multilevel submodel:
        p_repr <- suppressWarnings(do.call(project, c(
          list(object = refmods[[args_prj_i$tstsetup_ref]]),
          excl_nonargs(args_prj_i)
        )))
        info_i <- paste(backend, tstsetup)
        if (tstsetup %in% tstsetups_loop) {
          expect_gt(n_par, 0L, label = info_i)
          # A fitted submodel object from a worker (e.g., of lme4) can hold
          # class definitions that are not identical to those of this R
          # session, so compare the projected draws and the other components:
          draws_of <- function(prj) {
            suppressWarnings(suppressMessages(as.matrix(prj)))
          }
          expect_equal(draws_of(p_repr), draws_of(prjs[[tstsetup]]),
                       info = info_i)
          nms_other <- setdiff(names(p_repr), "outdmin")
          expect_equal(p_repr[nms_other], prjs[[tstsetup]][nms_other],
                       info = info_i)
        } else {
          expect_equal(p_repr, prjs[[tstsetup]], info = info_i)
        }
      }
    })
  }
})

# varsel() ----------------------------------------------------------------

test_that("varsel() in parallel gives the same results as sequentially", {
  skip_if_not(run_prll)
  skip_if_not(run_vs)
  skip_if(length(par_backends_tst) == 0)
  withr::local_options(list(projpred.parallel_proj_trigger = 0L))
  tstsetups <- grep("\\.glm\\.", names(vss), value = TRUE)
  for (backend in par_backends_tst) {
    local({
      local_par_workers(backend)
      for (tstsetup in tstsetups) {
        args_vs_i <- args_vs[[tstsetup]]
        vs_repr <- do.call(varsel, c(
          list(object = refmods[[args_vs_i$tstsetup_ref]]),
          excl_nonargs(args_vs_i)
        ))
        expect_equal(vs_repr, vss[[tstsetup]], info = paste(backend, tstsetup))
      }
    })
  }
})

# cv_varsel() -------------------------------------------------------------

test_that("cv_varsel() in parallel gives the same results as sequentially", {
  skip_if_not(run_prll)
  skip_if_not(run_cvvs)
  skip_if(length(par_backends_tst) == 0)
  withr::local_options(list(projpred.parallel_proj_trigger = 0L))
  tstsetups <- grep("\\.glm\\.", names(cvvss), value = TRUE)
  for (backend in par_backends_tst) {
    local({
      local_par_workers(backend)
      for (tstsetup in tstsetups) {
        args_cvvs_i <- args_cvvs[[tstsetup]]
        # Use suppressWarnings() because test_that() somehow redirects
        # stderr() and so throws warnings that projpred wants to capture
        # internally:
        cvvs_repr <- suppressWarnings(do.call(cv_varsel, c(
          list(object = refmods[[args_cvvs_i$tstsetup_ref]],
               cvfits = if (identical(args_cvvs_i$cv_method, "kfold")) {
                 cvfitss[[args_cvvs_i$tstsetup_ref]]
               } else {
                 refmods[[args_cvvs_i$tstsetup_ref]]$cvfits # should be `NULL`
               }),
          excl_nonargs(args_cvvs_i)
        )))
        expect_equal(cvvs_repr, cvvss[[tstsetup]],
                     info = paste(backend, tstsetup))
      }
    })
  }
})

# The four tests below do not need `run_prll`: the first two run the parallel
# loops in this R session (see local_record_par()), and the other two only
# read code and objects.

# Until `frame` exits, the parallel loops run in this R session, as one task
# of par_lapply(): par_backend() gives a backend without a check for workers,
# and par_lapply() records its calls in `calls` of the returned environment.
# For each call, it records the objects in the enclosing environment of the
# loop function (these objects go to the workers), whether the reference
# model there has no fit, whether the folds among the elements have reference
# models without fits and the predictions `eta_test`, and `.rng`.
local_record_par <- function(frame = parent.frame()) {
  rec <- new.env()
  rec$calls <- list()
  par_lapply_rec <- function(.X, .FUN, ..., .backend, .rng = FALSE,
                             .progress = FALSE) {
    dots <- list(...)
    fun <- if (is.null(dots$.fun)) .FUN else dots$.fun
    rec$calls <- c(rec$calls, list(list(
      env_nms = ls(environment(fun), all.names = TRUE),
      fit_null = is.null(environment(fun)$refmodel$fit),
      folds_slim = all(vapply(.X, function(x) {
        !is.list(x) || is.null(x$fold) ||
          (is.null(x$fold$refmodel$fit) && !is.null(x$fold$eta_test))
      }, logical(1))),
      rng = .rng
    )))
    par_chunk <- get("par_chunk", envir = asNamespace("projpred"))
    par_chunk(list(xs = .X, seeds = NULL),
              list(fun = .FUN, args = dots, opts = list(), globs = list()))
  }
  testthat::local_mocked_bindings(
    par_backend = function(...) "future", par_lapply = par_lapply_rec,
    .package = "projpred", .env = frame
  )
  rec
}

test_that(paste(
  "the fold function of the parallel cv_varsel() holds only the objects that",
  "it uses, and the reference models without their fits"
), {
  skip_if_not(run_cvvs)
  tstsetup_ref <- tstsetup_ref_par()
  skip_if(is.na(tstsetup_ref))
  rec <- local_record_par()
  for (cv_meth in c("LOO", "kfold")) {
    cvvs <- cvvs_par(tstsetup_ref, cv_meth)
    expect_true(inherits(cvvs, "vsel"), info = cv_meth)
  }
  # The calls of the fold loops, LOO first:
  calls_folds <- Filter(function(cl) "refmodel" %in% cl$env_nms, rec$calls)
  expect_length(calls_folds, 2)
  expect_false(any(c("psisloo", "loglik_forPSIS", "refdist_pred") %in%
                     calls_folds[[1]]$env_nms))
  expect_false("list_cv" %in% calls_folds[[2]]$env_nms)
  # The folds get the reference models without their fits (see
  # refmodel_for_folds()), and each fold gets an RNG stream:
  expect_false(is.null(refmods[[tstsetup_ref]]$fit))
  for (cl in calls_folds) {
    expect_true(cl$fit_null)
    expect_true(cl$folds_slim)
    expect_true(cl$rng)
  }
})

test_that(paste(
  "the parallel projection of single draws holds only the objects that it",
  "uses and gives the same results as the sequential one"
), {
  skip_if_not(run_prj)
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
  rec <- local_record_par()
  for (fun_nm in names(tstsetups)[!is.na(tstsetups)]) {
    tstsetup <- tstsetups[[fun_nm]]
    args_prj_i <- args_prj[[tstsetup]]
    rec$calls <- list()
    prj_i <- suppressWarnings(do.call(project, c(
      list(object = refmods[[args_prj_i$tstsetup_ref]]),
      excl_nonargs(args_prj_i)
    )))
    expect_equal(prj_i, prjs[[tstsetup]], info = tstsetup)
    expect_gt(length(rec$calls), 0)
    for (cl in rec$calls) {
      expect_setequal(cl$env_nms, nms_expected[[fun_nm]])
      expect_false(cl$rng)
    }
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
    # The workers find all other free variables only if they are in the
    # projpred namespace, its imports or base (not in the global
    # environment):
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
  "the parallel cv_varsel() gives the same results for all backends and",
  "numbers of workers, and its workers do not load the package of the",
  "reference model fit"
), {
  skip_if_not(run_prll)
  skip_if_not(run_cvvs)
  skip_if(length(par_backends_tst) == 0)
  tstsetups_ref <- vapply(pkg_nms, tstsetup_ref_par, character(1))
  skip_if(all(is.na(tstsetups_ref)))
  tstsetups_ref <- tstsetups_ref[!is.na(tstsetups_ref)]
  expect_same_cvvs <- function(cvvs1, cvvs2, label) {
    expect_identical(ranking(cvvs1)[["foldwise"]],
                     ranking(cvvs2)[["foldwise"]], label = label)
    expect_equal(cvvs1$summaries, cvvs2$summaries, label = label)
  }
  expect_no_fit_pkg <- function(ns, label) {
    for (pkg_fit in pkg_nms) {
      expect_false(any(vapply(ns, function(ns_i) pkg_fit %in% ns_i,
                              logical(1))),
                   label = paste(label, pkg_fit))
    }
  }
  # The first parallel result of each setting:
  cvvss_first <- list()
  for (backend in par_backends_tst) {
    for (n_wrk in n_wrks_tst) {
      local({
        ns_fun <- local_par_workers(backend, n_wrk)
        # The check comes after each run, so that it names the first CV
        # method that makes the workers load a package of a fit:
        for (tstsetup_ref in tstsetups_ref) {
          for (cv_meth in c("LOO", "kfold")) {
            setting <- paste(tstsetup_ref, cv_meth)
            label <- paste(backend, n_wrk, setting)
            cvvs <- cvvs_par(tstsetup_ref, cv_meth)
            expect_s3_class(cvvs, "vsel")
            expect_no_fit_pkg(ns_fun(), label)
            if (is.null(cvvss_first[[setting]])) {
              cvvss_first[[setting]] <<- cvvs
            } else {
              expect_same_cvvs(cvvs, cvvss_first[[setting]], label)
            }
          }
        }
      })
    }
  }
  # Without workers, `parallel = TRUE` gives the results of the workers:
  if (inherits(future::plan(), "sequential")) {
    for (tstsetup_ref in tstsetups_ref) {
      for (cv_meth in c("LOO", "kfold")) {
        setting <- paste(tstsetup_ref, cv_meth)
        expect_same_cvvs(cvvs_par(tstsetup_ref, cv_meth),
                         cvvss_first[[setting]], paste("no workers", setting))
      }
    }
  }
  # The sequential folds get the same RNG streams, so the sequential results
  # are the parallel ones:
  for (tstsetup_ref in tstsetups_ref) {
    for (cv_meth in c("LOO", "kfold")) {
      setting <- paste(tstsetup_ref, cv_meth)
      expect_same_cvvs(cvvs_par(tstsetup_ref, cv_meth, parallel = FALSE),
                       cvvss_first[[setting]], paste("sequential", setting))
    }
  }
  # The fits of setup.R have a formula object in their calls. Usually, the
  # call has the formula as an expression. rstanarm::kfold() then evaluates
  # it in the refit call, so the formulas of the fold reference models get an
  # environment that holds the fits:
  skip_if_not("rstanarm" %in% names(tstsetups_ref))
  tstsetup_ref <- tstsetups_ref[["rstanarm"]]
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
  for (backend in par_backends_tst) {
    local({
      ns_fun <- local_par_workers(backend)
      cvvs_expr <- suppressWarnings(cv_varsel(
        refm_expr, cv_method = "kfold", cvfits = cvfits_expr,
        validate_search = TRUE, nclusters = nclusters_tst,
        nclusters_pred = nclusters_pred_tst, nterms_max = nterms_max_tst,
        parallel = TRUE, seed = seed_tst, verbose = FALSE
      ))
      expect_s3_class(cvvs_expr, "vsel")
      expect_no_fit_pkg(
        ns_fun(), paste(backend, "rstanarm fit with formula expression, kfold")
      )
    })
  }
})
