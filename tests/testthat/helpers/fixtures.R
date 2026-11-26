# Fixtures --------------------------------------------------------------------

# The getters below (get_fit(), get_ref(), get_vs(), ...) return the objects
# that the tests share. The argument is the name of a test setup, that is, a
# name of the argument list of the object type (`args_ref` for get_ref(),
# `args_vs` for get_vs(), and so on). A getter builds the object at its first
# call and keeps it for the later calls. Each build has its own seed and a
# fixed context (see with_fixture_context_tst()), so an object does not depend
# on the order of the calls or on the test that calls the getter first, and the
# global random state does not change.
#
# A build that fails with a known upstream error (see is_known_setup_err())
# drops its test setup. A build that needs a dropped test setup drops its own
# test setup too. A getter call for a dropped test setup skips the rest of the
# test. A build that fails with another error (also in a build that it needs)
# keeps that error: each getter call for that test setup stops with it.

fixtures_cache <- new.env(parent = emptyenv())

# Evaluate `code` in the context of a build: the `projpred.*` options that
# setup.R sets (`opts_fixtures_tst`; other `projpred.*` options are unset) and
# the collation and language "C" (as in a test).
with_fixture_context_tst <- function(code) {
  nms_now <- names(options())
  nms_unset <- setdiff(nms_now[startsWith(nms_now, "projpred.")],
                       names(opts_fixtures_tst))
  opts <- c(opts_fixtures_tst,
            setNames(vector("list", length(nms_unset)), nms_unset))
  withr::with_options(opts, withr::with_collate(
    "C", withr::with_language("C", code)
  ))
}

# Get the object of type `type` for test setup `tstsetup`. `args` is the
# argument list of the type (its names are the valid test setups; its element
# `mod_nm` is used for is_known_setup_err()). `build` is a function of
# `tstsetup` that builds the object. `kfold` is passed to
# is_known_setup_err(). With `drop_known = FALSE`, a known upstream error is
# kept as an error too.
get_fixture_tst <- function(type, tstsetup, args, build, kfold = FALSE,
                            drop_known = TRUE) {
  stopifnot(is.character(tstsetup), length(tstsetup) == 1L, !is.na(tstsetup))
  if (!tstsetup %in% names(args)) {
    stop("Unknown test setup `", tstsetup, "` for object type \"", type,
         "\".")
  }
  key <- paste0(type, ":", tstsetup)
  if (!exists(key, envir = fixtures_cache, inherits = FALSE)) {
    obj <- tryCatch(
      with_fixture_context_tst(withr::with_seed(seed_tst, build(tstsetup))),
      skip_fixture_tst = function(cond) {
        structure(list(reason = paste0("needs a dropped test setup (",
                                       cond$reason, ")")),
                  class = "dropped_fixture_tst")
      },
      failed_fixture_tst_error = function(e) {
        structure(list(error = e), class = "failed_fixture_tst")
      },
      error = function(e) {
        msg <- conditionMessage(e)
        if (drop_known &&
            is_known_setup_err(msg, mod_nm = args[[tstsetup]]$mod_nm,
                               kfold = kfold)) {
          structure(list(reason = paste0("known upstream error: ", msg)),
                    class = "dropped_fixture_tst")
        } else {
          structure(list(error = e), class = "failed_fixture_tst")
        }
      }
    )
    assign(key, obj, envir = fixtures_cache)
  }
  obj <- get(key, envir = fixtures_cache, inherits = FALSE)
  if (inherits(obj, "dropped_fixture_tst")) {
    # A skip condition (as testthat::skip() makes it) with its own class, so
    # that a build can tell it from other skips:
    reason <- paste0(type, " `", tstsetup, "` dropped, ", obj$reason)
    stop(structure(
      list(message = paste0("Reason: ", reason), reason = reason),
      class = c("skip_fixture_tst", "skip", "condition")
    ))
  }
  if (inherits(obj, "failed_fixture_tst")) {
    # An error with its own class, so that a build that needs this object
    # keeps it as an error (and does not test it with is_known_setup_err()):
    stop(structure(
      list(message = paste0("The build of ", type, " `", tstsetup,
                            "` failed: ", conditionMessage(obj$error)),
           call = NULL),
      class = c("failed_fixture_tst_error", "error", "condition")
    ))
  }
  return(obj)
}

# Is the test setup `tstsetup` of object type `type` dropped? This builds the
# object if necessary, but does not skip.
is_dropped_tst <- function(type, tstsetup) {
  getter <- get(paste0("get_", type), mode = "function")
  tryCatch({
    getter(tstsetup)
    FALSE
  }, skip_fixture_tst = function(cond) TRUE)
}

# The names of the dropped test setups among `tstsetups` of object type
# `type`:
dropped_tst <- function(type, tstsetups) {
  Filter(function(tstsetup) is_dropped_tst(type, tstsetup), tstsetups)
}

# The first `n` test setups among `tstsetups` of object type `type` that are
# not dropped. The builds stop after `n` of them:
kept_tst <- function(type, tstsetups, n = Inf) {
  out <- character()
  for (tstsetup in tstsetups) {
    if (length(out) >= n) break
    if (!is_dropped_tst(type, tstsetup)) out <- c(out, tstsetup)
  }
  return(out)
}

## Getters --------------------------------------------------------------------

# A fit is never dropped (an error stops the tests that need it):
get_fit <- function(tstsetup) {
  get_fixture_tst("fit", tstsetup, args_fit, build_fit_tst,
                  drop_known = FALSE)
}

get_ref <- function(tstsetup) {
  get_fixture_tst("ref", tstsetup, args_ref, function(tstsetup) {
    args_ref_i <- args_ref[[tstsetup]]
    do.call(get_refmodel, c(
      list(object = get_fit(args_ref_i$tstsetup_fit)),
      excl_nonargs(args_ref_i)
    ))
  })
}

# The test setups of get_cvfits() are `tstsetups_ref_kfold` (names of
# `args_ref`):
get_cvfits <- function(tstsetup) {
  get_fixture_tst(
    "cvfits", tstsetup, args_ref[tstsetups_ref_kfold],
    function(tstsetup) {
      run_cvfun(object = get_ref(tstsetup), K = K_tst, seed = seed3_tst)
    },
    kfold = TRUE
  )
}

get_vs <- function(tstsetup) {
  get_fixture_tst("vs", tstsetup, args_vs, function(tstsetup) {
    args_vs_i <- args_vs[[tstsetup]]
    do.call(varsel, c(
      list(object = get_ref(args_vs_i$tstsetup_ref)),
      excl_nonargs(args_vs_i)
    ))
  })
}

get_cvvs <- function(tstsetup) {
  get_fixture_tst("cvvs", tstsetup, args_cvvs, function(tstsetup) {
    args_cvvs_i <- args_cvvs[[tstsetup]]
    refmod_i <- get_ref(args_cvvs_i$tstsetup_ref)
    if (identical(args_cvvs_i$cv_method, "kfold")) {
      cvfits_i <- get_cvfits(args_cvvs_i$tstsetup_ref)
    } else {
      cvfits_i <- refmod_i$cvfits # should be `NULL`
    }
    do.call(cv_varsel, c(
      list(object = refmod_i, cvfits = cvfits_i),
      excl_nonargs(args_cvvs_i)
    ))
  })
}

# Projections from the `refmodel` objects (get_prj()) and from the `vsel`
# objects (get_prj_vs(), get_prj_cvvs()):
get_prj <- function(tstsetup) {
  get_fixture_tst("prj", tstsetup, args_prj, function(tstsetup) {
    args_prj_i <- args_prj[[tstsetup]]
    do.call(project, c(
      list(object = get_ref(args_prj_i$tstsetup_ref)),
      excl_nonargs(args_prj_i)
    ))
  })
}

get_prj_vs <- function(tstsetup) {
  get_fixture_tst("prj_vs", tstsetup, args_prj_vs, function(tstsetup) {
    args_prj_vs_i <- args_prj_vs[[tstsetup]]
    do.call(project, c(
      list(object = get_vs(args_prj_vs_i$tstsetup_vsel)),
      excl_nonargs(args_prj_vs_i)
    ))
  })
}

get_prj_cvvs <- function(tstsetup) {
  get_fixture_tst("prj_cvvs", tstsetup, args_prj_cvvs, function(tstsetup) {
    args_prj_cvvs_i <- args_prj_cvvs[[tstsetup]]
    do.call(project, c(
      list(object = get_cvvs(args_prj_cvvs_i$tstsetup_vsel)),
      excl_nonargs(args_prj_cvvs_i)
    ))
  })
}

# Predictions from the projections. The test setups are those of the
# projections. `type` is "pl", "pe" or "pp", with the suffix of the projection
# type ("", "_vs" or "_cvvs"):
get_pred_tst <- function(type, tstsetup) {
  prj_type <- sub("^p[lep]", "prj", type)
  args <- switch(prj_type,
                 "prj" = args_prj,
                 "prj_vs" = args_prj_vs,
                 "prj_cvvs" = args_prj_cvvs,
                 stop("Unknown `prj_type`."))
  get_prj_type <- get(paste0("get_", prj_type), mode = "function")
  get_fixture_tst(type, tstsetup, args, function(tstsetup) {
    prj <- get_prj_type(tstsetup)
    switch(sub("_.*", "", type),
           "pl" = proj_linpred(prj, allow_nonconst_wdraws_prj = TRUE,
                               .seed = seed2_tst),
           "pe" = proj_epred(prj, allow_nonconst_wdraws_prj = TRUE,
                             .seed = seed2_tst),
           "pp" = proj_predict(prj, .seed = seed2_tst),
           stop("Unknown `type`."))
  })
}

get_pl <- function(tstsetup) get_pred_tst("pl", tstsetup)
get_pe <- function(tstsetup) get_pred_tst("pe", tstsetup)
get_pp <- function(tstsetup) get_pred_tst("pp", tstsetup)
get_pl_vs <- function(tstsetup) get_pred_tst("pl_vs", tstsetup)
get_pe_vs <- function(tstsetup) get_pred_tst("pe_vs", tstsetup)
get_pp_vs <- function(tstsetup) get_pred_tst("pp_vs", tstsetup)
get_pl_cvvs <- function(tstsetup) get_pred_tst("pl_cvvs", tstsetup)
get_pe_cvvs <- function(tstsetup) get_pred_tst("pe_cvvs", tstsetup)
get_pp_cvvs <- function(tstsetup) get_pred_tst("pp_cvvs", tstsetup)
