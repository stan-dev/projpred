# A customized version of vdiffr::expect_doppelganger() that is needed because
# in 'testthat' versions >= 3.3.0, testthat::expect_snapshot_file() (upon which
# vdiffr::expect_doppelganger() builds) errors if the snapshot file name has
# already been used.
expect_doppelganger_cust <- function(title,
                                     fig,
                                     test_filename = "methods_vsel",
                                     ...) {
  snap_nm_old <- testthat::test_path(
    "_snaps", test_filename, paste0(tolower(gsub("\\.|_", "-", title)), ".svg")
  )
  if (packageVersion("testthat") >= "3.3.0" && file.exists(snap_nm_old)) {
    variant_tmp <- paste0("tmp_doppelganger__",
                          format(Sys.time(), format = "%Y-%m-%d_%H-%M-%S"))
    vdiffr::expect_doppelganger(title, fig, variant = variant_tmp, ...)
    snap_nm_new <- sub("_snaps", file.path("_snaps", variant_tmp), snap_nm_old)
    expect_true(compare_file_text(old = snap_nm_old,
                                  new = snap_nm_new),
                info = title)
    file.remove(snap_nm_new)
    file.remove(testthat::test_path("_snaps", variant_tmp, test_filename))
    file.remove(testthat::test_path("_snaps", variant_tmp))
  } else {
    vdiffr::expect_doppelganger(title = title, fig = fig, ...)
  }
  return(invisible(TRUE))
}

# In edition 3, expect_warning() catches only one warning. This expectation
# catches all warnings whose message matches `regexp` (all warnings for
# `regexp = NULL`) and requires at least one of them. Other warnings are not
# caught. `regexp = NA` requires that there is no warning. Returns the value of
# `object` invisibly.
expect_warning_all <- function(object, regexp = NULL, ..., info = NULL) {
  if (identical(regexp, NA)) {
    expect_warning(val <- object, NA, info = info)
    return(invisible(val))
  }
  n_warn <- 0L
  val <- withCallingHandlers(object, warning = function(w) {
    if (is.null(regexp) || grepl(regexp, conditionMessage(w), ...)) {
      n_warn <<- n_warn + 1L
      invokeRestart("muffleWarning")
    }
  })
  expect(n_warn > 0L,
         if (is.null(regexp)) {
           "Expected a warning."
         } else {
           sprintf("Expected a warning that matches \"%s\".", regexp)
         },
         info = info, trace_env = parent.frame())
  return(invisible(val))
}

# Expects that all elements of the logical `object` are `TRUE`, as
# expect_true(all(object)) does (an empty `object` passes, `NA` fails). A
# failure message gives the number and the positions of the elements that are
# not `TRUE`: names if available, `[i, j, ...]` for a matrix or an array,
# otherwise indices. testthat::expect_all_true() has no `info`, stops for an
# empty `object` and fails for a matrix.
expect_all_true_tst <- function(object, info = NULL, label = NULL) {
  # rlang::enquo() evaluates `!!` at capture time, as in the testthat
  # expectations:
  act <- quasi_label(rlang::enquo(object), label)
  if (!is.logical(act$val)) {
    stop("`object` must be logical, not of class ", class(act$val)[1], ".")
  }
  val <- as.vector(act$val)
  not_true <- which(!(val %in% TRUE))
  if (length(not_true)) {
    nms <- names(act$val)
    dims <- dim(act$val)
    if (length(dims) > 1) {
      pos <- apply(arrayInd(not_true, dims), 1, function(idx) {
        paste0("[", paste(idx, collapse = ", "), "]")
      })
    } else if (!is.null(nms)) {
      pos <- ifelse(is.na(nms[not_true]) | nms[not_true] == "",
                    as.character(not_true), nms[not_true])
    } else {
      pos <- not_true
    }
    msg <- sprintf(
      "Expected all elements of %s to be TRUE. %d of %d are not TRUE, at %s%s.",
      act$lab, length(not_true), length(val),
      paste(head(pos, 10), collapse = ", "),
      if (length(not_true) > 10) ", ..." else ""
    )
  } else {
    msg <- ""
  }
  expect(length(not_true) == 0, msg, info = info, trace_env = parent.frame())
  return(invisible(act$val))
}

# The reference class objects in the slots `resp` and `pp` of an lme4 fit hold
# an external pointer, which differs between two fits with equal content.
# expect_equal() (edition 3) compares the pointers, so this replaces each lme4
# fit by a list of its slots, with the fields of the reference class objects
# except the pointer. `x` is a fit or a list that contains fits at any depth.
# The function does not look into other S4 objects, environments, closures or
# attributes.
drop_mer_ptrs <- function(x) {
  if (isS4(x)) {
    if (!methods::is(x, "merMod")) {
      return(x)
    }
    slot_nms <- methods::slotNames(x)
    out <- lapply(setNames(nm = slot_nms), function(slot_nm) {
      methods::slot(x, slot_nm)
    })
    for (slot_nm in intersect(c("resp", "pp"), slot_nms)) {
      obj <- out[[slot_nm]]
      fld_nms <- setdiff(names(obj$getRefClass()$fields()), "Ptr")
      out[[slot_nm]] <- lapply(setNames(nm = fld_nms), function(fld_nm) {
        obj$field(fld_nm)
      })
    }
    out$class_orig <- as.character(class(x))
    return(structure(out, class = "mer_noptr"))
  }
  if (is.list(x) && !is.data.frame(x)) {
    attrs <- attributes(x)
    x <- lapply(unclass(x), drop_mer_ptrs)
    attributes(x) <- attrs
  }
  return(x)
}
