# Monotonic terms ---------------------------------------------------------

# Basis for a monotonic term `mo(x)` (as in brms). For an ordered factor or an
# integer `x` with levels `l_0, ..., l_D` (for an integer `x`,
# `min(x), ..., max(x)`), column `k` (k = 1, ..., D) is the indicator of
# `x >= l_k`. With `j` the position of `x` in the levels (0 for `l_0`), brms
# models the effect of `mo(x)` as b * D * sum(zeta[seq_len(j)]) with `zeta` on
# the simplex. This is a linear combination of these columns, so the
# projection onto the unconstrained coefficients of the columns reproduces the
# reference model's effect. The submodels do not constrain the coefficients
# to a common sign. Argument `id` (brms: a shared simplex) has no effect here.
# Argument `mo_levels` gives the levels of the training data, for predictions
# for new data (see makepredictcall.mo_basis()).
mo <- function(x, id = NA, mo_levels = NULL) {
  if (is.null(mo_levels)) {
    if (is.factor(x)) {
      if (!is.ordered(x)) {
        stop("In mo() terms, the predictor must be an ordered factor or an ",
             "integer vector.")
      }
      mo_levels <- levels(x)
    } else {
      check_mo_int(x)
      if (all(is.na(x))) {
        stop("In mo() terms, the predictor must have non-missing values.")
      }
      mo_levels <- seq.int(as.integer(min(x, na.rm = TRUE)),
                           as.integer(max(x, na.rm = TRUE)))
    }
  }
  if (is.character(mo_levels)) {
    # An ordered factor in the training data:
    if (!is.factor(x) && !is.character(x)) {
      stop("In mo() terms, the predictor is an ordered factor in the ",
           "training data, so it must be a factor or a character vector in ",
           "new data.")
    }
    x_int <- match(as.character(x), mo_levels) - 1L
    if (anyNA(x_int[!is.na(x)])) {
      stop("In mo() terms, the predictor has levels which do not exist in ",
           "the training data.")
    }
  } else {
    # An integer in the training data:
    if (is.factor(x) || is.character(x)) {
      stop("In mo() terms, the predictor is an integer in the training data, ",
           "so it must be an integer in new data.")
    }
    check_mo_int(x)
    x_int <- as.integer(x) - mo_levels[1L]
    if (any(x_int < 0L | x_int > length(mo_levels) - 1L, na.rm = TRUE)) {
      stop("In mo() terms, the predictor has values outside the range of ",
           "the training data (in K-fold CV, the training data of a fold). ",
           "An ordered factor with fixed levels avoids this in K-fold CV.")
    }
  }
  n_cols <- length(mo_levels) - 1L
  if (n_cols < 1L) {
    stop("In mo() terms, the predictor must have at least two levels.")
  }
  out <- outer(x_int, seq_len(n_cols), ">=") * 1
  colnames(out) <- as.character(mo_levels[-1L])
  attr(out, "mo_levels") <- mo_levels
  class(out) <- c("mo_basis", "matrix")
  return(out)
}

# Check that `x` is a valid integer predictor of a mo() term.
check_mo_int <- function(x) {
  x_obs <- x[!is.na(x)]
  if (!is.numeric(x) || any(!is.finite(x_obs) | x_obs != round(x_obs))) {
    stop("In mo() terms, the predictor must be an ordered factor or an ",
         "integer vector.")
  }
  return(invisible(TRUE))
}

# Stores the levels of the training data in the "predvars" attribute of the
# terms. The function in the stored call is projpred's mo() with its namespace,
# because some fitters (mgcv::gam(), for example) set the environment of the
# terms to the global environment, where `mo` is brms::mo() or does not exist.
#' @noRd
#' @export
makepredictcall.mo_basis <- function(var, call) {
  mo_ns <- call(":::", as.name("projpred"), as.name("mo"))
  if (!identical(call[[1L]], quote(mo)) && !identical(call[[1L]], mo_ns)) {
    return(NextMethod())
  }
  call[[1L]] <- mo_ns
  call$mo_levels <- attr(var, "mo_levels")
  return(call)
}

# Does the expression `expr` contain a call of mo()?
has_mo_call <- function(expr) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  if (identical(expr[[1L]], quote(mo))) {
    return(TRUE)
  }
  return(any(vapply(as.list(expr)[-1L], has_mo_call, logical(1))))
}

# Is the expression `expr` a call of mo() (also with a namespace)?
is_mo_call <- function(expr) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  fun <- expr[[1L]]
  if (identical(fun, quote(mo))) {
    return(TRUE)
  }
  return(is.call(fun) && length(fun) == 3L &&
           (identical(fun[[1L]], as.name("::")) ||
              identical(fun[[1L]], as.name(":::"))) &&
           identical(fun[[3L]], quote(mo)))
}

# Column groups of the model matrix `mm` for the monotonicity constraint. `mm`
# must have the attribute "assign" (as from model.matrix()), and `trms` must be
# the terms of its model frame (with the attribute "predvars", see
# makepredictcall.mo_basis()). In a term with a mo() call, the columns with the
# same values of the other variables of the term form one group, for example
# all columns of `mo(x)` or of `mo(x):z`, or the columns of `mo(x):f` for one
# column of the factor `f`. (brms also gives each such column of a mo() term
# its own monotonic effect.) Groups with one column (a predictor with two
# levels) are left out because one coefficient is always monotonic. Returns a
# named list of the column indices in `mm`, ordered by the level of the mo()
# predictor. The names are the column names in which the mo() part is
# replaced by the mo() call, for example `mo(x)` or `mo(x):fb`.
mo_sign_groups <- function(trms, mm) {
  fcts <- attr(trms, "factors")
  if (length(fcts) == 0) {
    return(list())
  }
  var_nms <- rownames(fcts)
  in_term <- rowSums(fcts) > 0
  has_mo <- in_term & vapply(var_nms, function(nm) {
    has_mo_call(str2lang(nm))
  }, logical(1))
  if (!any(has_mo)) {
    return(list())
  }
  is_mo <- has_mo & vapply(var_nms, function(nm) {
    is_mo_call(str2lang(nm))
  }, logical(1))
  if (any(has_mo & !is_mo)) {
    stop("mo() calls inside other calls are not supported (variable ",
         paste0("`", var_nms[has_mo & !is_mo], "`", collapse = ", "), ").")
  }
  predvars <- attr(trms, "predvars")
  asgn <- attr(mm, "assign")
  if (is.null(predvars) || is.null(asgn)) {
    stop("Internal error: the attributes \"predvars\" or \"assign\" are ",
         "missing. Please notify the package maintainer.")
  }
  cnms <- colnames(mm)
  out <- list()
  for (t_idx in seq_len(ncol(fcts))) {
    mo_vars <- which(is_mo & fcts[, t_idx] > 0)
    if (length(mo_vars) == 0) {
      next
    }
    if (length(mo_vars) > 1) {
      stop("Terms with more than one mo() call are not supported with ",
           "monotonic projections (term `", colnames(fcts)[t_idx], "`). Use ",
           "`options(projpred.mo_monotone = FALSE)` for unconstrained ",
           "projections.")
    }
    var_nm <- var_nms[mo_vars]
    # The elements of the "predvars" call are the response and the variables,
    # in the order of the rows of the "factors" attribute:
    mo_levels <- eval(predvars[[mo_vars + 1L]]$mo_levels, baseenv())
    n_cols <- length(mo_levels) - 1L
    cols <- which(asgn == t_idx)
    if (is.null(mo_levels) || length(cols) %% n_cols != 0) {
      stop("Internal error: the columns of term `", colnames(fcts)[t_idx],
           "` do not match its mo() levels. Please notify the package ",
           "maintainer.")
    }
    if (n_cols < 2) {
      next
    }
    # The column name of an interaction is the `:`-separated list of the
    # column names of its variables. Find the part of the mo() variable (the
    # longest match, because the levels of an ordered factor may contain `:`)
    # and replace it by the mo() call:
    mo_parts <- paste0(var_nm, as.character(mo_levels[-1L]))
    lvl_lbl <- lapply(cnms[cols], function(cnm) {
      lbls <- vapply(mo_parts, function(part) {
        if (cnm == part) {
          return(var_nm)
        }
        part_c <- paste0(part, ":")
        if (startsWith(cnm, part_c)) {
          return(paste0(var_nm, ":", substring(cnm, nchar(part_c) + 1L)))
        }
        c_part <- paste0(":", part)
        if (endsWith(cnm, c_part)) {
          return(paste0(substr(cnm, 1L, nchar(cnm) - nchar(c_part)), ":",
                        var_nm))
        }
        pos <- regexpr(paste0(":", part, ":"), cnm, fixed = TRUE)
        if (pos > 0) {
          return(paste0(substr(cnm, 1L, pos), var_nm,
                        substring(cnm, pos + nchar(part) + 1L)))
        }
        return(NA_character_)
      }, character(1), USE.NAMES = FALSE)
      if (all(is.na(lbls))) {
        return(list(lvl = NA_integer_, lbl = NA_character_))
      }
      lvl <- which.max(ifelse(is.na(lbls), -1L, nchar(mo_parts)))
      return(list(lvl = lvl, lbl = lbls[lvl]))
    })
    lvls <- vapply(lvl_lbl, "[[", integer(1), "lvl")
    lbls <- vapply(lvl_lbl, "[[", character(1), "lbl")
    grps <- split(seq_along(cols), factor(lbls, levels = unique(lbls)))
    if (anyNA(lvls) || !all(vapply(grps, function(g) {
      identical(sort(lvls[g]), seq_len(n_cols))
    }, logical(1)))) {
      stop("Internal error: the columns of term `", colnames(fcts)[t_idx],
           "` could not be assigned to its mo() levels. Please notify the ",
           "package maintainer.")
    }
    out <- c(out, lapply(grps, function(g) cols[g][order(lvls[g])]))
  }
  return(out)
}

# If `formula` contains mo() terms, give it an environment in which `mo` is
# projpred's mo() (and not, for example, brms::mo()). The parent of this
# environment is the original one, so the other objects of the original
# environment stay visible.
mo_formula_env <- function(formula) {
  if (!has_mo_call(formula[[length(formula)]])) {
    return(formula)
  }
  env_old <- environment(formula) %||% globalenv()
  if (identical(get0("mo", envir = env_old, inherits = FALSE), mo)) {
    return(formula)
  }
  env_new <- new.env(parent = env_old)
  assign("mo", mo, envir = env_new)
  environment(formula) <- env_new
  return(formula)
}

# Throw an error for mo() terms in the cases that projpred does not support:
# in group-level terms and in additive multilevel models (the predictions of
# such submodels for new data would not use the levels of the training data).
# `fml_extractions` is the output of extract_terms_response().
check_mo_terms <- function(fml_extractions) {
  has_mo <- function(trms) {
    any(vapply(trms, function(trm) has_mo_call(str2lang(trm)), logical(1)))
  }
  grp_trms <- fml_extractions$group_terms
  if (length(grp_trms) > 0 && has_mo(grp_trms)) {
    stop("mo() terms are not supported in group-level terms.")
  }
  other_trms <- c(fml_extractions$individual_terms,
                  fml_extractions$interaction_terms)
  if (length(grp_trms) > 0 && length(fml_extractions$additive_terms) > 0 &&
      length(other_trms) > 0 && has_mo(other_trms)) {
    stop("mo() terms are not supported in additive multilevel models (models ",
         "with group-level and smooth terms).")
  }
  return(invisible(TRUE))
}
