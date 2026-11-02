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

# Column groups of the model matrix `mm` for the monotonicity constraint. `fr`
# is the model frame of `mm` (its terms with the attribute "predvars", see
# makepredictcall.mo_basis()), and `mm` must have the attributes "assign" and
# "contrasts" (as from model.matrix()). In a term with a mo() call, the
# columns with the same values of the other variables of the term form one
# group, for example all columns of `mo(x)` or of `mo(x):z`, or the columns of
# `mo(x):f` for one column of the factor `f`. As in brms, which gives each such
# column its own monotonic effect, the sum of two groups (for example the
# effect of `x` within one level of `f` in `mo(x) * f`) is not constrained.
# Groups with one column (a predictor with two levels) are left out because
# one coefficient is always monotonic. Returns a named list of the column
# indices in `mm`, ordered by the level of the mo() predictor. The names are
# the column names in which the mo() part is replaced by the mo() call, for
# example `mo(x)` or `mo(x):fb`.
mo_sign_groups <- function(fr, mm) {
  trms <- attr(fr, "terms")
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
  if (is.null(predvars) || is.null(asgn) ||
      !all(var_nms[in_term] %in% names(fr))) {
    stop("Internal error: the model frame or the model matrix does not have ",
         "the expected attributes. Please notify the package maintainer.")
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
    if (is.null(mo_levels)) {
      stop("Internal error: the levels of `", var_nm, "` are missing in the ",
           "terms. Please notify the package maintainer.")
    }
    n_cols <- length(mo_levels) - 1L
    if (n_cols < 2) {
      next
    }
    # model.matrix() builds the columns of a term as the product of the column
    # sets of its variables (in the order of the rows of the "factors"
    # attribute), with the first variable changing fastest. Thus the position
    # of a column gives its mo() level and the column of the other variables:
    vars_before <- which(fcts[, t_idx] > 0 & seq_along(var_nms) < mo_vars)
    n_before <- prod(vapply(vars_before, function(i) {
      n_cols_in_term(fr[[var_nms[i]]], coding = fcts[i, t_idx],
                     contr = attr(mm, "contrasts")[[var_nms[i]]])
    }, integer(1)))
    cols <- which(asgn == t_idx)
    if (length(cols) %% (n_before * n_cols) != 0) {
      stop("Internal error: the columns of term `", colnames(fcts)[t_idx],
           "` do not match its mo() levels. Please notify the package ",
           "maintainer.")
    }
    pos <- seq_along(cols) - 1L
    lvls <- (pos %/% n_before) %% n_cols + 1L
    rest <- pos %% n_before + n_before * (pos %/% (n_before * n_cols))
    grps <- lapply(split(seq_along(cols), rest), function(g) {
      cols[g][order(lvls[g])]
    })
    # The column name of an interaction is the `:`-separated list of the
    # column names of its variables. In the column of the first level, replace
    # the part of the mo() variable by the mo() call:
    mo_part <- paste0(var_nm, as.character(mo_levels[2L]))
    names(grps) <- vapply(grps, function(g) {
      mo_label(cnms[g[1L]], mo_part, var_nm)
    }, character(1), USE.NAMES = FALSE)
    out <- c(out, grps)
  }
  if (length(out) > 0) {
    names(out) <- make.unique(names(out))
  }
  return(out)
}

# The number of columns of the model frame variable `v` in a term of
# model.matrix(), with the coding `coding` from the "factors" attribute of the
# terms (1: contrasts, 2: indicators of all levels) and the contrasts `contr`
# (an element of the "contrasts" attribute of the model matrix, or `NULL`).
n_cols_in_term <- function(v, coding, contr) {
  dat <- list(.v = v)
  if (coding == 2) {
    return(ncol(model.matrix(~ 0 + .v, data = dat)))
  }
  contr_arg <- if (is.null(contr)) NULL else list(.v = contr)
  return(ncol(model.matrix(~ .v, data = dat, contrasts.arg = contr_arg)) - 1L)
}

# The column name `cnm` with the part `mo_part` (bounded by `:` or by the ends
# of `cnm`) replaced by `var_nm`.
mo_label <- function(cnm, mo_part, var_nm) {
  starts <- gregexpr(mo_part, cnm, fixed = TRUE)[[1L]]
  ends <- starts + nchar(mo_part) - 1L
  ok <- starts > 0 &
    (starts == 1L | substr(cnm, starts - 1L, starts - 1L) == ":") &
    (ends == nchar(cnm) | substr(cnm, ends + 1L, ends + 1L) == ":")
  if (!any(ok)) {
    stop("Internal error: the column name `", cnm, "` does not contain `",
         mo_part, "`. Please notify the package maintainer.")
  }
  s <- starts[ok][1L]
  return(paste0(substr(cnm, 1L, s - 1L), var_nm,
                substring(cnm, s + nchar(mo_part))))
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
