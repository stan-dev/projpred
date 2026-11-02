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
