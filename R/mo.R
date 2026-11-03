# Monotonic terms ---------------------------------------------------------

# Basis for a monotonic term `mo(x)` (as in brms). For an ordered factor or an
# integer `x` with levels `l_0, ..., l_D` (for an integer `x`,
# `min(x), ..., max(x)`), column `k` (k = 1, ..., D) is the indicator of
# `x >= l_k`. With `j` the position of `x` in the levels (0 for `l_0`), brms
# models the effect of `mo(x)` as b * D * sum(zeta[seq_len(j)]) with `zeta` on
# the simplex. This is a linear combination of these columns, so the
# projection onto the coefficients of the columns reproduces the reference
# model's effect. The submodel fitters constrain the coefficients to a common
# sign (see mo_sign_groups() and option `projpred.mo_monotone`). Argument `id`
# (brms: a shared simplex) has no effect here.
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
  if (is_mo_call(expr)) {
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
# example `mo(x)` or `mo(x):fb`. Attribute "keys" holds the names with the
# coding of the term appended; groups with the same key have the same meaning
# in different models (used for the signs from the full model). Terms with
# more than one mo() call give an error, or are left out if `skip_multi` is
# `TRUE`.
mo_sign_groups <- function(fr, mm, skip_multi = FALSE) {
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
    stop_mo_unsupported(paste0(
      "mo() calls inside other calls (variable ",
      paste0("`", var_nms[has_mo & !is_mo], "`", collapse = ", "), ")"
    ))
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
  keys <- character()
  for (t_idx in seq_len(ncol(fcts))) {
    mo_vars <- which(is_mo & fcts[, t_idx] > 0)
    if (length(mo_vars) == 0) {
      next
    }
    if (length(mo_vars) > 1) {
      if (skip_multi) {
        next
      }
      stop_mo_unsupported(paste0(
        "terms with more than one mo() call (term `", colnames(fcts)[t_idx],
        "`)"
      ))
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
    # The coding of the variables of the term (see above) distinguishes, for
    # example, `fb:mo(x)` in `f * mo(x)` (the difference to the first level
    # of `f`) from `fb:mo(x)` in `f + f:mo(x)` (the effect within level `b`):
    # (Only factors have columns that depend on the coding.)
    vars_term <- which(fcts[, t_idx] > 0)
    is_fac <- vapply(vars_term, function(i) {
      v <- fr[[var_nms[i]]]
      is.factor(v) || is.character(v) || is.logical(v)
    }, logical(1))
    term_coding <- paste(fcts[vars_term[is_fac], t_idx], collapse = "")
    keys <- c(keys, rep(term_coding, length(grps)))
    out <- c(out, grps)
  }
  if (length(out) > 0) {
    names(out) <- make.unique(names(out))
    attr(out, "keys") <- paste0(names(out), "#", keys)
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

# Monotonic projections ---------------------------------------------------

# Are the projections of mo() terms monotonic (global option
# `projpred.mo_monotone`)?
mo_monotone <- function() {
  return(isTRUE(getOption("projpred.mo_monotone", TRUE)))
}

# Does one of the predictor terms `terms_` (character) contain a mo() call?
has_mo_term <- function(terms_) {
  return(any(vapply(terms_, function(trm) {
    has_mo_call(str2lang(trm))
  }, logical(1))))
}

# Adds an environment for the signs from the full model (see mo_full_signs())
# to `p_ref` (output of get_refdist() or get_p_clust()) if the reference model
# has mo() terms. Copies of `p_ref` share this environment, and a new `p_ref`
# (other draws or clusters) gets a new one.
mo_add_cache <- function(p_ref, refmodel) {
  fml <- refmodel$formula
  if (inherits(fml, "formula") && has_mo_call(fml[[length(fml)]])) {
    p_ref$mo_cache <- new.env(parent = emptyenv())
  }
  return(p_ref)
}

# The signs of the sign groups (see mo_sign_groups()) in the unconstrained
# projection of `p_ref` onto the full model (all terms of the reference model
# formula). `args_divmin` are the arguments of the divergence minimizer for a
# submodel; their formula and data are replaced. Returns an S x G matrix (S
# projected draws or clusters, G groups) with the keys of the groups as column
# names. A sign is 1 if the sum of the coefficients of the group is
# nonnegative, otherwise -1. The result is cached in `p_ref$mo_cache`.
mo_full_signs <- function(p_ref, refmodel, args_divmin) {
  cache <- p_ref$mo_cache
  if (!is.null(cache) && !is.null(cache$signs)) {
    return(cache$signs)
  }
  data <- refmodel$fetch_data()
  terms_full <- setdiff(split_formula(refmodel$formula, data = data,
                                      add_lower_terms = FALSE), "1")
  subset <- subset_formula_and_data(formula = refmodel$formula,
                                    terms_ = terms_full, data = data,
                                    y = p_ref$mu)
  fml_full <- flatten_formula(subset$formula)
  if (formula_contains_additive_terms(fml_full)) {
    stop_mo_unsupported("reference models with smooth terms")
  }
  args_divmin$formula <- fml_full
  args_divmin$data <- subset$data
  args_divmin$projpred_mo_signs <- NULL
  outdmin <- tryCatch(
    do.call(refmodel$div_minimizer, args_divmin),
    error = function(e) {
      if (grepl("^Internal error", conditionMessage(e))) {
        stop(e)
      }
      warning("The projection onto the full model (for the signs of the ",
              "monotonicity constraint of mo() terms) failed with the error ",
              "\"", conditionMessage(e), "\". The signs come from the sums ",
              "of the unconstrained coefficients of the submodels instead.",
              call. = FALSE)
      return(NULL)
    }
  )
  if (is.null(outdmin)) {
    signs <- matrix(numeric(), nrow = NCOL(p_ref$mu), ncol = 0)
    if (!is.null(cache)) {
      cache$signs <- signs
    }
    return(signs)
  }
  # The sign groups of the population-level part of the full model:
  fml_fixed <- reformulas::nobars(fml_full)
  fml_rhs <- stats::as.formula(call("~", fml_fixed[[length(fml_fixed)]]),
                               env = environment(fml_full))
  fr <- model.frame(fml_rhs, data = subset$data, drop.unused.levels = TRUE)
  mm <- model.matrix(fml_rhs, data = fr)
  groups <- mo_sign_groups(fr, mm, skip_multi = TRUE)
  cnms <- colnames(mm)
  signs <- vapply(outdmin, function(fit) {
    b <- mo_fit_coefs(fit)
    vapply(groups, function(g) {
      idx <- match(cnms[g], names(b))
      if (anyNA(idx)) {
        stop("Internal error: the coefficients of the full model do not ",
             "match its model matrix. Please notify the package maintainer.")
      }
      mo_sign_of_sum(b[idx])
    }, numeric(1))
  }, numeric(length(groups)))
  # S x G (vapply() gives a G x S array, or a vector for G = 1):
  signs <- matrix(t(matrix(signs, nrow = length(groups))),
                  nrow = length(outdmin),
                  dimnames = list(NULL, attr(groups, "keys")))
  if (!is.null(cache)) {
    cache$signs <- signs
  }
  return(signs)
}

# The population-level coefficients of a submodel fit (a named vector).
# Columns that lme4 dropped (rank deficiency) get the coefficient 0.
mo_fit_coefs <- function(fit) {
  if (inherits(fit, "subfit")) {
    return(stats::setNames(as.vector(fit$beta), rownames(fit$beta)))
  }
  if (inherits(fit, "merMod")) {
    b <- lme4::fixef(fit, add.dropped = TRUE)
    b[is.na(b)] <- 0
    return(b)
  }
  if (inherits(fit, "lme")) {
    return(fit$coefficients$fixed)
  }
  return(stats::coef(fit))
}

# The sign (1 or -1) of the sum of the coefficients `b`.
mo_sign_of_sum <- function(b) {
  if (!all(is.finite(b))) {
    stop("A projection gave non-finite coefficients of a mo() term, so the ",
         "monotonicity constraint cannot be applied.")
  }
  return(if (sum(b) >= 0) 1 else -1)
}

# The sign groups (see mo_sign_groups()) as column indices of the model matrix
# without its intercept column. `fr` and `mm` as in mo_sign_groups().
mo_groups_x <- function(fr, mm) {
  groups <- mo_sign_groups(fr, mm)
  if (length(groups) == 0) {
    return(groups)
  }
  pos_x <- cumsum(colnames(mm) != "(Intercept)")
  keys <- attr(groups, "keys")
  groups <- lapply(groups, function(g) pos_x[g])
  attr(groups, "keys") <- keys
  return(groups)
}

# The scales (standard deviations, as in glm_ridge(), with observation weights
# `weights`) of the columns of `x`, for mo_binding().
mo_col_scale <- function(x, weights = NULL) {
  return(standardization(x, weights = weights %||% rep(1, nrow(x)))$scale)
}

# Is the monotonicity constraint binding for the coefficients `beta` (of the
# columns of `x`), that is, does a sign group have both positive and negative
# coefficients? The tolerance applies on the scale of the standardized
# columns (`scale`, see mo_col_scale()).
mo_binding <- function(beta, x, groups, weights = NULL,
                       scale = mo_col_scale(x, weights), tol = 1e-10) {
  if (length(groups) == 0) {
    return(FALSE)
  }
  b_std <- as.vector(beta) * scale
  return(any(vapply(groups, function(g) {
    any(b_std[g] > tol) && any(b_std[g] < -tol)
  }, logical(1))))
}

# The signs of the sign groups for the constrained fit: the sign from the full
# model (`signs_full`, a named vector with the keys of the groups as names),
# or the sign of the sum of the unconstrained coefficients `beta` if the full
# model has no group with the same key.
mo_cons_signs <- function(beta, groups, signs_full) {
  keys <- attr(groups, "keys")
  beta <- as.vector(beta)
  return(vapply(seq_along(groups), function(i) {
    s <- if (is.null(signs_full)) NA_real_ else signs_full[keys[i]]
    if (is.na(s)) {
      s <- mo_sign_of_sum(beta[groups[[i]]])
    }
    unname(s)
  }, numeric(1)))
}

# Additional arguments of glm_ridge() for the constrained fit, given the
# unconstrained fit `fit` with the model matrix `x` (without intercept
# column). The start is the unconstrained solution (on the scale of the
# standardized columns, as glm_ridge() expects it), unless `dot_args` (the
# other arguments of glm_ridge()) set the start or the standardization.
mo_cons_args <- function(fit, x, groups, signs_full, dot_args) {
  s_g <- mo_cons_signs(fit$beta, groups, signs_full)
  out <- list(cons_cols = unlist(groups, use.names = FALSE),
              cons_signs = rep(s_g, lengths(groups)))
  if (!any(c("beta_init", "beta0_init", "intercept", "normalize") %in%
           names(dot_args))) {
    transf <- standardization(
      x, weights = dot_args$weights %||% rep(1, nrow(x))
    )
    transf$scale[transf$scale == 0] <- 1
    beta <- as.vector(fit$beta)
    out$beta_init <- beta * transf$scale
    out$beta0_init <- fit$beta0 + sum(transf$shift * beta)
  }
  return(out)
}

# The number of mo() calls in the expression `expr`.
n_mo_calls <- function(expr) {
  if (!is.call(expr)) {
    return(0L)
  }
  if (is_mo_call(expr)) {
    return(1L)
  }
  return(sum(vapply(as.list(expr)[-1L], n_mo_calls, integer(1))))
}

# Throws an error that monotonic projections of mo() terms do not support
# `why`.
stop_mo_unsupported <- function(why) {
  stop("Monotonic projections of mo() terms do not support ", why, ". Use ",
       "`options(projpred.mo_monotone = FALSE)` for unconstrained projections.",
       call. = FALSE)
}

# The terms (language objects) of the predictor terms `terms_` (character;
# elements can be sums of terms, as in `search_terms`), without the
# intercept.
split_term_exprs <- function(terms_) {
  split_plus <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L) {
      return(c(split_plus(e[[2L]]), split_plus(e[[3L]])))
    }
    return(list(e))
  }
  exprs <- unlist(lapply(terms_, function(trm) split_plus(str2lang(trm))))
  return(Filter(function(e) !(is.numeric(e) && e %in% c(0, 1)), exprs))
}

# Does the term `expr` contain a mo() call inside another call (other than an
# interaction)?
has_nested_mo <- function(expr) {
  if (!is.call(expr) || is_mo_call(expr)) {
    return(FALSE)
  }
  if (identical(expr[[1L]], as.name(":"))) {
    return(any(vapply(as.list(expr)[-1L], has_nested_mo, logical(1))))
  }
  return(has_mo_call(expr))
}

# Throws an error at the start of a search or a projection for the cases that
# monotonic projections of mo() terms do not support. `terms_` are the
# predictor terms of the possible submodels (elements can be sums of terms, as
# in `search_terms`), and `method` is the search method (`NULL` for
# project()). For a search, the check is conservative: for example, it stops
# for a group term with a non-Gaussian family also if no submodel combines it
# with a mo() term.
check_mo_monotone <- function(refmodel, terms_, method = NULL) {
  if (!mo_monotone()) {
    return(invisible(TRUE))
  }
  exprs <- split_term_exprs(terms_)
  is_mo <- vapply(exprs, has_mo_call, logical(1))
  if (!any(is_mo)) {
    return(invisible(TRUE))
  }
  has_grp <- any(vapply(exprs, function(e) {
    grepl("|", paste(deparse(e), collapse = ""), fixed = TRUE)
  }, logical(1)))
  why <- if (refmodel$family$for_augdat) {
    "the augmented-data projection (the latent projection is supported)"
  } else if (!identical(refmodel$div_minimizer, divmin)) {
    "custom divergence minimizers"
  } else if (formula_contains_additive_terms(refmodel$formula)) {
    "reference models with smooth terms"
  } else if ((!is.null(method) || !has_grp) &&
             !identical(getOption("projpred.glm_fitter",
                                  "fit_glm_ridge_callback"),
                        "fit_glm_ridge_callback")) {
    paste("other submodel fitters than the default one (global option",
          "`projpred.glm_fitter`)")
  } else if (identical(method, "L1")) {
    "L1 search (forward search is supported)"
  } else if (has_grp && !mo_lmm_supported(refmodel$family)) {
    paste("multilevel submodels with a non-Gaussian family (the latent",
          "projection is supported)")
  } else if (any(vapply(exprs[is_mo], n_mo_calls, integer(1)) > 1L)) {
    "terms with more than one mo() call"
  } else if (any(vapply(exprs[is_mo], has_nested_mo, logical(1)))) {
    "mo() calls inside other calls"
  } else {
    NULL
  }
  if (!is.null(why)) {
    stop_mo_unsupported(why)
  }
  return(invisible(TRUE))
}

#' Monotonic effects of `mo()` terms in the parametrization of \pkg{brms}
#'
#' For a projection with `mo()` terms (see section "Monotonic terms" in
#' [init_refmodel()]), `mo_effects()` converts the projected coefficients of
#' each `mo()` term into the parametrization of \pkg{brms}: a coefficient `b`
#' and a simplex `zeta`, so that the effect of level \eqn{j} (counted from 0)
#' is \eqn{b D \sum_{k \le j} \zeta_k}{b * D * sum(zeta[1:j])} with \eqn{D}
#' the number of levels minus 1.
#'
#' @param object An object of class `projection` (returned by [project()],
#'   possibly as elements of a `list`).
#' @param ... Currently ignored.
#'
#' @details For the projected coefficients \eqn{\beta_1, \ldots, \beta_D}{beta_1,
#'   ..., beta_D} of the levels of a `mo()` term, \eqn{b = \sum_k \beta_k / D}{b
#'   = sum(beta) / D} and \eqn{\zeta_k = \beta_k / \sum_l \beta_l}{zeta_k =
#'   beta_k / sum(beta)}. If all coefficients are 0, `zeta` is `NA`. In the
#'   default monotonic projection (global option `projpred.mo_monotone`), all
#'   coefficients of a term have the same sign, so `zeta` is on the simplex.
#'   In an unconstrained projection, elements of `zeta` can be negative or
#'   larger than 1.
#'
#'   As in \pkg{brms}, an interaction of a `mo()` term with a factor has one
#'   monotonic effect per column of the factor (for example `mo(x):fb`), and
#'   an interaction with a numeric variable has one monotonic effect (for
#'   example `mo(x):z`). Terms with only two levels have no simplex and are
#'   left out.
#'
#' @return For a `projection` object, a [posterior::draws_matrix()] with one
#'   row per projected draw (weighted by the weights of the projected draws)
#'   and the columns `bsp_<term>` and `simo_<term>1[k]` (names as in
#'   \pkg{brms}, for example `bsp_mox` and `simo_mox1[1]` for `mo(x)`). For a
#'   `list` of `projection` objects, a `list` of such matrices.
#'
#' @export
mo_effects <- function(object, ...) {
  if (is.list(object) && !inherits(object, "projection") &&
      length(object) > 0 &&
      all(vapply(object, inherits, logical(1), what = "projection"))) {
    return(lapply(object, mo_effects, ...))
  }
  if (!inherits(object, "projection")) {
    stop("`object` must be of class `projection` (or a list of such objects).")
  }
  refm <- object$refmodel
  fml <- make_formula(object$predictor_terms, formula = refm$formula)
  fml_fixed <- reformulas::nobars(fml)
  fml_rhs <- stats::as.formula(call("~", fml_fixed[[length(fml_fixed)]]),
                               env = environment(fml))
  fr <- model.frame(fml_rhs, data = refm$fetch_data(),
                    drop.unused.levels = TRUE)
  mm <- model.matrix(fml_rhs, data = fr)
  groups <- mo_sign_groups(fr, mm, skip_multi = TRUE)
  if (length(groups) == 0) {
    stop("The submodel has no mo() terms with more than two levels.")
  }
  cnms <- colnames(mm)
  coefs <- do.call(rbind, lapply(object$outdmin, function(fit) {
    b <- mo_fit_coefs(fit)
    b[match(cnms, names(b))]
  }))
  out <- do.call(cbind, lapply(seq_along(groups), function(i) {
    g <- groups[[i]]
    beta <- coefs[, g, drop = FALSE]
    if (anyNA(beta)) {
      stop("Internal error: the coefficients of the submodel do not match ",
           "its model matrix. Please notify the package maintainer.")
    }
    beta_sum <- rowSums(beta)
    zeta <- beta / beta_sum
    zeta[beta_sum == 0, ] <- NA
    coef_nm <- brms_rename(names(groups)[i])
    out_i <- cbind(beta_sum / length(g), zeta)
    colnames(out_i) <- c(paste0("bsp_", coef_nm),
                         paste0("simo_", coef_nm, "1[", seq_along(g), "]"))
    out_i
  }))
  rownames(out) <- NULL
  attr(out, "wdraws_prj") <- object$wdraws_prj
  return(mat2drmat(out))
}

# The renaming of brms for parameter names (brms:::rename() with its default
# patterns).
brms_rename <- function(x) {
  pattern <- c(" ", "(", ")", "[", "]", ",", "\"", "'", "?", "+", "-", "*",
               "/", "^", "=", "$")
  replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ", "USD")
  for (i in seq_along(pattern)) {
    x <- gsub(pattern[i], replacement[i], x, fixed = TRUE)
  }
  return(x)
}

# Can the submodel fitters with family `family` apply the monotonicity
# constraint to multilevel submodels? (Only the Gaussian lmer() path.)
mo_lmm_supported <- function(family) {
  return(family$family == "gaussian" && family$link == "identity" &&
           !getOption("projpred.PQL", FALSE) &&
           getOption("projpred.gaussian_not_as_generalized", TRUE))
}

# The sign groups (see mo_sign_groups()) of the population-level part of the
# multilevel formula `formula`, as column indices of `X` (the population-level
# model matrix of lme4, which can lack columns because of rank deficiency).
mo_lmm_groups <- function(formula, data, X) {
  fml_fixed <- reformulas::nobars(formula)
  fr <- model.frame(fml_fixed, data = data, drop.unused.levels = TRUE)
  mm <- model.matrix(fml_fixed, data = fr)
  groups <- mo_sign_groups(fr, mm)
  if (length(groups) == 0) {
    return(groups)
  }
  keys <- attr(groups, "keys")
  groups <- lapply(groups, function(g) {
    idx <- match(colnames(mm)[g], colnames(X))
    idx[!is.na(idx)]
  })
  keep <- lengths(groups) >= 2
  groups <- groups[keep]
  attr(groups, "keys") <- keys[keep]
  return(groups)
}

# Gaussian linear multilevel submodel with the monotonicity constraint. `fit`
# is the unconstrained lmer() fit of `formula`, `lmer_args` are the other
# arguments of that lmer() call, and `signs_full` the signs from the full model
# (see mo_cons_signs()). If a sign group of the unconstrained fit has
# coefficients of both signs, the REML (or ML) criterion of lmer() is profiled
# over the variance parameters `theta`; for given `theta`, the
# population-level coefficients come from a generalized least-squares step
# with the sign bounds. (The REML criterion keeps the terms of the full
# population-level model matrix also for coefficients at a bound.) The
# solution is returned as an lmerMod object (lmer() with fixed `theta`, then
# the constrained coefficients and the matching conditional modes in its
# slots). Otherwise, `fit` is returned.
mo_lmer_cons <- function(fit, formula, lmer_args, signs_full) {
  args_lf <- lmer_args[intersect(names(lmer_args),
                                 methods::formalArgs(lme4::lFormula))]
  lf <- do.call(lme4::lFormula, c(list(formula = formula), args_lf))
  X <- lf$X
  beta_unc <- lme4::fixef(fit)
  if (!identical(names(beta_unc), colnames(X))) {
    stop("Internal error: the coefficients of the multilevel submodel do not ",
         "match its model matrix. Please notify the package maintainer.")
  }
  groups <- mo_lmm_groups(formula, lmer_args$data, X)
  w <- as.vector(stats::model.weights(lf$fr) %||% rep(1, nrow(X)))
  if (!mo_binding(beta_unc, X, groups, weights = w)) {
    return(fit)
  }
  s_g <- mo_cons_signs(beta_unc, groups, signs_full)
  cons <- as.integer(unlist(groups, use.names = FALSE))
  REML <- lme4::isREML(fit)
  n <- nrow(X)
  p <- ncol(X)
  y <- stats::model.response(lf$fr) -
    as.vector(stats::model.offset(lf$fr) %||% 0)
  sw <- sqrt(w)
  rt <- lf$reTrms
  Zt <- rt$Zt %*% Matrix::Diagonal(x = sw)
  # Flip the signs of the constrained columns so that all bounds are >= 0:
  sgn <- rep(1, p)
  sgn[cons] <- rep(s_g, lengths(groups))
  Xf <- sweep(X * sw, 2, sgn, "*")
  yw <- y * sw
  Lambdat <- rt$Lambdat
  # The constrained solution for given `theta`:
  gls <- function(theta) {
    Lambdat@x <- theta[rt$Lind]
    LtZt <- Lambdat %*% Zt
    L <- Matrix::Cholesky(Matrix::tcrossprod(LtZt), LDL = FALSE, Imult = 1)
    # V^{-1} a with V = I + Z Lambda Lambda' Z' (weighted):
    vinv <- function(a) {
      a <- as.matrix(a)
      a - as.matrix(Matrix::crossprod(
        LtZt, Matrix::solve(L, LtZt %*% a, system = "A")
      ))
    }
    qmat <- crossprod(Xf, vinv(Xf))
    qvec <- drop(crossprod(Xf, vinv(yw)))
    b0 <- tryCatch(solve(qmat, qvec), error = function(e) NULL)
    if (is.null(b0)) {
      # (A singular system, for example for very large `theta`.)
      return(list(crit = Inf))
    }
    bls <- bls_c(qmat, qvec, cons, b0)
    beta <- as.vector(bls$b)
    resid <- yw - drop(Xf %*% beta)
    r2 <- sum(resid * drop(vinv(resid)))
    ld_l <- as.numeric(Matrix::determinant(L, sqrt = FALSE)$modulus)
    ld_x <- as.numeric(determinant(qmat, logarithm = TRUE)$modulus)
    crit <- if (REML) {
      ld_l + ld_x + (n - p) * (1 + log(2 * pi * r2 / (n - p)))
    } else {
      ld_l + n * (1 + log(2 * pi * r2 / n))
    }
    if (!is.finite(crit)) {
      crit <- Inf
    }
    nlist(crit, beta = beta * sgn, r2, ld_l, ld_x, LtZt, L,
          bls_converged = bls$converged)
  }
  theta_unc <- lme4::getME(fit, "theta")
  opt <- lme4::nloptwrap(theta_unc, function(th) gls(th)$crit,
                         lower = rt$lower, upper = rep(Inf, length(rt$lower)),
                         control = list(xtol_abs = 1e-8, ftol_abs = 1e-10))
  sol <- gls(opt$par)
  if (!is.finite(sol$crit)) {
    stop("The monotonic projection onto a multilevel submodel failed (the ",
         "criterion is not finite at the optimum).")
  }
  beta <- unname(sol$beta)
  # The lmerMod object with fixed theta:
  ctrl <- lmer_args$control %||% lme4::lmerControl()
  ctrl["optimizer"] <- list(NULL)
  # (Singular fits are possible at the constrained optimum as in lmer(); the
  # internal call should not add a message.)
  ctrl$checkConv$check.conv.singular$action <- "ignore"
  args_lmer <- lmer_args[setdiff(names(lmer_args), c("control", "start"))]
  fit_c <- do.call(lme4::lmer, c(
    list(formula = formula, control = ctrl, start = list(theta = opt$par)),
    args_lmer
  ))
  # The conditional modes for the constrained coefficients:
  Lambdat@x <- opt$par[rt$Lind]
  u <- drop(as.matrix(Matrix::solve(
    sol$L, sol$LtZt %*% ((y - drop(X %*% beta)) * sw), system = "A"
  )))
  b <- drop(as.matrix(Matrix::crossprod(Lambdat, u)))
  gamma <- drop(X %*% beta) + drop(as.matrix(Matrix::crossprod(rt$Zt, b)))
  fit_c@beta <- beta
  fit_c@u <- u
  # The predictor and response modules are reference objects (shared with
  # the deviance function of `fit_c`), so update copies:
  pp <- fit_c@pp$copy()
  pp$setDelu(u - pp$u0)
  pp$setDelb(beta - pp$beta0)
  resp <- fit_c@resp$copy()
  resp$updateMu(gamma)
  fit_c@pp <- pp
  fit_c@resp <- resp
  wrss <- sum(w * (y - gamma)^2)
  cmp <- fit_c@devcomp$cmp
  cmp[["wrss"]] <- wrss
  cmp[["ussq"]] <- sum(u^2)
  cmp[["pwrss"]] <- wrss + sum(u^2)
  cmp[["ldL2"]] <- sol$ld_l
  cmp[["ldRX2"]] <- sol$ld_x
  cmp[["sigmaML"]] <- sqrt(cmp[["pwrss"]] / n)
  cmp[["sigmaREML"]] <- sqrt(cmp[["pwrss"]] / (n - p))
  cmp[[if (REML) "REML" else "dev"]] <- sol$crit - sum(log(w))
  fit_c@devcomp$cmp <- cmp
  optinfo <- fit@optinfo
  optinfo$optimizer <- "nloptwrap"
  optinfo$val <- opt$par
  optinfo$feval <- opt$feval
  optinfo$message <- opt$message
  optinfo$warnings <- list()
  optinfo$derivs <- NULL
  optinfo$conv <- list(opt = opt$conv, lme4 = list())
  if (!sol$bls_converged) {
    optinfo$warnings <- list(paste(
      "maximum number of iterations of the bounded least-squares step",
      "reached"
    ))
  }
  fit_c@optinfo <- optinfo
  return(fit_c)
}
