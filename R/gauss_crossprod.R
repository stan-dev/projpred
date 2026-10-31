# Gaussian submodels from cross-products --------------------------------------
#
# For a Gaussian submodel with the identity link, glm_ridge() solves a weighted
# ridge least-squares problem: the columns of the model matrix are centered and
# scaled with the observation weights, and the intercept is not penalized. Let
# Xc be the model matrix of all terms, centered at the weighted column means,
# W the observation weights, M the fitted values of the projected reference
# draws, sd the weighted column standard deviations (as in .weighted_sd()), and
# Lambda = regul * diag(sd^2). For the columns J of a submodel and draw s, the
# coefficients solve
#   (G[J, J] + Lambda[J, J]) beta_s = C[J, s],  G = Xc' W Xc,  C = Xc' W M,
# and the intercept is ybar_s - m[J]' beta_s, with ybar_s the weighted mean of
# M[, s] and m the weighted column means. With G[J, J] + Lambda[J, J] = R' R and
# y_s = R^-T C[J, s], let q_s = ||y_s||^2 + beta_s' Lambda[J, J] beta_s. Then
#   sum(W * (M[, s] - fit_s)^2) = sum(W * (M[, s] - ybar_s)^2) - q_s,
#   sum(W * (-2 * M[, s] * fit_s + fit_s^2)) = -ybar_s^2 * sum(W) - q_s.
# These give the projected dispersion (dis_gauss()) and the cross-entropy
# (ce_gauss()) without the n-dimensional fitted values.

# The condition for the cross-product computations: Gaussian family with the
# identity link (traditional or latent projection), the default divergence
# minimizer and submodel fitter, zero offsets, no argument for glm_ridge() other
# than `thresh` (not used for this family) in `divmin_args`, and a positive
# `regul`. Returns `regul` or `NULL`.
gauss_crossprod_regul <- function(refmodel, divmin_args) {
  fam <- refmodel$family
  glm_fitter <- getOption("projpred.glm_fitter", "fit_glm_ridge_callback")
  if (fam$family != "gaussian" || fam$link != "identity" ||
      isTRUE(fam$for_augdat) ||
      !identical(refmodel$div_minimizer, divmin) ||
      !identical(glm_fitter, "fit_glm_ridge_callback") ||
      any(refmodel$offset != 0)) {
    return(NULL)
  }
  # fit_glm_ridge_callback() uses `regul` and passes only the formal arguments
  # of glm_ridge() on:
  args_glm_ridge <- intersect(names(divmin_args), methods::formalArgs(glm_ridge))
  if (!all(args_glm_ridge %in% "thresh")) {
    return(NULL)
  }
  # A `regul` of zero or below zero can make the system singular. The draw-wise
  # path then gives an approximate solution:
  regul <- divmin_args$regul %||% 1e-4
  if (!isTRUE(regul > 0)) {
    return(NULL)
  }
  return(regul)
}

# The model matrix of the terms in `term_strings` (a string can consist of
# several terms, for example "x1 + x2 + x1:x2") for the training data of
# `refmodel`, and its cross-products. Returns `NULL` for group-level or smooth
# terms, a factor in an interaction term (its dummy coding depends on the other
# terms of a submodel), missing values, or a zero-variance column.
gauss_crossprod_design <- function(refmodel, term_strings, regul) {
  terms_all <- setdiff(unique(split_plus(term_strings)), "1")
  if (length(terms_all) == 0) {
    return(NULL)
  }
  trms <- extract_terms_response(make_formula(terms_all))
  if (length(trms$group_terms) > 0 || length(trms$additive_terms) > 0) {
    return(NULL)
  }
  # The environment of the reference model formula, as in the draw-wise
  # projections (for functions that the formula calls):
  tt <- terms(as.formula(paste("~", paste(terms_all, collapse = " + ")),
                         env = environment(refmodel$formula)))
  data <- refmodel$fetch_data()
  mf <- model.frame(tt, data = data, drop.unused.levels = TRUE)
  if (nrow(mf) != nrow(data)) {
    return(NULL)
  }
  fac <- attr(tt, "factors")
  is_categ <- vapply(rownames(fac), function(v) {
    is.factor(mf[[v]]) || is.character(mf[[v]]) || is.logical(mf[[v]])
  }, logical(1))
  if (any(attr(tt, "order") > 1 &
          colSums(fac[is_categ, , drop = FALSE]) > 0)) {
    return(NULL)
  }
  x <- model.matrix(tt, data = mf)
  assign_x <- attr(x, "assign")
  keep <- assign_x > 0
  x <- x[, keep, drop = FALSE]
  assign_x <- assign_x[keep]
  cols_by_term <- split(seq_along(assign_x),
                        norm_term_labels(attr(tt, "term.labels"))[assign_x])
  wobs <- refmodel$wobs
  sum_w <- sum(wobs)
  n <- nrow(x)
  x_means <- colSums(x * wobs) / sum_w
  x_c <- sweep(x, 2, x_means)
  # The weighted centered model matrix (for cross-products with the draws):
  x_cw <- x_c * wobs
  g_mat <- crossprod(x_c, x_cw)
  sd_x <- sqrt(n / (n - 1) * diag(g_mat) / sum_w)
  if (n < 2 || any(!is.finite(sd_x)) || any(sd_x == 0)) {
    return(NULL)
  }
  lambda_x <- regul * sd_x^2
  a_full <- g_mat
  diag(a_full) <- diag(a_full) + lambda_x
  # For new data, the terms of the model frame keep the data-dependent bases
  # of the training data (attribute "predvars", for example of poly()):
  return(nlist(tt_pred = attr(mf, "terms"), xlev = .getXlevels(tt, mf), keep,
               x_c, x_cw, x_means, cols_by_term, g_mat, a_full, lambda_x,
               wobs, sum_w))
}

# The term labels of a term string (see norm_term_labels()).
labels_gauss <- function(term_string) {
  if (identical(term_string, "1")) {
    return(character(0))
  }
  cached_by_strings("labels_gauss", term_string, function() {
    unique(norm_term_labels(split_plus(term_string)))
  })
}

# The performance evaluation of perf_eval() along the predictor ranking for
# Gaussian submodels, from cross-products (see the top of this file): for each
# size in `nterms`, the projected dispersion, the cross-entropy and the
# predictions for the test data, without the draw-wise submodel fits. Returns
# the list `out_by_size` of perf_eval() (elements `sub_summary`, `ce` and
# `dis_sub`), or `NULL` if a condition does not hold.
perf_eval_gauss <- function(search_path, nterms, refmodel, p_ref,
                            refmodel_fulldata, indices_test, newdata_test,
                            offset_test, wobs_test, y_test, y_oscale_test,
                            divmin_args) {
  regul <- gauss_crossprod_regul(refmodel, divmin_args)
  ranking <- search_path$predictor_ranking
  if (is.null(regul) || length(nterms) == 0 || any(nterms > length(ranking)) ||
      any(offset_test != 0)) {
    return(NULL)
  }
  des <- gauss_crossprod_design(refmodel, c("1", ranking), regul)
  if (is.null(des)) {
    return(NULL)
  }
  data_test <- newdata_test %||% refmodel_fulldata$fetch_data(obs = indices_test)
  x_test <- tryCatch({
    mf_test <- model.frame(des$tt_pred, data = data_test, xlev = des$xlev,
                           na.action = na.pass)
    model.matrix(des$tt_pred, data = mf_test)[, des$keep, drop = FALSE]
  }, error = function(e) NULL)
  if (is.null(x_test) || nrow(x_test) != nrow(data_test) || anyNA(x_test)) {
    return(NULL)
  }
  x_test_c <- sweep(x_test, 2, des$x_means)

  # Quantities of the projected reference draws:
  wobs <- des$wobs
  sum_w <- des$sum_w
  # Weighted means and sums of squares of the projected draws, in compiled code
  # without temporary observations x draws matrices. The weights scale the
  # (small) model matrix, not the draws.
  if (isTRUE(p_ref$moments_only)) {
    # One draw per cluster, without the observations x draws matrices of
    # `p_ref` (see get_p_clust()): the center of cluster j is
    # `ws[j] * draws_offs[, idx[j]]`, and its predictive variance is
    # `draws_predvar[j]` for all observations.
    draws <- p_ref$draws_offs
    idx <- p_ref$idx
    ws <- p_ref$ws
    mom <- wcol_moments_c(draws, idx - 1L, ws, wobs, sum_w)
    if (!all(is.finite(mom$ss))) {
      # A non-finite value in a projected draw: the full `p_ref` (with
      # predictive variances that depend on the observation) is necessary.
      return(NULL)
    }
    sum_wvar <- p_ref$draws_predvar * sum_w
    if (identical(idx, seq_len(ncol(draws)))) {
      c_mat <- crossprod(des$x_cw, draws)
    } else {
      # Thinned draws: the product only for the projected draws.
      c_mat <- crossprod(des$x_cw, draws[, idx, drop = FALSE])
    }
    c_mat <- c_mat * rep(ws, each = nrow(c_mat))
  } else {
    mu_ref <- unclass(p_ref$mu_offs)
    mom <- wcol_moments_c(mu_ref, seq_len(ncol(mu_ref)) - 1L,
                          rep(1, ncol(mu_ref)), wobs, sum_w)
    sum_wvar <- drop(crossprod(wobs, unclass(p_ref$var)))
    c_mat <- crossprod(des$x_cw, mu_ref)
  }
  ybar <- mom$mean
  sum_wr2 <- mom$ss
  nprj <- length(ybar)

  # Columns of the submodel of each size (the submodels are nested):
  labels_size <- lapply(seq_along(ranking), function(j) {
    unique(unlist(lapply(ranking[seq_len(j)], labels_gauss)))
  })
  if (!all(unlist(labels_size) %in% names(des$cols_by_term))) {
    return(NULL)
  }

  out_size <- function(cols, r_chol, y_mat) {
    if (length(cols) == 0) {
      quad <- rep(0, nprj)
      mu_test <- matrix(ybar, nrow = nrow(x_test), ncol = nprj, byrow = TRUE)
    } else {
      beta <- backsolve(r_chol, y_mat)
      quad <- colSums(y_mat^2) + colSums(des$lambda_x[cols] * beta^2)
      mu_test <- sweep(x_test_c[, cols, drop = FALSE] %*% beta, 2, ybar, "+")
    }
    rownames(mu_test) <- rownames(x_test)
    dis <- sqrt((sum_wvar + sum_wr2 - quad) / sum_w)
    ce <- weighted.mean((-ybar^2 * sum_w - quad) / sum_w, p_ref$wdraws_prj)
    sub_summary <- weighted_summary_means(
      y_wobs_test = data.frame(y = y_test, y_oscale = y_oscale_test,
                               wobs = wobs_test),
      data_aux_test = data_test,
      family = refmodel$family,
      wdraws = p_ref$wdraws_prj,
      mu = mu_test,
      dis = dis,
      cl_ref = p_ref$cl,
      wdraws_ref = p_ref$wdraws_orig
    )
    return(list(sub_summary = sub_summary, ce = ce, dis_sub = dis))
  }

  # Grow the Cholesky factor along the ranking (block Cholesky):
  out_by_size <- vector("list", length(nterms))
  cols <- integer(0)
  r_chol <- NULL
  y_mat <- NULL
  for (size_j in seq(0, max(nterms))) {
    if (size_j > 0) {
      cols_new <- setdiff(
        unlist(des$cols_by_term[labels_size[[size_j]]], use.names = FALSE),
        cols
      )
      if (length(cols_new) > 0) {
        a_nn <- des$a_full[cols_new, cols_new, drop = FALSE]
        rhs_n <- c_mat[cols_new, , drop = FALSE]
        if (length(cols) > 0) {
          r_cn <- forwardsolve(t(r_chol),
                               des$a_full[cols, cols_new, drop = FALSE])
          a_nn <- a_nn - crossprod(r_cn)
          rhs_n <- rhs_n - crossprod(r_cn, y_mat)
        }
        # A numerically singular system goes to the draw-wise path:
        r_n <- tryCatch(chol(a_nn), error = function(e) NULL)
        if (is.null(r_n)) {
          return(NULL)
        }
        y_n <- forwardsolve(t(r_n), rhs_n)
        if (length(cols) == 0) {
          r_chol <- r_n
          y_mat <- y_n
        } else {
          r_chol <- rbind(
            cbind(r_chol, r_cn),
            cbind(matrix(0, length(cols_new), length(cols)), r_n)
          )
          y_mat <- rbind(y_mat, y_n)
        }
        cols <- c(cols, cols_new)
      }
    }
    for (k in which(nterms == size_j)) {
      out_by_size[[k]] <- out_size(cols, r_chol, y_mat)
    }
  }
  return(out_by_size)
}
