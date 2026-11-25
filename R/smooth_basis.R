# Projected draws of smooth terms in the basis of the reference model fit
# ------------------------------------------------------------------------
#
# The additive submodels are mgcv::gam() fits in the basis of mgcv. brms and
# rstanarm use other bases for the same smooth terms: brms uses
# mgcv::smoothCon(absorb.cons = TRUE, modCon = 3, diagonal.penalty = TRUE) and
# mgcv::smooth2random(type = 2) (parameters `bs_*`, `s_*`, `sds_*`), rstanarm
# uses mgcv::jagam(diagonalize = TRUE) (parameters `s(x).1`, ...,
# `smooth_sd[...]`). Both bases span the same functions as the basis of mgcv
# (up to a constant, which goes into the intercept). A change of basis does not
# change a penalized fit, so as.matrix.projection() maps the coefficients of
# the submodel fits into the basis of the reference model fit:
#
# 1. B: the columns of the smooth term in the reference model fit (from
#    brms::standata() or the model matrix `x` of a `stanreg` fit), with column
#    blocks that have one SD parameter each.
# 2. The linear predictor of the smooth term in each draw (from
#    predict.gam(type = "lpmatrix")) is regressed on [1, B]. This gives the
#    coefficients in the target basis and a shift of the intercept.
# 3. The penalties hold in the basis of the setup of mgcv::gam() (`fit =
#    FALSE`), which can differ from the basis of predict.gam() (t2() terms).
#    With [1, B] T = X_setup, the penalty k in the target basis is
#    T^-T S_k T^-1. It must be c_bk * I on each block b and 0 elsewhere. Then
#    the SD of block b is sqrt(sig2 / kappa_b) with kappa_b = sum_k(sp_k *
#    c_bk), as in the mixed-model form of the smooth (prior N(0, SD^2 I) of the
#    coefficients of the block). A block without a penalty in the submodel (the
#    null space of a jagam() smooth) gets no SD column. If kappa_b is
#    negligible against the data information of the block (the fit does not
#    penalize the block, for example if the submodel reproduces the reference
#    model exactly), the SD is not identified and is `NA`.

# brms's default renaming of special characters in parameter names. Copied over
# from brms::rename() (brms 2.23.0) with the default `pattern` and
# `replacement`.
brms_rename <- function(x) {
  pattern <- c(" ", "(", ")", "[", "]", ",", "\"", "'", "?", "+", "-", "*",
               "/", "^", "=", "$")
  replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ", "USD")
  for (i in seq_along(pattern)) {
    x <- gsub(pattern[i], replacement[i], x, fixed = TRUE)
  }
  return(x)
}

# The bases of the smooth terms of the reference model fit for naming scheme
# `nm_scheme`, or `NULL` if the fit does not give them. The output is a list
# with one element per smooth term, named by the key that smooth_key() gives
# for the mgcv smooth of the same term. Each element is a list with:
# * `B`: the basis matrix (one row per observation of the reference model
#   data);
# * `blk`: the block of each column of `B` (0 for columns without an SD
#   parameter);
# * `coef_nms`: the parameter name of each column of `B`;
# * `sd_nms`: the name of the SD parameter of each block 1, 2, ....
smooth_targets <- function(fit, nm_scheme) {
  if (nm_scheme == "brms" && inherits(fit, "brmsfit")) {
    return(smooth_targets_brms(fit))
  }
  if (nm_scheme == "rstanarm" && inherits(fit, "stanreg") &&
      !is.null(fit$jam) && !is.null(fit$x)) {
    return(smooth_targets_rstanarm(fit))
  }
  return(NULL)
}

# Only the smooth terms of the linear predictor `mu` of a univariate model
# (Stan data `Xs`, `Zs_<j>_<k>`). A smooth term without `Xs` columns (a basis
# whose null space is only the constant) is not mapped; projpred does not
# accept the arguments of s() that give such a basis.
smooth_targets_brms <- function(fit) {
  sdata <- tryCatch(brms::standata(fit), error = function(e) NULL)
  Xs <- sdata[["Xs"]]
  smcols <- attr(Xs, "smcols")
  if (is.null(Xs) || is.null(smcols)) {
    return(NULL)
  }
  out <- list()
  for (j in seq_along(smcols)) {
    nb <- sdata[[paste0("nb_", j)]]
    if (!length(smcols[[j]]) || is.null(nb)) {
      next
    }
    Xs_nms <- colnames(Xs)[smcols[[j]]]
    lab <- sub("_1$", "", Xs_nms[1])
    Zs <- lapply(seq_len(nb), function(k) sdata[[paste0("Zs_", j, "_", k)]])
    if (any(vapply(Zs, is.null, logical(1)))) {
      next
    }
    nZs <- vapply(Zs, NCOL, integer(1))
    s_nms <- unlist(lapply(seq_len(nb), function(k) {
      paste0("s_", lab, "_", k, "[", seq_len(nZs[k]), "]")
    }))
    out[[lab]] <- list(
      B = unname(cbind(Xs[, smcols[[j]], drop = FALSE],
                       do.call(cbind, Zs))),
      blk = c(rep(0L, length(Xs_nms)), rep(seq_len(nb), nZs)),
      coef_nms = c(paste0("bs_", Xs_nms), s_nms),
      sd_nms = paste0("sds_", lab, "_", seq_len(nb))
    )
  }
  return(out)
}

smooth_targets_rstanarm <- function(fit) {
  X <- as.matrix(fit$x)
  sp_nms <- names(fit$jam$sp)
  out <- list()
  sp_start <- 1L
  for (sm in fit$jam$smooth) {
    nsp <- length(sm$rank)
    sp_idxs <- sp_start + seq_len(nsp) - 1L
    sp_start <- sp_start + nsp
    idxs <- sm$first.para:sm$last.para
    coef_nms <- fit$jam$term.names[idxs]
    if (sum(sm$rank) != length(idxs) || length(sp_nms) < max(sp_idxs) ||
        !all(coef_nms %in% colnames(X))) {
      next
    }
    out[[sm$label]] <- list(
      B = unname(X[, coef_nms, drop = FALSE]),
      blk = rep(seq_len(nsp), sm$rank),
      coef_nms = coef_nms,
      sd_nms = paste0("smooth_sd[", sp_nms[sp_idxs], "]")
    )
  }
  return(out)
}

# The key of an mgcv smooth for the lookup in the output of smooth_targets().
smooth_key <- function(sm, nm_scheme) {
  if (nm_scheme == "brms") {
    return(brms_rename(sm$label))
  }
  return(sm$label)
}

# Replace the columns of the smooth terms in `res` (the output of
# get_subparams() for the `gam` fits `fits`, one row per fit) by the
# parameters in the basis of the reference model fit. Each smooth term for
# which this is not possible keeps its columns in the basis of mgcv::gam(),
# with a warning. `refmodel` is the reference model object. `warn_no_basis`:
# give a warning if the reference model fit does not give the bases.
smooth_to_ref_basis <- function(res, fits, refmodel, nm_scheme, warn_no_basis) {
  fit_1 <- fits[[1]]
  if (!length(fit_1$smooth)) {
    return(res)
  }
  warn_mgcv_basis <- function(reason) {
    warning(reason, " The columns of the smooth terms are the coefficients ",
            "in the basis of mgcv::gam().", call. = FALSE)
  }
  targets <- smooth_targets(refmodel$fit, nm_scheme)
  if (is.null(targets)) {
    if (warn_no_basis) {
      warn_mgcv_basis(paste0(
        "The reference model fit does not give the bases of the smooth terms ",
        "for `nm_scheme = \"", nm_scheme, "\"`."
      ))
    }
    return(res)
  }
  data <- refmodel$fetch_data()
  if (!all(vapply(fits, function(fit) {
    NROW(fit$model) == nrow(data)
  }, logical(1))) ||
      !all(vapply(targets, function(tgt) {
        nrow(tgt$B) == nrow(data)
      }, logical(1)))) {
    warn_mgcv_basis(paste(
      "The submodel fits or the reference model fit do not have the",
      "observations of the reference model data."
    ))
    return(res)
  }
  if (!all(vapply(fits, gam_same_smooths, logical(1), fit_ref = fit_1))) {
    warn_mgcv_basis("The submodel fits of the projected draws differ.")
    return(res)
  }
  setup <- gam_setup_for_basis(fit_1, data)
  lp_mat <- tryCatch(
    stats::predict(fit_1, newdata = data, type = "lpmatrix"),
    error = function(e) NULL
  )
  if (is.null(setup) || is.null(lp_mat)) {
    warn_mgcv_basis("The setup of the submodel fit could not be rebuilt.")
    return(res)
  }
  coefs <- do.call(cbind, lapply(fits, stats::coef))
  icpt_nm <- mknms_icpt("(Intercept)", nm_scheme = nm_scheme)
  res_coef_nms <- colnames(res)
  if (nm_scheme == "brms") {
    icpt_nm <- paste0("b_", icpt_nm)
    res_coef_nms <- sub("^b_", "", res_coef_nms)
  }
  # The new columns of each mapped smooth term, by the name of its first column
  # in `res`:
  new_cols <- list()
  drop_nms <- character()
  for (j in seq_along(fit_1$smooth)) {
    sm <- fit_1$smooth[[j]]
    tgt <- targets[[smooth_key(sm, nm_scheme)]]
    idxs <- sm$first.para:sm$last.para
    res_nms <- colnames(res)[match(names(stats::coef(fit_1))[idxs],
                                   res_coef_nms)]
    mapped <- NULL
    if (!is.null(tgt) && !anyNA(res_nms) && icpt_nm %in% colnames(res)) {
      mapped <- smooth_map(
        lp_j = lp_mat[, idxs, drop = FALSE] %*% coefs[idxs, , drop = FALSE],
        tgt = tgt, setup = setup, j = j, fits = fits
      )
    }
    if (is.null(mapped)) {
      warning("The smooth term `", sm$label, "` of the submodel does not ",
              "match a smooth term of the reference model fit. Its columns ",
              "are the coefficients in the basis of mgcv::gam().",
              call. = FALSE)
      next
    }
    if (!is.null(mapped$sd_problem)) {
      warning("The smooth term `", sm$label, "` of the submodel ",
              mapped$sd_problem, call. = FALSE)
    }
    res[, icpt_nm] <- res[, icpt_nm] + mapped$icpt
    new_cols[[res_nms[1]]] <- mapped$cols
    drop_nms <- c(drop_nms, res_nms)
  }
  if (!length(new_cols)) {
    return(res)
  }
  out <- lapply(colnames(res), function(nm) {
    if (nm %in% names(new_cols)) {
      return(new_cols[[nm]])
    }
    if (nm %in% drop_nms) {
      return(NULL)
    }
    return(res[, nm, drop = FALSE])
  })
  out <- do.call(cbind, out)
  rownames(out) <- rownames(res)
  return(out)
}

# Whether the `gam` fit `fit` has the same smooth terms and coefficients as the
# `gam` fit `fit_ref`.
gam_same_smooths <- function(fit, fit_ref) {
  sm_info <- function(f) {
    lapply(f$smooth, function(sm) {
      list(sm$label, sm$first.para, sm$last.para, sm$first.sp, sm$last.sp)
    })
  }
  identical(sm_info(fit), sm_info(fit_ref)) &&
    identical(names(stats::coef(fit)), names(stats::coef(fit_ref)))
}

# The setup of mgcv::gam() (`fit = FALSE`) for the formula of the `gam` fit
# `fit` and the reference model data `data`, or `NULL`. The response does not
# matter for the basis and the penalties, so it is a constant.
gam_setup_for_basis <- function(fit, data) {
  fml <- fit$formula
  if (length(fml) != 3) {
    return(NULL)
  }
  resp_nm <- ".projpred_smooth_basis_y"
  fml[[2]] <- as.name(resp_nm)
  data[[resp_nm]] <- 0
  setup <- tryCatch(mgcv::gam(fml, data = data, fit = FALSE),
                    error = function(e) NULL)
  if (is.null(setup) || nrow(setup$X) != nrow(data) ||
      length(setup$smooth) != length(fit$smooth)) {
    return(NULL)
  }
  return(setup)
}

# Map smooth term `j` of the `gam` fits `fits` into the target basis `tgt` (an
# element of the output of smooth_targets()). `lp_j` is the linear predictor
# of the term (one column per fit), `setup` the output of
# gam_setup_for_basis(). Returns `NULL` if the bases do not match, else a list
# with the intercept shifts `icpt` (one per fit), the new columns `cols` (one
# row per fit) and `sd_problem` (`NULL` or the end of a warning message about
# the SD columns).
smooth_map <- function(lp_j, tgt, setup, j, fits, tol = 1e-8) {
  B1 <- cbind(1, tgt$B)
  qr_B1 <- qr(B1)
  if (qr_B1$rank < ncol(B1)) {
    return(NULL)
  }
  if (max(abs(qr.resid(qr_B1, lp_j))) > tol * max(abs(lp_j), 1)) {
    return(NULL)
  }
  theta <- qr.coef(qr_B1, lp_j)
  cols <- t(theta[-1, , drop = FALSE])
  colnames(cols) <- tgt$coef_nms
  sd_out <- smooth_sds(setup, j = j, fits = fits, qr_B1 = qr_B1, tgt = tgt,
                       tol = tol)
  if (!is.null(sd_out$sds)) {
    # Within the term, the SD columns come after the columns without an SD
    # (brms: `bs_*`) and before the other coefficients:
    n0 <- sum(tgt$blk == 0L)
    cols <- cbind(cols[, seq_len(n0), drop = FALSE], sd_out$sds,
                  cols[, n0 + seq_len(ncol(cols) - n0), drop = FALSE])
  }
  return(list(icpt = theta[1, ], cols = cols, sd_problem = sd_out$problem))
}

# The SD columns of smooth term `j`: a list with `sds` (one row per fit in
# `fits`, or `NULL`) and `problem` (`NULL` or the end of a warning message).
smooth_sds <- function(setup, j, fits, qr_B1, tgt, tol) {
  no_sds <- function(problem) list(sds = NULL, problem = problem)
  sm <- setup$smooth[[j]]
  idxs <- sm$first.para:sm$last.para
  X_j <- setup$X[, idxs, drop = FALSE]
  pen_problem <- paste(
    "has penalties that do not match the SD parameters of the reference",
    "model fit. Its SD columns are missing."
  )
  if (max(abs(qr.resid(qr_B1, X_j))) > tol * max(abs(X_j), 1)) {
    return(no_sds(pen_problem))
  }
  if (length(idxs) != ncol(tgt$B)) {
    return(no_sds(paste(
      "has fewer basis functions than in the reference model fit (for",
      "example because of an identifiability constraint of mgcv::gam() for a",
      "term that overlaps another one). Its SD columns are missing."
    )))
  }
  T_j <- qr.coef(qr_B1, X_j)[-1, , drop = FALSE]
  T_inv <- tryCatch(solve(T_j), error = function(e) NULL)
  if (is.null(T_inv) || is.null(sm$first.sp) || is.null(sm$last.sp) ||
      sm$last.sp < sm$first.sp) {
    return(no_sds(pen_problem))
  }
  sp_idxs <- sm$first.sp:sm$last.sp
  S_list <- setup$S[sp_idxs]
  if (!length(S_list) || anyNA(setup$off[sp_idxs])) {
    return(no_sds(pen_problem))
  }
  blk <- tgt$blk
  nblk <- length(tgt$sd_nms)
  # c_bk: the factor of penalty k on block b; has_pen: whether block b has a
  # penalty (relative to the largest entry of each penalty).
  c_bk <- matrix(0, nblk, length(S_list))
  has_pen <- rep(FALSE, nblk)
  same_blk <- outer(blk, blk, "==") & outer(blk > 0L, blk > 0L, "&")
  for (k in seq_along(S_list)) {
    S_k <- matrix(0, length(idxs), length(idxs))
    off_k <- setup$off[sp_idxs[k]] - sm$first.para
    S_idxs <- off_k + seq_len(nrow(S_list[[k]]))
    if (off_k < 0 || max(S_idxs) > length(idxs)) {
      return(no_sds(pen_problem))
    }
    S_k[S_idxs, S_idxs] <- S_list[[k]]
    P_k <- crossprod(T_inv, S_k %*% T_inv)
    P_max <- max(abs(P_k))
    if (P_max == 0) {
      next
    }
    if (max(abs(P_k[!same_blk])) > tol * P_max) {
      return(no_sds(pen_problem))
    }
    for (b in seq_len(nblk)) {
      bi <- which(blk == b)
      c_b <- mean(diag(P_k)[bi])
      if (max(abs(P_k[bi, bi] - c_b * diag(length(bi)))) > tol * P_max) {
        return(no_sds(pen_problem))
      }
      c_bk[b, k] <- c_b
      has_pen[b] <- has_pen[b] || abs(c_b) > tol * P_max
    }
  }
  if (!any(has_pen)) {
    return(no_sds(pen_problem))
  }
  # The penalty of a block is negligible if kappa_b is below `tol_ident` times
  # the mean diagonal element of B_b' W B_b (W: the working weights of the
  # fit). Then the SD is not identified. A scale of 0 gives no SD either.
  tol_ident <- 1e-6
  B <- tgt$B
  info_fun <- function(w) {
    vapply(seq_len(nblk), function(b) {
      mean(colSums(B[, blk == b, drop = FALSE]^2 * w))
    }, numeric(1))
  }
  info_const <- info_fun(1)
  sp_list <- lapply(fits, function(fit) {
    if (!is.null(fit$full.sp)) fit$full.sp else fit$sp
  })
  # Only a guard: projpred does not accept arguments of s() and t2() (such as
  # `id` or `sp`) that would give other smoothing parameters.
  if (any(lengths(sp_list) != length(setup$S))) {
    return(no_sds(pen_problem))
  }
  sds <- do.call(rbind, Map(function(fit, sp) {
    kappa_b <- drop(c_bk %*% sp[sp_idxs])
    w <- fit$weights
    info_b <- if (length(w) == nrow(B) && any(w != 1)) {
      info_fun(w)
    } else {
      info_const
    }
    sd_b <- sqrt(fit$sig2 / kappa_b)
    unident <- has_pen & (!(kappa_b > tol_ident * info_b) | !(fit$sig2 > 0))
    sd_b[unident] <- NA_real_
    sd_b
  }, fits, sp_list))
  sds <- sds[, has_pen, drop = FALSE]
  colnames(sds) <- tgt$sd_nms[has_pen]
  problem <- NULL
  if (anyNA(sds)) {
    n_na <- sum(rowSums(is.na(sds)) > 0)
    problem <- paste0(
      "has SD values that are not identified in ", n_na, " of ",
      length(fits), " projected draws, because the submodel fit does not ",
      "penalize the term (for example if the submodel reproduces the ",
      "reference model exactly). These SD values are `NA`."
    )
  }
  return(list(sds = sds, problem = problem))
}
