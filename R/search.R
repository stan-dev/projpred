# `return_outdmins = FALSE` tells the search that the caller does not use the
# submodel fits along the ranking (element `outdmins`), so
# search_forward_gauss() does not compute them. The draw-wise search computes
# them anyway.
search_forward <- function(p_ref, refmodel, nterms_max, verbose = TRUE,
                           search_terms, est_runtime = TRUE,
                           search_terms_was_null, return_outdmins = TRUE,
                           ...) {
  nterms_max_with_icpt <- nterms_max + 1L
  iq <- ceiling(quantile(seq_len(nterms_max_with_icpt), 1:10 / 10))
  if (is.null(search_terms)) {
    stop("Did not expect `search_terms` to be `NULL`. Please report this.")
  }

  # Gaussian submodels: rank the terms from cross-products (see
  # search_forward_gauss()); `NULL` if this case does not apply:
  if (getOption("projpred.search_gauss_crossprod", TRUE)) {
    out_gauss <- search_forward_gauss(
      p_ref = p_ref, refmodel = refmodel, nterms_max = nterms_max,
      verbose = verbose, search_terms = search_terms,
      return_outdmins = return_outdmins, ...
    )
    if (!is.null(out_gauss)) {
      return(out_gauss)
    }
  }

  chosen <- character()
  outdmins <- c()

  for (size in seq_len(nterms_max_with_icpt)) {
    # Determine candidate predictors for the current size:
    cands <- select_possible_terms_size(chosen, search_terms, size = size,
                                        use_term_info = TRUE)
    if (is.null(cands))
      next
    full_cands <- lapply(cands, function(cand) c(chosen, cand))

    # Perform the projections (for the intercept-only model, we measure the
    # runtime to estimate the total runtime for the search):
    if (size == 1 && est_runtime && getOption("projpred.mssg_time", TRUE)) {
      time_bef <- Sys.time()
    }
    submodls <- lapply(full_cands, proj_to_submodl, p_ref = p_ref,
                       refmodel = refmodel, verbose = verbose, ...)
    if (size == 1 && est_runtime && getOption("projpred.mssg_time", TRUE)) {
      time_aft <- Sys.time()
      dtime <- difftime(time_aft, time_bef, units = "secs")
      # Scale up to the remaining forward search where we have p * (p + 1) / 2
      # projections (with `p = nterms_max`) if `search_terms` is at its default:
      time_est_min <- dtime * nterms_max * (nterms_max + 1) / 2
      # Scale from the intercept-only submodel to GLM submodels (see PR #459 for
      # an empirical derivation of the factor):
      time_est_min <- time_est_min * 1.3
      # Initialize upper bound:
      time_est_max <- time_est_min
      # Adjust for multilevel and/or additive terms if necessary (see PR #459
      # for an empirical derivation of the factors):
      has_mlvl <- any(grepl("\\|", search_terms))
      has_smooth <- length(parse_additive_terms(search_terms)) > 0
      if (has_mlvl && !has_smooth) {
        time_est_max <- time_est_max * 26.9
      } else if (!has_mlvl && has_smooth) {
        time_est_max <- time_est_max * 9.8
      } else if (has_mlvl && has_smooth) {
        time_est_max <- time_est_max * 57.6
      }
      if (time_est_max > 3 * 60) {
        mssg_time_start <- "The remaining forward search is estimated to take "
        if (time_est_max == time_est_min) {
          mssg_time_est <- paste0(
            "ca. ", round(time_est_min / 60, 1), " minutes "
          )
        } else {
          mssg_time_est <- paste0(
            "between ca. ", round(time_est_min / 60, 1), " and ",
            round(time_est_max / 60, 1), " minutes "
          )
        }
        mssg_time_start_curr <- paste0(
          "(current time: ", format(time_aft, usetz = TRUE),
          "; estimated end time: "
        )
        if (time_est_max == time_est_min) {
          mssg_time_est_curr <- paste0(
            format(time_aft + time_est_min, usetz = TRUE), "). Note that ",
            "this is only a rough estimate."
          )
        } else {
          mssg_time_est_curr <- paste0(
            "between ", format(time_aft + time_est_min, usetz = TRUE), " and ",
            format(time_aft + time_est_max, usetz = TRUE), "). Note that ",
            "these are not guaranteed bounds but the bounds of a rough ",
            "interval estimate."
          )
        }
        mssg_time <- paste0(mssg_time_start, mssg_time_est,
                            mssg_time_start_curr, mssg_time_est_curr)
        if (has_mlvl) {
          # The empirically derived factor used above assumes that *all*
          # projections involve one multilevel term (and apart from that, the
          # empirical derivation used only group-level intercepts, not
          # group-level slopes):
          mssg_time <- paste0(
            mssg_time, " Furthermore, since there are multilevel predictor ",
            "terms, this estimate may be even more unreliable."
          )
        }
        if (has_smooth) {
          # The empirically derived factor used above assumes that *all*
          # projections involve one multilevel and one smooth term (and apart
          # from that, the empirical derivation used only group-level intercepts
          # and an s() smooth term, not more complex terms):
          mssg_time <- paste0(
            mssg_time, " Furthermore, since there are additive (\"smooth\") ",
            "predictor terms, this estimate may be even more unreliable."
          )
        }
        if (!search_terms_was_null) {
          # In case of custom `search_terms`, some model sizes may be skipped:
          mssg_time <- paste0(
            mssg_time, " Since argument `search_terms` differs from its ",
            "default, this estimate may be an overestimate."
          )
        }
        message(mssg_time)
      }
    }

    # Select best candidate:
    imin <- which.min(sapply(submodls, "[[", "ce"))
    chosen <- c(chosen, cands[imin])

    # Store `outdmin` (i.e., the object containing the projection results)
    # corresponding to the best candidate:
    outdmins <- c(outdmins, list(submodls[[imin]]$outdmin))

    # Verbose mode output:
    ct_chosen <- count_terms_chosen(chosen)
    if (verbose && ct_chosen %in% iq) {
      vtxt <- paste(names(iq)[max(which(ct_chosen == iq))], "of terms selected")
      if (verbose >= 2L) {
        vtxt <- paste0(vtxt, ": `~ ", paste(chosen, collapse = " + "), "`")
      }
      verb_out(vtxt)
    }

    # Free up some memory:
    rm(submodls)
    if (getOption("projpred.run_gc", FALSE)) {
      gc()
    }
  }

  # For `predictor_ranking`, `reduce_models(chosen)` used to be used instead of
  # `chosen`. However, `reduce_models(chosen)` and `chosen` should be identical
  # at this place because select_possible_terms_size() already avoids redundant
  # models. Thus, use `chosen` here because it matches `outdmins` (this
  # matching is necessary because later in perf_eval()'s `!refit_prj` case,
  # `outdmins` is indexed with integers which are based on `predictor_ranking`):
  return(nlist(predictor_ranking = setdiff(chosen, "1"), outdmins))
}

# For search terms that are each "1" or have one term label (see
# search_term_info()): the label of each term (`NA` for "1") and the string
# that select_possible_terms_size() returns for it. Returns `NULL` in other
# cases, or if the label of a returned string is not the label of its term.
single_label_info <- function(term_info) {
  if (is.null(term_info)) {
    return(NULL)
  }
  is_single <- vapply(term_info, function(info) {
    length(info$labels) <= 1L && !is.na(info$label_out)
  }, logical(1))
  if (!all(is_single)) {
    return(NULL)
  }
  lab <- vapply(term_info, function(info) {
    if (length(info$labels) == 1L) info$labels else NA_character_
  }, character(1), USE.NAMES = FALSE)
  out <- vapply(term_info, `[[`, character(1), "label_out",
                USE.NAMES = FALSE)
  lab_out <- lapply(out, labels_gauss)
  if (!identical(lengths(lab_out), as.integer(!is.na(lab))) ||
      !identical(as.character(unlist(lab_out)), lab[!is.na(lab)])) {
    return(NULL)
  }
  return(list(lab = lab, out = out))
}

# The candidates of one search step for the output `single` of
# single_label_info(): the result `cands` of select_possible_terms_size() and
# the index `idx` of the search term of each candidate, without formula code.
# The new term label of a candidate is `single$lab[idx]` (`NA`: no new label).
# With one label per term, the size of a candidate is the number of chosen
# labels, plus one for a new label, plus one for the intercept. `cands` is
# `NULL` if no candidate has the size. Returns `NULL` if a candidate is a term
# with a chosen label (select_possible_terms_size() builds a formula for it).
# In a search, this does not occur: each step adds at most one label, so
# after the intercept step, the number of chosen labels is at most `size - 2`.
cands_single_label <- function(single, labels_chosen, size) {
  in_chosen <- single$lab %in% labels_chosen
  n_new <- as.integer(!is.na(single$lab) & !in_chosen)
  is_valid <- n_new + length(labels_chosen) + 1L == size
  if (any(is_valid & in_chosen)) {
    return(NULL)
  }
  if (!any(is_valid)) {
    return(list(cands = NULL, idx = NULL))
  }
  idx <- which(is_valid)
  cands <- single$out[idx]
  keep <- !duplicated(cands)
  return(list(cands = cands[keep], idx = idx[keep]))
}

# Forward search for Gaussian submodels with the identity link, from the
# cross-products of the model matrix of all search terms. Let J be the
# model-matrix columns of a candidate. For draw s, the projection solves the
# weighted ridge problem of glm_ridge():
#   (G[J, J] + regul * diag(sd[J]^2)) beta_s = C[J, s].
# Here, G = Xc' W Xc and C = Xc' W M. Xc is the model matrix, centered at the
# weighted column means. W holds the observation weights, M the fitted values
# of the projected reference draws, and sd the weighted column standard
# deviations from .weighted_sd(). glm_ridge() centers and scales the columns
# in the same way and does not penalize the intercept. The part of ce_gauss()
# that differs between the candidates is
#   -beta_s' (G[J, J] + 2 * regul * diag(sd[J]^2)) beta_s / sum(W).
# Thus, the cost per candidate does not depend on the number of observations.
# This function gives only the ranking. The submodel fits along the ranking come
# from proj_to_submodl(), as in search_forward(). The function returns `NULL`
# if one of its conditions does not hold; then search_forward() continues with
# its own search.
search_forward_gauss <- function(p_ref, refmodel, nterms_max, verbose,
                                 search_terms, return_outdmins = TRUE, ...) {
  # The arguments for the divergence minimizer, as in proj_to_submodl():
  dot_args <- list(...)
  divmin_args <- dot_args$search_control %||%
    dot_args[setdiff(names(dot_args), "search_control")]
  regul <- gauss_crossprod_regul(refmodel, divmin_args)
  if (is.null(regul)) {
    return(NULL)
  }
  # Custom search terms are possible, because the candidates come from
  # select_possible_terms_size() (or cands_single_label(), which gives the same
  # result), as in search_forward():
  des <- gauss_crossprod_design(refmodel, search_terms, regul)
  if (is.null(des)) {
    return(NULL)
  }
  cols_by_term <- des$cols_by_term
  # The weights scale the (small) model matrix, not the draws:
  c_mat <- crossprod(des$x_cw, p_ref$mu)
  lambda_x <- des$lambda_x
  a_full <- des$a_full
  wdraws <- p_ref$wdraws_prj

  # The term labels of a candidate string. The same strings come back in every
  # step, so the labels are computed once. `[[` and assign() make a symbol from
  # the string, and a symbol has a limit of 10000 bytes, so long strings are
  # not cached:
  label_cache <- new.env(hash = TRUE, parent = emptyenv())
  labels_of <- function(cand) {
    if (nchar(cand, type = "bytes") >= 10000L) {
      return(labels_gauss(cand))
    }
    lbls <- label_cache[[cand]]
    if (is.null(lbls)) {
      lbls <- labels_gauss(cand)
      assign(cand, lbls, envir = label_cache)
    }
    return(lbls)
  }

  # The criterion is -sum_s wdraws_s * beta_s' (G[J, J] + 2 * Lambda) beta_s
  # (smaller is better, as for the cross-entropy). With
  # G[J, J] + Lambda = R' R and y_s = R^-T C[J, s], it equals
  # -sum_s wdraws_s * (||y_s||^2 + beta_s' Lambda beta_s). In one step, all
  # candidates share the columns of the chosen terms, so the Cholesky factor
  # `r_chosen` and `y_chosen` of these columns are computed once per step. Each
  # candidate extends them by its new columns `cols_new` (block Cholesky).
  crit_cand <- function(cols_new, cols_chosen, r_chosen, y_chosen) {
    quad <- 0
    r_cn <- NULL
    r_n <- NULL
    y_n <- NULL
    beta_n <- NULL
    if (length(cols_new) > 0) {
      a_nn <- a_full[cols_new, cols_new, drop = FALSE]
      rhs_n <- c_mat[cols_new, , drop = FALSE]
      if (length(cols_chosen) > 0) {
        r_cn <- forwardsolve(t(r_chosen),
                             a_full[cols_chosen, cols_new, drop = FALSE])
        a_nn <- a_nn - crossprod(r_cn)
        rhs_n <- rhs_n - crossprod(r_cn, y_chosen)
      }
      r_n <- tryCatch(chol(a_nn), error = function(e) NULL)
      if (is.null(r_n)) {
        return(list(crit = NA_real_))
      }
      y_n <- forwardsolve(t(r_n), rhs_n)
      beta_n <- backsolve(r_n, y_n)
      quad <- colSums(y_n^2) + colSums(lambda_x[cols_new] * beta_n^2)
    }
    if (length(cols_chosen) > 0) {
      y_c <- y_chosen
      if (!is.null(beta_n)) {
        y_c <- y_c - r_cn %*% beta_n
      }
      beta_c <- backsolve(r_chosen, y_c)
      quad <- quad + colSums(y_chosen^2) +
        colSums(lambda_x[cols_chosen] * beta_c^2)
    }
    return(list(crit = -sum(wdraws * quad), r_cn = r_cn, r_n = r_n, y_n = y_n))
  }

  # crit_cand() for many candidates that each add one column (`cols_single`)
  # at once, with matrix operations instead of one call per candidate. For
  # candidate j, beta_c = b0 - u_j beta_n[j, ]' with b0 = R^-1 y_chosen and
  # u_j = R^-1 r_cj, so the penalty term of the chosen columns is
  # sum(lambda b0^2) - 2 beta_n[j, ] p[j, ] + beta_n[j, ]^2 l[j], with
  # p = (u * lambda)' b0 and l = colSums(lambda u^2).
  crit_single <- function(cols_single, cols_chosen, r_chosen, y_chosen) {
    a_jj <- a_full[cbind(cols_single, cols_single)]
    rhs <- c_mat[cols_single, , drop = FALSE]
    if (length(cols_chosen) > 0) {
      r_cj <- forwardsolve(t(r_chosen),
                           a_full[cols_chosen, cols_single, drop = FALSE])
      a_jj <- a_jj - colSums(r_cj^2)
      rhs <- rhs - crossprod(r_cj, y_chosen)
    }
    if (!isTRUE(all(a_jj > 0))) {
      return(rep(NA_real_, length(cols_single)))
    }
    y_n <- rhs / sqrt(a_jj)
    beta_n <- y_n / sqrt(a_jj)
    quad <- y_n^2 + lambda_x[cols_single] * beta_n^2
    if (length(cols_chosen) > 0) {
      lambda_c <- lambda_x[cols_chosen]
      b0 <- backsolve(r_chosen, y_chosen)
      u <- backsolve(r_chosen, r_cj)
      q0 <- colSums(y_chosen^2) + colSums(lambda_c * b0^2)
      quad <- quad +
        matrix(q0, nrow = length(cols_single), ncol = length(q0), byrow = TRUE) -
        2 * beta_n * crossprod(u * lambda_c, b0) +
        colSums(lambda_c * u^2) * beta_n^2
    }
    return(-drop(quad %*% wdraws))
  }

  chosen <- character()
  labels_chosen <- character()
  cols_chosen <- integer(0)
  r_chosen <- NULL
  y_chosen <- NULL
  term_info <- search_term_info(search_terms)
  single <- single_label_info(term_info)
  if (!is.null(single)) {
    # The columns of each search term, looked up by name once for all steps:
    # `integer(0)` for "1" (no label), `NULL` if the model matrix has no
    # columns for the label.
    cols_term <- unname(cols_by_term[single$lab])
    cols_term[is.na(single$lab)] <- list(integer(0))
    term_has_cols <- !vapply(cols_term, is.null, logical(1))
  }
  for (size in seq_len(nterms_max + 1L)) {
    step <- NULL
    if (!is.null(single)) {
      step <- cands_single_label(single, labels_chosen, size)
    }
    if (!is.null(step)) {
      cands <- step$cands
      if (is.null(cands)) {
        next
      }
      if (!all(term_has_cols[step$idx])) {
        return(NULL)
      }
      cols_cands <- cols_term[step$idx]
    } else {
      cands <- select_possible_terms_size(chosen, search_terms, size = size,
                                          use_term_info = TRUE)
      if (is.null(cands)) {
        next
      }
      labels_new <- lapply(cands, function(cand) {
        setdiff(labels_of(cand), labels_chosen)
      })
      if (!all(unlist(labels_new) %in% names(cols_by_term))) {
        return(NULL)
      }
      cols_cands <- lapply(labels_new, function(lbls) {
        unlist(cols_by_term[lbls], use.names = FALSE)
      })
    }
    crits <- rep(NA_real_, length(cands))
    is_single <- lengths(cols_cands) == 1L
    if (any(is_single)) {
      crits[is_single] <- crit_single(unlist(cols_cands[is_single]),
                                      cols_chosen, r_chosen, y_chosen)
    }
    for (k in which(!is_single)) {
      crits[k] <- crit_cand(cols_cands[[k]], cols_chosen, r_chosen,
                            y_chosen)$crit
    }
    # A numerically singular system (for example, collinear columns with a
    # small `regul`) goes to the search of search_forward():
    if (anyNA(crits)) {
      return(NULL)
    }
    imin <- which.min(crits)
    chosen <- c(chosen, cands[imin])
    cols_new <- cols_cands[[imin]]
    if (length(cols_new) > 0) {
      # The factors of the chosen candidate, for the next step:
      best <- crit_cand(cols_new, cols_chosen, r_chosen, y_chosen)
      if (is.na(best$crit)) {
        return(NULL)
      }
      if (length(cols_chosen) == 0) {
        r_chosen <- best$r_n
        y_chosen <- best$y_n
      } else {
        r_chosen <- rbind(
          cbind(r_chosen, best$r_cn),
          cbind(matrix(0, length(cols_new), length(cols_chosen)), best$r_n)
        )
        y_chosen <- rbind(y_chosen, best$y_n)
      }
      cols_chosen <- c(cols_chosen, cols_new)
      label_new <- if (!is.null(step)) {
        single$lab[step$idx[imin]]
      } else {
        labels_new[[imin]]
      }
      labels_chosen <- c(labels_chosen, label_new)
    }
  }

  # The submodel fits along the ranking (as stored by search_forward()), if the
  # caller uses them:
  outdmins <- NULL
  if (return_outdmins) {
    outdmins <- lapply(seq_along(chosen), function(j) {
      proj_to_submodl(predictor_terms = chosen[seq_len(j)], p_ref = p_ref,
                      refmodel = refmodel, verbose = verbose, ...)$outdmin
    })
  }
  return(nlist(predictor_ranking = setdiff(chosen, "1"), outdmins))
}

#' Force search terms
#'
#' A helper function to construct the input for argument `search_terms` of
#' [varsel()] or [cv_varsel()] if certain predictor terms should be forced to be
#' selected first whereas other predictor terms are optional (i.e., they are
#' subject to the variable selection, but only after the inclusion of the
#' "forced" terms).
#'
#' @param forced_terms A character vector of predictor terms that should be
#'   selected first.
#' @param optional_terms A character vector of predictor terms that should be
#'   subject to the variable selection after the inclusion of the "forced"
#'   terms.
#'
#' @return A character vector that may be used as input for argument
#'   `search_terms` of [varsel()] or [cv_varsel()].
#'
#' @seealso [varsel()], [cv_varsel()]
#'
#' @examplesIf requireNamespace("rstanarm", quietly = TRUE)
#' # Data:
#' dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)
#'
#' # The `stanreg` fit which will be used as the reference model (with small
#' # values for `chains` and `iter`, but only for technical reasons in this
#' # example; this is not recommended in general):
#' fit <- rstanarm::stan_glm(
#'   y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
#'   QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
#' )
#'
#' # We will force X1 and X2 to be selected first:
#' search_terms_forced <- force_search_terms(
#'   forced_terms = paste0("X", 1:2),
#'   optional_terms = paste0("X", 3:5)
#' )
#'
#' # Run varsel() (here without cross-validation and with small values for
#' # `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the sake of
#' # speed in this example; this is not recommended in general):
#' vs <- varsel(fit, nclusters = 5, nclusters_pred = 10,
#'              search_terms = search_terms_forced, seed = 5555)
#' # Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
#' # and `?ranking` for possible post-processing functions.
#'
#' @export
force_search_terms <- function(forced_terms, optional_terms) {
  stopifnot(length(forced_terms) > 0)
  stopifnot(length(optional_terms) > 0)
  forced_terms <- paste(forced_terms, collapse = " + ")
  return(c(forced_terms, paste0(forced_terms, " + ", optional_terms)))
}

search_L1_surrogate <- function(p_ref, d_train, family, intercept, nterms_max,
                                penalty, search_control) {

  ## predictive mean and variance of the reference model (with parameters
  ## integrated out)
  mu <- p_ref$mu
  v <- p_ref$var

  if (NCOL(mu) > 1 || NCOL(v) > 1) {
    stop("Internal error: search_L1 received multiple draws. ",
         "Please report to the developers.")
  }

  ## L1-penalized projection (projection path).
  ## (Notice: here we use pmax = nterms_max+1 so that the computation gets
  ## carried until all the way down to the least regularization also for model
  ## size nterms_max)
  out_capt <- utils::capture.output(
    search <- glm_elnet(
      d_train$x, mu, family,
      lambda_min_ratio = search_control$lambda_min_ratio %||% 1e-5,
      nlambda = search_control$nlambda %||% 150, pmax = nterms_max + 1,
      pmax_strict = FALSE, weights = d_train$weights, intercept = intercept,
      obsvar = v, penalty = penalty, thresh = search_control$thresh %||% 1e-6
    )
  )
  out_capt <- unique(grep("[Ww]arning|bug", out_capt, value = TRUE))
  if (length(out_capt) > 0) {
    warning(paste(
      c(paste0("The following warnings have been thrown by projpred's ",
               "internal L1-search function:"),
        "---", out_capt, "---",
        paste0("It is recommended to inspect this in detail and (if ",
               "necessary) to adjust tuning parameters via argument ",
               "`search_control` (of varsel() or cv_varsel()).")),
      collapse = "\n"
    ))
  }

  ## sort the variables according to the order in which they enter the model in
  ## the L1-path
  entering_indices <- apply(search$beta != 0, 1, function(num) {
    which(num)[1] # na for those that did not enter
  })
  ## variables that entered at some point
  entered_variables <- c(seq_len(NCOL(d_train$x)))[!is.na(entering_indices)]
  ## variables that did not enter at any point
  notentered_variables <- c(seq_len(NCOL(d_train$x)))[is.na(entering_indices)]
  order_of_entered <- sort(entering_indices, index.return = TRUE)$ix
  order <- c(entered_variables[order_of_entered], notentered_variables)

  ## fetch the coefficients corresponding to those points at the search_path
  ## where new variable enters
  nvar <- length(order)
  n <- nrow(p_ref$mu)
  out <- list(
    alpha = rep(NA, nterms_max + 1),
    beta = matrix(0, nrow = nterms_max, ncol = nterms_max + 1),
    lambda = rep(NA, nterms_max + 1),
    w = matrix(NA, nrow = n, ncol = nterms_max + 1)
  )
  for (k in 0:nterms_max) {
    if (k == 0) {
      out$alpha[1] <- search$beta0[1]
      out$lambda[1] <- search$lambda[1]
      out$w[, 1] <- search$w[, 1]
    } else {
      ## find those points in the L1-path where only the k most relevant
      ## features can have nonzero coefficient, and then fetch their
      ## coefficients with least regularization
      ivar <- utils::tail(order, nvar - k)
      steps_k_var <- which(colSums(search$beta[ivar, , drop = FALSE] != 0) == 0)
      if (length(steps_k_var) > 0) {
        j <- utils::tail(steps_k_var, 1)
      } else {
        ## no steps where all the variables in set ivar would have zero
        ## coefficient (could be due to one or more of these variables having
        ## penalty = 0 so they are always in the model) so set the coefficients
        ## to be equal to the starting value
        j <- 1
      }
      out$alpha[k + 1] <- search$beta0[j]
      out$beta[1:k, k + 1] <- search$beta[order[1:k], j]
      out$lambda[k + 1] <- search$lambda[j]
      out$w[, k + 1] <- search$w[, j]
    }
  }

  out$predictor_ranking <- order[seq_len(nterms_max)]
  if (any(is.na(out$predictor_ranking)) &&
      length(entered_variables) < nterms_max) {
    if (length(setdiff(notentered_variables,
                       which(penalty == Inf))) > 0) {
      warning("Less than nterms_max variables entered L1-path. ",
              "Try reducing lambda_min_ratio. ")
    }
  }

  return(out)
}

search_L1 <- function(p_ref, refmodel, nterms_max, penalty, search_control) {
  if (nterms_max == 0) {
    stop("L1 search cannot be used for an empty (i.e. intercept-only) ",
         "full-model formula or `nterms_max = 0`.")
  }
  # Preparations:
  fr <- model.frame(refmodel$formula, data = refmodel$fetch_data(),
                    drop.unused.levels = TRUE)
  da_classes <- attr(attr(fr, "terms"), "dataClasses")
  nms_chr_fac <- names(da_classes)[da_classes %in% c("character", "factor")]
  resp_nm <- all.vars(attr(fr, "terms"))[attr(attr(fr, "terms"), "response")]
  nms_chr_fac <- setdiff(nms_chr_fac, resp_nm)
  if (length(nms_chr_fac) > 0) {
    xlvls <- lapply(setNames(nm = nms_chr_fac), function(nm_chr_fac) {
      levels(as.factor(fr[[nm_chr_fac]]))
    })
  } else {
    xlvls <- NULL
  }
  # TODO: In the following model.matrix() call, allow user-specified contrasts
  # to be passed to argument `contrasts.arg`. The `contrasts.arg` default
  # (`NULL`) uses `options("contrasts")` internally, but it might be more
  # convenient to let users specify contrasts directly. At that occasion,
  # contrasts should also be tested thoroughly (not done until now).
  x <- model.matrix(refmodel$formula, data = fr)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  ## it's important to keep the original order because that's the order
  ## in which lasso will estimate the parameters
  tt <- terms(refmodel$formula)
  terms_ <- attr(tt, "term.labels")
  search_path <- search_L1_surrogate(
    p_ref, nlist(x, weights = refmodel$wobs), refmodel$family,
    intercept = TRUE, ncol(x), penalty, search_control
  )

  predictor_ranking_orig <- collapse_ranked_predictors(
    path = colnames(x)[search_path$predictor_ranking],
    formula = refmodel$formula, data = fr
  )
  predictor_ranking <- utils::head(predictor_ranking_orig, nterms_max)
  # Place lower-order interaction terms before higher-order interaction terms,
  # but otherwise preserve the ranking:
  ia_orders <- sapply(gregexpr(":", predictor_ranking), function(greg_colon) {
    sum(greg_colon != -1)
  })
  ia_order_max <- max(ia_orders)
  for (ia_order in rev(seq_len(ia_order_max))) {
    ias <- predictor_ranking[ia_orders == ia_order]
    stopifnot(!any(duplicated(ias))) # safety measure for which.max()
    for (ia in ias) {
      ia_idx <- which.max(predictor_ranking == ia)
      if (ia_idx > nterms_max) break
      main_terms_ia <- strsplit(ia, ":")[[1]]
      ias_lower_split <- utils::combn(main_terms_ia, m = ia_order,
                                      simplify = FALSE)
      ias_lower <- lapply(ias_lower_split, all_ia_perms, is_split = TRUE)
      ias_lower <- unlist(ias_lower)
      ias_lower <- intersect(ias_lower, predictor_ranking_orig)
      prev_terms <- utils::head(predictor_ranking, ia_idx - 1L)
      has_lower_after <- !all(ias_lower %in% prev_terms)
      if (has_lower_after) {
        if (getOption("projpred.warn_L1_interactions", TRUE)) {
          warning("Interaction term `", ia, "` was selected before all ",
                  "corresponding lower-order interaction terms have been ",
                  "selected. This is a known deficiency of L1 search. Use ",
                  "forward search to avoid this. Now ranking the lower-order ",
                  "interaction terms before this interaction term.")
        }
        ias_lower <- setdiff(ias_lower, prev_terms)
        ias_lower <- ias_lower[order(match(ias_lower, predictor_ranking_orig))]
        new_head <- c(prev_terms, ias_lower, ia)
        predictor_ranking <- c(new_head, setdiff(predictor_ranking, new_head))
        predictor_ranking <- utils::head(predictor_ranking, nterms_max)
      }
    }
  }

  outdmins <- lapply(0:length(predictor_ranking), function(nterms) {
    if (nterms == 0) {
      formula <- make_formula(c("1"))
      beta <- NULL
      x <- x[, numeric(), drop = FALSE]
    } else {
      formula <- make_formula(predictor_ranking[seq_len(nterms)])
      variables <- unlist(lapply(
        predictor_ranking[seq_len(nterms)],
        function(term) {
          # TODO: In the following model.matrix() call, allow user-specified
          # contrasts to be passed to argument `contrasts.arg`. The
          # `contrasts.arg` default (`NULL`) uses `options("contrasts")`
          # internally, but it might be more convenient to let users specify
          # contrasts directly. At that occasion, contrasts should also be
          # tested thoroughly (not done until now).
          mm <- model.matrix(as.formula(paste("~ 1 +", term)), data = fr)
          return(setdiff(colnames(mm), "(Intercept)"))
        }
      ))
      indices <- match(variables, colnames(x)[search_path$predictor_ranking])
      indices <- indices[!is.na(indices)]
      beta <- search_path$beta[indices, max(indices) + 1, drop = FALSE]
      # Also reduce `x` (important for coef.subfit(), for example); note that
      # `x <- x[, variables, drop = FALSE]` should also be possible, but the
      # re-use of `colnames(x)` should provide another sanity check:
      x <- x[
        , colnames(x)[search_path$predictor_ranking[indices]], drop = FALSE
      ]
      # For consistency with fit_glm_ridge_callback():
      rownames(beta) <- colnames(x)
    }
    # Avoid model.frame.default()'s warning "variable '<...>' is not a factor"
    # when calling predict.subfit() later:
    xlvls <- xlvls[intersect(names(xlvls), all.vars(formula))]
    if (length(xlvls) == 0) {
      xlvls <- NULL
    }

    sub <- nlist(alpha = search_path$alpha[nterms + 1], beta,
                 w = search_path$w[, nterms + 1], formula, x, xlvls,
                 terms = attr(fr, "terms"))
    class(sub) <- "subfit"
    return(list(sub))
  })

  return(nlist(predictor_ranking, outdmins))
}
