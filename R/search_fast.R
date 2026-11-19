# Screened forward search (argument `search_screen` of varsel() and
# cv_varsel(), integer m, `NULL` = off). In each step, a cheap criterion ranks
# the candidates, and only the m best get an exact projection: the score
# statistic at the current submodel (search_forward_screen()), or for models
# without it, a latent screen (search_forward_latent()). Returns `NULL` if
# `search_screen` is `NULL` or no screen applies; then search_forward()
# continues with its own search (with a message in verbose mode).
search_forward_fast <- function(p_ref, refmodel, nterms_max, verbose,
                                search_terms, search_screen = NULL, ...) {
  if (is.null(search_screen)) {
    return(NULL)
  }
  out <- search_forward_screen(
    p_ref = p_ref, refmodel = refmodel, nterms_max = nterms_max,
    verbose = verbose, search_terms = search_terms,
    screen_m = as.integer(search_screen), ...
  )
  if (is.null(out)) {
    out <- search_forward_latent(
      p_ref = p_ref, refmodel = refmodel, nterms_max = nterms_max,
      verbose = verbose, search_terms = search_terms,
      screen_m = as.integer(search_screen), ...
    )
  }
  if (is.null(out)) {
    verb_out("Argument `search_screen` cannot be used for this model or ",
             "data. Running the plain forward search.", verbose = verbose)
  }
  return(out)
}

# Score-screened forward search for binomial (logit link) and Poisson (log
# link) submodels without group-level or smooth terms. For the current submodel
# A, projected draw s, the observation weights w and the new columns J of a
# candidate, the projection objective sum_i w_i (b(eta_i) - mu_ref_i eta_i) has
# the gradient g = -X_J' (w * (mu_ref - mu_A)) at beta_J = 0 and the efficient
# information
#   I = X_J' V X_J - X_J' V X_A (X_A' V X_A)^-1 X_A' V X_J,  V = diag(w b''),
# where X_A includes the intercept. One Newton step reduces the objective by
# g' I^-1 g / 2. The approximate gain of a candidate is this reduction,
# averaged over the draws with `wdraws_prj` and divided by sum(w) (the scale of
# the cross-entropy `ce`; only the order of the gains matters). The ridge
# penalty of the projection is ignored. The gains only select the candidates
# to project; the exact projections choose the term. The submodels are fitted
# to the reference means without the offsets, so `mu_ref` and `mu_A` are the
# means without the offsets, too.
search_forward_screen <- function(p_ref, refmodel, nterms_max, verbose,
                                  search_terms, screen_m, ...) {
  fam <- refmodel$family
  canon <- c(binomial = "logit", poisson = "log")
  if (!fam$family %in% names(canon) || fam$link != canon[[fam$family]] ||
      isTRUE(fam$for_augdat) || isTRUE(fam$for_latent) ||
      !identical(refmodel$div_minimizer, divmin) || anyNA(p_ref$mu) ||
      !isTRUE(screen_m >= 1L)) {
    return(NULL)
  }
  # The model matrix of all search terms (centered columns; the centering does
  # not change the score statistic of a submodel with an intercept). `regul`
  # only sets the unused `lambda_x`:
  des <- gauss_crossprod_design(refmodel, search_terms, regul = 1)
  if (is.null(des)) {
    return(NULL)
  }
  x_c <- des$x_c
  wobs <- des$wobs
  mu_ref <- unclass(p_ref$mu)
  wdraws <- p_ref$wdraws_prj
  step_cands <- fast_step_cands(search_terms, des$cols_by_term)

  # The approximate gains of the candidates with the new columns `cols_cands`
  # at the submodel fit `outdmin` with the columns `cols_chosen`:
  score_gains <- function(outdmin, cols_cands, cols_chosen) {
    mu_a <- refmodel$mu_fun(outdmin)
    if (!identical(dim(mu_a), dim(mu_ref)) || anyNA(mu_a)) {
      return(NULL)
    }
    return(screen_gains(x_c = x_c, mu_ref = mu_ref, mu_a = mu_a, wobs = wobs,
                        wdraws = wdraws, var_fun = fam$variance,
                        cols_cands = cols_cands, cols_chosen = cols_chosen))
  }

  iq <- ceiling(quantile(seq_len(nterms_max + 1L), 1:10 / 10))
  chosen <- character()
  labels_chosen <- character()
  cols_chosen <- integer(0)
  outdmins <- list()
  cur <- NULL
  for (size in seq_len(nterms_max + 1L)) {
    st <- step_cands(size, chosen, labels_chosen)
    if (is.null(st)) {
      return(NULL)
    }
    cands <- st$cands
    if (is.null(cands)) {
      next
    }
    gain <- NULL
    if (!is.null(cur) && length(cands) > 1L) {
      gain <- score_gains(cur$outdmin, st$cols, cols_chosen)
    }
    fit_idx <- if (is.null(gain)) {
      seq_along(cands)
    } else {
      utils::head(order(gain, decreasing = TRUE), screen_m)
    }
    submodls <- lapply(fit_idx, function(k) {
      proj_to_submodl(c(chosen, cands[k]), p_ref = p_ref, refmodel = refmodel,
                      verbose = verbose, ...)
    })
    ibest <- which.min(vapply(submodls, `[[`, numeric(1), "ce"))
    if (length(ibest) == 0) {
      return(NULL)
    }
    k <- fit_idx[ibest]
    cur <- submodls[[ibest]]
    rm(submodls)
    chosen <- c(chosen, cands[k])
    labels_chosen <- c(labels_chosen, st$labels_new[[k]])
    cols_chosen <- c(cols_chosen, st$cols[[k]])
    outdmins <- c(outdmins, list(cur$outdmin))

    # Verbose mode output, as in search_forward():
    ct_chosen <- count_terms_chosen(chosen)
    if (verbose && ct_chosen %in% iq) {
      vtxt <- paste(names(iq)[max(which(ct_chosen == iq))], "of terms selected")
      if (verbose >= 2L) {
        vtxt <- paste0(vtxt, ": `~ ", paste(chosen, collapse = " + "), "`")
      }
      verb_out(vtxt)
    }
  }
  return(nlist(predictor_ranking = setdiff(chosen, "1"), outdmins))
}

# The approximate gains of search_forward_screen() for the candidates with the
# new columns `cols_cands` (a list of indices into the columns of `x_c`) at a
# submodel with the columns `cols_chosen` and an intercept. `mu_ref` and `mu_a`
# are N x S matrices of the reference and submodel means, `wobs` the
# observation weights, `wdraws` the S draw weights and `var_fun` the variance
# function of the family. Returns `NULL` if X_A' V X_A is not positive definite
# for a draw (then all candidates of the step are projected).
screen_gains <- function(x_c, mu_ref, mu_a, wobs, wdraws, var_fun, cols_cands,
                         cols_chosen) {
  g_all <- crossprod(x_c, (mu_ref - mu_a) * wobs)
  v_all <- var_fun(mu_a) * wobs
  x_a <- cbind(1, x_c[, cols_chosen, drop = FALSE])
  is_single <- lengths(cols_cands) == 1L
  cols_single <- unlist(cols_cands[is_single])
  cols_all <- unique(unlist(cols_cands))
  # Column subsets outside the draw loop, so that the loop allocates only
  # N x ncol(x_a) per draw:
  x_cand <- x_c[, cols_all, drop = FALSE]
  x_sq_single <- x_c[, cols_single, drop = FALSE]^2
  pos <- setNames(seq_along(cols_all), cols_all)
  gain <- numeric(length(cols_cands))
  for (s in seq_along(wdraws)) {
    v <- v_all[, s]
    xav <- x_a * v
    r_a <- tryCatch(chol(crossprod(x_a, xav)), error = function(e) NULL)
    if (is.null(r_a)) {
      return(NULL)
    }
    # u' u = X_cand' V X_A (X_A' V X_A)^-1 X_A' V X_cand:
    u <- backsolve(r_a, crossprod(xav, x_cand), transpose = TRUE)
    gain_s <- numeric(length(cols_cands))
    if (length(cols_single) > 0) {
      info_tot <- drop(crossprod(v, x_sq_single))
      info <- info_tot -
        colSums(u[, pos[as.character(cols_single)], drop = FALSE]^2)
      # A column (almost) in the span of the chosen columns gets gain 0:
      gain_s[is_single] <- ifelse(info > 1e-10 * info_tot,
                                  g_all[cols_single, s]^2 / info, 0)
    }
    for (k in which(!is_single)) {
      cj <- cols_cands[[k]]
      if (length(cj) == 0) {
        next
      }
      uj <- u[, pos[as.character(cj)], drop = FALSE]
      x_j <- x_c[, cj, drop = FALSE]
      info <- crossprod(x_j, x_j * v) - crossprod(uj)
      gj <- g_all[cj, s]
      # A singular information (for example, a block in the span of the
      # chosen columns) gives gain 0:
      r_j <- tryCatch(chol(info), error = function(e) NULL)
      gain_s[k] <- if (is.null(r_j) ||
                       min(diag(r_j))^2 <= 1e-10 * max(diag(info))) {
        0
      } else {
        sum(backsolve(r_j, gj, transpose = TRUE)^2)
      }
    }
    gain <- gain + wdraws[s] * gain_s
  }
  return(gain / (2 * sum(wobs)))
}

# The candidates of one step for search terms without group-level or smooth
# terms: `cands` as from select_possible_terms_size() (`NULL` if no candidate
# has the size), the new term labels `labels_new` and the new model-matrix
# columns `cols` (indices into `cols_by_term`) of each candidate. Returns a
# function of (size, chosen, labels_chosen), which returns `NULL` if a new
# label has no columns.
fast_step_cands <- function(search_terms, cols_by_term) {
  function(size, chosen, labels_chosen) {
    cands <- select_possible_terms_size(chosen, search_terms, size = size,
                                        use_term_info = TRUE)
    if (is.null(cands)) {
      return(list(cands = NULL))
    }
    labels_new <- lapply(cands, function(cand) {
      setdiff(labels_gauss(cand), labels_chosen)
    })
    if (!all(unlist(labels_new) %in% names(cols_by_term))) {
      return(NULL)
    }
    cols <- lapply(labels_new, function(lbls) {
      unlist(cols_by_term[lbls], use.names = FALSE)
    })
    return(list(cands = cands, labels_new = labels_new, cols = cols))
  }
}
