search_forward <- function(p_ref, refmodel, nterms_max, verbose = TRUE, opt,
                           search_terms, ...) {
  nterms_max_with_icpt <- nterms_max + 1L
  iq <- ceiling(quantile(seq_len(nterms_max_with_icpt), 1:10 / 10))
  if (is.null(search_terms)) {
    stop("Did not expect `search_terms` to be `NULL`. Please report this.")
  }

  chosen <- character()
  outdmins <- c()

  for (size in seq_len(nterms_max_with_icpt)) {
    cands <- select_possible_terms_size(chosen, search_terms, size = size)
    if (is.null(cands))
      next
    full_cands <- lapply(cands, function(cand) c(chosen, cand))
    submodls <- lapply(full_cands, get_submodl_prj, p_ref = p_ref,
                       refmodel = refmodel, regul = opt$regul, ...)

    ## select best candidate
    imin <- which.min(sapply(submodls, "[[", "ce"))
    chosen <- c(chosen, cands[imin])

    ## append `outdmin`
    outdmins <- c(outdmins, list(submodls[[imin]]$outdmin))

    ct_chosen <- count_terms_chosen(chosen)
    if (verbose && ct_chosen %in% iq) {
      vtxt <- paste(names(iq)[max(which(ct_chosen == iq))], "of terms selected")
      if (getOption("projpred.extra_verbose", FALSE)) {
        vtxt <- paste0(vtxt, ": ", paste(chosen, collapse = " + "))
      }
      verb_out(vtxt)
    }
  }

  # For `solution_terms`, `reduce_models(chosen)` used to be used instead of
  # `chosen`. However, `reduce_models(chosen)` and `chosen` should be identical
  # at this place because select_possible_terms_size() already avoids redundant
  # models. Thus, use `chosen` here because it matches `outdmins` (this
  # matching is necessary because later in get_submodls()'s `!refit_prj` case,
  # `outdmins` is indexed with integers which are based on `solution_terms`):
  return(nlist(solution_terms = setdiff(chosen, "1"), outdmins))
}

search_L1_surrogate <- function(p_ref, d_train, family, intercept, nterms_max,
                                penalty, opt) {

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
  search <- glm_elnet(d_train$x, mu, family,
                      lambda_min_ratio = opt$lambda_min_ratio,
                      nlambda = opt$nlambda,
                      pmax = nterms_max + 1, pmax_strict = FALSE,
                      weights = d_train$weights,
                      intercept = intercept, obsvar = v, penalty = penalty,
                      thresh = opt$thresh)

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

  out$solution_terms <- order[seq_len(nterms_max)]
  if (any(is.na(out$solution_terms)) &&
      length(entered_variables) < nterms_max) {
    if (length(setdiff(notentered_variables,
                       which(penalty == Inf))) > 0) {
      warning("Less than nterms_max variables entered L1-path. ",
              "Try reducing lambda_min_ratio. ")
    }
  }

  return(out)
}

search_L1 <- function(p_ref, refmodel, nterms_max, penalty, opt) {
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
    intercept = TRUE, ncol(x), penalty, opt
  )

  solution_terms_orig <- collapse_ranked_predictors(
    path = colnames(x)[search_path$solution_terms], formula = refmodel$formula,
    data = fr
  )
  solution_terms <- utils::head(solution_terms_orig, nterms_max)
  # Check for interaction terms being selected before all involved main-effect
  # terms have been selected (and reorder `solution_terms` if that is the case):
  ia_terms <- grep(":", solution_terms, value = TRUE)
  stopifnot(!any(duplicated(ia_terms))) # safety measure for which.max()
  for (ia_term in ia_terms) {
    idx_ia <- which.max(solution_terms == ia_term)
    if (idx_ia > nterms_max) break
    main_terms_ia <- strsplit(ia_term, ":")[[1]]
    main_terms_ia <- intersect(main_terms_ia, solution_terms_orig)
    prev_terms <- utils::head(solution_terms, idx_ia - 1L)
    ia_sel_bef_main <- !all(main_terms_ia %in% prev_terms)
    if (ia_sel_bef_main) {
      if (getOption("projpred.warn_L1_interactions", TRUE)) {
        warning(
          "Interaction term `", ia_term, "` was selected before all involved ",
          "main-effect terms have been selected. This is a known deficiency ",
          "of L1 search. Use forward search to avoid this. Now modifying the ",
          "predictor ranking such that this interaction term comes after the ",
          "main-effect terms involved in it."
        )
      }
      main_terms_ia <- main_terms_ia[order(match(main_terms_ia,
                                                 solution_terms_orig))]
      new_head <- c(prev_terms, setdiff(main_terms_ia, prev_terms), ia_term)
      solution_terms <- c(new_head, setdiff(solution_terms, new_head))
      solution_terms <- utils::head(solution_terms, nterms_max)
    }
  }

  outdmins <- lapply(0:length(solution_terms), function(nterms) {
    if (nterms == 0) {
      formula <- make_formula(c("1"))
      beta <- NULL
      x <- x[, numeric(), drop = FALSE]
    } else {
      formula <- make_formula(solution_terms[seq_len(nterms)])
      variables <- unlist(lapply(
        solution_terms[seq_len(nterms)],
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
      indices <- match(variables, colnames(x)[search_path$solution_terms])
      indices <- indices[!is.na(indices)]
      beta <- search_path$beta[indices, max(indices) + 1, drop = FALSE]
      # Also reduce `x` (important for coef.subfit(), for example); note that
      # `x <- x[, variables, drop = FALSE]` should also be possible, but the
      # re-use of `colnames(x)` should provide another sanity check:
      x <- x[, colnames(x)[search_path$solution_terms[indices]], drop = FALSE]
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
                 w = search_path$w[, nterms + 1], formula, x, xlvls)
    class(sub) <- "subfit"
    return(list(sub))
  })

  return(nlist(solution_terms, outdmins))
}
