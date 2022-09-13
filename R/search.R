search_forward <- function(p_ref, refmodel, nterms_max, verbose = TRUE, opt,
                           search_terms = NULL, ...) {
  iq <- ceiling(quantile(seq_len(nterms_max), 1:10 / 10))
  if (is.null(search_terms)) {
    allterms <- split_formula(refmodel$formula, data = refmodel$fetch_data())
  } else {
    allterms <- search_terms
  }

  chosen <- character()
  total_terms <- count_terms_chosen(allterms)
  stop_search <- min(total_terms, nterms_max)
  submodls <- c()

  for (size in seq_len(stop_search)) {
    cands <- select_possible_terms_size(chosen, allterms, size = size)
    if (is.null(cands))
      next
    full_cands <- lapply(cands, function(cand) c(chosen, cand))
    subL <- lapply(full_cands, project_submodel,
                   p_ref = p_ref, refmodel = refmodel, regul = opt$regul, ...)

    ## select best candidate
    imin <- which.min(sapply(subL, "[[", "kl"))
    chosen <- c(chosen, cands[imin])

    ## append `submodl`
    submodls <- c(submodls, list(subL[[imin]]$submodl))

    if (verbose && count_terms_chosen(chosen) %in% iq) {
      print(paste0(names(iq)[max(which(count_terms_chosen(chosen) == iq))],
                   " of terms selected."))
    }
  }

  # For `solution_terms`, `reduce_models(chosen)` used to be used instead of
  # `chosen`. However, `reduce_models(chosen)` and `chosen` should be identical
  # at this place because select_possible_terms_size() already avoids redundant
  # models. Thus, use `chosen` here because it matches `submodls` (this
  # matching is necessary because later in .get_submodels()'s `!refit_prj` case,
  # `submodls` is indexed with integers which are based on `solution_terms`):
  return(nlist(solution_terms = setdiff(chosen, "1"), submodls))
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

  out$solution_terms <- order[1:nterms_max]
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
         "reference model.")
  }
  # TODO: In the following model.matrix() call, allow user-specified contrasts
  # to be passed to argument `contrasts.arg`. The `contrasts.arg` default
  # (`NULL`) uses `options("contrasts")` internally, but it might be more
  # convenient to let users specify contrasts directly. At that occasion,
  # contrasts should also be tested thoroughly (not done until now).
  x <- model.matrix(refmodel$formula, data = refmodel$fetch_data())
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  ## it's important to keep the original order because that's the order
  ## in which lasso will estimate the parameters
  tt <- terms(refmodel$formula)
  terms_ <- attr(tt, "term.labels")
  search_path <- search_L1_surrogate(
    p_ref, nlist(x, weights = refmodel$wobs), refmodel$family,
    refmodel$intercept, ncol(x), penalty, opt
  )
  solution_terms <- collapse_contrasts_solution_path(
    refmodel$formula, colnames(x)[search_path$solution_terms],
    refmodel$fetch_data()
  )
  submodls <- lapply(0:length(solution_terms), function(nterms) {
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
          mm <- model.matrix(as.formula(paste("~ 1 +", term)),
                             data = refmodel$fetch_data())
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
    }
    sub <- nlist(alpha = search_path$alpha[nterms + 1],
                 beta,
                 w = search_path$w[, nterms + 1],
                 formula,
                 x)
    class(sub) <- "subfit"
    return(list(sub))
  })
  return(list(solution_terms = solution_terms[seq_len(nterms_max)],
              submodls = submodls[seq_len(nterms_max + 1)]))
}
