search_rsens <- function(p_ref, refmodel, family, intercept, nterms_max,
                         verbose = TRUE, opt, search_terms = NULL,
                         increasing_order = TRUE) {
  ## 1. get possible terms
  ## 2. if there are no group terms, do rsens on the single variables level
  ## 3. if there are group terms, do rsens on interactions as well
  ## 4. order terms present in the formula according to our rules
  ndraws <- NCOL(p_ref$mu)

  formula <- refmodel$formula
  terms <- split_formula(formula,
    data = refmodel$fetch_data(),
    add_main_effects = FALSE
  )
  if (formula_contains_additive_terms(formula)) {
    stop("Rsens search for additive models is not implemented yet.")
  }
  ranks <- rankvars::rank(refmodel$fit, ndraws = ndraws, summary_type = "both")

  ordering <- bind_cols(
    as_data_frame(ranks$variables),
    as_data_frame(ranks$interactions)
  ) %>%
    as_data_frame() %>%
    gather() %>%
    filter(value > 0) %>%
    group_by(key) %>%
    summarise(value = mean(value)) %>%
    pivot_wider(names_from = "key", values_from = "value") %>%
    sort() %>%
    names()

  if (!formula_contains_group_terms(formula)) {
    return(list(solution_terms = ordering))
  }

  group_factors <- extract_terms_response(formula) %>%
    .$group_terms %>%
    strsplit(., "[ ]*\\|([^\\|]*\\||)[ ]*") %>%
    lapply(FUN = function(x) x[2]) %>%
    unlist() %>%
    unique()

  for (g in group_factors) {
    ordering <- gsub(paste0("^", g, "$"), paste0("(1 | ", g, ")"), ordering)
    ordering <- gsub(paste0("([\\w\\d.]+):", g),
      paste0("(\\1 | ", g, ")"),
      ordering,
      perl = TRUE
    )
  }

  return(list(solution_terms = ordering))
}

search_forward <- function(p_ref, refmodel, family, intercept, nterms_max,
                           verbose = TRUE, opt, search_terms = NULL,
                           increasing_order = TRUE) {
  ## initialize the forward selection
  ## proj performs the projection over draws
  projfun <- .get_proj_handle(refmodel, p_ref, family, opt$regul, intercept)

  formula <- refmodel$formula
  iq <- ceiling(quantile(seq_len(nterms_max), 1:10 / 10))
  if (is.null(search_terms)) {
    allterms <- split_formula(formula, data = refmodel$fetch_data())
  } else {
    allterms <- search_terms
  }

  chosen <- NULL
  total_terms <- count_terms_chosen(allterms)
  stop_search <- min(total_terms, nterms_max)
  submodels <- c()

  for (size in seq_len(stop_search)) {
    cands <- select_possible_terms_size(chosen, allterms, size = size)
    if (is.null(cands))
      break
    full_cands <- lapply(cands, function(cand) c(chosen, cand))
    sub <- sapply(full_cands, projfun)

    ## select best candidate
    imin <- which.min(sapply(seq_len(NCOL(sub)), function(i) {
      min(unlist(sub["kl", i]))
    }))
    chosen <- c(chosen, cands[imin])

    ## append submodels
    submodels <- c(submodels, sub["sub_fit", imin])

    if (verbose && length(chosen) %in% iq) {
      print(paste0(names(iq)[max(which(length(chosen) == iq))],
                   " of terms selected."))
    }
  }

  ## reduce chosen to a list of non-redundant accumulated models
  return(list(solution_terms = setdiff(reduce_models(chosen), "1"),
              sub_fits = submodels))
}

# copied over from search until we resolve the TODO below
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
    lambda_min_ratio = opt$lambda_min_ratio, nlambda = opt$nlambda,
    pmax = nterms_max + 1, pmax_strict = FALSE, offset = d_train$offset,
    weights = d_train$weights, intercept = intercept, obsvar = v,
    penalty = penalty, thresh = opt$thresh)

  ## sort the variables according to the order in which they enter the model in
  ## the L1-path
  entering_indices <- apply(search$beta != 0, 1, function(num)
    which(num)[1]) # na for those that did not enter
  ## variables that entered at some point
  entered_variables <-
    c(seq_len(NCOL(d_train$x)))[!is.na(entering_indices)]
  ## variables that did not enter at any point
  notentered_variables <-
    c(seq_len(NCOL(d_train$x)))[is.na(entering_indices)]
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

search_L1 <- function(p_ref, refmodel, family, intercept, nterms_max, penalty,
                      opt) {
  frame <- model.frame(refmodel$formula, refmodel$fetch_data())
  contrasts_arg <- get_contrasts_arg_list(
    refmodel$formula,
    refmodel$fetch_data()
  )
  x <- model.matrix(delete.intercept(refmodel$formula),
    data = frame,
    contrasts.arg = contrasts_arg
  )
  ## it's important to keep the original order because that's the order
  ## in which lasso will estimate the parameters
  tt <- terms(refmodel$formula)
  terms_ <- attr(tt, "term.labels")
  search_path <- search_L1_surrogate(
    p_ref, list(refmodel, x = x), family,
    intercept, ncol(x), penalty, opt
  )
  solution_terms <- collapse_contrasts_solution_path(
    refmodel$formula, colnames(x)[search_path$solution_terms],
    refmodel$fetch_data()
  )
  sub_fits <- lapply(0:length(solution_terms), function(nterms) {
    if (nterms == 0) {
      formula <- make_formula(c("1"))
      beta <- NULL
    } else {
      formula <- make_formula(solution_terms[seq_len(nterms)])
      variables <- unlist(lapply(
        solution_terms[seq_len(nterms)],
        function(term) {
          form <- as.formula(paste("~ 0 +", term))
          contrasts_arg <- get_contrasts_arg_list(
            form,
            refmodel$fetch_data()
          )
          return(colnames(model.matrix(form,
            data = refmodel$fetch_data(),
            contrasts.arg = contrasts_arg
          )))
        }
      ))
      indices <- match(variables, colnames(x)[search_path$solution_terms])
      beta <- search_path$beta[indices, max(indices) + 1, drop = FALSE]
    }
    sub <- nlist(
      alpha = search_path$alpha[nterms + 1],
      beta,
      w = search_path$w[, nterms + 1],
      formula,
      x
    )
    class(sub) <- "subfit"
    return(sub)
  })
  return(list(
    solution_terms = solution_terms[seq_len(nterms_max)],
    sub_fits = sub_fits[seq_len(nterms_max + 1)]
  ))
}
