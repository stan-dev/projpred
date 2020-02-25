search_forward_poc <- function(p_ref, refmodel, family_kl, intercept, nv_max,
                               verbose=TRUE, opt, groups=NULL, increasing_order=TRUE) {
  ## initialize the forward selection
  ## proj performs the projection over samples
  projfun <- .get_proj_handle_poc(family_kl, opt$regul)

  formula <- refmodel$formula
  i <- 1
  iq <- ceiling(quantile(1:nv_max, 1:10/10))
  if (is.null(groups))
    terms_ <- split_formula(formula)
  else
    terms_ <- groups
  if (increasing_order) {
    terms_ <- sort_submodels_by_size(unname(unlist(terms_)))
    order <- 1:length(terms_)
    current <- order[1]
    current_terms <- terms_[[current]]
  } else {
    order <- NULL
    current_terms <- terms_
  }

  chosen <- NULL
  total_variables <- count_terms_in_subformula(formula)
  submodels <- c()

  ## start adding terms one at a time
  while(count_variables_chosen(reduce_models(chosen)) < nv_max
        & count_variables_chosen(chosen) < total_variables) {

    notchosen <- setdiff(current_terms, chosen)

    ## if we have included all submodels in a class start with next class
    ## this only happens with multilevel or interaction models, where some terms may include
    ## more than one variable.
    ## In GLMs every term represents a single variable and therefore all terms are within
    ## size==1
    if (length(notchosen) == 0 & !is.null(order)) {
      if (current < order[length(order)]) {
        current <- current + 1
        current_terms <- terms_[[current]]
        notchosen <- setdiff(current_terms, chosen)

        already_selected <- lapply(notchosen, function(x)
          if (is_next_submodel_redundant(chosen, x)) x else NA)

        ## if redundant models add the terms to the list so we don't iterate forever
        chosen <- c(chosen, unname(unlist(already_selected[!is.na(already_selected)])))
      }
    }

    ## only add candidates that are not redundant with previous chosen submodels
    cands <- lapply(notchosen, function(x)
      if (is_next_submodel_redundant(chosen, x)) NA else c(chosen, x))

    ## remove already selected terms
    cands <- cands[!is.na(cands)]

    p_sub <- sapply(cands, projfun, p_ref, refmodel, intercept)

    imin <- which.min(p_sub['kl', ])
    chosen <- c(chosen, notchosen[imin])

    ## append submodels
    submodels <- c(submodels, p_sub['sub_fit', imin])

    if (verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))

    i <- i + 1
  }

  ## reduce chosen to a list of non-redundant accumulated models
  list(vind=reduce_models(chosen), sub_fits=submodels)
}

#' copied over from search until we resolve the TODO below
search_L1 <- function(p_ref, d_train, family, intercept, nv_max, penalty, opt) {

  ## predictive mean and variance of the reference model (with parameters integrated out)
  mu <- p_ref$mu
  v <- p_ref$var

  if (NCOL(mu) > 1 || NCOL(v) > 1)
    stop('Internal error: search_L1 received multiple draws. Please report to the developers.')

  ## L1-penalized projection (projection path).
  ## (Notice: here we use pmax = nv_max+1 so that the computation gets carried until all the way
  ## down to the least regularization also for model size nv_max)
  search <- glm_elnet(d_train$x, mu, family, lambda_min_ratio=opt$lambda_min_ratio, nlambda=opt$nlambda,
                      pmax=nv_max+1, pmax_strict=FALSE, offset=d_train$offset, weights=d_train$weights,
                      intercept=intercept, obsvar=v, penalty=penalty, thresh=opt$thresh)

  ## sort the variables according to the order in which they enter the model in the L1-path
  entering_indices <- apply(search$beta!=0, 1, function(num) which(num)[1]) # na for those that did not enter
  entered_variables <- c(1:NCOL(d_train$x))[!is.na(entering_indices)] # variables that entered at some point
  notentered_variables <- c(1:NCOL(d_train$x))[is.na(entering_indices)] # variables that did not enter at any point
  order_of_entered <- sort(entering_indices, index.return=TRUE)$ix
  order <- c(entered_variables[order_of_entered], notentered_variables)

  ## fetch the coefficients corresponding to those points at the searchpath where new variable enters
  nvar <- length(order)
  n <- nrow(p_ref$mu)
  out <- list(alpha=rep(NA, nv_max+1), beta=matrix(0, nrow=nv_max, ncol=nv_max+1),
              lambda=rep(NA, nv_max+1), w=matrix(NA, nrow=n, ncol=nv_max+1))
  for (k in 0:nv_max) {
    if (k == 0) {
      out$alpha[1] <- search$beta0[1]
      out$lambda[1] <- search$lambda[1]
      out$w[,1] <- search$w[,1]
    } else {
      ## find those points in the L1-path where only the k most relevant features can have nonzero
      ## coefficient, and then fetch their coefficients with least regularization
      ivar <- utils::tail(order, nvar-k)
      steps_k_var <- which(colSums(search$beta[ivar,,drop=F] != 0) == 0)
      if (length(steps_k_var) > 0)
        j <- utils::tail(steps_k_var, 1)
      else
        ## no steps where all the variables in set ivar would have zero coefficient (could be due
        ## to one or more of these variables having penalty = 0 so they are always in the model)
        ## so set the coefficients to be equal to the starting value
        j <- 1
      out$alpha[k+1] <- search$beta0[j]
      out$beta[1:k,k+1] <- search$beta[order[1:k],j]
      out$lambda[k+1] <- search$lambda[j]
      out$w[,k+1] <- search$w[,j]
    }
  }

  if (length(entered_variables) < nv_max)
    if (length(setdiff(notentered_variables, which(penalty == Inf))) > 0)
      warning('Less than nv_max variables entered L1-path. Try reducing lambda_min_ratio. ')

  out$vind <- order[1:nv_max]
  return(out)
}

search_L1_poc <- function(p_ref, refmodel, family, intercept, nv_max, penalty, opt) {
  terms_ <- split_formula(refmodel$formula)
  x <- model.matrix(refmodel$formula, refmodel$fetch_data())
  spath <- search_L1(p_ref, list(refmodel, x=x), family, intercept, nv_max, penalty, opt)
  ## TODO: check this? can we reduce the above line and this one to a single thing?
  ## extract the path from glmnet
  sub_fit <- glmnet::glmnet(x, p_ref$mu, family, intercept = intercept, dfmax = nv_max, lambda = penalty)
  return(list(vind=terms_[spath$vind - 1], sub_fit=sub_fit))
}
