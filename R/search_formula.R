search_forward_poc <- function(p_ref, refmodel, family_kl, intercept, nv_max,
                               verbose, opt, groups=NULL, increasing_order=TRUE) {
  ## initialize the forward selection
  ## proj performs the projection over samples
  projfun <- .get_proj_handle_poc(family_kl, opt$regul)
  i <- 1
  iq <- ceiling(quantile(1:nv_max, 1:10/10))
  if (is.null(groups))
    terms_ <- break_formula(refmodel$formula)
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
  total_variables <- count_terms_in_submodel(refmodel$formula)
  submodels <- c()

  ## start adding terms one at a time
  while(count_variables_chosen(refmodel, reduce_models(refmodel, chosen)) < nv_max
        & count_variables_chosen(refmodel, chosen) < total_variables) {

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
          if (is.redundant(refmodel, chosen, x)) x else NA)

        ## if redundant models add the terms to the list so we don't iterate forever
        chosen <- c(chosen, unname(unlist(already_selected[!is.na(already_selected)])))
      }
    }

    ## only add candidates that are not redundant with previous chosen submodels
    cands <- lapply(notchosen, function(x)
      if (is.redundant(refmodel, chosen, x)) NA else c(chosen, x))

    ## remove already selected terms
    cands <- cands[!is.na(cands)]

    p_sub <- sapply(cands, projfun, p_ref, refmodel, intercept)

    imin <- which.min(p_sub['kl', ])
    chosen <- c(chosen, notchosen[imin])

    ## append submodels
    submodels <- c(submodels, p_sub['sub_fit', imin])

    if(verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))

    i <- i + 1
  }

  ## reduce chosen to a list of non-redundant accumulated models
  list(vind=reduce_models(refmodel, chosen), sub_fits=submodels)
}

search_L1_poc <- function(p_ref, refmodel, family, intercept, nv_max, penalty, opt) {
  terms_ <- break_formula(refmodel$formula)
  x <- model.matrix(refmodel$formula, refmodel$fetch_data())
  spath <- search_L1(p_ref, list(refmodel, x=x), family, intercept, nv_max, penalty, opt)
  ## TODO: check this? can we reduce the above line and this one to a single thing?
  ## extract the path from glmnet
  sub_fit <- glmnet::glmnet(x, p_ref$mu, family, intercept = intercept, dfmax = nv_max, lambda = penalty)
  return(list(vind=terms_[spath$vind - 1], sub_fit=sub_fit))
}
