# Search heuristics
#

search_forward1 <- function(p_full, d_train, family, intercept, nv_max,
                           verbose, opt) {
  
  # predictive mean and variance of the reference model (with parameters integrated out)
  mu <- p_full$mu 
  v <- p_full$var
  
  if (NCOL(mu) > 1 || NCOL(v) > 1)
    stop('Internal error: search_forward1 received multiple draws. Please report to the developers.')
  
  # forward search
  search <- glm_forward(d_train$x, mu, family, lambda=opt$regul, offset=d_train$offset, weights=d_train$weights,
                        obsvar=v, intercept=intercept, pmax=nv_max)

  return(search$varorder)
}

search_forward <- function(p_full, d_train, family_kl, intercept, nv_max,
                           verbose, opt) {

  # initialize the forward selection
  # proj performs the projection over samples
  projfun <- .get_proj_handle(family_kl, opt$regul)
  i <- 1
  iq <- ceiling(quantile(1:nv_max, 1:10/10))
  cols <- 1:ncol(d_train$x)
  chosen <- NULL

  # start adding variables one at a time
  while(i <= nv_max) {

    notchosen <- setdiff(cols, chosen)
    cands <- lapply(notchosen, function(x) c(chosen, x))

    p_sub <- sapply(cands, projfun, p_full, d_train, intercept)

    imin <- which.min(p_sub['kl',])
    chosen <- c(chosen, notchosen[imin])

    if(verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))

    i <- i + 1
  }

  chosen
}



search_L1 <- function(p_full, d_train, family, intercept, nv_max, penalty, opt) {
  
  # predictive mean and variance of the reference model (with parameters integrated out)
  mu <- p_full$mu
  v <- p_full$var
  
  if (NCOL(mu) > 1 || NCOL(v) > 1)
    stop('Internal error: search_L1 received multiple draws. Please report to the developers.')
  
  # L1-penalized projection (projection path)
  search <- glm_elnet(d_train$x, mu, family, lambda_min_ratio=opt$lambda_min_ratio, nlambda=opt$nlambda,
                      pmax=nv_max, pmax_strict=FALSE, offset=d_train$offset, weights=d_train$weights,
                      intercept=intercept, obsvar=v, penalty=penalty)
  
  # sort the variables according to the order in which they enter the model in the L1-path
  entering_indices <- apply(search$beta!=0, 1, function(num) which(num)[1]) # na for those that did not enter
  entered_variables <- c(1:NCOL(d_train$x))[!is.na(entering_indices)] # variables that entered at some point
  notentered_variables <- c(1:NCOL(d_train$x))[is.na(entering_indices)] # variables that did not enter at any point
  order_of_entered <- sort(entering_indices, index.return=TRUE)$ix
  order <- c(entered_variables[order_of_entered], notentered_variables)
  
  
  if (length(entered_variables) < nv_max)
    if (length(setdiff(notentered_variables, which(penalty == Inf))) > 0)
	    warning('Less than nv_max variables entered L1-path. Try reducing lambda_min_ratio. ')
	
	return(order[1:nv_max])
}

