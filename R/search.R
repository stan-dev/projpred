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
  
  # L1-penalized projection (projection path).
  # (Notice: here we use pmax = nv_max+1 so that the computation gets carried until all the way 
  # down to the least regularization also for model size nv_max)
  search <- glm_elnet(d_train$x, mu, family, lambda_min_ratio=opt$lambda_min_ratio, nlambda=opt$nlambda,
                      pmax=nv_max+1, pmax_strict=FALSE, offset=d_train$offset, weights=d_train$weights,
                      intercept=intercept, obsvar=v, penalty=penalty, thresh=opt$thresh)
  
  # sort the variables according to the order in which they enter the model in the L1-path
  entering_indices <- apply(search$beta!=0, 1, function(num) which(num)[1]) # na for those that did not enter
  entered_variables <- c(1:NCOL(d_train$x))[!is.na(entering_indices)] # variables that entered at some point
  notentered_variables <- c(1:NCOL(d_train$x))[is.na(entering_indices)] # variables that did not enter at any point
  order_of_entered <- sort(entering_indices, index.return=TRUE)$ix
  order <- c(entered_variables[order_of_entered], notentered_variables)
  
  # fetch the coefficients corresponding to those points at the searchpath where new variable enters
  nvar <- length(order)
  n <- nrow(p_full$mu)
  out <- list(alpha=rep(NA, nv_max+1), beta=matrix(0, nrow=nv_max, ncol=nv_max+1), 
              w=matrix(NA, nrow=n, ncol=nv_max+1))
  for (k in 0:nv_max) {
    if (k == 0) {
      out$alpha[1] <- search$beta0[1]
      out$w[,1] <- search$w[,1]
    } else {
      # find those points in the L1-path where only the k most relevant features can have nonzero
      # coefficient, and then fetch their coefficients with least regularization
      ivar <- utils::tail(order, nvar-k)
      steps_k_var <- which(colSums(search$beta[ivar,,drop=F] != 0) == 0)
      if (length(steps_k_var) > 0) 
        j <- utils::tail(steps_k_var, 1)
      else 
        # no steps where all the variables in set ivar would have zero coefficient (could be due
        # to one or more of these variables having penalty = 0 so they are always in the model)
        # so set the coefficients to be equal to the starting value
        j <- 1
      out$alpha[k+1] <- search$beta0[j]
      out$beta[1:k,k+1] <- search$beta[order[1:k],j]
      out$w[,k+1] <- search$w[,j]
    }
  }
  
  if (length(entered_variables) < nv_max)
    if (length(setdiff(notentered_variables, which(penalty == Inf))) > 0)
	    warning('Less than nv_max variables entered L1-path. Try reducing lambda_min_ratio. ')
	
  
  out$vind <- order[1:nv_max]
	return(out)
}

