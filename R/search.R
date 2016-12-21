#' The forward selection algorithm
#'
#' p_full contains the parameters of the full model, d_train contains the data
#' and p_sub contains the parameters of the submodel. b0 is a guess for the
#' initial weight vector in the submodel. In addition, if a subset of the
#' samples is used for the variable selection, p_clust contains parameters
#' of the full model using only those samples.

fsel <- function(p_full, d_train, family_kl, intercept, nvmax, regul, coef_init, verbose) {

  # initialize the forward selection
  # proj performs the projection over samples
  proj <- .get_proj_handle(family_kl)
  i <- 1
  iq <- ceiling(quantile(1:nvmax, 1:10/10))
  cols <- 1:ncol(d_train$x)
  chosen <- NULL

  # start adding variables one at a time
  while(i <= nvmax) {

    notchosen <- setdiff(cols, chosen)
    cands <- lapply(notchosen, function(x) c(chosen, x))

    p_sub <- sapply(cands, proj, p_full, d_train, intercept, regul, coef_init)

    imin <- which.min(p_sub['kl',])
    chosen <- c(chosen, notchosen[imin])

    if(verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))

    i <- i + 1
  }

  chosen
}



search_L1 <- function(p_full, d_train, family, intercept, nvmax) {
	
	# prediction of full model (integrate over the uncertainty about the model parameters)
	mu <- p_full$mu %*% p_full$weights
	
	# create a grid of lambda values
	lambda_min_ratio <- 1e-3 # this should be small enough so that the computation does not stop before pmax
	nlam <- 100 # should be dense enough
	lambda <- lambda_grid(d_train$x, mu, family, alpha=1.0, eps=lambda_min_ratio, nlam=nlam)
	# add a few very small lambda values to proceed the variable selection almost up to the full model
	lambda <- c(lambda, rev(seq(log(1e-15), log(min(lambda)),  by=log(10))))
	
	
	# L1-penalized projection (projection path)
	search <- glm_elnet(d_train$x, mu, family, lambda=lambda, pmax=nvmax, pmax_strict=FALSE,
						offset=d_train$offset, weights=d_train$weights, intercept=intercept)
	
	# sort the variables according to the order in which they enter the model in the L1-path
	entering_indices <- apply(search$beta!=0, 1, function(num) which(num)[1])
	order <- sort(entering_indices, index.return=TRUE)$ix[1:nvmax]
	
	return(order)
}























