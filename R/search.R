#' The forward selection algorithm
#'
#' p_full contains the parameters of the full model, d_train contains the data
#' and p_sub contains the parameters of the submodel. b0 is a guess for the
#' initial weight vector in the submodel. In addition, if a subset of the
#' samples is used for the variable selection, p_clust contains parameters
#' of the full model using only those samples.

fsel <- function(p_full, d_train, b0, args) {

  # initialize the forward selection
  # proj performs the projection over samples
  proj <- .get_proj_handle(args$family_kl)
  i <- 1
  iq <- ceiling(quantile(1:args$nv, 1:10/10))
  cols <- 1:ncol(d_train$x)
  chosen <- NULL

  # if the full model has intercept, add it to the submodel 'as the first variable'
  if(args$intercept) {
    chosen <- 1
    i <- i + 1
  }

  # start adding variables one at a time
  while(i <= args$nv + args$intercept) {

    notchosen <- setdiff(cols, chosen)

    candidates <- lapply(notchosen, function(x) c(chosen, x))
    p_sub <- sapply(candidates, proj, p_full, d_train, b0, args)
    imin <- which.min(p_sub['kl',])
    chosen <- c(chosen, notchosen[imin])

    if(args$verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))

    i <- i + 1
  }

  chosen
}



search_L1 <- function(p_full, d_train, b0, args) {
	
	# prediction of full model (integrate over uncertainty about f)
	mu <- rowMeans(p_full$mu)
	
	# create a grid of lambda values
	lambda_min_ratio <- 1e-3 # this should be small enough so that the computation does not stop before pmax
	nlam <- 100 # should be dense enough
	lambda <- lambda_grid(d_train$x, mu, args$family_kl, alpha=1.0, eps=lambda_min_ratio, nlam=nlam)
	# add a few very small lambda values to proceed the variable selection almost up to the full model
	lambda <- c(lambda, rev(seq(log(1e-15), log(min(lambda)),  by=log(10))))
	
	
	# L1-penalized projection (projection path)
	pmax <- dim(d_train$x)[2] # TODO: THIS SHOULD BE GIVEN AS AN INPUT
	search <- glm_elnet(d_train$x, mu, args$family_kl, lambda=lambda, pmax=pmax, pmax_strict=FALSE,
						offset=d_train$offset, weights=d_train$weights)
	
	# sort the variables according to the order in which they enter the model in the L1-path
	entering_indices <- apply(search$beta!=0, 1, function(num) which(num)[1])
	order <- sort(entering_indices, index.return=TRUE)$ix
	
	print(entering_indices)
	return(order)
}























