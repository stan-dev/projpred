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
	
	# prediction of full model (integrated over uncertainty about f)
	mu <- rowMeans(family$linkinv(f))
	
	# L1-penalized projection (projection path)
	search <- glm_elnet(x,mu,family)
	
	if (relax) {
		# sort the variables according to the order in which they enter
		# the model in the L1-path
		order <- sort(apply(search$beta!=0,1,function(num) which(num)[1]), index.return=TRUE)$ix;
		
		# compute the relaxed (unregularized) projection for the selected models
		relax <- lapply(1:length(order), function(j){
			out <- glm_ridge(x[,order(1:j)],mu,family)
			return(list(beta=out$beta, beta0=out$beta0))
		})
		
	}
	
	stop('Not finished yet.')
	return(proj)
	
}























