#' Search heuristics
#'

search_forward <- function(p_full, d_train, family_kl, intercept, nv_max,
                           verbose) {

  # initialize the forward selection
  # proj performs the projection over samples
  proj <- .get_proj_handle(family_kl)
  i <- 1
  iq <- ceiling(quantile(1:nv_max, 1:10/10))
  cols <- 1:ncol(d_train$x)
  chosen <- NULL

  # start adding variables one at a time
  while(i <= nv_max) {

    notchosen <- setdiff(cols, chosen)
    cands <- lapply(notchosen, function(x) c(chosen, x))

    p_sub <- sapply(cands, proj, p_full, d_train, intercept)

    imin <- which.min(p_sub['kl',])
    chosen <- c(chosen, notchosen[imin])

    if(verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))

    i <- i + 1
  }

  chosen
}



search_L1 <- function(p_full, d_train, family, intercept, nv_max) {

    # prediction of full model (integrate over the uncertainty about the model parameters)
    mu <- p_full$mu %*% p_full$weights

    # create a grid of lambda values
    lambda_min_ratio <- 1e-3 # this is the relative value above which we use fairly dense grid
    nlam <- 100 # how many values in the grid before lambda_min_ratio
    lambda <- lambda_grid(d_train$x, mu, family, alpha=1.0, eps=lambda_min_ratio, nlam=nlam)

    # add a sparser grid of very small lambda values to proceed the variable selection almost up to the full model
    lambda_min <- 1e-15 # ultimate minimum for lambda
    lambda <- c(lambda, rev(seq(log(lambda_min), log(min(lambda)),  by=log(10))))


    # L1-penalized projection (projection path)
    search <- glm_elnet(d_train$x, mu, family, lambda=lambda, pmax=nv_max, pmax_strict=FALSE,
                        offset=d_train$offset, weights=d_train$weights, intercept=intercept)

    # sort the variables according to the order in which they enter the model in the L1-path
    entering_indices <- apply(search$beta!=0, 1, function(num) which(num)[1])
    entered_variables <- c(1:NCOL(d_train$x))[!is.na(entering_indices)] # variables that entered at some point
    order_of_entered <- sort(entering_indices, index.return=TRUE)$ix[1:nv_max]

    return(entered_variables[order_of_entered])
}























