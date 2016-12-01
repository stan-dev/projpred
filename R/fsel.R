#' The forward selection algorithm
#'
#' p_full contains the parameters of the full model, d_train contains the data
#' and p_sub contains the parameters of the submodel. b0 is a guess for the
#' initial weight vector in the submodel. In addition, if a subset of the
#' samples is used for the variable selection, p_clust contains parameters
#' of the full model using only those samples.

fsel <- function(p_full, d_train, family_kl, intercept, nv, regul, coef_init,
                 verbose) {

  # initialize the forward selection
  # proj performs the projection over samples
  proj <- .get_proj_handle(family_kl)
  i <- 1
  iq <- ceiling(quantile(1:nv, 1:10/10))
  cols <- 1:ncol(d_train$x)
  chosen <- NULL

  # start adding variables one at a time
  while(i <= nv) {

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
