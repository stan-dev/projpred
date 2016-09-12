#' The forward selection algorithm
#'
#' p_full contains the parameters of the full model, d_train contains the data
#' and p_sub contains the parameters of the submodel. b0 is a guess for the
#' initial weight vector in the submodel. In addition, if a subset of the
#' samples is used for the variable selection, p_clust contains parameters
#' of the full model using only those samples.

fsel <- function(p_full, d_train, d_test, p_clust = NULL, b0, args) {

  # initialize the forward selection
  # proj performs the projection over samples
  proj <- .get_proj_handle(args$family_kl)
  i <- 1
  iq <- ceiling(quantile(1:args$nv, 1:10/10))
  cols <- 1:ncol(d_train$x)
  chosen <- NULL
  kl <- rep(NA_real_, args$nv)
  mu <- vector('list', args$nv)
  lppd <- vector('list', args$nv)
  p_fsel <- if(args$clust) p_clust else p_full

  # if the full model has intercept, add it to the submodel 'as the first variable'
  if(args$intercept) {
    chosen <- 1
    p_sub <- proj(NULL, chosen, p_full, d_train, b0, args)
    kl[i] <- p_sub$kl
    sub_stats_temp <- .summary_stats(d_test, chosen, p_sub, args)
    mu[[i]] <- sub_stats_temp$mu
    lppd[[i]] <- sub_stats_temp$lppd
    i <- i + 1
  }

  # start adding variables one at a time
  while(i <= args$nv) {

    notchosen <- setdiff(cols, chosen)

    p_sub <- sapply(notchosen, proj, chosen, p_fsel, d_train, b0, args)
    imin <- which.min(p_sub['kl',])

    if(p_sub[['kl',imin]] == Inf) {
      warning(paste0('Numerical problems in the projection for a submodel of size ',
                     i, '. Ending forward selection.'))
      return(list(chosen = chosen[1:(i-1)], kl = kl[1:(i-1)],
                  mu = mu[1:(i-1)], lppd = lppd[1:(i-1)]))
    }

    # p_sel is the selected submodel
    p_sel <- if(args$clust) proj(notchosen[imin], chosen, p_full, d_train, b0, args) else p_sub[,imin]
    chosen <- c(chosen, notchosen[imin])

    kl[i] <- p_sel$kl
    sub_stats_temp <- .summary_stats(d_test, chosen, p_sel, args)
    mu[[i]] <- sub_stats_temp$mu
    lppd[[i]] <- sub_stats_temp$lppd

    if(args$verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))
    i <- i + 1
  }

  list(chosen = chosen, kl = kl, mu = mu, lppd = lppd)
}
