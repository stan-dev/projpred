#' The forward selection algorithm
#'
#' p_full contains the parameters of the full model, d_train contains the data
#' and p_sub contains the parameters of the submodel. b0 is a guess for the
#' initial weight vector in the submodel. In addition, if a subset for the
#' variable selection, p_clust contains parameters of the full model using
#' only those samples.


fsel <- function(p_full, d_train, d_test, p_clust = NULL, b0, args) {

  # for gaussian, precompute some values
  if(args$family_kl$family == 'gaussian') {
    d_train$covX <- crossprod(d_train$x)
    p_full$x_mu <- crossprod(d_train$x, p_full$mu)
    p_full$dis2 <- p_full$dis^2
    if(args$clust) {
      p_clust$x_mu <- crossprod(d_train$x,p_clust$mu)
      p_clust$dis2 <- p_clust$dis^2
    }
  }

  # initialize the forward selection
  # proj performs the projection over samples
  proj <- .get_proj_handle(args$family_kl$family)
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
    mu_temp <- args$family_kl$linkinv(d_test$x[,chosen]%*%p_sub$b)
    mu[[i]] <- rowMeans(mu_temp)
    lppd[[i]] <- apply(args$family_kl$ll_fun(mu_temp, p_sub$dis, d_test$y, d_test$w),
                       1, log_mean_exp)
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
    mu_temp <- args$family_kl$linkinv(d_test$x[,chosen]%*%p_sel$b + d_test$offset)
    lppd[[i]] <- apply(args$family_kl$ll_fun(mu_temp, p_sel$dis, d_test$y, d_test$w),
                       1, log_mean_exp)
    mu[[i]] <- rowMeans(mu_temp)

    if(args$verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))
    i <- i + 1
  }

  list(chosen = chosen, kl = kl, mu = mu, lppd = lppd)
}

# function handle for the projection over samples. Gaussian case
# uses analytical solution to do the projection over samples.
.get_proj_handle <- function(family) {

  # Use analytical solution for gaussian as it is a lot faster
  if(family == 'gaussian') {
    function(v_ind, chosen, p_full, d_train, b0, args) {
      v_inds <- c(chosen, v_ind)
      # for forward selection, these measures are precalculated
      covx <- if(is.null(d_train$covX)) crossprod(d_train$x[,v_inds, drop = F]) else d_train$covX[v_inds, v_inds, drop = F]
      x_mu <- if(is.null(p_full$x_mu)) crossprod(d_train$x[, v_inds, drop = F], p_full$mu) else p_full$x_mu[v_inds, , drop = F]
      dis2 <- if(is.null(p_full$dis2)) p_full$dis^2 else p_full$dis2

      # check if covariance matrix is invertible if it seems possible that it might not be
      if(args$rank_x - 2  <= length(v_inds)) {
        if(rankMatrix(covx) < length(v_inds)) return(list(b = NA, dis = NA, kl = Inf))
      }

      # Solution for the gaussian case
      p_sub <- list(b = solve(covx, x_mu))
      p_sub$dis <- sqrt(dis2 + colMeans((p_full$mu - d_train$x[, v_inds, drop = F]%*%p_sub$b)^2))
      p_sub$kl <- weighted.mean(log(p_sub$dis) - log(p_full$dis), p_full$cluster_w)
      p_sub
    }

  } else {
    function(v_ind, chosen, p_full, d_train, b0, args) {
      v_inds <- c(chosen, v_ind)

      # check if covariance matrix is invertible if it seems possible that it might not be
      # preferably this could be removed if NR could be guaranteed not to fail.
      if(args$rank_x - 2  <= length(v_inds)) {
        if(rankMatrix(d_train$x[, v_inds, drop =F]) < length(v_inds))
          return(list(b = NA, dis = NA, kl = Inf))
      }

      # perform the projection over samples
      res <- sapply(1:ncol(p_full$mu), function(s_ind) {
        NR(list(mu = p_full$mu[, s_ind, drop = F], dis = p_full$dis[s_ind]),
           list(x = d_train$x[, v_inds, drop = F], w = d_train$w, offset = d_train$offset),
           b0[v_inds,], args$family_kl)
      })

      # weight the results by sample weights (that are all 1 unless p_clust is used)
      p_sub <- list(kl = weighted.mean(unlist(res['kl',]), p_full$cluster_w),
                b = do.call(cbind, res['b',]))
      if('dis' %in% rownames(res)) p_sub$dis <- unlist(res['dis',])
      p_sub
    }
  }
}
