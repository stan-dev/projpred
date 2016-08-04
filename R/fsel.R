
fsel <- function(p, d, d_test, p_means = NULL, b0, args) {

  # use analytical solution for gaussian
  # and precompute some values
  if(args$family$family == 'gaussian') {
    d$covX <- crossprod(d$x)
    p$x_mu <- crossprod(d$x,p$mu)
    p$dis2 <- p$dis^2
    if(args$avg) {
      p_means$x_mu <- crossprod(d$x,p_means$mu)
      p_means$dis2 <- p_means$dis^2
    }
    kl_ind <- function(ind, chosen, p, d, b0, family) {
      varinds <- c(chosen, ind)
      q <- list(b = solve(d$covX[varinds, varinds, drop = F], p$x_mu[varinds, , drop = F]))
      q$dis <- sqrt(p$dis2 + colMeans((p$mu - d$x[, varinds, drop = F]%*%q$b)^2))
      q$kl <- sum(log(q$dis) - log(p$dis))/length(q$dis)
      q
    }
  } else {
    # function to calculate kl over samples
    kl_ind <- function(ind, chosen, p, d, b0, family) {
      varinds <- c(ind, chosen)
      s <- ncol(p$mu)
      res <- sapply(1:s, function(sind) {
        NR(list(mu = p$mu[, sind, drop = F], dis = p$dis[sind]),
           list(x = d$x[, varinds, drop = F], w = d$w), b0[varinds,], family)
      })
      q <- list(kl = sum(unlist(res['kl',]))/s, b = do.call(cbind, res['b',]))
      if('dis' %in% rownames(res)) q$dis <- unlist(res['dis',])
      q
    }
  }

  kl <- rep(NA_real_, args$d)
  test_stats <- matrix(NA_real_, nrow = args$d, ncol = 4, dimnames = list(NULL, c('mse', 'mlpd', 'r2', 'pctcorr')))

  i <- 1

  if(args$intercept) {
    chosen <- 1
    # calculate KL divergence with just intercept
    q <- kl_ind(NULL, chosen, p, d, b0, args$family)
    kl[i] <- q$kl
    if(is.list(d_test)) test_stats[i,] <- test_stat(q, chosen, d_test, args$family)
    #if(args$verbose) print(paste0(i, " of ", args$d, " variables selected."))
    i <- i + 1
  } else {
    chosen <- NULL
  }

  cols <- 1:ncol(d$x)
  notchosen <- setdiff(cols, chosen)
  if(args$verbose) iq <- ceiling(quantile(1:args$d, 1:10/10))

  # start adding variables one at a time
  while(i <= args$d) {

    if(args$avg) {
      q <- sapply(notchosen, kl_ind, chosen, p_means, d, b0, args$family)
    } else {
      q <- sapply(notchosen, kl_ind, chosen, p, d, b0, args$family)
    }

    imin <- which.min(q['kl',])
    chosen <- c(chosen, notchosen[imin])

    if(args$avg) {
      q <- kl_ind(NULL, chosen, p, d, b0, args$family)
      kl[i] <- q$kl
      if(is.list(d_test)) test_stats[i,] <- test_stat(q, chosen, d_test, args$family)
    } else {
      kl[i] <- q[['kl',imin]]
      if(is.list(d_test)) test_stats[i,] <- test_stat(q[, imin], chosen, d_test, args$family)
    }

    if(args$verbose && i %in% iq)
      print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))
    notchosen <- setdiff(cols, chosen)
    i <- i + 1
  }

  res <- list(chosen = chosen, kl = kl)
  if(is.list(d_test)) {
    res$mse <- test_stats[,'mse']
    res$mlpd <- test_stats[,'mlpd']
    if(args$family$family == 'gaussian') res$r2 <- test_stats[,'r2']
    if(args$family$family == 'binomial') res$pctcorr <- test_stats[,'pctcorr']
  }

  res
}
