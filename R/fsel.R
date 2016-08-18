
fsel <- function(p, d, d_test, p_means = NULL, b0, args) {


  # for gaussian, precompute some values
  if(args$family_kl$family == 'gaussian') {
    d$covX <- crossprod(d$x)
    p$x_mu <- crossprod(d$x, p$mu)
    p$dis2 <- p$dis^2
    if(args$avg) {
      p_means$x_mu <- crossprod(d$x,p_means$mu)
      p_means$dis2 <- p_means$dis^2
    }
  }

  # initialize forward selection
  # proj calculates
  proj <- .get_proj_handle(args$family_kl$family)
  i <- 1
  iq <- ceiling(quantile(1:args$nv, 1:10/10))
  cols <- 1:ncol(d$x)
  chosen <- NULL
  kl <- rep(NA_real_, args$nv)
  mu <- vector('list', args$nv)
  lppd <- vector('list', args$nv)
  p_fsel <- if(args$avg) p_means else p

  # if the full model has intercept, add it to the submodel 'as first variable'
  if(args$intercept) {
    chosen <- 1
    q <- proj(NULL, chosen, p_fsel, d, b0, args)
    kl[i] <- q$kl
    mu_temp <- args$family_kl$linkinv(d_test$x[,chosen]%*%q$b)
    mu[[i]] <- rowMeans(mu_temp)
    lppd[[i]] <- apply(args$family_kl$ll_fun(mu_temp, q$dis, d_test$y, d_test$w), 1, log_mean_exp)
    i <- i + 1
  }

  # start adding variables one at a time
  while(i <= args$nv) {

    notchosen <- setdiff(cols, chosen)

    q <- sapply(notchosen, proj, chosen, p_fsel, d, b0, args)
    kl_temp <- unlist(q['kl',])
    if(min(kl_temp) == Inf) {
      warning(paste0('Numerical problems in the projection for a submodel of size ', i, '. Ending forward selection.'))
      return(list(chosen = chosen[1:(i-1)], kl = kl[1:(i-1)], mu = mu[1:(i-1)], lppd = lppd[1:(i-1)]))
    }
    imin <- which.min(kl_temp)

    q_sel <- if(args$avg) proj(notchosen[imin], chosen, p_fsel, d, b0, args) else q[,imin]
    chosen <- c(chosen, notchosen[imin])

    kl[i] <- q_sel$kl
    mu_temp <- args$family_kl$linkinv(d_test$x[,chosen]%*%q_sel$b + d_test$offset)
    lppd[[i]] <- apply(args$family_kl$ll_fun(mu_temp, q_sel$dis, d_test$y, d_test$w), 1, log_mean_exp)
    mu[[i]] <- rowMeans(mu_temp)

    if(args$verbose && i %in% iq) print(paste0(names(iq)[max(which(i == iq))], " of variables selected."))
    i <- i + 1
  }

  list(chosen = chosen, kl = kl, mu = mu, lppd = lppd)
}

.get_proj_handle <- function(family) {

  # use analytical solution for gaussian
  if(family == 'gaussian') {
    function(v_ind, chosen, p, d, b0, args) {
      v_inds <- c(chosen, v_ind)
      covx <- d$covX[v_inds, v_inds, drop = F]

      if(ncol(d$x) - args$rank_x + length(v_inds) > ncol(covx)) {
        if(rankMatrix(covx) < length(v_inds)) return(list(b = NA, dis = NA, kl = Inf))
      }

      q <- list(b = solve(covx, p$x_mu[v_inds, , drop = F]))
      q$dis <- sqrt(p$dis2 + colMeans((p$mu - d$x[, v_inds, drop = F]%*%q$b)^2))
      q$kl <- sum(log(q$dis) - log(p$dis))/args$ns
      q
    }

  } else {
    proj <- function(v_ind, chosen, p, d, b0, args) {
      v_inds <- c(chosen, v_ind)

      if(ncol(d$x) - args$rank_x + length(v_inds) > ncol(d$x[, v_inds, drop =F])) {
        if(rankMatrix(d$x[, v_inds, drop =F]) < length(v_inds)) return(list(b = NA, dis = NA, kl = Inf))
      }

      res <- sapply(1:args$ns, function(s_ind) {
        NR(list(mu = p$mu[, s_ind, drop = F], dis = p$dis[s_ind]),
           list(x = d$x[, v_inds, drop = F], w = d$w, offset = d$offset), b0[v_inds,], args$family_kl)
      })

      q <- list(kl = sum(unlist(res['kl',]))/args$ns, b = do.call(cbind, res['b',]))
      if('dis' %in% rownames(res)) q$dis <- unlist(res['dis',])
      q
    }
  }
}
