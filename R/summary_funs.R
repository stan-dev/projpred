# Functions for calculating mu and lppd (ge)
# - .get_sub/full_summaries + .weighted_summary_means
# Functions for calculating mse, mlpd, etc. with (and without) bootstrapping
# - .bootstrap/calc_stats,

.get_sub_summaries <- function(submodels, d_test, family_kl) {

  res <- lapply(submodels, function(model) {
  	vind <- model$vind
    if(NROW(model$beta) == 0) {
      xt <- matrix(0, nrow = length(d_test$weights), ncol = 0)
    } else if(!is.matrix(d_test$x)) {
      xt <- matrix(d_test$x[vind], nrow = 1)
    } else {
      xt <- d_test$x[, vind, drop = F]
    }
    mu <- family_kl$mu_fun(xt, model$alpha, model$beta, d_test$offset)
    .weighted_summary_means(d_test, family_kl, model$weights, mu, model$dis)
  })
}


# Calculates weighted means of mu and lppd given samples of
# mu and dis, the full model and the data.
.weighted_summary_means <- function(d_test, family_kl, wsample, mu, dis) {

  loglik <- family_kl$ll_fun(mu, dis, matrix(d_test$y,nrow=NROW(mu)), d_test$weights)
  if (length(loglik) == 1) {
    # one observation, one sample
    list(mu = mu, lppd = loglik)
  } else if (is.null(dim(loglik))){
    # loglik is a vector, but not sure if it means one observation with many samples, or vice versa?
    stop('Internal error encountered: loglik is a vector, but should be a scalar or matrix')
  } else {
    # mu is a matrix, so apply weighted sum over the samples
    list(mu = c(mu %*% wsample),
         lppd = apply(loglik, 1, log_weighted_mean_exp, wsample))
  }
  
}



# .tabulate_stats_old <- function(varsel, alpha = 0.05, nfeat_baseline=Inf) {
#   #
#   # return a table of summary statistics with columns:
#   #
#   #  data: type of data ('sel','train','loo','kfold')
#   #  size: number of features in the submodel (Inf indicates the full/reference model)
#   #  delta: whether the value indicates the difference to the full model (=TRUE) or the actual value (=FALSE)
#   #  statistic: name of statistic ('kl', 'mlpd', 'mse', 'r2', ...)
#   #  value: (mean) value of the statistic
#   #  lq: lower credible bound for the statistic
#   #  uq: upper credible bound for the statistic
#   #  se: standard error for the statistic
#   #
#   # nfeat_baseline indicates the size of the model with respect to which the differences
#   # are computed (Inf means reference model)
#   
#   n <- length(varsel$d_test$y)
#   
#   # compute the baseline statistics
#   if (nfeat_baseline == Inf) {
#     summ_baseline <- varsel$summaries$full
#   } else {
#     summ_baseline <- varsel$summaries$sub[[nfeat_baseline+1]]
#   }
#   stats_bs_pw <- .pointwise_stats(summ_baseline$mu, summ_baseline$lppd, varsel$d_test, varsel$family_kl)
#   n_notna <- colSums(!is.na(stats_bs_pw)) # how many of the pointwise difference stats are non-NA
#   stat_names <- names(stats_bs_pw)
#   nstats <- length(stat_names)
#   
#   # pointwise weights for the baseline statistics
#   if (!is.null(summ_baseline$w))
#     w <- summ_baseline$w/sum(summ_baseline$w) # pointwise weights
#   else
#     w <- sapply(n_notna, function(nna) rep(1/nna, n))
#   if (NCOL(w)==1)
#     # repeat w for convenience of the upcoming calculations
#     w <- matrix(rep(w, nstats), ncol=nstats)
#   
#   # mean and se of baseline stats, these might be needed later on
#   m_bs <- colSums(w * stats_bs_pw, na.rm = T) 
#   se_bs <- sapply(1:nstats, function(i) weighted.sd(stats_bs_pw[,i], w[,i], na.rm=T)) / sqrt(n_notna)
#   
#   # compute statistics for the reference model
#   summ_ref <- varsel$summaries$full
#   stats_ref_pw <- .pointwise_stats(summ_ref$mu, summ_ref$lppd, varsel$d_test, varsel$family_kl)
#   pw_diff <- stats_ref_pw - stats_bs_pw
#   m_ref <- colMeans(stats_ref_pw, na.rm = T) # means
#   se_ref <- sqrt( apply(stats_ref_pw, 2, 'var') / n ) # standard errors
#   n_notna <- colSums(!is.na(pw_diff)) # how many of the pointwise difference stats are non-NA
#   m_diff <- colSums(w * pw_diff, na.rm = T) # relative mean
#   se_diff <- sapply(1:nstats, function(i) weighted.sd(pw_diff[,i], w[,i], na.rm=T)) / sqrt(n_notna) # relative standard errors
#   row1 <- data.frame(data = varsel$d_test$type, size=Inf, delta=F, statistic=stat_names, value=m_ref, 
#                      lq=qnorm(alpha/2, mean=m_ref, sd=se_ref), uq=qnorm(1-alpha/2, mean=m_ref, sd=se_ref),
#                      se=se_ref, row.names=1:nstats)
#   row2 <- data.frame(data = varsel$d_test$type, size=Inf, delta=T, statistic=stat_names, value=m_diff, 
#                      lq=qnorm(alpha/2, mean=m_diff, sd=se_diff), uq=qnorm(1-alpha/2, mean=m_diff, sd=se_diff),
#                      se=se_diff, row.names=1:nstats)
#   rows_stat <- rbind(row1,row2)
#   
#   
#   # compute statistics for the submodels by looping over the model sizes
#   
#   for (k in seq_along(varsel$summaries$sub)) {
#   	
#     # summaries for the submodel
#     summ_k <- varsel$summaries$sub[[k]]
#     stats_pw <- .pointwise_stats(summ_k$mu, summ_k$lppd, varsel$d_test, varsel$family_kl)
#     n_notna <- colSums(!is.na(stats_pw)) # how many of the pointwise stats are non-NA
#     stat_names <- names(stats_pw)
#     
#     # pointwise weights
#     if (!is.null(summ_k$w))
#       w <- summ_k$w/sum(summ_k$w) # pointwise weights
#     else
#     	w <- sapply(n_notna, function(nna) rep(1/nna, n))
#     if (NCOL(w)==1)
#     	# repeat w for convenience of the upcoming calculations
#     	w <- matrix(rep(w, nstats), ncol=nstats)
#     
#     # relative to the baseline model
#     if (all(is.na(stats_bs_pw))) {
#       # if the baseline model stats are all NA, set the differences to NA as well
#       m_diff <- rep(NA, ncol(stats_pw))
#       se_diff <- rep(NA, ncol(stats_pw))
#     } else {
#       pw_diff <- stats_pw - stats_bs_pw
#       m_diff <- colSums(w * pw_diff, na.rm = T) # means
#       se_diff <- sapply(1:nstats, function(i) weighted.sd(pw_diff[,i], w[,i], na.rm=T)) / sqrt(n_notna) # standard errors
#       # se_diff <- sqrt(colSums(w * pw_diff^2, na.rm = T) - m_diff^2) / sqrt(n_notna) # alternative way, but numerically unstable
#     }
#     row1 <- data.frame(data = varsel$d_test$type, size=k-1, delta=T, statistic=stat_names, value=m_diff, 
#                        lq=qnorm(alpha/2, mean=m_diff, sd=se_diff), uq=qnorm(1-alpha/2, mean=m_diff, sd=se_diff),
#                        se=se_diff, row.names=1:nstats)
#     
#     # actual value
#     if (all(!is.na(stats_bs_pw)) && all(n_notna < n)) {
#       # case where the statistics for the reference model have been computed for all the data points,
#       # but for the submodels using part of the data only, so compute the results for the submodels
#       # as the "difference to the reference model + the result for the reference model" 
#       m <- m_diff + m_bs
#       m_bs0 <- colSums(w * stats_bs_pw, na.rm=T)
#       # m0 <- m_diff + m_bs0 
#       se <- sqrt(se_diff^2 + se_bs^2)
#       # covterm <- (colSums(w * stats_pw*stats_bs_pw, na.rm = T) - m0*m_bs0) / n_notna # covariance uncertainty between submodel and reference model uncertainties
#       # se <- sqrt(se_diff^2 - se_bs^2 + 2*covterm) # Var(A) = Var(A-B) - Var(B) + 2*Cov(A,B), A = submodel statistic, B = refmodel statistic
#     } else {
#       m <- colSums(w * stats_pw, na.rm = T) # means
#       se <- sapply(1:nstats, function(i) weighted.sd(stats_pw[,i], w[,i], na.rm=T)) / sqrt(n_notna) # standard errors
#       # se <- sqrt(colSums(w * stats_pw^2, na.rm = T) - m^2) / sqrt(n_notna) # alternative way, but numerically unstable..
#     }
#     row2 <- data.frame(data = varsel$d_test$type, size=k-1, delta=F, statistic=stat_names, value=m, 
#                        lq=qnorm(alpha/2, mean=m, sd=se), uq=qnorm(1-alpha/2, mean=m, sd=se), 
#                        se=se, row.names=1:nstats)
#     
#     rows_stat <- rbind(rows_stat, row1, row2)
#   }
#   
#   # kl-values for the submodels and the reference model (for which kl=0 and indicated by size=Inf)
#   rows_kl <- data.frame(data = 'sel', size = c(seq_along(varsel$kl)-1, Inf), delta = F,
#                       statistic = 'kl',  value = c(varsel$kl, 0), lq = NA, uq = NA, se = NA)
#   
#   return(rbind(rows_kl, rows_stat))
# }






.tabulate_stats <- function(varsel, stats, alpha = 0.05, nfeat_baseline=NULL) {
  
  stat_tab <- data.frame()
  
  summ_ref <- varsel$summaries$full
  summ_sub <- varsel$summaries$sub
  
  # fetch the mu and lppd for the baseline model
  if (is.null(nfeat_baseline)) {
    # no baseline model, i.e, compute the statistics on the actual (non-relative) scale
    mu.bs <- NULL
    lppd.bs <- NULL
    delta <- FALSE
  } else {
    if (nfeat_baseline == Inf)
      summ.bs <- summ_ref
    else
      summ.bs <- summ_sub[[nfeat_baseline+1]]
    mu.bs <- summ.bs$mu
    lppd.bs <- summ.bs$lppd
    delta <- TRUE
  }
  
  for (s in seq_along(stats)) {
    
    stat <- stats[s]
    
    # reference model statistics
    summ <- summ_ref
    res <- get_stat(summ$mu, summ$lppd, varsel$d_test, varsel$family_kl, stat, 
                    mu.bs=mu.bs, lppd.bs=lppd.bs, weights=summ$w, alpha=alpha)
    row <- data.frame(data = varsel$d_test$type, size=Inf, delta=delta, statistic=stat, 
                      value=res$value, lq=res$lq, uq=res$uq, se=res$se)
    stat_tab <- rbind(stat_tab, row)
    
    # submodel statistics
    for (k in seq_along(varsel$summaries$sub)) {
      
      summ <- summ_sub[[k]]
      if (delta == F && sum(!is.na(summ_ref$mu)) > sum(!is.na(summ$mu))) {
        # special case (subsampling loo): reference model summaries computed for more points
        # than for the submodel, so utilize the reference model results to get more accurate 
        # statistic fot the submodel on the actual scale 
        res_ref <- get_stat(summ_ref$mu, summ_ref$lppd, varsel$d_test, varsel$family_kl, stat, 
                            mu.bs=NULL, lppd.bs=NULL, weights=summ_ref$w, alpha=alpha)
        res_diff <- get_stat(summ$mu, summ$lppd, varsel$d_test, varsel$family_kl, stat, 
                             mu.bs=summ_ref$mu, lppd.bs=summ_ref$lppd, weights=summ$w, alpha=alpha)
        val <- res_ref$value+res_diff$value
        val.se <- sqrt(res_ref$se^2+res_diff$se^2)
        lq <- qnorm(alpha/2, mean=val, sd=val.se)
        uq <- qnorm(1-alpha/2, mean=val, sd=val.se)
        row <- data.frame(data = varsel$d_test$type, size=k-1, delta=delta, statistic=stat, 
                          value=val, lq=lq, uq=uq, se=val.se)
      } else {
        # normal case
        res <- get_stat(summ$mu, summ$lppd, varsel$d_test, varsel$family_kl, stat, 
                        mu.bs=mu.bs, lppd.bs=lppd.bs, weights=summ$w, alpha=alpha)
        row <- data.frame(data = varsel$d_test$type, size=k-1, delta=delta, statistic=stat, 
                          value=res$value, lq=res$lq, uq=res$uq, se=res$se)
      }
      stat_tab <- rbind(stat_tab, row)
    }
  }
  
  stat_tab
}




get_stat <- function(mu, lppd, d_test, family, stat, mu.bs=NULL, lppd.bs=NULL, 
                     weights=NULL, alpha=0.1, seed=1208499, B=2000) {
  #
  # Calculates given statistic stat with standard error and confidence bounds.
  # mu.bs and lppd.bs are the pointwise mu and lppd for another model that is
  # used as a baseline for computing the difference in the given statistic,
  # for example the relative elpd. If these arguments are not given (NULL) then
  # the actual (non-relative) value is computed.
  
  n <- length(mu)
  
  if (stat %in% c('mlpd','elpd'))
    n_notna <- sum(!is.na(lppd))
  else
    n_notna <- sum(!is.na(mu))
  
  if (is.null(weights))
    # set default weights if not given
    weights <- rep(1/n_notna, n)
  # ensure the weights sum to n_notna
  weights <- n_notna*weights/sum(weights) 
  
  
  if (stat == 'mlpd') {
    
    if (!is.null(lppd.bs)) {
      value <- mean((lppd-lppd.bs)*weights, na.rm=T)
      value.se <- weighted.sd(lppd-lppd.bs, weights, na.rm=T) / sqrt(n_notna)
    } else {
      value <- mean(lppd*weights, na.rm=T)
      value.se <- weighted.sd(lppd, weights, na.rm=T) / sqrt(n_notna)
    }
    
  } else if (stat == 'elpd') {
    
    if (!is.null(lppd.bs)) {
      value <- sum((lppd-lppd.bs)*weights, na.rm=T)
      value.se <- weighted.sd(lppd-lppd.bs, weights, na.rm=T) / sqrt(n_notna) * n_notna
    } else {
      value <- sum(lppd*weights, na.rm=T)
      value.se <- weighted.sd(lppd, weights, na.rm=T) / sqrt(n_notna) * n_notna
    }
    
  } else if (stat == 'mse') {
    
    y <- d_test$y
    if (!is.null(mu.bs)) {
      value <- mean(weights*((mu-y)^2 - (mu.bs-y)^2), na.rm=T)
      value.se <- weighted.sd((mu-y)^2 - (mu.bs-y)^2, weights, na.rm=T) / sqrt(n_notna)
    } else {
      value <- mean(weights*(mu-y)^2, na.rm=T)
      value.se <- weighted.sd((mu-y)^2, weights, na.rm=T) / sqrt(n_notna)
    }
    
  } else if (stat == 'rmse') {
    
    y <- d_test$y
    if (!is.null(mu.bs)) {
      mu.bs[is.na(mu)] <- NA # make sure the relative rmse is computed using only those points for which mu is not NA
      value <- sqrt(mean(weights*(mu-y)^2, na.rm=T)) -  sqrt(mean(weights*(mu.bs-y)^2, na.rm=T))
      value.bootstrap1 <- bootstrap((mu-y)^2, function(resid2) sqrt(mean(weights*resid2, na.rm=T)), b=B, seed=seed)
      value.bootstrap2 <- bootstrap((mu.bs-y)^2, function(resid2) sqrt(mean(weights*resid2, na.rm=T)), b=B, seed=seed)
      value.se <- sd(value.bootstrap1-value.bootstrap2)
    } else {
      value <- sqrt(mean(weights*(mu-y)^2, na.rm=T))
      value.bootstrap <- bootstrap((mu-y)^2, function(resid2) sqrt(mean(weights*resid2, na.rm=T)), b=B)
      value.se <- sd(value.bootstrap)
    }
    
  } else if (stat == 'acc' || stat == 'pctcorr') {
    
    if (family$family != 'binomial')
      stop('Statistic ', stat, ' available only for binomial family.')
    y <- d_test$y
    if (!is.null(mu.bs)) {
      value <- mean( weights*((round(mu)==y) - (round(mu.bs)==y)), na.rm=T )
      value.se <- weighted.sd((round(mu)==y) - (round(mu.bs)==y), weights, na.rm=T) / sqrt(n_notna)
    } else {
      value <- mean(weights*(round(mu) == y), na.rm=T)
      value.se <- weighted.sd(round(mu) == y, weights, na.rm=T) / sqrt(n_notna)
    }
    
  } else {
    stop('Unknown statistic: ', stat)
  }
  
  lq <- qnorm(alpha/2, mean=value, sd=value.se)
  uq <- qnorm(1-alpha/2, mean=value, sd=value.se)
  
  return(list(value=value, se=value.se, lq=lq, uq=uq))
}



# .pointwise_stats <- function(mu, lppd, d_test, family, sample_weights=NULL) {
#   #
#   # calculate the pointwise statistics based on pointwise mu and lppd
#   #
#   stats <- list()
#   stats$mlpd <- lppd
#   stats$elpd <- lppd*length(mu)
#   stats$mse <- (d_test$y-mu)^2
#   stats$rmse <- sqrt(stats$mse)
#   
#   if(family$family == 'gaussian') {
#     stats$r2 <- 1 - stats$mse/mean((d_test$y-mean(d_test$y))^2)
#   }
#   if(family$family == 'binomial' && all(d_test$weights %in% c(0,1))) {
#     stats$pctcorr <- round(mu) == d_test$y
#     stats$acc <- stats$pctcorr
#   }
#   
#   as.data.frame(stats)
# }

.is_util <- function(stat) {
  # a simple function to determine whether a given statistic (string) is
  # a utility (we want to maximize) or loss (we want to minimize)
  recognized_stats <- c('elpd','mlpd','acc','pctcorr','r2','mse','rmse')
  if (!(stat %in% recognized_stats))
    stop(sprintf('Internal error: non-recognized statistic \'%s\'', stat))
  if (stat %in% c('rmse','mse'))
    return(F)
  else
    return(T)
}

.get_nfeat_baseline <- function(object, baseline, stat) {
  # get model size that is used as a baseline in comparisons.
  # baseline is one of 'best' or 'ref', stat is the statistic according to which
  # the selection is done
  if (baseline == 'best') {
    # find number of features that maximizes the utility (or minimizes the loss)
    tab <- .tabulate_stats(object, stat)
    stats_table <- subset(tab, tab$size != Inf)
    # tab <- .tabulate_stats(object)
    # stats_table <- subset(tab, tab$delta == F & tab$statistic == stat & tab$size != Inf)
    optfun <- ifelse(.is_util(stat), which.max, which.min)
    nfeat_baseline <- stats_table$size[optfun(stats_table$value)]
  } else {
    # use reference model
    nfeat_baseline <- Inf
  }
  return(nfeat_baseline)
}
