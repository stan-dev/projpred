

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
# mu and dis, and the test data.
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





.tabulate_stats <- function(varsel, stats, alpha = 0.05, nfeat_baseline=NULL) {
  
  #
  # Calculates the desired statistics, their standard errors and credible bounds with given
  # credible level alpha based on the variable selection information. If nfeat_baseline
  # is given, then compute the statistics relative to the baseline model with that size
  # (nfeat_baseline = Inf means reference model).
  
  stat_tab <- data.frame()
  
  summ_ref <- varsel$summaries$ref
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
      mu.bs[is.na(mu)] <- NA # make sure the relative rmse is computed using only those points for which 
      mu[is.na(mu.bs)] <- NA # both mu and mu.bs are non-NA
      value <- sqrt(mean(weights*(mu-y)^2, na.rm=T)) -  sqrt(mean(weights*(mu.bs-y)^2, na.rm=T))
      value.bootstrap1 <- bootstrap((mu-y)^2, function(resid2) sqrt(mean(weights*resid2, na.rm=T)), b=B, seed=seed)
      value.bootstrap2 <- bootstrap((mu.bs-y)^2, function(resid2) sqrt(mean(weights*resid2, na.rm=T)), b=B, seed=seed)
      value.se <- sd(value.bootstrap1-value.bootstrap2)
    } else {
      value <- sqrt(mean(weights*(mu-y)^2, na.rm=T))
      value.bootstrap <- bootstrap((mu-y)^2, function(resid2) sqrt(mean(weights*resid2, na.rm=T)), b=B, seed=seed)
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
    
  } else if (stat == 'auc') {

    if (family$family != 'binomial')
      stop('Statistic ', stat, ' available only for binomial family.')
    y <- d_test$y
    auc.data <- cbind(y, mu, weights)[!is.na(mu), ]
    if (!is.null(mu.bs)) {
      auc.data.bs <- cbind(y, mu.bs, weights)[!is.na(mu), ]
      value <- auc(auc.data) - auc(auc.data.bs)
      value.bootstrap1 <- bootstrap(auc.data, auc, b=B, seed=seed)
      value.bootstrap2 <- bootstrap(auc.data.bs, auc, b=B, seed=seed)
      value.se <- sd(value.bootstrap1 - value.bootstrap2)
    } else {
      value <- auc(auc.data)
      value.bootstrap <- bootstrap(auc.data, auc, b=B, seed=seed)
      value.se <- sd(value.bootstrap)
    }

  } else {
    stop('Unknown statistic: ', stat)
  }
  
  lq <- qnorm(alpha/2, mean=value, sd=value.se)
  uq <- qnorm(1-alpha/2, mean=value, sd=value.se)
  
  return(list(value=value, se=value.se, lq=lq, uq=uq))
}



.is_util <- function(stat) {
  # a simple function to determine whether a given statistic (string) is
  # a utility (we want to maximize) or loss (we want to minimize)
  # by the time we get here, stat should have already been validated

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
