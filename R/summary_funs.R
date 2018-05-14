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

.get_full_summaries <- function(p_full, d_test, coef_full, family_kl, intercept) {
  mu <- family_kl$mu_fun(d_test$x, coef_full$alpha, coef_full$beta, d_test$offset)
  .weighted_summary_means(d_test, family_kl, p_full$weights, mu, p_full$dis)
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


.bbweights <- function(N,B) {
  # generate Bayesian bootstrap weights, N = original sample size,
  # B = number of bootstrap samples
  bbw <- matrix(rgamma(N*B, 1), ncol = N)
  bbw <- bbw/rowSums(bbw)
  return(bbw)
}


.tabulate_stats <- function(varsel, alpha = 0.05) {
  #
  # return a table of summary statistics with columns:
  #
  #  data: type of data ('sel','train','loo','kfold')
  #  size: number of features in the submodel (Inf indicates the full/reference model)
  #  delta: whether the value indicates the difference to the full model (=TRUE) or the actual value (=FALSE)
  #  statistic: name of statistic ('kl', 'mlpd', 'mse', 'r2', ...)
  #  value: (mean) value of the statistic
  #  lq: lower credible bound for the statistic
  #  uq: upper credible bound for the statistic
  #  se: standard error for the statistic
  
  n <- length(varsel$d_test$y)
  
  # compute statistics for the full model
  summ_ref <- varsel$summaries$full
  stats_ref_pw <- .pointwise_stats(summ_ref$mu, summ_ref$lppd, varsel$d_test, varsel$family_kl)
  stat_names <- names(stats_ref_pw)
  nstats <- length(stat_names)
  m_ref <- colMeans(stats_ref_pw, na.rm = T) # means
  se_ref <- sqrt( apply(stats_ref_pw, 2, 'var') / n ) # standard errors
  row1 <- data.frame(data = varsel$d_test$type, size=Inf, delta=F, statistic=stat_names, value=m_ref, 
                     lq=qnorm(alpha/2, mean=m_ref, sd=se_ref), uq=qnorm(1-alpha/2, mean=m_ref, sd=se_ref),
                     se=se_ref, row.names=1:nstats)
  row2 <- data.frame(data = varsel$d_test$type, size=Inf, delta=T, statistic=stat_names, value=rep(0,nstats), 
                     lq=rep(0,nstats), uq=rep(0,nstats), se=rep(0,nstats), row.names=1:nstats)
  rows_stat <- rbind(row1,row2)
  
  
  # compute statistics for the submodels by looping over the model sizes
  
  # rows_stat <- data.frame()
  for (k in seq_along(varsel$summaries$sub)) {
    
    # summaries for the submodel
    summ_k <- varsel$summaries$sub[[k]]
    stats_pw <- .pointwise_stats(summ_k$mu, summ_k$lppd, varsel$d_test, varsel$family_kl)
    n_notna <- colSums(!is.na(stats_pw)) # how many of the pointwise stats are non-NA
    stat_names <- names(stats_pw)
    
    # pointwise weights
    if (!is.null(summ_k$w))
      w <- summ_k$w/sum(summ_k$w) # pointwise weights
    else
      w <- rep(1/n_notna, n)
    
    # relative to the reference model
    if (all(is.na(stats_ref_pw))) {
      # if the reference model stats are all NA, set the differences to NA as well
      m_diff <- rep(NA, ncol(stats_pw))
      se_diff <- rep(NA, ncol(stats_pw))
    } else {
      pw_diff <- stats_pw - stats_ref_pw
      m_diff <- colSums(w * pw_diff, na.rm = T) # means
      se_diff <- sqrt(colSums(w * pw_diff^2, na.rm = T) - m_diff^2) / sqrt(n_notna) # standard errors
    }
    row1 <- data.frame(data = varsel$d_test$type, size=k-1, delta=T, statistic=stat_names, value=m_diff, 
                       lq=qnorm(alpha/2, mean=m_diff, sd=se_diff), uq=qnorm(1-alpha/2, mean=m_diff, sd=se_diff),
                       se=se_diff, row.names=1:nstats)
    
    # actual value
    if (all(!is.na(stats_ref_pw)) && all(n_notna < n)) {
      # case where the statistics for the reference model have been computed for all the data points,
      # but for the submodels using part of the data only, so compute the results for the submodels
      # as the "difference to the reference model + the result for the reference model" 
      m <- m_diff + m_ref
      m_ref0 <- colSums(w * stats_ref_pw, na.rm=T) # mean for reference model within those points non-NA for submodel
      covterm <- (colSums(w * stats_pw*stats_ref_pw, na.rm = T) - m*m_ref0) / n_notna # covariance uncertainty between submodel and reference model uncertainties
      se <- sqrt(se_diff^2 - se_ref^2 + 2*covterm) # Var(A) = Var(A-B) - Var(B) + 2*Cov(A,B)
    } else {
      m <- colSums(w * stats_pw, na.rm = T) # means
      se <- sqrt(colSums(w * stats_pw^2, na.rm = T) - m^2) / sqrt(n_notna) # standard errors
    }
    row2 <- data.frame(data = varsel$d_test$type, size=k-1, delta=F, statistic=stat_names, value=m, 
                       lq=qnorm(alpha/2, mean=m, sd=se), uq=qnorm(1-alpha/2, mean=m, sd=se), 
                       se=se, row.names=1:nstats)
    
    rows_stat <- rbind(rows_stat, row1, row2)
  }
  
  # kl-values for the submodels and the reference model (for which kl=0 and indicated by size=Inf)
  rows_kl <- data.frame(data = 'sel', size = c(seq_along(varsel$kl)-1, Inf), delta = F,
                      statistic = 'kl',  value = c(varsel$kl, 0), lq = NA, uq = NA, se = NA)
  
  return(rbind(rows_kl, rows_stat))
}


.bootstrap_stats <- function(varsel, n_boot = 1000, alpha = 0.05) {
  #
  # Note: this function is deprecated and not used, see tabulate_stats instead.
  #
  # return a table of summary statistics with columns:
  #
  #  data: type of data ('sel','train','loo','kfold')
  #  size: number of features in the submodel
  #  delta: whether the value indicates the difference to the full model (=TRUE) or the actual value (=FALSE)
  #  statistic: name of statistic ('kl', 'mlpd', 'mse', 'r2', ...)
  #  value: (mean) value of the statistic
  #  lq: lower credible bound for the statistic
  #  uq: upper credible bound for the statistic
  
  n <- length(varsel$d_test$y)
  equal_weights <- matrix(1/n, 1, n)

  b_weights <- .bbweights(n,n_boot)

  # this function computes the average statistics given the 
  # (bootstrap) weights w for the observations
  hf <- function(stat, w) {
    .calc_stats(stat$mu, stat$lppd, varsel$d_test, varsel$family_kl, w)
  }

  # compute for each number of features k the bootstrapped statistics
  sub <- lapply(varsel$summaries$sub, function(stat_k) {
    list(stats = hf(stat_k, equal_weights),  boot = hf(stat_k, b_weights))
  })

  full <- list(stats = hf(varsel$summaries$full, equal_weights),
               boot = hf(varsel$summaries$full, b_weights))

  # apply over submodel sizes
  quantiles <- mapply(function(sub, size) {
    # apply over different stats
    mapply(function(name, sub_boot, full_boot, sub_statistic, full_statistic) {
      qs <- quantile(sub_boot, c(alpha/2, 1-alpha/2), na.rm=T)
      qs_delta <- quantile(sub_boot - full_boot, c(alpha/2, 1-alpha/2), na.rm=T)

      data.frame(size = rep(size, 2), delta = c(F, T), statistic = rep(name, 2),
                 value = c(sub_statistic, sub_statistic - full_statistic),
                 lq = c(qs[1], qs_delta[1]), uq = c(qs[2], qs_delta[2]))
    }, names(sub$boot), sub$boot, full$boot, sub$stats, full$stats, SIMPLIFY = F)
  }, sub, seq_along(varsel$summaries$sub) - 1, SIMPLIFY = F)

  quantiles_arr <- do.call(rbind, c(unlist(quantiles, recursive = F),
                                    make.row.names = F))

  stat_arr <- cbind(data = varsel$d_test$type, quantiles_arr)
  kl_arr <- data.frame(data = 'sel', size = seq_along(varsel$kl)-1, delta = F,
                       statistic = 'kl',  value = varsel$kl, lq = NA, uq = NA)

  rbind(kl_arr, stat_arr)
}






.calc_stats <- function(mu, lppd, d_test, family, sample_weights) {
  #
  # Note: this function is not used.
  #
  # calculate the average of the statistics based on pointwise mu and lppd,
  # assuming the observations are given weights sample_weights (which can be used
  # for bootstrapping the statistics)
  arr <- list(mlpd = lppd, elpd = lppd*length(mu), mse = (d_test$y-mu)^2)

  if(family$family == 'binomial' && all(d_test$weights %in% c(0,1))) {
    arr$pctcorr <- round(mu) == d_test$y
  }

  avg_ <- function(x) c(sample_weights%*%x)
  stats <- lapply(arr, avg_)

  if(family$family == 'gaussian') {
    stats$r2 <- 1 - stats$mse/avg_((d_test$y-mean(d_test$y))^2)
  }

  stats
}


.pointwise_stats <- function(mu, lppd, d_test, family, sample_weights=NULL) {
  #
  # calculate the pointwise statistics based on pointwise mu and lppd
  #
  stats <- list()
  stats$mlpd <- lppd
  stats$elpd <- lppd*length(mu)
  stats$mse <- (d_test$y-mu)^2
  stats$rmse <- sqrt(stats$mse)
  
  if(family$family == 'gaussian') {
    stats$r2 <- 1 - stats$mse/mean((d_test$y-mean(d_test$y))^2)
  }
  if(family$family == 'binomial' && all(d_test$weights %in% c(0,1))) {
    stats$pctcorr <- round(mu) == d_test$y
    stats$acc <- stats$pctcorr
  }
  
  as.data.frame(stats)
}

