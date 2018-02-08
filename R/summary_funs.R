# Functions for calculating mu and lppd (ge)
# - .get_sub/full_summaries + .weighted_summary_means
# Functions for calculating mse, mlpd, etc. with (and without) bootstrapping
# - .bootstrap/calc_stats,

.get_sub_summaries <- function(submodels, d_test, family_kl) {

  res <- lapply(submodels, function(submodels) {
  	vind <- submodels$vind
    if(NROW(submodels$beta) == 0) {
      xt <- matrix(0, nrow = length(d_test$weights), ncol = 0)
    } else if(!is.matrix(d_test$x)) {
      xt <- matrix(d_test$x[vind], nrow = 1)
    } else {
      xt <- d_test$x[, vind, drop = F]
    }

    mu <- family_kl$mu_fun(xt, submodels$alpha, submodels$beta, d_test$offset)
    .weighted_summary_means(d_test, family_kl, submodels$weights, mu, submodels$dis)
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
  #  size: number of features in the submodel
  #  delta: whether the value indicates the difference to the full model (=TRUE) or the actual value (=FALSE)
  #  statistic: name of statistic ('kl', 'mlpd', 'mse', 'r2', ...)
  #  value: (mean) value of the statistic
  #  lq: lower credible bound for the statistic
  #  uq: upper credible bound for the statistic
  
  n <- length(varsel$d_test$y)
  
  
  # compute statistics for the full model
  summ_ref <- varsel$summaries$full
  stats_ref_pw <- .pointwise_stats(summ_ref$mu, summ_ref$lppd, varsel$d_test, varsel$family_kl)
  stats_names <- names(stats_ref_pw)
  
  # compute statistics for the submodels
  
  # loop over model sizes
  rows_stat <- data.frame()
  for (k in seq_along(varsel$summaries$sub)) {
    
    summ_k <- varsel$summaries$sub[[k]]
    stats_pw <- .pointwise_stats(summ_k$mu, summ_k$lppd, varsel$d_test, varsel$family_kl)
    stat_names <- names(stats_pw)
    
    # actual value
    m <- colMeans(stats_pw) # means
    se <- sqrt(apply(stats_pw,2,'var') / n) # standard errors DIVIDE BY NUMBER OF NON-NAS!!
    row1 <- data.frame(data = varsel$d_test$type, size=k-1, delta=F, statistic=stat_names, value=m, 
                       lq=qnorm(alpha/2, mean=m, sd=se), uq=qnorm(1-alpha/2, mean=m, sd=se), 
                       row.names=1:length(m))
    
    # relative to the reference model
    m <- colMeans(stats_pw-stats_ref_pw) # means
    se <- sqrt(apply(stats_pw-stats_ref_pw,2,'var') / n) # standard errors DIVIDE BY NUMBER OF NON-NAS!!
    row2 <- data.frame(data = varsel$d_test$type, size=k-1, delta=T, statistic=stat_names, value=m, 
                       lq=qnorm(alpha/2, mean=m, sd=se), uq=qnorm(1-alpha/2, mean=m, sd=se),
                       row.names=1:length(m))
    
    rows_stat <- rbind(rows_stat, row1, row2)
  }
  
  # kl-value
  rows_kl <- data.frame(data = 'sel', size = seq_along(varsel$kl)-1, delta = F,
                      statistic = 'kl',  value = varsel$kl, lq = NA, uq = NA)
  
  return(rbind(rows_kl, rows_stat))
}


.bootstrap_stats <- function(varsel, n_boot = 1000, alpha = 0.05) {
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
  # calculate the average of the statistics based on pointwise mu and lppd,
  # assuming the observations are given weights sample_weights (which can be used
  # for bootstrapping the statistics)
  arr <- list(mlpd = lppd, mse = (d_test$y-mu)^2)

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
  # calculate the pointwise statistics based on pointwise mu and lppd
  stats <- list()
  stats$mlpd <- lppd
  stats$mse <- (d_test$y-mu)^2
  
  if(family$family == 'gaussian') {
    stats$r2 <- 1 - stats$mse/mean((d_test$y-mean(d_test$y))^2)
  }
  if(family$family == 'binomial' && all(d_test$weights %in% c(0,1))) {
    stats$pctcorr <- round(mu) == d_test$y
  }
  
  # avg_ <- function(x) c(sample_weights%*%x)
  # stats <- lapply(arr, avg_)
  
  # if(family$family == 'gaussian') {
  #   stats$r2 <- 1 - stats$mse/avg_((d_test$y-mean(d_test$y))^2)
  # }
  as.data.frame(stats)
  # stats
}

