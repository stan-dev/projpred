#' Functions for calculating mu and lppd (ge)
#' - .get_sub/full_summaries + .weighted_summary_means
#' Functions for calculating mse, mlpd, etc. with (and without) bootstrapping
#' - .bootstrap/calc_stats,

.get_sub_summaries <- function(chosen, d_test, p_sub, family_kl) {

  res <- lapply(p_sub, function(p_sub) {
  	ind <- 1:NROW(p_sub$beta)
    if(NROW(p_sub$beta) == 0) {
      xt <- matrix(0, nrow = length(d_test$weights), ncol = 0)
    } else if(!is.matrix(d_test$x)) {
      xt <- matrix(d_test$x[ind], nrow = 1)
    } else {
      xt <- d_test$x[, chosen[ind], drop = F]
    }

    mu <- family_kl$mu_fun(xt, p_sub$alpha, p_sub$beta, d_test$offset)
    .weighted_summary_means(d_test, family_kl, p_sub$weights, mu, p_sub$dis)
  })
}

.get_full_summaries <- function(p_full, d_test, coef_full, family_kl, intercept) {
  mu <- family_kl$mu_fun(d_test$x, coef_full$alpha, coef_full$beta, d_test$offset)
  .weighted_summary_means(d_test, family_kl, p_full$weights, mu, p_full$dis)
}

# Calculates weighted means of mu, dis, kl and lppd given samples of
# mu and dis, the full model and the data.
.weighted_summary_means <- function(d_test, family_kl, wsample, mu, dis) {

  loglik <- family_kl$ll_fun(mu, dis, matrix(d_test$y,nrow=NROW(mu)))
  if (length(loglik) == 1) {
      # one observation, one sample
      list(mu = mu, lppd = loglik)
  } else if (is.null(dim(loglik))){
      # loglik is a vector, but not sure if it means one observation with many samples, or vice versa?
      stop('loglik is a vector, but should be a scalar or matrix')
  } else {
      # mu is a matrix, so apply weighted sum over the samples
      list(mu = c(mu %*% wsample),
           lppd = apply(loglik, 1, log_weighted_mean_exp, wsample))
  }
  
}

.bootstrap_stats <- function(varsel, n_boot = 1000, alpha = 0.05) {

  n <- length(varsel$d_test$y)
  equal_weights <- matrix(1/n, 1, n)

  b_weights <- matrix(rexp(n * n_boot, 1), ncol = n)
  b_weights <- b_weights/rowSums(b_weights)

  hf <- function(x, w) {
    .calc_stats(x$mu, x$lppd, varsel$d_test, varsel$family_kl, w)
  }

  sub <- lapply(varsel$summaries$sub, function(x) {
    list(stats = hf(x, equal_weights),  boot = hf(x, b_weights))
  })

  full <- list(stats = hf(varsel$summaries$full, equal_weights),
               boot = hf(varsel$summaries$full, b_weights))

  # apply over submodel sizes
  quantiles <- mapply(function(sub, size) {
    # apply over different stats
    mapply(function(name, sub_boot, full_boot, sub_statistic, full_statistic) {
      qs <- quantile(sub_boot, c(alpha/2, 1-alpha/2))
      qs_delta <- quantile(sub_boot - full_boot, c(alpha/2, 1-alpha/2))

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

