# get bootstrapped 95%-intervals for the estimates
.bootstrap_stats <- function(sub_list, full_list, d_test, family_kl, b_weights,
                             intercept, alpha = 0.05) {

  equal_weights <- matrix(1/NROW(d_test$x), 1, NROW(d_test$x))

  res_sub <- lapply(sub_list, function(x) {
    list(stat = .calc_summaries(x$mu, x$lppd, d_test, family_kl, equal_weights),
         boot = .calc_summaries(x$mu, x$lppd, d_test, family_kl, b_weights))
  })

  res_full <- list(stat = .calc_summaries(full_list$mu, full_list$lppd, d_test,
                                          family_kl, equal_weights),
                   boot = .calc_summaries(full_list$mu, full_list$lppd, d_test,
                                          family_kl, b_weights))

  # get the quantiles from the bootstrap samples
  res_quantiles <- lapply(res_sub, function(res_sub) {
    mapply(function(name, size, boot, boot_full, stat, stat_full) {
      qs <- quantile(boot, c(alpha/2, 1-alpha/2))
      qs_delta <- quantile(boot - boot_full, c(alpha/2, 1-alpha/2))

      data.frame(size = rep(size, 2), delta = c(F, T), summary = rep(name, 2),
                 value = c(stat, stat - stat_full),
                 lq = c(qs[1], qs_delta[1]), uq = c(qs[2], qs_delta[2]))
    }, names(res_sub$boot), seq_along(sub_list)-1, res_sub$boot,
    res_full$boot, res_sub$stat, res_full$stat, SIMPLIFY = F)
  })

  # rbind the elements into one data.frame
  do.call(rbind, c(unlist(res_quantiles, recursive = F), make.row.names = F))
}

.gen_bootstrap_ws <- function(n_obs, n_boot = 1000) {
  b_weights <- matrix(rexp(n_obs * n_boot, 1), ncol = n_obs)
  b_weights/rowSums(b_weights)
}

.calc_summaries <- function(mu, lppd, d_test, family_kl, sample_weights) {
  arr <- list(mlpd = lppd, mse = (d_test$y-mu)^2)
  if(family_kl$family == 'binomial') {
    arr$pctcorr <- (round(d_test$weights*mu) == round(d_test$weights*d_test$y))
  }
  avg_ <- function(x) c(sample_weights%*%x)
  res <- lapply(arr, avg_)

  # check the correctness of this
  mlpd_null <- -2*(family_kl$dev.resids(
    d_test$y, weighted.mean(d_test$y,d_test$weights), d_test$weights))

  res$r2 <- 1 - res$mlpd/avg_(mlpd_null)

  res
}
