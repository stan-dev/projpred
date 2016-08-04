test_stat <- function(params, chosen, d_test, family) {
  if(!is.list(d_test)) return(rep(NA_real_, 4))
  mu <- family$linkinv(d_test$x[, chosen, drop = F]%*%params$b)
  mu_mean <- rowMeans(mu)*d_test$w
  ss_res <- sum((d_test$y - mu_mean)^2)
  lppds <- family$ll_fun(mu, params$dis, d_test$y, d_test$w)
  n <- length(mu_mean)

  c(mse = ss_res/n,
    mlpd = mean(apply(as.matrix(lppds), 2, log_mean_exp)),
    r2 = ifelse(family$family == 'gaussian', 1 - sum((d_test$y - mu_mean)^2)/sum((d_test$y - mean(d_test$y))^2), NA),
    pctcorr = ifelse(family$family == 'binomial', sum(round(mu_mean) == d_test$y)/n, NA))
}

# log likelihoods
log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x))) - log(length(x))
}
