.get_sub_summaries_poc <- function(submodels, test_points, refmodel, family_kl, groups=NULL) {

  has_group_features <- !is.null(groups)
  res <- lapply(submodels, function(model) {
  	vind <- model$vind
    if (length(vind) == 0)
      vind <- c("1")
    sub_fit <- model$sub_fit
    mu <- family_kl$mu_fun(sub_fit, obs=test_points)

    weights <- refmodel$wobs[test_points]
    y <- refmodel$y[test_points]
    y_test <- list(y=y, weights=weights)

    .weighted_summary_means_poc(y_test, family_kl, model$weights, matrix(mu, NROW(y), NCOL(mu)), model$dis)
  })
}

.weighted_summary_means_poc <- function(y_test, family_kl, wsample, mu, dis) {

  loglik <- family_kl$ll_fun(mu, dis, matrix(y_test$y, nrow=NROW(mu)), y_test$weights)
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
