.get_sub_summaries <- function(chosen, nv, d_train, d_test, p_full, family_kl, intercept) {

  # TODO FILL IN HOW TO SELECT WHICH SUMMARIES ARE COMPUTED (AND IMPLEMENT MSE, ACCURACY ETC.)

  # project onto the given model sizes
  psub <- .get_submodels(chosen, nv, family_kl, p_full, d_train, intercept)

  # compute the summaries on the test data for each sub model
  summaries <- lapply(1:length(nv),
                      function(j) {
                        if (nv[j] == 0)
                          ind <- integer(length=0) # empty
                        else
                          ind <- chosen[1:nv[j]]
                        if (is.null(dim(d_test$x)))
                          # only one test point
                          xt <- matrix(d_test$x[ind], nrow=1)
                        else
                          xt <- d_test$x[,ind,drop=F]

                        mu <- family_kl$mu_fun(xt, psub[[j]]$alpha, psub[[j]]$beta, d_test$offset)
                        loglik <- family_kl$ll_fun(mu, psub[[j]]$dis, d_test$y)
                        lppd <- apply(loglik, 1, log_weighted_mean_exp, p_full$weights)

                        return(list(lppd = lppd))
                      })

  return(summaries)

}

.get_sub_summaries2 <- function(chosen, d_train, d_test, p_full, family_kl,
                                intercept) {
  submodels <- .get_submodels(chosen, c(0, seq_along(chosen)), family_kl,
                              p_full, d_train, intercept)

  res <- lapply(submodels, function(p_sub) {
    if(NROW(p_sub$beta) == 0) {
      xt <- matrix(0, nrow = length(d_test$weights), ncol = 0)
    } else if(!is.matrix(d_test$x)) {
      xt <- matrix(d_test$x[inds], nrow = 1)
    } else {
      xt <- d_test$x[, chosen[1:NROW(p_sub$beta)], drop = F]
    }

    mu <- family_kl$mu_fun(xt, p_sub$alpha, p_sub$beta, d_test$offset)
    .weighted_summary_means(d_test, family_kl, p_full, mu, p_sub$dis)
  })
}

.get_full_summaries <- function(data, p_full, coef_full, family_kl) {
  mu <- family_kl$mu_fun(data$x, coef_full$alpha, coef_full$beta, data$offset)
  .weighted_summary_means(data, family_kl, p_full, mu, p_full$dis)
}

# Calculates weighted means of mu, dis, kl and lppd given samples of
# mu and dis, the full model and the data.
.weighted_summary_means <- function(data, family_kl, p_full, mu, dis) {
  loglik <- family_kl$ll_fun(mu, dis, data$y)
  kl <- family_kl$kl(p_full, data, list(mu = mu, dis = dis))
  avg_ <- function(x) c(x%*%p_full$weights)

  list(mu = avg_(mu), dis = avg_(dis), kl = avg_(kl),
       lppd = apply(loglik, 1, log_weighted_mean_exp, p_full$weights))
}
