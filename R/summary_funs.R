#' Functions for calculating the summary statistics
#' - .get_sub/full_summaries return mu, and lppd
#' - .calc/bootstrap_metrics return mse, mlpd, r2 (and pctcorr)
#' - .get_kl_array returns kl in an array similar to bootstrap*

.get_sub_summaries <- function(chosen, nv, d_train, d_test, p_full, family_kl, intercept) {

  # TODO FILL IN HOW TO SELECT WHICH SUMMARIES ARE COMPUTED (AND IMPLEMENT MSE, ACCURACY ETC.)
  # replace with get_sub_summaries2?

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

.get_sub_summaries2 <- function(chosen, p_full, data, p_sub, family_kl, intercept) {

  res <- lapply(p_sub, function(p_sub) {
    if(NROW(p_sub$beta) == 0) {
      xt <- matrix(0, nrow = length(data$weights), ncol = 0)
    } else if(!is.matrix(data$x)) {
      xt <- matrix(data$x[inds], nrow = 1)
    } else {
      xt <- data$x[, chosen[1:NROW(p_sub$beta)], drop = F]
    }

    mu <- family_kl$mu_fun(xt, p_sub$alpha, p_sub$beta, data$offset,
                           intercept)
    .weighted_summary_means(data, family_kl, p_full, mu, p_sub$dis)
  })
}


.get_full_summaries <- function(p_full, data, coef_full, family_kl, intercept) {
  mu <- family_kl$mu_fun(data$x, coef_full$alpha, coef_full$beta, data$offset,
                         intercept)
  .weighted_summary_means(data, family_kl, p_full, mu, p_full$dis)
}


# Calculates weighted means of mu, dis, kl and lppd given samples of
# mu and dis, the full model and the data.
.weighted_summary_means <- function(data, family_kl, p_full, mu, dis) {
  loglik <- family_kl$ll_fun(mu, dis, data$y)
  list(mu = c(mu%*%p_full$weights),
       lppd = apply(loglik, 1, log_weighted_mean_exp, p_full$weights))
}


.bootstrap_metrics <- function(sub_summaries, full_summaries, data, family_kl,
                               intercept, is_test, b_weights, alpha = 0.05) {

  equal_weights <- matrix(1/NROW(data$x), 1, NROW(data$x))

  sub <- lapply(sub_summaries, function(x) {
    list(metrics = .calc_metrics(x$mu, x$lppd, data, family_kl, equal_weights),
         boot = .calc_metrics(x$mu, x$lppd, data, family_kl, b_weights))
  })

  full <- list(metrics = .calc_metrics(full_summaries$mu, full_summaries$lppd,
                                       data, family_kl, equal_weights),
               boot = .calc_metrics(full_summaries$mu, full_summaries$lppd,
                                    data, family_kl, b_weights))

  # apply over submodel sizes
  quantiles <- mapply(function(sub, size) {
    # apply over different metrics
    mapply(function(name, sub_boot, full_boot, sub_metric, full_metric) {
      qs <- quantile(sub_boot, c(alpha/2, 1-alpha/2))
      qs_delta <- quantile(sub_boot - full_boot, c(alpha/2, 1-alpha/2))

      data.frame(size = rep(size, 2), delta = c(F, T), metric = rep(name, 2),
                 value = c(sub_metric, sub_metric - full_metric),
                 lq = c(qs[1], qs_delta[1]), uq = c(qs[2], qs_delta[2]))
      }, names(sub$boot), sub$boot, full$boot, sub$metrics, full$metrics, SIMPLIFY = F)
    }, sub, seq_along(sub_summaries) - 1, SIMPLIFY = F)

  # rbind the elements into one data.frame and add a column which indicates
  # whether the summaries are calculated from test or training data.
  cbind(data = ifelse(is_test, 'test', 'train'),
        do.call(rbind, c(unlist(quantiles, recursive = F), make.row.names = F)))
}


.get_bootstrap_ws <- function(n_obs, n_boot = 1000) {
  b_weights <- matrix(rexp(n_obs * n_boot, 1), ncol = n_obs)
  b_weights/rowSums(b_weights)
}


.calc_metrics <- function(mu, lppd, data, family_kl, sample_weights) {
  arr <- list(mlpd = lppd, mse = (data$y-mu)^2)
  if(family_kl$family == 'binomial') {
    arr$pctcorr <- (round(data$weights*mu) == round(data$weights*data$y))
  }
  avg_ <- function(x) c(sample_weights%*%x)
  metrics <- lapply(arr, avg_)

  # check the correctness of this
  mlpd_null <- -2*(family_kl$dev.resids(
    data$y, weighted.mean(data$y,data$weights), data$weights))

  metrics$r2 <- 1 - metrics$mlpd/avg_(mlpd_null)

  metrics
}


.get_kl_array <- function(p_sub) {
  data.frame(data = 'sel', size = seq_along(p_sub)-1, delta = F, metric = 'kl',
             value = sapply(p_sub, function(x) x$kl), lq = NA, uq = NA)
}
