.onAttach <- function(...) {
  ver <- utils::packageVersion("glmproj")
  packageStartupMessage("This is glmproj version ", ver)
}

# log likelihoods
log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x))) - log(length(x))
}

bootstrap_summaries <- function(mu, lppd, d, family_kl, b_weights = NULL, n_boot = 1000) {

  n <- length(d$y)
  if(is.null(b_weights)) {
    b_weights <- matrix(rexp(n * n_boot, 1), ncol = n)
    b_weights <- b_weights/rowSums(b_weights)
  }

  stats <- list(elpd = mean(lppd), mse = mean((d$y-mu)^2))
  boot_stats <- list(elpd = tcrossprod(lppd, b_weights), mse = tcrossprod((d$y-mu)^2, b_weights))
  # mcfadden pseudo r2
  lppd_null <- family_kl$ll_fun(mean(d$y), var(d$y)*(n-1)/n, d$y, d$w)
  stats$r2 <- 1 - stats$elpd/mean(lppd_null)
  boot_stats$r2 <- 1 - boot_stats$elpd/tcrossprod(lppd_null, b_weights)
  if(family_kl$family == 'binomial') {
    boot_stats$pctcorr <- tcrossprod(round(d$w*mu) == d$y, b_weights)
    stats$pctcorr <- mean(round(d$w*mu) == d$y)
  }

  aa <- do.call(rbind, c(mapply(function(boot_x, x, name) {
    data.frame(summary = name, value = x, lq = quantile(boot_x, 0.025), uq = quantile(boot_x, 0.975))
  }, boot_stats, stats, names(stats), SIMPLIFY = F), make.row.names = F))
}


extract_params <- function(fit) {
  res <- list()
  res$x <- unname(get_x(fit))
  attr(res$x, 'assign') <- NULL
  res$intercept <- attr(fit$terms,'intercept')

  e <- extract(fit$stanfit)
  res$b <- e$beta
  if(res$intercept) res$b <- unname(cbind(e$alpha, res$b))
  res$b <- t(res$b)
  res$dis <- unname(e$dispersion)
  if(is.null(res$dis)) res$dis <- NA
  res$offset <- fit$offset
  if(is.null(res$offset)) res$offset <- rep(0, nobs(fit))

  # get weights from the fit

  # get y
  y <- unname(get_y(fit))
  res$w <- unname(weights(fit))
  if(length(res$w) == 0) res$w <- rep(1, nobs(fit))
  if(NCOL(y) == 1) {
    res$y <- y
  } else {
    res$y <- y[, 1]
    res$w <- res$w * (y[, 1] + y[, 2])
  }

  res
}

summary_stats <- function(mus, lppds, d, nvs, family_kl, b_weights, eval_data) {
  do.call(rbind, c(
    mapply(function(mu, lppd, nv, d, family_kl, b_weights) {
      cbind(data = eval_data, nvar = nv, bootstrap_summaries(mu, lppd, d, family_kl, b_weights))
    }, mus, lppds, nvs, MoreArgs=list(d, family_kl, b_weights), SIMPLIFY = F),
    make.row.names = F))
}
