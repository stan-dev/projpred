.onAttach <- function(...) {
  ver <- utils::packageVersion("glmproj")
  packageStartupMessage("This is glmproj version ", ver)
}

# from rstanarm
log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x))) - log(length(x))
}

# from rstanarm
`%ORifNULL%` <- function(a, b) if (is.null(a)) b else a

# extract all important 'information' from a stanreg-object
# that is needed in the variable selection
.extract_vars <- function(fit) {
  e <- extract(fit$stanfit)

  res <- list(x = unname(get_x(fit)),
              b = t(unname(cbind(e$alpha,e$beta))),
              dis = unname(e$dispersion) %ORifNULL% rep(1, nobs(fit)),
              offset = fit$offset %ORifNULL% rep(0, nobs(fit)),
              intercept = attr(fit$terms,'intercept'))
  attr(res$x, 'assign') <- NULL

  y <- unname(get_y(fit))
  if(NCOL(y) == 1) {
    res$w <- if(length(weights(fit))) unname(weights(fit)) else rep(1, nobs(fit))
    res$y <- y
  } else {
    res$w <- rowSums(y)
    res$y <- y[, 1] / res$w
  }

  res
}

# get bootstrapped 95%-intervals for the estimates
.bootstrap_stats <- function(mus, lppds, d_test, nvs, family_kl, b_weights, eval_data) {
  nv_ind <- length(nvs)
  do.call(rbind, c(
    mapply(function(mu, lppd, nv, mu_full, lppd_full, d_test, family_kl, b_weights) {
      cbind(data = eval_data, nvar = nv,
            .bootstrap_helper(mu, lppd, mu_full, lppd_full, d_test, family_kl, b_weights))
      }, mus, lppds, nvs, SIMPLIFY = F,
      MoreArgs = list(mus[[nv_ind]], lppds[[nv_ind]], d_test, family_kl, b_weights)),
    make.row.names = F))
}

.bootstrap_helper <- function(mu, lppd, mu_full, lppd_full, d_test, family_kl, b_weights = NULL, n_boot = 1000) {

  n <- length(d_test$y)
  if(is.null(b_weights)) {
    b_weights <- matrix(rexp(n * n_boot, 1), ncol = n)
    b_weights <- b_weights/rowSums(b_weights)
  }

  arr <- rbind(elpd = lppd, mse = (d_test$y-mu)^2,
               elpd_delta = lppd-lppd_full,
               mse_delta = (mu+mu_full-2*d_test$y)*(mu-mu_full)) #=(y-mu)^2-(y-mu_full)^2
  # stats are the actual values, boot_stats are the bootstrap samples
  stats <- as.list(rowMeans(arr))
  boot_stats <- lapply(seq_len(length(stats)),
                       function(ind, tc) tc[ind, ], tcrossprod(arr, b_weights))
  names(boot_stats) <- names(stats)

  # mcfadden pseudo r2
  lppd_null <- family_kl$ll_fun(mean(d_test$y), var(d_test$y)*(n-1)/n, d_test$y, d_test$w)
  lppd_null_m <- mean(lppd_null)
  lppd_null_boot <- tcrossprod(lppd_null, b_weights)
  stats$r2 <- 1 - stats$elpd/lppd_null_m
  boot_stats$r2 <- 1 - boot_stats$elpd/lppd_null_boot
  stats$r2_delta <- stats$elpd_delta/lppd_null_m
  boot_stats$r2_delta <- boot_stats$elpd_delta/lppd_null_boot

  if(family_kl$family == 'binomial') {
    stats$pctcorr <- mean(round(d_test$w*mu) == d_test$w*d_test$y)
    stats$pctcorr_delta <- stats$pctcorr - mean(round(d_test$w*mu_full) == d_test$w*d_test$y)
    boot_stats$pctcorr <- tcrossprod(round(d_test$w*mu) == d_test$w*d_test$y, b_weights)
    boot_stats$pctcorr_delta <- tcrossprod((round(d_test$w*mu) == d_test$w*d_test$y) -
                                             (round(d_test$w*mu_full) == d_test$w*d_test$y), b_weights)
  }

  do.call(rbind, c(
    mapply(function(boot_x, x, name) {
      data.frame(delta = grepl('_delta', name), summary = gsub('_delta', '', name),
                 value = x, lq = quantile(boot_x, 0.025), uq = quantile(boot_x, 0.975))
    }, boot_stats, stats, names(stats), SIMPLIFY = F), make.row.names = F))
}

# initialize arguments to their default values if they are not specified
.init_args <- function(args, vars, fam) {
  res <- list(
    ns_total = ncol(vars$b),
    rank_x = rankMatrix(vars$x),
    ns = min(args$ns %ORifNULL% 400, ncol(vars$b)),
    nc = args$nc %ORifNULL% 0,
    intercept = vars$intercept,
    verbose = args$verbose %ORifNULL% FALSE,
    cv = args$cv %ORifNULL% F,
    family_kl = kl_helpers(fam)
  )
  res$clust <- res$nc > 0
  res$nv <- min(ncol(vars$x) - 1, args$nv, res$cv_nv, res$rank_x)

  res
}

# perform clustering over the samples
.get_p_clust <- function(mu, dis, args, cl = NULL) {
  # calculate the means of the variables
  cl <- cl %ORifNULL% kmeans(t(mu), args$nc, iter.max = 50)
  p <- list(mu = unname(t(cl$centers)),
            dis = sapply(1:args$nc, function(cl_ind, dis, cl_assign) {
              sqrt(mean(dis[which(cl_assign == cl_ind)]^2))
            }, dis, cl$cluster),
            cluster_w = cl$size)
  list(cl = cl, p = p)
}
