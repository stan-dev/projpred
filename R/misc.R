.onAttach <- function(...) {
  ver <- utils::packageVersion("glmproj")
  packageStartupMessage("This is glmproj version ", ver)
}

# from rstanarm
log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x))) - log(length(x))
}
is.stanreg <- function(x) inherits(x, "stanreg")

# Updated version of the kfold function in the rstanarm-package

#' @export
kfold <- function (x, K = 10, save_fits = FALSE)
{
  #validate_stanreg_object(x)
  #if (!used.sampling(x))
  #  STOP_sampling_only("kfold")
  #stopifnot(!is.null(x$data), K > 1, nrow(x$data) >= K)
  d <- x$data
  N <- nrow(d)
  wts <- x[["weights"]]
  perm <- sample.int(N)
  idx <- ceiling(seq(from = 1, to = N, length.out = K + 1))
  bin <- .bincode(perm, breaks = idx, right = FALSE, include.lowest = TRUE)
  lppds <- list()
  fits <- array(list(), c(K, 2), list(NULL, c('fit','omitted')))
  for (k in 1:K) {
    message("Fitting model ", k, " out of ", K)
    omitted <- which(bin == k)
    fit_k <- update(object = x, data = d[-omitted, ], weights = if (length(wts))
      wts[-omitted]
      else NULL, refresh = 0)
    lppds[[k]] <- log_lik(fit_k, newdata = d[omitted, ])
    if(save_fits) fits[k,] <- list(fit = fit_k, omitted = omitted)
  }
  elpds <- unlist(lapply(lppds, function(x) {
    apply(x, 2, log_mean_exp)
  }))
  out <- list(elpd_kfold = sum(elpds),
              se_elpd_kfold = sqrt(N * var(elpds)),
              pointwise = cbind(elpd_kfold = elpds))
  if(save_fits) out$fits <- fits
  structure(out, class = c("kfold", "loo"), K = K)
}

# check if the fit object is suitable for variable selection
.validate_for_varsel <- function(fit) {
  if(!is.stanreg(fit))
    stop('Object is not a stanreg object')

  if(!(gsub('rstanarm::', '', fit$call[1]) %in% c("stan_glm", "stan_lm")))
    stop('Only \'stan_lm\' and \'stan_glm\' are supported.')

  families <- c('gaussian','binomial','poisson')
  if(!(family(fit)$family %in% families))
    stop(paste0('Only the following families are supported:\n',
                paste(families, collapse = ', '), '.'))

  if(NCOL(get_x(fit)) < 4)
    stop('Not enought explanatory variables for variable selection')
}

# from rstanarm
`%ORifNULL%` <- function(a, b) if (is.null(a)) b else a

# extract all important 'information' from a stanreg object for variable selection
.extract_vars <- function(fit) {
  e <- extract(fit$stanfit)
  dis_name <- switch(family(fit)$family, 'gaussian' = 'sigma', 'Gamma' = 'shape',
                     'dispersion')
  res <- list(x = unname(get_x(fit)),
              b = t(unname(cbind(drop(e$alpha), drop(e$beta)))),
              dis = unname(e[[dis_name]]) %ORifNULL% rep(1, nrow(e$beta)),
              offset = fit$offset %ORifNULL% rep(0, nobs(fit)),
              intercept = attr(fit$terms,'intercept') %ORifNULL% F)

  # undo the random permutation to make results reproducible
  perm_inv <- c(mapply(function(p, i) order(p) + i*length(p),
                       fit$stanfit@sim$permutation, 1:fit$stanfit@sim$chains - 1))
  res$b <- res$b[, perm_inv]
  res$dis <- res$dis[perm_inv]

  res$x <- res$x[, as.logical(attr(res$x, 'assign'))]
  attr(res$x, 'assign') <- NULL

  y <- unname(get_y(fit))
  if(NCOL(y) == 1) {
    res$weights <- if(length(weights(fit))) unname(weights(fit)) else rep(1, nobs(fit))
    res$y <- y
  } else {
    res$weights <- rowSums(y)
    res$y <- y[, 1] / res$weights
  }

  res
}

# initialize arguments to their default values if they are not specified
.init_args <- function(args, vars, fam) {
  res <- list(
    ns_total = ncol(vars$b), # number of samples available
    rank_x = rankMatrix(vars$x), # nv is set to <= rank_x
    ns = min(args$ns %ORifNULL% 400, ncol(vars$b)),
    nc = min(args$nc %ORifNULL% 0, 40), # number of clusters, if samples are clustered
    n_boot = args$n_boot %ORifNULL% 1000, # bootstrap sample size
    intercept = vars$intercept %ORifNULL% F,
    verbose = args$verbose %ORifNULL% F,
    cv = args$cv %ORifNULL% F, # was function called from cv_varsel?
    regul = args$regul %ORifNULL% 1e-12, # small regul as in Dupuis & Robert
    max_it = args$max_it %ORifNULL% 300, # max IRLS steps
    epsilon = args$epsilon %ORifNULL% 1e-8, # used to determine if IRLS has converged
    family_kl = kl_helpers(fam)
  )
  res$clust <- res$nc > 0
  if(!is.null(args$nc) && args$nc > res$nc)
    print(paste0('Setting the number of clusters to ', res$nc, '.'))

  if(!is.null(args$ns) && args$ns > res$ns)
    print(paste0('Setting the number of samples to ', res$ns, '.'))

  res$nv <- min(ncol(vars$x) - res$intercept, args$nv, res$rank_x)
  if(!is.null(args$nv) && args$nv > res$nv)
    print(paste0('Setting the max number of variables
                 in the projection to ', res$nv, '.'))

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
            cluster_w = cl$size/sum(cl$size))
  list(cl = cl, p = p)
}

# function handle for the projection over samples. Gaussian case
# uses analytical solution to do the projection over samples.
.get_proj_handle <- function(family_kl) {

  # Use analytical solution for gaussian as it is a lot faster
  if(family_kl$family == 'gaussian' && family_kl$link == 'identity') {
    function(chosen, p_full, d_train, b0, args) {
      w <- sqrt(d_train$weights)

      regulvec <- c((1-args$intercept)*args$regul, rep(args$regul, length(chosen) - 1))
      regulmat <- diag(regulvec, length(regulvec), length(regulvec))
      # Solution for the gaussian case (with l2-regularization)
      p_sub <- list(b = solve(crossprod(w*d_train$x[,chosen, drop = F]) + regulmat,
                              crossprod(w*d_train$x[,chosen, drop = F], w*p_full$mu)))
      p_sub$dis <- sqrt(colMeans(d_train$weights*(
        p_full$mu - d_train$x[, chosen, drop = F]%*%p_sub$b)^2) + p_full$dis^2)
      p_sub$kl <- weighted.mean(log(p_sub$dis) - log(p_full$dis) + colSums(p_sub$b^2*regulvec), p_full$cluster_w)
      p_sub
    }

  } else {
    function(chosen, p_full, d_train, b0, args) {

      # perform the projection over samples
      res <- sapply(1:ncol(p_full$mu), function(s_ind) {
        IRLS(list(mu = p_full$mu[, s_ind, drop = F], dis = p_full$dis[s_ind]),
             list(x = d_train$x[, chosen, drop = F], weights = d_train$weights,
                  offset = d_train$offset), b0[chosen,], args)
      })

      # weight the results by sample weights (that are all 1 unless p_clust is used)
      p_sub <- list(kl = weighted.mean(unlist(res['kl',]), p_full$cluster_w),
                    b = do.call(cbind, res['b',]))
      if('dis' %in% rownames(res)) p_sub$dis <- unlist(res['dis',])
      p_sub
    }
  }
}

# calculate everything that needs to be saved from the submodel
.summary_stats <- function(chosen, d_train, d_test, p_full, b0, args) {

  projfun <- .get_proj_handle(args$family_kl)

  p_sub <- sapply(seq_along(chosen), function(x) {
    res <- projfun(chosen[1:x], p_full, d_train, b0, args)
    res$mu <- args$family_kl$linkinv(d_test$x[,chosen[1:x]]%*%res$b + d_test$offset)
    res$b <- NULL
    res
  })

  kl_list <- c(unlist(p_sub['kl',]), 0)

  mu_full <- args$family_kl$linkinv(d_test$x%*%p_full$b+ d_test$offset)
  mu_list <- c(p_sub['mu',], list(mu_full))

  dis_full <- p_full$dis %ORifNULL% rep(1, length(d_test$y))
  dis_list <- lapply(c(p_sub['dis',], list(dis_full)), function(x) {
    matrix(rep(x, each = length(d_test$y)), ncol = NCOL(mu_full))
  })

  lppd_list <- mapply(function(mu, dis) {
    apply(args$family_kl$ll_fun(mu, dis, d_test$y, d_test$weights), 1, log_mean_exp)
  }, mu_list, dis_list, SIMPLIFY = F)

  list(kl = kl_list,
       mu = lapply(mu_list, rowMeans),
       dis = lapply(dis_list, rowMeans),
       lppd = lppd_list)
}

# get bootstrapped 95%-intervals for the estimates
.bootstrap_stats <- function(mu_all, lppd_all, nv, d_test, family_kl, b_weights, data) {

  # calculate the bootstrap samples
  res_boot <- mapply(function(mu, lppd, nv, d_test, family_kl, b_weights) {
    c(.bootstrap_helper(mu, lppd, d_test, family_kl, b_weights), nv = nv)
  }, mu_all, lppd_all, nv, MoreArgs = list(d_test, family_kl, b_weights), SIMPLIFY = F)

  # get the quantiles from the bootstrap samples
  res_quantiles <- lapply(res_boot, function(res_boot, res_full, data) {
    mapply(function(size, name, stat, boot_stat, stat_full, boot_stat_full, nv, data) {
      qs <- quantile(boot_stat, c(0.025, 0.975))
      qs_delta <- quantile(boot_stat - boot_stat_full, c(0.025, 0.975))
      data.frame(data = data, size = size, delta = c(F, T),
                 summary = rep(name, 2), value = c(stat, stat - stat_full),
                 lq = c(qs[1], qs_delta[1]), uq = c(qs[2], qs_delta[2]))
    }, res_boot$nv, names(res_boot$stats), res_boot$stats, res_boot$boot_stats,
    res_full$stats, res_full$boot_stats, MoreArgs = list(nv, data), SIMPLIFY = F)
  }, res_boot[[length(res_boot)]], data)

  # rbind the elements into one data.frame
  do.call(rbind, c(unlist(res_quantiles, recursive = F), make.row.names = F))
}

.bootstrap_helper <- function(mu, lppd, d_test, family_kl, b_weights) {

  y <- d_test$y
  weights <- d_test$weights
  n <- length(y)
  arr <- rbind(mlpd = lppd, mse = (y-mu)^2)
  # stats are the actual values, boot_stats are the bootstrap samples
  stats <- as.list(rowMeans(arr))
  boot_stats <- lapply(seq_len(length(stats)),
                       function(ind, tc) tc[ind, ], tcrossprod(arr, b_weights))
  names(boot_stats) <- names(stats)

  # McFadden's pseudo r2
  lppd_null <- -2*family_kl$dev.resids(y, sum(weights * y)/sum(weights), weights)
  stats$r2 <- 1 - sum(lppd)/sum(lppd_null)
  boot_stats$r2 <- drop(1 - (b_weights%*%lppd)/(b_weights%*%lppd_null))

  if(family_kl$family == 'binomial') {
    stats$pctcorr <- mean(round(weights*mu) == weights*y)
    boot_stats$pctcorr <- drop(b_weights%*%(round(weights*mu) == weights*y))
  }

  list(stats = stats, boot_stats = boot_stats)
}

.gen_bootstrap_ws <- function(n_obs, n_boot) {
  b_weights <- matrix(rexp(n_obs * n_boot, 1), ncol = n_obs)
  b_weights/rowSums(b_weights)
}

.varsel_errors <- function(e) {
  if(grepl('computationally singular', e$message)) {
    stop(paste(
      'Numerical problems with inverting the covariance matrix. Possibly a',
      'problem with the convergence of the stan model?, If not, consider adding',
      'a small value to the diagonal elements by setting e.g. regul = 1e-10 or',
      'stopping the variable selection early by setting the variable nv accordingly.'
    ))
  } else {
    stop(e$message)
  }
}

