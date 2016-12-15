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

#
log_weighted_mean_exp <- function(x, w) {
  x <- x + log(w)
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}
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

  # undo the random permutation to make results reproducible
  perm_inv <- c(mapply(function(p, i) order(p) + i*length(p),
                       fit$stanfit@sim$permutation,1:fit$stanfit@sim$chains-1))
  res <- list(
    x = unname(get_x(fit)),
    alpha = unname(drop(e$alpha %ORifNULL% rep(0, NROW(e$beta))))[perm_inv],
    beta = t(unname(drop(e$beta)))[, perm_inv],
    dis = unname(e[['dispersion']]) %ORifNULL% rep(1, nrow(e$beta))[perm_inv],
    offset = fit$offset %ORifNULL% rep(0, nobs(fit)),
    intercept = attr(fit$terms,'intercept') %ORifNULL% F)

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
    ns_total = ncol(vars$beta), # number of samples available
    rank_x = rankMatrix(vars$x), # nv is set to <= rank_x
    ns = min(args$ns %ORifNULL% 400, ncol(vars$beta)),
    nc = min(args$nc %ORifNULL% 0, 40), # number of clusters, if samples are clustered
    n_boot = args$n_boot %ORifNULL% 1000, # bootstrap sample size
    intercept = vars$intercept %ORifNULL% F,
    verbose = args$verbose %ORifNULL% F,
    cv = args$cv %ORifNULL% F, # was function called from cv_varsel?
    regul = args$regul %ORifNULL% 1e-12 # small regul as in Dupuis & Robert
  )
  res$clust <- res$nc > 0
  if(!is.null(args$nc) && args$nc > res$nc)
    print(paste0('Setting the number of clusters to ', res$nc, '.'))

  if(!is.null(args$ns) && args$ns > res$ns)
    print(paste0('Setting the number of samples to ', res$ns, '.'))

  res$nv <- min(ncol(vars$x), args$nv, res$rank_x)
  if(!is.null(args$nv) && args$nv > res$nv)
    print(paste0('Setting the max number of variables
                 in the projection to ', res$nv, '.'))

  res
}

# perform clustering over the samples
.get_p_clust <- function(mu, dis, nc=1, wsample=rep(1,dim(mu)[2]), cl = NULL) {

  # TODO
  # THIS FUNCTION WORKS CURRENTLY ONLY FOR GAUSSIAN FAMILY.
  # SHOULD TAKE FAMILY AS AN INPUT AND ACT ACCORDINGLY

  # cluster the mu-samples if no clustering provided
  cl <- cl %ORifNULL% kmeans(t(mu), nc, iter.max = 50)

  if (typeof(cl)=='double') {
      # only cluster-indices provided, so create the list and put them there
      cl <- list(cluster=cl)
  }

  # (re)compute the cluster centers, because they may be different from the ones
  # returned by kmeans if the samples have differing weights
  nc <- max(cl$cluster) # number of clusters (assumes labeling 1,...,nc)
  centers <- matrix(0, nrow=nc, ncol=dim(mu)[1])
  wcluster <- rep(0,nc) # cluster weights
  for (j in 1:nc) {
      ind <- which(cl$cluster==j)
      ws <- wsample[ind]/sum(wsample[ind]) # normalized sample weights within the cluster
      centers[j,] <- mu[,ind] %*% ws
      wcluster[j] <- sum(wsample[ind]) # unnormalized weight for the jth cluster
  }
  cl$centers <- centers
  wcluster <- wcluster/sum(wcluster)

  # compute the dispersion parameters for each cluster
  disps <- sapply(1:nc,
                  function(cl_ind) {
                    ind <- which(cl$cluster== cl_ind)
                    ws <- wsample[ind]/sum(wsample[ind]) # normalized sample weights within the cluster
                    if (length(ind) > 1) {
                        mu_mean <- mu[,ind] %*% ws
                        mu_var <- mu[,ind]^2 %*% ws - mu_mean^2
                        sqrt( sum(ws*dis[ind]^2) + mean(mu_var) )
                    } else
                        sqrt(sum(ws*dis[ind]^2))
                  })

  # combine the results
  p <- list(mu = unname(t(cl$centers)),
            dis = disps,
            weights = wcluster)
  list(cl = cl, p = p)
}



.split_coef <- function(b, intercept) {
  if(intercept) {
    list(alpha = b[1, ], beta = b[-1, , drop = F])
  } else {
    list(alpha = rep(0, NCOL(b)), beta = b)
  }
}





# calculate everything that needs to be saved from the submodel
.summary_stats <- function(chosen, d_train, d_test, p_full, family_kl,
                           intercept, regul, coef_init, coef_full) {

  projfun <- .get_proj_handle(family_kl)

  # helper for re-evaluating dis and kl on test data
  mu_sub_kl <- function(mu, p_full, d_test) {
    ns <- NCOL(p_full$mu)
    dis <- family_kl$dis_fun(p_full, d_test, list(mu = mu))
    kl <- sapply(1:NCOL(p_full$mu), function(x) {
      family_kl$kl(p_full = list(mu = p_full$mu[,x], dis = p_full$dis[x]),
                   d_test, p_sub = list(mu = mu[,x], dis = dis[x]))
    })
    list(kl = kl, mu = mu, dis = dis)
  }

  p_sub <- sapply(seq_along(chosen), function(x) {
    res <- projfun(chosen[1:x], p_full, d_train, intercept, regul, coef_init)
    mu_sub <- family_kl$mu_fun(d_test$x[,chosen[1:x]], res$alpha,
                               res$beta[1:x,], d_test$offset, intercept)
    mu_sub_kl(mu_sub, p_full, d_test)
  })

  # do the projection also without any variables
  if(intercept) {
    p_null <- projfun(0, p_full, d_train, 1, regul, coef_init)
    mu_null <- family_kl$mu_fun(d_test$x[,0], p_null$alpha, p_null$beta[0,],
                                      d_test$offset, intercept)
  } else {
    mu_null <- family_kl$linkinv(matrix(0, NROW(p_full$mu), NCOL(p_full$mu)))
  }
  # add null model to the array of submodels
  p_sub <- cbind(mu_sub_kl(mu_null, p_full, d_test), p_sub)

  kl_list <- unname(unlist(p_sub['kl',]))

  mu_full <- family_kl$linkinv(d_test$offset +
    cbind(1, d_test$x)%*%rbind(coef_full$alpha, coef_full$beta))

  dis_rep <- function(x) {
    matrix(rep(x, each = length(d_test$y)), ncol = NCOL(mu_full))
  }
  lppd_fun <- function(mu, dis) {
    apply(family_kl$ll_fun(mu, dis, d_test$y, d_test$weights), 1,
          log_weighted_mean_exp, p_full$weights)
  }
  lppd_sub <- mapply(lppd_fun, p_sub['mu',], lapply(p_sub['dis',], dis_rep),
                     SIMPLIFY = F)
  # should somehow have dis calculated with test-data?
  lppd_full <- lppd_fun(mu_full, dis_rep(p_full$dis))

  # helper to avg mu and kl by sample weights
  avg_ <- function(x) c(x%*%p_full$weights)

  list(kl = sapply(p_sub['kl',], avg_),
       sub = list(mu = lapply(p_sub['mu',], avg_), lppd = lppd_sub),
       full = list(mu = avg_(mu_full), lppd = lppd_full))
}

# get bootstrapped 95%-intervals for the estimates
.bootstrap_stats <- function(stats, nv_list, d_test, family_kl, b_weights,
                             eval_data, intercept) {

  # calculate the bootstrap samples
  bh_fun <- function(mu, lppd, nv) {
    c(.bootstrap_helper(mu, lppd, d_test, family_kl, b_weights), nv = nv)
  }

  res_sub <- mapply(bh_fun, stats$sub$mu, stats$sub$lppd, nv_list, SIMPLIFY = F)
  res_full <- bh_fun(stats$full$mu, stats$full$lppd, NCOL(d_test$x) + intercept)

  # get the quantiles from the bootstrap samples
  res_quantiles <- lapply(res_sub, function(res_sub) {
    mapply(function(name, boot_stats, boot_stats_full, stats, stats_full){
      qs <- quantile(boot_stats, c(0.025, 0.975))
      qs_delta <- quantile(boot_stats - boot_stats_full, c(0.025, 0.975))

      data.frame(data = eval_data, size = res_sub$nv,
                 delta = c(F, T), summary = rep(name, 2),
                 value = c(stats, stats - stats_full),
                 lq = c(qs[1], qs_delta[1]), uq = c(qs[2], qs_delta[2]))
    }, names(res_sub$boot_stats), res_sub$boot_stats, res_full$boot_stats,
    res_sub$stats, res_full$stats, SIMPLIFY = F)
  })

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

  list(stats = unlist(stats), boot_stats = boot_stats)
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

