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

log_weighted_mean_exp <- function(x, w) {
  x <- x + log(w)
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}
# Updated version of the kfold function in the rstanarm-package

#' @export
kfold_ <- function (x, K = 10, save_fits = FALSE)
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
    dis = unname(e[['dispersion']]) %ORifNULL% rep(NA, nrow(e$beta))[perm_inv],
    offset = fit$offset %ORifNULL% rep(0, nobs(fit)),
    intercept = attr(fit$terms,'intercept') %ORifNULL% 0)

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

.split_coef <- function(b, intercept) {
  if(intercept) {
    list(alpha = b[1, ], beta = b[-1, , drop = F])
  } else {
    list(alpha = rep(0, NCOL(b)), beta = b)
  }
}

.varsel_errors <- function(e) {
  if(grepl('computationally singular', e$message)) {
    stop(paste(
      'Numerical problems with inverting the covariance matrix. Possibly a',
      'problem with the convergence of the stan model?, If not, consider',
      'stopping the selection early by setting the variable nv_max accordingly.'
    ))
  } else {
    stop(e$message)
  }
}

