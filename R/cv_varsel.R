#' @export cv_varsel cv_varsel.stanreg
cv_varsel <- function(fit, fits = NA, ...) {
  UseMethod('cv_varsel')
}

cv_varsel.stanreg <- function(fit, fits = NA, ...) {

  if(!is.list(fits)) fits <- cv_fit(fit)
  k <- nrow(fits)
  verbose <- ifelse(is.null(list(...)$verbose), F, list(...)$verbose)
  params <- extract_params(fit)
  family_kl <- kl_helpers(family(fit))
  if(ncol(params$x) < 2)
    stop('Data must have at least 2 features.')
  if(!(family_kl$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(family_kl$family, 'family not yet supported.'))
  cv_nv <- min(sapply(fits[,'d_test'], function(d) rankMatrix(d$x) - 1))

  msgs <- paste('Forward selection for the', c('full model.', paste0('fold number ', 1:k,'/',k,'.')))

  # perform forward selection
  sel <- mapply(function(fit, d_test, msg, verbose) {
    if(verbose) print(msg)
    varsel(fit, d_test, ..., cv = T, cv_nv = cv_nv)
  }, c(list(fit = fit), fits[,'fit']), c(list(d_test = NA), fits[,'d_test']), msgs, MoreArgs = list(verbose = verbose))

  # combine cross-validation results
  combcv <- function(x) as.list(data.frame(apply(x, 1, function(x) do.call(c, x))))
  d_cv <- combcv(simplify2array(fits[,'d_test'])[c('y','w','offset'),])
  mu_cv <- combcv(simplify2array(sel['mu',-1]))
  lppd_cv <- combcv(simplify2array(sel['lppd',-1]))

  # evaluate performance on test data and
  # use bayesian bootstrap to get 5% credible intervals
  n <- length(d_cv$y)
  nvs <- c(1:(length(sel[['mu',2]])-1), ncol(params$x)) - params$intercept
  n_boot <- 1000
  b_weights <- matrix(rexp(n * n_boot, 1), ncol = n)
  b_weights <- b_weights/rowSums(b_weights)
  test_stats <- summary_stats(mu_cv, lppd_cv, d_cv, nvs, family_kl, b_weights, 'test')

  # find out how many of cross-validated forward selection iterations select
  # the same variables as the forward selection that uses all the data.
  sub_chosen <- do.call(cbind, sel['chosen',-1])
  res <- list(chosen = sel[['chosen',1]])
  res$pctch <- mapply(function(var_ind, ind, arr, k) sum(arr[1:ind, ] == var_ind)/k,
                      res$chosen, seq_along(res$chosen), MoreArgs = list(sub_chosen, k))

  res$stats <- rbind(sel[['stats',1]], test_stats, make.row.names = F)
  res$family <- family(fit)

  structure(res, class = 'varsel')
}
