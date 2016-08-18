#' @export cv_fit cv_fit.stanreg
cv_fit <- function(fit, k = 10) {
  UseMethod('cv_fit')
}

cv_fit.stanreg <- function(fit, k = 10) {
  n <- nobs(fit)
  if(k > n) stop('Not enough observations for cross-validation.')

  t(sapply(1:k, function(i, fit, k, n, d) {
    print(paste0('Fitting model ', i, ' out of ', k))
    i_test <- seq(i, n, k)
    d_test <- list(x = d$x[i_test,], y = d$y[i_test],
                   w = d$w[i_test], offset = d$offset[i_test])
    list(fit = update(fit, subset = !(1:n %in% i_test), refresh = 0),
         d_test = d_test)
  }, fit, k, n, extract_params(fit)))
}

