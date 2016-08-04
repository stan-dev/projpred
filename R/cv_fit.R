#' @export

cv_fit <- function(fit, k = 10) {
  stopifnot(!is.null(fit$data), nrow(fit$data) >= k)
  n <- nobs(fit)
  params <- extract_params(fit)

  t(sapply(1:k, function(i, fit, k, n, d) {
    print(paste0('Fitting model ', i, ' out of ', k))
    i_test <- seq(i, n, k)
    list(fit = update(fit, subset = setdiff(1:n, i_test), refresh = 0),
         d_test = list(x = d$x[i_test,], y = d$y[i_test], w = d$w[i_test]))
  }, fit, k, n, params))
}

