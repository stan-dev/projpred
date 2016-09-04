#' K-fold cross-validation for models fitted with rstanarm.
#'
#'Do a 10-fold cross-validation for a \link[=stanreg-objects]{stanreg}-object.
#'
#' @param \code{fit} A \link[=stanreg-objects]{stanreg}-object.
#' @param \code{k} Number of cv-folds. Defaults to 10.
#'
#' @return An array with each row being one fold and first column containing
#' the fitted stan-model and the second column containing the corresponding testset.
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fits <- cv_fit(fit)
#' # do cross-validated variable selection
#' vars <- cv_varsel(fit, fits)
#' }
#'

#' @export
cv_fit <- function(fit, k = 10) {
  UseMethod('cv_fit')
}

#' @export
cv_fit.stanreg <- function(fit, k = 10) {
  n <- nobs(fit)
  if(k > n) stop('Not enough observations for cross-validation.')
  rperm <- sample(1:n)

  t(sapply(1:k, function(i, fit, k, n, d_full, rperm) {
    print(paste0('Fitting model ', i, ' out of ', k))
    i_test <- rperm[seq(i, n, k)]
    d_test <- list(x = d_full$x[i_test,], y = d_full$y[i_test],
                   w = d_full$w[i_test], offset = d_full$offset[i_test])
    list(fit = update(fit, subset = !(1:n %in% i_test), refresh = 0),
         d_test = d_test)
  }, fit, k, n, .extract_vars(fit), rperm))
}

