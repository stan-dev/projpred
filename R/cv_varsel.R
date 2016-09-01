#' Variable selection for generalized linear models with cross-validation
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit A \link[=stanreg-objects]{stanreg} object.
#' @param fits An array with cross-validated stanfits and the respective
#' test datasets returned by \link[=stanreg-objects]{cv_fit}(fit).
#' If not provided, \link[=stanreg-objects]{cv_fit}(fit) is called to
#' get the array.
#' @param ... Optional arguments. Possible arguments and their defaults are:
#' \describe{
#'  \item{\code{ns = min(400, [number of samples])}}{
#'    Number of samples used in the variable selection.
#'    Cannot be larger than the number of samples.}
#'  \item{\code{nc = 0}}{
#'    If nonzero, samples are clustered and the cluster centers are
#'    used in the variable selection instead of the actual samples.
#'  }
#'  \item{\code{nv = min(ncol(x) - 1, rankMatrix(x))}}{
#'    Maximum number of features to be used in the projection (incl. intercept).
#'    Cannot be larger than \code{min(ncol(x) - 1, rankMatrix(x))}.}
#'  \item{\code{verbose = FALSE}}{
#'    If \code{verbose = TRUE}, prints information about the progress of the
#'    variable selection.}
#' }
#'
#' @return A list with class \code{'varsel'} containing the following elements:
#' \describe{
#'  \item{\code{chosen}}{The order in which the features were added to the submodel.}
#'  \item{\code{pctch}}{Percentage of cross-validation runs that included the given
#'    variable to a model of given size.}
#'  \item{\code{stats}}{An array with statistics of the submodel performance.}
#'  \item{\code{family}}{A \code{\link{family}}-object.}
#' }
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fits <- cv_fit(fit)
#' vars <- cv_varsel(fit, fits)
#' plot(vars)
#' }
#'

#' @export
cv_varsel <- function(fit, fits = NA, ...) {
  UseMethod('cv_varsel')
}

#' @export
cv_varsel.stanreg <- function(fit, fits = NA, ...) {

  if(!is.list(fits)) fits <- cv_fit(fit)
  k <- nrow(fits)
  verbose <- ifelse(is.null(list(...)$verbose), F, list(...)$verbose)
  vars <- .extract_vars(fit)
  family_kl <- kl_helpers(family(fit))
  if(ncol(vars$x) < 2)
    stop('Data must have at least 2 features.')
  if(!(family_kl$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(family_kl$family, 'family not yet supported.'))
  # max number of variables
  cv_nv <- min(sapply(c(list(fit), fits[,'fit']), function(f) rankMatrix(get_x(f))))

  msgs <- paste('Forward selection for the',
                c('full model.', paste0('fold number ', 1:k,'/',k,'.')))

  # perform forward selection
  sel <- mapply(function(fit, d_test, msg, cv_nv, verbose) {
    if(verbose) print(msg)
    varsel(fit, d_test, ..., cv = T, cv_nv = cv_nv)
  }, c(list(fit = fit), fits[,'fit']),
  c(list(d_test = NA), fits[,'d_test']), msgs, MoreArgs = list(cv_nv, verbose))

  # combine cross-validation results
  combcv <- function(x) as.list(data.frame(apply(x, 1, function(x) do.call(c, x))))
  d_cv <- combcv(simplify2array(fits[,'d_test'])[c('y','w','offset'),])
  mu_cv <- combcv(simplify2array(sel['mu',-1]))
  lppd_cv <- combcv(simplify2array(sel['lppd',-1]))

  # evaluate performance on test data and
  # use bayesian bootstrap to get 5% credible intervals
  n <- length(d_cv$y)
  nvs <- c(1:(length(sel[['mu',2]])-1), ncol(vars$x)) - vars$intercept
  n_boot <- 1000
  b_weights <- matrix(rexp(n * n_boot, 1), ncol = n)
  b_weights <- b_weights/rowSums(b_weights)
  test_stats <- .bootstrap_stats(mu_cv, lppd_cv, d_cv, nvs, family_kl, b_weights, 'test')

  # find out how many of cross-validated forward selection iterations select
  # the same variables as the forward selection that uses all the data.
  sub_chosen <- do.call(cbind, sel['chosen',-1])
  res <- list(chosen = sel[['chosen',1]])
  res$pctch <- mapply(function(var_ind, ind, arr, k) sum(arr[1:ind, ] == var_ind)/k,
                      res$chosen, seq_along(res$chosen), MoreArgs = list(sub_chosen, k))

  res$stats <- rbind(sel[['stats',1]], test_stats, make.row.names = F)
  res$family <- family(fit)
  if(!is.null(list(...)$nc) && list(...)$nc > 0) res$cl <- sel[['cl',1]]

  structure(res, class = 'varsel')
}
