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
#'  \item{\code{nc = min(0, ns/4}}{
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
cv_varsel <- function(fit, fits = NULL, ...) {
  UseMethod('cv_varsel')
}

#' @export
cv_varsel.stanreg <- function(fit, fits = NULL, ...) {

  pfv <- .prob_for_varsel(fit)
  if(!is.null(pfv)) stop(pfv)
  if(is.null(fits)) fits <- cv_fit(fit) # change to use kfold
  if(!all(apply(fits, 1, function(fitrow, fit) {
    is.null(.prob_for_varsel(fitrow$fit)) && is.vector(fitrow$ind_test) &&
      max(fitrow$ind_test) <= nobs(fit) && all(fitrow$ind_test > 0)}, fit)))
    stop('fits does not have the correct form.')

  k <- nrow(fits)
  vars <- .extract_vars(fit)
  args <- .init_args(c(list(...), cv = T), vars, family(fit))

  d_test <- lapply(fits[,'ind_test'], function(inds, d_full) {
    list(x = d_full$x[inds,], y = d_full$y[inds],
         weights = d_full$weights[inds], offset = d_full$offset[inds])
  }, vars)

  # max number of variables to be projected
  args$nv <- min(c(sapply(c(list(fit), fits[,'fit']),
                          function(f) rankMatrix(get_x(f))), args$nv))

  msgs <- paste('Forward selection for the',
                c('full model.', paste0('fold number ', 1:k,'/',k,'.')))
  # perform the forward selection
  sel <- mapply(function(fit, d_test, msg, args) {
    if(args$verbose) print(msg)
    do.call(varsel, c(list(fit = fit, d_test = d_test), args))
  }, c(list(full = fit), fits[,'fit']), c(list(full = NA), d_test),
  msgs, MoreArgs = list(args))

  # combine cross validated results
  combcv <- function(x) as.list(data.frame(apply(x, 1, function(x) do.call(c, x))))
  d_cv <- combcv(simplify2array(d_test)[c('y', 'weights', 'offset'),])
  mu_cv <- combcv(simplify2array(sel['mu',-1]))
  lppd_cv <- combcv(simplify2array(sel['lppd',-1]))

  # evaluate performance on test data and
  # use bayesian bootstrap to get 95% credible intervals
  b_weights <- .gen_bootstrap_ws(length(d_cv$y), args$n_boot)
  nv <- c(1:(length(sel[['mu',2]])-1), nrow(vars$b)) - args$intercept
  stats <- rbind(sel[['stats',1]],
    .bootstrap_stats(sel[['mu',1]], sel[['lppd',1]], nv, vars, args$family_kl, b_weights, 'train'),
    .bootstrap_stats(mu_cv, lppd_cv, nv, d_cv, args$family_kl, b_weights, 'test'), make.row.names = F)

  # find out how many of cross-validated forward selection iterations select
  # the same variables as the forward selection that uses all the data.
  chosen_full <- sel[['chosen',1]]
  pctch <- mapply(function(var_ind, ind, arr, k) sum(arr[1:ind, ] == var_ind)/k,
                  chosen_full, seq_along(chosen_full),
                  MoreArgs = list(do.call(cbind, sel['chosen',-1]), k))

  res <- list(chosen = chosen_full, pctch = pctch, stats = stats, family = family(fit))
  if(args$clust) res$cl <- sel[['cl',1]]

  structure(res, class = 'varsel')
}
