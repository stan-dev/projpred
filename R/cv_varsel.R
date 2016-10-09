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
#'  \item{\code{ns = min(400, [number of draws])}}{
#'    Number of draws used in the variable selection.
#'    Cannot be larger than the number of draws in the full model.}
#'  \item{\code{nc = 0}}{
#'    If nonzero, a clustering with \code{nc} clusters is performed for
#'    the draws and the cluster centers are used in the variable selection
#'    instead of the actual draws.}
#'  \item{\code{nv = min(ncol(x) - 1, rankMatrix(x))}}{
#'    Maximum number of variables to be used in the projection (incl. intercept).
#'    Cannot be larger than \code{min(ncol(x) - 1, rankMatrix(x))}.}
#'  \item{\code{verbose = FALSE}}{
#'    If \code{verbose = TRUE}, prints information about the progress of the
#'    variable selection.}
#' }
#'
#' @return A list with class \code{'varsel'} containing the following elements:
#' \describe{
#'  \item{\code{chosen}}{The order in which the variables were added to the submodel.}
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
cv_varsel.stanreg <- function(fit, k_fold = NULL, ...) {

  .validate_for_varsel(fit)
  if(is.null(k_fold)) k_fold <- glmproj::kfold(fit, save_fits = T)

  if(!all(apply(k_fold$fits, 1, function(fits, fit) {
    .validate_for_varsel(fits$fit)
    is.vector(fits$omitted) && max(fits$omitted) <= nobs(fit) && all(fits$omitted > 0)
  }, fit))) stop('fits does not have the correct form.')

  k <- attr(k_fold, 'K')
  vars <- .extract_vars(fit)
  args <- .init_args(c(list(...), cv = T), vars, family(fit))

  d_test <- lapply(k_fold$fits[,'omitted'], function(omitted, d_full) {
    list(x = d_full$x[omitted,], y = d_full$y[omitted],
         weights = d_full$weights[omitted], offset = d_full$offset[omitted])
  }, vars)

  # max number of variables to be projected
  args$nv <- min(c(sapply(k_fold$fits[,'fit'], function(fit)
    rankMatrix(get_x(fit))), args$nv))

  msgs <- paste('Forward selection for the',
                c('full model.', paste0('fold number ', 1:k,'/',k,'.')))
  # perform the forward selection
  sel <- mapply(function(fit, d_test, msg, args) {
    print(msg)
    do.call(varsel, c(list(fit = fit, d_test = d_test), args))
  }, c(list(full = fit), k_fold$fits[,'fit']), c(list(full = NA), d_test),
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
  # the same variables as the forward selection with all the data.
  chosen_full <- sel[['chosen',1]]
  pctch <- mapply(function(var_ind, ind, arr, k) sum(arr[1:min(ind+0, nrow(arr)), ] == var_ind)/k,
                  chosen_full, seq_along(chosen_full),
                  MoreArgs = list(do.call(cbind, sel['chosen',-1]), k))

  res <- list(chosen = chosen_full, pctch = pctch, stats = stats, family = family(fit))
  if(args$clust) res$cl <- sel[['cl',1]]

  structure(res, class = 'varsel')
}
