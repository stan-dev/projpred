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
#' @return The original \link[=stanreg-objects]{stanreg} object augmented with an element 'varsel',
#' which is a list containing the following elements:
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
#' fits <- kfold(fit)
#' fit_v <- cv_varsel(fit, fits)
#' plot_varsel(fit_v)
#' }
#'

#' @export
cv_varsel <- function(fit, fits = NULL, method = 'L1', ns = 400L,
                      nv_max = NULL, intercept = NULL, verbose = F, ...) {
  UseMethod('cv_varsel')
}

#' @export
cv_varsel.stanreg <- function(fit, k_fold = NULL, method = 'L1', ns = 400L,
                              nv_max = NULL, intercept = NULL, verbose = F,
                              K = NULL, ...) {

  .validate_for_varsel(fit)
  if(is.null(k_fold)) {
    print(paste('k_fold not provided, performing 10-fold cross-validation',
                'for the stan model.'))
    #k_fold <- glmproj::kfold_(fit, save_fits = T)
    if(is.null(K)) K <- 10
    k_fold <- kfold_(fit, K, save_fits = T)
  }

  if(!all(apply(k_fold$fits, 1, function(fits, fit) {
    .validate_for_varsel(fits$fit)
    is.vector(fits$omitted) && max(fits$omitted) <= nobs(fit) && all(fits$omitted > 0)
  }, fit))) stop('k_fold does not have the correct form.')

  k <- attr(k_fold, 'K')
  msgs <- paste0(method, ' search for the fold number ', 1:k, '/', k, '.')

  vars <- .extract_vars(fit)
  vars_cv <- lapply(k_fold$fits[,'fit'], function(fit) .extract_vars(fit))
  family_kl <- kl_helpers(family(fit))
  if(is.null(intercept)) intercept <- vars$intercept
  if(is.null(nv_max) || nv_max > NROW(vars$beta)) nv_max <- NROW(vars$beta)

  d_test <- lapply(k_fold$fits[,'omitted'], function(omitted, d_full) {
    list(x = d_full$x[omitted,], y = d_full$y[omitted],
         weights = d_full$weights[omitted], offset = d_full$offset[omitted])
  }, vars)

  e <- get_data_and_parameters(vars, NA, intercept, ns, family_kl)
  e_cv <- mapply(function(vars, d_test) {
    get_data_and_parameters(vars, d_test, intercept, ns, family_kl)
  }, vars_cv, d_test, SIMPLIFY = F)

  paste(method, 'search for the full model.')
  chosen <- select(method, e$p_full, e$d_train, family_kl, intercept, nv_max,
                   verbose)
  chosen_cv <- mapply(function(e, msg) {
    print(msg)
    select(method, e$p_full, e$d_train, family_kl, intercept, nv_max, verbose)
  }, e_cv, msgs, SIMPLIFY = F)

  p_sub <- .get_submodels(chosen, c(0, seq_along(chosen)), family_kl, e$p_full,
                          e$d_train, intercept)
  p_sub_cv <- mapply(function(chosen, e) {
    .get_submodels(chosen, c(0, seq_along(chosen)), family_kl, e$p_full,
                   e$d_train, intercept)
  }, chosen_cv, e_cv, SIMPLIFY = F)

  sub_summaries <- .get_sub_summaries2(p_sub, chosen, e$d_eval, e$p_full,
                                       family_kl, intercept)
  full_summaries <- .get_full_summaries(e$d_eval, e$p_full, e$coef_full,
                                        family_kl, intercept)

  # extract and combine mu and lppd from the list
  hf <- function(x) as.list(do.call(rbind, x))

  sub_summaries_cv <- apply(
    mapply(function(p_sub, e) {
      lapply(.get_sub_summaries2(p_sub, chosen, e$d_eval, e$p_full,
                                 family_kl, intercept), data.frame)
    }, p_sub_cv, e_cv),
  1, hf)

  full_summaries_cv <- hf(lapply(e_cv, function(e) {
    data.frame(.get_full_summaries(e$d_eval, e$p_full, e$coef_full,
                                   family_kl, intercept))
  }))

  d_cv <- hf(lapply(d_test, function(d) {
    data.frame(d[c('y', 'weights', 'offset')])}))
  d_cv$x <- do.call(rbind, lapply(d_test, function(d) d$x))

  # evaluate performance on test data and
  # use bayesian bootstrap to get 95% credible intervals
  b_weights <- .get_bootstrap_ws(NROW(e$d_eval$x))
  metrics <- .bootstrap_metrics(sub_summaries, full_summaries, e$d_eval,
                               family_kl, intercept, b_weights, F)
  metrics_cv <- .bootstrap_metrics(sub_summaries_cv, full_summaries_cv, d_cv,
                                   family_kl, intercept, b_weights, T)
  kl <- .get_kl_array(p_sub)

  # find out how many of cross-validated forward selection iterations select
  # the same variables as the forward selection with all the data.
  chosen_array <- do.call(rbind, chosen_cv)
  pctch <- sapply(seq_along(chosen), function(ind) {
    sum(chosen_array[, 1:ind] == chosen[ind])/k
  })
  names(pctch) <- chosen

  fit$varsel <- list(chosen = chosen, pctch = pctch,
                     metrics = rbind(kl, metrics, metrics_cv))
  fit
}


