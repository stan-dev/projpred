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
cv_varsel <- function(fit, method = 'L1', cv_method = 'kfold', ns = 400L,
                      nv_max = NULL, intercept = NULL, verbose = F,
                      K = NULL, k_fold = NULL, ...) {
  UseMethod('cv_varsel')
}

#' @export
cv_varsel.stanreg <- function(fit,  method = 'L1', cv_method = 'kfold', ns = 400L,
                              nv_max = NULL, intercept = NULL, verbose = F,
                              K = NULL, k_fold = NULL, ...) {

  .validate_for_varsel(fit)
  vars <- .extract_vars(fit)
  if(is.null(intercept)) intercept <- vars$intercept
  if(is.null(nv_max) || nv_max > NROW(vars$beta)) nv_max <- NROW(vars$beta)

  print(paste('Performing', method, 'search for the full model.'))
  sel <- varsel(fit, NULL, method, ns, nv_max, intercept, verbose)$varsel

  if(tolower(cv_method) == 'kfold') {
    sel_cv <- kfold_varsel(fit, method, ns, nv_max, intercept, verbose, vars, K, k_fold)
  } else {
    stop(sprintf('Unknown cross-valdation method: %s.', method))
  }

  # find out how many of cross-validated iterations select
  # the same variables as the selection with all the data.
  pctch <- sapply(seq_along(sel$chosen), function(ind, chosen_array) {
    sum(chosen_array[, 1:ind] == sel$chosen[ind])/NROW(chosen_array)
  }, do.call(rbind, sel_cv$chosen_cv))

  fit$proj <- NULL
  fit$varsel <- c(sel[c('chosen', 'kl', 'family_kl')],
                  sel_cv[c('d_test', 'summaries')],
                  list(pctch = pctch))
  fit
}

kfold_varsel <- function(fit, method, ns, nv_max, intercept, verbose, vars,
                         K, k_fold) {
  # returns:
  #  - list of crossvalidated paths (chosen_cv),
  #  - list (d_test) with test outputs y, test weights and data type (string)
  #  - list with submodel and full model summaries

  # Construct the kfold-objects. The resulting list contains an element 'fits',
  # which is a K x 2 dimensional array. Each row corresponds to one of the K
  # folds. First column contains the rstanarm-objects and the second column
  # the indices of the omitted observations (aka test data).
  if(is.null(k_fold)) {
    if(is.null(K)) K <- 10
    print(paste0('k_fold not provided, performing ', K,
                 '-fold cross-validation for the stan model.'))
    k_fold <- kfold_(fit, K = K, save_fits = T)
  }
  family_kl <- kl_helpers(family(fit))

  # check that the fit-objects are valid for variable selection
  if(!all(apply(k_fold$fits, 1, function(fits, fit) {
    .validate_for_varsel(fits$fit)
    is.vector(fits$omitted) && max(fits$omitted) <= nobs(fit) && all(fits$omitted > 0)
  }, fit))) stop('k_fold does not have the correct form.')
  K <- attr(k_fold, 'K')

  # extract variables from each fit-object (stan-samples, x, y, etc.)
  # to a list of size K
  vars_cv <- lapply(k_fold$fits[,'fit'], .extract_vars)

  # List of size K with test data for each fold (note that vars is from
  # the full model, not from the cross-validated models).
  d_test <- lapply(k_fold$fits[,'omitted'], function(omitted) {
    list(x = vars$x[omitted,], y = vars$y[omitted],
         weights = vars$weights[omitted], offset = vars$offset[omitted])
  })

  # List of K elements, each containing d_train, p_full, etc. corresponding
  # to the corresponding fold.
  e_cv <- mapply(function(vars, d_test) {
    .get_data_and_parameters(vars, d_test, intercept, ns, family_kl)
  }, vars_cv, d_test, SIMPLIFY = F)

  # List of K elements, each a list of the variables selected for the
  # corresponding fold.
  msgs <- paste0(method, ' search for the fold number ', 1:K, '/', K, '.')
  chosen_cv <- mapply(function(e, msg) {
    print(msg)
    select(method, e$p_full, e$d_train, family_kl, intercept, nv_max, verbose)
  }, e_cv, msgs, SIMPLIFY = F)

  # Construct p_sub for each fold using .get_submodels.
  p_sub_cv <- mapply(function(chosen, e) {
    .get_submodels(chosen, c(0, seq_along(chosen)), family_kl, e$p_full,
                   e$d_train, intercept)
  }, chosen_cv, e_cv, SIMPLIFY = F)

  # Helper function extract and combine mu and lppd from K lists with each
  # n/K of the elements to one list with n elements
  hf <- function(x) as.list(do.call(rbind, x))

  # Apply some magic to manipulate the structure of the list so that instead of
  # list with K sub_summaries each containing n/K mu:s and lppd:s, we have only
  # one sub_summary-list that contains with all n mu:s and lppd:s.
  sub_cv <- apply(
    mapply(function(p_sub, e, chosen) {
      lapply(.get_sub_summaries(chosen, e$p_full, e$d_test, p_sub,
                                family_kl, intercept), data.frame)
    }, p_sub_cv, e_cv, chosen_cv),
    1, hf)

  # Same for the full model.
  full_cv <- hf(lapply(e_cv, function(e) {
    data.frame(.get_full_summaries(e$p_full, e$d_test, e$coef_full,
                                   family_kl, intercept))
  }))

  # Combine also the K separate test data sets into one list
  # with n y's and weights's.
  d_cv <- hf(lapply(d_test, function(d) {
    data.frame(d[c('y', 'weights')])}))

  list(chosen_cv = chosen_cv, d_test = c(d_cv, type = 'kfold'),
       summaries = list(sub = sub_cv, full = full_cv))
}

