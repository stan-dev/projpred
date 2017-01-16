#' Variable selection for generalized linear models
#'
#' Perform the projection predictive variable selection for a generalized
#' linear model fitted with rstanarm.
#' @param fit A \link[=stanreg-objects]{stanreg} object.
#' @param d_test A test dataset, which is used to evaluate model performance.
#' If not provided, training data is used.
#' @param method The method used in the variable selection. Possible options are
#' \code{'L1'} for L1-search and \code{'forward'} for forward selection.
#' Defaults to \code{'L1'}.
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
#'  \item{\code{stats}}{An array with statistics of the submodel performance.}
#'  \item{\code{family}}{A \code{\link{family}}-object.}
#' }
#'
#' @examples
#' \dontrun{
#' ### Usage with stanreg objects
#' fit <- stan_glm(y~x, binomial())
#' fit_v <- varsel(fit)
#' plot_varsel(fit_v)
#' }
#'

#' @export
varsel <- function(fit, d_test = NULL, method = 'L1', ns = 400L,
                   nv_max = NULL, intercept = NULL, verbose = F, ...) {
    UseMethod('varsel')
}

#' @export
varsel.stanreg <- function(fit, d_test = NULL, method = 'L1', ns = NULL, nc = NULL,
                           nv_max = NULL, intercept = NULL, verbose = F, ...) {
	
  # TODO, IMPLEMENT SENSIBILITY CHECKS FOR NS AND NC (INTO MISC.R) AND CALL THEM
  # TODO, FIGURE OUT HOW TO HANDLE THE TEST PREDICTIONS FOR FULL MODEL WHEN COEF_FULL ARE NOT AVAILABLE
	
  if (is.null(ns) && is.null(nc))
  	# by default, use one cluster for selection
  	nc <- 1
  
  .validate_for_varsel(fit)
  vars <- .extract_vars(fit)
  family_kl <- kl_helpers(family(fit))
  
  if(is.null(intercept)) intercept <- vars$intercept
  if(is.null(nv_max) || nv_max > NROW(vars$beta)) nv_max <- NROW(vars$beta)

  # training and test data
  d_train <- .get_traindata(fit)
  if (is.null(d_test)) {
  	d_test <- d_train
  	d_type <- 'train'
  } else {
  	d_test <- .fill_offset_and_weights(d_test)
  	d_type <- 'test'
  }
  
  # the reference distribution (or fit) used for selection
  p_full <- .get_refdist(fit, ns=ns, nc=nc)

  # perform the selection
  chosen <- select(method, p_full, d_train, family_kl, intercept, nv_max, verbose)

  # statistics for the selected submodels
  p_sub <- .get_submodels(chosen, c(0, seq_along(chosen)), family_kl, p_full, d_train, intercept)
  sub <- .get_sub_summaries(chosen, p_full, d_test, p_sub, family_kl, intercept)
  
  # TODO, FIGURE OUT HOW TO HANDLE THE TEST PREDICTIONS FOR FULL MODEL WHEN COEF_FULL ARE NOT AVAILABLE
  full <- NULL
  # full <- .get_full_summaries(p_full, d_test, list(alpha=vars$alpha, beta=vars$beta), family_kl, intercept)

  # ensure that after the new selection, there are no old projections in the fit structure
  fit$proj <- NULL
  
  # store the relevant fields into fit
  fit$varsel <- list(chosen = chosen,
                     kl = sapply(p_sub, function(x) x$kl),
                     d_test = c(d_test[c('y','weights')], type = d_type),
                     summaries = list(sub = sub, full = full),
                     family_kl = family_kl)

  fit
}

select <- function(method, p_full, d_train, family_kl, intercept, nv_max,
                   verbose) {
  #
  # Auxiliary function, performs variable selection with the given method,
  # and returns the variable ordering.
  #
  if (tolower(method) == 'l1') {
    chosen <- search_L1(p_full, d_train, family_kl, intercept, nv_max)
  } else if (tolower(method) == 'forward') {
    tryCatch(chosen <- search_forward(p_full, d_train, family_kl, intercept,
                                      nv_max, verbose),
             'error' = .varsel_errors)
  } else {
    stop(sprintf('Unknown search method: %s.', method))
  }
  return(chosen)
}




